"""
隔膜边界模型模块 (Collapsed 1D Separator Boundary Model)
=======================================================

物理背景
--------
本模块将隔膜 (separator) 简化为一维有限体积平板模型,
作为阴极电解质的阻抗/传输边界条件。

在完整的 Li 金属 | 隔膜 | NMC532 阴极半电池中:
- 隔膜位于阳极 (Li 金属) 和阴极 (NMC532) 之间
- 隔膜内充满电解液, 传导 Li+ 离子
- 隔膜厚度 ~25 µm (Khan et al. 2021)

本模块捕捉高倍率下的主要效应:
1. 欧姆压降 (Ohmic drop): I * L_sep / kappa_eff
2. 盐耗尽 (Salt depletion): 高电流下阴极侧 c_e 降低
3. Li 箔界面极化 (Li-foil BV overpotential): 可选

简化假设:
- 准稳态: 隔膜内浓度分布线性 (无积累)
- 一维: 仅沿厚度方向 (x) 变化
- 恒定孔隙率和 Bruggeman 迁曲度

模型参数来自 Khan et al. (2021) Table II。

对应 DERIVATION.md 章节: §4.1, §4.3

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from pnmcathode.physics.electrolyte import (
    electrolyte_diffusion_coefficient,
    electrolyte_ionic_conductivity,
)
from pnmcathode.physics.reaction import F, R


@dataclass(frozen=True)
class SeparatorParams:
    """
    一维隔膜边界模型参数。

    Parameters
    ----------
    enabled : bool
        是否启用隔膜模型。False 时返回零损失 (无限大电解质储液库)。
    thickness : float
        隔膜厚度 L_sep [m], 默认 25 µm。
    porosity : float
        隔膜孔隙率 epsilon_sep, 默认 0.39。
    bruggeman : float
        Bruggeman 迁曲度指数, 默认 1.5。
        有效传输系数: D_eff = D * epsilon^b, kappa_eff = kappa * epsilon^b。
    t_plus : float
        Li+ 迁移数 (transference number), 默认 0.363。
        物理含义: 电流中 Li+ 承载的份额。
    c_ref : float
        参考 (初始/储液库) 电解质浓度 [mol/m³], 默认 1200。
    c_floor : float
        最小允许浓度 [mol/m³], 防止数值奇异点, 默认 1。
    include_concentration_overpotential : bool
        是否包含浓度过电位 (Nernst 项)。
    include_li_foil_bv : bool
        是否包含 Li 箔 Butler-Volmer 过电位。
    i0_foil : float
        Li 箔交换电流密度 [A/m²], 默认 19。
    alpha_foil : float
        Li 箔 BV 对称因子, 默认 0.5。
    """

    enabled: bool = False
    thickness: float = 25e-6        # L_sep [m] 隔膜厚度
    porosity: float = 0.39          # epsilon_sep 隔膜孔隙率
    bruggeman: float = 1.5          # Bruggeman 迁曲度指数
    t_plus: float = 0.363           # Li+ 迁移数
    c_ref: float = 1200.0           # 参考浓度 [mol/m³]
    c_floor: float = 1.0            # 最小浓度 [mol/m³]
    include_concentration_overpotential: bool = True  # 浓度过电位开关
    include_li_foil_bv: bool = False  # Li 箔 BV 开关
    i0_foil: float = 19.0           # Li 箔交换电流密度 [A/m²]
    alpha_foil: float = 0.5         # Li 箔 BV 对称因子


@dataclass(frozen=True)
class SeparatorState:
    """
    隔膜边界计算结果。

    Attributes
    ----------
    c_cathode : float
        阴极侧电解质浓度 [mol/m³]。
    phi_e_cathode : float
        阴极侧电解质电位 [V] (含隔膜压降)。
    dphi_ohm : float
        欧姆压降 [V]。
    dphi_conc : float
        浓度过电位 [V] (Nernst 项)。
    eta_li : float
        Li 箔 BV 过电位 [V]。
    """

    c_cathode: float       # 阴极侧电解质浓度 [mol/m³]
    phi_e_cathode: float   # 阴极侧电解质电位 [V]
    dphi_ohm: float        # 欧姆压降 [V]
    dphi_conc: float       # 浓度过电位 [V]
    eta_li: float          # Li 箔 BV 过电位 [V]


def separator_boundary(
    I_app: float,
    T: float,
    params: SeparatorParams,
    n_iter: int = 3,
) -> SeparatorState:
    """
    计算隔膜传输模型的阴极侧电解质边界条件。

    算法:
    1. 准稳态迭代求解隔膜内平均浓度 (n_iter 次)
    2. 计算欧姆压降: dphi_ohm = I * L_sep / kappa_eff
    3. 计算浓度过电位 (Nernst 项): dphi_conc = (RT/F)(1-t+) ln(c_ref/c_cathode)
    4. 可选: Li 箔 BV 过电位

    Parameters
    ----------
    I_app : float
        施加电流密度 [A/m²] (阳极约定)。放电: I_app < 0。
    T : float
        温度 [K]。
    params : SeparatorParams
        隔膜模型参数。
    n_iter : int
        准稳态迭代次数 (默认 3)。

    Returns
    -------
    SeparatorState
        阴极侧浓度、电位和各项损失。
    """
    # 未启用时返回零损失 (无限大储液库)
    if not params.enabled:
        return SeparatorState(
            c_cathode=params.c_ref,
            phi_e_cathode=0.0,
            dphi_ohm=0.0,
            dphi_conc=0.0,
            eta_li=0.0,
        )

    I_dis = abs(I_app)  # 放电电流密度 (正值) [A/m²]
    c_sep_mean = params.c_ref  # 隔膜内平均浓度 [mol/m³]

    # ===== 准稳态迭代求解隔膜内平均浓度 =====
    # 隔膜内盐浓度梯度由 Li+ 迁移和扩散平衡决定:
    #   dc/dx = I * (1 - t+) / (F * D_eff)
    # 阴极侧浓度 = c_ref - dc_sep
    for _ in range(n_iter):
        # 有效扩散系数 [m²/s]: D_eff = D_e * epsilon^b (Bruggeman)
        D_eff = electrolyte_diffusion_coefficient(c_sep_mean, T) * params.porosity ** params.bruggeman
        # 有效离子电导 [S/m]: kappa_eff = kappa * epsilon^b
        k_eff = electrolyte_ionic_conductivity(c_sep_mean, T) * params.porosity ** params.bruggeman

        # 浓度降 [mol/m³]: dc = I * (1-t+) * L / (F * D_eff)
        dc_sep = I_dis * (1.0 - params.t_plus) * params.thickness / (F * D_eff)
        c_sep_cathode = max(params.c_floor, params.c_ref - dc_sep)
        # 用平均浓度更新传输系数 (准稳态近似)
        c_sep_mean = 0.5 * (params.c_ref + c_sep_cathode)

    # ===== 欧姆压降 [V] =====
    # dphi_ohm = I * L_sep / kappa_eff
    dphi_ohm = I_dis * params.thickness / k_eff

    # ===== 浓度过电位 (Nernst 项) [V] =====
    # dphi_conc = (RT/F) * (1-t+) * ln(c_ref / c_cathode)
    c_cathode_safe = max(c_sep_cathode, params.c_floor)
    if params.include_concentration_overpotential:
        dphi_conc = (R * T / F) * (1.0 - params.t_plus) * math.log(params.c_ref / c_cathode_safe)
    else:
        dphi_conc = 0.0

    # ===== Li 箔 BV 过电位 [V] =====
    # eta_li = (2RT/F) * asinh(I / (2*i0))
    # 禁用时 = 理想 Li 参考电极 (无极化)
    if params.include_li_foil_bv and I_dis > 0.0:
        eta_li = (2.0 * R * T / F) * math.asinh(I_dis / (2.0 * params.i0_foil))
    else:
        eta_li = 0.0

    # 阴极侧电解质参考电位 [V] (放电时为负值)
    phi_e_cathode = -(dphi_ohm + dphi_conc + eta_li)

    return SeparatorState(
        c_cathode=c_cathode_safe,
        phi_e_cathode=phi_e_cathode,
        dphi_ohm=dphi_ohm,
        dphi_conc=dphi_conc,
        eta_li=eta_li,
    )
