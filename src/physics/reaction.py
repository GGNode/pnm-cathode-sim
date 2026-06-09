"""
Butler-Volmer 电极反应动力学模块 (Reaction Kinetics)
=====================================================

物理背景
--------
本模块实现电解质/NMC 界面处的嵌锂/脱锂反应动力学。

对应嵌锂反应 (intercalation reaction):

    Li_s  <->  Li+_e + e-_s

即固相中的锂原子 ↔ 电解质中的锂离子 + 固相中的电子。

对应的 Butler-Volmer 方程给出界面反应电流密度:

    i_r = i0 * [exp(alpha_a * F * eta / (R*T)) - exp(-alpha_c * F * eta / (R*T))]

其中过电位 (overpotential) 定义为:

    eta = phi_s - phi_e - U_eq(c_s, c_e)

符号约定 (Sign Convention) — 阳极约定 (Anodic Convention):
---------------------------------------------------------
- eta > 0  → i_r > 0: 阳极反应 (脱锂/氧化, deintercalation/oxidation)
- eta < 0  → i_r < 0: 阴极反应 (嵌锂/还原, intercalation/reduction)
- eta = 0  → i_r = 0: 平衡态 (equilibrium)

在放电 (discharge/lithiation) 过程中:
- 阴极发生还原反应, i_r < 0
- 电解质中的 Li+ 被消耗 (c_e 减小)
- 固相中的 Li 浓度增加 (c_s 增大)

交换电流密度 (Exchange Current Density):
---------------------------------------
    i0 = F * k0 * ce^alpha_a * (cs_max - cs)^alpha_a * cs^alpha_c

其中:
- k0: 反应速率常数 [m^(2.5) / (mol^0.5 · s)]
- ce: 电解质 Li+ 浓度 [mol/m³]
- cs: 固相表面锂浓度 [mol/m³]
- cs_max: NMC 最大锂浓度 [mol/m³] (NMC532: 48900)
- alpha_a, alpha_c: 阳极/阴极传递系数 (通常各为 0.5)

数值稳定性
----------
指数项 exp(alpha * F * eta / (RT)) 在大过电位时可能溢出。
298 K 时, RT/F ≈ 25.7 mV, eta = 1V 对应指数约 19.5。
本模块使用 np.clip 将指数参数限制在 [-500, 500] 范围内。

对应 DERIVATION.md 章节: §1.3, §2.5, §5.2, §6.1

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np

# ===== 物理常数 =====
F = 96485.3329   # 法拉第常数 (Faraday constant) [C/mol]
                 # 定义: 1 mol 电子所携带的电荷量
R = 8.314462     # 摄氏气体常数 (Gas constant) [J/(mol·K)]


def butler_volmer(
    i0: float,
    eta: float | np.ndarray,
    T: float = 298.15,
    alpha_a: float = 0.5,
    alpha_c: float = 0.5,
) -> float | np.ndarray:
    """
    计算 Butler-Volmer 反应电流密度。

    对应 DERIVATION.md §2.5 的 Butler-Volmer 动力学方程:

        i_r = i0 * [exp(alpha_a * F * eta / (R*T))
                   - exp(-alpha_c * F * eta / (R*T))]

    其中:
        - 第一项 exp(alpha_a * F * eta / (R*T)) 为阳极 (氧化) 分支
        - 第二项 exp(-alpha_c * F * eta / (R*T)) 为阴极 (还原) 分支

    符号约定 (Anodic Convention):
        - eta > 0 → i_r > 0 (脱锂, deintercalation)
        - eta < 0 → i_r < 0 (嵌锂, lithiation / discharge)

    当 alpha_a = alpha_c = 0.5 时, 方程可简化为:
        i_r = 2 * i0 * sinh(F * eta / (2*R*T))

    Parameters
    ----------
    i0 : float
        交换电流密度 (exchange current density) [A/m²]。
        反映反应的内在速率, 由 exchange_current_density() 计算。
    eta : float or array
        过电位 (overpotential) [V]。
        eta = phi_s - phi_e - U_eq。
        正值 → 阳极脱锂, 负值 → 阴极嵌锂。
    T : float
        温度 [K], 默认 298.15 K (25°C)。
    alpha_a : float
        阳极传递系数 (anodic transfer coefficient), 无量纲。
        物理含义: 反应活化能中电场对阳极反应的贡献比例。
        典型值: 0.5 (对称势垒)。
    alpha_c : float
        阴极传递系数 (cathodic transfer coefficient), 无量纲。
        典型值: 0.5。
        注意: alpha_a + alpha_c = 1 (对于单电子转移反应)。

    Returns
    -------
    I_rxn : float or array
        反应电流密度 [A/m²]。
        正值 = 阳极 (脱锂), 负值 = 阴极 (嵌锂)。

    Notes
    -----
    数值稳定性: 使用 np.clip 将指数参数限制在 [-500, 500]。
    exp(500) ≈ 1.4e217 在 float64 范围内; exp(-500) ≈ 0。
    参见 DERIVATION.md §6.1 "Exponential Overflow in Butler-Volmer"。
    """
    eta = np.asarray(eta, dtype=float)

    # f = F/(R*T) 是热电压的倒数 [V^{-1}]
    # 在 298.15 K 时, f ≈ 38.92 V^{-1}, 即 RT/F ≈ 25.7 mV
    f = F / (R * T)

    # 计算指数参数并裁剪以防止溢出
    # arg_a = alpha_a * F * eta / (R*T): 阳极分支指数
    # arg_c = -alpha_c * F * eta / (R*T): 阴极分支指数
    # 裁剪到 [-500, 500] 保证 exp() 不溢出
    arg_a = np.clip(alpha_a * f * eta, -500.0, 500.0)
    arg_c = np.clip(-alpha_c * f * eta, -500.0, 500.0)

    # Butler-Volmer 方程: i_r = i0 * [exp(arg_a) - exp(arg_c)]
    result = i0 * (np.exp(arg_a) - np.exp(arg_c))

    # 如果输入是标量, 返回标量
    return float(result) if result.ndim == 0 else result


def exchange_current_density(
    k0: float,
    ce: float,
    cs: float,
    cs_max: float = 48900.0,
    alpha_a: float = 0.5,
    alpha_c: float = 0.5,
) -> float:
    """
    计算 Li 嵌入反应的交换电流密度。

    对应 DERIVATION.md §2.5 的交换电流密度模型:

        i0 = F * k0 * ce^alpha_a * (cs_max - cs)^alpha_a * cs^alpha_c

    物理含义:
    - ce^alpha_a: 电解质 Li+ 浓度对反应速率的贡献
      (更高的 Li+ 浓度 → 更快的反应)
    - (cs_max - cs)^alpha_a: 固相空位浓度的贡献
      (更多空位 → 更容易嵌锂)
    - cs^alpha_c: 固相锂浓度的贡献
      (更多固相锂 → 更容易脱锂)
    - F: 法拉第常数, 将摩尔速率转换为电流密度

    指数 gamma_e = alpha_a, gamma_v = alpha_a, gamma_s = alpha_c
    对应 DERIVATION.md §2.5 中的 "common Li-ion choice"。

    单位分析:
        [F] * [k0] * [ce]^alpha_a * [cs_max-cs]^alpha_a * [cs]^alpha_c
        = [C/mol] * [m^(2.5)/(mol^0.5·s)] * [mol/m³]^(alpha_a+alpha_a+alpha_c)
        = [C/mol] * [m^(2.5)/(mol^0.5·s)] * [mol/m³]^1
        (当 alpha_a + alpha_c = 1, 且 alpha_a = 0.5 时)
        = [A/m²]

    Parameters
    ----------
    k0 : float
        反应速率常数 (rate constant) [m^(2.5) / (mol^0.5 · s)]。
        反映反应的内在动力学速率。
    ce : float
        电解质 Li+ 浓度 [mol/m³]。典型值: 1000 (0C) ~ 1200 (1M)。
    cs : float
        固相表面锂浓度 [mol/m³]。范围: (0, cs_max)。
    cs_max : float
        NMC 最大锂浓度 [mol/m³]。NMC532: 48900 mol/m³。
    alpha_a : float
        阳极传递系数。
    alpha_c : float
        阴极传递系数。

    Returns
    -------
    i0 : float
        交换电流密度 [A/m²]。始终为正值。

    Notes
    -----
    为避免数值奇异点, ce 和 cs 被限制在远离零的小正数。
    参见 DERIVATION.md §5.2 "concentration derivatives are singular at
    c_e = 0, c_s = 0, and c_s = c_s,max"。
    """
    # 避免 ce = 0 导致的奇异点 (i0 会为零或导数发散)
    ce = max(ce, 1e-10)
    # 限制 cs 在 (0, cs_max) 之间, 避免 cs=0 或 cs=cs_max 的奇异点
    cs = np.clip(cs, 1e-10, cs_max - 1e-10)

    # 交换电流密度公式
    # 参见 DERIVATION.md §2.5 Eq. (i0)
    return k0 * F * (ce ** alpha_a) * ((cs_max - cs) ** alpha_a) * (cs ** alpha_c)
