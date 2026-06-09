"""
NMC532 开路电位 (Open-Circuit Voltage) 模块
=============================================

物理背景
--------
本模块提供 NMC532 (LiNi₀.₅Mn₀.₃Co₀.₂O₂) 正极活性材料的平衡电位
(开路电压, OCV) 曲线。

开路电压 U(x) 是固相锂嵌入度 x = c_s / c_s_max 的函数:
- x = 0: 完全脱锂 (delithiated), U ≈ 4.28 V vs Li/Li+
- x = 0.5: 半嵌锂, U ≈ 3.73 V
- x = 1: 完全嵌锂 (lithiated), U ≈ 3.53 V

OCV 曲线在电池模型中的作用
--------------------------
在 Butler-Volmer 方程中, OCV 决定了平衡态:

    eta = phi_s - phi_e - U(x)

其中:
- phi_s: 固相电位 [V]
- phi_e: 电解质电位 [V]
- U(x): 平衡电位 vs Li/Li+ [V]
- eta: 过电位 [V]

当 eta = 0 时, 系统处于局部平衡: phi_s - phi_e = U(x)。

多项式拟合
----------
使用 6 阶多项式拟合实验数据:

    U(x) = a0 + a1*x + a2*x² + a3*x³ + a4*x⁴ + a5*x⁵ + a6*x⁶

多项式系数通过最小二乘法拟合 NMC532 半电池实验数据获得。

注意: 本模块中的 OCV 已经是相对于 Li/Li+ 参考电极的值,
因此在计算过电位时不需要额外的 Nernst 修正。
参见 DERIVATION.md §2.5: "If the empirical curve is already measured
vs Li/Li+ at reference electrolyte concentration, then electrolyte
concentration should not be double-counted."

导数 dU/dx
----------
Newton-Raphson 求解器需要 dU/dc_s 来构建 Jacobian:

    dU/dc_s = dU/dx * dx/dc_s = dU/dx * (1/c_s_max)

本模块提供 dU/dx, 调用者需乘以 1/c_s_max 得到 dU/dc_s。

对应 DERIVATION.md 章节: §2.5, §5.2, §7.1

参考文献
--------
Bernardi & Gojkovic, J. Electrochem. Soc. 166(12), A2485 (2019).
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np

# ===== NMC532 OCV 多项式系数 =====
# Khan et al. Eq. 2.19, corrected against the cited NMC532/NCM523 OCP
# source (Verma et al., JES 164 A3380, Eq. 5). Khan's rendered PDF drops
# the leading "3" in the soc^8 coefficient, printing -5520.41099; that
# transcription gives U(0.5) = 121 V. The cited source prints -35520.41099.
# The independent variable is the lithium state of charge/stoichiometric
# fraction in NMC532, i.e. soc = c_s / c_s,max, and the fit is calibrated
# over the measured GITT voltage window rather than the whole [0, 1] domain.
_OCV_COEFFS = [
    5744.862289,
    -35520.41099,
    95714.29862,
    -147364.5514,
    142718.3782,
    -90095.81521,
    37061.41195,
    -9578.599274,
    1409.309503,
    -85.31153081,
]


def nmc532_ocv(soc: float | np.ndarray) -> float | np.ndarray:
    """
    计算 NMC532 的开路电压 (OCV)。

    对应 DERIVATION.md §2.5: "For NMC532, an empirical U(x) with
    x = c_s/c_s,max is usually preferable."

    Implements Khan et al. Eq. 2.19 with the missing digit in the ``soc^8``
    coefficient restored from the cited OCP source. The paper writes the
    independent variable as ``soc_i``; the source paper defines this as the
    state of charge of lithium in NMC523/NMC532.

    Parameters
    ----------
    soc : float or array
        Lithium state-of-charge/stoichiometric fraction used in Eq. 2.19,
        soc = c_s / c_s,max.

    Returns
    -------
    U_eq : float or array
        平衡电位 [V], 相对于 Li/Li+ 参考电极。
        单位: 伏特 (Volt)。

    Notes
    -----
    The value is not clipped so this routine remains an exact transcription
    of the corrected Eq. 2.19. The polynomial is only physically meaningful
    over the calibrated GITT stoichiometry window.
    """
    soc = np.asarray(soc, dtype=float)
    U = np.polyval(_OCV_COEFFS, soc)
    U = U - 0.0003 * np.exp(7.657 * (soc**115))

    # 标量输入返回标量
    return float(U) if U.ndim == 0 else U


def ocv_derivative(soc: float | np.ndarray) -> float | np.ndarray:
    """
    计算 OCV 对嵌锂度 x 的导数 dU/dx, 用于 Newton-Raphson 耦合。

    对应 DERIVATION.md §5.2 中的浓度导数:
        d i/d c_s = (d i0/d c_s)(E_a - E_c) - B * (dU/d c_s)

    其中 dU/dc_s = dU/dx * (1/c_s_max)。

    物理意义:
    - dU/dx < 0: OCV 随嵌锂度增加而降低 (正常工作范围)
    - |dU/dx| 在 OCV 曲线平坦区域较小, 在陡峭区域较大
    - dU/dx 用于构建 Jacobian 矩阵, 实现 Newton-Raphson 迭代的收敛

    Parameters
    ----------
    x : float or array
        嵌锂度 (lithiation degree), x = c_s / c_s_max, 范围 [0, 1]。

    Returns
    -------
    dUdx : float or array
        OCV 对 x 的导数 [V]。
        通常为负值 (OCV 随 x 增加而降低)。

    Notes
    -----
    调用者需乘以 1/c_s_max 得到 dU/dc_s [V·m³/mol]:
        dU/dc_s = dUdx / c_s_max
    """
    soc = np.asarray(soc, dtype=float)
    deriv_coeffs = np.polyder(_OCV_COEFFS)
    dU = np.polyval(deriv_coeffs, soc)
    dU = dU - 0.0003 * np.exp(7.657 * (soc**115)) * 7.657 * 115.0 * soc**114

    return float(dU) if dU.ndim == 0 else dU
