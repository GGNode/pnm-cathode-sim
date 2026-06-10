"""
NMC532 开路电位 (Open-Circuit Voltage) 模块
=============================================

使用 Khan et al. (2021) Eq. 2.19 多项式拟合。
多项式在 GITT 校准窗口内有效; x<0.25 时外推到非物理负值,
因此加 clip 到 [3.0, 4.5]V。

论文: J. Electrochem. Soc. 168 (2021) 070534
来源: PAPER_REFERENCE.md

对应 DERIVATION.md 章节: §2.5, §5.2, §7.1
"""

import numpy as np

# ===== Khan et al. Eq. 2.19 多项式系数 =====
# 注意: soc^8 系数为 -35520.41099 (论文勘误: 原文打印 -5520.41099)
_KHAN_COEFFS = [
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

    Khan et al. Eq. 2.19 多项式 + exp 修正。
    Clip 到 [3.0, 4.5]V 防止 x<0.25 时多项式外推到负值。

    Parameters
    ----------
    soc : float or array
        Lithium stoichiometric fraction, soc = c_s / c_s,max.

    Returns
    -------
    U_eq : float or array
        平衡电位 [V], 相对于 Li/Li+ 参考电极.

    验证值:
    - U(0.35) ≈ 4.30V
    - U(0.50) ≈ 4.00V
    - U(0.90) ≈ 3.67V
    - U(1.00) ≈ 2.94V → clip to 3.0V
    """
    soc = np.asarray(soc, dtype=float)
    U = np.polyval(_KHAN_COEFFS, soc)
    U = U - 0.0003 * np.exp(7.657 * (soc**115))
    # Clip: 多项式在 x<0.25 产生非物理负值, x>0.95 急剧下降
    U = np.clip(U, 3.0, 4.5)
    return float(U) if U.ndim == 0 else U


def ocv_derivative(soc: float | np.ndarray) -> float | np.ndarray:
    """
    计算 OCV 对嵌锂度 x 的导数 dU/dx, 用于 Newton-Raphson 耦合。

    Parameters
    ----------
    soc : float or array
        Lithium stoichiometric fraction.

    Returns
    -------
    dU_dx : float or array
        dU/dx [V].
    """
    soc = np.asarray(soc, dtype=float)
    deriv_coeffs = np.polyder(_KHAN_COEFFS)
    dU = np.polyval(deriv_coeffs, soc)
    # exp correction derivative
    dU = dU - 0.0003 * 7.657 * 115 * (soc**114) * np.exp(7.657 * (soc**115))
    return float(dU) if dU.ndim == 0 else dU
