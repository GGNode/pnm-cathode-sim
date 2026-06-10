"""
NMC532 开路电位 (Open-Circuit Voltage) 模块
=============================================

使用基于文献的查找表 + 线性插值。

数据来源: Xiang et al., J. Power Sources 241 (2013) 582-588
(NMC532 vs Li/Li+, GITT measurement at 25°C)

对应 DERIVATION.md 章节: §2.5, §5.2, §7.1
"""

import numpy as np

# ===== NMC532 OCV 查找表 (文献数据) =====
# x = c_s / c_s_max (锂嵌入度, 0=完全脱锂, 1=完全嵌锂)
# U = 开路电压 vs Li/Li+ [V]
#
# 数据来源: Xiang et al. (2013) Fig. 3, GITT measurement
# 物理: 脱锂越深 (x 越小), OCV 越高
_OCV_X = np.array([
    0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35,
    0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
    0.80, 0.85, 0.90, 0.95, 1.00,
])
_OCV_U = np.array([
    4.28, 4.15, 4.03, 3.95, 3.87, 3.81, 3.75, 3.71,
    3.68, 3.65, 3.63, 3.61, 3.59, 3.57, 3.55, 3.50,
    3.47, 3.40, 3.30, 3.10, 2.50,
])


def nmc532_ocv(soc: float | np.ndarray) -> float | np.ndarray:
    """
    计算 NMC532 的开路电压 (OCV)。

    使用查找表 + 线性插值覆盖完整的 [0, 1] 范围。
    对应 DERIVATION.md §2.5。

    Parameters
    ----------
    soc : float or array
        Lithium stoichiometric fraction, soc = c_s / c_s,max.

    Returns
    -------
    U_eq : float or array
        平衡电位 [V], 相对于 Li/Li+ 参考电极。
    """
    soc = np.asarray(soc, dtype=float)
    U = np.interp(soc, _OCV_X, _OCV_U)
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
    dU_dx = np.gradient(_OCV_U, _OCV_X)
    dU = np.interp(soc, _OCV_X, dU_dx)
    return float(dU) if dU.ndim == 0 else dU
