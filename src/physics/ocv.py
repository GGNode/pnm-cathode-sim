"""
NMC532 开路电位 (Open-Circuit Voltage) 模块
=============================================

使用查找表 + 线性插值。Khan et al. Eq. 2.19 多项式仅在 x∈[0.3, 0.95]
校准有效; x<0.3 时多项式产生非物理负值, x>0.95 时急剧下降。
查找表覆盖完整 [0, 1] 范围。

数据来源:
- x ∈ [0.30, 0.95]: Khan et al. Eq. 2.19 多项式 (校准范围)
- x ∈ [0.00, 0.30]: 基于 NMC532 高电压平台的文献外推
  (Xiang et al., J. Power Sources 241 (2013) 582-588)
- x ∈ [0.95, 1.00]: 多项式端点修正

对应 DERIVATION.md 章节: §2.5, §5.2, §7.1
"""

import numpy as np

# ===== NMC532 OCV 查找表 =====
# x = c_s / c_s_max (锂嵌入度, 0=完全脱锂, 1=完全嵌锂)
# U = 开路电压 vs Li/Li+ [V]
#
# 物理: 脱锂越深 (x 越小), OCV 越高 (更多能量可释放)
_OCV_X = np.array([
    0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35,
    0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
    0.80, 0.85, 0.90, 0.95, 1.00,
])
_OCV_U = np.array([
    4.50, 4.42, 4.35, 4.30, 4.38, 4.38, 4.38, 4.30,  # x=0.00-0.35
    4.19, 4.09, 4.00, 3.91, 3.84, 3.80, 3.77, 3.74,  # x=0.40-0.75
    3.72, 3.70, 3.67, 3.64, 3.00,                      # x=0.80-1.00
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
