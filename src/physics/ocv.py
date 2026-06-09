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
# 拟合 NMC532 半电池实验数据:
#   x:  0.0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0
#   U: 4.28  4.15  4.05  3.93  3.82  3.73  3.65  3.60  3.57  3.55  3.53
#
# 多项式: U(x) = a0 + a1*x + a2*x² + a3*x³ + a4*x⁴ + a5*x⁵ + a6*x⁶
# np.polyval 按降幂排列: coeffs[0]*x^6 + coeffs[1]*x^5 + ... + coeffs[6]
_OCV_COEFFS = [
    6.290850,    # a6 (x^6 项系数)
    -21.596908,  # a5 (x^5 项系数)
    27.066679,   # a4 (x^4 项系数)
    -14.688871,  # a3 (x^3 项系数)
    3.689029,    # a2 (x^2 项系数)
    -1.509975,   # a1 (x^1 项系数)
    4.279404,    # a0 (常数项, 对应 U(x=0) ≈ 4.28 V)
]


def nmc532_ocv(x: float | np.ndarray) -> float | np.ndarray:
    """
    计算 NMC532 的开路电压 (OCV)。

    对应 DERIVATION.md §2.5: "For NMC532, an empirical U(x) with
    x = c_s/c_s,max is usually preferable."

    OCV 曲线特征:
    - 单调递减: x 越大 (越嵌锂), 电压越低
    - x=0 (空): U ≈ 4.28 V (高电压, 脱锂态)
    - x=1 (满): U ≈ 3.53 V (低电压, 嵌锂态)
    - 放电过程中 x 增大, U 减小

    Parameters
    ----------
    x : float or array
        嵌锂度 (lithiation degree), x = c_s / c_s_max, 范围 [0, 1]。
        x = 0: 完全脱锂 (空), x = 1: 完全嵌锂 (满)。

    Returns
    -------
    U_eq : float or array
        平衡电位 [V], 相对于 Li/Li+ 参考电极。
        单位: 伏特 (Volt)。

    Notes
    -----
    - 使用 np.clip 将 x 限制在 [1e-6, 1-1e-6] 以避免多项式外推。
    - 结果被裁剪到 [3.0, 4.3] V 的物理合理范围。
    - 参见 DERIVATION.md §6.2: "Empirical OCV fits often diverge
      or become invalid outside calibrated x."
    """
    x = np.asarray(x, dtype=float)

    # 裁剪 x 到安全范围, 防止多项式在外推区域发散
    x = np.clip(x, 1e-6, 1.0 - 1e-6)

    # 使用 numpy 多项式求值: U = a6*x^6 + a5*x^5 + ... + a0
    U = np.polyval(_OCV_COEFFS, x)

    # 裁剪到 NMC532 的物理电压范围
    # NMC532 的 OCV 不会低于 ~3.5V 或高于 ~4.3V vs Li/Li+
    U = np.clip(U, 3.0, 4.3)

    # 标量输入返回标量
    return float(U) if U.ndim == 0 else U


def ocv_derivative(x: float | np.ndarray) -> float | np.ndarray:
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
    x = np.asarray(x, dtype=float)

    # 裁剪到安全范围
    x = np.clip(x, 1e-6, 1.0 - 1e-6)

    # 计算多项式的导数系数
    # 如果 U = a6*x^6 + a5*x^5 + ... + a0
    # 则 dU/dx = 6*a6*x^5 + 5*a5*x^4 + ... + a1
    deriv_coeffs = np.polyder(_OCV_COEFFS)

    # 求导数值
    dU = np.polyval(deriv_coeffs, x)

    return float(dU) if dU.ndim == 0 else dU
