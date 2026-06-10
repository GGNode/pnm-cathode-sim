"""
电解质传输物性模块 (Electrolyte Transport Properties)
=====================================================

物理背景
--------
本模块实现锂离子电池电解液的传输物性相关性 (correlations)。

电解液是 Li+ 离子在阴极内传输的介质。关键传输物性:
1. Li+ 扩散系数 D_e(c_e, T) [m²/s]: 控制浓度扩散通量
2. 离子电导率 kappa(c_e, T) [S/m]: 控制电流传输

这两个物性都依赖于电解质浓度 c_e 和温度 T。
模型中使用 Khan et al. (2021) Table II 的经验相关性。

浓度单位转换
-----------
论文 Table II 中的相关性使用 mol/L 作为浓度单位,
而模型状态变量使用 mol/m³。
    1 mol/L = 1000 mol/m³

本模块内部处理单位转换, 调用方只需传入 mol/m³。

扩散系数 (Diffusion Coefficient)
---------------------------------
Table II 给出:
    D_Li+ = 10^(-4.43 - 54.0/(T - 229 - 5c2) - 0.22c2)
其中 c2 为 mol/L。

原始公式给出 cm²/s, 乘以 1e-4 转换为 m²/s。
典型值: ~1e-10 m²/s (1M, 298K)。

离子电导率 (Ionic Conductivity)
-------------------------------
Table II footnote (2):
    kappa = c2 * (-10.5 + 0.0740*T - 6.96e-5*T² + ...)²
原始公式给出 mS/cm, 乘以 0.1 转换为 S/m。
典型值: ~1 S/m (1M, 298K)。

对应 DERIVATION.md 章节: §2.2, §3.1

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

from __future__ import annotations

import numpy as np


def _mol_m3_to_mol_l(c_e: float | np.ndarray) -> np.ndarray:
    """
    将电解质浓度从 mol/m³ 转换为 mol/L。

    论文 Table II 的相关性使用 mol/L, 而模型状态使用 mol/m³。
    转换: c [mol/L] = c [mol/m³] / 1000

    Parameters
    ----------
    c_e : float or ndarray
        电解质浓度 [mol/m³]。

    Returns
    -------
    c2 : ndarray
        电解质浓度 [mol/L]。
    """
    return np.asarray(c_e, dtype=float) / 1000.0


def electrolyte_diffusion_coefficient(
    c_e: float | np.ndarray,
    T: float = 303.0,
) -> float | np.ndarray:
    """
    计算 Li+ 在电解液中的扩散系数 D_e(c_e, T)。

    Khan et al. Table II 经验公式:
        D_Li+ = 10^(-4.43 - 54.0/(T - 229 - 5c2) - 0.22c2)

    其中 c2 为 mol/L 浓度。原始公式给出 cm²/s,
    乘以 1e-4 转换为 m²/s。

    Parameters
    ----------
    c_e : float or ndarray
        电解质 Li+ 浓度 [mol/m³]。典型值: 1000~1200。
    T : float
        温度 [K], 默认 303 K (30°C)。

    Returns
    -------
    D_e : float or ndarray
        Li+ 扩散系数 [m²/s]。
        典型值: ~1e-10 m²/s (1M, 298K)。
    """
    c2 = _mol_m3_to_mol_l(c_e)
    # 指数: -4.43 - 54/(T - 229 - 5c) - 0.22c
    exponent = -4.43 - 54.0 / (T - 229.0 - 5.0 * c2) - 0.22 * c2
    value_cm2_s = 10.0 ** exponent
    # cm²/s → m²/s (×1e-4)
    value = value_cm2_s * 1e-4
    return float(value) if value.ndim == 0 else value


def electrolyte_ionic_conductivity(
    c_e: float | np.ndarray,
    T: float = 303.0,
) -> float | np.ndarray:
    """
    计算电解液离子电导率 kappa(c_e, T)。

    Khan et al. Table II footnote (2) 经验公式:
        kappa = c2 * inner²
        inner = -10.5 + 0.0740*T - 6.96e-5*T²
                + 0.668*c2 - 0.0178*c2*T + 2.80e-5*c2*T²
                + 0.494*c2² - 8.86e-4*c2²*T

    其中 c2 为 mol/L。原始公式给出 mS/cm,
    乘以 0.1 转换为 S/m。

    Parameters
    ----------
    c_e : float or ndarray
        电解质 Li+ 浓度 [mol/m³]。
    T : float
        温度 [K], 默认 303 K。

    Returns
    -------
    kappa : float or ndarray
        离子电导率 [S/m]。
        典型值: ~1 S/m (1M, 298K)。
    """
    c2 = _mol_m3_to_mol_l(c_e)
    inner = (
        -10.5
        + 0.0740 * T
        - 6.96e-5 * T**2
        + 0.668 * c2
        - 0.0178 * c2 * T
        + 2.80e-5 * c2 * T**2
        + 0.494 * c2**2
        - 8.86e-4 * c2**2 * T
    )
    value_mS_cm = c2 * inner**2
    # mS/cm → S/m (×0.1)
    value = value_mS_cm * 0.1
    return float(value) if value.ndim == 0 else value
