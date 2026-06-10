"""
固相锂扩散模块 (Solid-Phase Li Diffusion)
==========================================

物理背景
--------
本模块实现 NMC532 活性材料颗粒内部的锂扩散。

NMC532 颗粒中的锂扩散遵循 Fick 定律:

    ∂c_s/∂t = ∇·(D_s(c_s, T) ∇c_s)

其中:
- c_s: 固相锂浓度 [mol/m³]
- D_s: 锂扩散系数 [m²/s], 依赖于 c_s 和 T

在球坐标系中 (假设球形颗粒):

    ∂c_s/∂t = (1/r²) * ∂/∂r (r² * D_s * ∂c_s/∂r)

扩散系数
--------
NMC532 的锂扩散系数 D_s 与浓度和温度有关:

    D_s(c_s, T) = D_ref * exp(-E_a/R * (1/T - 1/T_ref)) * f(x)

其中 x = c_s / c_s_max, f(x) 描述浓度依赖性。

典型值: D_ref ~ 1e-14 m²/s (298 K, 50% SoC)
        E_a ~ 30 kJ/mol (活化能)
        文献报道范围: 1e-16 ~ 1e-12 m²/s

有限体积离散化
--------------
将球形颗粒离散为 N 个同心壳层 (有限体积):

    Shell i: 内半径 r_i = i*dr, 外半径 r_(i+1) = (i+1)*dr
    壳层体积: V_i = (4/3)*π*(r_(i+1)³ - r_i³)
    壳层厚度: dr = R_p / N

壳层之间的扩散通量使用两点通量近似:

    J_{i→i+1} = -D_s * (c_s[i+1] - c_s[i]) / dr  [mol/(m²·s)]

界面面积: A = 4*π*r²

表面边界条件:
    flux_surface = -I_rxn / F  [mol/(m²·s)]
    正值 = 锂进入颗粒 (嵌锂), 负值 = 锂离开颗粒 (脱锂)

质量守恒
--------
有限体积离散化保证质量守恒:

    Σ_i (dc_s[i]/dt * V_i) = flux_surface * A_surface

即: 颗粒内总锂量的变化率 = 表面通量 × 表面积

符号约定
--------
- flux_surface > 0: 锂进入颗粒 (嵌锂, lithiation)
  对应放电过程中的阴极反应 (i_r < 0)
- flux_surface < 0: 锂离开颗粒 (脱锂, deintercalation)

    flux_surface = -I_rxn / F

当 I_rxn < 0 (阴极反应) 时, flux_surface > 0 (嵌锂)。

对应 DERIVATION.md 章节: §2.3, §3.4, §7.5

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np


def nmc532_diffusion_coefficient(c_s: float, T: float = 303.0) -> float:
    """
    计算 NMC532 锂扩散系数 D_s(c_s, T)。

    对应 DERIVATION.md §2.3: "Solid lithium in NMC is neutral
    intercalated lithium. Its flux is Fickian: N_solid = -D_s(c_s,T) grad c_s."

    Table II 脚注 (1) 给出了浓度依赖的 NMC532 扩散系数,
    以 soc 的 10 为底多项式表示:

        D_s = 10^(-2319soc^10 + 6642soc^9 - 5269soc^8 - 3319soc^7
                 + 10038soc^6 - 9806soc^5 + 5817soc^4 - 2286soc^3
                 + 575.3soc^2 - 83.16soc - 9.292)

    ``T`` 参数仅为 API 兼容性而保留, 该经验公式不含温度项。

    Parameters
    ----------
    c_s : float
        固相锂浓度 [mol/m³]。范围: (0, c_s_max=48900)。
    T : float
        温度 [K]。保留以兼容调用方；论文中该固相公式不含温度项。

    Returns
    -------
    D_s : float
        锂扩散系数 [m²/s]。
        典型值: ~1e-15 m²/s near 50% SoC.

    Notes
    -----
    参见 DERIVATION.md §3.4: "For nonlinear D_s(c_s,T), K_mn^s must
    be evaluated at t^{n+1} in a fully implicit Newton solve."
    """
    c_s_max = 48900.0
    soc = np.clip(c_s / c_s_max, 1e-9, 1.0 - 1e-9)
    exponent = (
        -2319.0 * soc**10
        + 6642.0 * soc**9
        - 5269.0 * soc**8
        - 3319.0 * soc**7
        + 10038.0 * soc**6
        - 9806.0 * soc**5
        + 5817.0 * soc**4
        - 2286.0 * soc**3
        + 575.3 * soc**2
        - 83.16 * soc
        - 9.292
    )
    return float(10.0 ** exponent)


def discretize_spherical_particle(R_p: float, N: int) -> tuple[float, np.ndarray]:
    """
    将球形颗粒离散为 N 个同心壳层 (有限体积法)。

    对应 DERIVATION.md §7.5 中的固体扩散离散化。

    壳层几何:
    - Shell 0 (核心):  r ∈ [0, dr]
    - Shell 1:          r ∈ [dr, 2*dr]
    - ...
    - Shell N-1 (表面): r ∈ [(N-1)*dr, N*dr=R_p]

    每个壳层的体积:
        V_i = (4/3)*π*(r_outer³ - r_inner³)
            = (4/3)*π*((i+1)³ - i³)*dr³

    所有壳层体积之和等于球的总体积:
        Σ V_i = (4/3)*π*R_p³

    Parameters
    ----------
    R_p : float
        颗粒半径 (particle radius) [m]。
        NMC532 典型值: 2.5~5 μm (2.5e-6 ~ 5e-6 m)。
    N : int
        壳层数量。更多的壳层 → 更高的空间分辨率, 但计算成本更高。
        典型值: 10~50。

    Returns
    -------
    dr : float
        壳层厚度 [m]。dr = R_p / N。
    volumes : ndarray of shape (N,)
        每个壳层的体积 [m³]。
        单位: m³。
    """
    dr = R_p / N  # 壳层厚度 [m]

    # 每个壳层的内外半径
    # Shell i: 内半径 = i*dr, 外半径 = (i+1)*dr
    r_inner = np.arange(N) * dr
    r_outer = r_inner + dr

    # 球壳体积: V = (4/3)*π*(r_outer³ - r_inner³)
    volumes = (4.0 / 3.0) * np.pi * (r_outer**3 - r_inner**3)

    return dr, volumes


def solid_diffusion_rhs(
    c_s: np.ndarray,
    R_p: float,
    D_s: float,
    flux_surface: float,
    N: int,
) -> np.ndarray:
    """
    计算固相扩散 ODE 系统的右端项 (RHS)。

    求解: d(c_s)/dt = RHS

    RHS 通过对球形扩散方程的有限体积离散化获得:

        ∂c_s/∂t = (1/r²) * ∂/∂r (r² * D_s * ∂c_s/∂r)

    离散化后, 对于内部壳层 i 和 i+1 之间的界面:
        J_{i→i+1} = -D_s * (c_s[i+1] - c_s[i]) / dr  [mol/(m²·s)]
        A_{face} = 4*π*r_face²  [m²]

    壳层 i 的浓度变化率:
        dc_s[i]/dt += J * A / V_i  (来自内侧界面)
        dc_s[i]/dt -= J * A / V_i  (来自外侧界面)

    表面边界条件 (壳层 N-1 的外侧):
        dc_s[N-1]/dt += flux_surface * A_surface / V_last

    符号约定:
    - flux_surface > 0: 锂进入颗粒 (嵌锂)
    - flux_surface < 0: 锂离开颗粒 (脱锂)
    - flux_surface = -I_rxn / F

    对应 DERIVATION.md §2.3, §3.4

    Parameters
    ----------
    c_s : ndarray of shape (N,)
        各壳层的锂浓度 [mol/m³]。
        c_s[0] = 核心浓度, c_s[N-1] = 表面浓度。
    R_p : float
        颗粒半径 [m]。
    D_s : float
        扩散系数 [m²/s]。假设颗粒内均匀 (不随位置变化)。
        注意: 实际 D_s 可能随 c_s 变化, 此处使用单一值。
    flux_surface : float
        颗粒表面的锂通量 [mol/(m²·s)]。
        正值 = 进入颗粒 (嵌锂), 负值 = 离开颗粒 (脱锂)。
        与 Butler-Volmer 反应电流的关系: flux_surface = -I_rxn / F。
    N : int
        壳层数量。

    Returns
    -------
    dc_dt : ndarray of shape (N,)
        各壳层浓度的时间导数 [mol/(m³·s)]。
        dc_dt[i] = dc_s[i]/dt。

    Notes
    -----
    质量守恒检验:
        Σ_i (dc_dt[i] * V_i) = flux_surface * A_surface
    即: 颗粒总锂量变化率 = 表面通量 × 表面积。
    """
    dr = R_p / N  # 壳层厚度 [m]
    dc_dt = np.zeros(N)

    # ===== 内部界面 (壳层 i 和壳层 i+1 之间) =====
    for i in range(N - 1):
        # 界面位置: r = (i+1) * dr
        r_face = (i + 1) * dr

        # 浓度梯度 (一阶差分近似)
        # dc/dr ≈ (c_s[i+1] - c_s[i]) / dr
        dcdr = (c_s[i + 1] - c_s[i]) / dr

        # Fick 定律: J = -D * dc/dr [mol/(m²·s)]
        # 负号表示扩散方向与浓度梯度方向相反
        flux = -D_s * dcdr

        # 界面面积: A = 4*π*r² (球面)
        area = 4.0 * np.pi * r_face**2

        # 壳层 i 的体积: V_i = (4/3)*π*(r_outer³ - r_inner³)
        r_i_inner = i * dr
        r_i_outer = (i + 1) * dr
        vol_i = (4.0 / 3.0) * np.pi * (r_i_outer**3 - r_i_inner**3)

        # 壳层 i+1 的体积
        r_ip1_inner = (i + 1) * dr
        r_ip1_outer = (i + 2) * dr
        vol_ip1 = (4.0 / 3.0) * np.pi * (r_ip1_outer**3 - r_ip1_inner**3)

        # 有限体积守恒: flux 进入壳层 i → dc_dt[i] 增加
        # flux 方向: 从 i+1 流向 i (当 c_s[i+1] > c_s[i] 时 flux > 0)
        dc_dt[i] += flux * area / vol_i      # 壳层 i 获得通量
        dc_dt[i + 1] -= flux * area / vol_ip1  # 壳层 i+1 失去通量

    # ===== 表面边界条件 (壳层 N-1 的外侧) =====
    # 表面通量 flux_surface 从外部进入最外层壳层
    r_surface = R_p  # 颗粒表面半径
    area_surface = 4.0 * np.pi * r_surface**2  # 表面积

    # 最外层壳层的体积
    vol_last = (4.0 / 3.0) * np.pi * (R_p**3 - (R_p - dr)**3)

    # 表面通量对最外层壳层的贡献
    # flux_surface > 0 (嵌锂) → dc_dt[N-1] 增加
    dc_dt[N - 1] += flux_surface * area_surface / vol_last

    return dc_dt
