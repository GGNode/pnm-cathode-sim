"""
单孔放电仿真器 (Single-Pore Discharge Simulator)
=================================================

物理背景
--------
本模块实现单个 NMC532 球形颗粒在充分搅拌电解质中的恒流放电仿真。
这是 Phase 1 验证测试用例, 用于隔离验证以下物理子系统:
    1. Butler-Volmer 界面反应动力学
    2. 固相锂扩散 (球坐标 Fick 扩散)
    3. NMC532 OCV 曲线

简化的物理假设 (相比完整 PNM):
    - 电解质浓度恒定: c_e = const (无传输限制)
    - 固相电位均匀: phi_s = const (无电子电阻)
    - 仅考虑固相扩散 + BV 表面反应
    - 颗粒为单一孔隙节点, 不涉及网络拓扑

这等价于 Newman 型单粒子模型 (Single Particle Model, SPM) 的核心部分。

对应 DERIVATION.md 章节: §7.1 (0C 平衡), §7.2 (小电流极限), §7.5 (固相扩散限制)

球形颗粒扩散方程
-----------------
在球坐标系中, 固相锂扩散方程为:

    ∂c_s/∂t = (1/r²) ∂/∂r [ r² D_s ∂c_s/∂r ]

边界条件:
    - r = 0 (球心): ∂c_s/∂r = 0          (对称性)
    - r = R_p (表面): -D_s ∂c_s/∂r = J_in  (表面通量)

其中 J_in = -I_rxn / F 为进入颗粒的锂摩尔通量 [mol/(m²·s)]。
    I_rxn < 0 (放电/嵌锂) → J_in > 0 (锂进入颗粒)

离散化采用壳层模型 (Shell Model):
    将球体分为 N_shell 个同心球壳, 每个壳层的浓度均匀。
    使用有限体积法 (FVM) 离散扩散方程。
    参见 src/physics/solid.py 中的 discretize_spherical_particle()。

恒流放电算法
-----------
对于恒流放电 (galvanostatic discharge):
    1. 给定 C-rate, 计算施加电流 I_app (阳极约定, 负值)
    2. 在每个时间步:
       a. 由表面浓度 c_s,surf 计算 U_eq = OCV(c_s,surf / cs_max)
       b. 由 BV 方程求解过电位 eta, 使 BV(i0, eta) = I_app
       c. 电池电压 V = phi_s - phi_e = eta + U_eq  (因 phi_e = 0)
       d. 计算表面通量 J_in = -I_rxn / F
       e. 前向欧拉 (Forward Euler) 推进固相浓度

注意: 这里使用半隐式格式 — BV 反应隐式求解 (通过 bisection),
固相扩散显式推进 (Forward Euler)。

C-rate 定义与电流密度转换
-------------------------
1C = 1 小时完全放电:
    I_1C = F * cs_max * R_p / 3 / 3600   [A/m²]

这是一个近似: 假设所有固相锂在 1 小时内均匀提取。
对于 C-rate 放电:
    I_app = -C_rate * I_1C   [A/m²]  (阳极约定, 放电为负)

符号约定
-------
阳极约定 (Anodic Convention):
    - I_rxn > 0: 脱锂 (阳极, 氧化)
    - I_rxn < 0: 嵌锂 (阴极, 还原, 放电)
    - eta > 0: 过电位驱动脱锂
    - eta < 0: 过电位驱动嵌锂

电池电压:
    V_cell = phi_s - phi_e(separator)
    对于单孔: phi_e = 0 (参考电位), phi_s = eta + U_eq

对应 DERIVATION.md 章节: §1.3, §1.4, §2.5, §7.1-7.5

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np

from pnmcathode.physics.ocv import nmc532_ocv
from pnmcathode.physics.reaction import butler_volmer, exchange_current_density, F, R
from pnmcathode.physics.solid import (
    nmc532_diffusion_coefficient,
    discretize_spherical_particle,
    solid_diffusion_rhs,
)


class SinglePoreDischarge:
    """
    单个 NMC532 球形颗粒的恒流放电仿真器。

    物理模型:
    - 充分搅拌的电解质 (c_e = const, 无浓度梯度)
    - 无电子电阻 (phi_s = uniform)
    - Butler-Volmer 表面反应动力学
    - 球坐标固相 Fick 扩散

    使用方法:
        sim = SinglePoreDischarge(R_p=5e-6, C_rate=1.0)
        result = sim.run(cutoff_V=3.0)
        # result['voltage'], result['capacity'] 等
    """

    def __init__(
        self,
        R_p: float = 5e-6,
        N_shell: int = 10,
        c_s_init: float = 0.5,
        c_e: float = 1200.0,
        T: float = 303.0,
        C_rate: float = 1.0,
        k0: float = 1e-10,
        cs_max: float = 48900.0,
    ):
        """
        初始化单孔放电仿真器。

        Parameters
        ----------
        R_p : float
            颗粒半径 [m], 默认 5 µm (典型 NMC532 二次颗粒)。
        N_shell : int
            扩散壳层数, 默认 10。越多越精确, 但计算量增大。
        c_s_init : float
            初始锂化度 x = c_s / c_s_max [0,1]。
            0 = 完全脱锂, 1 = 完全嵌锂。
            默认 0.5 (半充电状态)。
        c_e : float
            电解质 Li+ 浓度 [mol/m³], 默认 1200 (约 1M LiPF6)。
        T : float
            温度 [K], 默认 298.15 K (25°C)。
        C_rate : float
            放电倍率, 默认 1.0C。
            1C = 1 小时完全放电, 0.5C = 2 小时, 2C = 30 分钟。
        k0 : float
            BV 反应速率常数 [m^(2.5) / (mol^0.5 · s)]。
        cs_max : float
            NMC532 最大锂浓度 [mol/m³], 默认 48900。
        """
        self.R_p = R_p
        self.N_shell = N_shell
        self.c_e = c_e
        self.T = T
        self.C_rate = C_rate
        self.k0 = k0
        self.cs_max = cs_max

        # ===== 颗粒离散化 =====
        # discretize_spherical_particle 返回:
        #   dr: 壳层厚度 [m]
        #   volumes: 各壳层体积 [m³] (按球壳体积公式计算)
        self.dr, self.volumes = discretize_spherical_particle(R_p, N_shell)
        self.total_volume = sum(self.volumes)

        # ===== 初始条件 =====
        # c_s_init 是锂化度 (0~1), 乘以 cs_max 得浓度 [mol/m³]
        self.c_s_init = c_s_init * cs_max
        self.surface_area = 4.0 * np.pi * R_p**2  # 颗粒表面积 [m²]

        # ===== 1C 电流密度的估算 =====
        # 理论: 对于均匀提取, 1 小时内提取所有锂:
        #   Q = F * cs_max * V_particle   [C]
        #   I = Q / (3600 * A_surface)   [A/m²]
        #   I_1C = F * cs_max * R_p / 3 / 3600   [A/m²]
        #
        # 推导: V_particle = (4/3)*pi*R_p^3, A_surface = 4*pi*R_p^2
        #   V/A = R_p/3
        #   I_1C = F * cs_max * (V/A) / 3600 = F * cs_max * R_p / 3 / 3600
        self.I_1C = F * cs_max * R_p / 3.0 / 3600.0  # A/m² of surface

    def run(
        self,
        cutoff_V: float = 3.0,
        dt: float = 5.0,
        t_max: float = 40000.0,
    ) -> dict:
        """
        运行恒流放电仿真。

        算法流程:
        1. 初始化均匀浓度 c_s = c_s_init
        2. 计算施加电流 I_app = -C_rate * I_1C (阳极约定)
        3. 循环直到截止电压或 t_max:
           a. 由表面浓度计算 U_eq = OCV(x_surface)
           b. 由 exchange_current_density() 计算 i0
           c. 用二分法求解 BV 方程: BV(i0, eta) = I_app → eta
           d. 电池电压 V = eta + U_eq (因 phi_e = 0, phi_s = eta + U_eq)
           e. 计算表面通量 flux_surface = -I_rxn / F
           f. 固相扩散 RHS (球坐标离散)
           g. 前向欧拉推进: c_s = c_s + dt * rhs
           h. 裁剪 c_s 到物理范围

        对应 DERIVATION.md §7.5 (固相扩散限制)。

        Parameters
        ----------
        cutoff_V : float
            截止电压 [V], 默认 3.0 V。
            NMC532 的有效放电范围约 2.5-4.2 V vs Li/Li+。
        dt : float
            时间步 [s], 默认 5.0。
        t_max : float
            最大仿真时间 [s], 默认 40000 (约 11 小时)。

        Returns
        -------
        result : dict
            time: 时间数组 [s]
            voltage: 电压数组 [V]
            capacity: 容量数组 [A·h/m²]
            c_s_surface: 表面浓度数组 [mol/m³]
            c_s_bulk: 体平均浓度数组 [mol/m³]
        """
        N = self.N_shell
        # 初始化固相浓度 (均匀)
        c_s = np.ones(N) * self.c_s_init

        # 施加电流密度 (阳极约定: 放电为负)
        I_app = -self.C_rate * self.I_1C  # A/m² (negative = cathodic = lithiation)

        times = [0.0]
        voltages = []
        capacities = []
        c_s_surfaces = []
        c_s_bulks = []

        t = 0.0
        Q_delivered = 0.0  # 累积交付容量 [A·h/m²]

        # ===== 初始状态 (t=0, 开路) =====
        cs_surface = c_s[-1]                # 最外层壳层 = 表面
        x_surface = cs_surface / self.cs_max  # 锂化度 x = c_s / cs_max
        U_eq = nmc532_ocv(x_surface)        # 平衡电位 vs Li/Li+ [V]
        # 开路时 V = U_eq (无过电位)
        voltages.append(U_eq)
        c_s_surfaces.append(cs_surface)
        c_s_bulks.append(np.mean(c_s))
        capacities.append(0.0)

        # ===== 主时间循环 =====
        while t < t_max:
            # --- Step 1: 求解表面条件 ---
            # 恒流约束: BV(i0, eta) = I_app
            # 求解 eta (过电位), 然后计算电压

            cs_surface = c_s[-1]                    # 表面浓度 [mol/m³]
            x_surface = cs_surface / self.cs_max    # 锂化度
            U_eq = nmc532_ocv(x_surface)            # 平衡电位 [V]

            # 交换电流密度 [A/m²]
            i0 = exchange_current_density(
                k0=self.k0, ce=self.c_e, cs=cs_surface, cs_max=self.cs_max,
            )

            # ===== 求解 BV 方程中的 eta =====
            # BV: i0 * [exp(a_a*F*eta/(RT)) - exp(-a_c*F*eta/(RT))] = I_app
            #
            # 对于对称 alpha=0.5, 可解析求解:
            #   i0 * 2*sinh(0.5*F*eta/(RT)) = I_app
            #   eta = (2*RT/F) * arcsinh(I_app / (2*i0))
            #
            # 但为了一般性 (任意 alpha), 使用二分法求解。
            eta = self._solve_eta_bisection(i0, I_app)

            # ===== 电池电压 =====
            # V_cell = phi_s - phi_e(separator)
            # 对于单孔: phi_e = 0 (参考), phi_s = eta + U_eq
            # 因此 V_cell = eta + U_eq
            #
            # 对应 DERIVATION.md §4.5: V_cell = phi_s,collector - phi_s,Li
            V_cell = eta + U_eq  # phi_s

            # 截止电压检查
            if V_cell < cutoff_V:
                break

            # --- Step 2: 计算表面通量 ---
            # BV 反应电流 [A/m²]
            I_rxn = butler_volmer(i0=i0, eta=eta, T=self.T)
            # 锂进入颗粒的摩尔通量 [mol/(m²·s)]
            # 放电: I_rxn < 0 (还原) → flux_surface > 0 (Li 进入颗粒)
            #   NMC + Li+ + e- → Li-NMC
            # flux_surface = -I_rxn / F
            flux_surface = -I_rxn / F

            # --- Step 3: 更新固相浓度 ---
            # 半隐式格式: RHS 用当前 c_s, 时间推进用 Forward Euler
            # solid_diffusion_rhs() 计算球坐标离散的扩散 RHS
            #   rhs_i = (1/r²) ∂/∂r [r² D_s ∂c_s/∂r]
            D_s = nmc532_diffusion_coefficient(np.mean(c_s), self.T)
            rhs = solid_diffusion_rhs(c_s, self.R_p, D_s, flux_surface, N)
            c_s = c_s + dt * rhs

            # 裁剪到物理范围: c_s ∈ (0, cs_max)
            # 参见 DERIVATION.md §6.2 "Concentration Bounds"
            c_s = np.clip(c_s, 1.0, self.cs_max - 1.0)

            # --- Step 4: 累积容量 ---
            # Q += |I_rxn| * A_surface * dt / 3600  [A·h/m²]
            Q_delivered += abs(I_rxn) * self.surface_area * dt / 3600.0  # Ah

            t += dt
            times.append(t)
            voltages.append(V_cell)
            c_s_surfaces.append(c_s[-1])
            c_s_bulks.append(np.mean(c_s))
            capacities.append(Q_delivered)

        return {
            "time": np.array(times),
            "voltage": np.array(voltages),
            "capacity": np.array(capacities),
            "c_s_surface": np.array(c_s_surfaces),
            "c_s_bulk": np.array(c_s_bulks),
        }

    def _solve_eta_bisection(
        self, i0: float, I_target: float, tol: float = 1e-8,
    ) -> float:
        """
        用二分法求解 BV 方程中的过电位 eta。

        求解: butler_volmer(i0, eta) = I_target
        即: i0 * [exp(a_a*F*eta/(RT)) - exp(-a_c*F*eta/(RT))] = I_target

        BV 方程关于 eta 单调递增:
            eta > 0 → I_rxn > 0 (阳极)
            eta < 0 → I_rxn < 0 (阴极)

        因此对于给定的 I_target, 存在唯一的 eta。

        二分法原理:
            1. 选择初始区间 [eta_low, eta_high] 使得 f(eta_low)*f(eta_high) < 0
            2. 取中点 eta_mid, 检查 f(eta_mid) 的符号
            3. 缩小区间直到 |f(eta_mid)| < tol

        Parameters
        ----------
        i0 : float
            交换电流密度 [A/m²]。
        I_target : float
            目标反应电流 [A/m²] (等于 I_app)。
        tol : float
            收敛容差 [A/m²], 默认 1e-8。

        Returns
        -------
        eta : float
            过电位 [V]。
        """
        # 初始搜索区间
        # 阴极 (I_target < 0): eta < 0, 搜索 [-0.5, 0]
        # 阳极 (I_target > 0): eta > 0, 搜索 [0, 0.5]
        eta_low = -0.5 if I_target < 0 else -0.01
        eta_high = 0.01 if I_target < 0 else 0.5

        # 扩展区间直到包含根 (f_low * f_high < 0)
        for _ in range(20):
            f_low = butler_volmer(i0, eta_low, self.T) - I_target
            f_high = butler_volmer(i0, eta_high, self.T) - I_target
            if f_low * f_high < 0:
                break
            eta_low *= 2
            eta_high *= 2

        # 二分迭代
        for _ in range(100):
            eta_mid = 0.5 * (eta_low + eta_high)
            f_mid = butler_volmer(i0, eta_mid, self.T) - I_target
            if abs(f_mid) < tol:
                return eta_mid
            if f_mid < 0:
                eta_low = eta_mid
            else:
                eta_high = eta_mid

        return 0.5 * (eta_low + eta_high)
