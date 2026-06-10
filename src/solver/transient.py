"""
瞬态放电求解器 (Transient Discharge Solver)
============================================

物理背景
--------
本模块实现孔隙网络模型的瞬态放电仿真。它将稳态电位求解器 (SteadyStateSolver)
与时间相关的浓度场演化耦合在一起, 形成完整的时域仿真框架。

核心物理过程:
1. 在每个时间步, 用稳态求解器计算当前浓度下的电位分布 (phi_e, phi_s)
2. 由 Butler-Volmer 方程计算各界面反应速率
3. 用 backward Euler 更新电解质浓度 c_e 和固相浓度 c_s

控制方程 (Governing Equations)
-----------------------------
1. 电解质浓度 (Electrolyte Salt Conservation):
       V_i * (c_e^{n+1} - c_e^n) / dt = L_ce * c_e^{n+1} + V * S_e
   其中 L_ce 是扩散拉普拉斯矩阵, S_e 是反应源项。
   对应 DERIVATION.md §3.2。

2. 固相浓度 (Solid Lithium Conservation):
       V_m * (c_s^{n+1} - c_s^n) / dt = L_cs * c_s^{n+1} + V * S_s
   其中 L_cs 是固相扩散拉普拉斯矩阵 (含浓度依赖的 D_s)。
   对应 DERIVATION.md §3.4。

3. 电位方程在每个时间步内通过稳态求解器求解 (准稳态假设):
   - 电解质电位: div(kappa_eff * grad phi_e) = a_s * i_F  (§2.2)
   - 固相电位: div(sigma_eff * grad phi_s) = -a_s * i_F  (§2.4)

时间离散 (Temporal Discretization)
---------------------------------
采用 backward Euler (后向欧拉) 格式:

    V * (c^{n+1} - c^n) / dt = L * c^{n+1} + source

整理为线性系统:
    (V/dt * I - L) * c^{n+1} = (V/dt) * c^n + V * source

其中:
    M = V/dt * I - L  是刚度矩阵 (stiffness matrix)
    rhs = (V/dt) * c^n + V * source  是右端向量

Backward Euler 是一阶隐式格式, 无条件稳定 (对线性扩散方程)。
参见 DERIVATION.md §5.1 和 §6.6。

自适应时间步进 (Adaptive Time Stepping)
--------------------------------------
采用基于以下准则的自适应策略:
1. Newton 迭代次数: 收敛快 → 增大 dt, 收敛慢 → 减小 dt
2. 电压跳变: 过大的电压变化 → 减小 dt (截断电压附近尤其重要)
3. 收敛失败 → 回滚状态, 缩小 dt 重试

特征时间尺度:
    tau_e ~ L^2 / D_e,eff    (电解质扩散时间)
    tau_s ~ l^2 / D_s         (固相扩散时间)
    tau_rxn = 3600 / C_rate   (反应时间, 对应完全放电时长)

参见 DERIVATION.md §6.6 "Time Step Selection"。

C-rate 转换
-----------
放电 C-rate 定义: 1C = 1 小时内完全放电。
    I = C_rate * Q_remaining / 3600  [A]
    I_app = -I / A_cell  [A/m²]  (阳极约定: 放电电流为负)

参见 DERIVATION.md §1.4 "Current and Voltage Sign Convention"。

对应 DERIVATION.md 章节: §1.4, §3.2, §3.4, §5.1, §6.6, §6.7

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op
from scipy import sparse
from scipy.sparse.linalg import spsolve

from src.physics.electrolyte import electrolyte_diffusion_coefficient
from src.physics.ocv import nmc532_ocv
from src.physics.reaction import butler_volmer, exchange_current_density, F, R
from src.physics.separator import SeparatorParams, separator_boundary
from src.physics.solid import nmc532_diffusion_coefficient
from src.solver.steady import SteadyStateSolver


class TransientSolver:
    """
    瞬态放电求解器 — 耦合电位场和浓度场的时间推进。

    在每个时间步内:
    1. 用 SteadyStateSolver 求解当前浓度下的电位分布 (phi_e, phi_s)
    2. 由 Butler-Volmer 方程计算各 e/NMC 界面的反应电流 i_r
    3. 将反应源项加入浓度方程, 用 backward Euler 推进 c_e, c_s

    属性 (Attributes)
    ----------------
    c_e : ndarray, shape (Np,)
        全局电解质 Li+ 浓度 [mol/m³]。初始值 1200 mol/m³。
    c_s : ndarray, shape (Np,)
        全局固相锂浓度 [mol/m³]。NMC 节点初始值 24450 mol/m³ (x≈0.5)。
    L_ce : sparse matrix, shape (n_e, n_e)
        电解质扩散拉普拉斯矩阵 [m³/s]。K_ik = D_e,eff * A / L。
    L_cs : sparse matrix, shape (n_nmc, n_nmc)
        固相扩散拉普拉斯矩阵 [m³/s]。含浓度依赖的 D_s(c_s)。
    """

    def __init__(
        self,
        net: op.network.Cubic,
        T: float = 303.0,
        k0: float = 1e-10,
        geometric_area: float | None = None,
        separator: SeparatorParams | None = None,
    ):
        """
        初始化瞬态求解器。

        Parameters
        ----------
        net : openpnm.network.Cubic
            孔隙网络对象, 包含几何 (pore.volume, throat.area 等)
            和相标识 (pore.electrolyte, pore.nmc, pore.cbd) 属性。
        T : float
            温度 [K], 默认 298.15 K (25°C)。
        k0 : float
            Butler-Volmer 反应速率常数 [m^(2.5) / (mol^0.5 · s)]。
            默认 5e-10。
        geometric_area : float or None
            Projected geometric electrode area [m²]. When provided, it is used
            for C-rate current-density conversion and current boundary scaling
            instead of the network-derived collector throat area.
        """
        self.net = net
        self.T = T
        self.k0 = k0
        self.geometric_area = geometric_area
        self.separator = separator or SeparatorParams(enabled=False)
        self.Np = net.Np        # 总孔隙数
        self.Nt = net.Nt        # 总喉道数
        self.conns = net["throat.conns"]    # 喉道连接关系, shape (Nt, 2)
        self.coords = net["pore.coords"]    # 孔隙坐标, shape (Np, 3)

        # ===== 全局浓度场 =====
        # c_e: 电解质 Li+ 浓度 [mol/m³], 初始 1200 (约 1M LiPF6)
        # c_s: 固相锂浓度 [mol/m³], 初始 24450 (x = c_s/cs_max ≈ 0.5)
        self.c_e = np.ones(self.Np) * 1200.0
        self.c_s = np.ones(self.Np) * 24450.0
        self.c_e_reservoir = 1200.0

        # ===== 相标识掩码 (Phase Masks) =====
        # e_mask: 电解质孔隙 (含 separator 电解质)
        # nmc_mask: NMC 活性材料孔隙 (可嵌锂)
        # cbd_mask: 导电剂 (Carbon Black Domain) 孔隙 (仅导电子, 不储锂)
        # solid_mask: 所有固相孔隙 (NMC + CBD)
        self.e_mask = net["pore.electrolyte"]
        self.nmc_mask = net["pore.nmc"]
        self.cbd_mask = net["pore.cbd"]
        self.solid_mask = self.nmc_mask | self.cbd_mask

        # ===== 局部索引和映射 =====
        # 将全局节点索引映射到各自子矩阵的局部索引
        # e_indices: 电解质节点的全局索引
        # s_indices: 固相节点 (NMC+CBD) 的全局索引
        # nmc_indices: NMC 节点的全局索引
        self.e_indices = np.where(self.e_mask)[0]
        self.s_indices = np.where(self.solid_mask)[0]
        self.nmc_indices = np.where(self.nmc_mask)[0]
        self.n_e = len(self.e_indices)      # 电解质节点数
        self.n_s = len(self.s_indices)      # 固相节点总数
        self.n_nmc = len(self.nmc_indices)  # NMC 节点数 (参与反应)
        # 映射字典: 全局索引 → 局部索引
        self.e_map = {g: l for l, g in enumerate(self.e_indices)}
        self.s_map = {g: l for l, g in enumerate(self.s_indices)}
        self.nmc_map = {g: l for l, g in enumerate(self.nmc_indices)}

        # 构建扩散拉普拉斯矩阵和界面数据
        self._build_diffusion_matrices()

        # 稳态求解器 (处理电位方程和 BV 动力学)
        self._steady = SteadyStateSolver(net, T=T, k0=k0, separator=self.separator)
        if geometric_area is not None:
            if geometric_area <= 0:
                raise ValueError("geometric_area must be positive when provided")
            self._steady._cc_area = float(geometric_area)
        # 过滤: 只保留与稳态求解器一致的反应界面
        self._filter_reactive_interfaces()

    def set_concentration(self, c_e: float = 1200.0, c_s: float = 24450.0):
        """
        设置初始/重置浓度场。

        Parameters
        ----------
        c_e : float
            电解质 Li+ 浓度 [mol/m³]。
        c_s : float
            固相锂浓度 [mol/m³] (仅影响 NMC 节点)。
        """
        self.c_e[:] = c_e
        self.c_s[:] = c_s
        self.c_e_reservoir = c_e
        self._steady.set_concentration(c_e, c_s)
        self._build_electrolyte_diffusion_matrix()
        self._build_solid_diffusion_matrix()

    @property
    def collector_area(self) -> float:
        """
        集流体有效面积 [m²]。

        用于将电流密度 [A/m²] 转换为总电流 [A]:
            I_total = I_app * A_cell

        对应 DERIVATION.md §4.2: 积分集流体边界条件。
        """
        return self._steady._cc_area

    def discharge_capacity_coulombs(self, cs_max: float = 48900.0) -> float:
        """
        计算当前状态下的剩余嵌锂容量 [C] (库仑)。

        物理含义:
            Q_remaining = F * sum_m [ (c_s,max - c_s,m) * V_m ]

        其中求和仅对 NMC 节点, V_m 为控制体体积。
        这是理论上还能嵌入的锂对应的电荷量。

        Parameters
        ----------
        cs_max : float
            NMC 最大锂浓度 [mol/m³], 默认 48900 (NMC532)。

        Returns
        -------
        Q : float
            剩余容量 [C] (库仑)。
        """
        # vacancies = cs_max - c_s (空位浓度), 取 max 保证非负
        vacancies = np.maximum(cs_max - self.c_s[self.nmc_indices], 0.0)
        return float(F * np.sum(vacancies * self._vol_s))

    def current_density_for_c_rate(
        self,
        C_rate: float,
        cs_max: float = 48900.0,
        geometric_area: float | None = None,
    ) -> float:
        """
        将正的放电 C-rate 转换为阳极约定的电流密度 [A/m²]。

        转换公式:
            I = C_rate * Q_remaining / 3600     [A]
            I_app = -I / A_cell                  [A/m²]

        符号约定 (阳极约定):
            - I_app < 0: 放电 (阴极还原, 嵌锂)
            - I_app > 0: 充电 (阳极氧化, 脱锂)

        对应 DERIVATION.md §1.4:
            "a reported discharge magnitude I_dis > 0 corresponds to I_app = -I_dis"

        Parameters
        ----------
        C_rate : float
            放电倍率, 必须为正值。1C = 1 小时完全放电。
        cs_max : float
            NMC 最大锂浓度 [mol/m³]。
        geometric_area : float or None
            Optional projected geometric area [m²] to use for this conversion.
            Defaults to the solver's configured collector area.

        Returns
        -------
        I_app : float
            阳极约定电流密度 [A/m²] (放电时为负值)。
        """
        if C_rate < 0:
            raise ValueError("C_rate must be non-negative; use positive values for discharge")
        # Q_remaining [C], 除以 3600 得 [A·h], 再乘 C_rate 得 [A]
        capacity_c = self.discharge_capacity_coulombs(cs_max=cs_max)
        current_a = C_rate * capacity_c / 3600.0
        area = self.collector_area if geometric_area is None else float(geometric_area)
        if area <= 0:
            raise ValueError("current-density area must be positive")
        # 阳极约定: 放电电流为负
        return -current_a / area

    def characteristic_dt(
        self,
        I_app: float = 0.0,
        C_rate: float | None = None,
        safety: float = 0.05,
        dt_min: float = 1e-6,
        dt_max: float = 300.0,
    ) -> float:
        """
        估计初始时间步长 [s], 基于扩散和反应特征时间尺度。

        三个特征时间:
            tau_e = L^2 / D_e,eff          电解质扩散时间
            tau_s = l_particle^2 / D_s     固相扩散时间
            tau_rxn = Q / |I|              反应时间 (完全放电时长)

        最终 dt = safety * min(tau_e, tau_s, tau_rxn)

        safety (安全因子, 默认 0.05) 保证初始时间步足够小,
        避免 Newton 迭代不收敛或浓度过冲。

        对应 DERIVATION.md §6.6 "Time Step Selection"。

        Parameters
        ----------
        I_app : float
            应用电流密度 [A/m²]。
        C_rate : float or None
            放电 C-rate (可选)。
        safety : float
            安全因子, 默认 0.05。
        dt_min : float
            最小时间步 [s], 默认 1e-6。
        dt_max : float
            最大时间步 [s], 默认 300。

        Returns
        -------
        dt : float
            估计的初始时间步 [s], 裁剪到 [dt_min, dt_max]。
        """
        # ===== 电解质扩散时间 tau_e =====
        # L = 阴极厚度 (x 方向最大跨度)
        # D_e,eff = D_e * epsilon^1.5 (Bruggeman 修正, porosity=0.39)
        x = self.coords[:, 0]
        length = max(float(x.max() - x.min()), 1e-12)
        c_e_values = self.c_e[self.e_mask]
        c_e_mean = float(np.mean(c_e_values)) if c_e_values.size else 1200.0
        D_e = electrolyte_diffusion_coefficient(c_e_mean, self.T) * 0.39**1.5  # [m²/s]
        tau_e = length**2 / D_e

        # ===== 固相扩散时间 tau_s =====
        # 使用 NMC 节点间最小间距作为特征长度
        # D_s 取 NMC 节点的中位数 (浓度依赖)
        if self.n_nmc:
            ds_vals = np.array([
                nmc532_diffusion_coefficient(self.c_s[g], self.T)
                for g in self.nmc_indices
            ])
            D_s = max(float(np.nanmedian(ds_vals)), 1e-30)
            center_dist = np.linalg.norm(
                self.coords[self.conns[:, 1]] - self.coords[self.conns[:, 0]],
                axis=1,
            )
            spacing = max(float(np.min(center_dist)), 1e-12)
            tau_s = spacing**2 / D_s
        else:
            tau_s = np.inf

        # ===== 反应时间 tau_rxn =====
        # tau_rxn = Q_remaining / |I_total|
        # 对于 C-rate: tau_rxn = 3600 / C_rate [s]
        if C_rate is not None and C_rate > 0:
            tau_rxn = 3600.0 / C_rate
        elif abs(I_app) > 0:
            capacity_c = self.discharge_capacity_coulombs()
            current_a = abs(I_app) * self.collector_area
            tau_rxn = capacity_c / current_a if current_a > 0 else np.inf
        else:
            tau_rxn = np.inf

        # 取有限正值的最小时间尺度, 乘以安全因子
        times = [tau for tau in (tau_e, tau_s, tau_rxn) if np.isfinite(tau) and tau > 0]
        dt = safety * min(times) if times else 1.0
        return float(np.clip(dt, dt_min, dt_max))

    def _snapshot_state(self) -> tuple[np.ndarray, np.ndarray]:
        """
        保存当前浓度场快照 (用于自适应步进的状态回滚)。

        Returns
        -------
        snapshot : tuple of ndarray
            (c_e_copy, c_s_copy) — 浓度场的深拷贝。
        """
        return self.c_e.copy(), self.c_s.copy()

    def _restore_state(self, snapshot: tuple[np.ndarray, np.ndarray]) -> None:
        """
        从快照恢复浓度场 (用于时间步拒绝后的状态回滚)。

        Parameters
        ----------
        snapshot : tuple of ndarray
            由 _snapshot_state() 生成的快照。
        """
        self.c_e[:] = snapshot[0]
        self.c_s[:] = snapshot[1]
        self._build_electrolyte_diffusion_matrix()
        self._build_solid_diffusion_matrix()

    def _build_electrolyte_diffusion_matrix(self):
        """Build electrolyte diffusion Laplacian using D_e(c_e)."""
        conns = self.conns
        A = self.net["throat.area"]     # 喉道截面积 [m²]
        L = self.net["throat.length"]   # 喉道长度 [m]

        er, ec, ev = [], [], []  # 稀疏矩阵的行、列、值
        self._ee_throats = []    # 记录 (local_i, local_j, conductance)
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            # 只处理两端都是电解质的喉道
            if self.e_mask[g1] and self.e_mask[g2]:
                l1, l2 = self.e_map[g1], self.e_map[g2]
                D_e1 = electrolyte_diffusion_coefficient(self.c_e[g1], self.T)
                D_e2 = electrolyte_diffusion_coefficient(self.c_e[g2], self.T)
                D_e_avg = 2.0 * D_e1 * D_e2 / (D_e1 + D_e2 + 1e-30)
                De_eff = D_e_avg * 0.39**1.5
                # 扩散传导率: K = D_eff * A / L  [m³/s]
                g_val = De_eff * A[t] / L[t]
                # 填充 2×2 块: [l1,l2]=+K, [l2,l1]=+K, [l1,l1]=-K, [l2,l2]=-K
                er += [l1, l2, l1, l2]; ec += [l2, l1, l1, l2]
                ev += [g_val, g_val, -g_val, -g_val]
                self._ee_throats.append((l1, l2, g_val))
        self.L_ce = sparse.csr_matrix((ev, (er, ec)), shape=(self.n_e, self.n_e))

    def _build_solid_diffusion_matrix(self):
        """Build solid diffusion Laplacian using D_s(c_s)."""
        conns = self.conns
        A = self.net["throat.area"]
        L = self.net["throat.length"]

        # ===== 固相扩散拉普拉斯矩阵 (n_nmc × n_nmc) =====
        # 仅 NMC-NMC 喉道 (CBD 不储锂, 不参与 c_s 方程)
        # D_s 依赖于 c_s, 使用调和平均保证通量连续性
        sr, sc, sv = [], [], []
        self._ss_throats = []
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            if self.nmc_mask[g1] and self.nmc_mask[g2]:
                l1, l2 = self.nmc_map[g1], self.nmc_map[g2]
                # 调和平均: D_avg = 2*D1*D2/(D1+D2)
                # 当 D1=D2 时退化为 D1; 当 D1<<D2 时约等于 2*D1
                D_s1 = nmc532_diffusion_coefficient(self.c_s[g1], self.T)
                D_s2 = nmc532_diffusion_coefficient(self.c_s[g2], self.T)
                D_s_avg = 2.0 * D_s1 * D_s2 / (D_s1 + D_s2 + 1e-30)
                g_val = D_s_avg * A[t] / L[t]
                sr += [l1, l2, l1, l2]; sc += [l2, l1, l1, l2]
                sv += [g_val, g_val, -g_val, -g_val]
                self._ss_throats.append((l1, l2, g_val))
        self.L_cs = sparse.csr_matrix((sv, (sr, sc)), shape=(self.n_nmc, self.n_nmc))

    def _build_diffusion_matrices(self):
        """
        构建电解质和固相的扩散拉普拉斯矩阵。

        对每个 e-e (电解质-电解质) 喉道:
            K_ik = D_e,eff(c_e) * A_ik / L_ik     [m³/s]

        拉普拉斯矩阵 L_ce 的构建:
            L[i,k] = +K_ik  (非对角, 正)
            L[i,i] = -sum_k K_ik  (对角, 负)

        这样 L * c 给出扩散通量:
            (L * c)_i = sum_k K_ik * (c_k - c_i)   [mol/s]

        对应 DERIVATION.md §3.1 和 §3.2。

        固相扩散类似, 但使用 NMC-NMC 喉道,
        且 D_s 依赖于浓度, 采用调和平均:

            D_s,avg = 2 * D_s1 * D_s2 / (D_s1 + D_s2)

        这是两点通量近似的标准做法, 保证通量连续性。
        """
        self._build_electrolyte_diffusion_matrix()
        self._build_solid_diffusion_matrix()

        conns = self.conns
        A = self.net["throat.area"]

        # ===== 孔隙/控制体体积 =====
        self._vol_e = self.net["pore.volume"][self.e_indices]    # 电解质控制体 [m³]
        self._vol_s = self.net["pore.volume"][self.nmc_indices]  # NMC 控制体 [m³]

        # ===== e/NMC 反应界面数据 =====
        # 每个界面记录: (throat_idx, e_local, nmc_local, A_intf, e_global, nmc_global)
        # A_intf: 界面面积 [m²], 用于计算反应源项 i_r * A / F
        self._interfaces = []
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            if self.e_mask[g1] and self.nmc_mask[g2]:
                self._interfaces.append((
                    t, self.e_map[g1], self.nmc_map[g2], A[t], g1, g2,
                ))
            elif self.nmc_mask[g1] and self.e_mask[g2]:
                self._interfaces.append((
                    t, self.e_map[g2], self.nmc_map[g1], A[t], g2, g1,
                ))

    def _filter_reactive_interfaces(self):
        """
        过滤反应界面: 只保留稳态求解器识别的活性界面。

        原因: 预处理阶段可能已将某些不连通的 e/NMC 界面标记为非活性。
        参见 DERIVATION.md §6.4 "Disconnected Pores or Solids"。
        """
        reactive_pairs = {
            (e_g, nmc_g) for e_g, nmc_g, _area in self._steady.reactive_interfaces
            if nmc_g in self.nmc_map
        }
        self._interfaces = [
            interface for interface in self._interfaces
            if (interface[4], interface[5]) in reactive_pairs
        ]

    def step(self, dt: float, I_app: float = 0.0, C_rate: float | None = None) -> dict:
        """
        推进一个时间步。

        算法流程:
        1. 若指定 C_rate, 转换为 I_app (阳极约定)
        2. 用稳态求解器计算当前浓度下的电位分布
        3. 对每个 e/NMC 界面计算 BV 反应电流 i_r
        4. 将反应源项分配到各控制体
        5. 组装线性系统, backward Euler 求解 c_e^{n+1}, c_s^{n+1}
        6. 裁剪浓度到物理范围, 更新全局数组

        浓度方程 (Backward Euler):
            (V/dt * I - L) * c^{n+1} = (V/dt) * c^n + V * source

        其中 source 来自反应:
            c_e source: +i_r * A_intf / (F * V_e)  (阳极产生 Li+)
            c_s source: -i_r * A_intf / (F * V_s)  (阳极消耗固相 Li)

        对应 DERIVATION.md §3.2, §3.4, §3.6。

        Parameters
        ----------
        dt : float
            时间步 [s]。
        I_app : float
            应用电流密度 [A/m²] (阳极约定, 放电为负)。
        C_rate : float or None
            可选的放电 C-rate。若提供, 转换为 I_app。

        Returns
        -------
        result : dict
            phi_e: 电解质电位 [V]
            phi_s: 固相电位 [V]
            c_e: 电解质浓度 [mol/m³]
            c_s: 固相浓度 [mol/m³]
            voltage: 电池电压 [V]
            I_rxn: 各界面反应电流 [A/m²]
            eta: 各界面过电位 [V]
            reaction_mol_rate: 总反应摩尔速率 [mol/s]
            I_app: 实际使用的电流密度 [A/m²]
            dt: 使用的时间步 [s]
            iterations: Newton 迭代次数
            converged: 是否收敛
        """
        # C-rate → I_app 转换
        if C_rate is not None:
            I_app = self.current_density_for_c_rate(C_rate)

        # ===== Step 1: 求解电位场 =====
        # 在当前浓度下, 求解电解质和固相的电位方程
        self._steady.set_concentration(self.c_e, self.c_s)
        pot = self._steady.solve(I_app=I_app, tol=1e-8, max_iter=200)

        # 提取局部电位值 (用于计算 BV 过电位)
        phi_e_local = np.array([pot["phi_e"][g] for g in self.e_indices])
        phi_s_local = np.array([pot["phi_s"][g] for g in self.nmc_indices])

        # ===== Step 2: 计算反应速率和浓度源项 =====
        # 对每个 e/NMC 界面:
        #   eta = phi_s - phi_e - U_eq(c_s)     过电位 [V]
        #   i_r = BV(i0, eta)                    反应电流密度 [A/m²]
        #   flux = i_r * A_intf / F              反应摩尔速率 [mol/s]
        #
        # 符号约定 (阳极约定, DERIVATION.md §1.3):
        #   i_r > 0: 脱锂 (Li 从固相进入电解质)
        #   i_r < 0: 嵌锂/放电 (Li 从电解质进入固相)
        dce_dt = np.zeros(self.n_e)     # 电解质浓度变化率 [mol/(m³·s)]
        dcs_dt = np.zeros(self.n_nmc)   # 固相浓度变化率 [mol/(m³·s)]
        eta_all = np.zeros(self.Nt)     # 各喉道的过电位
        reaction_mol_rate = 0.0         # 总反应摩尔速率 [mol/s]

        for t, le, ls, A_intf, e_g, nmc_g in self._interfaces:
            cs = self.c_s[nmc_g]
            ce = self.c_e[e_g]
            # U_eq: NMC 平衡电位 vs Li/Li+ [V]
            U_eq = nmc532_ocv(cs / 48900.0)
            # i0: 交换电流密度 [A/m²]
            i0 = exchange_current_density(self.k0, ce, cs, 48900.0)

            # 获取电位值 (若求解器失败则用安全默认值)
            pe = phi_e_local[le] if le < len(phi_e_local) else 0.0
            ps = phi_s_local[ls] if ls < len(phi_s_local) else U_eq
            if np.isnan(pe):
                pe = 0.0
            if np.isnan(ps):
                ps = U_eq

            # 过电位: eta = phi_s - phi_e - U_eq
            # 参见 DERIVATION.md §2.5
            eta = ps - pe - U_eq
            # Butler-Volmer 反应电流密度 [A/m²]
            I_rxn = butler_volmer(i0, eta, self.T)
            eta_all[t] = eta

            # ===== 浓度源项 =====
            # flux = I_rxn * A_intf / F  [mol/s]
            # 阳极约定: I_rxn > 0 → 电解质 Li 增加, 固相 Li 减少
            #           I_rxn < 0 → 电解质 Li 减少, 固相 Li 增加
            flux = I_rxn * A_intf / F  # mol/s
            reaction_mol_rate += flux
            # 电解质: dc_e/dt += flux / V_e  (Li+ 由反应产生/消耗)
            dce_dt[le] += flux / self._vol_e[le] if self._vol_e[le] > 0 else 0.0
            # 固相: dc_s/dt -= flux / V_s  (固相 Li 由反应消耗/产生)
            dcs_dt[ls] -= flux / self._vol_s[ls] if self._vol_s[ls] > 0 else 0.0

        # ===== Step 3: 提取局部浓度 =====
        c_e_local = self.c_e[self.e_indices]
        c_s_local = self.c_s[self.nmc_indices]

        # ===== Step 4: Backward Euler 时间推进 =====
        # 离散方程: V*(c_new - c_old)/dt = L*c_new + V*source
        # 整理: (V/dt*I - L)*c_new = (V/dt)*c_old + V*source
        #
        # M = V/dt * I - L  (刚度矩阵)
        # rhs = (V/dt) * c_old + V * source

        # 电解质浓度方程
        M_e = sparse.diags(self._vol_e / dt) - self.L_ce
        rhs_e = (self._vol_e / dt) * c_e_local + self._vol_e * dce_dt
        # 隔膜/Li 侧连接电解液储液库: c_e = c_e_reservoir。
        # 参见 DERIVATION.md §4.1，固定浓度边界补充放电消耗的 Li+。
        if self.n_e and len(self._steady.sep_e) > 0:
            if self.separator.enabled:
                c_bc = separator_boundary(I_app, self.T, self.separator).c_cathode
            else:
                c_bc = self.c_e_reservoir
            M_e = M_e.tolil()
            for i in self._steady.sep_e:
                M_e[i, :] = 0.0
                M_e[i, i] = 1.0
                rhs_e[i] = c_bc
            M_e = M_e.tocsr()
        # 固相浓度方程
        M_s = sparse.diags(self._vol_s / dt) - self.L_cs
        rhs_s = (self._vol_s / dt) * c_s_local + self._vol_s * dcs_dt

        # 求解线性系统
        c_e_new = spsolve(M_e.tocsr(), rhs_e) if self.n_e else c_e_local
        c_s_new = spsolve(M_s.tocsr(), rhs_s) if self.n_nmc else c_s_local

        # ===== Step 5: 裁剪到物理范围 =====
        # c_e > 0: 电解质浓度必须为正
        # c_s < cs_max: 固相浓度不能超过最大值
        # 参见 DERIVATION.md §6.2 "Concentration Bounds"
        c_e_new = np.clip(c_e_new, 1.0, 10000.0)
        c_s_new = np.clip(c_s_new, 1.0, 48899.0)

        # ===== Step 6: 更新全局浓度数组 =====
        self.c_e[self.e_indices] = c_e_new
        self.c_s[self.nmc_indices] = c_s_new
        self._build_electrolyte_diffusion_matrix()
        self._build_solid_diffusion_matrix()

        # NOTE: Voltage is from start-of-step potentials (before concentration
        # update).  V-Q alignment is handled by recording capacity at step START
        # in the caller, matching this voltage to the same time instant.
        return {
            "phi_e": pot["phi_e"],
            "phi_s": pot["phi_s"],
            "c_e": self.c_e.copy(),
            "c_s": self.c_s.copy(),
            "voltage": pot["voltage"],
            "I_rxn": pot["I_rxn"],
            "eta": eta_all,
            "reaction_mol_rate": reaction_mol_rate,
            "I_app": I_app,
            "dt": dt,
            "iterations": pot["iterations"],
            "converged": pot["converged"],
        }

    def adaptive_step(
        self,
        dt: float,
        I_app: float,
        previous_voltage: float | None = None,
        voltage_jump_limit: float = 0.05,
        growth_factor: float = 1.25,
        shrink_factor: float = 0.5,
        fast_iterations: int = 50,
        dt_min: float = 1e-6,
        dt_max: float = 300.0,
        max_retries: int = 24,
    ) -> tuple[dict, float, float]:
        """
        自适应时间步推进 — 带拒绝和增长策略。

        算法:
        1. 尝试 step(trial_dt)
        2. 若收敛失败、电压 NaN、或电压跳变过大 → 拒绝, 回滚, 缩小 dt
        3. 若收敛且 Newton 迭代少 → 增大 dt (乘以 growth_factor)
        4. 若电压跳变接近限制 → 缩小 dt (乘以 shrink_factor)

        拒绝准则:
            - not converged
            - voltage is NaN/Inf
            - |V_new - V_old| > voltage_jump_limit

        对应 DERIVATION.md §6.6 "Time Step Selection"。

        Parameters
        ----------
        dt : float
            建议的时间步 [s]。
        I_app : float
            应用电流密度 [A/m²]。
        previous_voltage : float or None
            上一步的电压 [V], 用于检测电压跳变。
        voltage_jump_limit : float
            允许的最大电压跳变 [V], 默认 0.05。
        growth_factor : float
            dt 增长因子, 默认 1.25 (每次增长 25%)。
        shrink_factor : float
            dt 收缩因子, 默认 0.5 (每次减半)。
        fast_iterations : int
            Newton 迭代低于此值时增大 dt, 默认 50。
        dt_min : float
            最小时间步 [s]。
        dt_max : float
            最大时间步 [s]。
        max_retries : int
            最大重试次数, 默认 24。

        Returns
        -------
        result : dict
            时间步结果。
        dt_used : float
            实际使用的时间步 [s]。
        dt_next : float
            建议的下一步时间步 [s]。
        """
        trial_dt = float(np.clip(dt, dt_min, dt_max))
        last_result = None

        for _ in range(max_retries):
            # 保存状态快照 (用于回滚)
            snapshot = self._snapshot_state()
            result = self.step(trial_dt, I_app=I_app)
            last_result = result

            # ===== 拒绝判断 =====
            voltage = result["voltage"]
            jump = 0.0 if previous_voltage is None else abs(voltage - previous_voltage)
            reject = (
                (not result["converged"])       # Newton 未收敛
                or (not np.isfinite(voltage))    # 电压无效
                or jump > voltage_jump_limit     # 电压跳变过大
            )

            if reject and trial_dt > dt_min:
                # 回滚到快照, 缩小时间步重试
                self._restore_state(snapshot)
                trial_dt = max(trial_dt * shrink_factor, dt_min)
                continue

            # ===== 步长增长/收缩策略 =====
            dt_next = trial_dt
            if result["converged"] and result["iterations"] <= fast_iterations:
                # 收敛快, 增大下一步时间步
                dt_next = min(trial_dt * growth_factor, dt_max)
            elif jump > 0.5 * voltage_jump_limit:
                # 电压跳变接近限制, 缩小下一步时间步
                dt_next = max(trial_dt * shrink_factor, dt_min)

            result["dt"] = trial_dt
            result["dt_next"] = dt_next
            result["voltage_jump"] = jump
            return result, trial_dt, dt_next

        # 所有重试均失败
        voltage = None if last_result is None else last_result.get("voltage")
        raise RuntimeError(
            f"adaptive step failed to converge after {max_retries} retries; "
            f"last voltage={voltage}"
        )

    def run_discharge(
        self,
        C_rate: float,
        cutoff_voltage: float = 2.5,
        dt: float | None = None,
        dt_min: float = 1e-6,
        dt_max: float = 300.0,
        max_steps: int = 200,
        voltage_jump_limit: float = 0.05,
    ) -> dict:
        """
        运行完整的恒流放电仿真, 带自适应时间步进和截止电压检测。

        放电过程:
        1. 从开路状态 (I=0) 开始, 记录初始 OCV
        2. 逐步推进, 记录 V(t) 和 Q(t)
        3. 当电压降至 cutoff_voltage 时停止
        4. 使用线性插值确定精确的截止时刻

        截止电压插值:
            当 V_old > cutoff >= V_new 时, 用线性插值估计精确截止时刻:
                frac = (V_old - cutoff) / (V_old - V_new)
                t_cutoff = t_old + frac * dt

        对应 DERIVATION.md §6.6 和 §7.2。

        Parameters
        ----------
        C_rate : float
            放电 C-rate, 必须为正值。1C = 1 小时完全放电。
        cutoff_voltage : float
            截止电压 [V], 默认 2.5 V。
        dt : float or None
            初始时间步 [s]。若 None, 用 characteristic_dt() 估计。
        dt_min : float
            最小时间步 [s]。
        dt_max : float
            最大时间步 [s]。
        max_steps : int
            最大时间步数, 默认 200。
        voltage_jump_limit : float
            允许的最大电压跳变 [V]。

        Returns
        -------
        result : dict
            time: 时间数组 [s]
            capacity: 容量数组 [A·s/m²]
            capacity_Ah_m2: 容量数组 [A·h/m²]
            voltage: 电压数组 [V]
            dt: 各步时间步长 [s]
            iterations: 各步 Newton 迭代次数
            I_app: 应用电流密度 [A/m²]
            C_rate: 放电 C-rate
            cutoff_voltage: 截止电压 [V]
            reached_cutoff: 是否达到截止电压
        """
        if C_rate <= 0:
            raise ValueError("C_rate must be positive for discharge")

        # C-rate → 电流密度 (阳极约定, 放电 I_app < 0)
        I_app = self.current_density_for_c_rate(C_rate)
        # 初始时间步估计
        step_dt = dt if dt is not None else self.characteristic_dt(
            I_app=I_app, C_rate=C_rate, dt_min=dt_min, dt_max=dt_max,
        )

        # ===== 初始状态 (开路) =====
        self._steady.set_concentration(self.c_e, self.c_s)
        initial = self._steady.solve(I_app=0.0)
        previous_voltage = initial["voltage"]
        times = [0.0]
        capacities = [0.0]
        voltages = [previous_voltage]
        dts = [0.0]
        iterations = [initial["iterations"]]

        # 检查初始 OCV 是否已低于截止电压
        done = previous_voltage <= cutoff_voltage
        for _ in range(max_steps):
            if done:
                break

            # 自适应时间步推进
            result, dt_used, step_dt = self.adaptive_step(
                step_dt,
                I_app=I_app,
                previous_voltage=previous_voltage,
                voltage_jump_limit=voltage_jump_limit,
                dt_min=dt_min,
                dt_max=dt_max,
            )

            t_new = times[-1] + dt_used
            # 累积容量: Q = |I_app| * t [A·s/m²]
            q_new = abs(I_app) * t_new
            v_new = result["voltage"]

            # ===== 截止电压检测与线性插值 =====
            if previous_voltage > cutoff_voltage >= v_new:
                # 线性插值: 精确到截止时刻
                frac = (previous_voltage - cutoff_voltage) / max(previous_voltage - v_new, 1e-30)
                times.append(times[-1] + frac * dt_used)
                capacities.append(capacities[-1] + frac * (q_new - capacities[-1]))
                voltages.append(cutoff_voltage)
                dts.append(frac * dt_used)
                iterations.append(result["iterations"])
                done = True
                break

            # 记录数据点
            times.append(t_new)
            capacities.append(q_new)
            voltages.append(v_new)
            dts.append(dt_used)
            iterations.append(result["iterations"])
            previous_voltage = v_new

        return {
            "time": np.asarray(times),
            "capacity": np.asarray(capacities),
            "capacity_Ah_m2": np.asarray(capacities) / 3600.0,  # [A·s/m²] → [A·h/m²]
            "voltage": np.asarray(voltages),
            "dt": np.asarray(dts),
            "iterations": np.asarray(iterations),
            "I_app": I_app,
            "C_rate": C_rate,
            "cutoff_voltage": cutoff_voltage,
            "reached_cutoff": done,
        }
