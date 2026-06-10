"""
稳态电位求解器模块 (Steady-State Potential Solver)
====================================================

物理背景
--------
本模块在孔网络上求解耦合的电解质电位 (phi_e) 和固相电位 (phi_s) 分布。

核心方程 (对应 DERIVATION.md §3.3, §3.5):

1. 电解质电荷守恒 (Ohm's law):
    F_{phi_e,i} = Σ_k G_ik^e (phi_k^e - phi_i^e) + Σ_r i_r A_r = 0

2. 固相电荷守恒 (Ohm's law):
    F_{phi_s,m} = Σ_n G_mn^s (phi_n^s - phi_m^s) - Σ_r i_r A_r + F_{BC,m}^s = 0

其中:
- G_ik^e = kappa * A_ik / L_ik: 电解质离子电导 [S]
- G_mn^s = sigma * A_mn / L_mn: 固相电子电导 [S]
- i_r: Butler-Volmer 反应电流密度 [A/m²]
- A_r: 反应界面面积 [m²]
- F_{BC,m}^s: 集流体边界电流贡献

线性化求解策略
--------------
由于 Butler-Volmer 方程的非线性, 采用逐次超松弛 (SOR) 迭代:

1. 在当前电位下计算 BV 的线性化电导:
    g_bv = dI/deta * A_r
    其中 dI/deta = i0 * F/(RT) * (alpha_a * exp(arg_a) + alpha_c * exp(arg_c))

2. 将 BV 反应视为电导耦合项:
    电解质方程: (L_e - diag(g_bv_e)) @ phi_e + g_bv @ phi_s = rhs_e
    固相方程:   (L_s - diag(g_bv_s)) @ phi_s + g_bv @ phi_e = rhs_s

3. 求解线性系统, 更新电位

符号约定 (Sign Convention)
--------------------------
对应 DERIVATION.md §1.3, §1.4:
- 反应电流 i_r: 阢极约定, 正值 = 脱锂 (anodic)
- I_app > 0: 阳极电流, 正极充电/脱锂
- I_app < 0: 阴极电流, 放电/嵌锂
- V_cell = phi_s(collector) - phi_e(separator)

边界条件
--------
对应 DERIVATION.md §4.1, §4.2:
- 隔膜端 (x=0): phi_e = 0 V (Dirichlet, 电解质电位参考点)
- 集流体端 (x=L): I_app 施加在固相上 (Neumann)
- 电解质在集流体端: 无通量 (∂phi_e/∂n = 0)
- 固相在隔膜端: 无电子电流 (∂phi_s/∂n = 0)

连通性检查
----------
对应 DERIVATION.md §6.4 "Disconnected Pores or Solids":
- 仅保留与隔膜端连通的电解质分量
- 仅保留与集流体端连通且有活性界面的固相分量
- 未连通的节点被固定 (Dirichlet) 以避免奇异矩阵

对应 DERIVATION.md 章节: §2.2, §2.4, §3.3, §3.5, §4, §6.4

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

from src.physics.ocv import nmc532_ocv
from src.physics.reaction import butler_volmer, exchange_current_density, F, R


class SteadyStateSolver:
    """
    孔网络上的稳态电位求解器。

    求解耦合的电解质/固相电位分布, 使用 Butler-Volmer 反应耦合。
    采用线性化 + SOR 迭代策略处理 BV 非线性。

    对应 DERIVATION.md §3.3 (电解质电位) 和 §3.5 (固相电位)。
    """

    def __init__(self, net: op.network.Cubic, T: float = 298.15, k0: float = 5e-10):
        """
        初始化稳态求解器。

        Parameters
        ----------
        net : openpnm.network.Cubic
            孔网络 (由 generator.create_cathode_network 创建)。
        T : float
            温度 [K]。默认 298.15 K (25°C)。
        k0 : float
            Butler-Volmer 反应速率常数。
        """
        self.net = net
        self.T = T
        self.k0 = k0
        self.Np = net.Np  # 总节点数
        self.Nt = net.Nt  # 总喉道数
        self.conns = net["throat.conns"]  # 连接矩阵 [Nt, 2]
        self.coords = net["pore.coords"]  # 节点坐标 [Np, 3]

        # 初始浓度 (均匀)
        self.c_e = np.ones(self.Np) * 1200.0  # 电解质 Li+ 浓度 [mol/m³]
        self.c_s = np.ones(self.Np) * 24450.0  # 固相 Li 浓度 [mol/m³] (x=0.5)

        # 相掩码 (布尔数组)
        self.e_mask = net["pore.electrolyte"]  # 电解质节点
        self.nmc_mask = net["pore.nmc"]        # NMC 节点
        self.cbd_mask = net["pore.cbd"]        # CBD 节点
        self.solid_mask = self.nmc_mask | self.cbd_mask  # 固相节点 (NMC + CBD)

        # 局部索引映射: 全局索引 ↔ 局部索引
        # e_indices[i] = 全局索引, e_map[全局] = 局部索引
        self.e_indices = np.where(self.e_mask)[0]  # 电解质节点的全局索引
        self.s_indices = np.where(self.solid_mask)[0]  # 固相节点的全局索引
        self.n_e = len(self.e_indices)  # 电解质节点数
        self.n_s = len(self.s_indices)  # 固相节点数
        self.e_map = {g: l for l, g in enumerate(self.e_indices)}  # 全局→局部
        self.s_map = {g: l for l, g in enumerate(self.s_indices)}  # 全局→局部

        # 计算电导率和构建矩阵
        self._compute_conductances()
        self._build_matrices()
        self._phi_e_guess: np.ndarray | None = None
        self._phi_s_guess: np.ndarray | None = None

    def set_concentration(self, c_e: float = 1200.0, c_s: float = 24450.0):
        """设置均匀浓度场。"""
        self.c_e[:] = c_e
        self.c_s[:] = c_s

    def _compute_conductances(self):
        """
        计算喉道电导。

        对应 DERIVATION.md §3.1:
        - 电解质离子电导: G_ik^e = kappa * A_ik / L_ik [S]
        - 固相电子电导: G_mn^s = sigma_mn * A_mn / L_mn [S]

        电导率:
        - kappa(c_e) = F² * D_e * c_e / (RT) [S/m] (Nernst-Einstein)
          参见 DERIVATION.md §2.1
        - sigma_nmc = 0.01 S/m (NMC532 电子导电率)
        - sigma_cbd = 760 S/m (CBD 电子导电率)

        调和平均: 对于 NMC-CBD 喉道, 使用调和平均:
            sigma_avg = 2 * sig1 * sig2 / (sig1 + sig2)
        """
        net = self.net
        A = net["throat.area"]     # 喉道截面积 [m²]
        L = net["throat.length"]   # 喉道长度 [m]

        # 电解质参数
        D_e = 2.0e-10  # 电解质扩散系数 [m²/s]
        # Nernst-Einstein 电导率: kappa = F² * D_e * c_e / (RT)
        # 在 c_e=1200 mol/m³, T=298.15 K 时
        kappa = 1200.0 * F**2 * D_e / (R * self.T)  # [S/m]

        # 固相导电率
        sigma_nmc, sigma_cbd = 0.01, 760.0  # [S/m]

        # ===== 电解质离子电导 G_e =====
        # 仅在两个端点都是电解质节点时才非零
        p1, p2 = self.conns[:, 0], self.conns[:, 1]
        self.G_e = np.where(self.e_mask[p1] & self.e_mask[p2], kappa * A / L, 0.0)

        # ===== 固相电子电导 G_s =====
        # 仅在两个端点都是固相节点时才非零
        s1, s2 = self.solid_mask[p1], self.solid_mask[p2]
        # 根据节点相标记选择导电率
        sig1 = np.where(self.nmc_mask[p1], sigma_nmc, np.where(self.cbd_mask[p1], sigma_cbd, 0.0))
        sig2 = np.where(self.nmc_mask[p2], sigma_nmc, np.where(self.cbd_mask[p2], sigma_cbd, 0.0))
        # 调和平均 (当材料不同时)
        # 对应 DERIVATION.md §3.1: "use harmonic or series averaging when
        # material properties differ across adjacent half-throats"
        sig_avg = 2.0 * sig1 * sig2 / (sig1 + sig2 + 1e-30)
        self.G_s = np.where(s1 & s2, sig_avg * A / L, 0.0)

        # ===== 电解质/NMC 反应界面 =====
        # 对应 DERIVATION.md §1.1: "Electrolyte/active-material interfaces:
        # reactions r in R, connecting one electrolyte node i(r) to one active solid node m(r)"
        self.interfaces = []
        for t in range(self.Nt):
            i, j = int(self.conns[t, 0]), int(self.conns[t, 1])
            if self.e_mask[i] and self.nmc_mask[j]:
                self.interfaces.append((i, j, A[t]))  # (电解质节点, NMC节点, 面积)
            elif self.nmc_mask[i] and self.e_mask[j]:
                self.interfaces.append((j, i, A[t]))

    def _build_matrices(self):
        """
        构建 Laplacian 矩阵和边界条件。

        电解质 Laplacian 矩阵 L_e:
            L_e[i,k] = G_ik^e (i≠k, 连接)
            L_e[i,i] = -Σ_k G_ik^e (对角线)
        对应 DERIVATION.md §3.3: F_{phi_e,i} = Σ_k G_ik (phi_k - phi_i)

        固相 Laplacian 矩阵 L_s:
            类似结构, 使用固相电导 G_s。

        边界条件:
        - 隔膜端 (x=0): 电解质 Dirichlet (phi_e = 0)
        - 集流体端 (x=L): 固相 Neumann (I_app)
        """
        n_e, n_s, conns = self.n_e, self.n_s, self.conns

        # ===== 电解质 Laplacian L_e =====
        # 对应 DERIVATION.md §3.3: 电解质电位残差
        # 矩阵元素: L_e[l1,l2] = G_e, L_e[l1,l1] -= G_e
        er, ec, ev = [], [], []
        for t in range(self.Nt):
            if self.G_e[t] == 0:
                continue
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            l1, l2 = self.e_map[g1], self.e_map[g2]
            # 非对角: +G (流入), 对角: -G (流出)
            er += [l1, l2, l1, l2]; ec += [l2, l1, l1, l2]
            ev += [self.G_e[t], self.G_e[t], -self.G_e[t], -self.G_e[t]]
        self.L_e = sparse.csr_matrix((ev, (er, ec)), shape=(n_e, n_e))

        # ===== 固相 Laplacian L_s =====
        # 对应 DERIVATION.md §3.5: 固相电位残差
        sr, sc, sv = [], [], []
        for t in range(self.Nt):
            if self.G_s[t] == 0:
                continue
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            l1, l2 = self.s_map[g1], self.s_map[g2]
            sr += [l1, l2, l1, l2]; sc += [l2, l1, l1, l2]
            sv += [self.G_s[t], self.G_s[t], -self.G_s[t], -self.G_s[t]]
        self.L_s = sparse.csr_matrix((sv, (sr, sc)), shape=(n_s, n_s))

        # ===== 边界节点集合 =====
        # 对应 DERIVATION.md §4.1 (隔膜端) 和 §4.2 (集流体端)
        x = self.coords[:, 0]  # x 坐标
        x_min, x_max = x.min(), x.max()
        span = x_max - x_min if x_max > x_min else 1.0

        # 隔膜端电解质节点 (x ≈ x_min, 占总跨度的 15%)
        self.sep_e = np.array([self.e_map[g] for g in self.e_indices if (x[g] - x_min) < span * 0.15])
        # 隔膜端固相节点
        self.sep_s = np.array([self.s_map[g] for g in self.s_indices if (x[g] - x_min) < span * 0.15])
        # 集流体端固相节点 (x ≈ x_max)
        self.cc_s = np.array([self.s_map[g] for g in self.s_indices if (x_max - x[g]) < span * 0.15])

        # 确保至少有一个边界节点
        if len(self.sep_s) == 0 and n_s > 0:
            self.sep_s = np.array([int(np.argmin(x[self.s_indices] - x_min))])
        if len(self.cc_s) == 0 and n_s > 0:
            self.cc_s = np.array([int(np.argmax(x[self.s_indices]))])

        # 连通性检查
        self._mark_connected_components()
        self._filter_reactive_interfaces()

        # 集流体面积 (用于将 I_app [A/m²] 转换为 [A])
        cc_g = set(self.s_indices[self.active_cc_s])
        self._cc_area = 0.0
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            if g1 in cc_g or g2 in cc_g:
                self._cc_area += self.net["throat.area"][t]
        if self._cc_area == 0:
            self._cc_area = 1.0

    def _mark_connected_components(self):
        """
        标记可支持反应的离子/电子分量。

        对应 DERIVATION.md §6.4 "Disconnected Pores or Solids":
        - 仅保留与隔膜端连通的电解质分量
        - 仅保留与集流体端连通且有活性界面的固相分量

        未连通的节点不能参与电化学反应 (无离子/电子通路)。
        """
        # 电解质连通性: 哪些电解质节点与隔膜端 Dirichlet BC 连通?
        if self.n_e:
            n_comp_e, labels_e = connected_components(self.L_e != 0, directed=False)
            sep_components = set(labels_e[self.sep_e]) if len(self.sep_e) else set()
            # active_e[i] = True 表示节点 i 与隔膜端连通
            self.active_e = np.array([label in sep_components for label in labels_e])
        else:
            self.active_e = np.array([], dtype=bool)

        # 固相连通性: 哪些固相节点与集流体端连通且有活性界面?
        if self.n_s:
            n_comp_s, labels_s = connected_components(self.L_s != 0, directed=False)
            cc_components = set(labels_s[self.cc_s]) if len(self.cc_s) else set()

            # 检查哪些固相分量有活性的电解质/NMC 界面
            interface_components = set()
            for e_g, s_g, _area in self.interfaces:
                e_l = self.e_map[e_g]
                s_l = self.s_map[s_g]
                if self.active_e[e_l] and labels_s[s_l] in cc_components:
                    interface_components.add(labels_s[s_l])

            # 活性固相分量 = 与集流体连通 且 有活性界面
            active_components = cc_components & interface_components
            self.active_s = np.array([label in active_components for label in labels_s])
            # 活性集流体节点
            self.active_cc_s = np.array([i for i in self.cc_s if self.active_s[i]], dtype=int)
        else:
            self.active_s = np.array([], dtype=bool)
            self.active_cc_s = np.array([], dtype=int)

    def _filter_reactive_interfaces(self):
        """
        过滤掉没有离子和电子通路的界面。

        对应 DERIVATION.md §6.4: "Verify every reactive interface has
        both ionic and electronic connectivity."
        """
        self.reactive_interfaces = []
        for e_g, s_g, area in self.interfaces:
            e_l = self.e_map[e_g]
            s_l = self.s_map[s_g]
            # 仅保留电解质和固相都活性的界面
            if self.active_e[e_l] and self.active_s[s_l]:
                self.reactive_interfaces.append((e_g, s_g, area))

        # 预计算界面数组 (用于向量化)
        self._ie_e = np.array(
            [self.e_map[ie[0]] for ie in self.reactive_interfaces],
            dtype=int,
        )  # 界面对应的电解质局部索引
        self._ie_s = np.array(
            [self.s_map[ie[1]] for ie in self.reactive_interfaces],
            dtype=int,
        )  # 界面对应的固相局部索引
        self._ie_A = np.array([ie[2] for ie in self.reactive_interfaces], dtype=float)
        self._ie_e_g = np.array([ie[0] for ie in self.reactive_interfaces], dtype=int)
        self._ie_s_g = np.array([ie[1] for ie in self.reactive_interfaces], dtype=int)

    def solve(self, I_app: float = 0.0, tol: float = 1e-12, max_iter: int = 100) -> dict:
        """
        使用逐次超松弛 (SOR) 求解线性化 BV 系统。

        对应 DERIVATION.md §5.1, §5.3, §5.6。

        每次迭代:
        1. 在当前电位下计算 BV 线性化电导 g_bv = dI/deta * A_r
        2. 构建包含 BV Jacobian 的耦合算子
        3. 求解线性系统

        Butler-Volmer 线性化 (对应 DERIVATION.md §5.2):
            i_r ≈ i_r^0 + (dI/deta) * (eta - eta^0)
            其中 dI/deta = i0 * F/(RT) * (alpha_a * E_a + alpha_c * E_c)

        Parameters
        ----------
        I_app : float
            施加电流密度 [A/m²]。
            I_app > 0: 阳极 (脱锂/充电)
            I_app < 0: 阴极 (嵌锂/放电)
        tol : float
            收敛容差 [V]。
        max_iter : int
            最大迭代次数。

        Returns
        -------
        result : dict
            - phi_e: 电解质电位 [V] (Np 数组, 非电解质节点为 NaN)
            - phi_s: 固相电位 [V] (Np 数组, 非固相节点为 NaN)
            - voltage: 电池电压 [V]
            - I_rxn: 反应电流密度 [A/m²] (Nt 数组)
            - iterations: 迭代次数
            - converged: 是否收敛
        """
        n_e, n_s = self.n_e, self.n_s
        n_intf = len(self.reactive_interfaces)

        # ===== Dirichlet 初始值 =====
        # 电解质: phi_e = 0 V (隔膜端参考电位)
        phi_e_bc = np.zeros(n_e)
        # 固相: phi_s = U(x) (平衡态, OCV)
        phi_s_bc = np.array([
            nmc532_ocv(self.c_s[self.s_indices[i]] / 48900.0) for i in range(n_s)
        ])

        phi_e = phi_e_bc.copy()
        phi_s = phi_s_bc.copy()

        # ===== 零电流快速返回 =====
        # 对应 DERIVATION.md §7.1: "Zero-Current Equilibrium (0C)"
        if abs(I_app) < tol:
            phi_e_full = np.full(self.Np, np.nan)
            phi_s_full = np.full(self.Np, np.nan)
            for g in range(self.Np):
                if self.e_mask[g]:
                    phi_e_full[g] = phi_e[self.e_map[g]]
                if self.solid_mask[g]:
                    phi_s_full[g] = phi_s[self.s_map[g]]

            # V_cell = phi_s(collector) - phi_e(separator)
            # 对应 DERIVATION.md §4.5
            voltage_nodes = self.active_cc_s if len(self.active_cc_s) else self.cc_s
            if len(voltage_nodes) > 0 and len(self.sep_e) > 0:
                V_cell = float(np.mean(phi_s[voltage_nodes])) - float(np.mean(phi_e[self.sep_e]))
            else:
                V_cell = 0.0

            return {
                "phi_e": phi_e_full,
                "phi_s": phi_s_full,
                "voltage": V_cell,
                "I_rxn": np.zeros(self.Nt),
                "iterations": 0,
                "converged": True,
            }

        if self._phi_e_guess is not None and self._phi_e_guess.shape == (n_e,):
            phi_e = self._phi_e_guess.copy()
        if self._phi_s_guess is not None and self._phi_s_guess.shape == (n_s,):
            phi_s = self._phi_s_guess.copy()

        # ===== SOR 迭代 =====
        converged = False
        iteration = 0
        for iteration in range(max_iter):
            # --- 步骤 1: 计算界面 BV 量 ---
            if n_intf > 0:
                phi_e_intf = phi_e[self._ie_e]  # 界面处电解质电位
                phi_s_intf = phi_s[self._ie_s]  # 界面处固相电位
                cs = self.c_s[self._ie_s_g]     # 界面处固相浓度
                ce = self.c_e[self._ie_e_g]     # 界面处电解质浓度
                x_j = cs / 48900.0              # 嵌锂度

                # 平衡电位
                U_eq = nmc532_ocv(x_j)
                # 交换电流密度
                i0 = np.array([
                    exchange_current_density(self.k0, ce[k], cs[k], 48900.0)
                    for k in range(n_intf)
                ])
                # 过电位: eta = phi_s - phi_e - U
                # 对应 DERIVATION.md §2.5
                eta = phi_s_intf - phi_e_intf - U_eq
                # BV 反应电流
                I_rxn = np.array([butler_volmer(i0[k], eta[k], self.T) for k in range(n_intf)])

                # 线性化电导: dI/deta
                # 对应 DERIVATION.md §5.2: B = d i/d eta = i0 f (alpha_a E_a + alpha_c E_c)
                f_val = F / (R * self.T)  # F/(RT) [V^{-1}]
                arg_a = np.clip(0.5 * f_val * eta, -500, 500)
                arg_c = np.clip(-0.5 * f_val * eta, -500, 500)
                dI_deta = i0 * (0.5 * f_val * np.exp(arg_a) + 0.5 * f_val * np.exp(arg_c))
                # g_bv: 类似电导的 BV 耦合项 [S = A/V]
                g_bv = dI_deta * self._ie_A  # [A/V * m² = S?]
            else:
                I_rxn = np.array([])
                g_bv = np.array([])

            # --- 步骤 2: 构建耦合系统 ---
            # 对应 DERIVATION.md §5.3: Jacobian Contributions from One Interface

            # 电解质块: (L_e - diag(g_bv_e)) @ phi_e + off-diag
            # 对应 DERIVATION.md §3.3: F_{phi_e} = Σ G(phi_k - phi_i) + Σ i_r A_r
            # 线性化后: i_r A_r ≈ g_bv * (phi_s - phi_e - U_eq)
            # 对 phi_e 的 Jacobian 贡献: -g_bv (对角线)
            # 对 phi_s 的 Jacobian 贡献: +g_bv (非对角线)

            # 汇总每个电解质节点的 g_bv
            g_bv_e = np.zeros(n_e)
            if n_intf > 0:
                np.add.at(g_bv_e, self._ie_e, g_bv)

            # 修改电解质 Laplacian: L_e - diag(g_bv_e)
            M_e = self.L_e.tolil()
            for i in range(n_e):
                M_e[i, i] -= g_bv_e[i]
            M_e = M_e.tocsr()

            # 电解质 RHS: -Σ(g_bv * (phi_s - U_eq))
            rhs_e = np.zeros(n_e)
            if n_intf > 0:
                bv_source = g_bv * (phi_s[self._ie_s] - U_eq)
                np.add.at(rhs_e, self._ie_e, -bv_source)

            # 固相块: (L_s + diag(g_bv_s)) @ phi_s - off-diag
            # 对应 DERIVATION.md §3.5: F_{phi_s} = Σ G(phi_n - phi_m) - Σ i_r A_r + F_BC
            # 线性化后: -i_r A_r ≈ -g_bv * (phi_s - phi_e - U_eq)
            # 对 phi_s 的 Jacobian 贡献: -g_bv (对角线)
            # 对 phi_e 的 Jacobian 贡献: +g_bv (非对角线)

            g_bv_s = np.zeros(n_s)
            if n_intf > 0:
                np.add.at(g_bv_s, self._ie_s, g_bv)

            # 修改固相 Laplacian: L_s - diag(g_bv_s)
            M_s = self.L_s.tolil()
            for i in range(n_s):
                M_s[i, i] -= g_bv_s[i]
            M_s = M_s.tocsr()

            # 固相 RHS: -Σ(g_bv * (phi_e + U_eq))
            rhs_s = np.zeros(n_s)
            if n_intf > 0:
                bv_sink = g_bv * (phi_e[self._ie_e] + U_eq)
                np.add.at(rhs_s, self._ie_s, -bv_sink)

            # 施加电流边界条件 (集流体端固相)
            # 对应 DERIVATION.md §4.2: F_{BC,m}^s = I_app * w_m
            # I_app 是电流密度 [A/m²], 需要乘以集流体面积
            if len(self.active_cc_s) > 0:
                rhs_s[self.active_cc_s] -= I_app * self._cc_area / len(self.active_cc_s)

            # ===== Dirichlet 边界条件 =====
            # 对应 DERIVATION.md §5.5: "replace the corresponding residual row by F_y = y - y_B = 0"

            # 电解质 Dirichlet: 隔膜端 phi_e = 0
            M_e = M_e.tolil()
            for i in self.sep_e:
                M_e[i, :] = 0; M_e[i, i] = 1.0
                rhs_e[i] = phi_e_bc[i]
            # 未连通的电解质节点也固定
            for i in range(n_e):
                if not self.active_e[i] or M_e[i, :].nnz == 0:
                    M_e[i, :] = 0
                    M_e[i, i] = 1.0; rhs_e[i] = 0.0
            M_e = M_e.tocsr()

            # 固相 Dirichlet: 未连通的节点固定到 OCV
            M_s = M_s.tolil()
            for i in range(n_s):
                if not self.active_s[i] or M_s[i, :].nnz == 0:
                    M_s[i, :] = 0
                    M_s[i, i] = 1.0; rhs_s[i] = phi_s_bc[i]
            M_s = M_s.tocsr()

            # ===== 求解线性系统 =====
            try:
                phi_e_new = spsolve(M_e, rhs_e)
            except Exception:
                phi_e_new = phi_e.copy()
            try:
                phi_s_new = spsolve(M_s, rhs_s)
            except Exception:
                phi_s_new = phi_s.copy()

            # ===== 欠松弛更新 (SOR) =====
            # 使用 alpha=0.5 的欠松弛提高稳定性
            alpha = 0.5
            de = np.max(np.abs(phi_e_new - phi_e))
            ds = np.max(np.abs(phi_s_new - phi_s))
            phi_e = (1 - alpha) * phi_e + alpha * phi_e_new
            phi_s = (1 - alpha) * phi_s + alpha * phi_s_new

            # 收敛检查
            if max(de, ds) < tol:
                converged = True
                break

        # ===== 构建全网络输出 =====
        if converged:
            self._phi_e_guess = phi_e.copy()
            self._phi_s_guess = phi_s.copy()

        phi_e_full = np.full(self.Np, np.nan)
        phi_s_full = np.full(self.Np, np.nan)
        for g in range(self.Np):
            if self.e_mask[g]:
                phi_e_full[g] = phi_e[self.e_map[g]]
            if self.solid_mask[g]:
                phi_s_full[g] = phi_s[self.s_map[g]]

        # 电池电压: V_cell = phi_s(collector) - phi_e(separator)
        # 对应 DERIVATION.md §4.5
        voltage_nodes = self.active_cc_s if len(self.active_cc_s) else self.cc_s
        if len(voltage_nodes) > 0 and len(self.sep_e) > 0:
            V_cell = float(np.mean(phi_s[voltage_nodes])) - float(np.mean(phi_e[self.sep_e]))
        else:
            V_cell = 0.0

        # 计算所有喉道的反应电流
        I_rxn_all = np.zeros(self.Nt)
        reactive_pairs = {(e_g, s_g) for e_g, s_g, _area in self.reactive_interfaces}
        for t in range(self.Nt):
            g1, g2 = int(self.conns[t, 0]), int(self.conns[t, 1])
            if self.e_mask[g1] and self.nmc_mask[g2] and (g1, g2) in reactive_pairs:
                le, lnmc = self.e_map[g1], self.s_map[g2]
                cs = self.c_s[g2]
                U_eq = nmc532_ocv(cs / 48900.0)
                i0 = exchange_current_density(self.k0, self.c_e[g1], cs, 48900.0)
                eta = phi_s[lnmc] - phi_e[le] - U_eq
                I_rxn_all[t] = butler_volmer(i0, eta, self.T)
            elif self.nmc_mask[g1] and self.e_mask[g2] and (g2, g1) in reactive_pairs:
                le, lnmc = self.e_map[g2], self.s_map[g1]
                cs = self.c_s[g1]
                U_eq = nmc532_ocv(cs / 48900.0)
                i0 = exchange_current_density(self.k0, self.c_e[g2], cs, 48900.0)
                eta = phi_s[lnmc] - phi_e[le] - U_eq
                I_rxn_all[t] = butler_volmer(i0, eta, self.T)

        return {
            "phi_e": phi_e_full,
            "phi_s": phi_s_full,
            "voltage": V_cell,
            "I_rxn": I_rxn_all,
            "iterations": iteration + 1,
            "converged": converged,
        }
