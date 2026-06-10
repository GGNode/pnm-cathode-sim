"""
孔网络生成器模块 (Pore Network Generator)
==========================================

物理背景
--------
本模块创建阴极电极的立方孔网络 (cubic pore network)。

孔网络模型 (PNM, Pore Network Model) 将多孔电极离散化为:
- 节点 (node/pore): 代表孔隙体 (pore body), 存储物质和电位
- 键 (bond/throat): 代表喉道 (pore throat), 连接相邻孔隙, 传导通量

三相标记
--------
阴极电极包含三相:
- Phase 0 (电解质, electrolyte): 孔隙空间, 充满电解液
    → 承载 c_e (电解质 Li+ 浓度) 和 phi_e (电解质电位)
- Phase 1 (NMC532, 活性材料): 锂嵌入/脱锂的活性位点
    → 承载 c_s (固相 Li 浓度) 和 phi_s (固相电位)
    → 与电解质通过 Butler-Volmer 反应耦合
- Phase 2 (CBD, 碳粘结剂域): 导电但不存储锂
    → 仅承载 phi_s (电子导电), c_s = 0

相标记通过随机分配实现, 使得:
- 电解质体积分数 ≈ porosity (孔隙率)
- CBD 体积分数 ≈ cbd_fraction
- NMC 体积分数 ≈ 1 - porosity - cbd_fraction

几何属性
--------
每个节点 (pore):
- pore.diameter: 孔径 [m], 均匀分布
- pore.volume: 孔体积 [m³], 按球体计算: V = π/6 * d³

每个键 (throat):
- throat.diameter: 喉道直径 [m]
- throat.area: 截面积 [m²], A = π/4 * d²
- throat.length: 喉道长度 [m], 等于中心距减去两端孔半径

Khan et al. (2021) 的网络规模:
- 1CAL: 4637 节点, 31427 键
- 3CAL: 3510 节点, 23126 键
- 孔隙率: ~35%

对应 DERIVATION.md 章节: §1.1, §3.1, §8

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op
from scipy import sparse
from scipy.sparse.csgraph import connected_components


def check_percolation(net: op.network.Cubic, phase: str = "electrolyte") -> bool:
    """Return True when a phase has a connected path from separator to collector.

    The x direction is treated as the electrode thickness direction: the
    minimum-x plane is the separator side, and the maximum-x plane is the
    current-collector side.
    """
    phase_key = phase.lower()
    if phase_key in {"electrolyte", "e"}:
        mask = np.asarray(net["pore.electrolyte"], dtype=bool)
    elif phase_key in {"solid", "s"}:
        mask = np.asarray(net["pore.nmc"], dtype=bool) | np.asarray(net["pore.cbd"], dtype=bool)
    elif phase_key == "nmc":
        mask = np.asarray(net["pore.nmc"], dtype=bool)
    elif phase_key == "cbd":
        mask = np.asarray(net["pore.cbd"], dtype=bool)
    else:
        raise ValueError("phase must be one of: electrolyte, solid, nmc, cbd")

    if not np.any(mask):
        return False

    coords = np.asarray(net["pore.coords"], dtype=float)
    x = coords[:, 0]
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    sep = mask & np.isclose(x, x_min)
    cc = mask & np.isclose(x, x_max)
    if not np.any(sep) or not np.any(cc):
        return False

    phase_pores = np.where(mask)[0]
    phase_map = np.full(net.Np, -1, dtype=int)
    phase_map[phase_pores] = np.arange(phase_pores.size)

    conns = np.asarray(net["throat.conns"], dtype=int)
    keep = mask[conns[:, 0]] & mask[conns[:, 1]]
    if not np.any(keep):
        return False

    local_conns = phase_map[conns[keep]]
    rows = np.r_[local_conns[:, 0], local_conns[:, 1]]
    cols = np.r_[local_conns[:, 1], local_conns[:, 0]]
    data = np.ones(rows.size, dtype=bool)
    graph = sparse.csr_matrix((data, (rows, cols)), shape=(phase_pores.size, phase_pores.size))
    _n_components, labels = connected_components(graph, directed=False)

    sep_components = set(labels[phase_map[np.where(sep)[0]]])
    cc_components = set(labels[phase_map[np.where(cc)[0]]])
    return bool(sep_components & cc_components)


def create_cathode_network(
    shape: list[int] = [10, 10, 10],
    spacing: float = 1e-5,
    porosity: float = 0.35,
    cbd_fraction: float = 0.10,
    seed: int | None = None,
    throat_scale: float = 1.0,
) -> op.network.Cubic:
    """
    创建具有三相阴极标记的立方孔网络。

    对应 DERIVATION.md §1.1 "Network Sets":
    - 电解质图: 节点 i ∈ E, 喉道 (i,k) ∈ T_e
    - 固相图: 活性材料/CBD 节点 m ∈ S, 固相接触 (m,n) ∈ T_s
    - 电解质/活性材料界面: 反应 r ∈ R

    三相随机分配策略:
    - 对每个节点生成 [0,1) 的随机数
    - 若 < porosity → 电解质 (Phase 0)
    - 若 ∈ [porosity, porosity+cbd_fraction) → CBD (Phase 2)
    - 否则 → NMC532 (Phase 1)

    Parameters
    ----------
    shape : list of int
        网络维度 [nx, ny, nz], 每个方向的节点数。
        例如 [10, 10, 10] → 1000 个节点。
        nx 方向对应电极厚度方向 (x=0 为隔膜端, x=L 为集流体端)。
    spacing : float
        节点间距 (pore-to-pore spacing) [m]。
        典型值: 1e-5 m (10 μm)。
        对应 DERIVATION.md §1.1 中的几何量。
    porosity : float
        目标电解质体积分数 (Phase 0)。
        典型值: 0.35 (35% 孔隙率)。
    cbd_fraction : float
        CBD 体积分数 (Phase 2)。
        典型值: 0.10 (10%)。
        NMC 体积分数 = 1 - porosity - cbd_fraction。
    seed : int or None
        随机种子, 用于结果复现。None 表示不固定种子。

    Returns
    -------
    net : openpnm.network.Cubic
        带有以下属性的孔网络:
        - pore.phase_label: 相标记 (0=电解质, 1=NMC, 2=CBD)
        - pore.electrolyte: 布尔数组, 电解质节点
        - pore.nmc: 布尔数组, NMC 节点
        - pore.cbd: 布尔数组, CBD 节点
        - pore.diameter: 孔径 [m]
        - pore.volume: 孔体积 [m³]
        - throat.diameter: 喉道直径 [m]
        - throat.area: 喉道截面积 [m²]
        - throat.length: 喉道长度 [m]
        - pore.cs_max: NMC 最大锂浓度 (仅 NMC 节点非零) [mol/m³]
        - pore.sigma_nmc: NMC 电导率 (仅 NMC 节点非零) [S/m]
        - pore.sigma_cbd: CBD 电导率 (仅 CBD 节点非零) [S/m]
    """
    rng = np.random.default_rng(seed)

    # 创建立方网络拓扑
    # OpenPNM 会自动创建 throat.conns (连接矩阵) 和 pore.coords (坐标)
    net = op.network.Cubic(shape=shape, spacing=spacing)

    # ===== 相标记分配 =====
    # 对应 DERIVATION.md §1.1 的三相: E (电解质), S (活性固相), CBD
    n_pores = net.Np
    labels = np.zeros(n_pores, dtype=int)

    # 随机分配: 基于目标体积分数
    rand_vals = rng.random(n_pores)
    labels[rand_vals < porosity] = 0                        # Phase 0: 电解质
    labels[(rand_vals >= porosity) & (rand_vals < porosity + cbd_fraction)] = 2  # Phase 2: CBD
    labels[rand_vals >= porosity + cbd_fraction] = 1         # Phase 1: NMC532

    # 存储相标记到网络属性
    net["pore.phase_label"] = labels
    net["pore.electrolyte"] = labels == 0  # 电解质布尔掩码
    net["pore.nmc"] = labels == 1          # NMC 布尔掩码
    net["pore.cbd"] = labels == 2          # CBD 布尔掩码

    # ===== 几何属性 =====
    # 孔径: 在 spacing 的 [0.8, 1.2] 倍范围内均匀分布
    # 使用窄分布避免 SoL 双峰 (小孔太快饱和)
    pore_diameter = rng.uniform(0.8, 1.2, n_pores) * spacing
    net["pore.diameter"] = pore_diameter
    # 孔体积: 按球体计算, V = π/6 * d³
    net["pore.volume"] = (np.pi / 6.0) * pore_diameter**3

    # 喉道属性
    n_throats = net.Nt
    # 喉道直径: 在 spacing*0.5 的 [0.2, 0.8] 倍范围内
    # throat_scale > 1 增大接触面积, 补偿合成网络 vs XCT 微结构的几何差异
    throat_diameter = rng.uniform(0.2, 0.8, n_throats) * spacing * 0.5 * throat_scale
    net["throat.diameter"] = throat_diameter
    # 喉道截面积: A = π/4 * d²
    net["throat.area"] = (np.pi / 4.0) * throat_diameter**2

    # 喉道长度: 两孔中心距减去两端孔半径
    # 对应 DERIVATION.md §3.1: L_ik^e 是喉道长度
    conns = net["throat.conns"]  # 连接矩阵: conns[t] = [pore1, pore2]
    p1_coords = net["pore.coords"][conns[:, 0]]  # 孔1 坐标
    p2_coords = net["pore.coords"][conns[:, 1]]  # 孔2 坐标
    p1_radius = pore_diameter[conns[:, 0]] / 2.0  # 孔1 半径
    p2_radius = pore_diameter[conns[:, 1]] / 2.0  # 孔2 半径
    # 中心距
    center_dist = np.linalg.norm(p2_coords - p1_coords, axis=1)
    # 喉道长度 = 中心距 - 两端孔半径
    throat_length = center_dist - p1_radius - p2_radius
    # 确保最小长度 (避免零长度喉道)
    throat_length = np.maximum(throat_length, spacing * 0.01)
    net["throat.length"] = throat_length

    # ===== 相特有属性 =====
    # NMC532 最大锂浓度: c_s_max = 48900 mol/m³
    # 仅 NMC 节点有此属性 (Phase 1)
    net["pore.cs_max"] = np.where(labels == 1, 48900.0, 0.0)
    # NMC 电子导电率: sigma_nmc ≈ 0.01 S/m (较低)
    net["pore.sigma_nmc"] = np.where(labels == 1, 0.01, 0.0)  # S/m
    # CBD 电子导电率: sigma_cbd ≈ 760 S/m (较高, 导电碳)
    net["pore.sigma_cbd"] = np.where(labels == 2, 760.0, 0.0)  # S/m

    return net
