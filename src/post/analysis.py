"""
后处理与分析工具 (Post-Processing and Analysis)
================================================

物理背景
--------
本模块提供孔隙网络电池模型的结构分析和放电性能分析工具。

功能分类:
1. 网络结构分析:
   - 孔径分布 (pore size distribution): 反映多孔电极的微观结构
   - 配位数 (coordination number): 每个孔隙的连接数, 影响传输路径
   - 锂化度 (state of lithiation): x = c_s / c_s_max, 反映荷电状态

2. 放电性能分析:
   - 放电曲线 (discharge curve): V(t) 和 V(Q) 关系
   - 容量计算: Q = |I_app| * t [A·s/m²]

孔径分布
--------
孔径分布是多孔电极的基本结构特征:
- 影响比表面积 (a_s = 3/R_p for spheres)
- 影响渗透率和传输特性
- 与 Bruggeman 迂曲度修正相关

参见 DERIVATION.md §3.1 关于几何参数的定义。

配位数
------
配位数 (coordination number) 是每个孔隙节点连接的喉道数:
- 影响网络的连通性和传输效率
- 典型的多孔电极配位数约 4-6 (三维)
- 配位数过低可能导致传输瓶颈

锂化度
------
锂化度 (state of lithiation) 定义为:
    x = c_s / c_s_max

其中:
    x = 0: 完全脱锂 (charged state)
    x = 1: 完全嵌锂 (discharged state)

OCV 曲线 U(x) 直接由锂化度决定:
    V_cell ≈ U(x_surface) - |eta| - IR_drop

放电曲线
--------
放电曲线 (discharge curve) 描述电池电压随放电进行的变化:
    V(t) = U(x(t)) - eta(t) - IR_drop(t)

放电过程中的电压降来源:
    1. 活化过电位 eta: BV 动力学阻力
    2. 欧姆压降 IR: 电解质和固相电阻
    3. 浓度过降 Delta_V_conc: c_e 和 c_s 的不均匀分布

容量计算:
    Q(t) = |I_app| * t   [A·s/m²]
    Q_Ah(t) = Q(t) / 3600   [A·h/m²]

对应 DERIVATION.md 章节: §3.1, §7.2, §7.3

参考文献
--------
Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op

from src.solver.transient import TransientSolver


def pore_size_distribution(
    net: op.network.Cubic, n_bins: int = 20,
) -> tuple[np.ndarray, np.ndarray]:
    """
    计算孔径分布直方图。

    孔径分布是多孔电极的基本结构表征:
    - 影响比表面积: a_s ≈ 3/R_p (球形孔隙)
    - 影响有效传输系数: D_eff = D * epsilon^b (Bruggeman)
    - 与反应面积密度相关

    Parameters
    ----------
    net : openpnm.network.Cubic
        孔隙网络对象, 需包含 'pore.diameter' 属性。
    n_bins : int
        直方图的分箱数, 默认 20。

    Returns
    -------
    bins : ndarray, shape (n_bins + 1,)
        分箱边界 [m]。
    hist : ndarray, shape (n_bins,)
        各箱的孔隙数 (计数)。
    """
    diameters = net["pore.diameter"]
    hist, bins = np.histogram(diameters, bins=n_bins)
    return bins, hist


def coordination_number(net: op.network.Cubic) -> np.ndarray:
    """
    计算每个孔隙的配位数 (连接数)。

    配位数 = 与该孔隙相连的喉道数。
    对于三维立方网络, 内部节点配位数为 6, 表面节点较少。

    物理意义:
    - 高配位数 → 更多传输路径 → 更好的传输性能
    - 低配位数 → 传输瓶颈 → 可能限制反应均匀性
    - 配位数为 0 的节点是孤立节点 (不参与传输)

    Parameters
    ----------
    net : openpnm.network.Cubic
        孔隙网络对象, 需包含 'throat.conns' 属性。

    Returns
    -------
    cn : ndarray, shape (Np,), dtype int
        每个孔隙的配位数。
    """
    conns = net["throat.conns"]
    cn = np.zeros(net.Np, dtype=int)
    # 对每个喉道的两端节点各加 1
    np.add.at(cn, conns[:, 0], 1)
    np.add.at(cn, conns[:, 1], 1)
    return cn


def state_of_lithiation(c_s: np.ndarray, net: op.network.Cubic) -> np.ndarray:
    """
    计算每个孔隙的锂化度 (state of lithiation)。

    定义: x = c_s / c_s_max

    其中:
        x = 0: 完全脱锂 (空位满, charged state)
        x = 1: 完全嵌锂 (锂满, discharged state)

    仅 NMC 节点有物理意义; 非 NMC 节点 (电解质, CBD) 返回 NaN。

    锂化度与 OCV 的关系:
        V_cell ≈ U(x) = OCV_curve(x)
    NMC532 的 OCV 曲线 x 从 0 (4.2V) 到 1 (2.5V)。

    Parameters
    ----------
    c_s : ndarray, shape (Np,)
        固相锂浓度 [mol/m³]。所有孔隙的值, 但仅 NMC 节点有意义。
    net : openpnm.network.Cubic
        孔隙网络对象, 需包含 'pore.nmc' 属性。

    Returns
    -------
    sol : ndarray, shape (Np,)
        锂化度 [0, 1]。非 NMC 节点为 NaN。
    """
    cs_max = 48900.0  # NMC532 最大锂浓度 [mol/m³]
    sol = np.full(net.Np, np.nan)
    nmc_mask = net["pore.nmc"]
    # 裁剪到 [0, 1] 防止数值越界
    sol[nmc_mask] = np.clip(c_s[nmc_mask] / cs_max, 0.0, 1.0)
    return sol


def discharge_curve(
    net: op.network.Cubic,
    I_app: float = -0.001,
    dt: float = 10.0,
    n_steps: int = 100,
    T: float = 298.15,
    k0: float = 5e-10,
    c_e_init: float = 1200.0,
    c_s_init: float = 24450.0,
) -> dict:
    """
    运行瞬态放电仿真并返回电压-容量曲线。

    这是一个便捷函数, 封装了 TransientSolver 的典型使用流程:
    1. 创建瞬态求解器
    2. 设置初始浓度
    3. 循环推进时间步
    4. 记录 V(t) 和 Q(t)

    放电曲线的典型特征:
    - 初始电压降: 活化过电位建立
    - 中间平台: OCV 曲线相对平坦的区域
    - 末期急剧下降: 浓度极化 (c_e → 0 或 c_s → cs_max)

    容量计算:
        Q(t) = |I_app| * t   [A·s/m²]
        这是面积比容量 (areal capacity), 表示单位集流体面积的放电量。

    对应 DERIVATION.md §7.2 (小电流极限), §7.3 (高倍率耗尽)。

    Parameters
    ----------
    net : openpnm.network.Cubic
        孔隙网络对象。
    I_app : float
        应用电流密度 [A/m²] (阳极约定, 放电为负)。
        默认 -0.001 A/m² (非常小的电流, 接近 OCV)。
    dt : float
        时间步 [s], 默认 10.0。
    n_steps : int
        时间步数, 默认 100。总仿真时长 = dt * n_steps。
    T : float
        温度 [K], 默认 298.15。
    k0 : float
        BV 速率常数 [m^(2.5) / (mol^0.5 · s)]。
    c_e_init : float
        初始电解质浓度 [mol/m³]。
    c_s_init : float
        初始固相浓度 [mol/m³]。默认 24450 (x ≈ 0.5)。

    Returns
    -------
    result : dict
        voltage: 电压数组 [V], shape (n_steps,)
        capacity: 容量数组 [A·s/m²], shape (n_steps,)
        time: 时间数组 [s], shape (n_steps,)
    """
    # 创建瞬态求解器并设置初始浓度
    solver = TransientSolver(net, T=T, k0=k0)
    solver.set_concentration(c_e=c_e_init, c_s=c_s_init)

    voltages = np.zeros(n_steps)
    times = np.zeros(n_steps)
    capacities = np.zeros(n_steps)

    for i in range(n_steps):
        # 推进一个时间步
        result = solver.step(dt=dt, I_app=I_app)
        voltages[i] = result["voltage"]
        times[i] = (i + 1) * dt
        # 容量 = |I| * t  [A·s/m²] (面积比容量)
        capacities[i] = abs(I_app) * times[i]

    return {
        "voltage": voltages,
        "capacity": capacities,
        "time": times,
    }
