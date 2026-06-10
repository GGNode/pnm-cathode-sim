"""孔网络电池模型可视化工具模块

本模块提供 PNM 模拟结果的可视化功能，包括：
- 放电曲线（V-Q 曲线）绘制
- 孔网络相分布图（电解液/NMC/CBD 三相空间分布）
- 浓度场空间分布图（电解液 Li⁺ 浓度或固相 Li 浓度）

可视化对应论文 Figure 5（放电曲线）和 Figure 6-8（空间分布）。

物理背景：
- 放电曲线是电池性能的核心表征：电压随容量（放电深度）下降
- 空间分布揭示了微结构异质性对性能的影响：
  · 高倍率下电解液浓度从隔膜到集流体递减（Li⁺ 耗尽）
  · 固相锂化状态(SoL)在小颗粒中更快达到饱和
  · 电流密度在集流体附近可能不均匀

使用方法：
    from pnmcathode.post.visualization import plot_discharge_curve, plot_concentration_field
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    plot_discharge_curve(result, ax=axes[0])
    plot_concentration_field(net, c_e, title="电解液浓度", ax=axes[1])
    plot_concentration_field(net, c_s, title="固相 Li 浓度", ax=axes[2])

依赖：
    - matplotlib（可选，缺少时 raise ImportError）
    - numpy

对应 DERIVATION.md：无直接对应章节（后处理/可视化）
"""

import numpy as np
import openpnm as op

# matplotlib 可选依赖 — 无 matplotlib 时仍可运行模拟，只是不能画图
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 — 导入以启用 3D 投影
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def plot_discharge_curve(result: dict, ax=None, **kwargs):
    """绘制放电曲线（电压 vs 容量）

    对应论文 Figure 5：不同 C-rate 下的 V-Q 曲线。
    放电曲线的特征：
    - 初始电压 ≈ OCV（开路电位，约 3.7V for NMC532 at x=0.5）
    - 中段电压缓慢下降（欧姆极化 + 活化极化）
    - 末段电压急剧下降（电解液 Li⁺ 耗尽 + 固相扩散限制）
    - 截止电压通常为 3.0V（论文）或 2.5V（安全余量）

    Parameters
    ----------
    result : dict
        放电模拟结果，来自 TransientSolver.run_discharge() 或 analysis.discharge_curve()。
        必须包含：
        - "capacity": 容量数组 [A·s/m²]（面积比容量）
        - "voltage": 电压数组 [V]（端电压）
    ax : matplotlib.axes.Axes or None
        绘图坐标轴。None 时自动创建新图。
    **kwargs
        传递给 ax.plot() 的额外参数（如 color, linewidth, label）。

    Returns
    -------
    ax : matplotlib.axes.Axes
        绘图坐标轴，可用于叠加多条曲线。

    示例：
        # 运行 0.2C 放电并画图
        result = solver.run_discharge(C_rate=0.2)
        plot_discharge_curve(result, label="0.2C", color="blue")

        # 叠加 1C 放电
        result_1c = solver.run_discharge(C_rate=1.0)
        plot_discharge_curve(result_1c, ax=ax, label="1C", color="red")
    """
    if not HAS_MPL:
        raise ImportError("matplotlib is required for plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(result["capacity"], result["voltage"], **kwargs)
    ax.set_xlabel("Capacity (A·s/m²)")  # 面积比容量 [A·s/m²] = [C/m²]
    ax.set_ylabel("Cell Voltage (V)")    # 端电压 [V] vs Li/Li+
    ax.set_title("Discharge Curve")
    ax.grid(True, alpha=0.3)
    return ax


def plot_phase_map(
    net: op.network.Cubic,
    property_name: str = "pore.phase_label",
    ax=None,
):
    """绘制孔网络相分布的 2D 投影图

    将三相孔网络（电解液=0, NMC=1, CBD=2）投影到 x-y 平面，
    用颜色区分不同相。对应论文 Figure 3（微结构可视化）。

    物理意义：
    - 电解液相（蓝色）：传输 Li⁺ 的通道
    - NMC 相（黄色）：存储 Li 的活性材料
    - CBD 相（绿色）：传导电子的碳粘结剂域
    - 三相的空间分布决定了局部反应活性和传输效率

    Parameters
    ----------
    net : openpnm.network.Cubic
        孔网络对象，包含 "pore.coords" 和相标签属性。
    property_name : str
        用于着色的孔属性名，默认 "pore.phase_label"。
        也可以是 "pore.volume"（孔体积）、"pore.diameter"（孔径）等。
    ax : matplotlib.axes.Axes or None
        绘图坐标轴。

    Returns
    -------
    ax : matplotlib.axes.Axes

    示例：
        plot_phase_map(net, "pore.phase_label")  # 三相分布
        plot_phase_map(net, "pore.volume")        # 孔体积分布
    """
    if not HAS_MPL:
        raise ImportError("matplotlib is required for plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    # 取孔坐标（3D）投影到 x-y 平面
    coords = net["pore.coords"]
    values = net[property_name]

    # 散点图：x-y 坐标，颜色映射属性值
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=values, cmap="viridis", s=20)
    ax.set_xlabel("x (m)")   # x 方向：隔膜→集流体
    ax.set_ylabel("y (m)")   # y 方向：横向
    ax.set_title(f"{property_name}")
    ax.set_aspect("equal")    # 等比例显示，保持物理几何关系
    plt.colorbar(sc, ax=ax)
    return ax


def plot_concentration_field(
    net: op.network.Cubic,
    c: np.ndarray,
    title: str = "Concentration",
    ax=None,
):
    """绘制浓度场的 2D 空间分布图

    对应论文 Figure 7-8：电解液浓度和固相锂化状态的空间分布。
    这些图揭示了高倍率下的空间异质性：
    - 电解液浓度：隔膜侧高（Li⁺ 源），集流体侧低（Li⁺ 耗尽）
    - 固相浓度：颗粒表面高（锂化快），核心低（扩散限制）
    - 小颗粒比大颗粒锂化更快（表面积/体积比更大）

    使用 "RdYlBu_r" 色图（红-黄-蓝反转）：
    - 红色 = 高浓度（电解液充足 / 固相锂化充分）
    - 蓝色 = 低浓度（电解液耗尽 / 固相贫锂）

    Parameters
    ----------
    net : openpnm.network.Cubic
        孔网络对象。
    c : ndarray, shape (Np,)
        每个孔节点的浓度值 [mol/m³]。
        对于电解液：c_e（典型范围 0-1200 mol/m³）
        对于固相：c_s（典型范围 0-48900 mol/m³）
        非活性孔（如 CBD）的值应为 NaN，绘图时自动跳过。
    title : str
        图标题，如 "电解液 Li⁺ 浓度" 或 "固相锂化状态"。
    ax : matplotlib.axes.Axes or None
        绘图坐标轴。

    Returns
    -------
    ax : matplotlib.axes.Axes

    示例：
        # 电解液浓度分布
        c_e_full = np.full(net.Np, np.nan)
        c_e_full[e_mask] = solver.c_e
        plot_concentration_field(net, c_e_full, "电解液 Li⁺ 浓度 [mol/m³]")

        # 固相锂化状态 (SoL = c_s / c_s_max)
        sol_full = np.full(net.Np, np.nan)
        sol_full[nmc_mask] = solver.c_s[nmc_mask] / 48900.0
        plot_concentration_field(net, sol_full, "固相锂化状态 (SoL)")
    """
    if not HAS_MPL:
        raise ImportError("matplotlib is required for plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    coords = net["pore.coords"]
    # 过滤 NaN 值（非活性孔，如 CBD 不参与浓度计算）
    valid = ~np.isnan(c)

    # 散点图：红-黄-蓝色图，高浓度=红，低浓度=蓝
    sc = ax.scatter(
        coords[valid, 0], coords[valid, 1],
        c=c[valid], cmap="RdYlBu_r", s=20,
    )
    ax.set_xlabel("x (m)")   # x 方向：隔膜→集流体
    ax.set_ylabel("y (m)")   # y 方向：横向
    ax.set_title(title)
    ax.set_aspect("equal")
    plt.colorbar(sc, ax=ax, label="mol/m³")  # 浓度色标 [mol/m³]
    return ax
