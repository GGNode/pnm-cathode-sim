"""后处理模块: 可视化与分析。

子模块:
- visualization: 放电曲线、相分布图、浓度场图
- analysis: 孔径分布、配位数、锂化度、放电曲线生成
"""

from pnmcathode.post.analysis import (
    coordination_number,
    discharge_curve,
    pore_size_distribution,
    state_of_lithiation,
)
from pnmcathode.post.visualization import (
    plot_concentration_field,
    plot_discharge_curve,
    plot_phase_map,
)

__all__ = [
    "coordination_number",
    "discharge_curve",
    "plot_concentration_field",
    "plot_discharge_curve",
    "plot_phase_map",
    "pore_size_distribution",
    "state_of_lithiation",
]
