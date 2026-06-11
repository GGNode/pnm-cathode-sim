"""
P2D 结果存储与输出 (Results Storage)
====================================

职责:
- 保存 transient snapshots。
- 提供 voltage、capacity、平均浓度、诊断信息。

对应 TASK_PHASE0.md §1.9
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from pnmcathode.p2d.state import P2DState

Array = np.ndarray


@dataclass
class P2DSnapshot:
    """单个时间步的快照。

    Attributes
    ----------
    time : float
        时间 [s]。
    state : P2DState
        物理状态。
    voltage : float
        电池电压 [V]。
    current_density : float
        电流密度 [A/m²]。
    diagnostics : dict
        诊断信息。
    """

    time: float
    state: P2DState
    voltage: float
    current_density: float
    diagnostics: dict[str, float] = field(default_factory=dict)


@dataclass
class P2DResult:
    """P2D 仿真结果。

    Attributes
    ----------
    snapshots : list of P2DSnapshot
        时间步快照列表。
    metadata : dict
        仿真元数据。
    """

    snapshots: list[P2DSnapshot]
    metadata: dict[str, float | str]

    @property
    def time(self) -> Array:
        """时间数组 [s]。"""
        return np.array([s.time for s in self.snapshots])

    @property
    def voltage(self) -> Array:
        """电压数组 [V]。"""
        return np.array([s.voltage for s in self.snapshots])

    @property
    def capacity_Ah_m2(self) -> Array:
        """累积容量 [A·h/m²]。

        Q = ∫ |I_app| dt
        """
        times = self.time
        if len(times) < 2:
            return np.zeros(len(times))
        I_app = abs(self.snapshots[0].current_density) if self.snapshots else 0.0
        return I_app * times / 3600.0

    def final_state(self) -> P2DState:
        """返回最终状态。"""
        return self.snapshots[-1].state

    def to_npz(self, path: str) -> None:
            """保存结果到 .npz 文件。

            Parameters
            ----------
            path : str
                文件路径。
            """
            np.savez(
                path,
                time=self.time,
                voltage=self.voltage,
                capacity=self.capacity_Ah_m2,
            )
