"""仿真结果容器与序列化辅助函数。"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


def _metadata_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


@dataclass
class DischargeResult:
    """恒流放电输出数组的类型化包装。"""

    time: np.ndarray
    voltage: np.ndarray
    capacity_Ah_m2: np.ndarray
    current_density: float
    c_rate: float
    reached_cutoff: bool
    spatial: dict[str, np.ndarray] | None = None
    metadata: dict[str, Any] | None = None

    @property
    def final_capacity(self) -> float:
        """最终面积比容量 [A h m^-2]。"""

        if self.capacity_Ah_m2.size == 0:
            return 0.0
        return float(self.capacity_Ah_m2[-1])

    @classmethod
    def from_solver_output(
        cls,
        raw: dict[str, Any],
        spatial: dict[str, np.ndarray] | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> "DischargeResult":
        """从 ``TransientSolver.run_discharge`` 输出构建结果对象。"""

        merged_metadata: dict[str, Any] = {
            "cutoff_voltage": raw.get("cutoff_voltage"),
        }
        if metadata:
            merged_metadata.update(metadata)
        return cls(
            time=np.asarray(raw["time"], dtype=float),
            voltage=np.asarray(raw["voltage"], dtype=float),
            capacity_Ah_m2=np.asarray(raw["capacity_Ah_m2"], dtype=float),
            current_density=float(raw["I_app"]),
            c_rate=float(raw["C_rate"]),
            reached_cutoff=bool(raw["reached_cutoff"]),
            spatial=spatial,
            metadata=merged_metadata,
        )

    def to_npz(self, path: str | Path) -> None:
        """将结果序列化为压缩 NPZ 文件。"""

        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)

        spatial = self.spatial or {}
        spatial_keys = sorted(spatial)
        metadata_json = json.dumps(
            self.metadata or {},
            default=_metadata_default,
            sort_keys=True,
        )
        payload: dict[str, Any] = {
            "time": np.asarray(self.time, dtype=float),
            "voltage": np.asarray(self.voltage, dtype=float),
            "capacity_Ah_m2": np.asarray(self.capacity_Ah_m2, dtype=float),
            "current_density": np.asarray(self.current_density, dtype=float),
            "c_rate": np.asarray(self.c_rate, dtype=float),
            "reached_cutoff": np.asarray(self.reached_cutoff, dtype=bool),
            "spatial_keys": np.asarray(spatial_keys, dtype=str),
            "metadata_json": np.asarray(metadata_json, dtype=str),
        }
        for key in spatial_keys:
            payload[f"spatial.{key}"] = np.asarray(spatial[key])

        np.savez_compressed(target, **payload)

    @classmethod
    def from_npz(cls, path: str | Path) -> "DischargeResult":
        """从 :meth:`to_npz` 写出的 NPZ 文件反序列化结果。"""

        with np.load(Path(path), allow_pickle=False) as data:
            spatial_keys = [str(key) for key in data["spatial_keys"].tolist()]
            spatial = {
                key: np.asarray(data[f"spatial.{key}"])
                for key in spatial_keys
            }
            metadata_json = str(data["metadata_json"].item())
            metadata = json.loads(metadata_json) if metadata_json else {}
            return cls(
                time=np.asarray(data["time"], dtype=float),
                voltage=np.asarray(data["voltage"], dtype=float),
                capacity_Ah_m2=np.asarray(data["capacity_Ah_m2"], dtype=float),
                current_density=float(data["current_density"].item()),
                c_rate=float(data["c_rate"].item()),
                reached_cutoff=bool(data["reached_cutoff"].item()),
                spatial=spatial or None,
                metadata=metadata,
            )
