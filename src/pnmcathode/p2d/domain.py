"""
P2D 模型的区域定义与有限体积网格 (Domain & Mesh)
=================================================

职责:
- 定义 Phase 0 的 separator 和 positive electrode 区域。
- 构建 x 方向宏观有限体积网格 (MacroMesh)。
- 构建球形颗粒径向有限体积网格 (ParticleMesh)。
- 提供 face 几何和区域索引。

对应 TASK_PHASE0.md §1.2
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

Array = np.ndarray
RegionName = Literal["negative", "separator", "positive"]


@dataclass(frozen=True)
class P2DRegion:
    """P2D 模型中的一个均匀区域 (separator 或 positive electrode)。

    属性:
        name: 区域名称 ("negative", "separator", "positive")
        x_left: 区域左边界坐标 [m]
        x_right: 区域右边界坐标 [m]
        n_cells: 该区域的有限体积单元数
        epsilon_e: 电解质体积分数 (孔隙率)
        epsilon_s: 固相活性材料体积分数 (仅 positive electrode)
        particle_radius: 颗粒半径 [m] (仅 positive electrode)
        specific_area: 外部给定的比表面积 [m²/m³]；None 时用球形颗粒公式
        bruggeman_e: 电解质 Bruggeman 指数
        bruggeman_s: 固相 Bruggeman 指数
        sigma_s: 固相电子电导率 [S/m]；None 表示无固相导电
    """

    name: RegionName
    x_left: float
    x_right: float
    n_cells: int
    epsilon_e: float
    epsilon_s: float = 0.0
    particle_radius: float | None = None
    specific_area: float | None = None
    bruggeman_e: float = 1.5
    bruggeman_s: float = 1.5
    sigma_s: float | None = None

    def __post_init__(self) -> None:
        if self.n_cells <= 0:
            raise ValueError("n_cells must be positive")
        if self.x_right <= self.x_left:
            raise ValueError("x_right must be greater than x_left")
        if not 0.0 < self.epsilon_e <= 1.0:
            raise ValueError("epsilon_e must be in (0, 1]")

    @property
    def length(self) -> float:
        """区域厚度 [m]。"""
        return self.x_right - self.x_left

    def area_density(self) -> float:
        """计算比表面积 a_s [m²/m³]。

        优先使用外部给定的 specific_area；
        否则使用球形颗粒公式: a_s = 3 * epsilon_s / R_p。

        对应 TASK_PHASE0.md 公式 2.1
        """
        if self.specific_area is not None:
            return self.specific_area
        if self.particle_radius is None or self.particle_radius <= 0:
            raise ValueError(
                "particle_radius must be positive when specific_area is not provided"
            )
        return 3.0 * self.epsilon_s / self.particle_radius


@dataclass(frozen=True)
class MacroMesh:
    """x 方向宏观有限体积网格。

    覆盖 separator 和 positive electrode 区域。
    区域顺序: [separator cells, positive electrode cells]。

    对应 TASK_PHASE0.md §4.1
    """

    regions: tuple[P2DRegion, ...]
    x_faces: Array
    x_centers: Array
    dx: Array
    volumes: Array
    face_areas: Array
    region_id: Array
    separator_cells: Array
    positive_cells: Array
    electrode_cells: Array

    @classmethod
    def from_regions(
        cls,
        regions: tuple[P2DRegion, ...],
        area: float = 1.0,
    ) -> "MacroMesh":
        """从区域列表构建宏观网格。

        Parameters
        ----------
        regions : tuple of P2DRegion
            区域列表，按从左到右顺序排列。
        area : float
            几何截面积 [m²]，默认 1.0 m²。

        Returns
        -------
        MacroMesh
        """
        if area <= 0:
            raise ValueError("area must be positive")

        n_x = sum(r.n_cells for r in regions)
        x_faces = np.zeros(n_x + 1)
        x_centers = np.zeros(n_x)
        dx_arr = np.zeros(n_x)
        volumes = np.zeros(n_x)
        face_areas = np.full(n_x + 1, area)
        region_id = np.zeros(n_x, dtype=int)

        idx = 0
        for rid, reg in enumerate(regions):
            faces = np.linspace(reg.x_left, reg.x_right, reg.n_cells + 1)
            for i in range(reg.n_cells):
                x_faces[idx] = faces[i]
                x_centers[idx] = 0.5 * (faces[i] + faces[i + 1])
                dx_arr[idx] = faces[i + 1] - faces[i]
                volumes[idx] = dx_arr[idx] * area
                region_id[idx] = rid
                idx += 1
        x_faces[n_x] = regions[-1].x_right

        # 构建区域索引
        sep_ids = [
            i for i, r in enumerate(regions) if r.name == "separator"
        ]
        pos_ids = [
            i for i, r in enumerate(regions) if r.name == "positive"
        ]

        separator_cells = np.array(
            [i for i in range(n_x) if region_id[i] in sep_ids], dtype=int
        )
        positive_cells = np.array(
            [i for i in range(n_x) if region_id[i] in pos_ids], dtype=int
        )
        electrode_cells = positive_cells.copy()

        return cls(
            regions=regions,
            x_faces=x_faces,
            x_centers=x_centers,
            dx=dx_arr,
            volumes=volumes,
            face_areas=face_areas,
            region_id=region_id,
            separator_cells=separator_cells,
            positive_cells=positive_cells,
            electrode_cells=electrode_cells,
        )

    @property
    def n_x(self) -> int:
        """总单元数。"""
        return len(self.dx)

    @property
    def n_pos(self) -> int:
        """positive electrode 单元数。"""
        return len(self.positive_cells)

    def region_for_cell(self, cell: int) -> P2DRegion:
        """返回指定单元所属的区域。"""
        return self.regions[self.region_id[cell]]

    def left_face(self, cell: int) -> int:
        """返回单元左侧面的索引。"""
        return cell

    def right_face(self, cell: int) -> int:
        """返回单元右侧面的索引。"""
        return cell + 1


@dataclass(frozen=True)
class ParticleMesh:
    """球形颗粒径向有限体积网格。

    每个 positive macro cell 共享同一网格结构。
    r=0 处的 face area 为 0，自然满足中心对称条件。

    对应 TASK_PHASE0.md §4.2
    """

    radius: float
    n_shells: int
    r_faces: Array
    r_centers: Array
    shell_volumes: Array
    face_areas: Array

    @classmethod
    def spherical(cls, radius: float, n_shells: int) -> "ParticleMesh":
        """构建均匀球壳网格。

        Parameters
        ----------
        radius : float
            颗粒半径 [m]。
        n_shells : int
            壳层数量。

        Returns
        -------
        ParticleMesh
        """
        if radius <= 0:
            raise ValueError("radius must be positive")
        if n_shells <= 0:
            raise ValueError("n_shells must be positive")

        r_faces = np.linspace(0.0, radius, n_shells + 1)
        r_centers = 0.5 * (r_faces[:-1] + r_faces[1:])

        # 壳层体积: V_i = (4/3) * pi * (r_{i+1}^3 - r_i^3)
        shell_volumes = (4.0 / 3.0) * np.pi * (
            r_faces[1:] ** 3 - r_faces[:-1] ** 3
        )

        # 面积: A_i = 4 * pi * r_i^2；r=0 处为 0 (中心对称)
        face_areas = 4.0 * np.pi * r_faces ** 2
        face_areas[0] = 0.0  # 中心对称: 无通量

        return cls(
            radius=radius,
            n_shells=n_shells,
            r_faces=r_faces,
            r_centers=r_centers,
            shell_volumes=shell_volumes,
            face_areas=face_areas,
        )
