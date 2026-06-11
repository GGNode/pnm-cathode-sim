"""
P2D 状态向量管理 (State Vector Management)
==========================================

职责:
- 定义 monolithic Newton 向量的内存布局。
- 将物理状态 (P2DState) 与 Newton 向量 (Array) 互转。
- 保证 separator 中没有 phi_s 和 c_s unknown。

状态向量布局:
    y = [c_e(0:n_x), phi_e(0:n_x), phi_s(0:n_pos), c_s(0:n_pos*n_r)]

对应 TASK_PHASE0.md §1.3, §4.3
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pnmcathode.p2d.domain import MacroMesh, ParticleMesh

Array = np.ndarray


@dataclass(frozen=True)
class StateLayout:
    """描述 Newton 向量中各物理变量的切片位置。

    Attributes
    ----------
    n_x : int
        宏观单元总数 (separator + positive)。
    n_pos : int
        positive electrode 单元数。
    n_r : int
        颗粒壳层数。
    c_e : slice
        电解质浓度在 Newton 向量中的切片。
    phi_e : slice
        电解质电位在 Newton 向量中的切片。
    phi_s : slice
        固相电位在 Newton 向量中的切片。
    c_s : slice
        固相浓度 (展平) 在 Newton 向量中的切片。
    positive_to_macro : ndarray
        positive local index → macro cell index 映射。
    macro_to_positive : ndarray
        macro cell index → positive local index 映射 (-1 表示非 positive)。
    """

    n_x: int
    n_pos: int
    n_r: int
    c_e: slice
    phi_e: slice
    phi_s: slice
    c_s: slice
    positive_to_macro: Array
    macro_to_positive: Array

    @property
    def size(self) -> int:
        """Newton 向量总长度。"""
        return 2 * self.n_x + self.n_pos + self.n_pos * self.n_r

    @classmethod
    def from_meshes(cls, macro: MacroMesh, particle: ParticleMesh) -> "StateLayout":
        """从宏观网格和颗粒网格构建状态布局。

        Parameters
        ----------
        macro : MacroMesh
            宏观有限体积网格。
        particle : ParticleMesh
            颗粒径向网格。

        Returns
        -------
        StateLayout
        """
        n_x = macro.n_x
        n_pos = macro.n_pos
        n_r = particle.n_shells

        i0 = 0
        c_e = slice(i0, i0 + n_x)
        i0 += n_x
        phi_e = slice(i0, i0 + n_x)
        i0 += n_x
        phi_s = slice(i0, i0 + n_pos)
        i0 += n_pos
        c_s = slice(i0, i0 + n_pos * n_r)

        positive_to_macro = macro.positive_cells.copy()
        macro_to_positive = np.full(n_x, -1, dtype=int)
        for p, j in enumerate(positive_to_macro):
            macro_to_positive[j] = p

        return cls(
            n_x=n_x,
            n_pos=n_pos,
            n_r=n_r,
            c_e=c_e,
            phi_e=phi_e,
            phi_s=phi_s,
            c_s=c_s,
            positive_to_macro=positive_to_macro,
            macro_to_positive=macro_to_positive,
        )


@dataclass
class P2DState:
    """P2D 模型的物理状态。

    Attributes
    ----------
    c_e : ndarray, shape (n_x,)
        电解质 Li+ 浓度 [mol/m³]。
    phi_e : ndarray, shape (n_x,)
        电解质电位 [V]。
    phi_s : ndarray, shape (n_pos,)
        固相电位 [V] (仅 positive electrode)。
    c_s : ndarray, shape (n_pos, n_r)
        固相锂浓度 [mol/m³] (仅 positive electrode)。
    time : float
        当前时间 [s]。
    """

    c_e: Array
    phi_e: Array
    phi_s: Array
    c_s: Array
    time: float = 0.0

    def copy(self) -> "P2DState":
        """返回深拷贝。"""
        return P2DState(
            c_e=self.c_e.copy(),
            phi_e=self.phi_e.copy(),
            phi_s=self.phi_s.copy(),
            c_s=self.c_s.copy(),
            time=self.time,
        )

    def surface_concentration(self) -> Array:
        """返回各 positive cell 的表面浓度 (最外层壳层) [mol/m³]。

        Returns
        -------
        c_s_surf : ndarray, shape (n_pos,)
        """
        return self.c_s[:, -1].copy()

    def average_solid_concentration(self, particle: ParticleMesh) -> Array:
        """计算各 positive cell 的体积加权平均固相浓度。

        对应 TASK_PHASE0.md 公式 2.2: c̄_s = 3/R_p³ ∫ c_s r² dr

        Parameters
        ----------
        particle : ParticleMesh
            颗粒径向网格。

        Returns
        -------
        c_s_avg : ndarray, shape (n_pos,)
        """
        # 壳层体积加权平均: c̄_s = sum(c_s_i * V_i) / sum(V_i)
        # sum(V_i) = (4/3)*pi*R_p^3
        total_volume = (4.0 / 3.0) * np.pi * particle.radius ** 3
        # c_s shape: (n_pos, n_r), shell_volumes shape: (n_r,)
        return np.dot(self.c_s, particle.shell_volumes) / total_volume


def pack_state(state: P2DState, layout: StateLayout) -> Array:
    """将物理状态打包为 Newton 向量。

    Parameters
    ----------
    state : P2DState
        物理状态。
    layout : StateLayout
        状态布局。

    Returns
    -------
    y : ndarray
        Newton 向量。
    """
    y = np.empty(layout.size)
    y[layout.c_e] = state.c_e
    y[layout.phi_e] = state.phi_e
    y[layout.phi_s] = state.phi_s
    y[layout.c_s] = state.c_s.ravel()
    return y


def unpack_state(y: Array, layout: StateLayout) -> P2DState:
    """将 Newton 向量解包为物理状态。

    Parameters
    ----------
    y : ndarray
        Newton 向量。
    layout : StateLayout
        状态布局。

    Returns
    -------
    P2DState
    """
    return P2DState(
        c_e=y[layout.c_e].copy(),
        phi_e=y[layout.phi_e].copy(),
        phi_s=y[layout.phi_s].copy(),
        c_s=y[layout.c_s].reshape(layout.n_pos, layout.n_r).copy(),
    )


def make_initial_state(
    macro: MacroMesh,
    particle: ParticleMesh,
    c_e0: float,
    soc0: float,
    c_s_max: float,
    phi_e0: float = 0.0,
) -> P2DState:
    """构建初始状态。

    Parameters
    ----------
    macro : MacroMesh
        宏观网格。
    particle : ParticleMesh
        颗粒网格。
    c_e0 : float
        初始电解质浓度 [mol/m³]。
    soc0 : float
        初始 SOC (0~1)。
    c_s_max : float
        固相最大浓度 [mol/m³]。
    phi_e0 : float
        初始电解质电位 [V]，默认 0.0。

    Returns
    -------
    P2DState
    """
    n_x = macro.n_x
    n_pos = macro.n_pos
    n_r = particle.n_shells

    c_s0 = soc0 * c_s_max
    # 初始 OCV 作为 phi_s 初始值
    from pnmcathode.physics.ocv import nmc532_ocv

    phi_s0 = float(nmc532_ocv(soc0))

    return P2DState(
        c_e=np.full(n_x, c_e0),
        phi_e=np.full(n_x, phi_e0),
        phi_s=np.full(n_pos, phi_s0),
        c_s=np.full((n_pos, n_r), c_s0),
        time=0.0,
    )
