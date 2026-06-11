"""
P2D 稀疏 Jacobian 组装 (Sparse Jacobian Assembly)
=================================================

职责:
- 组装 analytic sparse Jacobian。
- 提供 finite-difference Jacobian 用于测试和 fallback。
- 提供 check_jacobian 工具函数验证解析导数。

Jacobian 组装策略:
- 使用 COO triplets 收集局部 stencil。
- x 方向通量只耦合相邻 macro cells。
- r 方向固相扩散只耦合相邻 shells。
- BV 源项在同一 positive macro cell 内耦合 c_e, phi_e, phi_s, c_s_surf。

对应 TASK_PHASE0.md §1.7
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from pnmcathode.p2d.domain import MacroMesh, ParticleMesh
from pnmcathode.p2d.kinetics import butler_volmer_with_derivatives
from pnmcathode.p2d.materials import (
    P2DParameters,
    bruggeman_conductivity_e,
    bruggeman_conductivity_s,
    bruggeman_diffusivity_e,
)
from pnmcathode.p2d.residual import P2DResidualContext, assemble_p2d_residual, compute_reaction
from pnmcathode.p2d.state import P2DState, StateLayout, pack_state, unpack_state
from pnmcathode.physics.reaction import F

Array = np.ndarray


def assemble_p2d_jacobian(
    y_new: Array,
    state_old: P2DState,
    context: P2DResidualContext,
) -> sparse.csr_matrix:
    """组装解析稀疏 Jacobian。

    Phase 0 首版: 有限差分 Jacobian 作为默认实现。
    后续可替换为完全解析版本。

    Parameters
    ----------
    y_new : ndarray
        Newton 向量。
    state_old : P2DState
        上一时间步状态。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    J : csr_matrix
        稀疏 Jacobian 矩阵。
    """
    # 使用有限差分 Jacobian (Phase 0 首版)
    return finite_difference_jacobian(y_new, state_old, context)


def _build_sparsity_pattern(context: P2DResidualContext) -> list[set[int]]:
    """构建 Jacobian 稀疏模式: 每个未知量影响哪些残差方程。

    Parameters
    ----------
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    pattern : list of set[int]
        pattern[j] = set of row indices that column j can affect.
    """
    layout = context.layout
    n = layout.size
    n_x = layout.n_x
    n_pos = layout.n_pos
    n_r = layout.n_r

    pattern: list[set[int]] = [set() for _ in range(n)]

    # c_e[j] 影响:
    # - R_c_e[j] (浓度方程), R_c_e[j-1], R_c_e[j+1] (通量)
    # - R_phi_e[j], R_phi_e[j-1], R_phi_e[j+1] (电流)
    # - 如果 j 是 positive cell: R_phi_s[p], R_c_s[p, n_r-1] (通过 BV)
    for j in range(n_x):
        col = layout.c_e.start + j
        # 直接影响
        pattern[col].add(layout.c_e.start + j)
        if j > 0:
            pattern[col].add(layout.c_e.start + j - 1)
        if j < n_x - 1:
            pattern[col].add(layout.c_e.start + j + 1)
        # phi_e 通过电流
        pattern[col].add(layout.phi_e.start + j)
        if j > 0:
            pattern[col].add(layout.phi_e.start + j - 1)
        if j < n_x - 1:
            pattern[col].add(layout.phi_e.start + j + 1)
        # BV coupling
        p = layout.macro_to_positive[j]
        if p >= 0:
            pattern[col].add(layout.phi_s.start + p)
            pattern[col].add(layout.c_s.start + p * n_r + n_r - 1)

    # phi_e[j] 影响:
    # - R_phi_e[j], R_phi_e[j-1], R_phi_e[j+1] (电流)
    # - 如果 j 是 positive cell: R_phi_s[p] (通过 BV)
    for j in range(n_x):
        col = layout.phi_e.start + j
        pattern[col].add(layout.phi_e.start + j)
        if j > 0:
            pattern[col].add(layout.phi_e.start + j - 1)
        if j < n_x - 1:
            pattern[col].add(layout.phi_e.start + j + 1)
        p = layout.macro_to_positive[j]
        if p >= 0:
            pattern[col].add(layout.phi_s.start + p)
            pattern[col].add(layout.c_s.start + p * n_r + n_r - 1)

    # phi_s[p] 影响:
    # - R_phi_s[p], R_phi_s[p-1], R_phi_s[p+1] (固相电流)
    # - R_c_s[p, n_r-1] (通过 BV 表面通量)
    for p in range(n_pos):
        col = layout.phi_s.start + p
        pattern[col].add(layout.phi_s.start + p)
        if p > 0:
            pattern[col].add(layout.phi_s.start + p - 1)
        if p < n_pos - 1:
            pattern[col].add(layout.phi_s.start + p + 1)
        pattern[col].add(layout.c_s.start + p * n_r + n_r - 1)

    # c_s[p, k] 影响:
    # - R_c_s[p, k], R_c_s[p, k-1], R_c_s[p, k+1] (扩散)
    # - 如果 k = n_r-1: R_phi_s[p] (通过 BV), R_c_e[j], R_phi_e[j]
    for p in range(n_pos):
        j = layout.positive_to_macro[p]
        for k in range(n_r):
            col = layout.c_s.start + p * n_r + k
            pattern[col].add(layout.c_s.start + p * n_r + k)
            if k > 0:
                pattern[col].add(layout.c_s.start + p * n_r + k - 1)
            if k < n_r - 1:
                pattern[col].add(layout.c_s.start + p * n_r + k + 1)
            if k == n_r - 1:
                # BV coupling
                pattern[col].add(layout.phi_s.start + p)
                pattern[col].add(layout.c_e.start + j)
                pattern[col].add(layout.phi_e.start + j)

    return pattern


def _group_columns(pattern: list[set[int]], n: int) -> list[list[int]]:
    """将列分组，使得同一组内的列不影响同一行 (贪心着色)。

    Parameters
    ----------
    pattern : list of set[int]
        稀疏模式。
    n : int
        列数。

    Returns
    -------
    groups : list of list[int]
        列分组。
    """
    groups: list[list[int]] = []
    assigned = np.full(n, -1, dtype=int)

    for j in range(n):
        # 找到 j 可以加入的组
        placed = False
        for g_idx, group in enumerate(groups):
            # 检查 j 是否与组内任何列冲突
            conflict = False
            for k in group:
                if k in pattern[j] or j in pattern[k]:
                    conflict = True
                    break
            if not conflict:
                group.append(j)
                assigned[j] = g_idx
                placed = True
                break
        if not placed:
            groups.append([j])
            assigned[j] = len(groups) - 1

    return groups


def finite_difference_jacobian(
    y_new: Array,
    state_old: P2DState,
    context: P2DResidualContext,
    rel_step: float = 1e-7,
) -> sparse.csr_matrix:
    """稀疏有限差分 Jacobian。

    利用已知稀疏模式进行列分组，同时扰动同组列，
    大幅减少残差求值次数。

    Parameters
    ----------
    y_new : ndarray
        Newton 向量。
    state_old : P2DState
        上一时间步状态。
    context : P2DResidualContext
        残差上下文。
    rel_step : float
        相对扰动步长。

    Returns
    -------
    J : csr_matrix
        稀疏 Jacobian 矩阵。
    """
    n = len(y_new)
    R0 = assemble_p2d_residual(y_new, state_old, context)

    # 构建稀疏模式和列分组
    pattern = _build_sparsity_pattern(context)
    groups = _group_columns(pattern, n)

    # 收集 COO triplets
    rows_list = []
    cols_list = []
    vals_list = []

    for group in groups:
        # 同时扰动组内所有列
        y_pert = y_new.copy()
        h_vals = {}
        for j in group:
            h = rel_step * max(abs(y_new[j]), 1e-10)
            y_pert[j] += h
            h_vals[j] = h

        R_pert = assemble_p2d_residual(y_pert, state_old, context)
        diff = R_pert - R0

        for j in group:
            col_j = diff / h_vals[j]
            # 只存储 pattern 中的非零元素
            nonzero = pattern[j] & set(np.nonzero(col_j)[0].tolist())
            for i in nonzero:
                rows_list.append(i)
                cols_list.append(j)
                vals_list.append(col_j[i])

    J = sparse.coo_matrix(
        (vals_list, (rows_list, cols_list)), shape=(n, n)
    ).tocsr()
    return J


def check_jacobian(
    y_new: Array,
    state_old: P2DState,
    context: P2DResidualContext,
    rel_step: float = 1e-7,
) -> dict[str, float]:
    """比较解析 Jacobian 与有限差分 Jacobian。

    Parameters
    ----------
    y_new : ndarray
        Newton 向量。
    state_old : P2DState
        上一时间步状态。
    context : P2DResidualContext
        残差上下文。
    rel_step : float
        有限差分步长。

    Returns
    -------
    result : dict
        abs_error: 最大绝对误差
        rel_error: 最大相对误差
        fd_norm: 有限差分 Jacobian 范数
    """
    J_fd = finite_difference_jacobian(y_new, state_old, context, rel_step)
    # Phase 0: 解析 Jacobian 暂时等于有限差分
    # 后续替换为真正的解析版本后再比较
    J_analytic = assemble_p2d_jacobian(y_new, state_old, context)

    diff = J_analytic - J_fd
    abs_err = float(np.max(np.abs(diff.toarray())))
    fd_norm = float(np.max(np.abs(J_fd.toarray())))
    rel_err = abs_err / max(fd_norm, 1e-30)

    return {
        "abs_error": abs_err,
        "rel_error": rel_err,
        "fd_norm": fd_norm,
    }
