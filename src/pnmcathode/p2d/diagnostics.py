"""
P2D 诊断与守恒检查 (Diagnostics & Conservation Checks)
======================================================

职责:
- Phase 0 验证与 solver runtime diagnostics。
- 守恒检查: 质量、电荷、总锂量。
- 平衡态残差检查。

对应 TASK_PHASE0.md §1.10, §5.1
"""

from __future__ import annotations

import numpy as np

from pnmcathode.p2d.domain import ParticleMesh
from pnmcathode.p2d.materials import P2DParameters
from pnmcathode.p2d.residual import P2DFluxes, P2DResidualContext, compute_fluxes
from pnmcathode.p2d.state import P2DState, StateLayout
from pnmcathode.physics.ocv import nmc532_ocv
from pnmcathode.physics.reaction import F

Array = np.ndarray


def total_lithium_inventory(
    state: P2DState,
    context: P2DResidualContext,
) -> float:
    """计算系统中总锂量 [mol]。

    电解质: ∑ V_j * ε_e,j * c_e,j
    固相:   ∑ V_j * ε_s,j * c̄_s,j

    Parameters
    ----------
    state : P2DState
        当前状态。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    total : float
        总锂量 [mol]。
    """
    macro = context.macro
    layout = context.layout
    params = context.params
    particle = context.particle
    pos = params.positive

    # 电解质锂
    Li_e = 0.0
    for j in range(layout.n_x):
        reg = macro.region_for_cell(j)
        Li_e += macro.volumes[j] * reg.epsilon_e * state.c_e[j]

    # 固相锂
    Li_s = 0.0
    c_s_avg = state.average_solid_concentration(particle)
    for p in range(layout.n_pos):
        j = layout.positive_to_macro[p]
        Li_s += macro.volumes[j] * pos.epsilon_s * c_s_avg[p]

    return Li_e + Li_s


def current_conservation_error(
    fluxes: P2DFluxes,
    context: P2DResidualContext,
) -> float:
    """检查总电流守恒: max| i_e + i_s - I_app |。

    对应 TASK_PHASE0.md 公式 2.4

    Parameters
    ----------
    fluxes : P2DFluxes
        通量信息。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    error : float
        总电流守恒最大绝对误差。
    """
    I_app = context.current_density
    i_total = fluxes.ionic_current_faces + fluxes.solid_current_faces
    return float(np.max(np.abs(i_total - I_app)))


def mass_conservation_error(
    state_new: P2DState,
    state_old: P2DState,
    fluxes: P2DFluxes,
    context: P2DResidualContext,
) -> dict[str, float]:
    """检查质量守恒。

    电解质: ∑ V_j * ε_e * (c_e_new - c_e_old) / dt = ∑ 边界通量 + ∑ 源项
    固相:   V_j * ε_s * (c̄_s_new - c̄_s_old) / dt = -V_j * a_s * i_F / F

    Parameters
    ----------
    state_new : P2DState
        新状态。
    state_old : P2DState
        旧状态。
    fluxes : P2DFluxes
        通量信息。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    errors : dict
        electrolyte: 电解质质量守恒相对误差
        solid: 固相质量守恒最大绝对误差
    """
    macro = context.macro
    layout = context.layout
    params = context.params
    particle = context.particle
    pos = params.positive
    dt = context.dt
    t_plus = params.material.separator.t_plus

    # 电解质质量守恒
    delta_Li_e = 0.0
    for j in range(layout.n_x):
        reg = macro.region_for_cell(j)
        delta_Li_e += macro.volumes[j] * reg.epsilon_e * (state_new.c_e[j] - state_old.c_e[j]) / dt

    # 边界通量 (盐通量)
    boundary_flux = 0.0
    # 无通量边界 → boundary_flux = 0

    # 反应源项
    source_total = 0.0
    a_s = pos.area_density()
    for p in range(layout.n_pos):
        j = layout.positive_to_macro[p]
        source_total += macro.volumes[j] * (1.0 - t_plus) * a_s * fluxes.faradaic_current[j] / F

    expected = boundary_flux + source_total
    e_error = abs(delta_Li_e - expected) / max(abs(expected), 1e-30)

    # 固相质量守恒
    c_s_avg_new = state_new.average_solid_concentration(particle)
    c_s_avg_old = state_old.average_solid_concentration(particle)
    max_s_error = 0.0
    for p in range(layout.n_pos):
        j = layout.positive_to_macro[p]
        delta = macro.volumes[j] * pos.epsilon_s * (c_s_avg_new[p] - c_s_avg_old[p]) / dt
        source = -macro.volumes[j] * a_s * fluxes.faradaic_current[j] / F
        max_s_error = max(max_s_error, abs(delta - source))

    return {
        "electrolyte": float(e_error),
        "solid": float(max_s_error),
    }


def equilibrium_error(
    state: P2DState,
    context: P2DResidualContext,
) -> dict[str, float]:
    """检查平衡态残差。

    平衡态条件 (I_app = 0):
    - η = 0 → phi_s - phi_e = U(c_s_surf)
    - i_F = 0
    - c_e 均匀
    - c_s 均匀

    Parameters
    ----------
    state : P2DState
        当前状态。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    errors : dict
        eta_max: 最大过电位绝对值 [V]
        i_f_max: 最大 Faradaic 电流 [A/m²]
        ce_variation: c_e 相对变化
        cs_variation: c_s 相对变化
    """
    layout = context.layout
    params = context.params
    mat = params.material
    cs_max = mat.active.cs_max

    c_e_pos = state.c_e[layout.positive_to_macro]
    phi_e_pos = state.phi_e[layout.positive_to_macro]
    c_s_surf = state.surface_concentration()

    soc = np.clip(c_s_surf / cs_max, 1e-9, 1.0 - 1e-9)
    ocv = mat.ocv_model.value(soc)

    eta = state.phi_s - phi_e_pos - ocv
    i0 = np.array([
        mat.kinetics.k0 * F * max(c_e_pos[p], 1e-10) ** mat.kinetics.alpha_a
        * max(cs_max - c_s_surf[p], 1e-10) ** mat.kinetics.alpha_a
        * max(c_s_surf[p], 1e-10) ** mat.kinetics.alpha_c
        for p in range(layout.n_pos)
    ])

    # BV 电流
    f = F / (8.314462 * mat.temperature)
    arg_a = np.clip(mat.kinetics.alpha_a * f * eta, -500.0, 500.0)
    arg_c = np.clip(-mat.kinetics.alpha_c * f * eta, -500.0, 500.0)
    i_f = i0 * (np.exp(arg_a) - np.exp(arg_c))

    ce_range = float(np.max(c_e_pos) - np.min(c_e_pos)) / max(float(np.mean(c_e_pos)), 1e-10)
    cs_range = float(np.max(c_s_surf) - np.min(c_s_surf)) / max(float(np.mean(c_s_surf)), 1e-10)

    return {
        "eta_max": float(np.max(np.abs(eta))),
        "i_f_max": float(np.max(np.abs(i_f))),
        "ce_variation": ce_range,
        "cs_variation": cs_range,
    }
