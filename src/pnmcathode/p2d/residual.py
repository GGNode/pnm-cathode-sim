"""
P2D 单调残差组装 (Monolithic Residual Assembly)
===============================================

职责:
- 组装完整的 nonlinear residual (c_e, phi_e, phi_s, c_s)。
- 所有传输项使用有限体积通量差分。
- 实现 half-cell 边界条件与 potential gauge。
- 提供诊断通量信息。

符号约定:
- 阳极约定: i_F > 0 → 脱锂, i_F < 0 → 嵌锂 (放电)
- 放电: I_app < 0, i_F < 0, c_s 增加
- c_e residual: mol/m³/s (per-volume)
- phi_e, phi_s residual: A/m³ (per-volume)
- c_s residual: mol/m³/s (per-volume)

对应 TASK_PHASE0.md §1.6, §2
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from pnmcathode.p2d.domain import MacroMesh, ParticleMesh
from pnmcathode.p2d.kinetics import (
    ReactionRates,
    butler_volmer_with_derivatives,
    exchange_current_density_vec,
)
from pnmcathode.p2d.materials import (
    P2DParameters,
    bruggeman_conductivity_e,
    bruggeman_conductivity_s,
    bruggeman_diffusivity_e,
)
from pnmcathode.p2d.state import P2DState, StateLayout, pack_state, unpack_state
from pnmcathode.physics.ocv import nmc532_ocv, ocv_derivative
from pnmcathode.physics.reaction import F

Array = np.ndarray


@dataclass(frozen=True)
class P2DResidualContext:
    """残差组装所需的全部上下文信息。

    Attributes
    ----------
    macro : MacroMesh
        宏观有限体积网格。
    particle : ParticleMesh
        颗粒径向网格。
    layout : StateLayout
        状态向量布局。
    params : P2DParameters
        模型参数。
    dt : float
        时间步长 [s]。
    current_density : float
        应用电流密度 [A/m²] (阳极约定, 放电 < 0)。
    """

    macro: MacroMesh
    particle: ParticleMesh
    layout: StateLayout
    params: P2DParameters
    dt: float
    current_density: float


@dataclass(frozen=True)
class P2DFluxes:
    """诊断通量信息。

    Attributes
    ----------
    salt_flux_faces : ndarray
        电解质盐通量 [mol/(m²·s)]，shape (n_x+1,)。
    ionic_current_faces : ndarray
        电解质电流密度 [A/m²]，shape (n_x+1,)。
    solid_current_faces : ndarray
        固相电流密度 [A/m²]，shape (n_x+1,)。
    faradaic_current : ndarray
        Faradaic 电流密度 [A/m²]，shape (n_x,)。
    volumetric_reaction : ndarray
        体积反应速率 [mol/(m³·s)]，shape (n_x,)。
    """

    salt_flux_faces: Array
    ionic_current_faces: Array
    solid_current_faces: Array
    faradaic_current: Array
    volumetric_reaction: Array


def compute_reaction(
    state: P2DState,
    context: P2DResidualContext,
) -> ReactionRates:
    """计算各 positive cell 的 Butler-Volmer 反应速率及导数。

    Parameters
    ----------
    state : P2DState
        当前物理状态。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    ReactionRates
    """
    layout = context.layout
    params = context.params
    mat = params.material
    pos = params.positive
    n_pos = layout.n_pos

    # 提取 positive cell 的局部状态
    c_e_pos = state.c_e[layout.positive_to_macro]
    phi_e_pos = state.phi_e[layout.positive_to_macro]
    phi_s = state.phi_s
    c_s_surf = state.surface_concentration()

    # OCV 及其导数
    cs_max = mat.active.cs_max
    soc = np.clip(c_s_surf / cs_max, 1e-9, 1.0 - 1e-9)
    ocv = mat.ocv_model.value(soc)
    docv_dsoc = mat.ocv_model.derivative(soc)
    # dU/d(c_s) = dU/d(soc) / c_s_max
    docv_dc_s = docv_dsoc / cs_max

    return butler_volmer_with_derivatives(
        c_e=c_e_pos,
        c_s_surf=c_s_surf,
        phi_e=phi_e_pos,
        phi_s=phi_s,
        c_s_max=cs_max,
        ocv=ocv,
        docv_dc_s=docv_dc_s,
        kinetics=mat.kinetics,
        temperature=mat.temperature,
    )


def compute_fluxes(
    state: P2DState,
    context: P2DResidualContext,
) -> P2DFluxes:
    """计算所有通量和反应源项。

    Parameters
    ----------
    state : P2DState
        当前物理状态。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    P2DFluxes
    """
    macro = context.macro
    layout = context.layout
    params = context.params
    mat = params.material
    pos = params.positive
    sep = params.separator_region
    n_x = layout.n_x

    # 反应速率
    rates = compute_reaction(state, context)

    # 将 i_F 映射到 macro 网格
    i_f_macro = np.zeros(n_x)
    i_f_macro[layout.positive_to_macro] = rates.i_f

    # 比表面积
    a_s = pos.area_density()
    a_s_macro = np.zeros(n_x)
    a_s_macro[layout.positive_to_macro] = a_s

    # ===== 电解质物性 (cell-centered) =====
    c_e = state.c_e
    D_e_eff = np.empty(n_x)
    kappa_eff = np.empty(n_x)
    for j in range(n_x):
        reg = macro.region_for_cell(j)
        D_e_val = bruggeman_diffusivity_e(np.atleast_1d(c_e[j]), reg, mat)
        D_e_eff[j] = float(np.atleast_1d(D_e_val)[0])
        kappa_val = bruggeman_conductivity_e(np.atleast_1d(c_e[j]), reg, mat)
        kappa_eff[j] = float(np.atleast_1d(kappa_val)[0])

    # ===== 电解质扩散通量 (face-centered) =====
    salt_flux = np.zeros(n_x + 1)
    for j in range(1, n_x):
        D_face = 2.0 * D_e_eff[j - 1] * D_e_eff[j] / (D_e_eff[j - 1] + D_e_eff[j] + 1e-30)
        dx_face = macro.x_centers[j] - macro.x_centers[j - 1]
        salt_flux[j] = D_face * (c_e[j] - c_e[j - 1]) / dx_face

    # ===== 电解质电流 (face-centered) =====
    # Half-cell 左边界: Li metal reference → i_e[0] = I_app
    # 右边界 (current collector): i_e[n_x] = 0 (全部电流走固相)
    i_e = np.zeros(n_x + 1)
    i_e[0] = context.current_density  # Li metal/separator 左边界
    t_plus = mat.separator.t_plus
    c_floor = params.concentration_floor
    for j in range(1, n_x):
        k_face = 2.0 * kappa_eff[j - 1] * kappa_eff[j] / (kappa_eff[j - 1] + kappa_eff[j] + 1e-30)
        dx_face = macro.x_centers[j] - macro.x_centers[j - 1]
        # Ohmic
        ohmic = k_face * (state.phi_e[j] - state.phi_e[j - 1]) / dx_face
        # 浓差
        ce_avg = 0.5 * (c_e[j] + c_e[j - 1])
        ce_safe = max(ce_avg, c_floor)
        ln_ce = np.log(ce_safe)
        ce_j = max(c_e[j], c_floor)
        ce_jm1 = max(c_e[j - 1], c_floor)
        dln = np.log(ce_j) - np.log(ce_jm1)
        T = mat.temperature
        conc = 2.0 * 8.314462 * T * k_face / F * (1.0 - t_plus) * dln / dx_face
        i_e[j] = -ohmic + conc

    # ===== 固相电流 (face-centered) =====
    sigma_eff = bruggeman_conductivity_s(pos, mat)
    i_s = np.zeros(n_x + 1)
    for j in range(1, n_x):
        if layout.macro_to_positive[j - 1] >= 0 and layout.macro_to_positive[j] >= 0:
            dx_face = macro.x_centers[j] - macro.x_centers[j - 1]
            i_s[j] = -sigma_eff * (state.phi_s[layout.macro_to_positive[j]] - state.phi_s[layout.macro_to_positive[j - 1]]) / dx_face
    # 右边界: i_s[n_x] = I_app
    i_s[n_x] = context.current_density

    return P2DFluxes(
        salt_flux_faces=salt_flux,
        ionic_current_faces=i_e,
        solid_current_faces=i_s,
        faradaic_current=i_f_macro,
        volumetric_reaction=a_s_macro * i_f_macro / F,
    )


def assemble_p2d_residual(
    y_new: Array,
    state_old: P2DState,
    context: P2DResidualContext,
) -> Array:
    """组装 P2D 单调残差向量。

    残差方程:
        R_c_e:  ε_e * (c_e_new - c_e_old) / dt + (N[j] - N[j+1]) / dx - S = 0
        R_phi_e: (i_e[j+1] - i_e[j]) / dx + a_s * i_F = 0
        R_phi_s: (i_s[j+1] - i_s[j]) / dx - a_s * i_F = 0
        R_c_s:  (c_s_new - c_s_old) / dt + div(J_r) = 0

    Parameters
    ----------
    y_new : ndarray
        Newton 向量。
    state_old : P2DState
        上一时间步的状态。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    residual : ndarray
        残差向量，shape (layout.size,)。
    """
    layout = context.layout
    params = context.params
    mat = params.material
    pos = params.positive
    sep = params.separator_region
    macro = context.macro
    particle = context.particle
    dt = context.dt

    state_new = unpack_state(y_new, layout)
    n_x = layout.n_x
    n_pos = layout.n_pos
    n_r = layout.n_r

    residual = np.zeros(layout.size)

    # ===== 反应速率 =====
    rates = compute_reaction(state_new, context)
    i_f_macro = np.zeros(n_x)
    i_f_macro[layout.positive_to_macro] = rates.i_f

    a_s = pos.area_density()
    t_plus = mat.separator.t_plus
    cs_max = mat.active.cs_max
    c_floor = params.concentration_floor

    # ===== 电解质物性 (cell-centered) =====
    c_e = state_new.c_e
    D_e_eff = np.empty(n_x)
    kappa_eff = np.empty(n_x)
    for j in range(n_x):
        reg = macro.region_for_cell(j)
        D_e_val = bruggeman_diffusivity_e(np.atleast_1d(c_e[j]), reg, mat)
        D_e_eff[j] = float(np.atleast_1d(D_e_val)[0])
        kappa_val = bruggeman_conductivity_e(np.atleast_1d(c_e[j]), reg, mat)
        kappa_eff[j] = float(np.atleast_1d(kappa_val)[0])

    sigma_eff = bruggeman_conductivity_s(pos, mat)

    # ===== 1. 电解质浓度残差 =====
    # R = ε_e * (c_e_new - c_e_old)/dt + (N[j] - N[j+1])/dx - (1-t+)*a_s*i_F/F
    for j in range(n_x):
        reg = macro.region_for_cell(j)
        eps_e = reg.epsilon_e
        dx_j = macro.dx[j]

        # 盐通量
        N_left = 0.0  # 边界无通量
        N_right = 0.0
        if j > 0:
            D_face = 2.0 * D_e_eff[j - 1] * D_e_eff[j] / (D_e_eff[j - 1] + D_e_eff[j] + 1e-30)
            dx_face = macro.x_centers[j] - macro.x_centers[j - 1]
            N_left = D_face * (c_e[j] - c_e[j - 1]) / dx_face
        if j < n_x - 1:
            D_face = 2.0 * D_e_eff[j] * D_e_eff[j + 1] / (D_e_eff[j] + D_e_eff[j + 1] + 1e-30)
            dx_face = macro.x_centers[j + 1] - macro.x_centers[j]
            N_right = D_face * (c_e[j + 1] - c_e[j]) / dx_face

        # 源项 (仅 positive cells)
        source = 0.0
        if layout.macro_to_positive[j] >= 0:
            source = (1.0 - t_plus) * a_s * i_f_macro[j] / F

        residual[layout.c_e][j] = (
            eps_e * (state_new.c_e[j] - state_old.c_e[j]) / dt
            + (N_left - N_right) / dx_j
            - source
        )

    # ===== 2. 电解质电位残差 =====
    # R = (i_e[j+1] - i_e[j]) / dx + a_s * i_F
    # Half-cell 左边界: i_e[0] = I_app (Li metal reference)
    for j in range(n_x):
        reg = macro.region_for_cell(j)
        dx_j = macro.dx[j]

        # 电流
        if j == 0:
            i_e_left = context.current_density  # Li metal 左边界
        else:
            i_e_left = 0.0
        i_e_right = 0.0
        if j > 0:
            k_face = 2.0 * kappa_eff[j - 1] * kappa_eff[j] / (kappa_eff[j - 1] + kappa_eff[j] + 1e-30)
            dx_face = macro.x_centers[j] - macro.x_centers[j - 1]
            ohmic = k_face * (state_new.phi_e[j] - state_new.phi_e[j - 1]) / dx_face
            ce_avg = 0.5 * (c_e[j] + c_e[j - 1])
            ce_safe = max(ce_avg, c_floor)
            ce_j = max(c_e[j], c_floor)
            ce_jm1 = max(c_e[j - 1], c_floor)
            dln = np.log(ce_j) - np.log(ce_jm1)
            T = mat.temperature
            conc = 2.0 * 8.314462 * T * k_face / F * (1.0 - t_plus) * dln / dx_face
            i_e_left = -ohmic + conc
        if j < n_x - 1:
            k_face = 2.0 * kappa_eff[j] * kappa_eff[j + 1] / (kappa_eff[j] + kappa_eff[j + 1] + 1e-30)
            dx_face = macro.x_centers[j + 1] - macro.x_centers[j]
            ohmic = k_face * (state_new.phi_e[j + 1] - state_new.phi_e[j]) / dx_face
            ce_avg = 0.5 * (c_e[j + 1] + c_e[j])
            ce_safe = max(ce_avg, c_floor)
            ce_jp1 = max(c_e[j + 1], c_floor)
            ce_j = max(c_e[j], c_floor)
            dln = np.log(ce_jp1) - np.log(ce_j)
            T = mat.temperature
            conc = 2.0 * 8.314462 * T * k_face / F * (1.0 - t_plus) * dln / dx_face
            i_e_right = -ohmic + conc

        source = 0.0
        if layout.macro_to_positive[j] >= 0:
            source = a_s * i_f_macro[j]

        residual[layout.phi_e][j] = (i_e_right - i_e_left) / dx_j + source

    # Potential gauge: phi_e(0) = 0 (Dirichlet)
    residual[layout.phi_e][0] = state_new.phi_e[0]

    # ===== 3. 固相电位残差 =====
    # R = (i_s_right - i_s_left) / dx - a_s * i_F
    for p in range(n_pos):
        j = layout.positive_to_macro[p]
        dx_j = macro.dx[j]

        # i_s left
        i_s_left = 0.0
        if p > 0:
            j_left = layout.positive_to_macro[p - 1]
            dx_face = macro.x_centers[j] - macro.x_centers[j_left]
            i_s_left = -sigma_eff * (state_new.phi_s[p] - state_new.phi_s[p - 1]) / dx_face
        # else: separator interface, i_s = 0

        # i_s right
        i_s_right = 0.0
        if p < n_pos - 1:
            j_right = layout.positive_to_macro[p + 1]
            dx_face = macro.x_centers[j_right] - macro.x_centers[j]
            i_s_right = -sigma_eff * (state_new.phi_s[p + 1] - state_new.phi_s[p]) / dx_face
        else:
            # 右边界: i_s = I_app
            i_s_right = context.current_density

        residual[layout.phi_s][p] = (i_s_right - i_s_left) / dx_j - a_s * i_f_macro[j]

    # ===== 4. 固相颗粒扩散残差 =====
    r_faces = particle.r_faces
    shell_volumes = particle.shell_volumes
    face_areas = particle.face_areas

    for p in range(n_pos):
        cs_old = state_old.c_s[p, :]
        cs_new = state_new.c_s[p, :]
        j = layout.positive_to_macro[p]
        reg = macro.region_for_cell(j)
        eps_s = reg.epsilon_s

        # 扩散系数 (使用 shell 平均值)
        cs_avg = 0.5 * (cs_new[:-1] + cs_new[1:])
        D_s_faces = np.array([
            mat.active.diffusivity(float(c), mat.temperature) for c in cs_avg
        ])

        # 面通量: J = -D_s * (c_s[k+1] - c_s[k]) / dr
        dr = particle.radius / n_r
        J = np.zeros(n_r + 1)
        J[0] = 0.0  # 中心对称
        for k in range(1, n_r):
            J[k] = -D_s_faces[k - 1] * (cs_new[k] - cs_new[k - 1]) / dr

        # 表面通量: J[R_p] = i_F / F
        # 阳极约定下放电 i_F < 0，负的外向通量表示 Li 进入颗粒。
        J[n_r] = i_f_macro[j] / F

        # 通量散度
        flux_div = np.zeros(n_r)
        for k in range(n_r):
            flux_in = J[k] * face_areas[k]
            flux_out = J[k + 1] * face_areas[k + 1]
            flux_div[k] = (flux_out - flux_in) / shell_volumes[k]

        offset = layout.c_s.start + p * n_r
        residual[offset:offset + n_r] = (
            (cs_new - cs_old) / dt + flux_div
        )

    return residual


def residual_norm(
    residual: Array,
    context: P2DResidualContext,
) -> float:
    """计算归一化残差范数。

    按残差方程的物理尺度归一化:
    - c_e: 除以 c_e_ref / dt
    - phi_e, phi_s: 除以特征体积电流密度 [A/m³]
    - c_s: 除以 c_s_max / dt

    Parameters
    ----------
    residual : ndarray
        残差向量。
    context : P2DResidualContext
        残差上下文。

    Returns
    -------
    norm : float
        归一化残差的 RMS 范数。
    """
    layout = context.layout
    mat = context.params.material
    params = context.params

    c_e_ref = 1000.0  # mol/m³
    c_s_max = mat.active.cs_max
    dt = max(context.dt, 1e-30)

    # 电荷守恒残差单位为 A/m³，不能用电压尺度归一化。
    # 使用外加电流在最小控制体长度上的尺度，并用平衡态 BV
    # 交换电流体积源作为低电流工况下的数值尺度。
    dx_min = max(float(np.min(context.macro.dx)), 1e-30)
    applied_scale = abs(context.current_density) / dx_min
    c_e0 = mat.electrolyte.c_init
    c_s0 = 0.5 * c_s_max
    i0_ref = float(exchange_current_density_vec(
        np.array([c_e0]),
        np.array([c_s0]),
        c_s_max,
        mat.kinetics,
    ).i0[0])
    reaction_scale = params.positive.area_density() * i0_ref
    current_scale = max(applied_scale, reaction_scale, 1.0)

    scaled = residual.copy()
    scaled[layout.c_e] /= c_e_ref / dt
    scaled[layout.phi_e] /= current_scale
    scaled[layout.phi_s] /= current_scale
    scaled[layout.c_s] /= c_s_max / dt

    return float(np.sqrt(np.mean(scaled ** 2)))
