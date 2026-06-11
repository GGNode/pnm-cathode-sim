"""
Phase 1: 非零电流边界条件与总电流守恒测试 (Boundary Current)
=============================================================

验证:
- 内部 faces: i_e + i_s = I_app
- separator/positive interface: i_s = 0
- separator 内部: i_e = I_app
- cathode current collector: i_s = I_app, i_e = 0

对 I_app = -1, -10, -50 A/m² 各跑一个收敛时间步。

对应 TASK_PHASE1.md §3.4, §5
"""

import numpy as np
import pytest

from pnmcathode.config import Kinetics, SolverSettings
from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh
from pnmcathode.p2d.materials import from_config
from pnmcathode.p2d.residual import P2DResidualContext, compute_fluxes
from pnmcathode.p2d.solver import P2DSolver
from pnmcathode.p2d.state import StateLayout, make_initial_state
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


def _make_boundary_solver(n_pos=10, n_shells=5):
    """构建边界电流测试用求解器。"""
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    separator = separator_khan2021(enabled=True)
    kinetics = Kinetics(k0=1e-10)
    settings = SolverSettings(temperature=303.0, newton_tol=1e-8)

    pos = P2DRegion(
        name="positive",
        x_left=separator.thickness,
        x_right=separator.thickness + 129e-6,
        n_cells=n_pos,
        epsilon_e=0.368,
        epsilon_s=0.4928,
        particle_radius=5e-6,
        sigma_s=active.sigma,
        bruggeman_e=1.5,
        bruggeman_s=1.5,
    )

    params = from_config(
        active=active,
        electrolyte=electrolyte,
        kinetics=kinetics,
        separator=separator,
        settings=settings,
        positive=pos,
        particle_shells=n_shells,
        n_separator_cells=5,
    )
    return P2DSolver(params, settings)


def total_current_error(state, context):
    """计算内部 faces 总电流守恒误差。"""
    fluxes = compute_fluxes(state, context)
    i_total = fluxes.ionic_current_faces + fluxes.solid_current_faces
    interior = slice(1, context.layout.n_x)
    scale = max(abs(context.current_density), 1.0)
    return float(np.max(np.abs(i_total[interior] - context.current_density)) / scale)


class TestTotalCurrentConservation:
    """总电流守恒: i_e + i_s = I_app。"""

    @pytest.mark.parametrize("I_app", [-1.0, -10.0, -50.0])
    def test_interior_total_current(self, I_app):
        """内部 faces 总电流 = I_app。"""
        solver = _make_boundary_solver()
        state0 = solver.initial_state(soc0=0.5)

        # 跑一个收敛时间步
        state, report = solver.step(state0, dt=1.0, current_density=I_app)
        assert report.converged, f"I_app={I_app}: Newton 未收敛"

        ctx = P2DResidualContext(
            macro=solver.macro,
            particle=solver.particle,
            layout=solver.layout,
            params=solver.params,
            dt=1.0,
            current_density=I_app,
        )
        error = total_current_error(state, ctx)
        assert error <= 1e-8, (
            f"I_app={I_app}: 总电流守恒误差 {error:.2e} > 1e-8"
        )


class TestSeparatorPositiveInterface:
    """Separator/positive 界面条件。"""

    @pytest.mark.parametrize("I_app", [-1.0, -10.0, -50.0])
    def test_solid_current_at_interface(self, I_app):
        """separator/positive 界面处 i_s = 0。"""
        solver = _make_boundary_solver()
        state0 = solver.initial_state(soc0=0.5)
        state, report = solver.step(state0, dt=1.0, current_density=I_app)
        assert report.converged

        ctx = P2DResidualContext(
            macro=solver.macro,
            particle=solver.particle,
            layout=solver.layout,
            params=solver.params,
            dt=1.0,
            current_density=I_app,
        )
        fluxes = compute_fluxes(state, ctx)

        # separator/positive 界面: separator 最后一个 cell 的右 face
        n_sep = len(solver.macro.separator_cells)
        # face n_sep 是 separator/positive 界面
        i_s_interface = fluxes.solid_current_faces[n_sep]
        assert abs(i_s_interface) <= 1e-10, (
            f"I_app={I_app}: 界面 i_s = {i_s_interface:.2e}, 应为 0"
        )


class TestSeparatorIonicCurrent:
    """Separator 内部 ionic current。"""

    @pytest.mark.parametrize("I_app", [-1.0, -10.0, -50.0])
    def test_separator_i_e_equals_I_app(self, I_app):
        """separator 内部 i_e = I_app。"""
        solver = _make_boundary_solver()
        state0 = solver.initial_state(soc0=0.5)
        state, report = solver.step(state0, dt=1.0, current_density=I_app)
        assert report.converged

        ctx = P2DResidualContext(
            macro=solver.macro,
            particle=solver.particle,
            layout=solver.layout,
            params=solver.params,
            dt=1.0,
            current_density=I_app,
        )
        fluxes = compute_fluxes(state, ctx)

        n_sep = len(solver.macro.separator_cells)
        # separator 内部 faces: 1 到 n_sep-1
        for j in range(1, n_sep):
            i_e_j = fluxes.ionic_current_faces[j]
            rel_err = abs(i_e_j - I_app) / max(abs(I_app), 1.0)
            assert rel_err <= 1e-8, (
                f"I_app={I_app}: separator face {j} i_e = {i_e_j:.6f}, "
                f"相对误差 {rel_err:.2e}"
            )


class TestCathodeCurrentCollector:
    """Cathode current collector 右边界。"""

    @pytest.mark.parametrize("I_app", [-1.0, -10.0, -50.0])
    def test_solid_current_at_collector(self, I_app):
        """右边界 i_s = I_app。"""
        solver = _make_boundary_solver()
        state0 = solver.initial_state(soc0=0.5)
        state, report = solver.step(state0, dt=1.0, current_density=I_app)
        assert report.converged

        ctx = P2DResidualContext(
            macro=solver.macro,
            particle=solver.particle,
            layout=solver.layout,
            params=solver.params,
            dt=1.0,
            current_density=I_app,
        )
        fluxes = compute_fluxes(state, ctx)

        n_x = solver.layout.n_x
        i_s_right = fluxes.solid_current_faces[n_x]
        assert abs(i_s_right - I_app) <= 1e-10, (
            f"I_app={I_app}: 右边界 i_s = {i_s_right:.6f}, 应为 {I_app}"
        )

    @pytest.mark.parametrize("I_app", [-1.0, -10.0, -50.0])
    def test_ionic_current_at_collector(self, I_app):
        """右边界 i_e = 0。"""
        solver = _make_boundary_solver()
        state0 = solver.initial_state(soc0=0.5)
        state, report = solver.step(state0, dt=1.0, current_density=I_app)
        assert report.converged

        ctx = P2DResidualContext(
            macro=solver.macro,
            particle=solver.particle,
            layout=solver.layout,
            params=solver.params,
            dt=1.0,
            current_density=I_app,
        )
        fluxes = compute_fluxes(state, ctx)

        n_x = solver.layout.n_x
        i_e_right = fluxes.ionic_current_faces[n_x]
        assert abs(i_e_right) <= 1e-10, (
            f"I_app={I_app}: 右边界 i_e = {i_e_right:.2e}, 应为 0"
        )


class TestLeftBoundaryIonicCurrent:
    """Li metal 左边界 ionic current。"""

    @pytest.mark.parametrize("I_app", [-1.0, -10.0, -50.0])
    def test_left_i_e_equals_I_app(self, I_app):
        """左边界 i_e[0] = I_app。"""
        solver = _make_boundary_solver()
        state0 = solver.initial_state(soc0=0.5)
        state, report = solver.step(state0, dt=1.0, current_density=I_app)
        assert report.converged

        ctx = P2DResidualContext(
            macro=solver.macro,
            particle=solver.particle,
            layout=solver.layout,
            params=solver.params,
            dt=1.0,
            current_density=I_app,
        )
        fluxes = compute_fluxes(state, ctx)

        i_e_left = fluxes.ionic_current_faces[0]
        assert abs(i_e_left - I_app) <= 1e-10, (
            f"I_app={I_app}: 左边界 i_e = {i_e_left:.6f}, 应为 {I_app}"
        )

    @pytest.mark.parametrize("I_app", [-1.0, -10.0, -50.0])
    def test_left_solid_current_zero(self, I_app):
        """左边界 i_s[0] = 0 (separator 中无固相)。"""
        solver = _make_boundary_solver()
        state0 = solver.initial_state(soc0=0.5)
        state, report = solver.step(state0, dt=1.0, current_density=I_app)
        assert report.converged

        ctx = P2DResidualContext(
            macro=solver.macro,
            particle=solver.particle,
            layout=solver.layout,
            params=solver.params,
            dt=1.0,
            current_density=I_app,
        )
        fluxes = compute_fluxes(state, ctx)

        i_s_left = fluxes.solid_current_faces[0]
        assert abs(i_s_left) <= 1e-10, (
            f"I_app={I_app}: 左边界 i_s = {i_s_left:.2e}, 应为 0"
        )
