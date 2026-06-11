"""
Phase 1: 零电流平衡态与 OCV 一致性测试 (Half-Cell Equilibrium)
================================================================

验证:
- I=0 时浓度不漂移、过电位为零、BV 电流为零。
- 电压 = OCV(SOC)。
- 对 SOC = 0.35, 0.5, 0.75 均成立。

对应 TASK_PHASE1.md §3.2, §5.1
"""

import numpy as np
import pytest

from pnmcathode.config import Kinetics, SolverSettings
from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh
from pnmcathode.p2d.diagnostics import equilibrium_error
from pnmcathode.p2d.materials import from_config
from pnmcathode.p2d.residual import P2DResidualContext, compute_fluxes, compute_reaction
from pnmcathode.p2d.solver import P2DSolver
from pnmcathode.p2d.state import StateLayout, make_initial_state
from pnmcathode.physics.ocv import nmc532_ocv
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


def _make_equilibrium_solver(n_cells=10, n_shells=5):
    """构建平衡态测试用求解器。"""
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    separator = separator_khan2021(enabled=True)
    kinetics = Kinetics(k0=1e-10)
    settings = SolverSettings(temperature=303.0, newton_tol=1e-8)

    pos = P2DRegion(
        name="positive",
        x_left=separator.thickness,
        x_right=separator.thickness + 129e-6,
        n_cells=n_cells,
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


def _make_context(solver, dt=1.0, I_app=0.0):
    """构建残差上下文。"""
    return P2DResidualContext(
        macro=solver.macro,
        particle=solver.particle,
        layout=solver.layout,
        params=solver.params,
        dt=dt,
        current_density=I_app,
    )


class TestZeroCurrentEquilibrium:
    """零电流平衡态测试。

    对 SOC = 0.35, 0.5, 0.75 各跑 10 个步长:
    - max|c_e(t) - c_e(0)| / c_e0 <= 1e-10
    - max|c_s(t) - c_s(0)| / c_s_max <= 1e-10
    - max|eta| <= 1e-8 V
    - max|i_F| <= 1e-9 A/m²
    - |V - U(SOC)| <= 1e-8 V
    """

    @pytest.mark.parametrize("soc0", [0.35, 0.5, 0.75])
    def test_ce_stable(self, soc0):
        """电解质浓度不漂移。"""
        solver = _make_equilibrium_solver()
        state0 = solver.initial_state(soc0=soc0)
        c_e0 = state0.c_e.copy()

        state = state0.copy()
        for _ in range(10):
            state, report = solver.step(state, dt=1.0, current_density=0.0)
            assert report.converged, f"Newton 未收敛: {report}"

        c_e0_val = solver.params.material.electrolyte.c_init
        drift = np.max(np.abs(state.c_e - c_e0)) / c_e0_val
        assert drift <= 1e-10, (
            f"SOC={soc0}: c_e 漂移 {drift:.2e} > 1e-10"
        )

    @pytest.mark.parametrize("soc0", [0.35, 0.5, 0.75])
    def test_cs_stable(self, soc0):
        """固相浓度不漂移。"""
        solver = _make_equilibrium_solver()
        state0 = solver.initial_state(soc0=soc0)
        c_s0 = state0.c_s.copy()

        state = state0.copy()
        for _ in range(10):
            state, report = solver.step(state, dt=1.0, current_density=0.0)
            assert report.converged

        cs_max = solver.params.material.active.cs_max
        drift = np.max(np.abs(state.c_s - c_s0)) / cs_max
        assert drift <= 1e-10, (
            f"SOC={soc0}: c_s 漂移 {drift:.2e} > 1e-10"
        )

    @pytest.mark.parametrize("soc0", [0.35, 0.5, 0.75])
    def test_overpotential_zero(self, soc0):
        """过电位应为零。"""
        solver = _make_equilibrium_solver()
        state = solver.initial_state(soc0=soc0)

        # 跑 10 步使 Newton 充分收敛
        for _ in range(10):
            state, report = solver.step(state, dt=1.0, current_density=0.0)
            assert report.converged

        ctx = _make_context(solver)
        errors = equilibrium_error(state, ctx)
        assert errors["eta_max"] <= 1e-8, (
            f"SOC={soc0}: max|eta| = {errors['eta_max']:.2e} > 1e-8 V"
        )

    @pytest.mark.parametrize("soc0", [0.35, 0.5, 0.75])
    def test_faradaic_current_zero(self, soc0):
        """BV Faradaic 电流应为零。"""
        solver = _make_equilibrium_solver()
        state = solver.initial_state(soc0=soc0)

        for _ in range(10):
            state, report = solver.step(state, dt=1.0, current_density=0.0)
            assert report.converged

        ctx = _make_context(solver)
        errors = equilibrium_error(state, ctx)
        assert errors["i_f_max"] <= 1e-9, (
            f"SOC={soc0}: max|i_F| = {errors['i_f_max']:.2e} > 1e-9 A/m²"
        )

    @pytest.mark.parametrize("soc0", [0.35, 0.5, 0.75])
    def test_voltage_equals_ocv(self, soc0):
        """零电流电压应等于 OCV(SOC)。"""
        solver = _make_equilibrium_solver()
        state = solver.initial_state(soc0=soc0)

        for _ in range(10):
            state, report = solver.step(state, dt=1.0, current_density=0.0)
            assert report.converged

        V = solver.voltage(state)
        U_expected = float(nmc532_ocv(soc0))
        assert abs(V - U_expected) <= 1e-8, (
            f"SOC={soc0}: |V - U| = {abs(V - U_expected):.2e} > 1e-8 V "
            f"(V={V:.10f}, U={U_expected:.10f})"
        )


class TestOCVConsistency:
    """OCV 一致性测试。

    - solver 电压 = phi_s[-1] (gauge 下 phi_e[0]=0)。
    - 零电流时 V = U(SOC)。
    """

    def test_voltage_is_phi_s_last(self):
        """voltage() 应返回 phi_s[-1]。"""
        solver = _make_equilibrium_solver()
        state = solver.initial_state(soc0=0.5)
        V = solver.voltage(state)
        assert V == pytest.approx(float(state.phi_s[-1]))

    @pytest.mark.parametrize("soc0", [0.35, 0.5, 0.75])
    def test_voltage_equals_ocv(self, soc0):
        """零电流电压等于 OCV。"""
        solver = _make_equilibrium_solver()
        state = solver.initial_state(soc0=soc0)
        V = solver.voltage(state)
        U = float(nmc532_ocv(soc0))
        assert abs(V - U) <= 1e-8, (
            f"|V - U({soc0})| = {abs(V - U):.2e} V"
        )


class TestOCVModelAlignment:
    """P2D OCV 模型与 nmc532_ocv 对齐测试。"""

    def test_p2d_ocv_matches_nmc532(self):
        """from_config 构造的 OCV 模型与 nmc532_ocv 完全一致。"""
        solver = _make_equilibrium_solver()
        ocv_model = solver.params.material.ocv_model

        soc_grid = np.linspace(0.3, 0.9, 13)
        u_p2d = ocv_model.value(soc_grid)
        u_ref = nmc532_ocv(soc_grid)

        np.testing.assert_allclose(u_p2d, u_ref, atol=1e-12)

    def test_ocv_derivative_vs_finite_diff(self):
        """OCV derivative 与中心差分一致 (避开 clip 区间)。"""
        solver = _make_equilibrium_solver()
        ocv_model = solver.params.material.ocv_model

        soc_grid = np.linspace(0.35, 0.85, 11)
        du_analytic = ocv_model.derivative(soc_grid)

        h = 1e-6
        du_fd = (ocv_model.value(soc_grid + h) - ocv_model.value(soc_grid - h)) / (2 * h)

        np.testing.assert_allclose(du_analytic, du_fd, atol=1e-4)
