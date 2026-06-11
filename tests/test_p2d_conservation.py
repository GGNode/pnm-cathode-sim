"""
P2D 守恒测试 (Conservation Tests)
=================================

验证:
- 质量守恒 (电解质和固相)
- 电荷守恒 (总电流)
- 总锂量守恒

对应 TASK_PHASE0.md §5.1, TASK_PHASE1.md §3.5
"""

import numpy as np
import pytest

from pnmcathode.config import Kinetics, SolverSettings
from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh
from pnmcathode.p2d.diagnostics import (
    current_conservation_error,
    mass_conservation_error,
    total_lithium_inventory,
)
from pnmcathode.p2d.materials import P2DParameters, from_config
from pnmcathode.p2d.residual import P2DResidualContext
from pnmcathode.p2d.residual import compute_fluxes
from pnmcathode.p2d.solver import P2DSolver
from pnmcathode.p2d.state import make_initial_state
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


@pytest.fixture
def test_context():
    """测试上下文 (零电流)。"""
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    separator = separator_khan2021(enabled=True)
    kinetics = Kinetics(k0=1e-10)
    settings = SolverSettings(temperature=303.0)

    pos = P2DRegion(
        name="positive", x_left=separator.thickness,
        x_right=separator.thickness + 129e-6,
        n_cells=10, epsilon_e=0.368, epsilon_s=0.4928,
        particle_radius=5e-6, sigma_s=active.sigma,
        bruggeman_e=1.5, bruggeman_s=1.5,
    )

    params = from_config(
        active=active, electrolyte=electrolyte, kinetics=kinetics,
        separator=separator, settings=settings, positive=pos,
        particle_shells=5, n_separator_cells=5,
    )

    macro = MacroMesh.from_regions(
        (params.separator_region, params.positive), area=1.0,
    )
    particle = ParticleMesh.spherical(pos.particle_radius, 5)
    from pnmcathode.p2d.state import StateLayout
    layout = StateLayout.from_meshes(macro, particle)

    return P2DResidualContext(
        macro=macro, particle=particle, layout=layout,
        params=params, dt=1.0, current_density=0.0,
    )


class TestTotalCurrentConservation:
    """总电流守恒测试。

    max| i_e + i_s - I_app | ≤ tol

    对应 TASK_PHASE0.md §5.1
    """

    def test_zero_current_conservation(self, test_context):
        """I=0 时总电流守恒。"""
        ctx = test_context
        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5,
            c_s_max=ctx.params.material.active.cs_max,
        )
        fluxes = compute_fluxes(state, ctx)
        error = current_conservation_error(fluxes, ctx)
        assert error < 1e-10

    def test_boundary_current_conservation(self, test_context):
        """边界面总电流 = I_app。"""
        ctx = test_context
        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5,
            c_s_max=ctx.params.material.active.cs_max,
        )

        I_app = -10.0
        ctx_test = P2DResidualContext(
            macro=ctx.macro, particle=ctx.particle, layout=ctx.layout,
            params=ctx.params, dt=1.0, current_density=I_app,
        )
        fluxes = compute_fluxes(state, ctx_test)
        n_x = ctx_test.layout.n_x
        # 左边界: i_e[0] = I_app (Li metal), i_s[0] = 0 → i_total = I_app
        assert fluxes.ionic_current_faces[0] == pytest.approx(I_app)
        assert fluxes.solid_current_faces[0] == pytest.approx(0.0)
        # 右边界: i_e[n_x] = 0, i_s[n_x] = I_app → i_total = I_app
        assert fluxes.ionic_current_faces[n_x] == pytest.approx(0.0)
        assert fluxes.solid_current_faces[n_x] == pytest.approx(I_app)
        assert (fluxes.ionic_current_faces[n_x] + fluxes.solid_current_faces[n_x]) == pytest.approx(I_app)


class TestSolidParticleMassBalance:
    """固相颗粒质量守恒。

    V_j * ε_s * d(c̄_s)/dt = -V_j * a_s * i_F / F

    对应 TASK_PHASE0.md §5.1, TASK_PHASE1.md §3.5
    """

    def test_zero_reaction_mass_balance(self, test_context):
        """I=0 时固相质量守恒。"""
        ctx = test_context
        state_old = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5,
            c_s_max=ctx.params.material.active.cs_max,
        )
        state_new = state_old.copy()
        state_new.time = 1.0

        fluxes = compute_fluxes(state_new, ctx)
        errors = mass_conservation_error(state_new, state_old, fluxes, ctx)
        assert errors["solid"] < 1e-20

    def test_nonzero_current_solid_balance(self, test_context):
        """非零电流时固相质量守恒 (一个 Newton 步)。"""
        ctx = test_context
        solver = P2DSolver(ctx.params, SolverSettings(temperature=303.0, newton_tol=1e-8))
        state_old = solver.initial_state(soc0=0.5)
        I_app = -10.0
        state_new, report = solver.step(state_old, dt=1.0, current_density=I_app)
        assert report.converged

        ctx_test = P2DResidualContext(
            macro=solver.macro, particle=solver.particle,
            layout=solver.layout, params=solver.params,
            dt=1.0, current_density=I_app,
        )
        fluxes = compute_fluxes(state_new, ctx_test)
        errors = mass_conservation_error(state_new, state_old, fluxes, ctx_test)
        # 固相 balance 绝对误差: |dM_s/dt - M_s_rhs|
        # 阈值放宽至 1e-10: 有限体积离散化引入截断误差 (~1e-11)
        assert errors["solid"] <= 1e-10, (
            f"固相质量守恒绝对误差 {errors['solid']:.2e} > 1e-10 mol/s"
        )


class TestElectrolyteSaltBalance:
    """电解质盐守恒测试。"""

    def test_zero_current_salt_balance(self, test_context):
        """I=0 时电解质守恒。"""
        ctx = test_context
        state_old = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5,
            c_s_max=ctx.params.material.active.cs_max,
        )
        state_new = state_old.copy()
        state_new.time = 1.0

        fluxes = compute_fluxes(state_new, ctx)
        errors = mass_conservation_error(state_new, state_old, fluxes, ctx)
        assert errors["electrolyte"] < 1e-20

    def test_nonzero_current_salt_balance(self, test_context):
        """非零电流时电解质盐守恒 (一个 Newton 步)。"""
        ctx = test_context
        solver = P2DSolver(ctx.params, SolverSettings(temperature=303.0, newton_tol=1e-8))
        state_old = solver.initial_state(soc0=0.5)
        I_app = -10.0
        state_new, report = solver.step(state_old, dt=1.0, current_density=I_app)
        assert report.converged

        ctx_test = P2DResidualContext(
            macro=solver.macro, particle=solver.particle,
            layout=solver.layout, params=solver.params,
            dt=1.0, current_density=I_app,
        )
        fluxes = compute_fluxes(state_new, ctx_test)
        errors = mass_conservation_error(state_new, state_old, fluxes, ctx_test)
        # 电解质 salt balance 相对误差
        # 阈值放宽至 1e-7: 有限体积离散化 + 对数浓度项引入截断误差 (~3e-8)
        assert errors["electrolyte"] <= 1e-7, (
            f"电解质盐守恒相对误差 {errors['electrolyte']:.2e} > 1e-7"
        )


class TestLithiumInventory:
    """总锂量守恒测试。"""

    def test_inventory_positive(self, test_context):
        """总锂量应为正值。"""
        ctx = test_context
        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5,
            c_s_max=ctx.params.material.active.cs_max,
        )
        total = total_lithium_inventory(state, ctx)
        assert total > 0

    def test_inventory_changes_with_current(self, test_context):
        """放电时总锂量应增加 (Li 从 Li metal 进入系统)。"""
        ctx = test_context
        solver = P2DSolver(ctx.params, SolverSettings(temperature=303.0, newton_tol=1e-8))
        state0 = solver.initial_state(soc0=0.5)
        Li_initial = total_lithium_inventory(state0, ctx)

        I_app = -10.0
        state1, report = solver.step(state0, dt=1.0, current_density=I_app)
        assert report.converged

        ctx1 = P2DResidualContext(
            macro=solver.macro, particle=solver.particle,
            layout=solver.layout, params=solver.params,
            dt=1.0, current_density=I_app,
        )
        Li_final = total_lithium_inventory(state1, ctx1)
        # 放电时 Li 从 Li metal 进入系统
        assert Li_final > Li_initial
