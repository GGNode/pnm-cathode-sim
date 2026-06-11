"""
P2D 守恒测试 (Conservation Tests)
=================================

验证:
- 质量守恒 (电解质和固相)
- 电荷守恒 (总电流)
- 总锂量守恒

对应 TASK_PHASE0.md §5.1
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
from pnmcathode.p2d.state import make_initial_state
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


@pytest.fixture
def test_context():
    """测试上下文。"""
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    separator = separator_khan2021(enabled=True)
    kinetics = Kinetics(k0=5e-10)
    settings = SolverSettings(temperature=303.0)

    pos = P2DRegion(
        name="positive", x_left=separator.thickness,
        x_right=separator.thickness + 75e-6,
        n_cells=10, epsilon_e=0.35, epsilon_s=0.55,
        particle_radius=5e-6, sigma_s=active.sigma,
    )

    params = from_config(
        active=active, electrolyte=electrolyte, kinetics=kinetics,
        separator=separator, settings=settings, positive=pos,
        particle_shells=5,
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
        from pnmcathode.p2d.state import pack_state

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
        # 左边界: i_e[0] = 0, i_s[0] = 0 → i_total = 0 (gauge face)
        # 右边界: i_e[n_x] = 0, i_s[n_x] = I_app → i_total = I_app
        n_x = ctx_test.layout.n_x
        assert fluxes.ionic_current_faces[n_x] == pytest.approx(0.0)
        assert fluxes.solid_current_faces[n_x] == pytest.approx(I_app)
        assert (fluxes.ionic_current_faces[n_x] + fluxes.solid_current_faces[n_x]) == pytest.approx(I_app)


class TestSolidParticleMassBalance:
    """固相颗粒质量守恒。

    V_j * ε_s * d(c̄_s)/dt = -V_j * a_s * i_F / F

    对应 TASK_PHASE0.md §5.1
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
