"""
P2D 残差测试 (Residual Tests)
=============================

验证:
- 残差维度正确
- 平衡态残差为零
- 边界条件正确施加
- 符号一致性

对应 TASK_PHASE0.md §5.2
"""

import numpy as np
import pytest

from pnmcathode.config import ActiveMaterial, Electrolyte, Kinetics, Separator, SolverSettings
from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh
from pnmcathode.p2d.materials import P2DParameters, P2DMaterial, OCVModel, from_config
from pnmcathode.p2d.residual import P2DResidualContext, assemble_p2d_residual
from pnmcathode.p2d.state import P2DState, StateLayout, make_initial_state, pack_state
from pnmcathode.physics.ocv import nmc532_ocv, ocv_derivative
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


@pytest.fixture
def default_context():
    """默认测试上下文 (小网格)。"""
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
    layout = StateLayout.from_meshes(macro, particle)

    return P2DResidualContext(
        macro=macro, particle=particle, layout=layout,
        params=params, dt=1.0, current_density=0.0,
    )


class TestResidualDimension:
    """残差维度测试。"""

    def test_residual_size(self, default_context):
        ctx = default_context
        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5,
            c_s_max=ctx.params.material.active.cs_max,
        )
        y = pack_state(state, ctx.layout)
        R = assemble_p2d_residual(y, state, ctx)
        assert len(R) == ctx.layout.size

    def test_residual_finite(self, default_context):
        """平衡态残差全部有限。"""
        ctx = default_context
        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5,
            c_s_max=ctx.params.material.active.cs_max,
        )
        y = pack_state(state, ctx.layout)
        R = assemble_p2d_residual(y, state, ctx)
        assert np.all(np.isfinite(R))


class TestOpenCircuitEquilibrium:
    """平衡态残差测试 (I_app = 0)。

    对应 TASK_PHASE0.md §5.2:
    - η = 0
    - i_F = 0
    - 所有 concentration residual 为 0
    - voltage = OCV
    """

    def test_equilibrium_residual_zero(self, default_context):
        """I=0, 均匀 c_e, c_s, phi_e=0, phi_s=U(soc) → R=0。"""
        ctx = default_context
        layout = ctx.layout
        cs_max = ctx.params.material.active.cs_max
        c_e0 = 1200.0
        soc0 = 0.5

        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=c_e0, soc0=soc0, c_s_max=cs_max,
        )

        # 设置 phi_s 为精确的 OCV 值
        U = float(nmc532_ocv(soc0))
        state.phi_s[:] = U

        y = pack_state(state, layout)
        R = assemble_p2d_residual(y, state, ctx)

        # 浓度残差应为零 (无反应源项)
        np.testing.assert_allclose(R[layout.c_e], 0.0, atol=1e-20)

        # phi_e gauge: R[phi_e][0] = phi_e[0] = 0
        assert R[layout.phi_e][0] == pytest.approx(0.0)

        # phi_s 残差: 无源项 → 应为零
        np.testing.assert_allclose(R[layout.phi_s], 0.0, atol=1e-20)

        # c_s 残差: 无表面通量 → 应为零
        np.testing.assert_allclose(R[layout.c_s], 0.0, atol=1e-20)

    def test_equilibrium_phi_e_gauge(self, default_context):
        """phi_e[0] 应被 gauge 约束为 0。"""
        ctx = default_context
        layout = ctx.layout
        cs_max = ctx.params.material.active.cs_max

        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5, c_s_max=cs_max,
        )
        y = pack_state(state, layout)
        R = assemble_p2d_residual(y, state, ctx)

        # gauge 残差 = phi_e[0]
        assert R[layout.phi_e][0] == pytest.approx(0.0)


class TestBoundaryConditions:
    """边界条件测试。"""

    def test_no_flux_left_boundary(self, default_context):
        """左边界无通量: N_e(0) = 0, i_e(0) = 0。"""
        ctx = default_context
        layout = ctx.layout
        cs_max = ctx.params.material.active.cs_max

        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5, c_s_max=cs_max,
        )
        y = pack_state(state, layout)
        R = assemble_p2d_residual(y, state, ctx)

        # 左边界 c_e 残差: 只有时间导数项 (无通量)
        # 对于均匀状态, 时间导数也为 0
        assert abs(R[layout.c_e][0]) < 1e-20

    def test_solid_current_right_boundary(self, default_context):
        """右边界固相电流 = I_app。"""
        ctx = default_context
        layout = ctx.layout
        cs_max = ctx.params.material.active.cs_max

        state = make_initial_state(
            ctx.macro, ctx.particle, c_e0=1200.0, soc0=0.5, c_s_max=cs_max,
        )
        y = pack_state(state, layout)

        # 设置非零电流
        ctx_test = P2DResidualContext(
            macro=ctx.macro, particle=ctx.particle, layout=layout,
            params=ctx.params, dt=1.0, current_density=-10.0,
        )
        R = assemble_p2d_residual(y, state, ctx_test)

        # 固相残差应非零 (因为 I_app ≠ 0)
        assert not np.allclose(R[layout.phi_s], 0.0, atol=1e-20)
