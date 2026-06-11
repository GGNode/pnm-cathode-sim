"""
P2D 状态向量测试 (State Vector Tests)
=====================================

验证:
- StateLayout 构建与切片正确性
- pack/unpack 往返一致性
- 初始状态构造
- surface_concentration / average_solid_concentration

对应 TASK_PHASE0.md §5
"""

import numpy as np
import pytest

from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh
from pnmcathode.p2d.state import P2DState, StateLayout, make_initial_state, pack_state, unpack_state


@pytest.fixture
def meshes():
    sep = P2DRegion(
        name="separator", x_left=0, x_right=25e-6,
        n_cells=5, epsilon_e=0.39,
    )
    pos = P2DRegion(
        name="positive", x_left=25e-6, x_right=100e-6,
        n_cells=20, epsilon_e=0.35, epsilon_s=0.55,
        particle_radius=5e-6,
    )
    macro = MacroMesh.from_regions((sep, pos), area=1.0)
    particle = ParticleMesh.spherical(5e-6, 10)
    return macro, particle


@pytest.fixture
def layout(meshes):
    macro, particle = meshes
    return StateLayout.from_meshes(macro, particle)


class TestStateLayout:
    """StateLayout 测试。"""

    def test_size(self, layout):
        n_x = 25
        n_pos = 20
        n_r = 10
        expected = 2 * n_x + n_pos + n_pos * n_r
        assert layout.size == expected

    def test_slices(self, layout):
        n_x = 25
        n_pos = 20
        n_r = 10
        assert layout.c_e == slice(0, n_x)
        assert layout.phi_e == slice(n_x, 2 * n_x)
        assert layout.phi_s == slice(2 * n_x, 2 * n_x + n_pos)
        assert layout.c_s == slice(2 * n_x + n_pos, 2 * n_x + n_pos + n_pos * n_r)

    def test_positive_to_macro(self, layout):
        """positive cells 从 separator 之后开始。"""
        assert layout.positive_to_macro[0] == 5
        assert layout.positive_to_macro[-1] == 24

    def test_macro_to_positive(self, layout):
        """separator cells 映射到 -1。"""
        assert layout.macro_to_positive[0] == -1
        assert layout.macro_to_positive[4] == -1
        assert layout.macro_to_positive[5] == 0
        assert layout.macro_to_positive[24] == 19


class TestPackUnpack:
    """pack/unpack 往返测试。"""

    def test_roundtrip(self, layout, meshes):
        macro, particle = meshes
        state = make_initial_state(macro, particle, c_e0=1200.0, soc0=0.5, c_s_max=48900.0)
        y = pack_state(state, layout)
        state2 = unpack_state(y, layout)

        np.testing.assert_allclose(state.c_e, state2.c_e)
        np.testing.assert_allclose(state.phi_e, state2.phi_e)
        np.testing.assert_allclose(state.phi_s, state2.phi_s)
        np.testing.assert_allclose(state.c_s, state2.c_s)

    def test_y_size(self, layout, meshes):
        macro, particle = meshes
        state = make_initial_state(macro, particle, c_e0=1200.0, soc0=0.5, c_s_max=48900.0)
        y = pack_state(state, layout)
        assert len(y) == layout.size


class TestP2DState:
    """P2DState 方法测试。"""

    def test_surface_concentration(self, meshes):
        macro, particle = meshes
        state = make_initial_state(macro, particle, c_e0=1200.0, soc0=0.5, c_s_max=48900.0)
        c_s_surf = state.surface_concentration()
        assert len(c_s_surf) == 20
        np.testing.assert_allclose(c_s_surf, 0.5 * 48900.0)

    def test_average_solid_concentration_uniform(self, meshes):
        """均匀 c_s → 平均值等于 c_s 值。"""
        macro, particle = meshes
        state = make_initial_state(macro, particle, c_e0=1200.0, soc0=0.5, c_s_max=48900.0)
        c_s_avg = state.average_solid_concentration(particle)
        np.testing.assert_allclose(c_s_avg, 0.5 * 48900.0, rtol=1e-10)

    def test_copy_independence(self, meshes):
        macro, particle = meshes
        state = make_initial_state(macro, particle, c_e0=1200.0, soc0=0.5, c_s_max=48900.0)
        state2 = state.copy()
        state2.c_e[0] = 999.0
        assert state.c_e[0] == 1200.0


class TestInitialState:
    """初始状态测试。"""

    def test_c_e_uniform(self, meshes):
        macro, particle = meshes
        state = make_initial_state(macro, particle, c_e0=1200.0, soc0=0.5, c_s_max=48900.0)
        np.testing.assert_allclose(state.c_e, 1200.0)

    def test_phi_e_zero(self, meshes):
        macro, particle = meshes
        state = make_initial_state(macro, particle, c_e0=1200.0, soc0=0.5, c_s_max=48900.0)
        np.testing.assert_allclose(state.phi_e, 0.0)

    def test_phi_s_ocv(self, meshes):
        macro, particle = meshes
        from pnmcathode.physics.ocv import nmc532_ocv
        state = make_initial_state(macro, particle, c_e0=1200.0, soc0=0.5, c_s_max=48900.0)
        expected = float(nmc532_ocv(0.5))
        np.testing.assert_allclose(state.phi_s, expected)
