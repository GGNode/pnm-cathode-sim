"""
P2D 网格与区域定义测试 (Domain & Mesh Tests)
=============================================

验证:
- P2DRegion 属性计算
- MacroMesh 构建与区域索引
- ParticleMesh 球壳几何
- face averaging 与面积守恒

对应 TASK_PHASE0.md §5
"""

import numpy as np
import pytest

from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh


class TestP2DRegion:
    """P2DRegion 基本属性测试。"""

    def test_length(self):
        reg = P2DRegion(
            name="positive", x_left=25e-6, x_right=100e-6,
            n_cells=10, epsilon_e=0.35, epsilon_s=0.55,
            particle_radius=5e-6,
        )
        assert reg.length == pytest.approx(75e-6)

    def test_area_density_spherical(self):
        """球形颗粒: a_s = 3 * epsilon_s / R_p"""
        reg = P2DRegion(
            name="positive", x_left=0.0, x_right=75e-6,
            n_cells=10, epsilon_e=0.35, epsilon_s=0.55,
            particle_radius=5e-6,
        )
        expected = 3.0 * 0.55 / 5e-6
        assert reg.area_density() == pytest.approx(expected)

    def test_area_density_override(self):
        """外部给定 specific_area 优先。"""
        reg = P2DRegion(
            name="positive", x_left=0.0, x_right=75e-6,
            n_cells=10, epsilon_e=0.35, epsilon_s=0.55,
            particle_radius=5e-6, specific_area=1e5,
        )
        assert reg.area_density() == pytest.approx(1e5)

    def test_validation(self):
        """参数校验。"""
        with pytest.raises(ValueError, match="n_cells"):
            P2DRegion(name="positive", x_left=0, x_right=1, n_cells=0, epsilon_e=0.35)
        with pytest.raises(ValueError, match="x_right"):
            P2DRegion(name="positive", x_left=1, x_right=0, n_cells=10, epsilon_e=0.35)


class TestMacroMesh:
    """MacroMesh 构建与索引测试。"""

    @pytest.fixture
    def mesh(self):
        sep = P2DRegion(
            name="separator", x_left=0, x_right=25e-6,
            n_cells=5, epsilon_e=0.39,
        )
        pos = P2DRegion(
            name="positive", x_left=25e-6, x_right=100e-6,
            n_cells=20, epsilon_e=0.35, epsilon_s=0.55,
            particle_radius=5e-6,
        )
        return MacroMesh.from_regions((sep, pos), area=1.0)

    def test_total_cells(self, mesh):
        assert mesh.n_x == 25

    def test_positive_cells(self, mesh):
        assert mesh.n_pos == 20

    def test_separator_cells(self, mesh):
        assert len(mesh.separator_cells) == 5

    def test_x_faces_span(self, mesh):
        assert mesh.x_faces[0] == pytest.approx(0.0)
        assert mesh.x_faces[-1] == pytest.approx(100e-6)

    def test_dx_sum(self, mesh):
        """总厚度 = 区域长度之和。"""
        assert np.sum(mesh.dx) == pytest.approx(100e-6)

    def test_volumes(self, mesh):
        """体积 = dx * area。"""
        np.testing.assert_allclose(mesh.volumes, mesh.dx * 1.0)

    def test_region_for_cell(self, mesh):
        """区域索引正确。"""
        sep = mesh.regions[0]
        pos = mesh.regions[1]
        assert mesh.region_for_cell(0).name == sep.name
        assert mesh.region_for_cell(4).name == sep.name
        assert mesh.region_for_cell(5).name == pos.name
        assert mesh.region_for_cell(24).name == pos.name

    def test_face_indices(self, mesh):
        assert mesh.left_face(0) == 0
        assert mesh.right_face(0) == 1
        assert mesh.left_face(24) == 24
        assert mesh.right_face(24) == 25

    def test_centers_midpoint(self, mesh):
        """中心 = 面中点。"""
        for j in range(mesh.n_x):
            assert mesh.x_centers[j] == pytest.approx(
                0.5 * (mesh.x_faces[j] + mesh.x_faces[j + 1])
            )


class TestParticleMesh:
    """ParticleMesh 球壳几何测试。"""

    @pytest.fixture
    def pmesh(self):
        return ParticleMesh.spherical(radius=5e-6, n_shells=10)

    def test_shell_count(self, pmesh):
        assert pmesh.n_shells == 10

    def test_total_volume(self, pmesh):
        """壳层体积之和 = 球体积。"""
        V_total = (4.0 / 3.0) * np.pi * (5e-6) ** 3
        assert np.sum(pmesh.shell_volumes) == pytest.approx(V_total, rel=1e-10)

    def test_center_face_area_zero(self, pmesh):
        """r=0 处 face area = 0 (中心对称)。"""
        assert pmesh.face_areas[0] == pytest.approx(0.0)

    def test_surface_face_area(self, pmesh):
        """表面 face area = 4*pi*R_p^2。"""
        A_expected = 4.0 * np.pi * (5e-6) ** 2
        assert pmesh.face_areas[-1] == pytest.approx(A_expected)

    def test_r_faces_span(self, pmesh):
        assert pmesh.r_faces[0] == pytest.approx(0.0)
        assert pmesh.r_faces[-1] == pytest.approx(5e-6)

    def test_shell_volumes_positive(self, pmesh):
        assert np.all(pmesh.shell_volumes > 0)

    def test_validation(self):
        with pytest.raises(ValueError, match="radius"):
            ParticleMesh.spherical(0, 10)
        with pytest.raises(ValueError, match="n_shells"):
            ParticleMesh.spherical(5e-6, 0)
