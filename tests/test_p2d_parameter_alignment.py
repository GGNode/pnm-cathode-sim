"""
Phase 1: 参数快照与参数源一致性测试 (Parameter Alignment)
=========================================================

验证:
- Khan 2021 参数快照与论文 Table I/II 一致。
- 1CAL / 3CAL 电极参数正确。
- 参数源隔离断言通过。

对应 TASK_PHASE1.md §3.1, §6
"""

import pytest

from pnmcathode.p2d.benchmark import (
    P2DBenchmarkCase,
    P2DParameterSnapshot,
    assert_no_mixed_parameter_sources,
    make_khan_solver,
    make_parameter_snapshot,
)


def _make_case(electrode_kind: str) -> P2DBenchmarkCase:
    """构建测试用 benchmark case。"""
    return P2DBenchmarkCase(
        name=f"khan_{electrode_kind.lower()}_test",
        parameter_source="khan2021",
        electrode_kind=electrode_kind,
        c_rate=1.0,
        initial_soc=0.5,
        temperature=303.0,
        cutoff_voltage=3.0,
        n_positive_cells=10,
        n_separator_cells=5,
        n_particle_shells=5,
        dt_initial=1.0,
        t_final=10.0,
    )


class TestKhanParameterSnapshot:
    """Khan 2021 参数快照测试。"""

    def test_temperature(self):
        """温度应为 303.0 K。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert snap.temperature == 303.0

    def test_cs_max(self):
        """c_s_max 应为 48900.0 mol/m³。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert snap.c_s_max == 48900.0

    def test_ce_init(self):
        """c_e0 应为 1200.0 mol/m³。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert snap.c_e0 == 1200.0

    def test_separator_thickness(self):
        """隔膜厚度应为 25 μm。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert snap.separator_thickness == pytest.approx(25e-6)

    def test_t_plus(self):
        """Li+ 迁移数应为 0.363。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert snap.t_plus == pytest.approx(0.363)

    def test_k0(self):
        """反应速率常数应为 1e-10。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert snap.k0 == 1e-10


class TestElectrode1CAL:
    """1CAL 电极参数测试。"""

    @pytest.fixture
    def snap(self):
        solver = make_khan_solver(_make_case("1CAL"))
        return make_parameter_snapshot(solver, source="khan2021")

    def test_thickness(self, snap):
        """1CAL 厚度 = 129 μm。"""
        assert snap.positive_thickness == pytest.approx(129e-6)

    def test_porosity(self, snap):
        """1CAL 孔隙率 = 0.368。"""
        assert snap.positive_porosity == pytest.approx(0.368)

    def test_active_fraction(self, snap):
        """1CAL 活性材料体积分数 = 0.4928。"""
        assert snap.positive_active_fraction == pytest.approx(0.4928)


class TestElectrode3CAL:
    """3CAL 电极参数测试。"""

    @pytest.fixture
    def snap(self):
        solver = make_khan_solver(_make_case("3CAL"))
        return make_parameter_snapshot(solver, source="khan2021")

    def test_thickness(self, snap):
        """3CAL 厚度 = 88 μm。"""
        assert snap.positive_thickness == pytest.approx(88e-6)

    def test_porosity(self, snap):
        """3CAL 孔隙率 = 0.366。"""
        assert snap.positive_porosity == pytest.approx(0.366)

    def test_active_fraction(self, snap):
        """3CAL 活性材料体积分数 = 0.4974。"""
        assert snap.positive_active_fraction == pytest.approx(0.4974)


class TestBenchmarkMetadata:
    """Benchmark 输出 metadata 测试。"""

    def test_parameter_source_in_metadata(self):
        """metadata 中 parameter_source 应为 'khan2021' 或 'pybamm'。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert snap.source in {"khan2021", "pybamm"}


class TestParameterSourceIsolation:
    """参数源隔离断言测试。"""

    def test_khan_source_passes(self):
        """Khan 2021 参数快照应通过隔离断言。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert_no_mixed_parameter_sources(snap)

    def test_khan_ocv_reference(self):
        """Khan OCV reference 应为 'li_metal_empirical'。"""
        solver = make_khan_solver(_make_case("1CAL"))
        snap = make_parameter_snapshot(solver, source="khan2021")
        assert snap.ocv_reference == "li_metal_empirical"
