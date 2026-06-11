"""
Phase 1: Khan 2021 电压曲线趋势测试 (Voltage Trend)
====================================================

验证:
- 初始电压在合理范围 (3.8 ~ 4.2 V)。
- 倍率越高电压越低。
- 0.2C vs 1C, 1C vs 3C 中位电压差 >= 0.02 V。

算例: khan_1cal_{0.2c,0.5c,1c,3c}

对应 TASK_PHASE1.md §3.6
"""

import numpy as np
import pytest

from pnmcathode.p2d.benchmark import (
    P2DBenchmarkCase,
    extract_voltage_curve,
    run_p2d_benchmark,
)
from pnmcathode.physics.ocv import nmc532_ocv


def _make_khan_case(name: str, c_rate: float) -> P2DBenchmarkCase:
    """构建 Khan benchmark case (小网格加速测试)。"""
    return P2DBenchmarkCase(
        name=name,
        parameter_source="khan2021",
        electrode_kind="1CAL",
        c_rate=c_rate,
        initial_soc=0.35,
        temperature=303.0,
        cutoff_voltage=3.0,
        n_positive_cells=12,
        n_separator_cells=4,
        n_particle_shells=6,
        dt_initial=10.0,
        t_final=500.0,
    )


class TestInitialVoltage:
    """初始电压范围测试。"""

    @pytest.mark.parametrize("c_rate", [0.2, 1.0, 3.0])
    def test_initial_voltage_range(self, c_rate):
        """初始电压应在 OCV(0.35) 附近合理范围。"""
        case = _make_khan_case(f"test_{c_rate}c", c_rate)
        result = run_p2d_benchmark(case)
        V0 = result.voltage[0]
        U0 = float(nmc532_ocv(case.initial_soc))
        # 低倍率接近 OCV，高倍率有 IR drop
        # V(0) 应在 [U0 - 0.5, U0 + 0.05] 范围内
        assert U0 - 0.5 <= V0 <= U0 + 0.05, (
            f"{case.name}: V(0) = {V0:.4f} V, "
            f"OCV({case.initial_soc}) = {U0:.4f} V"
        )


class TestRateTrend:
    """倍率趋势测试: 高倍率 → 低电压。"""

    def test_voltage_decreases_with_rate(self):
        """同容量区间内: V_0.2C > V_1C > V_3C。"""
        cases = {
            rate: _make_khan_case(f"1cal_{rate}c", rate)
            for rate in [0.2, 1.0, 3.0]
        }
        results = {rate: run_p2d_benchmark(case) for rate, case in cases.items()}

        # 在共同容量区间内插值比较
        caps = {}
        volts = {}
        for rate, res in results.items():
            t, c, v = extract_voltage_curve(res)
            caps[rate] = c
            volts[rate] = v

        q_min = max(c.min() for c in caps.values())
        q_max = min(c.max() for c in caps.values())
        if q_max <= q_min:
            pytest.skip("容量区间无重叠")
        q_grid = np.linspace(q_min, q_max, 50)

        v_interp = {}
        for rate in [0.2, 1.0, 3.0]:
            v_interp[rate] = np.interp(q_grid, caps[rate], volts[rate])

        # 0.2C vs 1C
        diff_02_1 = np.median(v_interp[0.2] - v_interp[1.0])
        assert diff_02_1 >= 0.02, (
            f"median(V_0.2C - V_1C) = {diff_02_1:.4f} V < 0.02 V"
        )

        # 1C vs 3C
        diff_1_3 = np.median(v_interp[1.0] - v_interp[3.0])
        assert diff_1_3 >= 0.02, (
            f"median(V_1C - V_3C) = {diff_1_3:.4f} V < 0.02 V"
        )


class TestDischargeBehavior:
    """放电行为测试。"""

    def test_voltage_decreases_during_discharge(self):
        """放电过程中电压应下降。"""
        case = _make_khan_case("1cal_1c_trend", 1.0)
        result = run_p2d_benchmark(case)
        V = result.voltage
        assert V[-1] < V[0], (
            f"最终电压 {V[-1]:.4f} >= 初始电压 {V[0]:.4f}"
        )

    def test_capacity_positive(self):
        """放电容量应为正值。"""
        case = _make_khan_case("1cal_1c_cap", 1.0)
        result = run_p2d_benchmark(case)
        _, cap, _ = extract_voltage_curve(result)
        assert cap[-1] > 0, f"最终容量 {cap[-1]:.4f} <= 0"


class TestOCVAtZeroCurrent:
    """零电流时初始电压 = OCV。"""

    def test_zero_rate_voltage_is_ocv(self):
        """极低倍率下初始电压应接近 OCV(0.35)。"""
        case = _make_khan_case("1cal_0p01c", 0.01)
        result = run_p2d_benchmark(case)
        V0 = result.voltage[0]
        U0 = float(nmc532_ocv(0.35))
        assert abs(V0 - U0) < 0.01, (
            f"|V(0) - U(0.35)| = {abs(V0 - U0):.4f} V >= 0.01 V"
        )
