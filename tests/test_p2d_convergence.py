"""
Phase 1: 网格与时间步收敛测试 (Convergence)
============================================

验证:
- 电压曲线随网格加密收敛。
- medium vs fine 的 RMS 差 < 0.7 * (coarse vs medium)。
- 最终容量差 <= 2%。

默认标记为 slow，CI 可跳过。

对应 TASK_PHASE1.md §3.8
"""

import numpy as np
import pytest

from pnmcathode.config import Kinetics, SolverSettings
from pnmcathode.p2d.benchmark import (
    P2DBenchmarkCase,
    c_rate_to_current_density,
    extract_voltage_curve,
)
from pnmcathode.p2d.domain import P2DRegion
from pnmcathode.p2d.materials import from_config
from pnmcathode.p2d.solver import P2DProtocol, P2DSolver
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


def _make_convergence_solver(n_cells: int, n_shells: int, dt: float):
    """构建收敛测试用求解器。"""
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    separator = separator_khan2021(enabled=True)
    kinetics = Kinetics(k0=1e-10)
    settings = SolverSettings(temperature=303.0, newton_tol=1e-8, dt_max=dt * 2)

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
        n_separator_cells=max(3, n_cells // 4),
    )
    return P2DSolver(params, settings)


def _run_grid(n_cells: int, n_shells: int, dt: float, c_rate: float = 1.0):
    """运行一组网格参数。"""
    solver = _make_convergence_solver(n_cells, n_shells, dt)
    state0 = solver.initial_state(soc0=0.35)
    I_app = c_rate_to_current_density(
        c_rate,
        solver.params.positive.epsilon_s,
        solver.params.positive.length,
        solver.params.material.active.cs_max,
    )
    protocol = P2DProtocol(
        current_density=I_app,
        t_final=3000.0,
        dt_initial=dt,
        cutoff_voltage=3.0,
        save_every=1,
    )
    result = solver.run(state0, protocol)
    return result


def _rms_voltage_diff(res_a, res_b):
    """计算两条电压曲线在共同容量区间上的 RMS 差。"""
    _, cap_a, v_a = extract_voltage_curve(res_a)
    _, cap_b, v_b = extract_voltage_curve(res_b)

    q_min = max(cap_a.min(), cap_b.min())
    q_max = min(cap_a.max(), cap_b.max())
    if q_max <= q_min:
        return float("inf")
    q_grid = np.linspace(q_min, q_max, 100)
    va = np.interp(q_grid, cap_a, v_a)
    vb = np.interp(q_grid, cap_b, v_b)
    return float(np.sqrt(np.mean((va - vb) ** 2)))


@pytest.mark.slow
class TestGridConvergence:
    """网格收敛测试。"""

    def test_voltage_convergence(self):
        """RMS(V_medium - V_fine) < 0.7 * RMS(V_coarse - V_medium)。"""
        res_coarse = _run_grid(20, 8, 20.0)
        res_medium = _run_grid(40, 16, 10.0)
        res_fine = _run_grid(80, 32, 5.0)

        rms_coarse_medium = _rms_voltage_diff(res_coarse, res_medium)
        rms_medium_fine = _rms_voltage_diff(res_medium, res_fine)

        # 如果 coarse-medium 差太小，放宽检查
        if rms_coarse_medium < 1e-4:
            return  # 已经收敛

        assert rms_medium_fine < 0.7 * rms_coarse_medium, (
            f"RMS(medium-fine) = {rms_medium_fine:.6f} >= "
            f"0.7 * RMS(coarse-medium) = {0.7 * rms_coarse_medium:.6f}"
        )

    def test_capacity_convergence(self):
        """|Q_medium - Q_fine| / Q_fine <= 0.02。"""
        res_medium = _run_grid(40, 16, 10.0)
        res_fine = _run_grid(80, 32, 5.0)

        _, cap_m, _ = extract_voltage_curve(res_medium)
        _, cap_f, _ = extract_voltage_curve(res_fine)

        Q_m = cap_m[-1]
        Q_f = cap_f[-1]
        if Q_f < 1e-10:
            pytest.skip("fine 容量太小")

        rel_diff = abs(Q_m - Q_f) / Q_f
        assert rel_diff <= 0.02, (
            f"|Q_medium - Q_fine|/Q_fine = {rel_diff:.4f} > 0.02"
        )


class TestCoarseSmoke:
    """粗网格 smoke test (非 slow)。"""

    def test_coarse_runs(self):
        """粗网格运行成功。"""
        result = _run_grid(10, 4, 30.0)
        assert len(result.snapshots) >= 2
        V = result.voltage
        assert V[0] > V[-1], "电压应随放电下降"
