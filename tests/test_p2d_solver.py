"""
P2D 求解器测试 (Solver Tests)
=============================

验证:
- Newton 收敛
- 零电流平衡态不漂移
- 低倍率放电电压接近 OCV
- 高倍率 smoke test

对应 TASK_PHASE0.md §5.2, §5.5
"""

import numpy as np
import pytest

from pnmcathode.config import Kinetics, SolverSettings
from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh
from pnmcathode.p2d.materials import P2DParameters, from_config
from pnmcathode.p2d.solver import P2DProtocol, P2DSolver
from pnmcathode.p2d.state import make_initial_state
from pnmcathode.physics.ocv import nmc532_ocv
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


def make_test_solver(n_cells_pos=5, n_shells=3, k0=5e-10):
    """构建测试用求解器 (小网格加速测试)。"""
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    separator = separator_khan2021(enabled=True)
    kinetics = Kinetics(k0=k0)
    settings = SolverSettings(temperature=303.0, newton_tol=1e-8)

    pos = P2DRegion(
        name="positive", x_left=separator.thickness,
        x_right=separator.thickness + 75e-6,
        n_cells=n_cells_pos, epsilon_e=0.35, epsilon_s=0.55,
        particle_radius=5e-6, sigma_s=active.sigma,
    )

    params = from_config(
        active=active, electrolyte=electrolyte, kinetics=kinetics,
        separator=separator, settings=settings, positive=pos,
        particle_shells=n_shells,
    )
    return P2DSolver(params, settings)


class TestNewtonConvergence:
    """Newton 收敛测试。"""

    def test_single_step_zero_current(self):
        """I=0 单步 Newton 应收敛。"""
        solver = make_test_solver()
        state0 = solver.initial_state(soc0=0.5)
        state1, report = solver.step(state0, dt=1.0, current_density=0.0)
        assert report.converged, f"Newton did not converge: {report}"

    def test_single_step_small_current(self):
        """小电流单步应收敛。"""
        solver = make_test_solver()
        state0 = solver.initial_state(soc0=0.5)
        state1, report = solver.step(state0, dt=0.01, current_density=-0.1)
        assert report.converged, f"Newton did not converge: {report}"


class TestZeroCurrentStateDoesNotDrift:
    """零电流平衡态不漂移测试。

    对应 TASK_PHASE0.md §5.2:
    - I=0 时 c_e 不变, c_s 不变, phi_s-phi_e = U
    """

    def test_concentration_stable(self):
        """I=0 多步后浓度不变。"""
        solver = make_test_solver()
        state0 = solver.initial_state(soc0=0.5)
        c_e0 = state0.c_e.copy()
        c_s0 = state0.c_s.copy()

        state = state0.copy()
        for _ in range(5):
            state, report = solver.step(state, dt=1.0, current_density=0.0)
            assert report.converged

        np.testing.assert_allclose(state.c_e, c_e0, rtol=1e-6)
        np.testing.assert_allclose(state.c_s, c_s0, rtol=1e-6)

    def test_voltage_equals_ocv(self):
        """I=0 时电压等于 OCV。"""
        solver = make_test_solver()
        state0 = solver.initial_state(soc0=0.5)
        V = solver.voltage(state0)
        U_expected = float(nmc532_ocv(0.5))
        assert abs(V - U_expected) < 0.01  # 10 mV tolerance


class TestLowRateDischarge:
    """低倍率放电测试。

    对应 TASK_PHASE0.md §5.5:
    - 极小倍率放电电压接近 OCV
    """

    def test_low_rate_voltage_near_ocv(self):
        """极小电流放电，电压应接近 OCV。"""
        solver = make_test_solver(k0=1e-8)
        state0 = solver.initial_state(soc0=0.5)
        U0 = float(nmc532_ocv(0.5))

        # 极小电流
        I_small = -0.01  # A/m²，非常小
        state, report = solver.step(state0, dt=1.0, current_density=I_small)
        assert report.converged

        V = solver.voltage(state)
        # 低倍率下电压应接近 OCV (偏差 < 50 mV)
        assert abs(V - U0) < 0.05


class TestHighRateSmoke:
    """高倍率 smoke test。

    对应 TASK_PHASE0.md §5.5:
    - 高倍率不出现负浓度或 Newton silent failure
    """

    def test_high_rate_no_negative_concentration(self):
        """高倍率不产生负浓度。"""
        solver = make_test_solver(k0=1e-8)
        state0 = solver.initial_state(soc0=0.5)

        I_high = -100.0  # A/m²
        state, report = solver.step(state0, dt=0.001, current_density=I_high)
        assert report.converged
        assert np.all(state.c_e > 0)
        assert np.all(state.c_s > 0)
        assert np.all(state.c_s < 48900.0)


class TestRunProtocol:
    """恒流放电 run 测试。"""

    def test_run_basic(self):
        """基本 run 测试。"""
        solver = make_test_solver(k0=1e-8)
        state0 = solver.initial_state(soc0=0.5)
        protocol = P2DProtocol(
            current_density=-0.1,
            t_final=10.0,
            dt_initial=1.0,
        )
        result = solver.run(state0, protocol)
        assert len(result.snapshots) >= 2
        assert result.final_state().time > 0
