"""
Phase 1: 浓度分布形状测试 (Spatial Profiles)
=============================================

验证:
- 放电后 mean(c̄_s) > mean(c̄_s(0))。
- 高倍率下 c_e 梯度显著。
- c_e 保持正值。
- c_s 保持物理范围。
- 高倍率下 c_s_surf - c̄_s 幅值大于低倍率。

对应 TASK_PHASE1.md §3.7
"""

import numpy as np
import pytest

from pnmcathode.p2d.benchmark import (
    P2DBenchmarkCase,
    P2DProfile,
    extract_profile,
    make_khan_solver,
)
from pnmcathode.p2d.solver import P2DProtocol
from pnmcathode.p2d.residual import P2DResidualContext, compute_fluxes


def _run_to_sol(target_sol: float, c_rate: float, n_cells: int = 12):
    """运行放电到目标 SoL 或截止。"""
    from pnmcathode.materials.presets import (
        electrolyte_khan2021,
        nmc532_khan2021,
        separator_khan2021,
    )
    from pnmcathode.config import Kinetics, SolverSettings
    from pnmcathode.p2d.domain import P2DRegion
    from pnmcathode.p2d.materials import from_config
    from pnmcathode.p2d.benchmark import c_rate_to_current_density
    from pnmcathode.p2d.solver import P2DSolver

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
        particle_shells=6,
        n_separator_cells=4,
    )
    solver = P2DSolver(params, settings)

    state0 = solver.initial_state(soc0=0.35)
    I_app = c_rate_to_current_density(
        c_rate, pos.epsilon_s, pos.length, active.cs_max
    )

    # 计算目标时间: target_sol * F * eps_s * L * cs_max / |I_app|
    from pnmcathode.physics.reaction import F
    t_target = target_sol * F * pos.epsilon_s * pos.length * active.cs_max / abs(I_app)
    t_final = min(t_target * 1.2, 500.0)

    protocol = P2DProtocol(
        current_density=I_app,
        t_final=t_final,
        dt_initial=5.0,
        cutoff_voltage=3.0,
        save_every=1,
    )
    result = solver.run(state0, protocol)
    return solver, result


class TestSolidConcentrationIncrease:
    """放电后固相平均浓度增加。"""

    def test_mean_cs_increases(self):
        """放电后 mean(c̄_s) > mean(c̄_s(0))。"""
        solver, result = _run_to_sol(0.1, 1.0)
        profile0 = extract_profile(
            result.snapshots[0], solver.macro, solver.particle,
            solver.params.material.active.cs_max,
        )
        profile_f = extract_profile(
            result.snapshots[-1], solver.macro, solver.particle,
            solver.params.material.active.cs_max,
        )
        assert np.mean(profile_f.c_s_average) > np.mean(profile0.c_s_average), (
            "放电后 mean(c̄_s) 未增加"
        )


class TestElectrolyteGradient:
    """高倍率下电解质浓度梯度。"""

    def test_1c_gradient_significant(self):
        """1C 放电后 c_e 梯度 >= 0.05 * c_e0。"""
        solver, result = _run_to_sol(0.15, 1.0)
        profile = extract_profile(
            result.snapshots[-1], solver.macro, solver.particle,
            solver.params.material.active.cs_max,
        )
        c_e0 = solver.params.material.electrolyte.c_init
        # 只检查 separator + positive 中的 c_e
        ce_range = np.max(profile.c_e) - np.min(profile.c_e)
        assert ce_range >= 0.05 * c_e0, (
            f"1C: c_e 梯度 {ce_range:.1f} mol/m³ < 0.05 * c_e0 = {0.05*c_e0:.1f}"
        )

    def test_3c_gradient_larger_than_1c(self):
        """3C 梯度应大于 1C 梯度。"""
        solver1, result1 = _run_to_sol(0.1, 1.0, n_cells=15)
        solver3, result3 = _run_to_sol(0.1, 3.0, n_cells=15)

        p1 = extract_profile(
            result1.snapshots[-1], solver1.macro, solver1.particle,
            solver1.params.material.active.cs_max,
        )
        p3 = extract_profile(
            result3.snapshots[-1], solver3.macro, solver3.particle,
            solver3.params.material.active.cs_max,
        )
        grad1 = np.max(p1.c_e) - np.min(p1.c_e)
        grad3 = np.max(p3.c_e) - np.min(p3.c_e)
        assert grad3 > grad1, (
            f"3C 梯度 {grad3:.1f} <= 1C 梯度 {grad1:.1f}"
        )


class TestConcentrationPhysicalRange:
    """浓度保持物理范围。"""

    def test_ce_positive(self):
        """c_e 保持正值。"""
        solver, result = _run_to_sol(0.2, 3.0)
        for snap in result.snapshots:
            assert np.min(snap.state.c_e) >= 1.0, (
                f"t={snap.time:.1f}s: min(c_e) = {np.min(snap.state.c_e):.2f} < 1.0 mol/m³"
            )

    def test_cs_in_range(self):
        """0 < c_s < c_s_max。"""
        solver, result = _run_to_sol(0.2, 1.0)
        cs_max = solver.params.material.active.cs_max
        for snap in result.snapshots:
            assert np.all(snap.state.c_s > 0), (
                f"t={snap.time:.1f}s: 存在 c_s <= 0"
            )
            assert np.all(snap.state.c_s < cs_max), (
                f"t={snap.time:.1f}s: 存在 c_s >= c_s_max"
            )


class TestSurfaceAverageConcentrationGap:
    """高倍率下表面-平均浓度差大于低倍率。"""

    def test_3c_gap_larger_than_02c(self):
        """max|c_s_surf - c̄_s|_3C > max|c_s_surf - c̄_s|_0.2C。"""
        solver02, result02 = _run_to_sol(0.1, 0.2, n_cells=15)
        solver3, result3 = _run_to_sol(0.1, 3.0, n_cells=15)

        cs_max02 = solver02.params.material.active.cs_max
        cs_max3 = solver3.params.material.active.cs_max

        gap02 = 0.0
        for snap in result02.snapshots:
            surf = snap.state.surface_concentration()
            avg = snap.state.average_solid_concentration(solver02.particle)
            gap02 = max(gap02, float(np.max(np.abs(surf - avg))))

        gap3 = 0.0
        for snap in result3.snapshots:
            surf = snap.state.surface_concentration()
            avg = snap.state.average_solid_concentration(solver3.particle)
            gap3 = max(gap3, float(np.max(np.abs(surf - avg))))

        assert gap3 > gap02, (
            f"3C gap {gap3:.1f} <= 0.2C gap {gap02:.1f} mol/m³"
        )
