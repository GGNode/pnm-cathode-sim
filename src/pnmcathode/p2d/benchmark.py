"""
P2D 基准测试工具 (Benchmark Utilities)
=======================================

职责:
- 定义 benchmark case 数据结构与参数快照。
- C-rate → 电流密度换算。
- 构建 Khan 2021 参数的 P2DSolver。
- 提取电压曲线、浓度/电位分布。
- 保存可复现 benchmark 数据。
- 参数源隔离断言。

对应 TASK_PHASE1.md §4, §6
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from pnmcathode.config import Kinetics, SolverSettings
from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh
from pnmcathode.p2d.materials import P2DParameters, from_config
from pnmcathode.p2d.results import P2DResult, P2DSnapshot
from pnmcathode.p2d.solver import P2DProtocol, P2DSolver
from pnmcathode.p2d.state import P2DState
from pnmcathode.physics.reaction import F

Array = np.ndarray


# ===== 数据结构 =====


@dataclass(frozen=True)
class P2DBenchmarkCase:
    """Benchmark 算例定义。

    Attributes
    ----------
    name : str
        算例名称 (如 "khan_1cal_1c")。
    parameter_source : str
        参数来源 ("khan2021" 或 "pybamm")。
    electrode_kind : str
        电极类型 ("1CAL", "3CAL", "pybamm-default")。
    c_rate : float
        C 倍率。
    initial_soc : float
        初始 SOC。
    temperature : float
        温度 [K]。
    cutoff_voltage : float
        截止电压 [V]。
    n_positive_cells : int
        positive electrode 网格单元数。
    n_separator_cells : int
        separator 网格单元数。
    n_particle_shells : int
        颗粒径向壳层数。
    dt_initial : float
        初始时间步 [s]。
    t_final : float
        终止时间 [s]。
    expected_voltage_path : Path or None
        参考电压曲线文件路径。
    """

    name: str
    parameter_source: Literal["khan2021", "pybamm"]
    electrode_kind: Literal["1CAL", "3CAL", "pybamm-default"]
    c_rate: float
    initial_soc: float
    temperature: float
    cutoff_voltage: float
    n_positive_cells: int
    n_separator_cells: int
    n_particle_shells: int
    dt_initial: float
    t_final: float
    expected_voltage_path: Path | None = None


@dataclass(frozen=True)
class P2DParameterSnapshot:
    """参数快照，用于验证参数源一致性。

    Attributes
    ----------
    source : str
        参数来源 ("khan2021" 或 "pybamm")。
    temperature : float
        温度 [K]。
    separator_thickness : float
        隔膜厚度 [m]。
    separator_porosity : float
        隔膜孔隙率。
    positive_thickness : float
        positive electrode 厚度 [m]。
    positive_porosity : float
        positive electrode 孔隙率。
    positive_active_fraction : float
        活性材料体积分数。
    particle_radius : float
        颗粒半径 [m]。
    c_e0 : float
        初始电解质浓度 [mol/m³]。
    c_s_max : float
        固相最大浓度 [mol/m³]。
    k0 : float
        反应速率常数。
    t_plus : float
        Li+ 迁移数。
    sigma_s : float
        固相电子电导率 [S/m]。
    bruggeman_e : float
        电解质 Bruggeman 指数。
    bruggeman_s : float
        固相 Bruggeman 指数。
    ocv_reference : str
        OCV 模型参考标识。
    """

    source: str
    temperature: float
    separator_thickness: float
    separator_porosity: float
    positive_thickness: float
    positive_porosity: float
    positive_active_fraction: float
    particle_radius: float
    c_e0: float
    c_s_max: float
    k0: float
    t_plus: float
    sigma_s: float
    bruggeman_e: float
    bruggeman_s: float
    ocv_reference: str


@dataclass(frozen=True)
class P2DProfile:
    """空间分布数据。

    Attributes
    ----------
    x_macro : ndarray
        宏观网格中心坐标 [m]。
    x_positive : ndarray
        positive electrode 网格中心坐标 [m]。
    c_e : ndarray
        电解质浓度 [mol/m³]。
    phi_e : ndarray
        电解质电位 [V]。
    phi_s : ndarray
        固相电位 [V]。
    c_s_surface : ndarray
        表面固相浓度 [mol/m³]。
    c_s_average : ndarray
        平均固相浓度 [mol/m³]。
    sol_surface : ndarray
        表面 SoL (c_s_surf / c_s_max)。
    sol_average : ndarray
        平均 SoL。
    """

    x_macro: Array
    x_positive: Array
    c_e: Array
    phi_e: Array
    phi_s: Array
    c_s_surface: Array
    c_s_average: Array
    sol_surface: Array
    sol_average: Array


# ===== 工具函数 =====


def c_rate_to_current_density(
    c_rate: float,
    epsilon_s: float,
    thickness: float,
    cs_max: float,
    soc_window: float = 1.0,
) -> float:
    """C-rate → 电流密度换算 (阳极约定, 放电 < 0)。

    I_1C = F * epsilon_s * L_p * c_s_max * Delta_SOC / 3600

    Parameters
    ----------
    c_rate : float
        C 倍率。
    epsilon_s : float
        活性材料体积分数。
    thickness : float
        positive electrode 厚度 [m]。
    cs_max : float
        固相最大浓度 [mol/m³]。
    soc_window : float
        SOC 窗口，默认 1.0。

    Returns
    -------
    float
        电流密度 [A/m²] (放电为负)。
    """
    i_1c = F * epsilon_s * thickness * cs_max * soc_window / 3600.0
    return -c_rate * i_1c


def make_khan_solver(case: P2DBenchmarkCase) -> P2DSolver:
    """从 benchmark case 构建 Khan 2021 参数的 P2DSolver。

    Parameters
    ----------
    case : P2DBenchmarkCase
        基准算例定义。

    Returns
    -------
    P2DSolver
    """
    from pnmcathode.materials.presets import (
        electrolyte_khan2021,
        nmc532_khan2021,
        separator_khan2021,
    )

    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    separator = separator_khan2021(enabled=True)
    settings = SolverSettings(
        temperature=case.temperature,
        newton_tol=1e-8,
        dt_max=max(case.dt_initial, 1.0),
    )
    kinetics = Kinetics(k0=1e-10, alpha_a=0.5, alpha_c=0.5)

    if case.electrode_kind == "1CAL":
        thickness = 129e-6
        epsilon_e = 0.368
        epsilon_s = 0.4928
    elif case.electrode_kind == "3CAL":
        thickness = 88e-6
        epsilon_e = 0.366
        epsilon_s = 0.4974
    else:
        raise ValueError(f"不支持的 Khan 电极类型: {case.electrode_kind}")

    positive = P2DRegion(
        name="positive",
        x_left=separator.thickness,
        x_right=separator.thickness + thickness,
        n_cells=case.n_positive_cells,
        epsilon_e=epsilon_e,
        epsilon_s=epsilon_s,
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
        positive=positive,
        particle_shells=case.n_particle_shells,
        n_separator_cells=case.n_separator_cells,
    )
    return P2DSolver(params=params, settings=settings)


def run_p2d_benchmark(case: P2DBenchmarkCase) -> P2DResult:
    """运行 P2D benchmark 算例。

    Parameters
    ----------
    case : P2DBenchmarkCase
        基准算例定义。

    Returns
    -------
    P2DResult
    """
    solver = make_khan_solver(case)
    state0 = solver.initial_state(soc0=case.initial_soc)
    current_density = c_rate_to_current_density(
        c_rate=case.c_rate,
        epsilon_s=solver.params.positive.epsilon_s,
        thickness=solver.params.positive.length,
        cs_max=solver.params.material.active.cs_max,
    )
    protocol = P2DProtocol(
        current_density=current_density,
        t_final=case.t_final,
        dt_initial=case.dt_initial,
        cutoff_voltage=case.cutoff_voltage,
        save_every=1,
    )
    return solver.run(state0, protocol)


def make_parameter_snapshot(
    solver: P2DSolver,
    source: str = "khan2021",
) -> P2DParameterSnapshot:
    """从 P2DSolver 提取参数快照。

    Parameters
    ----------
    solver : P2DSolver
        P2D 求解器。
    source : str
        参数来源标识。

    Returns
    -------
    P2DParameterSnapshot
    """
    params = solver.params
    mat = params.material
    pos = params.positive
    sep = params.separator_region

    return P2DParameterSnapshot(
        source=source,
        temperature=mat.temperature,
        separator_thickness=sep.x_right - sep.x_left,
        separator_porosity=sep.epsilon_e,
        positive_thickness=pos.x_right - pos.x_left,
        positive_porosity=pos.epsilon_e,
        positive_active_fraction=pos.epsilon_s,
        particle_radius=pos.particle_radius,
        c_e0=mat.electrolyte.c_init,
        c_s_max=mat.active.cs_max,
        k0=mat.kinetics.k0,
        t_plus=mat.separator.t_plus,
        sigma_s=pos.sigma_s if pos.sigma_s is not None else 0.0,
        bruggeman_e=pos.bruggeman_e,
        bruggeman_s=pos.bruggeman_s,
        ocv_reference=mat.ocv_model.reference,
    )


def extract_voltage_curve(
    result: P2DResult,
) -> tuple[Array, Array, Array]:
    """提取电压曲线。

    Parameters
    ----------
    result : P2DResult
        仿真结果。

    Returns
    -------
    time_s : ndarray
        时间 [s]。
    capacity_ah_m2 : ndarray
        容量 [A·h/m²]。
    voltage_v : ndarray
        电压 [V]。
    """
    return result.time, result.capacity_Ah_m2, result.voltage


def compare_voltage_curve(
    candidate_capacity: Array,
    candidate_voltage: Array,
    reference_capacity: Array,
    reference_voltage: Array,
) -> dict[str, float]:
    """比较两条电压曲线。

    Parameters
    ----------
    candidate_capacity : ndarray
        候选曲线容量。
    candidate_voltage : ndarray
        候选曲线电压。
    reference_capacity : ndarray
        参考曲线容量。
    reference_voltage : ndarray
        参考曲线电压。

    Returns
    -------
    dict
        rms_v, max_abs_v, bias_v, q_min, q_max
    """
    q_min = max(float(candidate_capacity.min()), float(reference_capacity.min()))
    q_max = min(float(candidate_capacity.max()), float(reference_capacity.max()))
    if q_max <= q_min:
        return {
            "rms_v": float("nan"),
            "max_abs_v": float("nan"),
            "bias_v": float("nan"),
            "q_min_Ah_m2": q_min,
            "q_max_Ah_m2": q_max,
        }
    q_grid = np.linspace(q_min, q_max, 200)
    v_candidate = np.interp(q_grid, candidate_capacity, candidate_voltage)
    v_reference = np.interp(q_grid, reference_capacity, reference_voltage)
    err = v_candidate - v_reference
    return {
        "rms_v": float(np.sqrt(np.mean(err**2))),
        "max_abs_v": float(np.max(np.abs(err))),
        "bias_v": float(np.mean(err)),
        "q_min_Ah_m2": q_min,
        "q_max_Ah_m2": q_max,
    }


def extract_profile(
    snapshot: P2DSnapshot,
    macro: MacroMesh,
    particle: ParticleMesh,
    cs_max: float,
) -> P2DProfile:
    """从快照提取空间分布。

    Parameters
    ----------
    snapshot : P2DSnapshot
        时间步快照。
    macro : MacroMesh
        宏观网格。
    particle : ParticleMesh
        颗粒网格。
    cs_max : float
        固相最大浓度 [mol/m³]。

    Returns
    -------
    P2DProfile
    """
    state = snapshot.state
    pos = macro.positive_cells
    c_s_surface = state.surface_concentration()
    c_s_average = state.average_solid_concentration(particle)
    return P2DProfile(
        x_macro=macro.x_centers.copy(),
        x_positive=macro.x_centers[pos].copy(),
        c_e=state.c_e.copy(),
        phi_e=state.phi_e.copy(),
        phi_s=state.phi_s.copy(),
        c_s_surface=c_s_surface,
        c_s_average=c_s_average,
        sol_surface=c_s_surface / cs_max,
        sol_average=c_s_average / cs_max,
    )


def save_benchmark_npz(
    path: Path | str,
    result: P2DResult,
    solver: P2DSolver,
) -> None:
    """保存 benchmark 结果到 .npz 文件。

    Parameters
    ----------
    path : Path or str
        输出文件路径。
    result : P2DResult
        仿真结果。
    solver : P2DSolver
        求解器 (提供网格信息)。
    """
    final = result.final_state()
    np.savez(
        path,
        time=result.time,
        voltage=result.voltage,
        capacity_Ah_m2=result.capacity_Ah_m2,
        x=solver.macro.x_centers,
        x_positive=solver.macro.x_centers[solver.macro.positive_cells],
        c_e_final=final.c_e,
        phi_e_final=final.phi_e,
        phi_s_final=final.phi_s,
        c_s_surface_final=final.surface_concentration(),
        c_s_average_final=final.average_solid_concentration(solver.particle),
        current_density=result.snapshots[0].current_density,
    )


def assert_no_mixed_parameter_sources(snapshot: P2DParameterSnapshot) -> None:
    """断言参数源未混用。

    Parameters
    ----------
    snapshot : P2DParameterSnapshot
        参数快照。

    Raises
    ------
    AssertionError
        如果参数源混用。
    """
    if snapshot.source == "khan2021":
        assert snapshot.temperature == 303.0, (
            f"Khan 2021 温度应为 303.0 K, 实际 {snapshot.temperature}"
        )
        assert snapshot.c_s_max == 48900.0, (
            f"Khan 2021 c_s_max 应为 48900.0, 实际 {snapshot.c_s_max}"
        )
        assert snapshot.k0 == 1e-10, (
            f"Khan 2021 k0 应为 1e-10, 实际 {snapshot.k0}"
        )
        assert snapshot.ocv_reference == "li_metal_empirical", (
            f"Khan 2021 OCV reference 应为 'li_metal_empirical', "
            f"实际 '{snapshot.ocv_reference}'"
        )
    elif snapshot.source == "pybamm":
        assert snapshot.ocv_reference.startswith("pybamm:"), (
            f"PyBaMM OCV reference 应以 'pybamm:' 开头, "
            f"实际 '{snapshot.ocv_reference}'"
        )
    else:
        raise AssertionError(f"未知参数源: {snapshot.source}")
