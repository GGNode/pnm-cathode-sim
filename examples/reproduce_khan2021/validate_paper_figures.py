"""渲染论文图表并与验证仿真结果进行对比。

从论文 PDF 中提取 Figure 4/6/7, 在 10x10x10 合成网络上运行
0.2C/0.5C/1C/3C 放电仿真, 生成对比图并保存验证指标。
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

os.environ.setdefault("MPLCONFIGDIR", str((ROOT / "data/.matplotlib").resolve()))
os.environ.setdefault("MPLBACKEND", "Agg")

import fitz
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pnmcathode.network.generator import check_percolation, create_cathode_network
from pnmcathode.physics.ocv import nmc532_ocv
from pnmcathode.physics.separator import SeparatorParams
from pnmcathode.solver.transient import TransientSolver


PDF_PATH = ROOT / "data/khan2021_pnm_lib_cathode.pdf"
PAPER_DIR = ROOT / "data/paper_figures"
OUT_DIR = ROOT / "data/validation"

SHAPE = [10, 10, 10]
SPACING = 1e-5
A_GEOMETRIC = (10 * SPACING) ** 2
POROSITY = 0.368
CBD_FRACTION = 0.1392  # 论文 1CAL: 13.92% CBD, 49.28% NMC
THROAT_SCALE = 1.0  # 无全局喉道面积缩放
SEED = 42
T = 303.0
K0 = 1e-10
C_E0 = 1200.0
CUTOFF = 3.0
CS_MAX = 48900.0
INITIAL_SOL = 0.50
C_S0 = INITIAL_SOL * CS_MAX
NMC_DENSITY_G_M3 = 4.75e6
# 论文等效电流密度 (基于 Khan 1CAL 面积负载)
PAPER_MASS_LOADING_KG_M2 = 297.8e-3  # kg/m² (297.8 g/m²)
SPECIFIC_CAPACITY_C_KG = 178.0 * 3600.0  # 178 mAh/g → C/kg
I_1C_PAPER = PAPER_MASS_LOADING_KG_M2 * SPECIFIC_CAPACITY_C_KG / 3600.0  # A/m² ≈ 53.0

# --- C-rate / 容量基准 ---
CURRENT_BASIS = "network"       # "network" 或 "paper_areal"
CAPACITY_BASIS = "network_mass" # "network_mass" 或 "paper_mass"

# --- 隔膜默认参数 ---
SEPARATOR_ENABLED = False
SEPARATOR_PARAMS = SeparatorParams(enabled=SEPARATOR_ENABLED)

# 论文参考几何参数
PAPER_THICKNESS_UM = 75.0  # 1CAL 电极厚度 [um]


def _safe_name(crate: float) -> str:
    return f"{crate:g}".replace(".", "p")


def render_pdf_assets() -> dict[str, str]:
    PAPER_DIR.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(PDF_PATH)
    matrix = fitz.Matrix(200 / 72, 200 / 72)
    outputs: dict[str, str] = {}

    for page_no in [8, 9, 10, 11]:
        page = doc[page_no - 1]
        out = PAPER_DIR / f"page_{page_no:02d}_200dpi.png"
        page.get_pixmap(matrix=matrix, alpha=False).save(out)
        outputs[f"page_{page_no:02d}"] = str(out.relative_to(ROOT))

    # PDF 将每个目标图表嵌入为页面图像。Page 9 是 Figure 5;
    # Figure 6 和 Figure 7 在此 PDF 的第 10 和 11 页。
    figure_pages = {"figure4": 8, "figure6": 10, "figure7": 11}
    for figure, page_no in figure_pages.items():
        page = doc[page_no - 1]
        images = page.get_images(full=True)
        if not images:
            continue
        rect = page.get_image_rects(images[0][0])[0]
        out = PAPER_DIR / f"{figure}_crop_200dpi.png"
        page.get_pixmap(matrix=matrix, clip=rect, alpha=False).save(out)
        outputs[figure] = str(out.relative_to(ROOT))

    return outputs


def make_solver(separator: SeparatorParams | None = None) -> TransientSolver:
    net = create_cathode_network(
        shape=SHAPE,
        spacing=SPACING,
        porosity=POROSITY,
        cbd_fraction=CBD_FRACTION,
        seed=SEED,
        throat_scale=THROAT_SCALE,
    )
    percolation = {
        "electrolyte": check_percolation(net, phase="electrolyte"),
        "solid": check_percolation(net, phase="solid"),
    }
    for phase, ok in percolation.items():
        if not ok:
            print(
                f"WARNING: {phase} phase does not percolate "
                "from separator to collector",
                flush=True,
            )
    solver = TransientSolver(
        net, T=T, k0=K0, geometric_area=A_GEOMETRIC, separator=separator,
    )
    solver.set_concentration(c_e=C_E0, c_s=C_S0)
    return solver


def print_connectivity_summary(solver: TransientSolver) -> None:
    net = solver.net
    percolates_e = check_percolation(net, phase="electrolyte")
    percolates_s = check_percolation(net, phase="solid")
    active_e = int(np.count_nonzero(solver._steady.active_e))
    active_s = int(np.count_nonzero(solver._steady.active_s))
    reactive_interfaces = len(solver._steady.reactive_interfaces)

    print("network connectivity summary:", flush=True)
    print(f"  shape={SHAPE}, spacing={SPACING:.3g} m, geometric_area={A_GEOMETRIC:.6g} m2", flush=True)
    print(f"  pores={solver.Np}, throats={solver.Nt}", flush=True)
    print(
        f"  NMC pores={solver.n_nmc}, electrolyte pores={solver.n_e}, "
        f"solid pores={solver.n_s}",
        flush=True,
    )
    print(
        f"  percolation: electrolyte={percolates_e}, solid={percolates_s}",
        flush=True,
    )
    print(
        f"  active_e={active_e}/{solver.n_e}, active_s={active_s}/{solver.n_s}, "
        f"reactive_interfaces={reactive_interfaces}",
        flush=True,
    )


def nmc_mass_loading_g_m2(solver: TransientSolver) -> float:
    nmc_vol_m3 = float(np.sum(solver.net["pore.volume"][solver.nmc_indices]))
    return NMC_DENSITY_G_M3 * nmc_vol_m3 / A_GEOMETRIC


def _range(values: np.ndarray) -> tuple[float | None, float | None]:
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return None, None
    return float(np.min(finite)), float(np.max(finite))


def _fmt_range_mV(bounds: tuple[float | None, float | None]) -> str:
    lo, hi = bounds
    if lo is None or hi is None:
        return "n/a"
    return f"{1e3 * lo:.3f}..{1e3 * hi:.3f} mV (span={1e3 * (hi - lo):.3f} mV)"


def eta_values(solver: TransientSolver, potential: dict) -> np.ndarray:
    eta = []
    for e_g, s_g, _area in solver._steady.reactive_interfaces:
        pe = potential["phi_e"][e_g]
        ps = potential["phi_s"][s_g]
        cs = solver.c_s[s_g]
        if np.isfinite(pe) and np.isfinite(ps):
            eta.append(ps - pe - nmc532_ocv(cs / CS_MAX))
    return np.asarray(eta, dtype=float)


def print_diagnostics(
    solver: TransientSolver,
    potential: dict,
    crate: float,
    i_app: float,
    mass_loading: float,
) -> dict:
    steady_diag = solver._steady.diagnostics(potential)
    active_e = steady_diag["active_e"]
    active_e_fraction = active_e / max(solver.n_e, 1)
    phi_e_range = _range(potential["phi_e"])
    phi_s_range = _range(potential["phi_s"])
    eta_range = _range(eta_values(solver, potential))

    diagnostics = {
        "active_e_fraction": float(active_e_fraction),
        "active_e": int(active_e),
        "total_e": int(solver.n_e),
        "active_s": int(steady_diag["active_s"]),
        "reactive_interfaces": int(steady_diag["reactive_interfaces"]),
        "phi_e_range_V": phi_e_range,
        "phi_s_range_V": phi_s_range,
        "eta_range_V": eta_range,
        "mass_loading_g_m2": float(mass_loading),
        "current_density_A_m2": float(i_app),
        "geometric_area_m2": float(A_GEOMETRIC),
        "cc_area_m2": float(steady_diag["cc_area"]),
        "I_rxn": steady_diag["I_rxn"],
    }

    print(
        f"  diagnostics {crate:g}C: "
        f"active_e_fraction={active_e_fraction:.3f} ({active_e}/{solver.n_e}), "
        f"reactive_interfaces={steady_diag['reactive_interfaces']}, "
        f"phi_e_range={_fmt_range_mV(phi_e_range)}, "
        f"phi_s_range={_fmt_range_mV(phi_s_range)}, "
        f"eta_range={_fmt_range_mV(eta_range)}, "
        f"mass_loading={mass_loading:.3f} g/m2, "
        f"current_density={i_app:.6g} A/m2",
        flush=True,
    )
    if active_e_fraction < 0.5:
        print("  WARNING: LOW CONNECTIVITY active_e_fraction < 0.5", flush=True)
    return diagnostics


def add_gravimetric_capacity(
    result: dict,
    solver: TransientSolver,
    mass_loading: float | None = None,
) -> dict:
    if mass_loading is None:
        mass_loading = nmc_mass_loading_g_m2(solver)
    result = {key: np.asarray(value) if isinstance(value, list) else value for key, value in result.items()}
    result["mass_loading_g_m2"] = mass_loading
    # capacity_Ah_m2 [Ah/m²] / mass_loading_kg_m2 [kg/m²] = [Ah/kg] = [mAh/g]
    if CAPACITY_BASIS == "network_mass":
        result["capacity_mAh_g"] = result["capacity_Ah_m2"] / (mass_loading / 1000.0)
    else:
        result["capacity_mAh_g"] = result["capacity_Ah_m2"] / PAPER_MASS_LOADING_KG_M2
    return result


def accepted_step(
    solver: TransientSolver,
    dt: float,
    i_app: float,
    dt_min: float = 1e-2,
    max_retries: int = 16,
) -> tuple[dict, float]:
    """执行一个时间步, 仅在电位求解失败时用更小的 dt 重试。"""
    trial_dt = dt
    for _ in range(max_retries):
        snapshot = solver._snapshot_state()
        step = solver.step(dt=trial_dt, I_app=i_app)
        if step["converged"] and np.isfinite(step["voltage"]):
            return step, trial_dt
        sol = solver.c_s[solver.nmc_indices] / CS_MAX
        ce = solver.c_e[solver.e_indices]
        bounded_state = (
            np.isfinite(step["voltage"])
            and 0.0 < step["voltage"] < 5.0
            and float(np.min(sol)) >= 0.0
            and float(np.max(sol)) <= 1.0
            and float(np.min(ce)) >= 1.0
            and float(np.max(ce)) <= 10000.0
        )
        if bounded_state:
            return step, trial_dt
        solver._restore_state(snapshot)
        trial_dt = max(0.5 * trial_dt, dt_min)
    raise RuntimeError(f"could not find a converged validation step at I_app={i_app}")


def run_discharge(
    crate: float,
    separator: SeparatorParams | None = None,
) -> tuple[dict, TransientSolver]:
    solver = make_solver(separator=separator)
    if CURRENT_BASIS == "network":
        i_app = solver.current_density_for_c_rate(
            crate, geometric_area=A_GEOMETRIC,
        )
    else:
        i_app = -crate * I_1C_PAPER  # 论文等效电流密度
    dt = min(100.0, 50.0 / crate, solver.characteristic_dt(I_app=i_app, C_rate=crate))
    mass_loading = nmc_mass_loading_g_m2(solver)

    solver._steady.set_concentration(solver.c_e, solver.c_s)
    initial = solver._steady.solve(I_app=i_app)
    diagnostics = print_diagnostics(solver, initial, crate, i_app, mass_loading)
    times = [0.0]
    capacities = [0.0]
    voltages = [initial["voltage"]]
    dts = [0.0]
    iterations = [initial["iterations"]]
    reached_cutoff = False

    max_time = 1.25 * 3600.0 / crate
    max_steps = int(np.ceil(max_time / dt)) + 500
    # 步进前电压 (t_0 时刻, 第一步之前)
    v_prev = voltages[0]
    for _ in range(max_steps):
        step, dt_used = accepted_step(solver, dt, i_app)
        t_new = times[-1] + dt_used
        v_new = step["voltage"]  # 时间步开始时的电压 (浓度更新之前)

        # 时间步开始时的容量 (与电压时间对齐)
        q_start = abs(i_app) * times[-1]

        # 超过 max_time 时停止
        if t_new >= max_time and not reached_cutoff:
            # 记录时间终点的最终数据点
            q_end = abs(i_app) * t_new
            times.append(t_new)
            capacities.append(q_end)
            voltages.append(v_new)
            break

        if voltages[-1] > CUTOFF >= v_new:
            frac = (voltages[-1] - CUTOFF) / max(voltages[-1] - v_new, 1e-30)
            times.append(times[-1] + frac * dt_used)
            capacities.append(q_start + frac * abs(i_app) * dt_used)
            voltages.append(CUTOFF)
            dts.append(frac * dt_used)
            iterations.append(step["iterations"])
            reached_cutoff = True
            break

        # 记录: t_n 时刻的电压, t_n 时刻的容量 (步开始对齐)
        times.append(t_new)
        capacities.append(abs(i_app) * t_new)
        voltages.append(v_new)
        dts.append(dt_used)
        iterations.append(step["iterations"])

        if np.mean(solver.c_s[solver.nmc_indices] / CS_MAX) >= 0.999:
            break

    cutoff_limited = reached_cutoff
    if not cutoff_limited:
        print(
            f"  警告: {crate:g}C 容量受时间限制; "
            "运行未达到截止电压。",
            flush=True,
        )

    result = {
        "time": np.asarray(times),
        "capacity": np.asarray(capacities),
        "capacity_Ah_m2": np.asarray(capacities) / 3600.0,
        "voltage": np.asarray(voltages),
        "dt": np.asarray(dts),
        "iterations": np.asarray(iterations),
        "I_app": i_app,
        "C_rate": crate,
        "cutoff_voltage": CUTOFF,
        "reached_cutoff": reached_cutoff,
        "cutoff_limited": cutoff_limited,
        "converged_fraction": 1.0,
        "final_mean_sol": float(np.mean(solver.c_s[solver.nmc_indices] / CS_MAX)),
        "final_min_sol": float(np.min(solver.c_s[solver.nmc_indices] / CS_MAX)),
        "final_max_sol": float(np.max(solver.c_s[solver.nmc_indices] / CS_MAX)),
        "diagnostics": diagnostics,
        "current_basis": CURRENT_BASIS,
        "capacity_basis": CAPACITY_BASIS,
        "I_1C_A_m2": solver.discharge_capacity_coulombs() / 3600.0 / A_GEOMETRIC,
        "separator_enabled": solver.separator.enabled,
    }
    result = add_gravimetric_capacity(result, solver, mass_loading=mass_loading)
    out = OUT_DIR / f"discharge_{_safe_name(crate)}c.npz"
    np.savez(out, **result)
    return result, solver


def run_to_sol(crate: float, target_sol: float = 0.75) -> dict:
    solver = make_solver()
    if CURRENT_BASIS == "network":
        i_app = solver.current_density_for_c_rate(
            crate, geometric_area=A_GEOMETRIC,
        )
    else:
        i_app = -crate * I_1C_PAPER  # Paper-equivalent current density
    dt = min(100.0, 50.0 / crate, solver.characteristic_dt(I_app=i_app, C_rate=crate))

    solver._steady.set_concentration(solver.c_e, solver.c_s)
    initial = solver._steady.solve(I_app=0.0)
    previous_voltage = initial["voltage"]
    time_s = 0.0
    steps = 0
    mean_sol = float(np.mean(solver.c_s[solver.nmc_indices] / CS_MAX))

    max_steps = int(np.ceil(1.25 * 3600.0 / crate / dt)) + 500
    while mean_sol < target_sol and previous_voltage > CUTOFF and steps < max_steps:
        step, dt_used = accepted_step(solver, dt, i_app)
        time_s += dt_used
        steps += 1
        previous_voltage = step["voltage"]
        mean_sol = float(np.mean(solver.c_s[solver.nmc_indices] / CS_MAX))

    solver._steady.set_concentration(solver.c_e, solver.c_s)
    pot = solver._steady.solve(I_app=i_app)

    coords = solver.net["pore.coords"]
    result = {
        "coords": coords,
        "x": coords[:, 0],
        "e_mask": solver.e_mask,
        "nmc_mask": solver.nmc_mask,
        "solid_mask": solver.solid_mask,
        "c_e": solver.c_e.copy(),
        "phi_e": pot["phi_e"],
        "sol": np.clip(solver.c_s.copy() / CS_MAX, 0.0, 1.0),
        "phi_s": pot["phi_s"],
        "time_s": np.array(time_s),
        "steps": np.array(steps),
        "C_rate": np.array(crate),
        "target_sol": np.array(target_sol),
        "mean_sol": np.array(mean_sol),
        "voltage": np.array(pot["voltage"]),
        "reached_target": np.array(mean_sol >= target_sol),
        "converged_fraction": np.array(1.0),
    }
    out = OUT_DIR / f"spatial_{_safe_name(crate)}c_75sol.npz"
    np.savez(out, **result)
    return result


def _imread(path: Path) -> np.ndarray:
    return plt.imread(path)


def plot_vq_comparison(results: dict[float, dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.2), constrained_layout=True)
    axes[0].imshow(_imread(PAPER_DIR / "figure4_crop_200dpi.png"))
    axes[0].set_title("Paper Figure 4")
    axes[0].axis("off")

    colors = {0.2: "#1f77b4", 0.5: "#2ca02c", 1.0: "#d62728", 3.0: "#9467bd"}
    for crate, result in results.items():
        axes[1].plot(
            result["capacity_mAh_g"],
            result["voltage"],
            label=f"{crate:g}C",
            linewidth=2.0,
            color=colors.get(crate),
        )
    axes[1].axhline(CUTOFF, color="0.25", linestyle="--", linewidth=1.0, label="3.0 V cutoff")
    axes[1].set_xlabel("Capacity [mAh/g NMC, nominal density]")
    axes[1].set_ylabel("Voltage [V]")
    axes[1].set_title("This model")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(frameon=False)
    fig.savefig(OUT_DIR / "figure4_vq_comparison.png", dpi=200)
    plt.close(fig)


def _scatter_field(ax, spatial: dict, values: np.ndarray, mask: np.ndarray, title: str, cbar_label: str) -> None:
    coords = spatial["coords"]
    valid = mask & np.isfinite(values)
    sc = ax.scatter(coords[valid, 0] * 1e6, coords[valid, 1] * 1e6, c=values[valid], s=32, cmap="viridis")
    ax.set_title(title)
    ax.set_xlabel("x [um]")
    ax.set_ylabel("y [um]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.18)
    plt.colorbar(sc, ax=ax, label=cbar_label, fraction=0.046, pad=0.04)


def plot_spatial_comparison(spatial: dict, crate: float) -> None:
    fig = plt.figure(figsize=(13.0, 9.0), constrained_layout=True)
    subfigs = fig.subfigures(1, 2, width_ratios=[1.0, 1.05])
    ax_left = subfigs[0].subplots(1, 1)
    ax_left.imshow(_imread(PAPER_DIR / "figure6_crop_200dpi.png"))
    ax_left.set_title("Paper Figure 6 (3C, 75% SoL)")
    ax_left.axis("off")

    axes = subfigs[1].subplots(2, 2)
    _scatter_field(axes[0, 0], spatial, spatial["c_e"], spatial["e_mask"], "c_e", "mol/m3")
    _scatter_field(axes[0, 1], spatial, spatial["phi_e"], spatial["e_mask"], "phi_e", "V")
    _scatter_field(axes[1, 0], spatial, spatial["sol"], spatial["nmc_mask"], "SoL", "c_s/c_s,max")
    _scatter_field(axes[1, 1], spatial, spatial["phi_s"], spatial["solid_mask"], "phi_s", "V")
    subfigs[1].suptitle(f"This model ({crate:g}C, mean SoL={float(spatial['mean_sol']):.3f})")
    fig.savefig(OUT_DIR / f"figure6_spatial_comparison_{_safe_name(crate)}c.png", dpi=200)
    plt.close(fig)


def _plot_profile_axis(ax, spatial: dict, values: np.ndarray, mask: np.ndarray, label: str) -> None:
    x = spatial["x"] * 1e6
    valid = mask & np.isfinite(values)
    ax.scatter(x[valid], values[valid], s=28, alpha=0.8, label=label)
    if np.count_nonzero(valid) >= 2:
        order = np.argsort(x[valid])
        xs = x[valid][order]
        ys = values[valid][order]
        bins = np.linspace(xs.min(), xs.max(), min(8, len(xs)) + 1)
        bx, by = [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            sel = (xs >= lo) & (xs <= hi)
            if np.any(sel):
                bx.append(float(np.mean(xs[sel])))
                by.append(float(np.mean(ys[sel])))
        ax.plot(bx, by, linewidth=2.0)
    ax.set_xlabel("Distance from separator [um]")
    ax.grid(True, alpha=0.25)


def plot_profile_comparison(spatial_1c: dict, spatial_3c: dict) -> None:
    fig = plt.figure(figsize=(13.0, 8.6), constrained_layout=True)
    subfigs = fig.subfigures(1, 2, width_ratios=[1.0, 1.05])
    ax_left = subfigs[0].subplots(1, 1)
    ax_left.imshow(_imread(PAPER_DIR / "figure7_crop_200dpi.png"))
    ax_left.set_title("Paper Figure 7")
    ax_left.axis("off")

    axes = subfigs[1].subplots(3, 2)
    for col, (crate, spatial) in enumerate([(1.0, spatial_1c), (3.0, spatial_3c)]):
        axes[0, col].set_title(f"This model {crate:g}C")
        _plot_profile_axis(axes[0, col], spatial, spatial["c_e"], spatial["e_mask"], "c_e")
        axes[0, col].set_ylabel("c_e [mol/m3]")
        _plot_profile_axis(axes[1, col], spatial, spatial["phi_e"], spatial["e_mask"], "phi_e")
        axes[1, col].set_ylabel("phi_e [V]")
        _plot_profile_axis(axes[2, col], spatial, spatial["sol"], spatial["nmc_mask"], "SoL")
        axes[2, col].set_ylabel("SoL")
    fig.savefig(OUT_DIR / "figure7_profiles_comparison.png", dpi=200)
    plt.close(fig)


def end_to_end_delta(values: np.ndarray, x: np.ndarray, mask: np.ndarray) -> float | None:
    valid = mask & np.isfinite(values)
    if np.count_nonzero(valid) < 2:
        return None
    xv = x[valid]
    yv = values[valid]
    lo = xv <= np.quantile(xv, 0.25)
    hi = xv >= np.quantile(xv, 0.75)
    return float(np.mean(yv[hi]) - np.mean(yv[lo]))


def build_metrics(rendered: dict[str, str], discharges: dict[float, dict], spatials: dict[float, dict]) -> dict:
    # 几何参数披露
    thickness_modeled_um = (SHAPE[0] - 1) * SPACING * 1e6
    summary_solver = make_solver()
    model_mass_loading = nmc_mass_loading_g_m2(summary_solver)
    model_i_1c = summary_solver.discharge_capacity_coulombs() / 3600.0 / A_GEOMETRIC

    metrics: dict = {
        "rendered": rendered,
        "simulation_settings": {
            "shape": SHAPE,
            "spacing_m": SPACING,
            "porosity": POROSITY,
            "cbd_fraction": CBD_FRACTION,
            "seed": SEED,
            "cutoff_V": CUTOFF,
            "initial_sol": INITIAL_SOL,
            "nmc_density_g_m3_for_capacity": NMC_DENSITY_G_M3,
            "geometric_area_m2": A_GEOMETRIC,
            "current_basis": CURRENT_BASIS,
            "capacity_basis": CAPACITY_BASIS,
            "separator_enabled": SEPARATOR_ENABLED,
        },
        "geometry_disclosure": {
            "model_thickness_um": float(thickness_modeled_um),
            "paper_thickness_um": float(PAPER_THICKNESS_UM),
            "thickness_ratio_paper_over_model": float(PAPER_THICKNESS_UM / thickness_modeled_um),
            "model_mass_loading_g_m2": float(model_mass_loading),
            "paper_mass_loading_g_m2": float(PAPER_MASS_LOADING_KG_M2 * 1000.0),
            "mass_loading_ratio_paper_over_model": float(
                PAPER_MASS_LOADING_KG_M2 * 1000.0 / model_mass_loading
            ),
            "model_I_1C_A_m2": float(model_i_1c),
            "paper_I_1C_A_m2": float(I_1C_PAPER),
            "I_1C_ratio_paper_over_model": float(I_1C_PAPER / model_i_1c),
            "model_nodes": int(summary_solver.Np),
            "model_throats": int(summary_solver.Nt),
            "active_interfaces": int(len(summary_solver._steady.reactive_interfaces)),
            "note": (
                "10x10x10 合成网络不是 Khan 1CAL XCT 网络的定量几何匹配。"
                "除非几何和质量负载匹配, 否则对比仅为定性。"
            ),
        },
        "figure4": {},
        "spatial": {},
    }
    for crate, result in discharges.items():
        cap = result["capacity_mAh_g"]
        voltage = result["voltage"]
        metrics["figure4"][f"{crate:g}C"] = {
            "initial_voltage_V": float(voltage[0]),
            "final_voltage_V": float(voltage[-1]),
            "final_capacity_mAh_g": float(cap[-1]),
            "final_capacity_Ah_m2": float(result["capacity_Ah_m2"][-1]),
            "reached_cutoff": bool(result["reached_cutoff"]),
            "cutoff_limited": bool(result["cutoff_limited"]),
            "time_limited": not bool(result["cutoff_limited"]),
            "converged_fraction": float(result["converged_fraction"]),
            "final_mean_sol": float(result["final_mean_sol"]),
            "final_min_sol": float(result["final_min_sol"]),
            "final_max_sol": float(result["final_max_sol"]),
            "diagnostics": result["diagnostics"],
            "n_points": int(len(voltage)),
        }

    for crate, spatial in spatials.items():
        x = spatial["x"]
        metrics["spatial"][f"{crate:g}C"] = {
            "mean_sol": float(spatial["mean_sol"]),
            "voltage_V": float(spatial["voltage"]),
            "time_s": float(spatial["time_s"]),
            "reached_target": bool(spatial["reached_target"]),
            "converged_fraction": float(spatial["converged_fraction"]),
            "ce_delta_cc_minus_sep_mol_m3": end_to_end_delta(spatial["c_e"], x, spatial["e_mask"]),
            "phie_delta_cc_minus_sep_V": end_to_end_delta(spatial["phi_e"], x, spatial["e_mask"]),
            "sol_delta_cc_minus_sep": end_to_end_delta(spatial["sol"], x, spatial["nmc_mask"]),
            "phis_delta_cc_minus_sep_V": end_to_end_delta(spatial["phi_s"], x, spatial["solid_mask"]),
        }
    return metrics


def main() -> None:
    PAPER_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rendered = render_pdf_assets()
    summary_solver = make_solver(separator=SEPARATOR_PARAMS)
    print_connectivity_summary(summary_solver)

    # --- 几何参数披露 ---
    thickness_modeled_um = (SHAPE[0] - 1) * SPACING * 1e6
    model_mass = nmc_mass_loading_g_m2(summary_solver)
    paper_mass = PAPER_MASS_LOADING_KG_M2 * 1000.0
    model_i1c = summary_solver.discharge_capacity_coulombs() / 3600.0 / A_GEOMETRIC
    print("\ngeometry disclosure:", flush=True)
    print(f"  model thickness: {thickness_modeled_um:.1f} um "
          f"(paper 1CAL: {PAPER_THICKNESS_UM:.1f} um, "
          f"ratio={PAPER_THICKNESS_UM / thickness_modeled_um:.2f})", flush=True)
    print(f"  model mass loading: {model_mass:.2f} g/m2 "
          f"(paper: {paper_mass:.2f} g/m2, "
          f"ratio={paper_mass / model_mass:.2f})", flush=True)
    print(f"  model I_1C: {model_i1c:.4f} A/m2 "
          f"(paper: {I_1C_PAPER:.4f} A/m2)", flush=True)
    print(f"  current_basis={CURRENT_BASIS}, capacity_basis={CAPACITY_BASIS}", flush=True)
    print(f"  separator_enabled={SEPARATOR_ENABLED}", flush=True)

    discharges: dict[float, dict] = {}
    for crate in [0.2, 0.5, 1.0, 3.0]:
        print(f"\nrunning {crate:g}C discharge (basis={CURRENT_BASIS})", flush=True)
        result, _solver = run_discharge(crate, separator=SEPARATOR_PARAMS)
        discharges[crate] = result

    plot_vq_comparison(discharges)

    spatials: dict[float, dict] = {}
    for crate in [1.0, 3.0]:
        print(f"running {crate:g}C spatial snapshot", flush=True)
        spatial = run_to_sol(crate, target_sol=0.75)
        spatials[crate] = spatial
        plot_spatial_comparison(spatial, crate)

    plot_profile_comparison(spatials[1.0], spatials[3.0])

    metrics = build_metrics(rendered, discharges, spatials)
    metrics_path = OUT_DIR / "validation_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(f"saved metrics: {metrics_path}", flush=True)

    # --- 4 场景汇总 ---
    print_4scenario_summary()


def print_4scenario_summary() -> None:
    """运行 4 种场景并打印对比表。"""
    scenarios = [
        ("network, sep OFF", "network", "network_mass", SeparatorParams(enabled=False)),
        ("network, sep ON",  "network", "network_mass", SeparatorParams(enabled=True)),
        ("paper,   sep OFF", "paper_areal", "paper_mass", SeparatorParams(enabled=False)),
        ("paper,   sep ON",  "paper_areal", "paper_mass", SeparatorParams(enabled=True)),
    ]
    print("\n" + "=" * 80, flush=True)
    print("4-SCENARIO SUMMARY", flush=True)
    print("=" * 80, flush=True)

    for label, cur_basis, cap_basis, sep_params in scenarios:
        # 临时覆盖全局变量
        global CURRENT_BASIS, CAPACITY_BASIS
        CURRENT_BASIS = cur_basis
        CAPACITY_BASIS = cap_basis

        print(f"\n--- {label} ---", flush=True)
        for crate in [0.2, 3.0]:
            result, solver = run_discharge(crate, separator=sep_params)
            v0 = result["voltage"][0]
            vf = result["voltage"][-1]
            cap = result["capacity_mAh_g"][-1]
            cutoff = result["reached_cutoff"]
            mass = result["mass_loading_g_m2"]
            i1c = result["I_1C_A_m2"]
            i_app = result["I_app"]

            # 隔膜损失
            if sep_params.enabled:
                sep_state = solver._steady.separator_state
                sep_info = (
                    f"  sep ohmic={sep_state.dphi_ohm * 1e3:.2f} mV, "
                    f"conc={sep_state.dphi_conc * 1e3:.2f} mV, "
                    f"Li BV={sep_state.eta_li * 1e3:.2f} mV"
                )
            else:
                sep_info = "  (无隔膜模型)"

            # 从诊断信息获取 phi_e 跨度
            diag = result.get("diagnostics", {})
            phi_e_range = diag.get("phi_e_range_V", (None, None))
            if phi_e_range[0] is not None and phi_e_range[1] is not None:
                phi_e_span = phi_e_range[1] - phi_e_range[0]
            else:
                phi_e_span = 0.0

            print(
                f"  {crate:g}C: V0={v0:.4f}V, Vf={vf:.4f}V, "
                f"cap={cap:.1f} mAh/g, cutoff={cutoff}, "
                f"I_app={i_app:.4f} A/m2, I_1C={i1c:.4f} A/m2, "
                f"mass={mass:.2f} g/m2",
                flush=True,
            )
            print(
                f"  phi_e span={phi_e_span * 1e3:.2f} mV{sep_info}",
                flush=True,
            )

    # 恢复默认值
    CURRENT_BASIS = "network"
    CAPACITY_BASIS = "network_mass"


if __name__ == "__main__":
    main()
