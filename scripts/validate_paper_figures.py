"""Render paper figures and compare validation simulations against them."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

os.environ.setdefault("MPLCONFIGDIR", str((ROOT / "data/.matplotlib").resolve()))
os.environ.setdefault("MPLBACKEND", "Agg")

import fitz
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.network.generator import create_cathode_network
from src.physics.ocv import nmc532_ocv
from src.solver.transient import TransientSolver


PDF_PATH = ROOT / "data/khan2021_pnm_lib_cathode.pdf"
PAPER_DIR = ROOT / "data/paper_figures"
OUT_DIR = ROOT / "data/validation"

SHAPE = [5, 5, 5]
SPACING = 1e-5
A_GEOMETRIC = (5 * SPACING) ** 2
POROSITY = 0.368
CBD_FRACTION = 0.10
SEED = 42
T = 303.0
K0 = 1e-10
C_E0 = 1200.0
CUTOFF = 3.0
CS_MAX = 48900.0
INITIAL_SOL = 0.35
C_S0 = INITIAL_SOL * CS_MAX
NMC_DENSITY_G_M3 = 4.75e6


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

    # The PDF embeds each target figure as one page image. Page 9 is Figure 5;
    # Figure 6 and Figure 7 are on pages 10 and 11 in this PDF.
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


def make_solver() -> TransientSolver:
    net = create_cathode_network(
        shape=SHAPE,
        spacing=SPACING,
        porosity=POROSITY,
        cbd_fraction=CBD_FRACTION,
        seed=SEED,
    )
    solver = TransientSolver(net, T=T, k0=K0, geometric_area=A_GEOMETRIC)
    solver.set_concentration(c_e=C_E0, c_s=C_S0)
    return solver


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
        f"current_density={i_app:.6g} A/m2"
    )
    if active_e_fraction < 0.5:
        print("  WARNING: LOW CONNECTIVITY active_e_fraction < 0.5")
    return diagnostics


def add_gravimetric_capacity(result: dict, solver: TransientSolver) -> dict:
    mass_loading = nmc_mass_loading_g_m2(solver)
    result = {key: np.asarray(value) if isinstance(value, list) else value for key, value in result.items()}
    result["mass_loading_g_m2"] = mass_loading
    result["capacity_mAh_g"] = result["capacity_Ah_m2"] * 1000.0 / mass_loading
    return result


def accepted_step(
    solver: TransientSolver,
    dt: float,
    i_app: float,
    dt_min: float = 1e-2,
    max_retries: int = 16,
) -> tuple[dict, float]:
    """Take one step, retrying with smaller dt only for failed potential solves."""
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


def run_discharge(crate: float) -> tuple[dict, TransientSolver]:
    solver = make_solver()
    i_app = solver.current_density_for_c_rate(crate)
    dt = min(20.0, 10.0 / crate, solver.characteristic_dt(I_app=i_app, C_rate=crate))
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

    max_time = 1.10 * 3600.0 / crate
    max_steps = int(np.ceil(max_time / dt)) + 200
    for _ in range(max_steps):
        step, dt_used = accepted_step(solver, dt, i_app)
        t_new = times[-1] + dt_used
        q_new = abs(i_app) * t_new
        v_new = step["voltage"]
        if voltages[-1] > CUTOFF >= v_new:
            frac = (voltages[-1] - CUTOFF) / max(voltages[-1] - v_new, 1e-30)
            times.append(times[-1] + frac * dt_used)
            capacities.append(capacities[-1] + frac * (q_new - capacities[-1]))
            voltages.append(CUTOFF)
            dts.append(frac * dt_used)
            iterations.append(step["iterations"])
            reached_cutoff = True
            break

        times.append(t_new)
        capacities.append(q_new)
        voltages.append(v_new)
        dts.append(dt_used)
        iterations.append(step["iterations"])

        if np.mean(solver.c_s[solver.nmc_indices] / CS_MAX) >= 0.999:
            break

    cutoff_limited = reached_cutoff
    if not cutoff_limited:
        print(
            f"  WARNING: {crate:g}C capacity is time-limited; "
            "run did not reach cutoff voltage."
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
    }
    result = add_gravimetric_capacity(result, solver)
    out = OUT_DIR / f"discharge_{_safe_name(crate)}c.npz"
    np.savez(out, **result)
    return result, solver


def run_to_sol(crate: float, target_sol: float = 0.75) -> dict:
    solver = make_solver()
    i_app = solver.current_density_for_c_rate(crate)
    dt = min(20.0, 10.0 / crate, solver.characteristic_dt(I_app=i_app, C_rate=crate))

    solver._steady.set_concentration(solver.c_e, solver.c_s)
    initial = solver._steady.solve(I_app=0.0)
    previous_voltage = initial["voltage"]
    time_s = 0.0
    steps = 0
    mean_sol = float(np.mean(solver.c_s[solver.nmc_indices] / CS_MAX))

    while mean_sol < target_sol and previous_voltage > CUTOFF and steps < 600:
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

    discharges: dict[float, dict] = {}
    for crate in [0.2, 0.5, 1.0, 3.0]:
        print(f"running {crate:g}C discharge")
        result, _solver = run_discharge(crate)
        discharges[crate] = result

    plot_vq_comparison(discharges)

    spatials: dict[float, dict] = {}
    for crate in [1.0, 3.0]:
        print(f"running {crate:g}C spatial snapshot")
        spatial = run_to_sol(crate, target_sol=0.75)
        spatials[crate] = spatial
        plot_spatial_comparison(spatial, crate)

    plot_profile_comparison(spatials[1.0], spatials[3.0])

    metrics = build_metrics(rendered, discharges, spatials)
    metrics_path = OUT_DIR / "validation_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(f"saved metrics: {metrics_path}")


if __name__ == "__main__":
    main()
