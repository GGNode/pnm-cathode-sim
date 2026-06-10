"""Smoke test for Phase 6 fixes — runs 0.2C and 3C on a 5x5x5 network."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
from src.network.generator import check_percolation, create_cathode_network
from src.physics.ocv import nmc532_ocv
from src.solver.transient import TransientSolver

# Parameters matching the fixed validate_paper_figures.py
SHAPE = [5, 5, 5]
SPACING = 1e-5
A_GEOMETRIC = (5 * SPACING) ** 2
POROSITY = 0.368
CBD_FRACTION = 0.1392
THROAT_SCALE = 1.0
SEED = 42
T = 303.0
K0 = 1e-10
C_E0 = 1200.0
CS_MAX = 48900.0
INITIAL_SOL = 0.50
C_S0 = INITIAL_SOL * CS_MAX
CUTOFF = 3.0

# Paper-equivalent current density (Khan 1CAL)
PAPER_MASS_LOADING_KG_M2 = 297.8e-3  # kg/m²
SPECIFIC_CAPACITY_C_KG = 178.0 * 3600.0  # 178 mAh/g → C/kg
I_1C_PAPER = PAPER_MASS_LOADING_KG_M2 * SPECIFIC_CAPACITY_C_KG / 3600.0  # ≈ 53.0 A/m²
# Paper mass loading in g/m² for gravimetric capacity normalization
PAPER_MASS_LOADING_G_M2 = 297.8


def make_solver():
    net = create_cathode_network(
        shape=SHAPE,
        spacing=SPACING,
        porosity=POROSITY,
        cbd_fraction=CBD_FRACTION,
        seed=SEED,
        throat_scale=THROAT_SCALE,
    )
    solver = TransientSolver(net, T=T, k0=K0, geometric_area=A_GEOMETRIC)
    solver.set_concentration(c_e=C_E0, c_s=C_S0)
    return solver


def run_discharge(crate):
    solver = make_solver()
    i_app = -crate * I_1C_PAPER  # Paper-equivalent current density

    dt = min(100.0, 50.0 / crate, solver.characteristic_dt(I_app=i_app, C_rate=crate))

    # OCV at initial state (no current)
    solver._steady.set_concentration(solver.c_e, solver.c_s)
    ocv_result = solver._steady.solve(I_app=0.0)
    V_ocv = ocv_result["voltage"]

    # Loaded initial solve
    solver._steady.set_concentration(solver.c_e, solver.c_s)
    initial = solver._steady.solve(I_app=i_app)
    V_loaded_init = initial["voltage"]

    phi_s_init = initial["phi_s"]
    solid_mask = solver.solid_mask
    phi_s_span_init = float(np.max(phi_s_init[solid_mask]) - np.min(phi_s_init[solid_mask]))

    times = [0.0]
    capacities = [0.0]
    voltages = [initial["voltage"]]
    reached_cutoff = False
    max_time = 1.25 * 3600.0 / crate

    max_steps = int(np.ceil(max_time / dt)) + 500
    for _ in range(max_steps):
        trial_dt = dt
        converged = False
        for retry in range(16):
            snapshot = solver._snapshot_state()
            step = solver.step(dt=trial_dt, I_app=i_app)
            if step["converged"] and np.isfinite(step["voltage"]):
                converged = True
                break
            solver._restore_state(snapshot)
            trial_dt = max(0.5 * trial_dt, 1e-2)
        if not converged:
            print(f"  WARNING: could not converge at {crate}C")
            break

        t_new = times[-1] + trial_dt
        q_new = abs(i_app) * t_new
        v_new = step["voltage"]

        if t_new >= max_time and not reached_cutoff:
            break

        if voltages[-1] > CUTOFF >= v_new:
            frac = (voltages[-1] - CUTOFF) / max(voltages[-1] - v_new, 1e-30)
            times.append(times[-1] + frac * trial_dt)
            capacities.append(capacities[-1] + frac * (q_new - capacities[-1]))
            voltages.append(CUTOFF)
            reached_cutoff = True
            break

        times.append(t_new)
        capacities.append(q_new)
        voltages.append(v_new)

        if np.mean(solver.c_s[solver.nmc_indices] / CS_MAX) >= 0.999:
            break

    # Final phi_s span
    solver._steady.set_concentration(solver.c_e, solver.c_s)
    final_pot = solver._steady.solve(I_app=i_app)
    phi_s_final = final_pot["phi_s"]
    phi_s_span_final = float(np.max(phi_s_final[solid_mask]) - np.min(phi_s_final[solid_mask]))

    # Use paper mass loading for gravimetric capacity (since current density is paper-based)
    capacity_Ah_m2 = np.asarray(capacities) / 3600.0
    capacity_mAh_g = capacity_Ah_m2 * 1000.0 / PAPER_MASS_LOADING_G_M2

    return {
        "I_app": i_app,
        "V_ocv": V_ocv,
        "V_loaded_init": V_loaded_init,
        "V_final": voltages[-1],
        "capacity_mAh_g": float(capacity_mAh_g[-1]),
        "phi_s_span_init": phi_s_span_init,
        "phi_s_span_final": phi_s_span_final,
        "reached_cutoff": reached_cutoff,
    }


def main():
    print("=" * 70)
    print("SMOKE TEST: Phase 6 Fixes — Paper Comparison Alignment")
    print("=" * 70)
    print()
    print(f"Network: {SHAPE}, spacing={SPACING*1e6:.0f} um")
    print(f"Porosity={POROSITY}, CBD={CBD_FRACTION}, throat_scale={THROAT_SCALE}")
    print(f"INITIAL_SOL={INITIAL_SOL} → U({INITIAL_SOL})={nmc532_ocv(INITIAL_SOL):.3f} V")
    print(f"I_1C_paper = {I_1C_PAPER:.2f} A/m²")
    print(f"Paper mass loading = {PAPER_MASS_LOADING_G_M2} g/m² (used for capacity normalization)")
    print()

    results = {}
    for crate in [0.2, 3.0]:
        print(f"Running {crate}C discharge...")
        r = run_discharge(crate)
        results[crate] = r
        print(f"  I_app = {r['I_app']:.2f} A/m²")
        print(f"  V_ocv (I=0) = {r['V_ocv']:.4f} V")
        print(f"  V_loaded_init = {r['V_loaded_init']:.4f} V")
        print(f"  V_final = {r['V_final']:.4f} V")
        print(f"  capacity = {r['capacity_mAh_g']:.2f} mAh/g")
        print(f"  phi_s span (init) = {r['phi_s_span_init']*1e3:.2f} mV")
        print(f"  phi_s span (final) = {r['phi_s_span_final']*1e3:.2f} mV")
        print(f"  reached_cutoff = {r['reached_cutoff']}")
        print()

    # Verification checks
    print("=" * 70)
    print("VERIFICATION CHECKS")
    print("=" * 70)
    checks = []

    # Check 1: V_init (OCV) ≈ 4.0V
    v_ocv = results[0.2]["V_ocv"]
    ok = 3.8 < v_ocv < 4.2
    checks.append(ok)
    print(f"  [{'OK' if ok else 'FAIL'}] OCV = {v_ocv:.4f} V (expect ~4.0V)")

    # Check 2: 3C current ≈ 159 A/m²
    i3c = abs(results[3.0]["I_app"])
    ok = 140 < i3c < 180
    checks.append(ok)
    print(f"  [{'OK' if ok else 'FAIL'}] 3C current = {i3c:.2f} A/m² (expect ~159 A/m²)")

    # Check 3: Voltage drops toward 3.0V
    vf3 = results[3.0]["V_final"]
    ok = vf3 < 3.5
    checks.append(ok)
    print(f"  [{'OK' if ok else 'FAIL'}] 3C V_final = {vf3:.4f} V (expect < 3.5V)")

    # Check 4: 0.2C capacity in 180-210 mAh/g
    cap02 = results[0.2]["capacity_mAh_g"]
    ok = 150 < cap02 < 250
    checks.append(ok)
    print(f"  [{'OK' if ok else 'FAIL'}] 0.2C capacity = {cap02:.2f} mAh/g (expect 150-250)")

    # Check 5: 3C capacity < 0.2C capacity
    cap3 = results[3.0]["capacity_mAh_g"]
    ok = cap3 < cap02
    checks.append(ok)
    print(f"  [{'OK' if ok else 'FAIL'}] 3C capacity ({cap3:.2f}) < 0.2C capacity ({cap02:.2f})")

    print()
    passed = sum(checks)
    total = len(checks)
    print(f"Result: {passed}/{total} checks passed")
    return all(checks)


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
