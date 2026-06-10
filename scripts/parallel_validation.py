"""多 C-rate 并行放电验证。

使用 multiprocessing 并行运行 0.2C、0.5C、1C、3C 放电。
网络: 13x13x13, spacing=10μm, 匹配论文 1CAL 厚度 ~130μm。
"""

from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np


@dataclass
class RunConfig:
    """单次放电运行的配置参数。"""
    C_rate: float
    shape: list[int]
    spacing: float
    porosity: float
    cbd_fraction: float
    seed: int
    throat_scale: float
    T: float
    k0: float
    c_e0: float
    c_s0_sol: float
    cs_max: float
    cutoff_V: float
    current_basis: str  # "network" 或 "paper_areal"
    I_1C_paper: float  # paper_areal 模式使用
    separator_enabled: bool


def run_one(cfg: RunConfig) -> dict:
    """单个 C-rate 放电 — 在子进程中运行。"""
    from src.network.generator import create_cathode_network
    from src.physics.separator import SeparatorParams
    from src.solver.transient import TransientSolver

    net = create_cathode_network(
        shape=cfg.shape, spacing=cfg.spacing, porosity=cfg.porosity,
        cbd_fraction=cfg.cbd_fraction, seed=cfg.seed,
        throat_scale=cfg.throat_scale,
    )
    geo_area = (cfg.shape[0] * cfg.spacing) ** 2
    sep_params = SeparatorParams(enabled=cfg.separator_enabled) if cfg.separator_enabled else None
    solver = TransientSolver(net, T=cfg.T, k0=cfg.k0, geometric_area=geo_area, separator=sep_params)

    solver.set_concentration(c_e=cfg.c_e0, c_s=cfg.c_s0_sol * cfg.cs_max)

    if cfg.current_basis == "network":
        i_app = solver.current_density_for_c_rate(cfg.C_rate, geometric_area=geo_area)
        I_1C_init = abs(solver.current_density_for_c_rate(1.0, geometric_area=geo_area))
    else:
        i_app = -cfg.C_rate * cfg.I_1C_paper
        I_1C_init = cfg.I_1C_paper

    # 运行放电
    dt = min(5.0, 50.0 / cfg.C_rate, solver.characteristic_dt(I_app=i_app, C_rate=cfg.C_rate))
    max_time = 3.0 * 3600.0 / cfg.C_rate  # 3 倍标称放电时间
    max_steps = int(max_time / dt) + 2000

    pot0 = solver._steady.solve(I_app=i_app, tol=1e-8)
    v_init = pot0["voltage"]

    times = [0.0]
    voltages = [v_init]
    t0 = time.time()

    for step_i in range(max_steps):
        try:
            result, dt_used, _ = solver.adaptive_step(dt=dt, I_app=i_app, dt_max=30.0)
        except Exception:
            break

        t_new = times[-1] + dt_used
        v_new = result["voltage"]

        if t_new >= max_time:
            times.append(t_new)
            voltages.append(v_new)
            break

        if voltages[-1] > cfg.cutoff_V >= v_new:
            frac = (voltages[-1] - cfg.cutoff_V) / max(voltages[-1] - v_new, 1e-30)
            times.append(times[-1] + frac * dt_used)
            voltages.append(cfg.cutoff_V)
            break

        times.append(t_new)
        voltages.append(v_new)

        dt = min(dt * 1.2, 30.0)

        # 提前停止: 电压收敛 (最近 200 步无显著变化)
        if step_i > 200 and len(voltages) > 200:
            recent = voltages[-200:]
            if np.max(recent) - np.min(recent) < 0.005:  # <5mV variation
                break

    elapsed = time.time() - t0
    times = np.array(times)
    voltages = np.array(voltages)

    # 质量负载
    geo = solver.net["pore.volume"][solver.nmc_indices]
    nmc_vol = geo.sum()
    mass_loading = nmc_vol * 4.75e6 / geo_area  # g/m²

    # 容量 (mAh/g)
    i_dis = abs(i_app)
    cap_ah_m2 = i_dis * times / 3600.0
    cap_mah_g = cap_ah_m2 / (mass_loading / 1000.0)

    reached_cutoff = voltages[-1] <= cfg.cutoff_V + 0.01

    # 隔膜损失
    sep_ohm = sep_conc = 0.0
    if cfg.separator_enabled:
        from src.physics.separator import separator_boundary, SeparatorParams as SP
        sep = separator_boundary(i_app, cfg.T, SP(enabled=True))
        sep_ohm = sep.dphi_ohm * 1000  # mV
        sep_conc = sep.dphi_conc * 1000

    return {
        "C_rate": cfg.C_rate,
        "current_basis": cfg.current_basis,
        "I_app": float(i_app),
        "I_1C": float(I_1C_init),
        "V_init": float(voltages[0]),
        "V_final": float(voltages[-1]),
        "capacity_mAh_g": float(cap_mah_g[-1]),
        "reached_cutoff": bool(reached_cutoff),
        "steps": len(times),
        "elapsed_s": round(elapsed, 1),
        "mass_loading_g_m2": round(mass_loading, 1),
        "separator_ohm_mV": round(sep_ohm, 2),
        "separator_conc_mV": round(sep_conc, 2),
        "shape": cfg.shape,
        "thickness_um": (cfg.shape[0] - 1) * cfg.spacing * 1e6,
        # V-Q 曲线 (降采样)
        "V_curve": voltages[::max(1, len(voltages)//100)].tolist(),
        "Q_curve": cap_mah_g[::max(1, len(cap_mah_g)//100)].tolist(),
    }


def main():
    base = dict(
        shape=[13, 13, 13],
        spacing=1e-5,
        porosity=0.368,
        cbd_fraction=0.1392,
        seed=42,
        throat_scale=1.0,
        T=303.0,
        k0=1e-10,
        c_e0=1200.0,
        c_s0_sol=0.50,
        cs_max=48900.0,
        cutoff_V=3.0,
    )

    configs = []
    for C in [0.2, 0.5, 1.0, 3.0]:
        for basis in ["network", "paper_areal"]:
            for sep in [False, True]:
                cfg = RunConfig(
                    C_rate=C, current_basis=basis, separator_enabled=sep,
                    I_1C_paper=53.01, **base,
                )
                configs.append(cfg)

    print(f"Running {len(configs)} scenarios on 8 cores...")
    print()

    results = []
    t_total = time.time()

    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(run_one, cfg): cfg for cfg in configs}
        for future in as_completed(futures):
            cfg = futures[future]
            try:
                r = future.result()
                results.append(r)
                label = f"{r['current_basis']}+sep={'ON' if cfg.separator_enabled else 'OFF'}"
                print(f"  {cfg.C_rate}C {label:30s} V0={r['V_init']:.3f} Vf={r['V_final']:.3f} "
                      f"cap={r['capacity_mAh_g']:.1f} cutoff={r['reached_cutoff']} "
                      f"steps={r['steps']} time={r['elapsed_s']}s")
            except Exception as e:
                print(f"  {cfg.C_rate}C FAILED: {e}")

    print(f"\nTotal wall time: {time.time()-t_total:.1f}s")

    # 保存结果
    out = Path(ROOT / "data" / "validation")
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "parallel_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {out / 'parallel_results.json'}")


if __name__ == "__main__":
    main()
