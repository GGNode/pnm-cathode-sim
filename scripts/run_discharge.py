"""Run an adaptive transient cathode discharge simulation."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

os.environ.setdefault("MPLCONFIGDIR", str((ROOT / "data/.matplotlib").resolve()))

import numpy as np

from src.network.generator import create_cathode_network
from src.solver.transient import TransientSolver


def parse_shape(value: str) -> list[int]:
    parts = value.lower().replace("x", ",").split(",")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("shape must have three dimensions, e.g. 5x5x5")
    try:
        shape = [int(part) for part in parts]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("shape dimensions must be integers") from exc
    if any(dim <= 1 for dim in shape):
        raise argparse.ArgumentTypeError("shape dimensions must be greater than 1")
    return shape


def run_simulation(args: argparse.Namespace) -> dict:
    net = create_cathode_network(
        shape=args.shape,
        spacing=args.spacing,
        porosity=args.porosity,
        cbd_fraction=args.cbd_fraction,
        seed=args.seed,
    )
    solver = TransientSolver(net, T=args.temperature, k0=args.k0)
    solver.set_concentration(c_e=args.ce, c_s=args.cs)
    return solver.run_discharge(
        C_rate=args.crate,
        cutoff_voltage=args.cutoff,
        dt=args.dt,
        dt_min=args.dt_min,
        dt_max=args.dt_max,
        max_steps=args.max_steps,
        voltage_jump_limit=args.voltage_jump_limit,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crate", type=float, default=0.2, help="positive discharge C-rate")
    parser.add_argument("--cutoff", type=float, default=2.5, help="cutoff voltage [V]")
    parser.add_argument("--shape", type=parse_shape, default=[5, 5, 5], help="network shape, e.g. 5x5x5")
    parser.add_argument("--spacing", type=float, default=1e-5, help="network spacing [m]")
    parser.add_argument("--porosity", type=float, default=0.5, help="electrolyte pore fraction")
    parser.add_argument("--cbd-fraction", type=float, default=0.10, help="CBD pore fraction")
    parser.add_argument("--seed", type=int, default=42, help="random seed")
    parser.add_argument("--temperature", type=float, default=298.15, help="temperature [K]")
    parser.add_argument("--k0", type=float, default=5e-10, help="reaction rate constant")
    parser.add_argument("--ce", type=float, default=1200.0, help="initial electrolyte concentration [mol/m3]")
    parser.add_argument("--cs", type=float, default=24450.0, help="initial solid concentration [mol/m3]")
    parser.add_argument("--dt", type=float, default=None, help="initial time step [s]")
    parser.add_argument("--dt-min", type=float, default=1e-6, help="minimum adaptive time step [s]")
    parser.add_argument("--dt-max", type=float, default=300.0, help="maximum adaptive time step [s]")
    parser.add_argument("--max-steps", type=int, default=80, help="maximum accepted time steps")
    parser.add_argument("--voltage-jump-limit", type=float, default=0.05, help="maximum accepted voltage jump [V]")
    parser.add_argument("--output", type=Path, default=Path("data/discharge.npz"), help="NPZ output path")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.crate <= 0:
        raise SystemExit("--crate must be positive for discharge")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = run_simulation(args)
    np.savez(args.output, **result)

    print(f"C-rate: {result['C_rate']:.3g}C")
    print(f"I_app: {result['I_app']:.6g} A/m2")
    print(f"steps: {len(result['time']) - 1}")
    print(f"final time: {result['time'][-1]:.6g} s")
    print(f"final capacity: {result['capacity_Ah_m2'][-1]:.6g} Ah/m2")
    print(f"final voltage: {result['voltage'][-1]:.6g} V")
    print(f"reached cutoff: {bool(result['reached_cutoff'])}")
    print(f"saved: {args.output}")


if __name__ == "__main__":
    main()
