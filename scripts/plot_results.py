"""Generate a 0.2C voltage-capacity discharge plot."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

os.environ.setdefault("MPLCONFIGDIR", str((ROOT / "data/.matplotlib").resolve()))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.run_discharge import build_parser, run_simulation


def build_plot_parser() -> argparse.ArgumentParser:
    parser = build_parser()
    parser.description = __doc__
    parser.set_defaults(crate=0.2, output=ROOT / "data/discharge_0.2c.npz", max_steps=80)
    parser.add_argument("--plot-output", type=Path, default=ROOT / "data/v_q_0.2c.png")
    return parser


def main() -> None:
    args = build_plot_parser().parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.plot_output.parent.mkdir(parents=True, exist_ok=True)

    result = run_simulation(args)

    fig, ax = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)
    ax.plot(result["capacity_Ah_m2"], result["voltage"], color="#1f77b4", linewidth=2.0)
    ax.set_xlabel("Areal capacity [Ah/m2]")
    ax.set_ylabel("Cell voltage [V]")
    ax.set_title("0.2C discharge")
    ax.grid(True, alpha=0.3)
    fig.savefig(args.plot_output, dpi=200)
    plt.close(fig)

    print(f"saved plot: {args.plot_output}")


if __name__ == "__main__":
    main()
