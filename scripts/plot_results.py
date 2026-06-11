"""
放电曲线绘图脚本 (Discharge Curve Plotter)
==========================================

运行 0.2C 放电仿真并生成 V-Q 曲线图。
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# 将项目根目录加入 Python 路径
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# matplotlib 配置: 无头模式 + 缓存目录
os.environ.setdefault("MPLCONFIGDIR", str((ROOT / "data/.matplotlib").resolve()))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pnmcathode import Cathode, CathodeGeometry, DischargeProtocol, Simulation, SolverSettings


def build_plot_parser() -> argparse.ArgumentParser:
    """构建绘图脚本参数解析器。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.description = __doc__
    parser.add_argument("--c-rate", type=float, default=0.2, dest="c_rate")
    parser.add_argument("--max-steps", type=int, default=80)
    parser.add_argument("--cutoff-voltage", type=float, default=2.5)
    parser.add_argument("--output", type=Path, default=ROOT / "data/discharge_0.2c.npz")
    parser.add_argument("--plot-output", type=Path, default=ROOT / "data/v_q_0.2c.png")
    return parser


def run_simulation(args: argparse.Namespace):
    """使用 pnmcathode API 运行放电仿真。"""
    cathode = Cathode.cubic(
        geometry=CathodeGeometry(shape=(10, 10, 10), seed=42),
    )
    protocol = DischargeProtocol(
        c_rate=args.c_rate,
        cutoff_voltage=args.cutoff_voltage,
        max_steps=args.max_steps,
    )
    result = Simulation(cathode, protocol, settings=SolverSettings()).run()
    result.to_npz(args.output)
    return result


def main() -> None:
    """主入口: 运行仿真并绘制 V-Q 曲线。"""
    args = build_plot_parser().parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.plot_output.parent.mkdir(parents=True, exist_ok=True)

    # 运行放电仿真
    result = run_simulation(args)

    # 绘制 V-Q 曲线
    fig, ax = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)
    ax.plot(result.capacity_Ah_m2, result.voltage, color="#1f77b4", linewidth=2.0)
    ax.set_xlabel("Areal capacity [Ah/m²]")  # 面积比容量 [A·h/m²]
    ax.set_ylabel("Cell voltage [V]")          # 端电压 [V]
    ax.set_title("0.2C discharge")
    ax.grid(True, alpha=0.3)
    fig.savefig(args.plot_output, dpi=200)
    plt.close(fig)

    print(f"saved plot: {args.plot_output}")


if __name__ == "__main__":
    main()
