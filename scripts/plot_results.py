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

from scripts.run_discharge import build_parser, run_simulation


def build_plot_parser() -> argparse.ArgumentParser:
    """构建绘图专用参数解析器 (扩展 run_discharge 的参数)。"""
    parser = build_parser()
    parser.description = __doc__
    parser.set_defaults(crate=0.2, output=ROOT / "data/discharge_0.2c.npz", max_steps=80)
    parser.add_argument("--plot-output", type=Path, default=ROOT / "data/v_q_0.2c.png")
    return parser


def main() -> None:
    """主入口: 运行仿真并绘制 V-Q 曲线。"""
    args = build_plot_parser().parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.plot_output.parent.mkdir(parents=True, exist_ok=True)

    # 运行放电仿真
    result = run_simulation(args)

    # 绘制 V-Q 曲线
    fig, ax = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)
    ax.plot(result["capacity_Ah_m2"], result["voltage"], color="#1f77b4", linewidth=2.0)
    ax.set_xlabel("Areal capacity [Ah/m²]")  # 面积比容量 [A·h/m²]
    ax.set_ylabel("Cell voltage [V]")          # 端电压 [V]
    ax.set_title("0.2C discharge")
    ax.grid(True, alpha=0.3)
    fig.savefig(args.plot_output, dpi=200)
    plt.close(fig)

    print(f"saved plot: {args.plot_output}")


if __name__ == "__main__":
    main()
