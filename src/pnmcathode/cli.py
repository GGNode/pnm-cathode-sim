"""pnmcathode 命令行入口。"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from pnmcathode import Cathode, CathodeGeometry, DischargeProtocol, Simulation, SolverSettings


def build_parser() -> argparse.ArgumentParser:
    """构建放电仿真命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="pnmcathode-discharge",
        description="运行一个 pnmcathode 恒流放电仿真并保存 NPZ 结果。",
    )
    parser.add_argument("--c-rate", type=float, default=1.0, dest="c_rate")
    parser.add_argument("--max-steps", type=int, default=100)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cutoff-voltage", type=float, default=2.5)
    parser.add_argument("--dt", type=float, default=None)
    return parser


def run_from_args(args: argparse.Namespace) -> "DischargeResult":
    """根据命令行参数运行一次放电仿真。"""
    cathode = Cathode.cubic(geometry=CathodeGeometry(seed=42))
    protocol = DischargeProtocol(
        c_rate=args.c_rate,
        cutoff_voltage=args.cutoff_voltage,
        max_steps=args.max_steps,
    )
    settings = SolverSettings(dt=args.dt)
    result = Simulation(cathode, protocol, settings=settings).run()
    result.to_npz(args.output)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """命令行主函数。"""
    args = build_parser().parse_args(argv)
    result = run_from_args(args)
    print(
        "saved result: "
        f"{args.output} "
        f"(steps={len(result.time) - 1}, final_voltage={result.voltage[-1]:.6g} V)"
    )
    return 0


__all__ = ["build_parser", "main", "run_from_args"]
