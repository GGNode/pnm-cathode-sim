"""
阴极放电仿真脚本 (Cathode Discharge Simulation Runner)
======================================================

使用 PNM 瞬态求解器运行自适应恒流放电仿真。

用法示例:
    python scripts/run_discharge.py --crate 0.2 --shape 5x5x5
    python scripts/run_discharge.py --crate 1.0 --cutoff 3.0 --output data/1c.npz
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

# matplotlib 缓存目录
os.environ.setdefault("MPLCONFIGDIR", str((ROOT / "data/.matplotlib").resolve()))

import numpy as np

from src.network.generator import create_cathode_network
from src.solver.transient import TransientSolver


def parse_shape(value: str) -> list[int]:
    """解析网络形状参数, 如 '5x5x5' 或 '5,5,5'。"""
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
    """
    运行放电仿真。

    流程:
    1. 创建孔网络 (create_cathode_network)
    2. 创建瞬态求解器 (TransientSolver)
    3. 设置初始浓度
    4. 运行放电 (run_discharge)
    """
    # 创建三相孔网络
    net = create_cathode_network(
        shape=args.shape,           # 网络维度 [nx, ny, nz]
        spacing=args.spacing,       # 节点间距 [m]
        porosity=args.porosity,     # 电解质体积分数
        cbd_fraction=args.cbd_fraction,  # CBD 体积分数
        seed=args.seed,             # 随机种子
    )
    # 创建瞬态求解器
    solver = TransientSolver(net, T=args.temperature, k0=args.k0)
    # 设置初始浓度
    solver.set_concentration(c_e=args.ce, c_s=args.cs)
    # 运行放电
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
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crate", type=float, default=0.2, help="放电 C-rate (正值)")
    parser.add_argument("--cutoff", type=float, default=2.5, help="截止电压 [V]")
    parser.add_argument("--shape", type=parse_shape, default=[5, 5, 5], help="网络形状, 如 5x5x5")
    parser.add_argument("--spacing", type=float, default=1e-5, help="节点间距 [m]")
    parser.add_argument("--porosity", type=float, default=0.5, help="电解质孔隙率")
    parser.add_argument("--cbd-fraction", type=float, default=0.10, help="CBD 体积分数")
    parser.add_argument("--seed", type=int, default=42, help="随机种子")
    parser.add_argument("--temperature", type=float, default=298.15, help="温度 [K]")
    parser.add_argument("--k0", type=float, default=5e-10, help="BV 反应速率常数")
    parser.add_argument("--ce", type=float, default=1200.0, help="初始电解质浓度 [mol/m³]")
    parser.add_argument("--cs", type=float, default=24450.0, help="初始固相浓度 [mol/m³]")
    parser.add_argument("--dt", type=float, default=None, help="初始时间步 [s]")
    parser.add_argument("--dt-min", type=float, default=1e-6, help="最小自适应时间步 [s]")
    parser.add_argument("--dt-max", type=float, default=300.0, help="最大自适应时间步 [s]")
    parser.add_argument("--max-steps", type=int, default=80, help="最大时间步数")
    parser.add_argument("--voltage-jump-limit", type=float, default=0.05, help="最大允许电压跳变 [V]")
    parser.add_argument("--output", type=Path, default=Path("data/discharge.npz"), help="NPZ 输出路径")
    return parser


def main() -> None:
    """主入口: 解析参数、运行仿真、保存结果。"""
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
