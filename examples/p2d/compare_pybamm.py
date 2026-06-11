"""
PyBaMM DFN 参考数据生成器 (可选)。

此脚本需要在安装了 PyBaMM 的环境中运行:
    pip install pybamm
    python examples/p2d/compare_pybamm.py

生成:
    data/validation/p2d_pybamm_*.npz

若 PyBaMM 未安装，脚本会打印提示并跳过。

对应 TASK_PHASE1.md §2.2, §8
"""

from __future__ import annotations

import sys
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "validation"


def main() -> None:
    """生成 PyBaMM DFN 参考数据。"""
    try:
        import pybamm
    except ImportError:
        print("PyBaMM 未安装。跳过参考数据生成。")
        print("安装方法: pip install pybamm")
        sys.exit(0)

    print(f"PyBaMM 版本: {pybamm.__version__}")

    # 使用 Chen 2020 参数集作为参考
    parameter_set = "Chen2020"
    model = pybamm.lithium_ion.DFN()
    param = pybamm.ParameterValues(parameter_set)

    # 记录关键参数
    print(f"参数集: {parameter_set}")
    print("注意: 这是 full-cell DFN，与 half-cell P2D 不完全等价。")
    print("比较时只使用 positive-side 趋势，不作为硬性误差标准。")

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    for c_rate in [0.5, 1.0, 3.0]:
        print(f"\n运行 PyBaMM DFN {c_rate}C ...")

        sim = pybamm.Simulation(
            model,
            parameter_values=param,
            solver=pybamm.CasadiSolver(mode="fast"),
        )
        sim.solve([0, 3600 / c_rate * 1.5])  # 足够长以捕获完整放电

        # 提取 positive electrode 数据
        V = sim.solution["Terminal voltage [V]"].entries
        t = sim.solution["Time [s]"].entries
        I = sim.solution["Current [A]"].entries

        out_path = DATA_DIR / f"p2d_pybamm_{c_rate}c.npz"
        import numpy as np

        np.savez(
            out_path,
            time=t,
            voltage=V,
            current=I,
            parameter_set=parameter_set,
            pybamm_version=pybamm.__version__,
            model_kind="DFN",
            c_rate=c_rate,
        )
        print(f"  保存: {out_path}")


if __name__ == "__main__":
    main()
