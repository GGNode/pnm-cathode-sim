"""
Khan 2021 1CAL/3CAL P2D 基准算例运行器。

用法:
    python examples/p2d/run_khan_p2d.py

输出:
    data/validation/p2d_khan_1cal_1c.npz
    data/validation/p2d_khan_3cal_1c.npz

对应 TASK_PHASE1.md §8
"""

from __future__ import annotations

from pathlib import Path

from pnmcathode.p2d.benchmark import (
    P2DBenchmarkCase,
    make_khan_solver,
    make_parameter_snapshot,
    run_p2d_benchmark,
    save_benchmark_npz,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "validation"


def main() -> None:
    """运行 Khan 2021 1CAL 和 3CAL P2D benchmark。"""
    cases = [
        P2DBenchmarkCase(
            name="khan_1cal_1c",
            parameter_source="khan2021",
            electrode_kind="1CAL",
            c_rate=1.0,
            initial_soc=0.35,
            temperature=303.0,
            cutoff_voltage=3.0,
            n_positive_cells=30,
            n_separator_cells=8,
            n_particle_shells=15,
            dt_initial=5.0,
            t_final=5000.0,
        ),
        P2DBenchmarkCase(
            name="khan_3cal_1c",
            parameter_source="khan2021",
            electrode_kind="3CAL",
            c_rate=1.0,
            initial_soc=0.35,
            temperature=303.0,
            cutoff_voltage=3.0,
            n_positive_cells=30,
            n_separator_cells=8,
            n_particle_shells=15,
            dt_initial=5.0,
            t_final=5000.0,
        ),
    ]

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    for case in cases:
        print(f"\n{'='*60}")
        print(f"运行: {case.name}")
        print(f"  电极: {case.electrode_kind}")
        print(f"  C-rate: {case.c_rate}")
        print(f"  网格: {case.n_positive_cells} positive cells, "
              f"{case.n_separator_cells} separator cells, "
              f"{case.n_particle_shells} shells")

        solver = make_khan_solver(case)
        snap = make_parameter_snapshot(solver, source="khan2021")
        print(f"  参数源: {snap.source}")
        print(f"  温度: {snap.temperature} K")
        print(f"  c_s_max: {snap.c_s_max} mol/m³")

        result = run_p2d_benchmark(case)
        V0 = result.voltage[0]
        Vf = result.voltage[-1]
        print(f"  初始电压: {V0:.4f} V")
        print(f"  最终电压: {Vf:.4f} V")
        print(f"  步数: {result.metadata.get('steps', 'N/A')}")

        out_path = DATA_DIR / f"p2d_{case.name}.npz"
        save_benchmark_npz(out_path, result, solver)
        print(f"  保存: {out_path}")


if __name__ == "__main__":
    main()
