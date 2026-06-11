"""打包与 CLI smoke 测试。"""

from __future__ import annotations

import numpy as np


def test_import_pnmcathode_and_subpackages():
    """主包和所有子包应可导入。"""
    import pnmcathode
    import pnmcathode.materials
    import pnmcathode.network
    import pnmcathode.physics
    import pnmcathode.post
    import pnmcathode.solver

    assert pnmcathode.__version__ == "0.1.0"


def test_top_level_api_available():
    """顶层 API 应暴露核心对象。"""
    from pnmcathode import Cathode, DischargeProtocol, Simulation

    assert Cathode.__name__ == "Cathode"
    assert DischargeProtocol.__name__ == "DischargeProtocol"
    assert Simulation.__name__ == "Simulation"


def test_short_simulation_runs():
    """极短仿真应能端到端跑通。"""
    from pnmcathode import Cathode, CathodeGeometry, DischargeProtocol, Simulation, SolverSettings

    cathode = Cathode.cubic(geometry=CathodeGeometry(shape=(3, 3, 3), seed=3))
    result = Simulation(
        cathode,
        DischargeProtocol(c_rate=0.2, max_steps=2),
        settings=SolverSettings(dt=0.01, dt_max=0.01),
    ).run()

    assert len(result.time) == 3
    assert np.all(np.isfinite(result.voltage))
    assert result.final_capacity > 0.0


def test_cli_writes_npz(tmp_path):
    """CLI 入口应能运行并保存 NPZ 结果。"""
    from pnmcathode.cli import main
    from pnmcathode.results import DischargeResult

    output = tmp_path / "result.npz"
    exit_code = main(
        [
            "--c-rate",
            "0.2",
            "--max-steps",
            "1",
            "--dt",
            "0.01",
            "--output",
            str(output),
        ]
    )

    assert exit_code == 0
    assert output.exists()
    loaded = DischargeResult.from_npz(output)
    assert len(loaded.time) == 2
    assert loaded.c_rate == 0.2
