from __future__ import annotations

import numpy as np

import pnmcathode
import pnmcathode.analysis as analysis
import pnmcathode.plotting as plotting
from pnmcathode import (
    ActiveMaterial,
    Cathode,
    CathodeGeometry,
    ConductiveAdditive,
    DischargeProtocol,
    DischargeResult,
    Electrolyte,
    Kinetics,
    Separator,
    Simulation,
    SolverSettings,
)
from pnmcathode.materials import ActiveMaterial as MaterialActiveMaterial
from pnmcathode.materials.presets import (
    cbd_khan2021,
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


def constant_transport(_concentration: float, _temperature: float) -> float:
    return 1.0


def constant_ocv(_soc: float) -> float:
    return 4.0


def test_config_dataclasses_and_public_exports():
    active = ActiveMaterial(
        name="test active",
        cs_max=100.0,
        sigma=2.0,
        diffusivity=constant_transport,
        ocv=constant_ocv,
    )
    electrolyte = Electrolyte(
        name="test electrolyte",
        diffusivity=constant_transport,
        conductivity=constant_transport,
    )
    additive = ConductiveAdditive(name="test additive", sigma=5.0)
    kinetics = Kinetics(k0=1e-9)
    separator = Separator(enabled=True)

    assert pnmcathode.__version__ == "0.1.0"
    assert MaterialActiveMaterial is ActiveMaterial
    assert active.ocv(0.5) == 4.0
    assert electrolyte.c_init == 1200.0
    assert additive.sigma == 5.0
    assert kinetics.alpha_a == 0.5
    assert separator.enabled is True
    assert hasattr(analysis, "pore_size_distribution")
    assert hasattr(plotting, "plot_discharge_curve")


def test_cathode_cubic_creation():
    geometry = CathodeGeometry(shape=(3, 3, 3), porosity=0.45, seed=7)

    cathode = Cathode.cubic(geometry=geometry)

    assert cathode.network is not None
    assert cathode.network.Np == 27
    assert cathode.active_material.name == "NMC532 (Khan 2021)"
    assert cathode.electrolyte.name == "LiPF6 carbonate electrolyte (Khan 2021)"
    assert np.any(cathode.network["pore.electrolyte"])
    assert np.any(cathode.network["pore.nmc"])


def test_simulation_run_returns_discharge_result(monkeypatch):
    calls = {}

    class DummyTransientSolver:
        def __init__(self, net, T, k0, geometric_area=None, separator=None):
            calls["net"] = net
            calls["T"] = T
            calls["k0"] = k0
            calls["geometric_area"] = geometric_area
            calls["separator"] = separator

        def set_concentration(self, c_e, c_s):
            calls["c_e"] = c_e
            calls["c_s"] = c_s

        def run_discharge(
            self,
            C_rate,
            cutoff_voltage,
            dt,
            dt_min,
            dt_max,
            max_steps,
            voltage_jump_limit,
        ):
            calls["run"] = {
                "C_rate": C_rate,
                "cutoff_voltage": cutoff_voltage,
                "dt": dt,
                "dt_min": dt_min,
                "dt_max": dt_max,
                "max_steps": max_steps,
                "voltage_jump_limit": voltage_jump_limit,
            }
            return {
                "time": np.array([0.0, 1.0]),
                "capacity_Ah_m2": np.array([0.0, 0.01]),
                "voltage": np.array([4.0, 3.9]),
                "I_app": -36.0,
                "C_rate": C_rate,
                "cutoff_voltage": cutoff_voltage,
                "reached_cutoff": False,
            }

    monkeypatch.setattr("pnmcathode.simulation.TransientSolver", DummyTransientSolver)

    cathode = Cathode.cubic(
        geometry=CathodeGeometry(shape=(3, 3, 3), seed=1, geometric_area=1.2e-6)
    )
    protocol = DischargeProtocol(c_rate=1.0, cutoff_voltage=2.5, max_steps=3)
    settings = SolverSettings(temperature=303.0, dt=0.1)

    result = Simulation(cathode, protocol, settings=settings).run()

    assert isinstance(result, DischargeResult)
    assert result.final_capacity == 0.01
    assert result.current_density == -36.0
    assert result.c_rate == 1.0
    assert calls["T"] == 303.0
    assert calls["k0"] == 5e-10
    assert calls["geometric_area"] == 1.2e-6
    assert calls["c_e"] == 1200.0
    assert calls["c_s"] == 24450.0
    assert calls["run"]["max_steps"] == 3


def test_presets_load_correctly():
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    additive = cbd_khan2021()
    separator = separator_khan2021()

    assert active.cs_max == 48900.0
    assert active.sigma == 0.01
    assert callable(active.diffusivity)
    assert callable(active.ocv)
    assert electrolyte.c_init == 1200.0
    assert callable(electrolyte.diffusivity)
    assert callable(electrolyte.conductivity)
    assert additive.sigma == 760.0
    assert separator.enabled is True
    assert separator.thickness == 25e-6


def test_discharge_result_serialization_round_trip(tmp_path):
    path = tmp_path / "result.npz"
    result = DischargeResult(
        time=np.array([0.0, 2.0]),
        voltage=np.array([4.1, 3.8]),
        capacity_Ah_m2=np.array([0.0, 0.02]),
        current_density=-12.5,
        c_rate=0.5,
        reached_cutoff=True,
        spatial={"c_e": np.array([[1200.0, 1190.0]])},
        metadata={"case": "round_trip"},
    )

    result.to_npz(path)
    loaded = DischargeResult.from_npz(path)

    np.testing.assert_allclose(loaded.time, result.time)
    np.testing.assert_allclose(loaded.voltage, result.voltage)
    np.testing.assert_allclose(loaded.capacity_Ah_m2, result.capacity_Ah_m2)
    assert loaded.current_density == result.current_density
    assert loaded.c_rate == result.c_rate
    assert loaded.reached_cutoff is True
    assert loaded.metadata == {"case": "round_trip"}
    assert loaded.spatial is not None
    np.testing.assert_allclose(loaded.spatial["c_e"], result.spatial["c_e"])
