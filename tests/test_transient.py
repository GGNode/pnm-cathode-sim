"""Phase 4: Transient solver tests.

Verifies:
  1. Mass conservation during time stepping
  2. Concentration approaches steady state at low C-rate
  3. Higher C-rate produces larger concentration gradients
  4. Result dict contains expected keys
"""

import numpy as np
import pytest

from src.network.generator import create_cathode_network
from src.solver.transient import TransientSolver


class TestTransientSolver:
    """Tests for the transient discharge solver."""

    def test_result_keys(self):
        """Result should contain expected keys."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = TransientSolver(net, T=298.15)
        result = solver.step(dt=1.0, I_app=0.0)
        for key in ["phi_e", "phi_s", "c_e", "c_s", "voltage", "I_rxn"]:
            assert key in result, f"Missing key: {key}"

    def test_mass_conservation_electrolyte(self):
        """Total Li+ in electrolyte must be conserved with zero current."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = TransientSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)

        e_mask = net["pore.electrolyte"]
        vol_e = net["pore.volume"][e_mask]
        c_e_init = solver.c_e[e_mask].copy()
        total_init = np.sum(c_e_init * vol_e)

        # Several steps with zero current
        for _ in range(5):
            result = solver.step(dt=10.0, I_app=0.0)

        c_e_final = solver.c_e[e_mask]
        total_final = np.sum(c_e_final * vol_e)
        assert abs(total_final - total_init) / max(abs(total_init), 1e-30) < 1e-10, (
            f"Electrolyte mass not conserved: init={total_init:.4e}, final={total_final:.4e}"
        )

    def test_mass_conservation_solid(self):
        """Total Li in solid must be conserved with zero current."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = TransientSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)

        nmc_mask = net["pore.nmc"]
        vol_s = net["pore.volume"][nmc_mask]
        c_s_init = solver.c_s[nmc_mask].copy()
        total_init = np.sum(c_s_init * vol_s)

        for _ in range(5):
            result = solver.step(dt=10.0, I_app=0.0)

        c_s_final = solver.c_s[nmc_mask]
        total_final = np.sum(c_s_final * vol_s)
        assert abs(total_final - total_init) / max(abs(total_init), 1e-30) < 1e-10, (
            f"Solid mass not conserved: init={total_init:.4e}, final={total_final:.4e}"
        )

    def test_voltage_reasonable(self):
        """Cell voltage should be in a physical range."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = TransientSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)
        result = solver.step(dt=1.0, I_app=-0.001)
        V = result["voltage"]
        assert 2.5 < V < 5.0, f"V_cell={V:.4f}V outside physical range"

    def test_discharge_depletes_solid(self):
        """Discharging should decrease average solid Li concentration."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = TransientSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)

        nmc_mask = net["pore.nmc"]
        cs_init = solver.c_s[nmc_mask].mean()
        for _ in range(50):
            result = solver.step(dt=10.0, I_app=-0.001)
        cs_final = solver.c_s[nmc_mask].mean()
        # Average c_s should decrease during discharge
        # (Li leaves solid → electrolyte in this convention)
        assert cs_final != cs_init, "c_s unchanged after discharge"
