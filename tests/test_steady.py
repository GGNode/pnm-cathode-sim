"""Phase 3: Steady-state solver tests.

Verifies:
  1. Uniform potential in a single-phase network with no current
  2. Potential drop across network under applied current
  3. Butler-Volmer coupling at electrolyte/NMC interface
  4. Low C-rate V-Q curve near equilibrium
"""

import numpy as np
import pytest

from src.network.generator import create_cathode_network
from src.solver.steady import SteadyStateSolver


class TestSteadyStateSolver:
    """Tests for the steady-state potential solver."""

    def test_uniform_potential_no_current(self):
        """With no applied current, potential should be uniform."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = SteadyStateSolver(net, T=298.15)
        # Set uniform concentration
        solver.set_concentration(c_e=1200.0, c_s=24450.0)
        # Solve with zero applied current
        result = solver.solve(I_app=0.0)
        phi_e = result["phi_e"]
        phi_s = result["phi_s"]
        # All electrolyte pores should have same phi_e
        e_mask = net["pore.electrolyte"]
        if np.any(e_mask):
            phi_e_range = phi_e[e_mask].max() - phi_e[e_mask].min()
            assert phi_e_range < 1e-10, (
                f"phi_e range = {phi_e_range:.2e}, expected ~0"
            )
        # All NMC pores should have same phi_s
        nmc_mask = net["pore.nmc"]
        if np.any(nmc_mask):
            phi_s_range = phi_s[nmc_mask].max() - phi_s[nmc_mask].min()
            assert phi_s_range < 1e-10, (
                f"phi_s range = {phi_s_range:.2e}, expected ~0"
            )

    def test_potential_drop_with_current(self):
        """Applied current should create a potential gradient."""
        net = create_cathode_network(
            shape=[10, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = SteadyStateSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)
        result = solver.solve(I_app=-0.1)  # small cathodic current
        phi_e = result["phi_e"]
        # There should be a potential gradient in x-direction
        e_mask = net["pore.electrolyte"]
        if np.sum(e_mask) > 1:
            coords = net["pore.coords"][e_mask]
            x_coords = coords[:, 0]
            # Sort by x and check for gradient
            sorted_idx = np.argsort(x_coords)
            phi_sorted = phi_e[e_mask][sorted_idx]
            # At least some variation expected
            assert phi_sorted.max() - phi_sorted.min() > 1e-12

    def test_voltage_near_ocv_at_low_current(self):
        """At very low current, cell voltage should be near OCV."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = SteadyStateSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)
        result = solver.solve(I_app=-0.001)  # very small current
        V_cell = result["voltage"]
        # Should be near OCV at x=0.5 (~3.73V)
        from src.physics.ocv import nmc532_ocv
        U_eq = nmc532_ocv(0.5)
        assert abs(V_cell - U_eq) < 0.1, (
            f"V_cell={V_cell:.4f}V, expected near U_eq={U_eq:.4f}V"
        )

    def test_small_current_polarization_scales_linearly(self):
        """Low-current voltage polarization should be linear in current magnitude."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = SteadyStateSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)

        ocv = solver.solve(I_app=0.0)["voltage"]
        pol_low = ocv - solver.solve(I_app=-1e-4)["voltage"]
        pol_high = ocv - solver.solve(I_app=-2e-4)["voltage"]

        assert pol_low > 0.0
        np.testing.assert_allclose(pol_high / pol_low, 2.0, rtol=0.1)

    def test_result_keys(self):
        """Result should contain expected keys."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = SteadyStateSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)
        result = solver.solve(I_app=0.0)
        for key in ["phi_e", "phi_s", "voltage", "I_rxn"]:
            assert key in result, f"Missing key: {key}"
