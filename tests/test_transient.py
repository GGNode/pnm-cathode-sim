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
from src.physics.ocv import nmc532_ocv
from src.physics.reaction import F
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

    def test_discharge_lithiates_solid(self):
        """Cathodic discharge should increase average solid Li concentration."""
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
        assert cs_final > cs_init, "c_s did not increase during cathodic discharge"

    def test_reaction_lithium_conservation_at_half_c(self):
        """Li changes should match Faradaic source plus separator reservoir flux."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = TransientSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)

        e_mask = net["pore.electrolyte"]
        nmc_mask = net["pore.nmc"]
        vol_e = net["pore.volume"][e_mask]
        vol_s = net["pore.volume"][nmc_mask]
        I_app = solver.current_density_for_c_rate(0.5)
        dt = min(solver.characteristic_dt(I_app=I_app, C_rate=0.5), 2.0)

        for _ in range(5):
            c_e_old = solver.c_e[e_mask].copy()
            c_old = solver.c_s[nmc_mask].copy()
            result = solver.step(dt=dt, I_app=I_app)
            c_e_new = solver.c_e[e_mask].copy()
            c_new = solver.c_s[nmc_mask].copy()

            delta_e = np.sum(vol_e * (c_e_new - c_e_old))
            delta_s = np.sum(vol_s * (c_new - c_old))
            reaction_mol_rate = 0.0
            for t, *_ in solver._interfaces:
                reaction_mol_rate += result["I_rxn"][t] * net["throat.area"][t] / F
            rhs_e = dt * reaction_mol_rate
            rhs_s = -dt * reaction_mol_rate

            total_lithium = max(np.sum(vol_e * c_e_old) + np.sum(vol_s * c_old), 1e-30)
            err_e = abs(delta_e - rhs_e)
            err_s = abs(delta_s - rhs_s)
            assert err_e / total_lithium < 0.01, (
                f"electrolyte Li balance error too large: lhs={delta_e:.4e}, rhs={rhs_e:.4e}"
            )
            assert err_s / total_lithium < 0.01, (
                f"solid Li balance error too large: lhs={delta_s:.4e}, rhs={rhs_s:.4e}"
            )
            reservoir_flux = delta_e + delta_s
            assert reservoir_flux >= -1e-18, (
                f"separator reservoir should not remove Li during discharge: {reservoir_flux:.4e}"
            )
            sep_e = solver._steady.sep_e
            np.testing.assert_allclose(solver.c_e[solver.e_indices[sep_e]], 1200.0)

    def test_geometric_area_overrides_c_rate_area(self):
        """C-rate conversion and current BC should use configured geometric area."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        geometric_area = (5e-5) ** 2
        solver = TransientSolver(net, T=298.15, geometric_area=geometric_area)
        expected = -solver.discharge_capacity_coulombs() / 3600.0 / geometric_area

        assert solver.collector_area == geometric_area
        assert solver._steady._cc_area == geometric_area
        assert solver.current_density_for_c_rate(1.0) == pytest.approx(expected)

    def test_zero_current_limit_recovers_ocv(self):
        """At I_app=0, eta=0, I_rxn=0, and V_cell=U(x0)."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = TransientSolver(net, T=298.15)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)

        result = solver.step(dt=1.0, I_app=0.0)
        U_eq = nmc532_ocv(24450.0 / 48900.0)

        assert np.allclose(result["eta"], 0.0, atol=1e-10)
        assert np.allclose(result["I_rxn"], 0.0, atol=1e-12)
        assert abs(result["voltage"] - U_eq) < 1e-10
