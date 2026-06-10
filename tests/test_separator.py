"""Separator boundary model tests.

Verifies:
  1. Zero-current identity: no drops, c_cathode = c_ref
  2. Ohmic linearity: dphi_ohm(2I) / dphi_ohm(I) ~ 2
  3. Monotone concentration: c_cathode decreases with I, bounded by c_floor
  4. Voltage drop: separator-on < separator-off at nonzero current
"""

import numpy as np
import pytest

from src.network.generator import create_cathode_network
from src.physics.separator import SeparatorParams, separator_boundary
from src.solver.steady import SteadyStateSolver


class TestSeparatorBoundary:
    """Unit tests for separator_boundary() function."""

    def test_zero_current_identity(self):
        """At I_app=0, phi_e_cathode=0, c_cathode=c_ref, all drops zero."""
        params = SeparatorParams(enabled=True)
        state = separator_boundary(0.0, 303.0, params)
        assert state.phi_e_cathode == pytest.approx(0.0, abs=1e-15)
        assert state.c_cathode == pytest.approx(params.c_ref)
        assert state.dphi_ohm == pytest.approx(0.0, abs=1e-15)
        assert state.dphi_conc == pytest.approx(0.0, abs=1e-15)
        assert state.eta_li == pytest.approx(0.0, abs=1e-15)

    def test_ohmic_linearity(self):
        """Ohmic drop should be linear in current for small currents."""
        params = SeparatorParams(enabled=True, include_concentration_overpotential=False)
        I1, I2 = -0.1, -0.2
        s1 = separator_boundary(I1, 303.0, params)
        s2 = separator_boundary(I2, 303.0, params)
        assert s1.dphi_ohm > 0.0
        np.testing.assert_allclose(s2.dphi_ohm / s1.dphi_ohm, 2.0, rtol=0.01)

    def test_monotone_concentration(self):
        """Cathode-side concentration should decrease monotonically with |I|."""
        params = SeparatorParams(enabled=True)
        currents = [-0.1, -1.0, -10.0, -100.0]
        concs = [separator_boundary(I, 303.0, params).c_cathode for I in currents]
        for a, b in zip(concs, concs[1:]):
            assert b <= a, f"Concentration not decreasing: {a} -> {b}"
        # Must stay above floor
        for c in concs:
            assert c >= params.c_floor

    def test_disabled_returns_identity(self):
        """When disabled, separator_boundary should return identity values."""
        params = SeparatorParams(enabled=False)
        state = separator_boundary(-100.0, 303.0, params)
        assert state.phi_e_cathode == 0.0
        assert state.c_cathode == 1200.0
        assert state.dphi_ohm == 0.0
        assert state.dphi_conc == 0.0
        assert state.eta_li == 0.0


class TestSeparatorSteadySolver:
    """Integration tests: separator with SteadyStateSolver."""

    def _make_solver(self, separator=None):
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        solver = SteadyStateSolver(net, T=298.15, separator=separator)
        solver.set_concentration(c_e=1200.0, c_s=24450.0)
        return solver

    def test_separator_on_voltage_below_separator_off(self):
        """With separator enabled and nonzero discharge, V_cell should be lower."""
        sep_params = SeparatorParams(enabled=True)
        solver_off = self._make_solver(separator=None)
        solver_on = self._make_solver(separator=sep_params)

        V_off = solver_off.solve(I_app=-1.0)["voltage"]
        V_on = solver_on.solve(I_app=-1.0)["voltage"]
        assert V_on < V_off, (
            f"Separator-on voltage ({V_on:.4f}V) should be lower than "
            f"separator-off ({V_off:.4f}V)"
        )

    def test_separator_zero_current_matches(self):
        """At zero current, separator-on and separator-off should give same V."""
        sep_params = SeparatorParams(enabled=True)
        solver_off = self._make_solver(separator=None)
        solver_on = self._make_solver(separator=sep_params)

        V_off = solver_off.solve(I_app=0.0)["voltage"]
        V_on = solver_on.solve(I_app=0.0)["voltage"]
        assert abs(V_on - V_off) < 1e-10, (
            f"Zero-current mismatch: on={V_on:.6f}V vs off={V_off:.6f}V"
        )

    def test_backward_compatibility_default(self):
        """Default separator (disabled) must be identical to no separator."""
        solver_default = self._make_solver(separator=None)
        solver_explicit = self._make_solver(
            separator=SeparatorParams(enabled=False),
        )

        for I_app in [0.0, -0.001, -0.1]:
            r1 = solver_default.solve(I_app=I_app)
            r2 = solver_explicit.solve(I_app=I_app)
            assert abs(r1["voltage"] - r2["voltage"]) < 1e-12, (
                f"Voltage mismatch at I_app={I_app}: "
                f"default={r1['voltage']:.6f}, explicit={r2['voltage']:.6f}"
            )
