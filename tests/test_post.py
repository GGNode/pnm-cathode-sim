"""Phase 5: Post-processing tests.

Verifies:
  1. Structural analysis (pore size distribution, coordination number)
  2. Discharge curve computation
  3. State-of-lithiation calculation
"""

import numpy as np
import pytest

from src.network.generator import create_cathode_network
from src.post.analysis import (
    coordination_number,
    discharge_curve,
    pore_size_distribution,
    state_of_lithiation,
)


class TestStructuralAnalysis:
    """Tests for network structural analysis."""

    def test_pore_size_distribution(self):
        """PSD should return bins and counts."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        bins, hist = pore_size_distribution(net)
        assert len(bins) > 0
        assert len(hist) == len(bins) - 1
        assert np.all(hist >= 0)
        assert np.sum(hist) == net.Np

    def test_coordination_number(self):
        """Coordination number should be positive and reasonable."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        cn = coordination_number(net)
        assert cn.shape == (net.Np,)
        assert np.all(cn >= 0)
        # Cubic network: max coordination = 6
        assert np.max(cn) <= 6
        # Average should be reasonable (3-6 for cubic)
        assert 2 < np.mean(cn) < 7


class TestStateOfLithiation:
    """Tests for state-of-lithiation calculation."""

    def test_sol_uniform_initial(self):
        """SoL should be uniform at initial concentration."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        c_s = np.ones(net.Np) * 24450.0
        sol = state_of_lithiation(c_s, net)
        nmc_mask = net["pore.nmc"]
        sol_nmc = sol[nmc_mask]
        # All NMC pores should have SoL = 0.5
        np.testing.assert_allclose(sol_nmc, 0.5, atol=1e-10)

    def test_sol_range(self):
        """SoL should be in [0, 1] for NMC pores."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        c_s = np.random.default_rng(42).uniform(0, 48900, net.Np)
        sol = state_of_lithiation(c_s, net)
        nmc_mask = net["pore.nmc"]
        assert np.all(sol[nmc_mask] >= 0)
        assert np.all(sol[nmc_mask] <= 1)


class TestDischargeCurve:
    """Tests for discharge curve computation."""

    def test_discharge_curve_returns_voltage_capacity(self):
        """Should return voltage and capacity arrays."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        result = discharge_curve(
            net, I_app=-0.001, dt=10.0, n_steps=5, T=298.15,
        )
        assert "voltage" in result
        assert "capacity" in result
        assert len(result["voltage"]) == 5
        assert len(result["capacity"]) == 5

    def test_discharge_curve_accepts_positive_c_rate(self):
        """Positive discharge C-rate should be converted to signed I_app < 0."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        result = discharge_curve(
            net, I_app=None, C_rate=0.2, dt=10.0, n_steps=2, T=298.15,
        )
        assert result["I_app"] < 0.0
        assert result["C_rate"] == 0.2
        np.testing.assert_allclose(
            result["capacity"],
            abs(result["I_app"]) * result["time"],
        )
        np.testing.assert_allclose(result["capacity_Ah_m2"], result["capacity"] / 3600.0)

    def test_discharge_voltage_in_range(self):
        """All voltages should be in physical range."""
        net = create_cathode_network(
            shape=[5, 5, 5], spacing=1e-5, porosity=0.5, seed=42,
        )
        result = discharge_curve(
            net, I_app=-0.001, dt=10.0, n_steps=3, T=298.15,
        )
        assert np.all(result["voltage"] > 2.0)
        assert np.all(result["voltage"] < 5.0)
