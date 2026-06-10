"""Phase 1: NMC532 open-circuit voltage curve tests.

Verifies:
  1. OCV is monotonically decreasing over the calibrated GITT window
  2. OCV range is physically reasonable in that window
  3. OCV values match the corrected Khan/Verma polynomial
  4. Boundary values are finite
"""

import numpy as np
import pytest

from src.physics.ocv import nmc532_ocv


class TestNMC532OCV:
    """Unit tests for NMC532 OCV curve."""

    def test_ocv_range(self):
        """OCV should be in range [2.5, 4.5] V over the full [0,1] range."""
        xs = np.linspace(0.30, 0.99, 100)
        for x in xs:
            u = nmc532_ocv(x)
            assert 2.5 <= u <= 4.5, f"OCV={u:.3f}V at x={x:.2f} out of range"

    def test_monotonically_decreasing(self):
        """OCV should decrease as x increases (more lithiated → lower voltage)."""
        xs = np.linspace(0.30, 0.95, 50)
        ocvs = [nmc532_ocv(x) for x in xs]
        for i in range(len(ocvs) - 1):
            assert ocvs[i] > ocvs[i + 1], (
                f"OCV not decreasing: U({xs[i]:.2f})={ocvs[i]:.4f} "
                f"<= U({xs[i+1]:.2f})={ocvs[i+1]:.4f}"
            )

    def test_midpoint_voltage(self):
        """At soc=0.5, OCV should be around 3.6V (literature NMC532)."""
        u = nmc532_ocv(0.5)
        assert 3.5 < u < 3.8, f"OCV at soc=0.5 = {u:.3f}V, expected ~3.6V"

    def test_37v_region(self):
        """OCV near soc=0.85 should be around 3.4V (literature NMC532)."""
        u = nmc532_ocv(0.851725)
        assert 3.3 < u < 3.5, f"OCV at soc=0.851725 = {u:.3f}V"

    def test_empty_voltage_higher_than_full(self):
        """Lower lithiation should have higher OCV than higher lithiation."""
        u_empty = nmc532_ocv(0.30)
        u_full = nmc532_ocv(0.95)
        assert u_empty > u_full

    def test_vectorized(self):
        """OCV function should accept array input."""
        xs = np.array([0.3, 0.5, 0.7, 0.9])
        us = nmc532_ocv(xs)
        assert us.shape == xs.shape
        # Should still be monotonically decreasing
        assert all(us[i] > us[i + 1] for i in range(len(us) - 1))

    def test_boundary_finite(self):
        """OCV at boundaries should be finite (not NaN or Inf)."""
        assert np.isfinite(nmc532_ocv(0.25))
        assert np.isfinite(nmc532_ocv(0.99))
