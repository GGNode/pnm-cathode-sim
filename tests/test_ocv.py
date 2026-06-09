"""Phase 1: NMC532 open-circuit voltage curve tests.

Verifies:
  1. OCV is monotonically decreasing with x (lithiation degree)
  2. OCV range is physically reasonable (3.0–4.3 V for NMC532)
  3. OCV at x=0.5 is approximately 3.7 V
  4. Boundary values are finite
"""

import numpy as np
import pytest

from src.physics.ocv import nmc532_ocv


class TestNMC532OCV:
    """Unit tests for NMC532 OCV curve."""

    def test_ocv_range(self):
        """OCV should be in range [3.0, 4.3] V for x in [0.01, 0.99]."""
        xs = np.linspace(0.01, 0.99, 100)
        for x in xs:
            u = nmc532_ocv(x)
            assert 3.0 <= u <= 4.3, f"OCV={u:.3f}V at x={x:.2f} out of range"

    def test_monotonically_decreasing(self):
        """OCV should decrease as x increases (more lithiated → lower voltage)."""
        xs = np.linspace(0.05, 0.95, 50)
        ocvs = [nmc532_ocv(x) for x in xs]
        for i in range(len(ocvs) - 1):
            assert ocvs[i] > ocvs[i + 1], (
                f"OCV not decreasing: U({xs[i]:.2f})={ocvs[i]:.4f} "
                f"<= U({xs[i+1]:.2f})={ocvs[i+1]:.4f}"
            )

    def test_midpoint_voltage(self):
        """At x=0.5, OCV should be approximately 3.7 V for NMC532."""
        u = nmc532_ocv(0.5)
        assert 3.5 < u < 3.9, f"OCV at x=0.5 = {u:.3f}V, expected ~3.7V"

    def test_empty_voltage_higher_than_full(self):
        """Empty cathode (x→0) should have higher OCV than full (x→1)."""
        u_empty = nmc532_ocv(0.05)
        u_full = nmc532_ocv(0.95)
        assert u_empty > u_full

    def test_vectorized(self):
        """OCV function should accept array input."""
        xs = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        us = nmc532_ocv(xs)
        assert us.shape == xs.shape
        # Should still be monotonically decreasing
        assert all(us[i] > us[i + 1] for i in range(len(us) - 1))

    def test_boundary_finite(self):
        """OCV at boundaries should be finite (not NaN or Inf)."""
        assert np.isfinite(nmc532_ocv(0.01))
        assert np.isfinite(nmc532_ocv(0.99))
