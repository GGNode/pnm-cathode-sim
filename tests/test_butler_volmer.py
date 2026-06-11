"""Phase 1: Butler-Volmer kinetics unit tests.

Verifies:
  1. BV current matches analytic formula for given overpotential
  2. Symmetry: i(eta) = -i(-eta) for symmetric transfer coefficients
  3. Zero overpotential → zero current
  4. Tafel limiting behavior at large |eta|
  5. Numerical stability at large overpotentials (no overflow)
"""

import numpy as np
import pytest

from pnmcathode.physics.reaction import butler_volmer, exchange_current_density

F = 96485.3329
R = 8.314462
T = 298.15


class TestButlerVolmer:
    """Unit tests for Butler-Volmer equation."""

    def test_zero_overpotential_zero_current(self):
        """At eta=0, BV current must be zero."""
        i0 = 10.0  # A/m²
        result = butler_volmer(i0=i0, eta=0.0, T=T)
        assert abs(result) < 1e-15

    def test_symmetry(self):
        """For alpha_a == alpha_c == 0.5, i(eta) == -i(-eta)."""
        i0 = 10.0
        eta = 0.05  # V
        i_pos = butler_volmer(i0=i0, eta=eta, T=T)
        i_neg = butler_volmer(i0=i0, eta=-eta, T=T)
        assert abs(i_pos + i_neg) < 1e-10

    def test_analytic_formula(self):
        """BV current must match: i = i0 * [exp(alpha_a*F*eta/(RT)) - exp(-alpha_c*F*eta/(RT))]."""
        i0 = 5.0
        eta = 0.02
        alpha_a = 0.5
        alpha_c = 0.5

        expected = i0 * (
            np.exp(alpha_a * F * eta / (R * T))
            - np.exp(-alpha_c * F * eta / (R * T))
        )
        result = butler_volmer(i0=i0, eta=eta, T=T, alpha_a=alpha_a, alpha_c=alpha_c)
        assert abs(result - expected) < 1e-8

    def test_positive_overpotential_positive_current(self):
        """Positive overpotential (anodic) should give positive current."""
        i0 = 10.0
        result = butler_volmer(i0=i0, eta=0.1, T=T)
        assert result > 0

    def test_negative_overpotential_negative_current(self):
        """Negative overpotential (cathodic) should give negative current."""
        i0 = 10.0
        result = butler_volmer(i0=i0, eta=-0.1, T=T)
        assert result < 0

    def test_large_overpotential_no_overflow(self):
        """Large overpotentials must not cause numerical overflow."""
        i0 = 10.0
        # eta = 1V → exp(0.5*96485*1/(8.314*298)) ≈ exp(19.5) ≈ 3e8 — should be fine
        result = butler_volmer(i0=i0, eta=1.0, T=T)
        assert np.isfinite(result)
        # Even larger
        result = butler_volmer(i0=i0, eta=5.0, T=T)
        assert np.isfinite(result)

    def test_tafel_limit_large_positive_eta(self):
        """At large positive eta, BV ≈ i0 * exp(alpha_a*F*eta/(RT)) (Tafel)."""
        i0 = 10.0
        eta = 0.3  # Large enough for Tafel approximation
        bv = butler_volmer(i0=i0, eta=eta, T=T)
        tafel = i0 * np.exp(0.5 * F * eta / (R * T))
        # Should be within 1% at eta=0.3V
        assert abs(bv - tafel) / tafel < 0.01

    def test_vectorized_eta(self):
        """BV function should accept array of eta values."""
        i0 = 10.0
        etas = np.array([-0.1, -0.05, 0.0, 0.05, 0.1])
        results = butler_volmer(i0=i0, eta=etas, T=T)
        assert results.shape == etas.shape
        # Zero at eta=0
        assert abs(results[2]) < 1e-15
        # Symmetry
        assert abs(results[0] + results[4]) < 1e-10
        assert abs(results[1] + results[3]) < 1e-10


class TestExchangeCurrentDensity:
    """Unit tests for exchange current density."""

    def test_basic_calculation(self):
        """i0 = k0 * F * ce^alpha_a * (cs_max - cs)^alpha_a * cs^alpha_c."""
        k0 = 1e-10  # rate constant
        ce = 1200.0  # mol/m³
        cs = 24450.0  # mol/m³ (half lithiated)
        cs_max = 48900.0  # mol/m³
        alpha_a = 0.5
        alpha_c = 0.5

        expected = k0 * F * (ce ** alpha_a) * ((cs_max - cs) ** alpha_a) * (cs ** alpha_c)
        result = exchange_current_density(
            k0=k0, ce=ce, cs=cs, cs_max=cs_max,
            alpha_a=alpha_a, alpha_c=alpha_c,
        )
        assert abs(result - expected) < 1e-10

    def test_positive_value(self):
        """Exchange current density must always be positive."""
        result = exchange_current_density(
            k0=1e-10, ce=1200.0, cs=24450.0, cs_max=48900.0,
        )
        assert result > 0
