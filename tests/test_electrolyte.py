"""Unit tests for electrolyte transport correlations."""

import numpy as np
import pytest

from src.physics.electrolyte import (
    electrolyte_diffusion_coefficient,
    electrolyte_ionic_conductivity,
)


class TestElectrolyteDiffusion:
    """D_e(c_e, T) should return m^2/s."""

    def test_value_at_standard_conditions(self):
        """At c_e=1200 mol/m^3, T=303K: D_e ≈ 3.25e-10 m^2/s."""
        De = electrolyte_diffusion_coefficient(1200.0, 303.0)
        assert 1e-10 < De < 1e-9, f"D_e = {De:.4e}, expected ~3.25e-10 m^2/s"
        np.testing.assert_allclose(De, 3.25e-10, rtol=0.05)

    def test_units_are_m2_not_cm2(self):
        """Raw formula gives ~3.25e-6 (cm^2/s). Must convert to m^2/s."""
        De = electrolyte_diffusion_coefficient(1200.0, 303.0)
        # If someone forgets the conversion, De would be ~3.25e-6
        assert De < 1e-8, (
            f"D_e = {De:.4e} — looks like cm^2/s, not m^2/s. "
            "Missing 1e-4 conversion?"
        )

    def test_array_input(self):
        """Should handle array input."""
        De = electrolyte_diffusion_coefficient(
            np.array([600.0, 1200.0, 2000.0]), 303.0,
        )
        assert De.shape == (3,)
        assert np.all(De > 0)


class TestElectrolyteConductivity:
    """κ(c_e, T) should return S/m."""

    def test_value_at_standard_conditions(self):
        """At c_e=1200 mol/m^3, T=303K: κ ≈ 1.28 S/m."""
        kappa = electrolyte_ionic_conductivity(1200.0, 303.0)
        assert 0.5 < kappa < 5.0, f"κ = {kappa:.4f}, expected ~1.28 S/m"
        np.testing.assert_allclose(kappa, 1.284, rtol=0.05)

    def test_units_are_Sm_not_mS_cm(self):
        """Raw formula gives ~12.84 (mS/cm). Must convert to S/m."""
        kappa = electrolyte_ionic_conductivity(1200.0, 303.0)
        # If someone forgets the conversion, kappa would be ~12.84
        assert kappa < 5.0, (
            f"κ = {kappa:.4f} — looks like mS/cm, not S/m. "
            "Missing ×0.1 conversion?"
        )

    def test_positive_conductivity(self):
        """Conductivity must be positive for all reasonable concentrations."""
        for ce in [100.0, 500.0, 1200.0, 2000.0, 3000.0]:
            kappa = electrolyte_ionic_conductivity(ce, 303.0)
            assert kappa > 0, f"κ({ce}) = {kappa} ≤ 0"

    def test_array_input(self):
        """Should handle array input."""
        kappa = electrolyte_ionic_conductivity(
            np.array([600.0, 1200.0, 2000.0]), 303.0,
        )
        assert kappa.shape == (3,)
        assert np.all(kappa > 0)
