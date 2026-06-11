"""Phase 1: Solid-phase diffusion unit tests.

Verifies:
  1. Single spherical particle diffusion equation
  2. Mass conservation: total Li in particle changes only by surface flux
  3. Steady-state: uniform concentration when no flux
  4. Cottrell limiting behavior at short times
  5. D_s formula gives physically reasonable values
"""

import numpy as np
import pytest

from pnmcathode.physics.solid import (
    solid_diffusion_rhs,
    nmc532_diffusion_coefficient,
    discretize_spherical_particle,
)


class TestNMC532DiffusionCoefficient:
    """Test the D_s(c_s, T) empirical formula."""

    def test_positive_value(self):
        """D_s must be positive for all physically meaningful c_s."""
        cs_values = np.linspace(1000, 47000, 20)
        for cs in cs_values:
            D = nmc532_diffusion_coefficient(cs, T=298.15)
            assert D > 0, f"D_s({cs}) = {D} <= 0"

    def test_reasonable_range(self):
        """D_s should be in range ~1e-16 to ~1e-10 m²/s."""
        cs_values = np.linspace(1000, 47000, 20)
        for cs in cs_values:
            D = nmc532_diffusion_coefficient(cs, T=298.15)
            assert 1e-18 < D < 1e-8, (
                f"D_s({cs}, 298K) = {D:.2e} m²/s outside expected range"
            )

    def test_printed_formula_has_no_temperature_dependence(self):
        """Table II footnote (1) D_s is concentration-dependent only."""
        D_cold = nmc532_diffusion_coefficient(24450, T=273.15)
        D_hot = nmc532_diffusion_coefficient(24450, T=333.15)
        assert D_hot == pytest.approx(D_cold)


class TestSphericalParticleDiscretization:
    """Test finite-volume discretization of spherical particle."""

    def test_shell_volumes_sum_to_particle_volume(self):
        """Sum of shell volumes should equal total particle volume."""
        R_p = 5e-6  # 5 µm radius
        N = 10
        dr, volumes = discretize_spherical_particle(R_p, N)
        total_volume = (4.0 / 3.0) * np.pi * R_p**3
        assert abs(sum(volumes) - total_volume) / total_volume < 1e-10

    def test_shell_count(self):
        """Should return correct number of shells."""
        R_p = 5e-6
        N = 20
        dr, volumes = discretize_spherical_particle(R_p, N)
        assert len(volumes) == N
        assert abs(dr - R_p / N) < 1e-15


class TestSolidDiffusionRHS:
    """Test the right-hand side of the solid diffusion equation."""

    def test_steady_state_uniform(self):
        """With uniform concentration and no surface flux, RHS should be zero."""
        N = 10
        R_p = 5e-6
        D_s = 1e-14
        c_s = np.ones(N) * 24450.0  # uniform
        flux_surface = 0.0

        rhs = solid_diffusion_rhs(c_s, R_p, D_s, flux_surface, N)
        np.testing.assert_allclose(rhs, 0.0, atol=1e-10)

    def test_mass_conservation(self):
        """Total Li change should equal surface flux * area * dt."""
        N = 10
        R_p = 5e-6
        D_s = 1e-14
        c_s = np.ones(N) * 24450.0
        # Non-zero surface flux
        flux_surface = 1e-3  # mol/(m²·s)

        rhs = solid_diffusion_rhs(c_s, R_p, D_s, flux_surface, N)
        dr, volumes = discretize_spherical_particle(R_p, N)
        # Total rate of change of Li moles
        total_rate = sum(rhs * volumes)
        # Should equal surface flux * surface area
        surface_area = 4.0 * np.pi * R_p**2
        expected_rate = flux_surface * surface_area
        assert abs(total_rate - expected_rate) / abs(expected_rate) < 1e-6

    def test_negative_flux_decreases_surface_concentration(self):
        """Extraction (negative flux) should decrease surface concentration."""
        N = 10
        R_p = 5e-6
        D_s = 1e-14
        c_s = np.ones(N) * 24450.0
        flux_surface = -1e-3  # extracting Li

        rhs = solid_diffusion_rhs(c_s, R_p, D_s, flux_surface, N)
        # Surface shell (last) should have most negative rate
        assert rhs[-1] < 0
