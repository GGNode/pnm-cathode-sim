"""Electrolyte transport correlations from Khan et al. (2021)."""

from __future__ import annotations

import numpy as np


def _mol_m3_to_mol_l(c_e: float | np.ndarray) -> np.ndarray:
    """Convert electrolyte concentration from mol/m^3 to mol/L for Table II."""
    return np.asarray(c_e, dtype=float) / 1000.0


def electrolyte_diffusion_coefficient(
    c_e: float | np.ndarray,
    T: float = 303.0,
) -> float | np.ndarray:
    """Li+ diffusivity in electrolyte from Table II, in m^2/s.

    Table II prints:
        D_Li+ = 10^(-4.43 - 54.0/(T - 229 - 5c2) - 0.22c2)

    The concentration ``c2`` is used in mol/L in this empirical correlation,
    while the model state stores electrolyte concentration in mol/m^3.
    """
    c2 = _mol_m3_to_mol_l(c_e)
    exponent = -4.43 - 54.0 / (T - 229.0 - 5.0 * c2) - 0.22 * c2
    value = 10.0 ** exponent
    return float(value) if value.ndim == 0 else value


def electrolyte_ionic_conductivity(
    c_e: float | np.ndarray,
    T: float = 303.0,
) -> float | np.ndarray:
    """Electrolyte ionic conductivity correlation, in S/m.

    This is Table II footnote (2), which is printed under the AM diffusivity
    row but is the standard electrolyte conductivity expression.
    """
    c2 = _mol_m3_to_mol_l(c_e)
    inner = (
        -10.5
        + 0.0740 * T
        - 6.96e-5 * T**2
        + 0.668 * c2
        - 0.0178 * c2 * T
        + 2.80e-5 * c2 * T**2
        + 0.494 * c2**2
        - 8.86e-4 * c2**2 * T
    )
    value = c2 * inner**2
    return float(value) if value.ndim == 0 else value
