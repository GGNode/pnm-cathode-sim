"""Butler-Volmer kinetics for electrode reaction.

Implements the Butler-Volmer equation for the intercalation reaction
at the electrolyte/NMC interface:

    I_rxn = i0 * [exp(alpha_a * F * eta / (R*T)) - exp(-alpha_c * F * eta / (R*T))]

where eta = phi_s - phi_e - U_eq is the overpotential.

References:
    Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np

# Physical constants
F = 96485.3329   # Faraday constant [C/mol]
R = 8.314462     # Gas constant [J/(mol·K)]


def butler_volmer(
    i0: float,
    eta: float | np.ndarray,
    T: float = 298.15,
    alpha_a: float = 0.5,
    alpha_c: float = 0.5,
) -> float | np.ndarray:
    """Compute Butler-Volmer current density.

    Parameters
    ----------
    i0 : float
        Exchange current density [A/m²].
    eta : float or array
        Overpotential [V]. Positive = anodic (discharge).
    T : float
        Temperature [K].
    alpha_a : float
        Anodic transfer coefficient.
    alpha_c : float
        Cathodic transfer coefficient.

    Returns
    -------
    I_rxn : float or array
        Reaction current density [A/m²].
    """
    eta = np.asarray(eta, dtype=float)
    f = F / (R * T)

    # Clip exponent arguments to prevent overflow
    # exp(500) ≈ 1.4e217 is safe for float64; use 500 as limit
    arg_a = np.clip(alpha_a * f * eta, -500.0, 500.0)
    arg_c = np.clip(-alpha_c * f * eta, -500.0, 500.0)

    result = i0 * (np.exp(arg_a) - np.exp(arg_c))

    # Return scalar if input was scalar
    return float(result) if result.ndim == 0 else result


def exchange_current_density(
    k0: float,
    ce: float,
    cs: float,
    cs_max: float = 48900.0,
    alpha_a: float = 0.5,
    alpha_c: float = 0.5,
) -> float:
    """Compute exchange current density for Li intercalation.

    i0 = k0 * F * ce^alpha_a * (cs_max - cs)^alpha_a * cs^alpha_c

    Parameters
    ----------
    k0 : float
        Rate constant [m^(2.5) / (mol^0.5 · s)].
    ce : float
        Electrolyte Li+ concentration [mol/m³].
    cs : float
        Surface solid Li concentration [mol/m³].
    cs_max : float
        Maximum solid Li concentration [mol/m³].
    alpha_a, alpha_c : float
        Transfer coefficients.

    Returns
    -------
    i0 : float
        Exchange current density [A/m²].
    """
    ce = max(ce, 1e-10)  # avoid zero
    cs = np.clip(cs, 1e-10, cs_max - 1e-10)

    return k0 * F * (ce ** alpha_a) * ((cs_max - cs) ** alpha_a) * (cs ** alpha_c)
