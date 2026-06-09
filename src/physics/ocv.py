"""NMC532 open-circuit voltage (OCV) curve.

Empirical OCV(U_eq) vs. lithiation degree x = c_s / c_s_max for NMC532.

Polynomial fit (6th order) to published NMC532 half-cell data:
    U_eq(x) = a0 + a1*x + a2*x² + a3*x³ + a4*x⁴ + a5*x⁵ + a6*x⁶

where x = c_s / c_s_max (0 = empty/delithiated, 1 = fully lithiated).

The voltage decreases with x: high x → low voltage (more discharged).
Typical range: ~4.2V at x=0 to ~3.6V at x=1.

References:
    Bernardi & Gojkovic, J. Electrochem. Soc. 166(12), A2485 (2019).
    Used in Khan et al. (2021) PNM model.
"""

import numpy as np

# Polynomial coefficients for NMC532 OCV(x), x = c_s/c_s_max
# Fitted to published NMC532 half-cell data points:
#   x:  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0
#   U: 4.28 4.15 4.05 3.93 3.82 3.73 3.65 3.60 3.57 3.55 3.53
# Gives: U(0)≈4.28V, U(0.5)≈3.73V, U(1)≈3.53V
_OCV_COEFFS = [
    6.290850,    # a0
    -21.596908,  # a1
    27.066679,   # a2
    -14.688871,  # a3
    3.689029,    # a4
    -1.509975,   # a5
    4.279404,    # a6
]


def nmc532_ocv(x: float | np.ndarray) -> float | np.ndarray:
    """Compute NMC532 open-circuit voltage.

    Parameters
    ----------
    x : float or array
        Lithiation degree, x = c_s / c_s_max, in [0, 1].
        x=0: empty (delithiated), x=1: full (lithiated).

    Returns
    -------
    U_eq : float or array
        Equilibrium potential [V] vs. Li/Li+.
    """
    x = np.asarray(x, dtype=float)
    x = np.clip(x, 1e-6, 1.0 - 1e-6)

    # Evaluate 6th-order polynomial
    U = np.polyval(_OCV_COEFFS, x)

    # Clamp to physical range
    U = np.clip(U, 3.0, 4.3)

    return float(U) if U.ndim == 0 else U


def ocv_derivative(x: float | np.ndarray) -> float | np.ndarray:
    """Compute dU/dx for Newton-Raphson coupling.

    Parameters
    ----------
    x : float or array
        Lithiation degree.

    Returns
    -------
    dUdx : float or array
        Derivative of OCV w.r.t. x [V].
    """
    x = np.asarray(x, dtype=float)
    x = np.clip(x, 1e-6, 1.0 - 1e-6)

    # Derivative of polynomial
    deriv_coeffs = np.polyder(_OCV_COEFFS)
    dU = np.polyval(deriv_coeffs, x)

    return float(dU) if dU.ndim == 0 else dU
