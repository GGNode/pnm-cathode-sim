"""Solid-phase Li diffusion in NMC532 particles.

Implements Fick's law for Li diffusion in a spherical active material particle,
discretized using finite-volume method with uniform shells.

    ∂c_s/∂t = D_s * ∇²c_s

In spherical coordinates:
    ∂c_s/∂t = (1/r²) * ∂/∂r (r² * D_s * ∂c_s/∂r)

The NMC532 diffusion coefficient D_s depends on c_s and T
(see nmc532_diffusion_coefficient).

References:
    Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np


def nmc532_diffusion_coefficient(c_s: float, T: float = 298.15) -> float:
    """NMC532 Li diffusion coefficient D_s(c_s, T).

    Empirical formula from the literature for NMC532.
    The formula in Khan et al. is complex; here we use a simpler but
    physically representative expression that captures the key behavior:
    D_s varies with state of lithiation and temperature.

    D_s = D_ref * exp(-E_a/R * (1/T - 1/T_ref)) * f(c_s)

    Parameters
    ----------
    c_s : float
        Solid Li concentration [mol/m³].
    T : float
        Temperature [K].

    Returns
    -------
    D_s : float
        Diffusion coefficient [m²/s].
    """
    R = 8.314462  # J/(mol·K)

    # Reference values at 298.15 K
    D_ref = 1e-14        # m²/s — typical for NMC at 50% SoC
    E_a = 30000.0        # J/mol — activation energy
    T_ref = 298.15       # K

    # Concentration dependence: D_s has a minimum near x=0 and x=1
    # Use a smooth function that is ~1 at x=0.5 and decreases at extremes.
    # c_s_max for NMC532 = 48900 mol/m³
    c_s_max = 48900.0
    x = np.clip(c_s / c_s_max, 0.01, 0.99)

    # Empirical fit: D_s(x) has a U-shape or shallow dependence
    # Literature shows D_s ~ 1e-14 to 1e-12 for NMC
    f_cs = 0.1 + 0.9 * (4.0 * x * (1.0 - x))  # peaks at x=0.5, ~0.1 at extremes

    D_s = D_ref * np.exp(-E_a / R * (1.0 / T - 1.0 / T_ref)) * f_cs

    return D_s


def discretize_spherical_particle(R_p: float, N: int) -> tuple[float, np.ndarray]:
    """Discretize a spherical particle into N concentric shells (finite volume).

    Parameters
    ----------
    R_p : float
        Particle radius [m].
    N : int
        Number of shells.

    Returns
    -------
    dr : float
        Shell thickness [m].
    volumes : ndarray
        Volume of each shell [m³].
    """
    dr = R_p / N
    # Shell radii: r_i = (i + 0.5) * dr for i = 0..N-1
    r_inner = np.arange(N) * dr
    r_outer = r_inner + dr
    # Volume of spherical shell: V = (4/3)*pi*(r_outer³ - r_inner³)
    volumes = (4.0 / 3.0) * np.pi * (r_outer**3 - r_inner**3)
    return dr, volumes


def solid_diffusion_rhs(
    c_s: np.ndarray,
    R_p: float,
    D_s: float,
    flux_surface: float,
    N: int,
) -> np.ndarray:
    """Right-hand side of the solid diffusion ODE system.

    d(c_s)/dt = RHS, where RHS is computed from finite-volume discretization
    of spherical diffusion with flux boundary condition at the surface.

    Parameters
    ----------
    c_s : ndarray of shape (N,)
        Li concentration in each shell [mol/m³].
    R_p : float
        Particle radius [m].
    D_s : float
        Diffusion coefficient [m²/s] (assumed uniform in particle).
    flux_surface : float
        Li flux at particle surface [mol/(m²·s)].
        Positive = into particle (lithiation), negative = extraction.
    N : int
        Number of shells.

    Returns
    -------
    dc_dt : ndarray of shape (N,)
        Time derivative of concentration in each shell [mol/(m³·s)].
    """
    dr = R_p / N
    dc_dt = np.zeros(N)

    # Interior faces (between shell i and shell i+1)
    for i in range(N - 1):
        # Face at r = (i+1) * dr
        r_face = (i + 1) * dr
        # Concentration gradient at face
        dcdr = (c_s[i + 1] - c_s[i]) / dr
        # Flux through face: J = -D * dc/dr  [mol/(m²·s)]
        flux = -D_s * dcdr
        # Area of face: A = 4*pi*r²
        area = 4.0 * np.pi * r_face**2
        # Volume of shell i
        r_i_inner = i * dr
        r_i_outer = (i + 1) * dr
        vol_i = (4.0 / 3.0) * np.pi * (r_i_outer**3 - r_i_inner**3)
        # Volume of shell i+1
        r_ip1_inner = (i + 1) * dr
        r_ip1_outer = (i + 2) * dr
        vol_ip1 = (4.0 / 3.0) * np.pi * (r_ip1_outer**3 - r_ip1_inner**3)

        dc_dt[i] += flux * area / vol_i
        dc_dt[i + 1] -= flux * area / vol_ip1

    # Surface boundary (flux into outermost shell)
    r_surface = R_p
    area_surface = 4.0 * np.pi * r_surface**2
    vol_last = (4.0 / 3.0) * np.pi * (R_p**3 - (R_p - dr)**3)
    dc_dt[N - 1] += flux_surface * area_surface / vol_last

    return dc_dt
