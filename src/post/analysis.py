"""Post-processing and analysis tools for pore network battery model.

Provides structural analysis of the pore network and performance
analysis of discharge simulations.

References:
    Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op

from src.solver.transient import TransientSolver


def pore_size_distribution(
    net: op.network.Cubic, n_bins: int = 20,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute pore size distribution histogram.

    Parameters
    ----------
    net : openpnm.network.Cubic
        Pore network.
    n_bins : int
        Number of histogram bins.

    Returns
    -------
    bins : ndarray
        Bin edges (length n_bins + 1).
    hist : ndarray
        Counts per bin (length n_bins).
    """
    diameters = net["pore.diameter"]
    hist, bins = np.histogram(diameters, bins=n_bins)
    return bins, hist


def coordination_number(net: op.network.Cubic) -> np.ndarray:
    """Compute coordination number (number of connections) per pore.

    Parameters
    ----------
    net : openpnm.network.Cubic
        Pore network.

    Returns
    -------
    cn : ndarray of shape (Np,)
        Coordination number for each pore.
    """
    conns = net["throat.conns"]
    cn = np.zeros(net.Np, dtype=int)
    np.add.at(cn, conns[:, 0], 1)
    np.add.at(cn, conns[:, 1], 1)
    return cn


def state_of_lithiation(c_s: np.ndarray, net: op.network.Cubic) -> np.ndarray:
    """Compute state of lithiation x = c_s / c_s_max for each pore.

    Parameters
    ----------
    c_s : ndarray
        Solid Li concentration [mol/m³] for all pores.
    net : openpnm.network.Cubic
        Pore network.

    Returns
    -------
    sol : ndarray
        State of lithiation in [0, 1].  Non-NMC pores get NaN.
    """
    cs_max = 48900.0
    sol = np.full(net.Np, np.nan)
    nmc_mask = net["pore.nmc"]
    sol[nmc_mask] = np.clip(c_s[nmc_mask] / cs_max, 0.0, 1.0)
    return sol


def discharge_curve(
    net: op.network.Cubic,
    I_app: float = -0.001,
    dt: float = 10.0,
    n_steps: int = 100,
    T: float = 298.15,
    k0: float = 5e-10,
    c_e_init: float = 1200.0,
    c_s_init: float = 24450.0,
) -> dict:
    """Run a transient discharge and return voltage vs. capacity.

    Parameters
    ----------
    net : openpnm.network.Cubic
        Pore network.
    I_app : float
        Applied current density [A/m²].
    dt : float
        Time step [s].
    n_steps : int
        Number of time steps.
    T : float
        Temperature [K].
    k0 : float
        Rate constant.
    c_e_init, c_s_init : float
        Initial concentrations.

    Returns
    -------
    result : dict
        Keys: 'voltage' (n_steps array), 'capacity' (n_steps array, A·s/m²),
        'time' (n_steps array, s).
    """
    solver = TransientSolver(net, T=T, k0=k0)
    solver.set_concentration(c_e=c_e_init, c_s=c_s_init)

    voltages = np.zeros(n_steps)
    times = np.zeros(n_steps)
    capacities = np.zeros(n_steps)

    for i in range(n_steps):
        result = solver.step(dt=dt, I_app=I_app)
        voltages[i] = result["voltage"]
        times[i] = (i + 1) * dt
        # Capacity = I * t (A·s/m²)
        capacities[i] = abs(I_app) * times[i]

    return {
        "voltage": voltages,
        "capacity": capacities,
        "time": times,
    }
