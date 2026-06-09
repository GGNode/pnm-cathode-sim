"""Visualization tools for pore network battery model.

Provides 2D/3D plotting of network properties, concentration fields,
and discharge curves.

References:
    Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def plot_discharge_curve(result: dict, ax=None, **kwargs):
    """Plot voltage vs. capacity from discharge_curve() output.

    Parameters
    ----------
    result : dict
        Output from analysis.discharge_curve().
    ax : matplotlib.axes.Axes or None
        Axes to plot on.  If None, creates a new figure.
    **kwargs
        Passed to ax.plot().

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if not HAS_MPL:
        raise ImportError("matplotlib is required for plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(result["capacity"], result["voltage"], **kwargs)
    ax.set_xlabel("Capacity (A·s/m²)")
    ax.set_ylabel("Cell Voltage (V)")
    ax.set_title("Discharge Curve")
    ax.grid(True, alpha=0.3)
    return ax


def plot_phase_map(
    net: op.network.Cubic,
    property_name: str = "pore.phase_label",
    ax=None,
):
    """Plot 2D projection of network colored by a pore property.

    Parameters
    ----------
    net : openpnm.network.Cubic
        Pore network.
    property_name : str
        Name of pore property to color by.
    ax : matplotlib.axes.Axes or None
        Axes to plot on.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if not HAS_MPL:
        raise ImportError("matplotlib is required for plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    coords = net["pore.coords"]
    values = net[property_name]

    sc = ax.scatter(coords[:, 0], coords[:, 1], c=values, cmap="viridis", s=20)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"{property_name}")
    ax.set_aspect("equal")
    plt.colorbar(sc, ax=ax)
    return ax


def plot_concentration_field(
    net: op.network.Cubic,
    c: np.ndarray,
    title: str = "Concentration",
    ax=None,
):
    """Plot 2D projection of a concentration field.

    Parameters
    ----------
    net : openpnm.network.Cubic
        Pore network.
    c : ndarray
        Concentration values per pore (Np).
    title : str
        Plot title.
    ax : matplotlib.axes.Axes or None

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if not HAS_MPL:
        raise ImportError("matplotlib is required for plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    coords = net["pore.coords"]
    valid = ~np.isnan(c)

    sc = ax.scatter(
        coords[valid, 0], coords[valid, 1],
        c=c[valid], cmap="RdYlBu_r", s=20,
    )
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(title)
    ax.set_aspect("equal")
    plt.colorbar(sc, ax=ax, label="mol/m³")
    return ax
