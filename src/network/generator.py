"""Pore network generation for cathode electrode.

Creates a cubic pore network with three-phase labeling:
  - Phase 0: Electrolyte (pore space filled with electrolyte)
  - Phase 1: NMC532 (active material)
  - Phase 2: CBD (carbon binder domain)

The network is a regular cubic lattice where each node represents a pore body
and each bond represents a pore throat. Phase labels are assigned stochastically
to approximate the target porosity.

For the Khan et al. (2021) paper:
  - 1CAL: 4637 nodes, 31427 bonds
  - 3CAL: 3510 nodes, 23126 bonds
  - Porosity: ~35%

References:
    Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op


def create_cathode_network(
    shape: list[int] = [10, 10, 10],
    spacing: float = 1e-5,
    porosity: float = 0.35,
    cbd_fraction: float = 0.10,
    seed: int | None = None,
) -> op.network.Cubic:
    """Create a cubic pore network with three-phase cathode labeling.

    Parameters
    ----------
    shape : list of int
        Network dimensions [nx, ny, nz].
    spacing : float
        Pore-to-pore spacing [m].
    porosity : float
        Target electrolyte volume fraction (phase 0).
    cbd_fraction : float
        CBD volume fraction (phase 2). NMC fraction = 1 - porosity - cbd_fraction.
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    net : openpnm.network.Cubic
        Network with pore properties assigned.
    """
    rng = np.random.default_rng(seed)

    # Create cubic network
    net = op.network.Cubic(shape=shape, spacing=spacing)

    # --- Assign phase labels ---
    # Phase 0: electrolyte, Phase 1: NMC532, Phase 2: CBD
    n_pores = net.Np
    labels = np.zeros(n_pores, dtype=int)

    # Random assignment based on volume fractions
    rand_vals = rng.random(n_pores)
    labels[rand_vals < porosity] = 0                        # electrolyte
    labels[(rand_vals >= porosity) & (rand_vals < porosity + cbd_fraction)] = 2  # CBD
    labels[rand_vals >= porosity + cbd_fraction] = 1         # NMC532

    net["pore.phase_label"] = labels
    net["pore.electrolyte"] = labels == 0
    net["pore.nmc"] = labels == 1
    net["pore.cbd"] = labels == 2

    # --- Assign geometric properties ---
    # Pore diameter: uniform distribution around spacing
    pore_diameter = rng.uniform(0.5, 1.5, n_pores) * spacing
    net["pore.diameter"] = pore_diameter
    net["pore.volume"] = (np.pi / 6.0) * pore_diameter**3

    # Throat properties
    n_throats = net.Nt
    throat_diameter = rng.uniform(0.2, 0.8, n_throats) * spacing * 0.5
    net["throat.diameter"] = throat_diameter
    net["throat.area"] = (np.pi / 4.0) * throat_diameter**2

    # Throat length: distance between pore centers minus pore radii
    conns = net["throat.conns"]
    p1_coords = net["pore.coords"][conns[:, 0]]
    p2_coords = net["pore.coords"][conns[:, 1]]
    p1_radius = pore_diameter[conns[:, 0]] / 2.0
    p2_radius = pore_diameter[conns[:, 1]] / 2.0
    center_dist = np.linalg.norm(p2_coords - p1_coords, axis=1)
    throat_length = center_dist - p1_radius - p2_radius
    # Ensure minimum length
    throat_length = np.maximum(throat_length, spacing * 0.01)
    net["throat.length"] = throat_length

    # --- Assign phase-specific properties ---
    # NMC properties
    net["pore.cs_max"] = np.where(labels == 1, 48900.0, 0.0)
    net["pore.sigma_nmc"] = np.where(labels == 1, 0.01, 0.0)  # S/m
    net["pore.sigma_cbd"] = np.where(labels == 2, 760.0, 0.0)  # S/m

    return net
