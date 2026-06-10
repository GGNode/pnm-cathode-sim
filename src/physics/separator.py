"""
Collapsed 1D separator boundary model.

Represents the separator as a finite-volume slab with quasi-steady transport,
collapsed into an impedance/transport boundary condition for the cathode
electrolyte.  This captures the dominant high-rate effects (ohmic drop, salt
depletion, and optional Li-foil interfacial polarisation) without adding full
coupled separator grid nodes to the Newton unknowns.

Paper constants from Khan et al. (2021) Table II unless noted otherwise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from src.physics.electrolyte import (
    electrolyte_diffusion_coefficient,
    electrolyte_ionic_conductivity,
)
from src.physics.reaction import F, R


@dataclass(frozen=True)
class SeparatorParams:
    """Parameters for the collapsed 1-D separator boundary model."""

    enabled: bool = False
    thickness: float = 25e-6        # L_sep [m]
    porosity: float = 0.39          # epsilon_sep
    bruggeman: float = 1.5          # Bruggeman exponent
    t_plus: float = 0.363           # Li+ transference number
    c_ref: float = 1200.0           # reference concentration [mol/m3]
    c_floor: float = 1.0            # minimum allowed concentration [mol/m3]
    include_concentration_overpotential: bool = True
    include_li_foil_bv: bool = False
    i0_foil: float = 19.0           # Li foil exchange current density [A/m2]
    alpha_foil: float = 0.5         # Li foil BV symmetry factor


@dataclass(frozen=True)
class SeparatorState:
    """Result of :func:`separator_boundary`."""

    c_cathode: float
    phi_e_cathode: float
    dphi_ohm: float
    dphi_conc: float
    eta_li: float


def separator_boundary(
    I_app: float,
    T: float,
    params: SeparatorParams,
    n_iter: int = 3,
) -> SeparatorState:
    """Compute cathode-side electrolyte BC from separator transport model.

    Parameters
    ----------
    I_app : float
        Applied current density in anodic convention [A/m2].
        Discharge: I_app < 0.
    T : float
        Temperature [K].
    params : SeparatorParams
        Separator model parameters.
    n_iter : int
        Number of quasi-steady iterations for mean concentration (default 3).

    Returns
    -------
    SeparatorState
        Cathode-side concentration, potential, and individual loss terms.
    """
    if not params.enabled:
        return SeparatorState(
            c_cathode=params.c_ref,
            phi_e_cathode=0.0,
            dphi_ohm=0.0,
            dphi_conc=0.0,
            eta_li=0.0,
        )

    I_dis = abs(I_app)  # positive discharge current density
    c_sep_mean = params.c_ref

    # Quasi-steady iteration for mean separator concentration
    for _ in range(n_iter):
        D_eff = electrolyte_diffusion_coefficient(c_sep_mean, T) * params.porosity ** params.bruggeman
        k_eff = electrolyte_ionic_conductivity(c_sep_mean, T) * params.porosity ** params.bruggeman

        dc_sep = I_dis * (1.0 - params.t_plus) * params.thickness / (F * D_eff)
        c_sep_cathode = max(params.c_floor, params.c_ref - dc_sep)
        c_sep_mean = 0.5 * (params.c_ref + c_sep_cathode)

    # Ohmic potential drop across separator
    dphi_ohm = I_dis * params.thickness / k_eff

    # Concentration overpotential (Nernst term)
    c_cathode_safe = max(c_sep_cathode, params.c_floor)
    if params.include_concentration_overpotential:
        dphi_conc = (R * T / F) * (1.0 - params.t_plus) * math.log(params.c_ref / c_cathode_safe)
    else:
        dphi_conc = 0.0

    # Li foil BV overpotential (disabled = ideal Li reference)
    if params.include_li_foil_bv and I_dis > 0.0:
        eta_li = (2.0 * R * T / F) * math.asinh(I_dis / (2.0 * params.i0_foil))
    else:
        eta_li = 0.0

    # Cathode-side electrolyte reference potential (negative for discharge)
    phi_e_cathode = -(dphi_ohm + dphi_conc + eta_li)

    return SeparatorState(
        c_cathode=c_cathode_safe,
        phi_e_cathode=phi_e_cathode,
        dphi_ohm=dphi_ohm,
        dphi_conc=dphi_conc,
        eta_li=eta_li,
    )
