"""Named material presets used by published validation cases."""

from __future__ import annotations

from pnmcathode.config import ActiveMaterial, ConductiveAdditive, Electrolyte, Separator
from pnmcathode.physics.electrolyte import (
    electrolyte_diffusion_coefficient,
    electrolyte_ionic_conductivity,
)
from pnmcathode.physics.ocv import nmc532_ocv, ocv_derivative
from pnmcathode.physics.solid import nmc532_diffusion_coefficient


def nmc532_khan2021() -> ActiveMaterial:
    """Return the NMC532 active material model used for Khan et al. 2021."""

    return ActiveMaterial(
        name="NMC532 (Khan 2021)",
        cs_max=48900.0,
        sigma=0.01,
        diffusivity=nmc532_diffusion_coefficient,
        ocv=nmc532_ocv,
        ocv_derivative=ocv_derivative,
    )


def cbd_khan2021() -> ConductiveAdditive:
    """Return the conductive carbon/binder domain used for Khan et al. 2021."""

    return ConductiveAdditive(name="CBD (Khan 2021)", sigma=760.0)


def electrolyte_khan2021() -> Electrolyte:
    """Return the electrolyte transport model used for Khan et al. 2021."""

    return Electrolyte(
        name="LiPF6 carbonate electrolyte (Khan 2021)",
        c_init=1200.0,
        diffusivity=electrolyte_diffusion_coefficient,
        conductivity=electrolyte_ionic_conductivity,
        transference_number=0.363,
        bruggeman=1.5,
    )


def separator_khan2021(enabled: bool = True) -> Separator:
    """Return the separator boundary model used for Khan et al. 2021."""

    return Separator(
        enabled=enabled,
        thickness=25e-6,
        porosity=0.39,
        bruggeman=1.5,
        t_plus=0.363,
        c_ref=1200.0,
        c_floor=1.0,
        include_concentration_overpotential=True,
        include_li_foil_bv=False,
        i0_foil=19.0,
        alpha_foil=0.5,
    )


__all__ = [
    "nmc532_khan2021",
    "cbd_khan2021",
    "electrolyte_khan2021",
    "separator_khan2021",
]
