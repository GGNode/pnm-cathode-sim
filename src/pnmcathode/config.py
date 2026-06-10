"""Public configuration dataclasses for cathode simulations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


TransportFunction = Callable[[float, float], float]
OcvFunction = Callable[[float], float]


@dataclass(frozen=True)
class CathodeGeometry:
    """Synthetic cubic cathode network geometry and phase fractions."""

    shape: tuple[int, int, int] = (10, 10, 10)
    spacing: float = 1e-5
    porosity: float = 0.35
    cbd_fraction: float = 0.10
    seed: int | None = None
    throat_scale: float = 1.0
    geometric_area: float | None = None

    def __post_init__(self) -> None:
        if len(self.shape) != 3:
            raise ValueError("shape must contain exactly three dimensions")
        if any(int(dim) <= 1 for dim in self.shape):
            raise ValueError("all shape dimensions must be greater than 1")
        if self.spacing <= 0.0:
            raise ValueError("spacing must be positive")
        if not 0.0 < self.porosity < 1.0:
            raise ValueError("porosity must be between 0 and 1")
        if not 0.0 <= self.cbd_fraction < 1.0:
            raise ValueError("cbd_fraction must be between 0 and 1")
        if self.porosity + self.cbd_fraction >= 1.0:
            raise ValueError("porosity + cbd_fraction must be less than 1")
        if self.throat_scale <= 0.0:
            raise ValueError("throat_scale must be positive")
        if self.geometric_area is not None and self.geometric_area <= 0.0:
            raise ValueError("geometric_area must be positive when provided")


@dataclass(frozen=True)
class ActiveMaterial:
    """Active solid material model."""

    name: str
    cs_max: float
    sigma: float
    diffusivity: TransportFunction
    ocv: OcvFunction
    ocv_derivative: OcvFunction | None = None

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("name must be non-empty")
        if self.cs_max <= 0.0:
            raise ValueError("cs_max must be positive")
        if self.sigma < 0.0:
            raise ValueError("sigma must be non-negative")


@dataclass(frozen=True)
class Electrolyte:
    """Electrolyte transport model."""

    name: str
    diffusivity: TransportFunction
    conductivity: TransportFunction
    c_init: float = 1200.0
    transference_number: float = 0.363
    bruggeman: float = 1.5

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("name must be non-empty")
        if self.c_init <= 0.0:
            raise ValueError("c_init must be positive")
        if not 0.0 <= self.transference_number <= 1.0:
            raise ValueError("transference_number must be between 0 and 1")
        if self.bruggeman < 0.0:
            raise ValueError("bruggeman must be non-negative")


@dataclass(frozen=True)
class ConductiveAdditive:
    """Non-active electronically conductive phase."""

    name: str = "CBD"
    sigma: float = 760.0

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("name must be non-empty")
        if self.sigma < 0.0:
            raise ValueError("sigma must be non-negative")


@dataclass(frozen=True)
class Kinetics:
    """Butler-Volmer kinetic parameters."""

    k0: float = 5e-10
    alpha_a: float = 0.5
    alpha_c: float = 0.5

    def __post_init__(self) -> None:
        if self.k0 <= 0.0:
            raise ValueError("k0 must be positive")
        if self.alpha_a <= 0.0:
            raise ValueError("alpha_a must be positive")
        if self.alpha_c <= 0.0:
            raise ValueError("alpha_c must be positive")


@dataclass(frozen=True)
class Separator:
    """Collapsed one-dimensional separator boundary model parameters."""

    enabled: bool = False
    thickness: float = 25e-6
    porosity: float = 0.39
    bruggeman: float = 1.5
    t_plus: float = 0.363
    c_ref: float = 1200.0
    c_floor: float = 1.0
    include_concentration_overpotential: bool = True
    include_li_foil_bv: bool = False
    i0_foil: float = 19.0
    alpha_foil: float = 0.5

    def __post_init__(self) -> None:
        if self.thickness <= 0.0:
            raise ValueError("thickness must be positive")
        if not 0.0 < self.porosity <= 1.0:
            raise ValueError("porosity must be between 0 and 1")
        if self.bruggeman < 0.0:
            raise ValueError("bruggeman must be non-negative")
        if not 0.0 <= self.t_plus <= 1.0:
            raise ValueError("t_plus must be between 0 and 1")
        if self.c_ref <= 0.0:
            raise ValueError("c_ref must be positive")
        if self.c_floor <= 0.0:
            raise ValueError("c_floor must be positive")
        if self.i0_foil <= 0.0:
            raise ValueError("i0_foil must be positive")
        if self.alpha_foil <= 0.0:
            raise ValueError("alpha_foil must be positive")


@dataclass(frozen=True)
class DischargeProtocol:
    """Galvanostatic discharge protocol."""

    c_rate: float
    cutoff_voltage: float = 2.5
    initial_soc: float = 0.5
    max_steps: int = 1000
    save_spatial: bool = False
    spatial_interval: int = 1

    def __post_init__(self) -> None:
        if self.c_rate <= 0.0:
            raise ValueError("c_rate must be positive")
        if self.cutoff_voltage <= 0.0:
            raise ValueError("cutoff_voltage must be positive")
        if not 0.0 <= self.initial_soc <= 1.0:
            raise ValueError("initial_soc must be between 0 and 1")
        if self.max_steps <= 0:
            raise ValueError("max_steps must be positive")
        if self.spatial_interval <= 0:
            raise ValueError("spatial_interval must be positive")


@dataclass(frozen=True)
class SolverSettings:
    """Numerical solver controls."""

    temperature: float = 298.15
    dt: float | None = None
    dt_min: float = 1e-6
    dt_max: float = 300.0
    voltage_jump_limit: float = 0.05
    newton_tol: float = 1e-8
    newton_max_iter: int = 200

    def __post_init__(self) -> None:
        if self.temperature <= 0.0:
            raise ValueError("temperature must be positive")
        if self.dt is not None and self.dt <= 0.0:
            raise ValueError("dt must be positive when provided")
        if self.dt_min <= 0.0:
            raise ValueError("dt_min must be positive")
        if self.dt_max <= 0.0:
            raise ValueError("dt_max must be positive")
        if self.dt_min > self.dt_max:
            raise ValueError("dt_min must be less than or equal to dt_max")
        if self.voltage_jump_limit <= 0.0:
            raise ValueError("voltage_jump_limit must be positive")
        if self.newton_tol <= 0.0:
            raise ValueError("newton_tol must be positive")
        if self.newton_max_iter <= 0:
            raise ValueError("newton_max_iter must be positive")
