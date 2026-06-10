# Toolkit Execution Plan

This is the step-by-step implementation contract for converting the current Khan 2021 reproduction code into a pip-installable `pnmcathode` toolkit while preserving legacy `src.*` imports used by the existing tests.

Rules for Claude Code:

- Execute steps in order.
- Use `git mv` for moved files so history is preserved.
- Do not rewrite numerical kernels unless a step explicitly says to edit them.
- After moving kernel modules under `src/pnmcathode/`, replace their internal imports from `src.*` to `pnmcathode.*`.
- Keep compatibility wrappers under the legacy `src.*` package until all old imports are removed.

## Step 1: Package Shell

Goal: create the `src/pnmcathode/` package, move existing kernel modules into it, create legacy compatibility wrappers, create `pyproject.toml`, verify editable install, and verify current tests still pass.

### 1.1 Move Existing Kernel Packages

Run these commands from the repository root:

```bash
mkdir -p src/pnmcathode
git mv src/network src/pnmcathode/network
git mv src/physics src/pnmcathode/physics
git mv src/solver src/pnmcathode/solver
git mv src/post src/pnmcathode/post
```

Modify every moved Python file under `src/pnmcathode/` by replacing imports as follows:

```text
from src. -> from pnmcathode.
import src. -> import pnmcathode.
```

The expected affected files are:

- `src/pnmcathode/physics/separator.py`
- `src/pnmcathode/solver/transient.py`
- `src/pnmcathode/solver/steady.py`
- `src/pnmcathode/solver/single_pore.py`
- `src/pnmcathode/post/analysis.py`

### 1.2 Create Legacy Compatibility Wrapper Files

Create `src/__init__.py` with exactly:

```python
"""Compatibility package for legacy ``src.*`` imports."""

from pnmcathode import __version__

__all__ = ["__version__"]
```

Create `src/network/__init__.py` with exactly:

```python
"""Compatibility wrappers for legacy ``src.network`` imports."""

from pnmcathode.network.generator import check_percolation, create_cathode_network

__all__ = ["check_percolation", "create_cathode_network"]
```

Create `src/network/generator.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.network.generator``."""

from pnmcathode.network.generator import *  # noqa: F401,F403
```

Create `src/physics/__init__.py` with exactly:

```python
"""Compatibility wrappers for legacy ``src.physics`` imports."""
```

Create `src/physics/electrolyte.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.physics.electrolyte``."""

from pnmcathode.physics.electrolyte import *  # noqa: F401,F403
```

Create `src/physics/ocv.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.physics.ocv``."""

from pnmcathode.physics.ocv import *  # noqa: F401,F403
```

Create `src/physics/reaction.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.physics.reaction``."""

from pnmcathode.physics.reaction import *  # noqa: F401,F403
```

Create `src/physics/separator.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.physics.separator``."""

from pnmcathode.physics.separator import *  # noqa: F401,F403
```

Create `src/physics/solid.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.physics.solid``."""

from pnmcathode.physics.solid import *  # noqa: F401,F403
```

Create `src/solver/__init__.py` with exactly:

```python
"""Compatibility wrappers for legacy ``src.solver`` imports."""

from pnmcathode.solver.single_pore import SinglePoreDischarge
from pnmcathode.solver.steady import SteadyStateSolver
from pnmcathode.solver.transient import TransientSolver

__all__ = ["SinglePoreDischarge", "SteadyStateSolver", "TransientSolver"]
```

Create `src/solver/single_pore.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.solver.single_pore``."""

from pnmcathode.solver.single_pore import *  # noqa: F401,F403
```

Create `src/solver/steady.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.solver.steady``."""

from pnmcathode.solver.steady import *  # noqa: F401,F403
```

Create `src/solver/transient.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.solver.transient``."""

from pnmcathode.solver.transient import *  # noqa: F401,F403
```

Create `src/post/__init__.py` with exactly:

```python
"""Compatibility wrappers for legacy ``src.post`` imports."""
```

Create `src/post/analysis.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.post.analysis``."""

from pnmcathode.post.analysis import *  # noqa: F401,F403
```

Create `src/post/visualization.py` with exactly:

```python
"""Compatibility wrapper for ``pnmcathode.post.visualization``."""

from pnmcathode.post.visualization import *  # noqa: F401,F403
```

### 1.3 Create `pyproject.toml`

Create `pyproject.toml` now using the exact final content shown in Step 9. Step 9 is a final metadata check, not a second competing version of the file.

### 1.4 Verify Step 1

Run:

```bash
python -m pip install -e .
python -m pytest
```

Expected result:

- `import pnmcathode` works.
- Existing tests that import `src.*` still pass.

## Step 2: Config Dataclasses

Goal: create all public configuration dataclasses in `src/pnmcathode/config.py`.

Create `src/pnmcathode/config.py` with exactly:

```python
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
```

Add tests for this module in Step 8.

## Step 3: Presets

Goal: create Khan 2021 material presets.

Create directory `src/pnmcathode/materials/`.

Create `src/pnmcathode/materials/__init__.py` with exactly:

```python
"""Material configuration objects and presets."""

from pnmcathode.config import ActiveMaterial, ConductiveAdditive, Electrolyte, Separator

__all__ = ["ActiveMaterial", "ConductiveAdditive", "Electrolyte", "Separator"]
```

Create `src/pnmcathode/materials/presets.py` with exactly:

```python
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
```

Add tests for these presets in Step 8.

## Step 4: DischargeResult

Goal: create a typed result wrapper with NPZ serialization.

Create `src/pnmcathode/results.py` with exactly:

```python
"""Simulation result containers and serialization helpers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


def _metadata_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


@dataclass
class DischargeResult:
    """Typed wrapper around galvanostatic discharge output arrays."""

    time: np.ndarray
    voltage: np.ndarray
    capacity_Ah_m2: np.ndarray
    current_density: float
    c_rate: float
    reached_cutoff: bool
    spatial: dict[str, np.ndarray] | None = None
    metadata: dict[str, Any] | None = None

    @property
    def final_capacity(self) -> float:
        """Final areal capacity in A h m^-2."""

        if self.capacity_Ah_m2.size == 0:
            return 0.0
        return float(self.capacity_Ah_m2[-1])

    @classmethod
    def from_solver_output(
        cls,
        raw: dict[str, Any],
        spatial: dict[str, np.ndarray] | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> "DischargeResult":
        """Build a result object from ``TransientSolver.run_discharge`` output."""

        merged_metadata: dict[str, Any] = {
            "cutoff_voltage": raw.get("cutoff_voltage"),
        }
        if metadata:
            merged_metadata.update(metadata)
        return cls(
            time=np.asarray(raw["time"], dtype=float),
            voltage=np.asarray(raw["voltage"], dtype=float),
            capacity_Ah_m2=np.asarray(raw["capacity_Ah_m2"], dtype=float),
            current_density=float(raw["I_app"]),
            c_rate=float(raw["C_rate"]),
            reached_cutoff=bool(raw["reached_cutoff"]),
            spatial=spatial,
            metadata=merged_metadata,
        )

    def to_npz(self, path: str | Path) -> None:
        """Serialize the result to a compressed NPZ file."""

        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)

        spatial = self.spatial or {}
        spatial_keys = sorted(spatial)
        metadata_json = json.dumps(
            self.metadata or {},
            default=_metadata_default,
            sort_keys=True,
        )
        payload: dict[str, Any] = {
            "time": np.asarray(self.time, dtype=float),
            "voltage": np.asarray(self.voltage, dtype=float),
            "capacity_Ah_m2": np.asarray(self.capacity_Ah_m2, dtype=float),
            "current_density": np.asarray(self.current_density, dtype=float),
            "c_rate": np.asarray(self.c_rate, dtype=float),
            "reached_cutoff": np.asarray(self.reached_cutoff, dtype=bool),
            "spatial_keys": np.asarray(spatial_keys, dtype=str),
            "metadata_json": np.asarray(metadata_json, dtype=str),
        }
        for key in spatial_keys:
            payload[f"spatial.{key}"] = np.asarray(spatial[key])

        np.savez_compressed(target, **payload)

    @classmethod
    def from_npz(cls, path: str | Path) -> "DischargeResult":
        """Deserialize a result from an NPZ file written by :meth:`to_npz`."""

        with np.load(Path(path), allow_pickle=False) as data:
            spatial_keys = [str(key) for key in data["spatial_keys"].tolist()]
            spatial = {
                key: np.asarray(data[f"spatial.{key}"])
                for key in spatial_keys
            }
            metadata_json = str(data["metadata_json"].item())
            metadata = json.loads(metadata_json) if metadata_json else {}
            return cls(
                time=np.asarray(data["time"], dtype=float),
                voltage=np.asarray(data["voltage"], dtype=float),
                capacity_Ah_m2=np.asarray(data["capacity_Ah_m2"], dtype=float),
                current_density=float(data["current_density"].item()),
                c_rate=float(data["c_rate"].item()),
                reached_cutoff=bool(data["reached_cutoff"].item()),
                spatial=spatial or None,
                metadata=metadata,
            )
```

Add serialization tests in Step 8.

## Step 5: Simulation Wrapper

Goal: create the orchestration layer that builds a solver and returns `DischargeResult`.

Create `src/pnmcathode/simulation.py` with exactly:

```python
"""High-level simulation orchestration."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

from pnmcathode.cathode import Cathode
from pnmcathode.config import DischargeProtocol, Kinetics, Separator, SolverSettings
from pnmcathode.physics.separator import SeparatorParams
from pnmcathode.results import DischargeResult
from pnmcathode.solver.transient import TransientSolver


def _separator_params(separator: Separator) -> SeparatorParams:
    return SeparatorParams(
        enabled=separator.enabled,
        thickness=separator.thickness,
        porosity=separator.porosity,
        bruggeman=separator.bruggeman,
        t_plus=separator.t_plus,
        c_ref=separator.c_ref,
        c_floor=separator.c_floor,
        include_concentration_overpotential=separator.include_concentration_overpotential,
        include_li_foil_bv=separator.include_li_foil_bv,
        i0_foil=separator.i0_foil,
        alpha_foil=separator.alpha_foil,
    )


class Simulation:
    """High-level galvanostatic discharge simulation."""

    def __init__(
        self,
        cathode: Cathode,
        protocol: DischargeProtocol,
        kinetics: Kinetics | None = None,
        separator: Separator | None = None,
        settings: SolverSettings | None = None,
    ):
        self.cathode = cathode
        self.protocol = protocol
        self.kinetics = kinetics or Kinetics()
        self.separator = separator or Separator(enabled=False)
        self.settings = settings or SolverSettings()

    def _metadata(self) -> dict[str, Any]:
        return {
            "geometry": asdict(self.cathode.geometry),
            "active_material": self.cathode.active_material.name,
            "electrolyte": self.cathode.electrolyte.name,
            "conductive_additive": self.cathode.conductive_additive.name,
            "protocol": asdict(self.protocol),
            "kinetics": asdict(self.kinetics),
            "separator": asdict(self.separator),
            "settings": asdict(self.settings),
        }

    def run(self) -> DischargeResult:
        """Run the discharge simulation and return typed results."""

        net = self.cathode.ensure_network()
        solver = TransientSolver(
            net,
            T=self.settings.temperature,
            k0=self.kinetics.k0,
            geometric_area=self.cathode.geometry.geometric_area,
            separator=_separator_params(self.separator),
        )
        solver.set_concentration(
            c_e=self.cathode.electrolyte.c_init,
            c_s=self.protocol.initial_soc * self.cathode.active_material.cs_max,
        )
        raw = solver.run_discharge(
            C_rate=self.protocol.c_rate,
            cutoff_voltage=self.protocol.cutoff_voltage,
            dt=self.settings.dt,
            dt_min=self.settings.dt_min,
            dt_max=self.settings.dt_max,
            max_steps=self.protocol.max_steps,
            voltage_jump_limit=self.settings.voltage_jump_limit,
        )
        return DischargeResult.from_solver_output(raw, metadata=self._metadata())
```

Add tests for this wrapper in Step 8. Use monkeypatching in the unit test so the public API test is fast and does not depend on a full transient solve.

## Step 6: Cathode Class

Goal: create the public `Cathode` object that bundles geometry, materials, and the generated network.

Create `src/pnmcathode/cathode.py` with exactly:

```python
"""Cathode object and network factory methods."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pnmcathode.config import (
    ActiveMaterial,
    CathodeGeometry,
    ConductiveAdditive,
    Electrolyte,
)
from pnmcathode.network.generator import create_cathode_network


@dataclass
class Cathode:
    """Bundle cathode geometry, material models, and an OpenPNM network."""

    geometry: CathodeGeometry
    active_material: ActiveMaterial
    electrolyte: Electrolyte
    conductive_additive: ConductiveAdditive = field(default_factory=ConductiveAdditive)
    network: Any | None = None

    @classmethod
    def cubic(
        cls,
        geometry: CathodeGeometry | None = None,
        active_material: ActiveMaterial | None = None,
        electrolyte: Electrolyte | None = None,
        conductive_additive: ConductiveAdditive | None = None,
    ) -> "Cathode":
        """Create a cubic synthetic cathode network."""

        from pnmcathode.materials.presets import (
            cbd_khan2021,
            electrolyte_khan2021,
            nmc532_khan2021,
        )

        cathode = cls(
            geometry=geometry or CathodeGeometry(),
            active_material=active_material or nmc532_khan2021(),
            electrolyte=electrolyte or electrolyte_khan2021(),
            conductive_additive=conductive_additive or cbd_khan2021(),
        )
        cathode.build_network()
        return cathode

    def build_network(self) -> Any:
        """Generate and attach a network using the configured geometry."""

        net = create_cathode_network(
            shape=list(self.geometry.shape),
            spacing=self.geometry.spacing,
            porosity=self.geometry.porosity,
            cbd_fraction=self.geometry.cbd_fraction,
            seed=self.geometry.seed,
            throat_scale=self.geometry.throat_scale,
        )
        self._apply_material_properties(net)
        self.network = net
        return net

    def ensure_network(self) -> Any:
        """Return the attached network, generating it if needed."""

        if self.network is None:
            return self.build_network()
        return self.network

    def _apply_material_properties(self, net: Any) -> None:
        nmc_mask = np.asarray(net["pore.nmc"], dtype=bool)
        cbd_mask = np.asarray(net["pore.cbd"], dtype=bool)
        net["pore.cs_max"] = np.where(nmc_mask, self.active_material.cs_max, 0.0)
        net["pore.sigma_nmc"] = np.where(nmc_mask, self.active_material.sigma, 0.0)
        net["pore.sigma_cbd"] = np.where(cbd_mask, self.conductive_additive.sigma, 0.0)
```

Add tests for this class in Step 8.

## Step 7: Public API and Examples

Goal: define clean top-level exports and move paper-specific scripts into examples.

### 7.1 Create Public API Files

Create `src/pnmcathode/__init__.py` with exactly:

```python
"""Reusable pore-network cathode simulation toolkit."""

from pnmcathode.cathode import Cathode
from pnmcathode.config import (
    ActiveMaterial,
    CathodeGeometry,
    ConductiveAdditive,
    DischargeProtocol,
    Electrolyte,
    Kinetics,
    Separator,
    SolverSettings,
)
from pnmcathode.results import DischargeResult
from pnmcathode.simulation import Simulation

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "ActiveMaterial",
    "Cathode",
    "CathodeGeometry",
    "ConductiveAdditive",
    "DischargeProtocol",
    "DischargeResult",
    "Electrolyte",
    "Kinetics",
    "Separator",
    "Simulation",
    "SolverSettings",
]
```

Create `src/pnmcathode/analysis.py` with exactly:

```python
"""Public analysis helpers."""

from pnmcathode.post.analysis import *  # noqa: F401,F403
```

Create `src/pnmcathode/plotting.py` with exactly:

```python
"""Public plotting helpers."""

from pnmcathode.post.visualization import *  # noqa: F401,F403
```

### 7.2 Move Paper-Specific Scripts

Run:

```bash
mkdir -p examples/reproduce_khan2021
git mv scripts/run_discharge.py examples/reproduce_khan2021/run_discharge.py
git mv scripts/validate_paper_figures.py examples/reproduce_khan2021/validate_paper_figures.py
git mv scripts/parallel_validation.py examples/reproduce_khan2021/parallel_validation.py
```

In each moved example file:

- Replace `ROOT = Path(__file__).resolve().parents[1]` with `ROOT = Path(__file__).resolve().parents[2]`.
- Replace `from src.network.generator import ...` with `from pnmcathode.network.generator import ...`.
- Replace `from src.physics...` with `from pnmcathode.physics...`.
- Replace `from src.solver...` with `from pnmcathode.solver...`.

Expected imports after replacement:

`examples/reproduce_khan2021/run_discharge.py`:

```python
from pnmcathode.network.generator import create_cathode_network
from pnmcathode.solver.transient import TransientSolver
```

`examples/reproduce_khan2021/validate_paper_figures.py`:

```python
from pnmcathode.network.generator import check_percolation, create_cathode_network
from pnmcathode.physics.ocv import nmc532_ocv
from pnmcathode.physics.separator import SeparatorParams
from pnmcathode.solver.transient import TransientSolver
```

`examples/reproduce_khan2021/parallel_validation.py` has local imports inside worker functions. Update those local imports to:

```python
from pnmcathode.network.generator import create_cathode_network
from pnmcathode.physics.separator import SeparatorParams
from pnmcathode.solver.transient import TransientSolver
```

and:

```python
from pnmcathode.physics.separator import separator_boundary, SeparatorParams as SP
```

Do not delete validation data under `data/`.

## Step 8: Tests

Goal: add public API tests for the new toolkit layer.

Create `tests/test_toolkit_api.py` with exactly:

```python
from __future__ import annotations

import numpy as np

import pnmcathode
import pnmcathode.analysis as analysis
import pnmcathode.plotting as plotting
from pnmcathode import (
    ActiveMaterial,
    Cathode,
    CathodeGeometry,
    ConductiveAdditive,
    DischargeProtocol,
    DischargeResult,
    Electrolyte,
    Kinetics,
    Separator,
    Simulation,
    SolverSettings,
)
from pnmcathode.materials import ActiveMaterial as MaterialActiveMaterial
from pnmcathode.materials.presets import (
    cbd_khan2021,
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)


def constant_transport(_concentration: float, _temperature: float) -> float:
    return 1.0


def constant_ocv(_soc: float) -> float:
    return 4.0


def test_config_dataclasses_and_public_exports():
    active = ActiveMaterial(
        name="test active",
        cs_max=100.0,
        sigma=2.0,
        diffusivity=constant_transport,
        ocv=constant_ocv,
    )
    electrolyte = Electrolyte(
        name="test electrolyte",
        diffusivity=constant_transport,
        conductivity=constant_transport,
    )
    additive = ConductiveAdditive(name="test additive", sigma=5.0)
    kinetics = Kinetics(k0=1e-9)
    separator = Separator(enabled=True)

    assert pnmcathode.__version__ == "0.1.0"
    assert MaterialActiveMaterial is ActiveMaterial
    assert active.ocv(0.5) == 4.0
    assert electrolyte.c_init == 1200.0
    assert additive.sigma == 5.0
    assert kinetics.alpha_a == 0.5
    assert separator.enabled is True
    assert hasattr(analysis, "pore_size_distribution")
    assert hasattr(plotting, "plot_discharge_curve")


def test_cathode_cubic_creation():
    geometry = CathodeGeometry(shape=(3, 3, 3), porosity=0.45, seed=7)

    cathode = Cathode.cubic(geometry=geometry)

    assert cathode.network is not None
    assert cathode.network.Np == 27
    assert cathode.active_material.name == "NMC532 (Khan 2021)"
    assert cathode.electrolyte.name == "LiPF6 carbonate electrolyte (Khan 2021)"
    assert np.any(cathode.network["pore.electrolyte"])
    assert np.any(cathode.network["pore.nmc"])


def test_simulation_run_returns_discharge_result(monkeypatch):
    calls = {}

    class DummyTransientSolver:
        def __init__(self, net, T, k0, geometric_area=None, separator=None):
            calls["net"] = net
            calls["T"] = T
            calls["k0"] = k0
            calls["geometric_area"] = geometric_area
            calls["separator"] = separator

        def set_concentration(self, c_e, c_s):
            calls["c_e"] = c_e
            calls["c_s"] = c_s

        def run_discharge(
            self,
            C_rate,
            cutoff_voltage,
            dt,
            dt_min,
            dt_max,
            max_steps,
            voltage_jump_limit,
        ):
            calls["run"] = {
                "C_rate": C_rate,
                "cutoff_voltage": cutoff_voltage,
                "dt": dt,
                "dt_min": dt_min,
                "dt_max": dt_max,
                "max_steps": max_steps,
                "voltage_jump_limit": voltage_jump_limit,
            }
            return {
                "time": np.array([0.0, 1.0]),
                "capacity_Ah_m2": np.array([0.0, 0.01]),
                "voltage": np.array([4.0, 3.9]),
                "I_app": -36.0,
                "C_rate": C_rate,
                "cutoff_voltage": cutoff_voltage,
                "reached_cutoff": False,
            }

    monkeypatch.setattr("pnmcathode.simulation.TransientSolver", DummyTransientSolver)

    cathode = Cathode.cubic(
        geometry=CathodeGeometry(shape=(3, 3, 3), seed=1, geometric_area=1.2e-6)
    )
    protocol = DischargeProtocol(c_rate=1.0, cutoff_voltage=2.5, max_steps=3)
    settings = SolverSettings(temperature=303.0, dt=0.1)

    result = Simulation(cathode, protocol, settings=settings).run()

    assert isinstance(result, DischargeResult)
    assert result.final_capacity == 0.01
    assert result.current_density == -36.0
    assert result.c_rate == 1.0
    assert calls["T"] == 303.0
    assert calls["k0"] == 5e-10
    assert calls["geometric_area"] == 1.2e-6
    assert calls["c_e"] == 1200.0
    assert calls["c_s"] == 24450.0
    assert calls["run"]["max_steps"] == 3


def test_presets_load_correctly():
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    additive = cbd_khan2021()
    separator = separator_khan2021()

    assert active.cs_max == 48900.0
    assert active.sigma == 0.01
    assert callable(active.diffusivity)
    assert callable(active.ocv)
    assert electrolyte.c_init == 1200.0
    assert callable(electrolyte.diffusivity)
    assert callable(electrolyte.conductivity)
    assert additive.sigma == 760.0
    assert separator.enabled is True
    assert separator.thickness == 25e-6


def test_discharge_result_serialization_round_trip(tmp_path):
    path = tmp_path / "result.npz"
    result = DischargeResult(
        time=np.array([0.0, 2.0]),
        voltage=np.array([4.1, 3.8]),
        capacity_Ah_m2=np.array([0.0, 0.02]),
        current_density=-12.5,
        c_rate=0.5,
        reached_cutoff=True,
        spatial={"c_e": np.array([[1200.0, 1190.0]])},
        metadata={"case": "round_trip"},
    )

    result.to_npz(path)
    loaded = DischargeResult.from_npz(path)

    np.testing.assert_allclose(loaded.time, result.time)
    np.testing.assert_allclose(loaded.voltage, result.voltage)
    np.testing.assert_allclose(loaded.capacity_Ah_m2, result.capacity_Ah_m2)
    assert loaded.current_density == result.current_density
    assert loaded.c_rate == result.c_rate
    assert loaded.reached_cutoff is True
    assert loaded.metadata == {"case": "round_trip"}
    assert loaded.spatial is not None
    np.testing.assert_allclose(loaded.spatial["c_e"], result.spatial["c_e"])
```

Run:

```bash
python -m pytest tests/test_toolkit_api.py
```

## Step 9: Final `pyproject.toml`

Goal: ensure package metadata and package discovery are final.

The final `pyproject.toml` must be exactly:

```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "pnmcathode"
version = "0.1.0"
description = "Pore-network cathode simulation toolkit for lithium-ion battery cathodes"
readme = "README.md"
requires-python = ">=3.10"
license = { text = "MIT" }
authors = [
  { name = "PNM-LIB-Cathode contributors" }
]
keywords = [
  "battery",
  "cathode",
  "pore-network",
  "simulation",
  "electrochemistry"
]
classifiers = [
  "Development Status :: 3 - Alpha",
  "Intended Audience :: Science/Research",
  "License :: OSI Approved :: MIT License",
  "Programming Language :: Python :: 3",
  "Programming Language :: Python :: 3.10",
  "Programming Language :: Python :: 3.11",
  "Programming Language :: Python :: 3.12",
  "Topic :: Scientific/Engineering",
]
dependencies = [
  "numpy>=1.24",
  "scipy>=1.10",
  "matplotlib>=3.7",
  "openpnm>=4.0",
  "porespy>=2.0",
]

[project.optional-dependencies]
yaml = ["pyyaml>=6.0"]
test = ["pytest>=7.0"]
dev = ["pytest>=7.0", "pyyaml>=6.0"]

[project.urls]
Homepage = "https://github.com/pnm-lib/pnm-lib-cathode"
Repository = "https://github.com/pnm-lib/pnm-lib-cathode"

[tool.setuptools]
packages = [
  "pnmcathode",
  "pnmcathode.materials",
  "pnmcathode.network",
  "pnmcathode.physics",
  "pnmcathode.solver",
  "pnmcathode.post",
  "src",
  "src.network",
  "src.physics",
  "src.solver",
  "src.post",
]

[tool.setuptools.package-dir]
pnmcathode = "src/pnmcathode"
"pnmcathode.materials" = "src/pnmcathode/materials"
"pnmcathode.network" = "src/pnmcathode/network"
"pnmcathode.physics" = "src/pnmcathode/physics"
"pnmcathode.solver" = "src/pnmcathode/solver"
"pnmcathode.post" = "src/pnmcathode/post"
src = "src"
"src.network" = "src/network"
"src.physics" = "src/physics"
"src.solver" = "src/solver"
"src.post" = "src/post"

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["."]
```

## Step 10: Verify

Goal: run the full verification loop and fix import issues.

Run:

```bash
python -m pip install -e .[dev]
python - <<'PY'
import pnmcathode
from pnmcathode import Cathode, DischargeProtocol, Simulation
from pnmcathode.materials.presets import nmc532_khan2021, electrolyte_khan2021
from src.network.generator import create_cathode_network

print("pnmcathode", pnmcathode.__version__)
print(Cathode)
print(DischargeProtocol)
print(Simulation)
print(nmc532_khan2021().name)
print(electrolyte_khan2021().name)
print(create_cathode_network)
PY
python -m pytest
```

If imports fail:

1. Check that moved files under `src/pnmcathode/` no longer import from `src.*`.
2. Check that every legacy wrapper file in Step 1.2 exists.
3. Check that `src/pnmcathode/__init__.py` exports `__version__`.
4. Check that `pyproject.toml` package mappings match the actual file layout.

Completion criteria:

- `python -m pip install -e .[dev]` succeeds.
- `import pnmcathode` succeeds.
- `from src.network.generator import create_cathode_network` still succeeds.
- `python -m pytest` passes.
