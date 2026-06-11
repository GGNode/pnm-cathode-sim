# Toolkit Design Plan

This document describes how to turn the current Khan et al. 2021 reproduction code into a reusable, pip-installable Python toolkit for pore-network cathode simulation.

The core strategy is:

1. Keep the existing numerical kernels initially intact.
2. Add a stable public package API around configuration, materials, geometry, simulation, and results.
3. Move Khan-specific values into named presets and examples.
4. Gradually replace hardcoded material functions inside the solvers with injected model objects.

## A. Package Design

### Package Name Options

Recommended name: `pnmcathode`

Three viable options:

1. `pnmcathode`
   - Short, import-friendly, specific to pore-network cathode modeling.
   - Example: `from pnmcathode import Cathode, Simulation`

2. `batterypnm`
   - Broader, leaves room for anodes, full cells, and alternate porous-electrode models.
   - Example: `from batterypnm import Cathode, GalvanostaticDischarge`

3. `pnm_battery`
   - Explicit and readable, but slightly less idiomatic because of the underscore.
   - Example: `from pnm_battery import materials, run_discharge`

Use `pnmcathode` for the first package release because the current codebase is cathode-specific.

### Proposed Package Layout

```text
pnm-lib-cathode/
|-- pyproject.toml
|-- src/
|   `-- pnmcathode/
|       |-- __init__.py
|       |-- config.py
|       |-- network/
|       |   |-- __init__.py
|       |   `-- generator.py
|       |-- materials/
|       |   |-- __init__.py
|       |   |-- active.py
|       |   |-- electrolyte.py
|       |   |-- separator.py
|       |   `-- presets.py
|       |-- physics/
|       |   |-- __init__.py
|       |   |-- reaction.py
|       |   |-- ocv.py
|       |   |-- electrolyte.py
|       |   |-- solid.py
|       |   `-- separator.py
|       |-- solver/
|       |   |-- __init__.py
|       |   |-- steady.py
|       |   `-- transient.py
|       |-- simulation.py
|       |-- results.py
|       |-- analysis.py
|       |-- plotting.py
|       `-- io.py
|-- examples/
|   |-- simple_1c.py
|   |-- sweep_c_rates.py
|   |-- custom_models.py
|   `-- reproduce_khan2021/
|       |-- run_discharge.py
|       |-- validate_paper_figures.py
|       `-- parallel_validation.py
`-- tests/
```

During migration, keep the current `src/network`, `src/physics`, `src/solver`, and `src/post` modules as compatibility wrappers or move them under `src/pnmcathode` and leave wrapper imports behind.

### Top-Level Import API

The public API should expose common user workflows without requiring users to know the internal solver modules.

```python
from pnmcathode import (
    Cathode,
    CathodeGeometry,
    Simulation,
    DischargeProtocol,
    SolverSettings,
)
from pnmcathode.materials import ActiveMaterial, Electrolyte, Separator
from pnmcathode.materials.presets import nmc532_khan2021, electrolyte_khan2021
from pnmcathode.analysis import capacity, voltage_capacity_curve
from pnmcathode.plotting import plot_voltage_capacity
```

Lower-level imports remain available for advanced users:

```python
from pnmcathode.network import create_cathode_network, check_percolation
from pnmcathode.solver import TransientSolver, SteadyStateSolver
from pnmcathode.physics.reaction import butler_volmer, exchange_current_density
```

### Core User-Facing Objects

#### `CathodeGeometry`

Defines the synthetic network geometry and phase fractions.

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class CathodeGeometry:
    shape: tuple[int, int, int] = (10, 10, 10)
    spacing: float = 1e-5
    porosity: float = 0.35
    cbd_fraction: float = 0.10
    seed: int | None = None
    throat_scale: float = 1.0
    geometric_area: float | None = None
```

#### `ActiveMaterial`

Represents active-solid parameters and callbacks.

```python
from dataclasses import dataclass
from typing import Callable

@dataclass(frozen=True)
class ActiveMaterial:
    name: str
    cs_max: float
    sigma: float
    diffusivity: Callable[[float, float], float]
    ocv: Callable[[float], float]
    ocv_derivative: Callable[[float], float] | None = None
```

#### `Electrolyte`

Represents liquid transport and initial concentration.

```python
@dataclass(frozen=True)
class Electrolyte:
    name: str
    c_init: float = 1200.0
    diffusivity: Callable[[float, float], float]
    conductivity: Callable[[float, float], float]
    transference_number: float = 0.363
    bruggeman: float = 1.5
```

#### `ConductiveAdditive`

Represents the non-active solid phase.

```python
@dataclass(frozen=True)
class ConductiveAdditive:
    name: str = "CBD"
    sigma: float = 760.0
```

#### `Kinetics`

Represents reaction model parameters.

```python
@dataclass(frozen=True)
class Kinetics:
    k0: float = 5e-10
    alpha_a: float = 0.5
    alpha_c: float = 0.5
```

#### `Cathode`

Bundles geometry, material models, and the generated OpenPNM network.

```python
@dataclass
class Cathode:
    geometry: CathodeGeometry
    active_material: ActiveMaterial
    electrolyte: Electrolyte
    conductive_additive: ConductiveAdditive = ConductiveAdditive()
    network: object | None = None

    @classmethod
    def cubic(
        cls,
        geometry: CathodeGeometry | None = None,
        active_material: ActiveMaterial | None = None,
        electrolyte: Electrolyte | None = None,
        conductive_additive: ConductiveAdditive | None = None,
    ) -> "Cathode":
        ...
```

#### `DischargeProtocol`

Defines the experiment to run.

```python
@dataclass(frozen=True)
class DischargeProtocol:
    c_rate: float
    cutoff_voltage: float = 2.5
    initial_soc: float = 0.5
    max_steps: int = 1000
    save_spatial: bool = False
    spatial_interval: int = 1
```

#### `SolverSettings`

Collects numerical controls separately from physics.

```python
@dataclass(frozen=True)
class SolverSettings:
    temperature: float = 298.15
    dt: float | None = None
    dt_min: float = 1e-6
    dt_max: float = 300.0
    voltage_jump_limit: float = 0.05
    newton_tol: float = 1e-8
    newton_max_iter: int = 200
```

#### `Simulation`

Main orchestration layer.

```python
class Simulation:
    def __init__(
        self,
        cathode: Cathode,
        protocol: DischargeProtocol,
        kinetics: Kinetics | None = None,
        separator: Separator | None = None,
        settings: SolverSettings | None = None,
    ):
        ...

    def run(self) -> "DischargeResult":
        ...
```

#### `DischargeResult`

Typed wrapper around solver output arrays.

```python
@dataclass
class DischargeResult:
    time: np.ndarray
    voltage: np.ndarray
    capacity_Ah_m2: np.ndarray
    current_density: float
    c_rate: float
    reached_cutoff: bool
    spatial: dict[str, np.ndarray] | None = None
    metadata: dict[str, object] | None = None

    @property
    def final_capacity(self) -> float:
        return float(self.capacity_Ah_m2[-1])

    def to_npz(self, path: str | Path) -> None:
        ...

    @classmethod
    def from_npz(cls, path: str | Path) -> "DischargeResult":
        ...
```

### Simulation in Fewer Than 10 Lines

```python
from pnmcathode import Cathode, DischargeProtocol, Simulation
from pnmcathode.materials.presets import nmc532_khan2021, electrolyte_khan2021

cathode = Cathode.cubic(
    active_material=nmc532_khan2021(),
    electrolyte=electrolyte_khan2021(),
)
protocol = DischargeProtocol(c_rate=1.0, cutoff_voltage=2.5)
result = Simulation(cathode, protocol).run()
print(result.final_capacity)
```

## B. Refactoring Plan

### What Stays Mostly As-Is

Keep these components as numerical kernels in the first migration phase:

- `network/generator.py`
  - `create_cathode_network`
  - `check_percolation`
- `physics/reaction.py`
  - `butler_volmer`
  - `exchange_current_density`
  - constants `F`, `R`
- `physics/electrolyte.py`
  - Existing Khan electrolyte correlation functions, renamed or re-exported as presets.
- `physics/ocv.py`
  - Existing NMC532 OCV and derivative as preset model functions.
- `physics/separator.py`
  - `SeparatorParams` and `separator_boundary`, eventually renamed to public `Separator`.
- `solver/steady.py` and `solver/transient.py`
  - Keep the algorithms, but change constructors to accept model objects over time.
- `post/analysis.py` and `post/visualization.py`
  - Keep functions, but expose them under `pnmcathode.analysis` and `pnmcathode.plotting`.

### What Needs To Change

The current code has Khan/NMC532 assumptions embedded in low-level places. These need to become parameters or presets:

| Current Hardcoding | Location | Target Design |
|---|---|---|
| `cs_max = 48900.0` | `generator.py`, `transient.py`, `analysis.py`, reaction calls | `ActiveMaterial.cs_max` |
| `sigma_nmc = 0.01` | `generator.py`, `steady.py` | `ActiveMaterial.sigma` |
| `sigma_cbd = 760.0` | `generator.py`, `steady.py` | `ConductiveAdditive.sigma` |
| `nmc532_ocv` | `steady.py`, `transient.py` | `ActiveMaterial.ocv` |
| `nmc532_diffusion_coefficient` | `transient.py` | `ActiveMaterial.diffusivity` |
| electrolyte correlations | `steady.py`, `transient.py`, `separator.py` | `Electrolyte.diffusivity`, `Electrolyte.conductivity` |
| default `k0` values | scripts, solver constructors | `Kinetics.k0` |
| paper validation defaults | `scripts/` | `examples/reproduce_khan2021/` presets |
| dictionary result shape only | `TransientSolver.run_discharge` | `DischargeResult` with backward-compatible `raw` dict |

### Library Code vs Paper-Specific Code

Library code belongs under `src/pnmcathode/`:

- Network creation and percolation checks.
- General material dataclasses.
- General discharge protocols.
- General steady/transient solvers.
- Result objects, analysis, plotting, and IO.
- Presets that are useful beyond a single validation script.

Paper-specific code belongs under `examples/reproduce_khan2021/`:

- `validate_paper_figures.py`
- `parallel_validation.py`
- Figure extraction/comparison logic.
- Hardcoded Figure 4, 6, and 7 comparison settings.
- Khan-specific plotting style and validation metrics.

Paper reference files can remain in the repository but should not be installed as package data by default:

- `DERIVATION.md`
- `PAPER_REFERENCE.md`
- `AUDIT_FINAL.md`
- `data/paper_figures/`
- `data/validation/`

### Configuration System

Use Python dataclasses as the primary configuration API.

Reasons:

- Lightweight and explicit.
- Easy to type check.
- Easy to document.
- No required runtime dependency on Pydantic.
- Works naturally with scientific callables for OCV, diffusion, and conductivity.

Add optional YAML/JSON loading as IO helpers, not as the core API.

```python
from pnmcathode.io import load_simulation

sim = load_simulation("configs/nmc532_1c.yaml")
result = sim.run()
```

Example YAML should deserialize into the same dataclasses:

```yaml
cathode:
  geometry:
    shape: [10, 10, 10]
    spacing: 1.0e-5
    porosity: 0.35
    cbd_fraction: 0.10
    seed: 42
  active_material: nmc532_khan2021
  electrolyte: khan2021
kinetics:
  k0: 5.0e-10
separator:
  enabled: true
  thickness: 25.0e-6
  porosity: 0.39
protocol:
  c_rate: 1.0
  cutoff_voltage: 2.5
settings:
  temperature: 298.15
  max_steps: 500
```

Keep YAML support optional:

```text
pip install pnmcathode[yaml]
```

### Handling Khan-Specific Constants

Move Khan-specific constants into presets:

```python
# src/pnmcathode/materials/presets.py
from pnmcathode.materials import ActiveMaterial, Electrolyte, Separator, ConductiveAdditive
from pnmcathode.physics.ocv import nmc532_ocv, ocv_derivative
from pnmcathode.physics.solid import nmc532_diffusion_coefficient
from pnmcathode.physics.electrolyte import (
    electrolyte_diffusion_coefficient,
    electrolyte_ionic_conductivity,
)

def nmc532_khan2021() -> ActiveMaterial:
    return ActiveMaterial(
        name="NMC532 (Khan 2021)",
        cs_max=48900.0,
        sigma=0.01,
        diffusivity=nmc532_diffusion_coefficient,
        ocv=nmc532_ocv,
        ocv_derivative=ocv_derivative,
    )

def cbd_khan2021() -> ConductiveAdditive:
    return ConductiveAdditive(name="CBD (Khan 2021)", sigma=760.0)

def electrolyte_khan2021() -> Electrolyte:
    return Electrolyte(
        name="LiPF6 carbonate electrolyte (Khan 2021)",
        c_init=1200.0,
        diffusivity=electrolyte_diffusion_coefficient,
        conductivity=electrolyte_ionic_conductivity,
        transference_number=0.363,
        bruggeman=1.5,
    )

def separator_khan2021(enabled: bool = True) -> Separator:
    return Separator(
        enabled=enabled,
        thickness=25e-6,
        porosity=0.39,
        bruggeman=1.5,
        t_plus=0.363,
        c_ref=1200.0,
    )
```

The default convenience API may use Khan presets for backward compatibility, but documentation must call them presets, not universal defaults.

## C. API Examples

### 1. Simple: 1C Discharge on Default NMC532 Cathode

```python
from pnmcathode import Cathode, DischargeProtocol, Simulation
from pnmcathode.materials.presets import nmc532_khan2021, electrolyte_khan2021

cathode = Cathode.cubic(
    active_material=nmc532_khan2021(),
    electrolyte=electrolyte_khan2021(),
)

result = Simulation(
    cathode=cathode,
    protocol=DischargeProtocol(c_rate=1.0, cutoff_voltage=2.5),
).run()

print(f"Capacity: {result.final_capacity:.4f} Ah/m2")
print(f"Final voltage: {result.voltage[-1]:.3f} V")
```

### 2. Medium: Custom Cathode Geometry and Electrolyte, Sweep 0.2C-3C

```python
import matplotlib.pyplot as plt

from pnmcathode import Cathode, CathodeGeometry, DischargeProtocol, Simulation
from pnmcathode.materials import Electrolyte
from pnmcathode.materials.presets import nmc532_khan2021
from pnmcathode.plotting import plot_voltage_capacity

def custom_diffusivity(c_e: float, T: float) -> float:
    return 1.8e-10

def custom_conductivity(c_e: float, T: float) -> float:
    return 1.25

geometry = CathodeGeometry(
    shape=(12, 12, 12),
    spacing=8e-6,
    porosity=0.42,
    cbd_fraction=0.08,
    seed=7,
    throat_scale=1.2,
)

electrolyte = Electrolyte(
    name="constant-property electrolyte",
    c_init=1000.0,
    diffusivity=custom_diffusivity,
    conductivity=custom_conductivity,
    transference_number=0.38,
)

cathode = Cathode.cubic(
    geometry=geometry,
    active_material=nmc532_khan2021(),
    electrolyte=electrolyte,
)

fig, ax = plt.subplots()
for c_rate in [0.2, 0.5, 1.0, 2.0, 3.0]:
    result = Simulation(
        cathode=cathode,
        protocol=DischargeProtocol(c_rate=c_rate, cutoff_voltage=2.7),
    ).run()
    plot_voltage_capacity(result, ax=ax, label=f"{c_rate:g}C")

ax.set_xlabel("Capacity / Ah m$^{-2}$")
ax.set_ylabel("Voltage / V")
ax.legend()
plt.show()
```

Important implementation detail: each sweep run should start from a fresh solver state. `Simulation.run()` must either clone the cathode state internally or document that users should call `cathode.copy()` when reusing mutable networks.

### 3. Advanced: Custom OCV Curve, Custom Separator, Spatial Output

```python
import numpy as np

from pnmcathode import (
    Cathode,
    CathodeGeometry,
    DischargeProtocol,
    Kinetics,
    Simulation,
    SolverSettings,
)
from pnmcathode.materials import ActiveMaterial, Electrolyte, Separator
from pnmcathode.materials.presets import electrolyte_khan2021

def lfp_ocv(soc: float | np.ndarray) -> float | np.ndarray:
    soc = np.asarray(soc, dtype=float)
    return 3.42 - 0.08 * np.tanh((soc - 0.5) / 0.08)

def lfp_diffusivity(c_s: float, T: float) -> float:
    return 3.0e-15

lfp = ActiveMaterial(
    name="Example LFP",
    cs_max=22860.0,
    sigma=0.05,
    diffusivity=lfp_diffusivity,
    ocv=lfp_ocv,
)

separator = Separator(
    enabled=True,
    thickness=20e-6,
    porosity=0.45,
    bruggeman=1.5,
    t_plus=0.40,
    c_ref=1000.0,
    include_concentration_overpotential=True,
    include_li_foil_bv=False,
)

cathode = Cathode.cubic(
    geometry=CathodeGeometry(
        shape=(16, 8, 8),
        spacing=6e-6,
        porosity=0.38,
        cbd_fraction=0.12,
        seed=11,
    ),
    active_material=lfp,
    electrolyte=electrolyte_khan2021(),
)

result = Simulation(
    cathode=cathode,
    protocol=DischargeProtocol(
        c_rate=2.0,
        cutoff_voltage=2.8,
        initial_soc=0.35,
        save_spatial=True,
        spatial_interval=5,
    ),
    kinetics=Kinetics(k0=1.0e-10, alpha_a=0.5, alpha_c=0.5),
    separator=separator,
    settings=SolverSettings(temperature=303.15, dt_max=60.0),
).run()

result.to_npz("outputs/lfp_2c_spatial.npz")

final_sol = result.spatial["state_of_lithiation"][-1]
final_ce = result.spatial["c_e"][-1]
print(f"Final capacity: {result.final_capacity:.4f} Ah/m2")
print(f"Mean final SoL: {np.nanmean(final_sol):.3f}")
print(f"Minimum final electrolyte concentration: {np.nanmin(final_ce):.1f} mol/m3")
```

## D. Packaging

### `pyproject.toml` Structure

```toml
[build-system]
requires = ["hatchling>=1.25"]
build-backend = "hatchling.build"

[project]
name = "pnmcathode"
version = "0.1.0"
description = "Pore-network cathode discharge simulation toolkit"
readme = "README.md"
requires-python = ">=3.10"
license = { text = "MIT" }
authors = [
  { name = "PNM Cathode contributors" }
]
keywords = [
  "battery",
  "lithium-ion",
  "pore-network-model",
  "cathode",
  "simulation",
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
]

[project.optional-dependencies]
dev = [
  "pytest>=7",
  "pytest-cov>=4",
  "ruff>=0.5",
  "mypy>=1.8",
]
yaml = [
  "pyyaml>=6",
]
docs = [
  "mkdocs>=1.5",
  "mkdocs-material>=9",
]

[project.urls]
Homepage = "https://github.com/OWNER/pnm-lib-cathode"
Repository = "https://github.com/OWNER/pnm-lib-cathode"
Issues = "https://github.com/OWNER/pnm-lib-cathode/issues"

[project.scripts]
pnmcathode-discharge = "pnmcathode.cli:discharge_main"

[tool.hatch.build.targets.wheel]
packages = ["src/pnmcathode"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-q"

[tool.ruff]
line-length = 100
target-version = "py310"

[tool.ruff.lint]
select = ["E", "F", "I", "UP", "B"]

[tool.mypy]
python_version = "3.10"
warn_return_any = true
warn_unused_configs = true
ignore_missing_imports = true
```

### Include in the Package

Install these:

- `pnmcathode.network`
- `pnmcathode.materials`
- `pnmcathode.physics`
- `pnmcathode.solver`
- `pnmcathode.simulation`
- `pnmcathode.results`
- `pnmcathode.analysis`
- `pnmcathode.plotting`
- `pnmcathode.io`
- CLI entry points that are general, not paper-validation-specific.

Do not install these as runtime package modules:

- `scripts/validate_paper_figures.py`
- `scripts/parallel_validation.py`
- `data/paper_figures/`
- `data/validation/`
- temporary probe files

Move general runnable examples to `examples/`:

- `examples/simple_1c.py`
- `examples/sweep_c_rates.py`
- `examples/custom_models.py`

Move paper reproduction scripts to:

- `examples/reproduce_khan2021/run_discharge.py`
- `examples/reproduce_khan2021/validate_paper_figures.py`
- `examples/reproduce_khan2021/parallel_validation.py`

### Test Suite Handling

Keep all 88 existing tests and add compatibility tests during migration.

Test categories:

1. Existing kernel tests
   - Continue testing low-level physics, network generation, steady solver, transient solver, separator, and post-processing.

2. Public API tests
   - `tests/test_public_api.py`
   - Verify imports from `pnmcathode`.
   - Verify `Cathode.cubic()` builds a valid network.
   - Verify a short `Simulation(...).run()` returns `DischargeResult`.

3. Preset tests
   - `tests/test_material_presets.py`
   - Verify `nmc532_khan2021().cs_max == 48900.0`.
   - Verify preset OCV and electrolyte functions match existing functions.

4. Compatibility tests
   - Verify current imports keep working temporarily:

```python
def test_legacy_imports_still_work():
    from src.network.generator import create_cathode_network
    from src.solver.transient import TransientSolver

    assert create_cathode_network is not None
    assert TransientSolver is not None
```

5. Packaging smoke tests
   - Run `python -m build`.
   - Install the wheel into a clean virtual environment.
   - Run:

```bash
python - <<'PY'
from pnmcathode import Cathode, DischargeProtocol, Simulation
from pnmcathode.materials.presets import nmc532_khan2021, electrolyte_khan2021

cathode = Cathode.cubic(
    active_material=nmc532_khan2021(),
    electrolyte=electrolyte_khan2021(),
)
result = Simulation(cathode, DischargeProtocol(c_rate=0.2, max_steps=2)).run()
print(result.voltage[-1])
PY
```

## E. Migration Strategy

### Refactor Without Breaking Existing Tests

Use a compatibility-first migration.

#### Phase 1: Create Package Shell

Create `src/pnmcathode/` and copy or move modules under it.

Option A, lowest risk:

- Keep existing `src/network`, `src/physics`, `src/solver`, `src/post`.
- Add `src/pnmcathode/` modules that import from the existing modules.
- Add `pyproject.toml`.
- Add public dataclasses and presets.

Example compatibility wrapper:

```python
# src/pnmcathode/network/generator.py
from src.network.generator import check_percolation, create_cathode_network

__all__ = ["check_percolation", "create_cathode_network"]
```

Then update tests to exercise both old and new imports.

#### Phase 2: Add Public Simulation Layer

Implement:

- `pnmcathode.config`
- `pnmcathode.materials`
- `pnmcathode.materials.presets`
- `pnmcathode.results`
- `pnmcathode.simulation`
- `pnmcathode.analysis`
- `pnmcathode.plotting`

At this phase, `Simulation.run()` can still call the existing `TransientSolver`:

```python
class Simulation:
    def run(self) -> DischargeResult:
        net = self.cathode.network or self.cathode.build_network()

        solver = TransientSolver(
            net,
            T=self.settings.temperature,
            k0=self.kinetics.k0,
            geometric_area=self.cathode.geometry.geometric_area,
            separator=self.separator.to_solver_params(),
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
            save_spatial=self.protocol.save_spatial,
            spatial_interval=self.protocol.spatial_interval,
        )
        return DischargeResult.from_solver_dict(raw, metadata=self.metadata())
```

This is not fully material-generic yet, but it creates the user-facing contract.

#### Phase 3: Inject Material Models Into Solvers

Change solver constructors from scalar parameters to model objects while preserving old arguments:

```python
class TransientSolver:
    def __init__(
        self,
        net,
        T: float = 298.15,
        k0: float | None = None,
        active_material: ActiveMaterial | None = None,
        electrolyte: Electrolyte | None = None,
        kinetics: Kinetics | None = None,
        separator: SeparatorParams | None = None,
        geometric_area: float | None = None,
    ):
        self.active_material = active_material or nmc532_khan2021()
        self.electrolyte = electrolyte or electrolyte_khan2021()
        self.kinetics = kinetics or Kinetics(k0=5e-10 if k0 is None else k0)
```

Then replace hardcoded calls:

```python
# before
U_eq = nmc532_ocv(cs / 48900.0)
i0 = exchange_current_density(self.k0, ce, cs, 48900.0)
D_s1 = nmc532_diffusion_coefficient(self.c_s[g1], self.T)

# after
cs_max = self.active_material.cs_max
U_eq = self.active_material.ocv(cs / cs_max)
i0 = exchange_current_density(
    self.kinetics.k0,
    ce,
    cs,
    cs_max,
    alpha_a=self.kinetics.alpha_a,
    alpha_c=self.kinetics.alpha_c,
)
D_s1 = self.active_material.diffusivity(self.c_s[g1], self.T)
```

Replace electrolyte transport:

```python
# before
D_e = electrolyte_diffusion_coefficient(c_e_mean, self.T)
kappa = electrolyte_ionic_conductivity(c_e_mean, self.T)

# after
D_e = self.electrolyte.diffusivity(c_e_mean, self.T)
kappa = self.electrolyte.conductivity(c_e_mean, self.T)
```

Replace conductivity:

```python
# before
sigma_nmc, sigma_cbd = 0.01, 760.0

# after
sigma_nmc = self.active_material.sigma
sigma_cbd = self.conductive_additive.sigma
```

#### Phase 4: Move Paper Scripts

Move scripts into examples and update imports.

Before:

```python
from src.network.generator import create_cathode_network
from src.solver.transient import TransientSolver
```

After:

```python
from pnmcathode import Cathode, CathodeGeometry, DischargeProtocol, Simulation
from pnmcathode.materials.presets import (
    cbd_khan2021,
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)
```

Keep a general CLI:

```bash
pnmcathode-discharge --crate 1.0 --shape 10x10x10 --preset khan2021 --output data/1c.npz
```

#### Phase 5: Documentation and Release

Update:

- `README.md` quickstart to use `pnmcathode`.
- `examples/README.md`.
- API reference docstrings.
- Changelog.

Release checklist:

```bash
python -m pip install -e ".[dev,yaml]"
python -m pytest
python -m build
python -m pip install dist/pnmcathode-0.1.0-py3-none-any.whl --force-reinstall
pnmcathode-discharge --crate 0.2 --shape 5x5x5 --max-steps 2
```

### Recommended Implementation Order

1. Add `pyproject.toml` and package metadata.
2. Add `src/pnmcathode/__init__.py` and import wrappers.
3. Add dataclasses in `config.py` and `materials/`.
4. Add Khan presets.
5. Add `DischargeResult`.
6. Add `Simulation` orchestration wrapper.
7. Add tests for public imports, presets, and a short simulation.
8. Update README quickstart.
9. Move general scripts to `examples/`.
10. Move paper-validation scripts to `examples/reproduce_khan2021/`.
11. Inject material models into `SteadyStateSolver` and `TransientSolver`.
12. Update tests to validate custom OCV, custom electrolyte, and custom separator.
13. Deprecate legacy `src.*` imports with warnings.
14. In a later release, remove legacy wrappers if desired.

### Estimated Effort

| Phase | Scope | Estimate |
|---|---|---:|
| 1 | Package shell, pyproject, wrappers | 0.5-1 day |
| 2 | Dataclasses, presets, result object, simulation wrapper | 1-2 days |
| 3 | Tests for public API and packaging smoke tests | 0.5-1 day |
| 4 | README and examples reorganization | 0.5 day |
| 5 | Material-model injection into solvers | 2-4 days |
| 6 | CLI cleanup and YAML loader | 1 day |
| 7 | Documentation polish and release checklist | 0.5-1 day |

Total: about 6-10 engineering days for a clean `0.1.0` release, depending on how much solver generalization is included before the first package cut.

### First Release Boundary

For `0.1.0`, prioritize a reliable package over a perfect abstraction.

Must have:

- `pip install -e .` works.
- `from pnmcathode import Cathode, Simulation` works.
- Existing tests still pass.
- A short public simulation returns a `DischargeResult`.
- Khan settings are represented as presets, not scattered script constants.
- Examples show simple, sweep, and custom-model workflows.

Can wait until `0.2.0`:

- Full YAML schema validation.
- Complete replacement of all hardcoded material calls inside solver internals.
- Non-cubic imported XCT network support.
- Full-cell simulation.
- Thermal coupling.
- Degradation models.
