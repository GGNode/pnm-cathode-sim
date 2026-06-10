# PNM-LIB-Cathode

**Pore Network Model for Li-ion Battery Cathode Discharge**

Reproduction of: Khan ZA, Elkamel A, Gostick JT. *"Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes."* J. Electrochem. Soc. 168(7), 070534 (2021).

## Overview

This project implements a transient pore network model (PNM) that simulates galvanostatic discharge of an NMC532 lithium-ion battery cathode. The model discretizes the porous electrode into a graph of pore nodes (storing species and potentials) connected by throat bonds (conducting fluxes), and couples:

- **Electrolyte salt transport** — diffusion and migration of Li⁺ in the liquid phase
- **Solid-phase lithium diffusion** — Fickian diffusion in spherical NMC532 particles
- **Butler-Volmer kinetics** — interfacial reaction at electrolyte/NMC boundaries
- **Charge conservation** — Ohmic conduction in both electrolyte and solid phases
- **Separator boundary model** — collapsed 1D separator with ohmic/concentration losses

## Project Structure

```
pnm-lib-cathode/
├── src/                          # Core library
│   ├── __init__.py
│   ├── network/
│   │   └── generator.py          # Cubic pore network with 3-phase labeling
│   ├── physics/
│   │   ├── electrolyte.py        # D_e(c_e,T) and kappa(c_e,T) correlations
│   │   ├── ocv.py                # NMC532 OCV polynomial (Khan Eq. 2.19)
│   │   ├── reaction.py           # Butler-Volmer kinetics and exchange current
│   │   ├── separator.py          # Collapsed 1D separator boundary model
│   │   └── solid.py              # Spherical solid diffusion (FVM shells)
│   ├── solver/
│   │   ├── steady.py             # Coupled phi_e/phi_s Newton-Raphson solver
│   │   ├── transient.py          # Adaptive backward Euler discharge solver
│   │   └── single_pore.py        # Single-particle model (verification)
│   └── post/
│       ├── analysis.py           # Pore size distribution, SoL, discharge curves
│       └── visualization.py      # Matplotlib plotting utilities
├── scripts/
│   ├── run_discharge.py          # CLI entry point for single discharge
│   ├── validate_paper_figures.py # Paper figure comparison (Figures 4, 6, 7)
│   ├── parallel_validation.py    # Multi-C-rate parallel validation
│   └── plot_results.py           # Quick V-Q plot generator
├── tests/                        # pytest test suite (88 tests)
├── data/
│   ├── paper_figures/            # Cropped paper figure PNGs
│   └── validation/               # Pre-computed validation results
├── DERIVATION.md                 # Full mathematical derivations
├── PAPER_REFERENCE.md            # Paper data extraction
└── AUDIT_FINAL.md                # Final audit report
```

## Installation

```bash
# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

**Dependencies:** numpy, scipy, matplotlib, openpnm (≥4.0), pytest

## Quick Start

Run a 0.2C discharge on a 5×5×5 network in 5 lines:

```python
from src.network.generator import create_cathode_network
from src.solver.transient import TransientSolver

net = create_cathode_network(shape=[5, 5, 5], spacing=1e-5, porosity=0.35, seed=42)
solver = TransientSolver(net, T=298.15, k0=5e-10)
result = solver.run_discharge(C_rate=0.2, cutoff_voltage=3.0)

print(f"Capacity: {result['capacity_Ah_m2'][-1]:.3f} Ah/m²")
print(f"Final voltage: {result['voltage'][-1]:.3f} V")
```

Or from the command line:

```bash
python scripts/run_discharge.py --crate 0.2 --shape 5x5x5 --cutoff 3.0
```

## Usage Guide

### Single Discharge (`scripts/run_discharge.py`)

```bash
python scripts/run_discharge.py --crate 1.0 --shape 10x10x10 --cutoff 2.5 --output data/1c.npz
```

Key options: `--crate`, `--cutoff`, `--shape`, `--spacing`, `--porosity`, `--k0`, `--temperature`, `--output`

### Paper Validation (`scripts/validate_paper_figures.py`)

Runs 0.2C/0.5C/1C/3C discharges on a 10×10×10 network, generates comparison plots against paper Figures 4, 6, and 7:

```bash
python scripts/validate_paper_figures.py
```

### Parallel Validation (`scripts/parallel_validation.py`)

Runs 4 C-rates × 4 scenarios (paper/network current basis, separator on/off) in parallel on a 13×13×13 network:

```bash
python scripts/parallel_validation.py
```

## Physics Model Summary

### Governing Equations

| Equation | Description | Reference |
|----------|-------------|-----------|
| `∂c_e/∂t = ∇·(D_e,eff ∇c_e) + S_e` | Electrolyte salt conservation | §3.2 |
| `∇·(κ_eff ∇φ_e) = a_s·i_F` | Electrolyte charge conservation | §2.2, §3.3 |
| `∂c_s/∂t = ∇·(D_s ∇c_s)` | Solid lithium diffusion (spherical) | §2.3, §3.4 |
| `∇·(σ_eff ∇φ_s) = -a_s·i_F` | Solid charge conservation | §2.4, §3.5 |
| `i_r = i0·[exp(α_a f η) - exp(-α_c f η)]` | Butler-Volmer kinetics | §2.5 |
| `η = φ_s - φ_e - U_eq(c_s)` | Overpotential | §2.5 |

### Key Material Properties (Khan et al. Table II)

| Property | Symbol | Value | Unit |
|----------|--------|-------|------|
| Max Li concentration (NMC532) | c_s,max | 48,900 | mol/m³ |
| NMC conductivity | σ_nmc | 0.01 | S/m |
| CBD conductivity | σ_cbd | 760 | S/m |
| Rate constant | k0 | 5×10⁻¹⁰ | m²·⁵/(mol⁰·⁵·s) |
| Exchange current density | i0 | ~10⁻³ | A/m² |
| Li diffusion in NMC | D_s | ~10⁻¹⁵ | m²/s |
| Electrolyte diffusivity | D_e | ~10⁻¹⁰ | m²/s |

### Discretization

- **Network**: Cubic lattice with 3-phase random labeling (electrolyte/NMC/CBD)
- **Solid diffusion**: N-shell finite volume method on spherical particles
- **Potentials**: Coupled Newton-Raphson with backtracking line search
- **Time stepping**: Backward Euler with adaptive dt (growth/shrink/reject)

## Configuration Reference

All tunable parameters (with defaults):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `shape` | [5,5,5] | Network dimensions [nx, ny, nz] |
| `spacing` | 1e-5 m | Pore-to-pore spacing |
| `porosity` | 0.35 | Electrolyte volume fraction |
| `cbd_fraction` | 0.10 | CBD volume fraction |
| `T` | 298.15 K | Temperature |
| `k0` | 5e-10 | BV rate constant |
| `c_e_init` | 1200 mol/m³ | Initial electrolyte concentration |
| `c_s_init` | 24450 mol/m³ | Initial solid concentration (x≈0.5) |
| `cutoff_voltage` | 2.5 V | Discharge cutoff voltage |
| `dt_min` | 1e-6 s | Minimum adaptive time step |
| `dt_max` | 300 s | Maximum adaptive time step |

## Test Suite

```bash
# Run all tests
python -m pytest tests/ -q

# Run with verbose output
python -m pytest tests/ -v

# Run specific module
python -m pytest tests/test_transient.py -v
```

**Coverage (88 tests):**
- Butler-Volmer kinetics and exchange current density
- OCV polynomial and derivative
- Electrolyte transport correlations
- Solid diffusion (spherical FVM)
- Network generation and percolation
- Steady-state Newton-Raphson solver
- Transient solver with adaptive stepping
- Single-pore discharge verification
- Separator boundary model
- Post-processing analysis tools

## Known Limitations

1. **Synthetic network vs XCT**: The paper uses XCT-resolved cathode microstructure (1CAL: 4637 nodes, 31427 throats). This model uses synthetic cubic networks with random phase labeling. The geometry differs significantly, leading to:
   - Different mass loading and active area
   - Different phase connectivity and percolation
   - Quantitative capacity/polarization mismatch (especially at high C-rates)

2. **Mass loading gap**: Synthetic networks produce lower areal capacity than XCT networks due to different volume-to-area ratios and active site connectivity.

3. **Separator model**: Collapsed 1D quasi-steady model — adequate for mV-scale losses but not a full 2D/3D separator simulation.

4. **Isothermal**: Temperature is uniform and constant (no thermal coupling).

5. **No degradation**: SEI growth, particle cracking, and capacity fade are not modeled.

## Citation

If you use this code, please cite:

> Khan ZA, Elkamel A, Gostick JT. "Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes." *J. Electrochem. Soc.* 168(7), 070534 (2021).

## References

- [DERIVATION.md](DERIVATION.md) — Full mathematical derivations
- [PAPER_REFERENCE.md](PAPER_REFERENCE.md) — Paper data extraction
- [AUDIT_FINAL.md](AUDIT_FINAL.md) — Final audit report
