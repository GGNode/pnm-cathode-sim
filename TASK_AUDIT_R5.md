# Paper Audit: NMC532 Pore Network Model — Visual + Code Review

## Your Role

You are auditing a pore network model (PNM) implementation against its source paper. Your job is to be CRITICAL and THOROUGH — find every discrepancy, every missing feature, every parameter mismatch.

## Context

Project directory: `/VOLUMES/1TB/projects/pnm-lib-cathode`
Paper PDF: The paper is referenced in `DERIVATION.md`. Key figures have been extracted to `data/paper_figures/`.

## What to Do

### Step 1: Read the paper content
- Read `DERIVATION.md` (1325 lines) — this is the full mathematical derivation extracted from the paper
- Read `PAPER_EXTRACTION.md` — additional paper content
- Read `AUDIT.md` — previous audit (OUTDATED, based on old code)

### Step 2: Examine the paper figures visually
- Look at `data/paper_figures/figure4_crop_200dpi.png` — V-Q discharge curves
- Look at `data/paper_figures/figure6_crop_200dpi.png` — spatial distributions
- Look at `data/paper_figures/figure7_crop_200dpi.png` — pore profiles
- Look at `data/paper_figures/page_08_200dpi.png` through `page_11_200dpi.png` — full paper pages

Use the vision_analyze tool to examine each image. Describe what you see in detail:
- What are the axis labels and units?
- What C-rates are shown?
- What are the capacity values at each C-rate?
- What spatial patterns are visible?

### Step 3: Read ALL source code
Read every file in `src/`:
- `src/physics/ocv.py` — OCV curve
- `src/physics/reaction.py` — Butler-Volmer kinetics
- `src/physics/electrolyte.py` — electrolyte properties
- `src/physics/solid.py` — solid diffusion
- `src/network/generator.py` — network generation
- `src/solver/steady.py` — Newton-Raphson potential solver
- `src/solver/transient.py` — transient discharge solver
- `src/solver/single_pore.py` — single-pore validation
- `src/post/analysis.py` — post-processing
- `src/post/visualization.py` — plotting

### Step 4: Read tests
- All files in `tests/`

### Step 5: Run validation
```bash
cd /VOLUMES/1TB/projects/pnm-lib-cathode
.venv/bin/pytest -v 2>&1 | tail -30
```

Run a quick 5×5×5 discharge test:
```python
import sys; sys.path.insert(0, 'src')
from network.generator import create_cathode_network
from solver.transient import TransientSolver
import numpy as np

CS_MAX = 48900.0
net = create_cathode_network(shape=[5,5,5], spacing=1e-5, porosity=0.5, cbd_fraction=0.1, seed=42)
for C_rate in [0.2, 1.0, 3.0]:
    solver = TransientSolver(net, T=298.15, k0=5e-10)
    solver.set_concentration(c_e=1200.0, c_s=0.35*CS_MAX)
    I_app = solver.current_density_for_c_rate(C_rate)
    dt = 20.0
    for i in range(500):
        result = solver.step(dt=dt, I_app=I_app)
        if result['voltage'] <= 3.0:
            cap = abs(I_app) * (i+1) * dt / 3600.0
            print(f'{C_rate}C: cutoff at step {i+1}, cap_Ah_m2={cap:.2f}, V={result["voltage"]:.3f}')
            break
    else:
        cap = abs(I_app) * 500 * dt / 3600.0
        print(f'{C_rate}C: no cutoff after 500 steps, cap_Ah_m2={cap:.2f}, V={result["voltage"]:.3f}')
```

### Step 6: Produce the audit report

Write the report to `AUDIT_R5.md` with these sections:

#### A. Paper Figure Analysis
For each figure (4, 6, 7), describe:
- What the paper shows (axes, values, trends)
- What our simulation produces (from validation data)
- Exact quantitative comparison (capacities, voltages, spatial distributions)

#### B. Parameter Verification Table
For EVERY parameter in the paper's Table II (or equivalent), compare:
| Parameter | Paper Value | Our Code Value | File:Line | Match? | Impact |

#### C. Physics Implementation Audit
For each physical equation in DERIVATION.md:
| Equation | DERIVATION Section | Implemented? | File:Line | Correct? | Notes |

#### D. Missing Features
List everything the paper has that our code doesn't:
- Features, analyses, visualizations
- For each: difficulty estimate (low/medium/high) and impact on results

#### E. Discrepancy Analysis
For every difference between paper results and our results:
- Root cause (parameter? geometry? physics? numerics?)
- Suggested fix
- Priority (P0/P1/P2)

#### F. Next Steps — Detailed Implementation Plan
For each item above, provide:
- Concrete task description
- Files to modify
- Estimated effort (hours)
- Dependencies on other tasks
- Verification criteria

Order the tasks by priority and dependency.

## Output Format

Write the full report to `/VOLUMES/1TB/projects/pnm-lib-cathode/AUDIT_R5.md`.

Be SPECIFIC. Cite file paths and line numbers. Show exact numbers. Don't say "approximately" when you can give an exact value. This audit will be used to plan the next phase of work.
