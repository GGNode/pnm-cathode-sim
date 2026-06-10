# TASK: Deep Audit — PNM Cathode Discharge Simulation

## Goal
Audit the PNM cathode discharge simulation code for remaining physics bugs, parameter errors, and design issues that prevent matching the Khan et al. (2021) paper results.

## Current State
- **Tests**: 81 passed, 1 skipped
- **Network**: 10×10×10 synthetic cubic, 3-phase (electrolyte/NMC/CBD)
- **Solver**: Coupled Newton-Raphson for phi_e/phi_s, backward Euler for c_e/c_s

### Current Results (10×10×10, throat_scale=2.0, k0=1e-10, porosity=0.368)
| C-rate | V_init | V_final | phi_s span | phi_e span | Capacity | Cutoff hit? |
|--------|--------|---------|-----------|-----------|----------|------------|
| 0.2C | 4.286V | ~4.28V | 20 mV | 0.7 mV | 224* | ❌ (time limit) |
| 0.5C | 4.266V | ~4.26V | 49 mV | 1.9 mV | 224* | ❌ (time limit) |
| 1C | 4.219V | ~4.21V | 94 mV | 3.9 mV | 224* | ❌ (time limit) |
| 3C | 4.159V | ~4.14V | 236 mV | 14 mV | 224* | ❌ (time limit) |

*Time-limited (max_time=50000s), not cutoff-limited

### Paper Reference (Khan 2021, 1CAL)
| C-rate | Capacity | V_init | Notes |
|--------|----------|--------|-------|
| 0.2C | ~195 mAh/g | ~4.0V | Hits 3.0V cutoff |
| 0.5C | ~185 mAh/g | ~3.95V | Hits 3.0V cutoff |
| 1C | ~175 mAh/g | ~3.9V | Hits 3.0V cutoff |
| 3C | ~150 mAh/g | ~3.85V | Hits 3.0V cutoff |

## Critical Questions for Audit

### Q1: Why doesn't the voltage drop to 3.0V cutoff?
This is THE central mystery. At 3C, V_init=4.159V but after 224 mAh/g of capacity, V is still ~4.14V. The voltage barely moves!

In the paper, 3C starts at ~3.85V and drops to 3.0V after ~150 mAh/g. Our 3C starts higher (4.159V) but never drops.

Possible causes to investigate:
1. **Activation overpotential too small**: BV kinetics may not produce enough η to drop voltage
2. **phi_s span not translating to V drop**: The 236mV phi_s span at 3C is measured across ALL solid nodes. The terminal voltage is V = max(phi_s) - min(phi_e) - separator_drop. Is this right?
3. **Terminal voltage calculation bug**: How exactly is V_terminal computed? Check `transient.py` terminal voltage formula.
4. **C-rate conversion**: Is the applied current correct? Check `current_density_for_c_rate()`.
5. **NMC volume fraction**: With porosity=0.368 and CBD=0.1392, NMC fraction = 0.4928. The mass load is 127.8 g/m² vs paper's 297.8 g/m². This means our network has ~2.3× less NMC per unit area, so the same C-rate gives a lower current density. This could explain why voltage doesn't drop!

### Q2: Is the C-rate conversion correct?
Paper: 178 mAh/g × C_rate = current per gram of NMC
Our code: `current_density_for_c_rate()` computes I from `discharge_capacity_coulombs()` which sums over NMC node volumes.

Check: Does the code correctly compute the geometric current density [A/m²] from the C-rate? Or does it use a wrong area normalization?

### Q3: Is the terminal voltage formula correct?
In PNM, terminal voltage should be:
V = max(phi_s at collector) - min(phi_e at separator) - separator_voltage_drop

Or more precisely:
V = phi_s(collector_BC) - phi_e(separator_BC)

Check how V_terminal is actually computed in the code.

### Q4: Are the diffusion matrices correct?
- `_build_electrolyte_diffusion_matrix()`: Uses D_e(c_e) with correct Bruggeman?
- `_build_solid_diffusion_matrix()`: Uses D_s(c_s) with SoC-dependent formula?
- Are the conductance formulas G = D * A / L correct for the finite-volume scheme?

### Q5: Is the steady-state solver correctly coupling phi_e and phi_s?
- The Newton-Raphson couples phi_e and phi_s through BV kinetics
- Check: Does the Jacobian correctly include d(BV)/d(phi_e) and d(BV)/d(phi_s) terms?
- Check: Are the boundary conditions correct? (collector for phi_s, separator for phi_e)

### Q6: What about the throat_scale=2.0 hack?
The synthetic network has smaller throats than XCT. throat_scale=2.0 doubles throat area. Is this:
- Applied in the right place?
- Applied consistently to both electrolyte and solid throats?
- Physically reasonable?

### Q7: Pore radius distribution
Current: mean=0.5μm, std=0.1μm (narrow), scale=[0.8, 1.2]
Paper (1CAL): 0.9-4.5μm range (small), 4.5-8.0μm (medium), 8.0-12μm (large)
Our mean is 10× smaller than paper! This could affect throat conductances.

## Files to Review
1. `src/solver/transient.py` — TransientSolver, terminal voltage, C-rate, time stepping
2. `src/solver/steady.py` — SteadyStateSolver, Newton-Raphson, Jacobian, BCs
3. `src/physics/electrolyte.py` — D_e(c_e) and κ(c_e) correlations
4. `src/physics/solid.py` — D_s(c_s) SoC-dependent diffusion
5. `src/physics/ocv.py` — Khan OCV polynomial
6. `src/physics/reaction.py` — BV kinetics, exchange current density
7. `src/network/generator.py` — Network generation, phase assignment, geometry
8. `scripts/validate_paper_figures.py` — Validation script parameters

## Deliverables
Write `AUDIT_R6.md` with:
1. **Bug list**: Any physics/implementation bugs found, with exact file:line references
2. **Root cause analysis**: Why voltage doesn't drop to cutoff
3. **Parameter gap analysis**: Remaining parameter mismatches vs paper
4. **Execution plan**: Prioritized fix list with expected impact on each item
5. **Recommended parameter changes**: Specific values to try

Focus on ROOT CAUSES, not symptoms. The question is "why doesn't V drop?" not "how to make V drop".
