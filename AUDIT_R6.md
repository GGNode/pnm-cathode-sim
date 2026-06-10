# AUDIT R6: Deep Audit - PNM Cathode Discharge Simulation

Date: 2026-06-10

Scope: `src/solver/transient.py`, `src/solver/steady.py`, `src/physics/*`, `src/network/generator.py`, and `scripts/validate_paper_figures.py`.

## Executive Finding

The voltage does not drop to the Khan 3.0 V cutoff because the simulation being compared is not applying the same physical problem as the paper. The dominant root causes are:

1. The synthetic 10x10x10 electrode has only about `127.8 g/m2` NMC loading, versus Khan 1CAL `297.8 g/m2`. Because C-rate is converted from the simulated active volume, 1C is `22.9 A/m2` instead of the paper-equivalent `~53.0 A/m2`. The 3C current is therefore `~68.7 A/m2` instead of `~159 A/m2`.
2. The validation starts at `INITIAL_SOL = 0.35`, where the current Khan OCV implementation gives `U = 4.300 V`. The paper curves start near 4.0 V. With this OCV function, `x ~= 0.50` is the 4.0 V initial state.
3. Transport polarization is suppressed by synthetic geometry. For the current 10x10x10 setup, 58% of throats are clipped to `0.1 um` length, and `throat_scale=2.0` quadruples throat/interface area. This produces very high ionic/electronic conductance and too little electrolyte drop.
4. `TransientSolver.step()` returns voltage from the pre-update concentration state while recording capacity at the post-update time. This can make V-Q curves look flatter or delayed, especially with large time steps or when a run stops at the time/SoL limit.
5. The model omits the paper's separator voltage model and electrolyte concentration-potential terms, so even with matching current density it will underpredict high-rate voltage loss.

## Bug List

### P0 - V-Q State Is Time-Misaligned

`src/solver/transient.py:564-567` solves potentials using the old concentrations. The code then updates `c_e`/`c_s` at `src/solver/transient.py:659-663`, but returns `pot["voltage"]` at `src/solver/transient.py:670`.

This means each reported point pairs end-of-step capacity with start-of-step voltage. `run_discharge()` records `q_new = |I| * t_new` at `src/solver/transient.py:879-882`, so the plotted voltage is stale relative to capacity. Fix by either:

- reporting voltage at `t_n` with capacity at `t_n`, or
- recomputing the potential after concentration update before returning/reporting the end-of-step voltage.

Expected impact: fixes misleadingly flat/delayed V-Q curves and cutoff detection near the end of discharge. It will not by itself create the missing paper-scale polarization.

### P0 - C-rate Is Correct Internally, Wrong for Paper Comparison

`current_density_for_c_rate()` uses simulated remaining capacity:

- capacity from NMC vacancies: `src/solver/transient.py:218-240`
- current density: `src/solver/transient.py:280-286`

For the 10x10x10 setup with `A = (10 * 1e-5)^2`, `porosity=0.368`, `cbd=0.1392`, `throat_scale=2`, and `INITIAL_SOL=0.35`, the probe gives:

| Quantity | Current model | Khan 1CAL |
|---|---:|---:|
| NMC mass loading | `127.8 g/m2` | `297.8 g/m2` |
| 1C current density | `22.9 A/m2` | `178 mAh/g * 297.8 g/m2 = 53.0 A/m2` |
| 3C current density | `68.7 A/m2` | `159 A/m2` |

So the C-rate conversion is not a sign or area bug when `geometric_area` is supplied, but it is not paper-equivalent because the simulated electrode loading is too low. Without `geometric_area`, `_cc_area` is a summed throat area (`src/solver/steady.py:289-297`), which is not a projected cell area and is unsafe for paper validation.

Expected impact: using paper areal loading/current should increase all polarization terms by about 2.3x before transport-geometry effects.

### P0 - Initial Stoichiometry Does Not Match Paper Starting Voltage

`scripts/validate_paper_figures.py:45-46` sets:

```python
INITIAL_SOL = 0.35
C_S0 = INITIAL_SOL * CS_MAX
```

But the current OCV gives:

- `U(0.35) = 4.300 V`
- `U(0.50) = 3.996 V`

The task's initial voltages (`4.286 V` at 0.2C, `4.159 V` at 3C) are consistent with starting at `x=0.35`. The paper reference starts near `~4.0 V`, which corresponds to `x ~= 0.5` in this code.

Expected impact: changing only `INITIAL_SOL` to 0.5 lowers the initial voltage by about 300 mV, but it also reduces simulated vacancy capacity unless the specific-capacity/current normalization is recalibrated. Do not tune this independently from the capacity basis.

### P1 - Synthetic Geometry Overestimates Conductance

`src/network/generator.py:195` sets pore diameters to `0.8-1.2 * spacing`. With `spacing=10 um`, pore diameters are `8-12 um`. `src/network/generator.py:219-221` computes throat length as center distance minus pore radii, then clips it to `0.01 * spacing = 0.1 um`.

For the audited setup:

- pore diameter: `8.00-12.00 um`, mean `10.03 um`
- throat diameter with `throat_scale=2`: `2.00-8.00 um`, mean `4.96 um`
- throat length: mean `0.36 um`
- clipped throat fraction: `57.6%`

This makes `G = property * A / L` far too large for many throats. The paper formulas use series pore/throat conductances (`PAPER_REFERENCE.md:94-120`, `122-129`), while the code uses only throat `A/L` in `src/solver/steady.py:159`, `src/solver/steady.py:218`, `src/solver/transient.py:417`, and `src/solver/transient.py:444`.

Expected impact: fixing geometry/series conductance should materially increase `phi_e` and solid/contact losses. This is a major reason the task reports `phi_e` span of only `14 mV` at 3C.

### P1 - `throat_scale=2.0` Is an Area Quadrupling, Not an Area Doubling

`src/network/generator.py:204` multiplies throat diameter by `throat_scale`; `src/network/generator.py:207` computes area from diameter squared. Therefore `throat_scale=2.0` makes throat area and reaction interface area 4x larger.

That same area is used for electrolyte conductance, solid conductance, solid/electrolyte diffusion, and BV interface current:

- electrolyte conductance: `src/solver/steady.py:159`
- solid conductance: `src/solver/steady.py:218`
- reaction interface list: `src/solver/steady.py:223-229`
- transient reaction interface list: `src/solver/transient.py:487-496`

Expected impact: reducing/removing this hack increases ohmic/activation losses. If the intent is 2x area, use `throat_scale=sqrt(2)`, not `2.0`.

### P1 - Separator and Electrolyte Migration Physics Are Missing

The paper includes separator concentration/potential equations and transference effects (`PAPER_REFERENCE.md:178-190`). The current code pins separator-side electrolyte concentration at the reservoir in `src/solver/transient.py:635-643` and pins separator electrolyte potential to zero in `src/solver/steady.py:469-475` and `src/solver/steady.py:597-606`.

The electrolyte charge equation is purely Ohmic:

- residual: `src/solver/steady.py:543-548`
- no concentration-potential term in the current
- no separator ohmic/concentration drop added to terminal voltage

Expected impact: high-rate voltage is too high, especially 1C/3C. The paper's 3C electrolyte potential variation is hundreds of mV; the audited 10x10x10 setup gives only about 14 mV.

### P1 - OCV Clipping Masks the 3.0 V Tail

`src/physics/ocv.py:59-60` clips OCV to `[3.0, 4.5]`. The raw corrected polynomial gives `U(1.0) ~= 2.94 V`, but the model cannot thermodynamically go below 3.0 V. With weak polarization, cutoff can be delayed or only reached by overpotential.

Expected impact: small at high current if polarization is correct, but important for low-rate cutoff behavior.

### P2 - Validation Stop/Diagnostics Use Unweighted Mean SoL

`scripts/validate_paper_figures.py:293` stops on the unweighted mean of NMC node stoichiometry. Metrics at `scripts/validate_paper_figures.py:317-319` also report unweighted mean/min/max. Because pore volumes vary, this is not the mass-average SoL and can disagree with delivered capacity.

Expected impact: misleading diagnostics and stop conditions. Use `sum(c_s * V) / (cs_max * sum(V))`.

## Answers to Critical Questions

### Q1: Why does voltage not drop?

Because the run has too little applied paper-equivalent current, starts at too high an OCV, has unrealistically large conductances/interface areas, omits separator/migration drops, and reports voltage one concentration step behind capacity.

The 3C initial solve is internally consistent: with the 10x10x10 network, `I_app ~= -68.7 A/m2`, terminal voltage is `4.159 V`, `phi_s` span is `236 mV`, `phi_e` span is `14 mV`, and total Faradaic current matches the imposed current. The problem is that paper 3C should be closer to `-159 A/m2`, plus separator/electrolyte concentration losses.

### Q2: Is C-rate conversion correct?

Correct for the simulated active volume when `geometric_area` is set. Incorrect as a paper comparison unless the simulated mass loading matches the paper or the current density is overridden from paper areal capacity.

Recommended paper-equivalent formula:

```text
I_1C = 178 Ah/kg * mass_loading_kg_m2 / 1 h
     = 178 mAh/g * 297.8 g/m2 / 1000
     ~= 53.0 A/m2
```

Use `I_app = -C_rate * 53.0 A/m2` for Khan 1CAL validation, or construct a network whose NMC loading is `~297.8 g/m2`.

### Q3: Is terminal voltage formula correct?

Mostly yes in principle. `src/solver/steady.py:809-812` computes:

```python
V_cell = mean(phi_s[active collector solid]) - mean(phi_e[separator electrolyte])
```

That matches the half-cell definition if `phi_e(separator)` is the Li reference and the separator drop is explicitly zero. It is not the same as `max(phi_s) - min(phi_e)`, and the full `phi_s` span across all solid nodes should not be interpreted as terminal voltage loss.

The missing piece is separator voltage/concentration drop, not the collector-minus-separator sign.

### Q4: Are diffusion matrices correct?

The signs of the Laplacians are consistent with `L*c = sum K(c_neighbor - c_i)`. The weak points are physical fidelity:

- electrolyte diffusion uses fixed `0.39**1.5` at `src/solver/transient.py:415`, not the actual cathode porosity and possibly double-counts tortuosity already represented by the network;
- electrolyte transport omits migration/upwind terms from the paper;
- conductance uses one throat `A/L`, not pore-throat-pore series resistance;
- solid diffusion is explicit in `D_s(c_s)` per step, not fully implicit in concentration.

### Q5: Is the steady solver coupling correct?

The BV residual/Jacobian signs are correct:

- electrolyte residual gets `+i_rxn A`: `src/solver/steady.py:544-548`
- solid residual gets `-i_rxn A`: `src/solver/steady.py:550-553`
- Jacobian terms follow `eta = phi_s - phi_e - U`: `src/solver/steady.py:563-590`

The comment at `src/solver/steady.py:36` says `J_se = -g_bv`, but the code uses `+g_bv`, which is correct for the implemented residual. This is a documentation bug, not a physics bug.

Boundary conditions are reasonable for the simplified model: separator `phi_e=0`, collector current on solid. They are incomplete relative to Khan because the separator model is collapsed to a pinned boundary.

### Q6: What about `throat_scale=2.0`?

It is applied globally and consistently, but it is not physically controlled. It multiplies diameters, so areas become 4x. It also increases reaction area, not only transport area. This lowers ohmic and activation polarization simultaneously, which directly works against voltage drop.

### Q7: Pore radius distribution?

The task description says mean radius `0.5 um`, but the current generator does not do that. With `spacing=10 um`, `src/network/generator.py:195` produces pore radii of `4-6 um`. This overlaps neighboring pores and triggers throat-length clipping. The issue is not pores being too small; it is that the simple cubic geometry creates many near-zero throat lengths and too-low mass loading for the paper thickness.

## Parameter Gap Analysis

| Parameter | Current validation | Khan 1CAL | Impact |
|---|---:|---:|---|
| Thickness | coordinate span `90 um` | `129 um` | Lower mass loading and shorter transport path |
| NMC loading | `127.8 g/m2` for 10x10x10 probe | `297.8 g/m2` | Current density at a paper C-rate is 2.3x too low |
| Initial SoL | `0.35` | voltage suggests `~0.5` under current OCV | Initial voltage too high by ~0.3 V |
| 1C current density | `22.9 A/m2` | `~53.0 A/m2` | Too little polarization |
| Throat scale | diameter x2, area x4 | XCT-derived network | Suppresses ohmic and BV losses |
| Separator | pinned reservoir/reference | finite 25 um separator, t+ effects | Missing voltage and concentration losses |
| Electrolyte charge | Ohmic only | includes concentration effects | Underpredicts high-rate phi_e loss |
| OCV lower bound | clipped at 3.0 V | raw curve reaches below 3.0 near full lithiation | Can mask cutoff tail |

## Execution Plan

1. Fix V-Q time alignment in `TransientSolver.step()` and validation plotting. Recompute/report end-of-step voltage after concentration update, or store start-of-step capacity with start voltage.
2. Add a validation mode that uses paper current density directly: `I_1C = 53.0 A/m2` for 1CAL. Do not call this "C-rate from synthetic network capacity".
3. Reconcile initial stoichiometry and capacity basis. Try `INITIAL_SOL=0.50` for voltage matching, but pair it with a paper-specific capacity/current normalization.
4. Remove `throat_scale=2.0` for the baseline. If a throat correction is retained, calibrate conductance and reaction area separately; do not scale every throat area globally.
5. Replace throat-only conductance with pore-throat-pore series conductance and avoid artificial `0.1 um` throat lengths dominating the network.
6. Add separator voltage drop at minimum as an analytic 1D correction, then implement the paper's separator concentration/potential equations.
7. Add migration/concentration-potential terms to electrolyte current and the paper's electrolyte concentration transport terms.
8. Change validation diagnostics to volume-weighted SoL and add a consistency check: `delivered_capacity` vs `sum(V_s * delta c_s) * F`.

## Recommended Parameters to Try

For a paper-comparison run before geometry refactor:

```python
INITIAL_SOL = 0.50
K0 = 1e-10
T = 303.0
POROSITY = 0.368
CBD_FRACTION = 0.1392
THROAT_SCALE = 1.0
I_1C_PAPER_1CAL = 53.0  # A/m2
I_app = -C_rate * I_1C_PAPER_1CAL
CUTOFF = 3.0
```

For a geometry-consistent 1CAL-like synthetic network:

```python
spacing = 398e-9 or another voxel/network-calibrated length scale
thickness ~= 129e-6
geometric_area chosen so NMC mass_loading ~= 297.8 g/m2
use XCT-derived pore/throat radii if available
use no global throat_scale until conductance calibration is measured
```

Expected impact order:

1. Paper current density and initial SoL: fixes starting voltage/current scale.
2. Remove area quadrupling and throat clipping: increases internal polarization.
3. Separator/migration model: required for the 3C drop and spatial phi_e profiles.
4. V-Q time alignment: makes cutoff and capacity reporting trustworthy.

