# Implementation Audit Against DERIVATION.md

Baseline verification: `.venv/bin/pytest -q` passes with `65 passed, 1 skipped`.

## Section 10 Checklist

| Item | Status | Evidence |
|---|---|---|
| Positive current convention documented | FAIL | `DERIVATION.md:92-101` defines `I_app > 0` as anodic charge and discharge as `I_app < 0`. `src/physics/reaction.py:35` says positive overpotential is anodic but incorrectly labels it discharge. `src/solver/steady.py:220` says positive `I_app` is discharge, contradicting the derivation. `src/solver/single_pore.py:110-111` correctly uses negative discharge current. |
| Discharge C-rate converted to `I_app < 0` | WARN | The single-pore model converts positive C-rate to negative current at `src/solver/single_pore.py:110-111`. The network discharge API accepts only raw `I_app` in `src/post/analysis.py:82-90`; there is no proper C-rate-to-current calculation from electrode capacity. |
| `eta = phi_s - phi_e - U` | PASS | Reaction documentation uses this form at `src/physics/reaction.py:8`; steady solve uses it at `src/solver/steady.py:164`, `src/solver/steady.py:289`, and `src/solver/steady.py:296`; transient uses it at `src/solver/transient.py:163`; single-pore uses it conceptually at `src/solver/single_pore.py:137-157`. |
| Electrolyte and solid reaction source signs opposite | FAIL | Derivation requires `F_c_e += -iA/F` and `F_c_s += +iA/F` in residual form, so state updates should use `dc_e/dt = +iA/(FV_e)` and `dc_s/dt = -iA/(FV_s)` (`DERIVATION.md:497-500`). The transient solver applies `dce_dt -= flux/V` and `dcs_dt += flux/V` at `src/solver/transient.py:169-171`, which reverses the derivation. The comments at `src/solver/transient.py:166-168` also incorrectly say `I_rxn > 0` is discharge. |
| `V_cell` references Li/Li+ correctly | PASS | With separator electrolyte pinned to zero, `V_cell = phi_s,collector - phi_e,separator` is allowed by `DERIVATION.md:651-657`. The implementation computes this at `src/solver/steady.py:274-276`. The single-pore solver reports `phi_s = U + eta` with `phi_e = 0` at `src/solver/single_pore.py:155-157`. |
| Exactly one potential gauge pinned | FAIL | Electrolyte separator nodes are pinned to zero at `src/solver/steady.py:225-229`, which is the intended reference boundary. The solid phase is also pinned at separator-side solid nodes at `src/solver/steady.py:236-239`, contradicting `DERIVATION.md:916-923` for galvanostatic operation. |
| Disconnected components removed or deactivated | FAIL | The derivation requires keeping only electrolyte components connected to the separator and reactive solid components connected to the collector (`DERIVATION.md:933-938`). The implementation only pins all-zero matrix rows at `src/solver/steady.py:230-233` and `src/solver/steady.py:240-242`; connected-but-floating components are not removed or deactivated. |
| BV exponentials stabilized | PASS | BV exponent arguments are clipped at `src/physics/reaction.py:51-56`. The steady solver applies matching clipping for its BV slope at `src/solver/steady.py:167-170`. |
| Concentration bounds enforced | WARN | `exchange_current_density` clips local values at `src/physics/reaction.py:92-95`, OCV clips lithiation at `src/physics/ocv.py:50-57`, and transient state updates clamp concentrations at `src/solver/transient.py:186-188`. This is bound-preserving after the step, not a rejected/damped solve as recommended by `DERIVATION.md:908-914`. |
| Analytical Jacobian signs consistent with residual signs | WARN | There is no full residual/Jacobian assembly despite the architecture target in `DERIVATION.md:1198-1229`. The steady solver uses a linearized BV conductance `g_bv = dI/deta A` at `src/solver/steady.py:167-171`; its electrolyte sign at `src/solver/steady.py:188-198` is consistent, and the solid diagonal code at `src/solver/steady.py:209-218` is consistent with `F_phi_s = L_s phi_s - iA + F_BC`, although nearby comments are inconsistent. |
| `kappa(c)` and `i0(c)` derivatives included if implicit | WARN | `kappa` is evaluated as a constant from 1200 mol/m3 at `src/solver/steady.py:57-62`, not implicitly as a function of `c_e`. `i0(c)` is evaluated from current concentrations at `src/solver/steady.py:156-163` and `src/solver/transient.py:149-164`, but concentration equations are explicit/semi-implicit (`src/solver/transient.py:180-184`), so the missing concentration derivatives are not currently used in a Newton system. |
| SI units throughout | WARN | Geometry uses meters and m3/m2 at `src/network/generator.py:27`, `src/network/generator.py:77-83`; concentrations are mol/m3 at `src/solver/steady.py:31-32`; constants are SI at `src/physics/reaction.py:16-18`. Capacity outputs are mixed: network post-processing reports `A·s/m²` at `src/post/analysis.py:114-129`, while single-pore reports Ah at `src/solver/single_pore.py:179-180`. |
| Bruggeman/tortuosity applied once | WARN | Electrolyte diffusion applies `0.39**1.5` at `src/solver/transient.py:73`. Electrolyte conductivity in steady solve applies no Bruggeman correction at `src/solver/steady.py:57-62`. There is no central transport builder documenting whether stochastic pore geometry already accounts for tortuosity, so double-counting cannot be ruled out globally. |
| `0C` recovers OCV with zero reaction current | WARN | The test suite checks near-OCV at low current at `tests/test_steady.py:67-80` and zero-current concentration conservation at `tests/test_transient.py:30-73`. There is no explicit assertion that `I_rxn == 0`, `eta == 0`, and `V_cell == U` at `I_app=0` as required by `DERIVATION.md:1029-1048`. |
| Small-current polarization scales linearly | FAIL | `tests/test_steady.py:67-80` checks only one low-current case. No test varies current magnitude and verifies linear scaling as required by `DERIVATION.md:1058-1079`. |
| High-rate discharge depletes electrolyte toward collector | FAIL | No implementation or test verifies the high-rate electrolyte depletion shape in `DERIVATION.md:1081-1123`. The transient test at `tests/test_transient.py:86-101` checks only that solid concentration changes and does not check electrolyte depletion location. |

## Additional Required Checks

| Check | Status | Evidence |
|---|---|---|
| Sign conventions: `eta = phi_s - phi_e - U` consistent | PASS | See checklist `eta` item above. |
| `V_cell` formula uses collector `phi_s`, not local `phi_s - phi_e` | PASS | `src/solver/steady.py:274-276` uses collector solid potential relative to separator electrolyte, and single-pore uses terminal `phi_s` at `src/solver/single_pore.py:155-157`. |
| Reaction source signs opposite in electrolyte vs solid | FAIL | Signs are opposite but reversed relative to the derivation in transient concentration updates at `src/solver/transient.py:169-171`. Potential reaction terms are opposite in the steady solver linearization at `src/solver/steady.py:188-218`. |
| Boundary conditions: exactly one gauge and disconnected handling | FAIL | Solid separator nodes are pinned under galvanostatic operation at `src/solver/steady.py:236-239`; disconnected connected components are not deactivated, only zero rows are pinned at `src/solver/steady.py:230-242`. |
| Jacobian derivatives consistent with residuals | WARN | Potential BV derivative signs are locally consistent in `src/solver/steady.py:167-218`, but the code has no full analytical concentration Jacobian or residual module. |
| Units SI throughout | WARN | Core physics is SI, but capacity reporting mixes `A·s/m²` and Ah as noted above. |
| Bruggeman applied once | WARN | Bruggeman appears only in electrolyte diffusion (`src/solver/transient.py:73`), not electrolyte conductivity (`src/solver/steady.py:57-62`), and no central policy prevents double-counting. |
| `0C` limit recovers OCV | WARN | Intended behavior is partially tested, but no explicit zero-current reaction/voltage test exists. |
| Mass conservation verified | FAIL | Existing tests verify zero-current mass conservation only (`tests/test_transient.py:30-73`). No test audits total lithium conservation at each nonzero-current time step as required by `DERIVATION.md:975-1006` and `DERIVATION.md:1300-1303`. |

## Fix Plan

1. Correct current-sign documentation and transient reaction concentration signs.
2. Remove the extra galvanostatic solid potential pin and handle disconnected ionic/electronic components explicitly.
3. Add conservation and limit tests for 0C, small-current scaling, and nonzero-current lithium conservation.
4. Add C-rate-driven network discharge, adaptive time stepping, and a plotting script for a 0.2C V-Q curve.
