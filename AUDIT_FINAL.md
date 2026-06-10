# Final Audit - Khan 2021 PNM Cathode Reproduction

## Verification

- Requested command `python -m pytest tests/ -v` was run and failed before collection because the default interpreter lacks `numpy`.
- Equivalent project-env run passed: `MPLCONFIGDIR="$PWD/data/.matplotlib" .venv/bin/python -m pytest tests/ -v`
  - Result: `88 passed, 1 skipped`.

## A. Remaining Code Bugs That Would Produce Wrong Physics

- No core steady/transient solver bug found that would explain the reported sep-off discharge results.
- Butler-Volmer sign convention is internally consistent:
  - `eta = phi_s - phi_e - U_eq`
  - `I_rxn < 0` for discharge/lithiation
  - electrolyte source decreases and solid source increases during discharge.
- Coupled Newton residual/Jacobian signs are consistent with anodic-current convention.
- Concentration-dependent `D_e(c_e)`, `kappa(c_e)`, and `D_s(c_s)` are wired into the relevant matrices.
- Real issue found in `scripts/parallel_validation.py`, not the core solver:
  - Separator-on cases mutate `solver._steady.separator` after construction instead of passing `separator=SeparatorParams(enabled=True)` into `TransientSolver`.
  - This enables separator voltage BC in the steady solve but leaves `TransientSolver.separator.enabled == False`, so the transient separator concentration BC is not applied.
  - Impact is small for current results because separator concentration/ohmic losses are only a few mV, but sep-on validation numbers are not fully self-consistent.
- Minor reporting issue in `scripts/parallel_validation.py`:
  - Returned `I_1C` is recomputed after discharge from remaining capacity, so it can report final remaining-current basis, not the initial 1C basis used to run the case.

## B. Separator Model Integration

- Core integration is correct when separator parameters are passed through the constructor.
- Steady solver:
  - Computes separator state from `I_app`.
  - Applies cathode-side electrolyte Dirichlet potential `phi_e_sep < 0` during discharge.
  - Reports cell voltage as `phi_s(collector)` when separator is enabled, because Li/Li+ reference is outside the separator at 0 V.
- Transient solver:
  - Shares the separator object with `SteadyStateSolver` when constructed normally.
  - Applies separator cathode-side concentration as the electrolyte boundary concentration.
- Tests cover zero-current identity, separator-on voltage drop, monotonic separator concentration depletion, and disabled/default backward compatibility.
- Caveat: `parallel_validation.py` should construct `TransientSolver(..., separator=SeparatorParams(enabled=True))` for sep-on runs.

## C. Physical Reasonableness of 3C = 30.6 mAh/g

- Yes, it is physically reasonable for the current synthetic network geometry.
- It is not a quantitative Khan 1CAL reproduction.
- The early 3C cutoff is consistent with:
  - synthetic network geometry rather than XCT-resolved 1CAL topology,
  - lower/changed mass loading and active interfacial connectivity,
  - large local polarization at high current,
  - many active sites becoming transport-limited before full cathode utilization.
- Separator losses are too small to rescue the gap:
  - saved parallel results show only about `5.3 mV` ohmic and `2.0 mV` concentration drop at 3C.
  - sep-on 3C capacity shifts only from about `30.6` to `29.7 mAh/g`.
- Low-C time-limited capacities should not be interpreted as true cutoff capacities; they are run-duration artifacts when voltage never reaches 3.0 V.

## D. Must-Do Changes / Readiness

- No must-do core solver change remains for a publishable/demonstrable method-development scaffold.
- Before presenting separator-on validation as final, fix `scripts/parallel_validation.py` separator construction and `I_1C` reporting.
- Before claiming quantitative reproduction of Khan Figure 4, geometry/mass loading must be matched to the paper XCT network; current results should be framed as qualitative method validation.

## E. Single Most Impactful Improvement

- Replace or calibrate the synthetic cubic network with a geometry-matched network:
  - match Khan 1CAL mass loading,
  - match electrode thickness/projected area,
  - match active interfacial area and phase connectivity,
  - preferably use the XCT-derived topology or a calibrated surrogate.
- This would dominate over further separator tuning, because the separator contributes only mV-scale changes while geometry controls the large capacity and high-rate polarization mismatch.
