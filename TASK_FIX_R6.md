# TASK: Phase 6 Fixes — Paper Comparison Alignment

Based on AUDIT_R6.md findings, implement the following fixes in priority order.

## Fix 1: Paper-Equivalent Current Density (P0)

**Problem**: The validation uses `current_density_for_c_rate()` which derives current from the synthetic network's NMC volume. This gives 1C = 22.9 A/m² instead of paper's ~53.0 A/m².

**Fix in `scripts/validate_paper_figures.py`**:
- Override the C-rate-to-current conversion with paper areal loading:
  ```python
  PAPER_MASS_LOADING = 297.8e-3  # kg/m² (297.8 g/m²)
  SPECIFIC_CAPACITY = 178.0 * 3600  # Ah/kg → C/kg (178 mAh/g)
  I_1C_PAPER = PAPER_MASS_LOADING * SPECIFIC_CAPACITY / 3600  # A/m² ≈ 53.0
  ```
- In the discharge loop, compute `I_app = -C_rate * I_1C_PAPER` directly instead of calling `solver.current_density_for_c_rate(C_rate)`.
- Keep the solver's internal `current_density_for_c_rate()` unchanged (it's correct for its own capacity basis).

## Fix 2: Initial Stoichiometry (P0)

**Problem**: `INITIAL_SOL = 0.35` gives `U(0.35) = 4.30V`. Paper starts at ~4.0V which corresponds to `x ≈ 0.50`.

**Fix in `scripts/validate_paper_figures.py`**:
- Change `INITIAL_SOL = 0.35` → `INITIAL_SOL = 0.50`
- This gives `U(0.50) ≈ 3.996V` matching paper's starting voltage.

## Fix 3: Remove throat_scale Hack (P0)

**Problem**: `throat_scale=2.0` multiplies diameter → area becomes 4×, not 2×. This quadruples conductance AND reaction interface area simultaneously.

**Fix in `scripts/validate_paper_figures.py`**:
- Change `THROAT_SCALE = 2.0` → `THROAT_SCALE = 1.0`

**Also fix in `src/network/generator.py`**:
- Verify that throat_scale=1.0 means no modification (pass-through).

## Fix 4: V-Q Time Alignment (P0)

**Problem**: `TransientSolver.step()` returns voltage computed BEFORE the concentration update, but the capacity is recorded AFTER. This makes V-Q curves flatter than reality.

**Fix in `src/solver/transient.py`**:
- In the `step()` method (around line 564-670), the voltage should be reported at the same time as the capacity. Two options:
  - Option A: Record capacity at start-of-step (before update) to match start-of-step voltage.
  - Option B: Recompute potential after concentration update (expensive but more accurate).
- **Preferred**: Option A — store `t_n` capacity before the update, and report that with the `t_n` voltage. Then after the update, the next step will report `t_{n+1}` voltage with `t_{n+1}` capacity.
- In `run_discharge()` (around line 879-882), the `q_delivered` should be `|I| * t_n` (step start time), not `|I| * t_new` (step end time). Or equivalently, record voltage at step end after re-solving.

**Minimal fix**: After the concentration update at line 659-663, add a second steady-state solve with the new concentrations and return THAT voltage. This is the physically correct approach — the voltage at the end of a time step should reflect the end-of-step state.

```python
# After c_e/c_s update (line ~663), re-solve potentials:
pot_new = self._steady.solve(I_app=I_app, c_e_new=self.c_e, c_s_new=self.c_s)
# Return pot_new["voltage"] instead of pot["voltage"]
```

If that's too expensive, at minimum fix the capacity recording to match the voltage time.

## Fix 5: Correct Initial Voltage in Report

**Problem**: The validate script's `get_initial_voltage()` computes OCV at INITIAL_SOL. With the new INITIAL_SOL=0.50, this should give ~4.0V.

**Fix**: Just verify after Fix 2 that V_init reports correctly.

## Verification

After all fixes:
1. `python -m pytest tests/ -q` — all tests pass
2. Quick smoke test: run 0.2C and 3C with the new parameters, verify:
   - V_init ≈ 4.0V (not 4.3V)
   - 3C current ≈ 159 A/m² (not 68.7)
   - Voltage actually drops toward 3.0V
   - 0.2C capacity is in the range 180-210 mAh/g
   - 3C capacity < 0.2C capacity
3. If validate_paper_figures.py runs within timeout, generate new figures

## Constraints
- Do NOT modify the physics equations (electrolyte.py, solid.py, reaction.py, ocv.py) — those are correct.
- Do NOT modify the steady-state solver's Newton-Raphson — it's verified correct.
- Focus only on the validation script parameters and the V-Q alignment bug.
- Keep all existing tests passing.
