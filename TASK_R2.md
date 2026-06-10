# PNM Fix Round 2: Remove Overly Aggressive Flux Limiter

## What happened in Round 1

Codex fixed:
- Transient reaction signs (correct)
- Separator reservoir BC (correct)  
- Steady solver warm-start (correct)
- C-rate API in analysis.py (correct)
- New tests (correct)

But also added a **flux limiter** in transient.py (around line 570-620) that is WRONG.

## The Problem

The flux limiter caps each interface flux so a single step cannot overfill/deplete a control volume. But it uses the CURRENT concentration as reference, not the full range:

```python
max_rate = max(48900.0 - cs, 0.0) * vol / dt
```

At initial SoL=0.5 (cs=24450), this allows only (48900-24450)*vol/dt flux — about half the theoretical rate. At 0.2C with dt=20s, this severely limits the reaction. At 1C with dt=10s, the cap is more generous (same flux limit but smaller dt means more steps fit in the same time).

Result: **capacity ordering is inverted** (0.2C=69.6 mAh/g, 1C=122.6 mAh/g). Physics says 0.2C should give ~178 mAh/g.

## What to Do

### Step 1: Remove the flux limiter entirely

In `src/solver/transient.py`, find the flux limiter block that Codex added (the block with `interface_fluxes: list`, `flux_arr`, `le_arr`, `ls_arr`, `lithiation_rate`, `delithiation_rate`, `max_rate`, etc.). Remove it entirely. 

Replace with just the original per-interface loop:

```python
for t, le, ls, A_intf, e_g, nmc_g in self._interfaces:
    # ... existing BV calculation ...
    flux = I_rxn * A_intf / F  # mol/s
    reaction_mol_rate += flux
    dce_dt[le] += flux / self._vol_e[le] if self._vol_e[le] > 0 else 0.0
    dcs_dt[ls] -= flux / self._vol_s[ls] if self._vol_s[ls] > 0 else 0.0
```

The transient solver ALREADY clamps concentrations after the implicit diffusion solve (line ~658):
```python
self.c_s = np.clip(self.c_s, 1.0, cs_max - 1.0)
```

This is sufficient. No flux limiter needed.

### Step 2: Keep everything else from Round 1

Do NOT touch:
- The separator reservoir BC (`c_e_reservoir` pinning at sep_e nodes)
- The steady solver warm-start (`_phi_e_guess`, `_phi_s_guess`)
- The C-rate API in analysis.py
- The new tests in test_steady.py and test_post.py
- The lithium conservation test improvements

### Step 3: Run tests

```bash
cd /VOLUMES/1TB/projects/pnm-lib-cathode
.venv/bin/pytest -q
```

### Step 4: Run validation

```bash
.venv/bin/python scripts/validate_paper_figures.py
```

Check that:
1. Capacity ordering: 0.2C > 0.5C > 1C > 3C (correct physical ordering)
2. 0.2C capacity is close to 178 mAh/g (paper value)
3. No negative voltages
4. All C-rates reach cutoff

Report the validation metrics after fixing.
