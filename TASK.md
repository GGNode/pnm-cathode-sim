# PNM Cathode Solver Fix — Full Audit Remediation

## Context

This is a pore network model (PNM) for Li metal | separator | NMC532 cathode half-cell.
The implementation has passed unit tests (65/66) but **fails validation against the reference paper**:

| C-rate | Paper capacity (mAh/g) | Current simulation | Issue |
|--------|----------------------|-------------------|-------|
| 0.2C | ~178 | 62.2 | Too low |
| 0.5C | ~170 | 9.1 | Way too low, should be close to 0.2C |
| 1C | ~155 | 10.5 | Capacity ordering inverted |
| 3C | ~120 | 13.0 | Capacity ordering inverted |

Additionally, spatial validation shows **voltage divergence** (-17V at 1C, -53V at 3C).

The root causes are identified in AUDIT.md. Your job is to fix ALL of them.

## Reference Documents

- `DERIVATION.md` — Full mathematical derivation with sign conventions, BCs, and verification criteria
- `AUDIT.md` — Checklist of implementation vs derivation, with fix plan

## Phase 1: Fix Sign Conventions (CRITICAL)

### 1a. Transient reaction concentration signs

**File:** `src/solver/transient.py` ~lines 166-171

The derivation (§8.2, lines 497-500) requires:
- `dc_e/dt = +i_a/(F*V_e)` (electrolyte gains Li during discharge)
- `dc_s/dt = -i_a/(F*V_s)` (solid loses Li during discharge)

Current code has these **reversed**. Fix the signs AND the comments (which incorrectly say `I_rxn > 0` is discharge).

### 1b. Current sign documentation

**File:** `src/solver/steady.py` ~line 220

The derivation (§1.3) defines `I_app < 0` for discharge. The code comment at line 220 says positive `I_app` is discharge. Fix the comment to match the derivation.

**File:** `src/physics/reaction.py` ~line 35

Comment says positive overpotential is discharge; derivation says it's anodic (charge). Fix the comment.

## Phase 2: Fix Boundary Conditions (CRITICAL)

### 2a. Remove extra solid potential pin

**File:** `src/solver/steady.py` ~lines 236-239

Under galvanostatic operation, the derivation (§15.3, lines 916-923) says:
- Pin electrolyte potential at separator nodes to 0 (gauge)
- Do NOT pin solid potential at separator

Current code pins solid at separator-side nodes too. Remove this extra pin.

### 2b. Handle disconnected components

**File:** `src/solver/steady.py` ~lines 230-242

The derivation (§15.5, lines 933-938) requires:
- Keep only electrolyte components connected to the separator
- Keep only solid components connected to the collector
- Deactivate disconnected nodes (set reaction rate to zero, remove from linear system)

Currently the code only pins zero rows, which is insufficient. Implement proper connected-component detection and deactivation.

## Phase 3: Fix Capacity Calculation and C-rate Conversion

### 3a. C-rate to current conversion

**File:** `src/post/analysis.py` ~lines 82-90

The network discharge API accepts raw `I_app` but has no proper C-rate-to-current conversion from electrode capacity. The single-pore model does this correctly at `src/solver/single_pore.py:110-111`. Add equivalent functionality for the network model.

### 3b. Capacity calculation

Verify that `src/post/analysis.py` computes capacity correctly. The current result of 62 mAh/g at 0.2C vs paper's 178 mAh/g suggests either:
- The current applied is too high (C-rate conversion bug)
- The capacity integration is wrong
- The active material mass/volume accounting is wrong

Investigate and fix.

### 3c. Capacity ordering

Higher C-rate should give LOWER capacity (less time before cutoff). Currently 3C > 0.2C which is physically wrong. This is likely a downstream effect of fixing 3a/3b, but verify.

## Phase 4: Add Missing Tests

### 4a. Zero-current OCV test

Add a test that at `I_app = 0`:
- `I_rxn == 0` everywhere
- `eta == 0` everywhere  
- `V_cell == U(SoL)` (OCV)

### 4b. Small-current linear scaling

Add a test that varies current magnitude at low rates and verifies polarization scales linearly.

### 4c. Lithium conservation at nonzero current

Add a test that audits total lithium conservation at each time step during nonzero-current discharge, per DERIVATION.md §13 (lines 975-1006).

## Phase 5: Validation Rerun

After all fixes, rerun `scripts/validate_against_paper.py` and report:
1. V-Q curves for 0.2C, 0.5C, 1C, 3C vs paper Figure 4
2. Spatial distributions at 3C vs paper Figure 6
3. Capacity values vs paper Table III

## Constraints

- Read DERIVATION.md thoroughly before making any changes — it defines ALL sign conventions
- Keep all existing tests passing (run `pytest` after each phase)
- Use SI units throughout
- Do not change the network generator or OCV curve — those are correct
- Commit after each phase with descriptive message
- If a fix is uncertain, add a `# TODO(pnm-fix):` comment explaining the uncertainty

## Verification

After each phase:
```bash
cd /VOLUMES/1TB/projects/pnm-lib-cathode
.venv/bin/pytest -q
```

After all phases:
```bash
.venv/bin/python scripts/validate_against_paper.py
```
