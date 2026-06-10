# Task: Comprehensive Review of PNM Cathode Simulation

## Context
We are reproducing Khan et al. (2021) "Pore Network Modelling of Lithium-Ion Cells" (J. Electrochem. Soc. 168, 070534) using OpenPNM. The model simulates NMC532 cathode discharge in a 2032 coin cell.

We just completed Phase 5A (parameter calibration). Need your expert review of the current state, root cause analysis of remaining issues, and revised implementation plan.

## Current Validation Results (5×5×5, seed=42, ε=0.368, T=303K)

| C-rate | Our V_init | Paper V_init | Our Capacity | Paper Capacity | Status |
|--------|-----------|-------------|-------------|----------------|--------|
| 0.2C   | 4.30V     | ~4.05V      | 200 mAh/g   | ~195           | ✅ close but V_init wrong |
| 0.5C   | 4.30V     | ~4.05V      | 205 mAh/g   | ~185           | ❌ 10% high |
| 1C     | 4.30V     | ~4.00V      | 214 mAh/g   | ~175           | ❌ 22% high |
| 3C     | 4.30V     | ~3.85V      | 9.5 mAh/g   | ~150           | ❌❌ catastrophically low |

## Critical Issues to Analyze

### Issue 1: V_init = 4.30V but Paper shows 4.05V at 0.2C

Our OCV function at x=0.35 (initial SoL) gives 4.30V. But the paper Figure 4 shows all curves starting at ~4.05V. This 0.25V offset is wrong.

**Khan polynomial (Eq. 2.19) at x=0.35:**
```python
nmc532_ocv(0.35) → 4.30V
```

**Paper Figure 4 initial voltage:** ~4.05V at 0.2C, ~3.85V at 3C (before IR drop)

Possible causes:
1. Khan polynomial is defined for SOC ∈ [0,1] not x ∈ [0,1]. SOC=0 might correspond to x=0.25 (fully discharged), SOC=1 to x=1.0 (fully charged). Then initial state x=0.35 → SOC = (0.35-0.25)/(1.0-0.25) = 0.133.
2. The polynomial was fit to Figure 4 data but with different x-axis mapping.
3. The polynomial might use a different convention where SOC is normalized differently.

**Key data:** If we compute U(0.35) = 4.30V but the paper says 4.05V, there's a mapping error between x (stoichiometry) and the polynomial's input variable.

### Issue 2: 3C gives only 9.5 mAh/g (paper: ~150)

After Phase 5A fixes, 3C hits cutoff at step 40 with only 9.5 mAh/g. The voltage drops from 4.30V to 3.0V almost immediately.

Possible causes:
1. **V_init wrong → fast cutoff**: If initial voltage is 4.30V (should be ~3.85V), the IR drop at 3C only needs to drop 1.3V instead of 0.85V to reach cutoff. With η_act ≈ 0.3-0.5V + η_ohmic ≈ 0.2-0.4V, the total overpotential could easily exceed 1.3V.
2. **Overpotentials too large**: k0=1e-10 (correct value) gives large activation overpotentials. At 3C this dominates.
3. **OCV function is flat/stiff at x=0.35**: dU/dx might be very small, meaning even a slight SoL change gives huge voltage change.

### Issue 3: 0.2C and 0.5C capacities too close (200 vs 205)

In the paper, 0.2C gives ~195 and 1C gives ~175 — there should be clear differentiation between C-rates. Our 0.2C/0.5C/1C give 200/205/214 which shows almost no C-rate dependence.

Possible causes:
1. The 5×5×5 network is too coarse — transport limitations (which cause capacity loss at high C) are insufficient.
2. The discharge might be stopping too early (voltage drops below cutoff due to wrong OCV, not due to actual capacity depletion).

### Issue 4: Spatial heterogeneity in 5×5×5

At 3C snapshot:
- SoL range: [0.35, 0.77] (separator to collector varies by 0.13)
- phi_s gradient: 0.39V across electrode
- c_e almost uniform (delta ~ 1e-12)

The c_e being perfectly uniform suggests electrolyte transport has NO effect. This is either:
1. Correct for 5×5×5 (too small for meaningful gradients)
2. De/D_eff too large

## Code Parameters After Phase 5A

| Parameter | Code Value | Paper Value | Status |
|-----------|-----------|-------------|--------|
| k0 | 1e-10 m/s | Table II: 1e-10 | ✅ correct |
| ε (porosity) | 0.368 | Table I: 0.368 | ✅ correct |
| T | 303K | 30°C (303K) | ✅ correct |
| OCV | Khan Eq.2.19 polynomial | Figure 4: V_init≈4.05 | ❌ wrong mapping |
| D_e | 1.37e-10 m²/s | Eq.2.21: concentration-dependent | ❌ constant, not concentration-dependent |
| κ | 0.655 S/m | Eq.2.22: concentration-dependent | ❌ constant, not concentration-dependent |
| D_s | 1e-14 m²/s | Table II: 1e-14 | ✅ correct |
| R_p | 5.22e-6 m | Table I: 5.22μm | ✅ correct |
| cs_max | 48900 mol/m³ | Table I: 48900 | ✅ correct |
| L_cathode | ~100μm | Table I: 105μm | ⚠️ network-dependent |

## OCV Polynomial Details

```python
# Khan et al. Eq. 2.19 (corrected soc^8 coefficient from -5520 to -35520)
# U = sum(coeffs[i] * soc^(9-i)) - 0.0003*exp(7.657*soc^115)
# Coefficients:
[5744.862289, -35520.41099, 95714.29862, -147364.5514, 142718.3782,
 -90095.81521, 37061.41195, -9578.599274, 1409.309503, -85.31153081]

# Evaluations:
U(0.35) = 4.30V  ← code uses this as initial
U(0.50) = 4.00V
U(0.90) = 3.67V
U(1.00) = 2.94V → clip to 3.0V
```

**The question**: What does the paper mean by "SOC" in the polynomial? Is SOC the same as x (Li stoichiometric fraction), or is SOC = (x - x_min)/(x_max - x_min)?

If SOC = (x - 0.25)/(1.0 - 0.25), then initial x=0.35 → SOC = 0.133, and U(0.133) might give ~4.05V. Check this.

## Questions for Codex

1. **Root cause of V_init mismatch**: Is the polynomial input supposed to be SOC (0-1 normalized) rather than x (Li stoichiometry)? What's the correct mapping?

2. **3C catastrophic failure**: Is this purely due to V_init being 0.25V too high (so IR drop easily exceeds the 1.3V headroom), or are the overpotentials genuinely too large?

3. **Missing concentration-dependent transport**: How important is De(c_e) and κ(c_e) for matching Figure 4? Would constant De/κ with correct V_init be sufficient?

4. **Network size**: Is 5×5×5 sufficient for Figure 4 reproduction, or do we need 20×20×20?

5. **Revised Phase plan**: Given the above, what's the correct priority order for Phase 5B/5C/5D? Please provide a concrete step-by-step plan with expected impact on each metric.

## Key Source Files
- `src/physics/ocv.py` — OCV polynomial (72 lines)
- `src/solver/transient.py` — transient solver (868 lines)
- `src/solver/steady.py` — Newton-Raphson steady solver (762 lines)
- `src/solver/single_pore.py` — single-pore BV + solid diffusion (212 lines)
- `src/physics/solid.py` — NMC diffusion coefficient (290 lines)
- `src/network/generator.py` — synthetic network creation
- `scripts/validate_paper_figures.py` — validation script

## Reference Documents
- `PAPER_REFERENCE.md` — complete paper data (tables, equations, figure values)
- `AUDIT_R5.md` — detailed parameter audit with file:line references
