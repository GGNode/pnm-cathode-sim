# TASK: Final Gap Analysis & Implementation Plan

## Context

We have a working PNM cathode discharge simulation (Khan et al. 2021 paper reproduction). Current state:
- 81 tests pass
- Core physics is correct: BV kinetics, Newton-Raphson coupling, D_e/κ/D_s concentration-dependent, OCV (Khan polynomial)
- Recent fixes: paper-equivalent current density, initial SoL=0.50, throat_scale removed, V-Q time alignment
- Remaining gap: separator model missing, network geometry doesn't match paper

## Current Results (10×10×10, throat_scale=1.0, k0=1e-10, SoL=0.50)

### With paper-equivalent current (I_1C = 53 A/m²)
| C-rate | V_init | V_final | phi_e span | phi_s span | Capacity | Cutoff? |
|--------|--------|---------|-----------|-----------|----------|---------|
| 0.2C | 3.896V | 3.355V | 8 mV | 158 mV | 222.5 | ❌ (time) |
| 3C | 3.547V | ~3.0V | 167 mV | 1377 mV | ~18.5 | ✅ |

### Paper reference (Khan 1CAL)
| C-rate | V_init | Capacity | Cutoff? |
|--------|--------|----------|---------|
| 0.2C | ~4.0V | ~195 | ✅ |
| 3C | ~3.85V | ~150 | ✅ |

### Key issue
Synthetic network mass loading = 127.8 g/m² vs paper's 297.8 g/m² (2.3× gap). Using paper current density on thinner electrode → overpolarization. Using network's own C-rate → underpolarization (previous behavior).

## My Assessment of Remaining Gaps

### MUST FIX (hard science without it)
1. **Separator model**: Only missing piece that changes physics conclusions. At 3C, separator contributes 200-300mV. Without it, high-rate V-Q shape is wrong. Needs at minimum: 1D separator Ohmic drop + concentration gradient.

### SHOULD FIX (engineering, not physics)  
2. **C-rate conversion approach**: Need to decide: use paper current (overpolarized on thin network) or network's own current (underpolarized without throat_scale)? Recommend: use network's own mass loading for C-rate, since we can't make the synthetic network match paper's mass loading without changing grid size/spacing.
3. **Network sizing**: 10×10×10 @ 10μm = 100μm vs paper 129μm. Could go to 13×13×13 but slow. Recommend: keep 10×10×10, accept the thickness difference.

### NICE TO HAVE (don't block on these)
4. pore-throat-pore series conductance
5. Migration/concentration-potential terms in electrolyte charge equation  
6. Better pore radius distribution matching XCT statistics

## Questions for Codex

1. **Separator model design**: Design a minimal but physically correct 1D separator model that can be added to the existing solver. The paper uses Eqs 2.21-2.24 (separator diffusion + potential + Li foil BV). Propose the simplest implementation that captures the dominant effect (Ohmic drop + concentration gradient).

2. **C-rate strategy**: Verify my assessment. Should we use network's own mass loading for C-rate conversion (self-consistent but thin electrode) or paper's mass loading (correct current but overpolarized)? What would produce more publishable/meaningful results?

3. **Is the current physics "good enough" for a paper-level reproduction?** Be honest. If someone reviewed this, what would they flag as a deal-breaker vs an acceptable simplification?

4. **Concrete execution plan**: For each MUST FIX and SHOULD FIX item, specify:
   - Exact files to modify
   - What to add/change (pseudocode level)
   - Expected impact on results
   - Estimated implementation effort

5. **Final validation strategy**: What test/scenario should we run at the end to demonstrate the implementation is correct?

Write your analysis to `PLAN_FINAL.md` in the project root.
