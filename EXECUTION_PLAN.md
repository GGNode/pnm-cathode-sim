# PNM Cathode Simulation — Revised Execution Plan

> **Last updated:** Round 5 (Codex GPT-5.5 review)
> **Codex verdict:** OCV mapping was correct. Core blockers are network pathology, current scaling, and validation methodology.

---

## Critical Finding: V_init = 4.30V IS CORRECT

Codex confirmed: the paper Figure 4 1CAL 0.2C PNM curve starts at ~4.3V. Our `U(0.35) = 4.30V` matches the paper. The "4.05V" from vision analysis was the under-load voltage (not OCV). **Do NOT remap OCV input.**

---

## Current Status Summary (5×5×5, post-Phase-5A)

| Metric | 0.2C | 0.5C | 1C | 3C | Paper 1CAL |
|--------|------|------|-----|-----|-----------|
| V_init (unloaded) | 4.30V | 4.30V | 4.30V | 4.30V | ~4.30V ✓ |
| V_init (loaded) | 4.13V | 4.02V | 3.90V | 3.52V | ~4.05V |
| Capacity | 200* | 205* | 214* | 9.5† | 195/185/175/150 |
| Reached cutoff? | ❌ | ❌ | ❌ | ✅ | all reach cutoff |
| Reactive interfaces | 15 | 15 | 15 | 15 | ~1000+ |

\* = time-limited (not cutoff), not comparable to paper capacity
† = network pathology (5/43 active electrolyte nodes, 1.03V phi_s drop)

---

## Root Cause Analysis

### Problem 1: Network is too small and disconnected (BLOCKING)
- 5×5×5 = 125 nodes → only **5/43 electrolyte nodes active**, **15 reactive interfaces**
- phi_e is **exactly flat** (no electrolyte transport dynamics at all)
- phi_s spans **1.03V** at 3C (paper: a few mV) → excessive solid ohmic drop
- **No percolation check** — synthetic random network may have disconnected paths

### Problem 2: Current density scaling is wrong (BLOCKING)
- Synthetic mass loading: **640 g/m²** vs paper 1CAL: **297.8 g/m²**
- Our 1C current density: **115 A/m²** vs paper: **53 A/m²**
- Current is ~2× too high → overpotentials scale accordingly
- Cause: `collector_area` uses throat-sum area, not geometric projected area

### Problem 3: Validation methodology is flawed
- Script records **unloaded OCV** as V_init (I_app=0) — should report loaded V_init
- Low C-rate runs **never reach cutoff** (time-limited) — capacities are NOT comparable to paper
- No diagnostics: active node fraction, reactive interface count, phi ranges, eta ranges

### Problem 4: Transport parameters are constant (secondary)
- D_e = constant 1.37e-10 m²/s (paper: concentration-dependent, Eq. 2.21)
- κ = Nernst-Einstein with constant D_e (paper: polynomial, Eq. 2.22)
- BUT: with flat phi_e, transport dynamics have no effect regardless
- **Priority: fix network first, then transport matters**

---

## Revised Phase Plan

### Phase 5B: Validation Cleanup + Diagnostics (1-2h) — NEXT

**Goal:** Fix validation methodology so we can trust the numbers.

1. **Add network diagnostics** to validate script:
   - `active_e_fraction` = active electrolyte nodes / total electrolyte nodes
   - `reactive_interfaces` count
   - `phi_e_range`, `phi_s_range` (mV)
   - `eta_range` (mV)
   - Mass loading (g/m²) and current density (A/m²)
   - Print "⚠️ LOW CONNECTIVITY" warning if active_e < 50%

2. **Report loaded V_init** (not unloaded OCV):
   - Record first step voltage (after I_app applied), not I_app=0 voltage

3. **Fix capacity comparison**:
   - Mark "time-limited" vs "cutoff-limited" explicitly
   - Only compare against paper when reached_cutoff=True

4. **Geometric collector area**:
   - Use projected electrode area (A_cell = π × r²) for current density
   - Or match paper's areal current density directly

### Phase 5C: Network Scaling (2-4h)

**Goal:** Get a percolating network with enough nodes for transport dynamics.

1. **Scale to 20×20×20** (8000 nodes):
   - Target: hundreds of reactive interfaces, >50% active electrolyte
   - Newton solver already converges at this size (verified Round 4)

2. **Percolation validation**:
   - Check both electrolyte and solid phases percolate from separator to collector
   - Remove disconnected clusters
   - Report connectivity diagnostics

3. **Match paper's network statistics**:
   - Paper 1CAL: 1032 pores, 2936 throats
   - Paper 3CAL: 4637 pores, 14808 throats
   - Our 20³: ~8000 pores, ~24000 throats (comparable to 3CAL)

4. **Target metrics**:
   - phi_s range < 50mV at 0.2C (paper: a few mV)
   - >100 reactive interfaces (paper: ~1000+)
   - 0.2C capacity: cutoff-limited, ~195 mAh/g

### Phase 5D: Current + Transport Scaling (2-4h)

**Goal:** Match paper's current density and electrolyte dynamics.

1. **Fix current density**:
   - Match paper's areal current: 1C ≈ 53 A/m² (1CAL) or adjust by mass loading ratio
   - Rerun all C-rates with correct current

2. **Wire concentration-dependent transport**:
   - `D_e(c_e)` from paper Eq. 2.21 (verify units: cm²/s → m²/s)
   - `κ(c_e)` from paper Eq. 2.22 (mS/cm → S/m)
   - Rebuild electrolyte conductances each time step

3. **Separator coupling**:
   - Replace fixed cathode-side reservoir with proper separator boundary condition
   - Implement ionic current flux at separator interface

### Phase 5E: Figure Comparison (1-2h)

**Goal:** Generate publication-quality comparison plots.

1. Figure 4: Voltage vs Capacity at 0.2C/0.5C/1C/3C overlay on paper data
2. Figure 6: SoL spatial distribution at specific capacity points
3. Figure 7: phi_e spatial distribution at high C-rate

---

## Decision Log

| Decision | Rationale |
|----------|-----------|
| **Do NOT remap OCV** | Codex confirmed V_init=4.30V matches paper's PNM curve |
| **Network before transport** | phi_e flat → transport dynamics meaningless at 5×5×5 |
| **Diagnostics before fixes** | Can't fix what you can't measure |
| **Geometric area, not throat-sum** | Match paper's reported mass loading and current density |
| **20×20×20 as target** | Already verified convergence; comparable to paper's 3CAL |

---

## Risk Register

| Risk | Impact | Mitigation |
|------|--------|-----------|
| 20×20×20 memory | 324MB for pore arrays, ~1GB total | Feasible on 16GB system |
| 20×20×20 solve time | Newton: ~5s per step (Round 4) | Acceptable for full discharge |
| Percolation failure | Network may not percolate | Use percolation-aware generator |
| D_e unit error | Codex flagged: electrolyte.py may have cm²/s vs m²/s issue | Verify before wiring in |
