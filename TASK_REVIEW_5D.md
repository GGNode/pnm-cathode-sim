# Codex Review: Current State & Next Steps

## Context
We are reproducing Khan et al. (2021) PNM cathode simulation. After Phase 5A-C, we have:
- 10×10×10 synthetic network with percolation validation
- Newton-Raphson steady solver + transient discharge
- Loaded V_init reporting, diagnostics, geometric area

Please review the current code and validation results, then provide a prioritized plan for the next phase.

## Current Results (10×10×10, ε=0.368, T=303K, k0=1e-10)

| C-rate | V_init (loaded) | V_final | Capacity | Status | phi_e span | phi_s span |
|--------|----------------|---------|----------|--------|------------|------------|
| 0.2C   | 4.126V         | 3.29V   | 200 mAh/g* | time-limited | 14 mV | 154 mV |
| 0.5C   | 4.023V         | —       | —        | time-limited | 39 mV | 352 mV |
| 1C     | 3.901V         | —       | —        | time-limited | 82 mV | 641 mV |
| 3C     | 3.860V         | 3.00V   | 78 mAh/g | cutoff | 237 mV | **1668 mV** |

*0.2C is time-limited (never reaches cutoff), so 200 mAh/g is not the true capacity.

**Paper reference (1CAL XCT network):**
- 0.2C: ~195 mAh/g, 1C: ~175 mAh/g, 3C: ~150 mAh/g
- All reach cutoff at 3.0V
- phi_s varies by only a few mV at low C

## Network Stats
- 1000 pores, 2700 throats
- NMC: 531, Electrolyte: 371, Solid: 629
- active_e: 279/371 (75%), reactive_interfaces: 750
- Both phases percolate ✓
- mass_loading: 167.7 g/m² (paper 1CAL: 297.8 g/m²)
- geometric_area: 1e-8 m²

## Key Issues to Review

### Issue 1: phi_s spans 1.67V at 3C (paper: few mV)
This is the biggest red flag. Current solid conductivities:
```python
sigma_nmc = 0.01 S/m   # line 170, steady.py
sigma_cbd = 760.0 S/m  # line 170, steady.py
```

Questions:
- Is sigma_nmc=0.01 S/m correct for NMC532? Literature values range from 0.01 to 10 S/m.
- The effective composite conductivity depends on the volume fraction and connectivity. With 53% NMC + 10% CBD, what's the expected effective conductivity?
- Should we use a Bruggeman or percolation-based effective conductivity instead of raw material values?

### Issue 2: Electrolyte transport is constant, not concentration-dependent
Currently (steady.py line 164-167):
```python
D_e = 2.0e-10  # constant
kappa = 1200.0 * F**2 * D_e / (R * self.T)  # Nernst-Einstein, constant
```

The paper has:
- D_e(c_e) = 10^(-4.43 - 54/(T-229-5c2) - 0.22c2) [Eq. 2.21]
- κ(c_e) polynomial [Eq. 2.22]

The `electrolyte.py` module already implements these but they're not wired in.

Questions:
- Should we wire in electrolyte.py before or after fixing solid conductivity?
- The electrolyte.py `electrolyte_diffusion_coefficient` returns values — verify the units are correct (m²/s, not cm²/s).

### Issue 3: Low-rate runs never reach cutoff
0.2C/0.5C/1C are time-limited. The discharge runs for max_time = 1.25 * 3600/C_rate seconds but voltage never drops to 3.0V.

Possible causes:
- The time limit is too short
- The network is too small (no transport limitations at low C)
- The cutoff should be reached but the solver stops too early

### Issue 4: Current density at 3C
3C gives I_app = -90.22 A/m². Paper 1CAL at 3C would be ~160 A/m² (3 × 53 A/m²). Our current is lower because mass loading is lower (167.7 vs 297.8 g/m²).

Is this a fair comparison? Should we match the paper's areal current density or use our network's actual mass loading?

### Issue 5: 0.2C capacity of 200 mAh/g is time-limited
This means the discharge didn't complete. The true capacity would be higher. But the paper's 0.2C is ~195 mAh/g and reaches cutoff. Why doesn't ours reach cutoff?

Possible: our phi_s span of 154 mV at 0.2C means some nodes have phi_s = 4.15V (near OCV) while others have phi_s = 4.30V. The voltage reported is the average at the collector, which stays above 3.0V because the network is small and well-connected at low C.

## Code to Review

### src/solver/steady.py — _compute_conductances()
Lines 142-190. Review the conductance model:
- Is the Nernst-Einstein kappa correct?
- Is the harmonic mean for mixed NMC-CBD throats correct?
- Should conductances be recomputed each time step (they're currently computed once in __init__)?

### src/solver/transient.py — step() method
Review the transient coupling:
- Are the source terms for c_e and c_s correct?
- Is the backward Euler implementation stable?
- Are the concentration bounds (clip) correct?

### src/physics/electrolyte.py
Review the transport correlations:
- Are the units correct (m²/s, S/m)?
- Is the T-229 denominator safe (T must be > 229+5c2)?

### src/network/generator.py — create_cathode_network()
Review the network generation:
- Is the random phase assignment correct for the paper's electrode composition?
- Are pore/throat sizes reasonable?
- Is the spacing=1e-5 m appropriate?

## Questions for Codex

1. **What is the single highest-impact fix?** Is it solid conductivity, electrolyte transport, or network scaling?

2. **Is sigma_nmc=0.01 S/m physically correct?** If not, what should it be? Should we use a composite effective conductivity?

3. **Should we wire in electrolyte.py now?** Or fix solid conductivity first?

4. **Is the 0.2C time-limited behavior expected?** Or does it indicate a bug?

5. **What's the correct comparison methodology?** Should we match paper's areal current density or use our network's mass loading?

6. **Please provide a concrete Phase 5D plan** with specific code changes, expected impact, and verification steps.

## Reference
- `PAPER_REFERENCE.md` — paper data (Table I/II, equations)
- `AUDIT_R5.md` — parameter audit
- `EXECUTION_PLAN.md` — current plan (needs updating)
