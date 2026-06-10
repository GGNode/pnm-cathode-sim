# Code Review: Newton-Raphson Steady-State Solver — Current Conservation Bug

## Context

This is a pore network model (PNM) for Li-NMC532 cathode. The steady-state solver was recently refactored from SOR to coupled Newton-Raphson.

## The Bug

After Newton convergence, current conservation is violated:

```
sum(I_rxn * A_interface) / (I_app * A_collector) = 0.30  (should be 1.0)
```

The solver converges (residual < tol) but the converged solution doesn't satisfy the physical constraint that all applied current must go through interfacial reactions.

## Key Files

- `src/solver/steady.py` — Newton-Raphson solver (just refactored)
- `src/physics/reaction.py` — Butler-Volmer kinetics
- `src/physics/ocv.py` — OCV lookup table (just refactored)
- `DERIVATION.md` — Full mathematical derivation with sign conventions

## What to Investigate

1. Read `DERIVATION.md` §3-§5 to understand the correct residual equations
2. Read the Newton solver in `steady.py` and check:
   - Is the collector BC (I_app * A_cc) correctly applied to the solid residual?
   - Are the BV conductance signs correct in the Jacobian blocks?
   - Does the Dirichlet BC for phi_e at separator properly allow current flow?
   - Is there a missing term coupling the reaction to the collector current?
3. Check if `_build_coupled_jacobian` or the residual assembly has a sign error
4. Verify that `sum(F_s) = 0` implies current conservation (sum of reaction = applied)

## Diagnostic Script

Run this after reading the code to reproduce the bug:

```python
import sys; sys.path.insert(0, 'src')
from network.generator import create_cathode_network
from solver.transient import TransientSolver
import numpy as np

net = create_cathode_network(shape=[5,5,5], spacing=1e-5, porosity=0.5, cbd_fraction=0.1, seed=42)
solver = TransientSolver(net, T=298.15, k0=5e-10)
solver.set_concentration(c_e=1200.0, c_s=0.35*48900.0)
I_app = solver.current_density_for_c_rate(0.2)

pot = solver._steady.solve(I_app=I_app)
interfaces = solver._steady.reactive_interfaces
I_rxn = pot['I_rxn']

total_rxn = sum(I_rxn[k] * A for k, (_, _, A) in enumerate(interfaces))
cc_current = I_app * solver._steady._cc_area
print(f'Applied: {cc_current:.6e} A, Reaction: {total_rxn:.6e} A, Ratio: {total_rxn/cc_current:.4f}')
```

Expected: ratio ≈ 1.0. Actual: ratio ≈ 0.30.

## Deliverable

Provide a detailed analysis of the root cause and a concrete fix (code diff or description of what to change).
