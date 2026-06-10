# PNM Fix Round 3: Newton-Raphson Solver for Steady State

## Problem

The current steady-state solver uses SOR (Successive Over-Relaxation) with alpha=0.5, solving phi_e and phi_s decoupled. This diverges on networks larger than 5×5×5=125 nodes. The 10×10×10 network diverges at step 7 (V jumps to -85V).

The root cause: decoupled SOR cannot handle the strong phi_e-phi_s coupling through Butler-Volmer kinetics. When the BV reaction is fast (high i0), small changes in eta = phi_s - phi_e - U cause large reaction currents, which the decoupled solver can't track.

## Goal

Replace the SOR solver in `src/solver/steady.py` with a coupled Newton-Raphson solver. Keep the same class interface (`solve(I_app, ...)`) and output format.

## Architecture

The unknown vector is `x = [phi_e; phi_s]` of length `n_e + n_s`.

The residual vector is:
```
F_e[i] = Σ_k G_ik^e (phi_k^e - phi_i^e) + Σ_r i_rxn * A_r     (electrolyte charge conservation)
F_s[m] = Σ_n G_mn^s (phi_n^s - phi_m^s) - Σ_r i_rxn * A_r + BC  (solid charge conservation)
```

where `i_rxn` is the Butler-Volmer current at each reactive interface.

The Jacobian J = dF/dx is a 2×2 block matrix:
```
J = [ J_ee  J_es ]
    [ J_se  J_ss ]
```

where:
- J_ee = dF_e/dphi_e = L_e - diag(g_bv_e)  (electrolyte Laplacian minus BV conductance)
- J_es = dF_e/dphi_s = +g_bv at interface entries
- J_se = dF_s/dphi_e = -g_bv at interface entries  
- J_ss = dF_s/dphi_s = L_s - diag(g_bv_s)

g_bv = d(i_rxn * A)/d(eta) is the linearized BV conductance (already computed in the current code).

## Implementation Steps

### Step 1: Build the coupled Jacobian

In `src/solver/steady.py`, inside the `solve()` method, replace the SOR loop with:

```python
def _build_coupled_jacobian(self, phi_e, phi_s, I_app):
    """Build the full (n_e+n_s) x (n_e+n_s) Jacobian and residual."""
    n_e = len(phi_e)
    n_s = len(phi_s)
    N = n_e + n_s
    
    # Compute BV quantities at interfaces
    # (same as current code: i0, eta, I_rxn, g_bv)
    
    # Build residual F
    F = np.zeros(N)
    # F[0:n_e] = L_e @ phi_e + reaction_terms_e + BC_e
    # F[n_e:N] = L_s @ phi_s + reaction_terms_s + BC_s
    
    # Build Jacobian J (sparse)
    # J[0:n_e, 0:n_e] = L_e - diag(g_bv_e)  (J_ee)
    # J[0:n_e, n_e:N] = +g_bv entries        (J_es)
    # J[n_e:N, 0:n_e] = -g_bv entries        (J_se)
    # J[n_e:N, n_e:N] = L_s - diag(g_bv_s)  (J_ss)
    
    # Apply Dirichlet BCs: replace rows for pinned nodes
    # sep_e nodes: F[i] = phi_e[i] - 0, J[i,:] = unit row
    # inactive nodes: similar
    
    return J, F
```

### Step 2: Newton iteration

```python
# Newton-Raphson loop
x = np.concatenate([phi_e, phi_s])
for iteration in range(max_iter):
    phi_e, phi_s = x[:n_e], x[n_e:]
    J, F = self._build_coupled_jacobian(phi_e, phi_s, I_app)
    
    # Solve J @ dx = -F
    dx = spsolve(J, -F)
    
    # Line search (optional but recommended for stability)
    alpha = 1.0
    x_new = x + alpha * dx
    
    # Check convergence
    if np.max(np.abs(dx)) < tol:
        converged = True
        break
    
    x = x_new
```

### Step 3: Line search for robustness

Add a simple backtracking line search:
```python
# If residual increases, halve the step
alpha = 1.0
for _ in range(10):
    x_trial = x + alpha * dx
    phi_e_t, phi_s_t = x_trial[:n_e], x_trial[n_e:]
    _, F_trial = self._build_coupled_jacobian(phi_e_t, phi_s_t, I_app)
    if np.linalg.norm(F_trial) < np.linalg.norm(F) * 1.1:
        break
    alpha *= 0.5
```

### Step 4: Preserve the warm-start

Keep the `_phi_e_guess` and `_phi_s_guess` for warm-starting between transient steps.

### Step 5: Keep the zero-current fast return

The `abs(I_app) < tol` early return (returning OCV) should stay as-is.

## Key Constraints

1. **Same interface**: `solve(I_app, tol, max_iter)` must return the same dict format: `phi_e`, `phi_s`, `voltage`, `I_rxn`, `iterations`, `converged`
2. **Same BCs**: electrolyte separator Dirichlet (phi_e=0), solid collector current BC, inactive node pinning
3. **Keep the existing Laplacian matrices** (`self.L_e`, `self.L_s`) — they're correct
4. **Keep the conductance computation** (`self.G_e`, `self.G_s`) — they're correct
5. **Sparse format**: use `scipy.sparse.lil_matrix` for assembly, convert to `csr` for solve
6. **Don't change**: `transient.py`, `reaction.py`, `ocv.py`, `generator.py`, `analysis.py`

## Verification

After implementing:

1. Run existing tests: `.venv/bin/pytest -q`
2. Test convergence on 5×5×5: should converge in ~5-10 Newton iterations
3. Test convergence on 10×10×10: should converge (this is where SOR failed)
4. Test convergence on 20×20×20: should converge (stretch goal)
5. Run a quick 0.2C discharge on 10×10×10 for 50 steps, verify voltage stays physical (2.5-4.5V)

```bash
cd /VOLUMES/1TB/projects/pnm-lib-cathode
.venv/bin/pytest -q
```

Then test:
```python
import sys; sys.path.insert(0, 'src')
from network.generator import create_cathode_network
from solver.transient import TransientSolver
import numpy as np

for shape in [[5,5,5], [10,10,10], [20,20,20]]:
    net = create_cathode_network(shape=shape, spacing=1e-5, porosity=0.5, cbd_fraction=0.1, seed=42)
    solver = TransientSolver(net, T=298.15, k0=5e-10)
    solver.set_concentration(c_e=1200.0, c_s=0.35*48900)
    I_app = solver.current_density_for_c_rate(0.2)
    
    # Test 50 steps
    dt = 20.0
    ok = True
    for i in range(50):
        result = solver.step(dt=dt, I_app=I_app)
        if not result['converged'] or not np.isfinite(result['voltage']):
            print(f'{shape}: FAILED at step {i}, V={result["voltage"]}')
            ok = False
            break
    if ok:
        mean_sol = np.mean(solver.c_s[solver.nmc_indices] / 48900.0)
        print(f'{shape}: OK, V={result["voltage"]:.3f}, mean_SoL={mean_sol:.3f}, iters={result["iterations"]}')
```
