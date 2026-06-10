# Final Gap Analysis and Implementation Plan

## Executive assessment

The core cathode PNM is now internally coherent: Butler-Volmer signs, coupled
solid/electrolyte potentials, concentration-dependent transport, OCV, initial
SoL, and end-of-step voltage alignment are in the right shape. The remaining
scientific gap is the separator/reference-electrode boundary. The current code
pins separator-side cathode electrolyte to an ideal `c_e = 1200 mol/m3`,
`phi_e = 0 V` reservoir, so it removes separator ohmic loss, separator salt
depletion, and Li-metal interfacial polarization from the voltage.

The geometry/mass-loading mismatch is a separate interpretation problem. A
10x10x10 synthetic network at 10 um spacing is not the paper's XCT network and
does not have the paper's areal loading. It can still support meaningful
method validation and qualitative rate-trend reproduction, but it should not be
presented as a quantitative reproduction of Khan Figure 4 unless the geometry
or scaling basis is made paper-consistent.

## 1. Separator model design

### Recommended minimal model

Implement a 1D separator boundary model as a finite-volume slab collapsed into
an impedance/transport boundary condition. This captures the dominant high-rate
effect without adding a full coupled separator grid to the Newton unknowns.

Use paper constants:

- `L_sep = 25e-6 m`
- `epsilon_sep = 0.39`
- `bruggeman_sep = 1.5`
- `t_plus = 0.363`
- `i0_foil = 19 A/m2`
- `alpha_foil = 0.5`
- `c_ref = 1200 mol/m3`

For each time step and applied current density:

```python
I_dis = abs(I_app)                    # positive discharge current density
c_sep_mean = c_ref                    # first pass; optionally iterate 2-3 times

for _ in range(n_separator_iterations):
    D_eff = D_e(c_sep_mean, T) * epsilon_sep**bruggeman_sep
    k_eff = kappa(c_sep_mean, T) * epsilon_sep**bruggeman_sep

    # Quasi-steady separator salt gradient, Li side at c_ref.
    dc_sep = I_dis * (1 - t_plus) * L_sep / (F * D_eff)
    c_sep_li = c_ref
    c_sep_cathode = max(c_floor, c_sep_li - dc_sep)
    c_sep_mean = 0.5 * (c_sep_li + c_sep_cathode)

# Potential drop from Li/reference side to cathode separator face.
dphi_ohm = I_dis * L_sep / k_eff
dphi_conc = (R * T / F) * (1 - t_plus) * log(c_sep_li / c_sep_cathode)

# Li foil BV overpotential, zero for disabled ideal Li reference.
eta_li_loss = (2 * R * T / F) * asinh(I_dis / (2 * i0_foil))

# Cathode-side electrolyte reference. For discharge this is negative.
phi_e_cathode_sep = -(dphi_ohm + dphi_conc + eta_li_loss)
```

Then replace the current cathode electrolyte Dirichlet value `phi_e = 0` at
`sep_e` with `phi_e = phi_e_cathode_sep`. The reported cell voltage should be
the cathode current-collector solid potential relative to the Li metal
reference, so with the Li reference fixed at zero:

```python
V_cell = mean(phi_s[collector_nodes])
```

Do not subtract `mean(phi_e[sep_e])` once the separator boundary is active;
that would remove the separator drop from the terminal voltage.

For electrolyte concentration, replace the existing direct reservoir clamp at
the cathode separator nodes with either:

1. Minimal first implementation: Dirichlet cathode-side separator concentration
   `c_e[sep_e] = c_sep_cathode`. This is simple and stable.
2. Better second implementation: Robin flux through separator:
   `N_sep = D_eff * A_cell / L_sep * (c_ref - c_sep_cathode_node)`, added to
   the separator-face electrolyte concentration residual. This avoids forcing
   all boundary pores to exactly the same concentration.

Start with option 1. It is enough to capture the dominant high-rate
concentration loss and is much easier to validate.

### Why not full separator Newton coupling first?

Khan Eqs. 2.21-2.24 can be represented by explicit 1D separator nodes coupled
to the cathode network and Li foil BV. That is more faithful, but it expands
the steady-state unknown vector, changes residual assembly, and adds another
nonlinear interface. The collapsed boundary model gives the main voltage and
salt-gradient effect with a small, auditable change. Once validated, it can be
upgraded to finite-volume separator cells without changing the cathode model.

## 2. C-rate strategy

Your assessment is mostly correct, with one important distinction:

- For a self-consistent synthetic-network study, use the network's own
  discharge capacity and projected geometric area for C-rate conversion:
  `TransientSolver.current_density_for_c_rate(C_rate)`. Normalize capacity by
  the network's own NMC mass loading. This is the meaningful default because
  the imposed current matches the amount of active material actually present.

- For quantitative Khan Figure 4 reproduction, paper current density only makes
  sense if the geometry and mass loading also match the paper. Applying
  `I_1C = 53 A/m2` to a network with 127.8 g/m2 effectively overdrives the
  synthetic electrode by about `297.8 / 127.8 = 2.33x` on a mass basis. That is
  useful as a stress test, not a fair C-rate comparison.

Recommended reporting:

1. Default validation mode: `current_basis = "network"`.
2. Optional sensitivity mode: `current_basis = "paper_areal"`.
3. Always report `mass_loading_g_m2`, `I_1C_A_m2`, and `current_basis` in
   validation metrics.

For publishable/meaningful results, use network C-rate for the synthetic
network and state that the comparison is qualitative unless the network is
resized or reconstructed to match the paper. If the goal is "paper-level
reproduction", build a geometry/mass-loading matched network and then use
paper-equivalent current.

## 3. Is the current physics paper-level?

Not yet for a paper-level reproduction claim.

Likely deal-breakers in review:

- Missing separator/reference model at high C-rate. This directly changes the
  voltage curve and accessible capacity.
- Synthetic network does not match XCT morphology, node/bond counts, thickness,
  interfacial area, or areal loading.
- Current-density basis is ambiguous if paper current is used on a thinner
  synthetic electrode.
- Electrolyte migration/concentration-potential terms from the paper are still
  simplified in the cathode network charge and species equations.

Acceptable simplifications if clearly disclosed:

- 10x10x10 synthetic network for method development and regression testing.
- Collapsed 1D separator boundary instead of explicit separator cells.
- Network-derived C-rate with network-derived mass normalization.
- No pore-throat-pore series conductance initially, if transport coefficients
  and throat geometry are stated as effective conductances.

The current model is good enough for an internal reproduction scaffold and for
qualitative physics exploration. It is not yet good enough to claim quantitative
reproduction of Khan 1CAL/3CAL curves.

## 4. Concrete execution plan

### MUST FIX: Separator boundary model

Files to modify:

- `src/physics/electrolyte.py`
- Add `src/physics/separator.py`
- `src/solver/steady.py`
- `src/solver/transient.py`
- `tests/test_steady.py`
- `tests/test_transient.py`
- `scripts/validate_paper_figures.py`
- `DERIVATION.md` or a short model note section in `README.md`

Add `src/physics/separator.py`:

```python
@dataclass(frozen=True)
class SeparatorParams:
    enabled: bool = False
    thickness: float = 25e-6
    porosity: float = 0.39
    bruggeman: float = 1.5
    t_plus: float = 0.363
    c_ref: float = 1200.0
    c_floor: float = 1.0
    include_concentration_overpotential: bool = True
    include_li_foil_bv: bool = False
    i0_foil: float = 19.0
    alpha_foil: float = 0.5

def separator_boundary(I_app, T, params):
    # returns c_cathode, phi_e_cathode, drops dict
```

Modify `SteadyStateSolver.__init__`:

```python
def __init__(..., separator: SeparatorParams | None = None):
    self.separator = separator or SeparatorParams(enabled=False)
    self.separator_state = separator_boundary(0.0, T, self.separator)
```

Modify `SteadyStateSolver.solve(I_app=...)` and residual/Jacobian boundary
rows:

```python
sep = separator_boundary(I_app, self.T, self.separator)
self.separator_state = sep
for i in self.sep_e:
    res[i] = phi_e[i] - sep.phi_e_cathode
```

Modify voltage reporting:

```python
if self.separator.enabled:
    V_cell = mean(phi_s[voltage_nodes])
else:
    V_cell = mean(phi_s[voltage_nodes]) - mean(phi_e[self.sep_e])
```

Modify `TransientSolver.__init__` to accept and pass `separator`:

```python
def __init__(..., separator: SeparatorParams | None = None):
    self.separator = separator or SeparatorParams(enabled=False)
    self._steady = SteadyStateSolver(net, T=T, k0=k0, separator=self.separator)
```

Modify electrolyte concentration BC in `TransientSolver.step`:

```python
if self.separator.enabled:
    c_bc = separator_boundary(I_app, self.T, self.separator).c_cathode
else:
    c_bc = self.c_e_reservoir

for i in self._steady.sep_e:
    M_e[i, :] = 0.0
    M_e[i, i] = 1.0
    rhs_e[i] = c_bc
```

Expected impact:

- 0.2C: small downward voltage shift, likely tens of mV or less.
- 3C: significant downward shift from separator ohmic plus concentration loss.
  This should steepen the high-rate curve and can change cutoff timing.
- The effect will not fix the paper-current overdrive problem on its own. With
  paper current on a thin network, 3C may become even more cutoff-limited.

Effort:

- Implementation: 0.5-1 day.
- Tests and validation rerun: 0.5 day.

### SHOULD FIX: C-rate basis and validation modes

Files to modify:

- `scripts/validate_paper_figures.py`
- `scripts/run_discharge.py`
- `src/solver/transient.py` only if adding a formal enum/API for current basis
- `tests/test_transient.py`

Change validation script constants:

```python
CURRENT_BASIS = "network"  # "network" or "paper_areal"
CAPACITY_BASIS = "network_mass"  # "network_mass" or "paper_mass"
```

Change current selection:

```python
if CURRENT_BASIS == "network":
    i_app = solver.current_density_for_c_rate(crate, geometric_area=A_GEOMETRIC)
elif CURRENT_BASIS == "paper_areal":
    i_app = -crate * I_1C_PAPER
```

Change capacity normalization:

```python
if CAPACITY_BASIS == "network_mass":
    result["capacity_mAh_g"] = result["capacity_Ah_m2"] / (mass_loading_g_m2 / 1000)
else:
    result["capacity_mAh_g"] = result["capacity_Ah_m2"] / PAPER_MASS_LOADING_KG_M2
```

Expected impact:

- Network current basis should raise 3C initial voltage relative to the current
  paper-current stress test and make the curves internally interpretable.
- Paper-current mode remains available to show why unmatched areal loading
  overpolarizes the synthetic electrode.

Effort:

- Implementation: 1-2 hours.
- Tests and validation: 1-2 hours.

### SHOULD FIX: Network sizing and geometry disclosure

Files to modify:

- `scripts/validate_paper_figures.py`
- `src/network/generator.py` only if adding helper metadata
- `data/validation/validation_metrics.json` generated by rerun
- `README.md` or `PLAN_FINAL.md` follow-up documentation

Do not immediately switch to 13x13x13 as the default unless runtime remains
acceptable. Instead:

```python
SHAPE = [10, 10, 10]
SPACING = 1e-5
THICKNESS_MODELED_UM = (SHAPE[0] - 1) * SPACING * 1e6  # coordinate span
```

Report:

- model thickness
- paper 1CAL thickness
- model mass loading
- paper mass loading
- ratio of paper/model mass loading
- node/bond counts and active interface count

Expected impact:

- No physics change.
- Prevents overclaiming and makes the current-basis decision auditable.

Effort:

- 1 hour.

### NICE TO HAVE: Pore-throat-pore series conductance

Files to modify:

- `src/network/generator.py` or new `src/physics/conductance.py`
- `src/solver/steady.py`
- `src/solver/transient.py`
- `tests/test_electrolyte.py`
- `tests/test_network.py`

Replace `G = property * throat.area / throat.length` with:

```python
R_total = L_half_i / (property_i * A_i) \
        + L_throat / (property_t * A_t) \
        + L_half_j / (property_j * A_j)
G = 1.0 / R_total
```

Expected impact:

- More defensible transport resistance.
- Possible additional polarization if pore bodies are a meaningful part of the
  path length.

Effort:

- 0.5-1 day.

### NICE TO HAVE: Migration and concentration-potential terms

Files to modify:

- `src/solver/steady.py`
- `src/solver/transient.py`
- `tests/test_steady.py`
- `tests/test_transient.py`

Add electrolyte current contribution:

```python
i_e = -kappa_eff grad(phi_e) - kappa_D_eff grad(log(c_e))
```

and update concentration transport to include migration/upwind terms if aiming
closer to Khan Eq. 2.2 and Eq. 2.7.

Expected impact:

- Most visible at high C-rate and large concentration gradients.
- Improves fidelity but increases Jacobian and stability complexity.

Effort:

- 1-2 days for a robust implementation.

## 5. Final validation strategy

### Unit tests

Add separator-specific tests:

1. `separator_boundary(I_app=0)` returns `phi_e_cathode = 0`,
   `c_cathode = c_ref`, and zero drops.
2. Separator ohmic drop is linear in small current:
   `drop(2I) / drop(I) ~= 2`.
3. Separator cathode concentration decreases monotonically with discharge
   current and stays bounded by `c_floor`.
4. With separator enabled and nonzero discharge current, `SteadyStateSolver`
   reports lower terminal voltage than the same case with separator disabled.
5. With separator disabled, existing zero-current OCV and mass-conservation
   tests remain unchanged.

### Analytical smoke scenario

Run a small 1D-like network or the existing 5x5x5 test network with BV nearly
linear at small current. Verify:

```text
V_no_sep - V_with_sep ~= I_dis * L_sep / kappa_eff
```

when concentration and Li foil terms are disabled. Then enable concentration
drop and check the additional voltage loss has the expected sign and increases
with current.

### End-to-end validation

Run four scenarios and save them separately:

1. Network current, separator off.
2. Network current, separator on.
3. Paper areal current, separator off.
4. Paper areal current, separator on.

For each, report:

- initial voltage
- cutoff capacity
- whether cutoff was reached
- `phi_e` span
- separator ohmic/concentration/Li-foil drops
- mass loading
- `I_1C`
- current basis

The final "correct" curve for this synthetic network should be scenario 2.
Scenario 4 should be labeled as a stress test against paper areal current, not
as a fair 3C reproduction.

### Acceptance criteria

- All unit tests pass.
- Separator-off results match current behavior within tolerance.
- Separator-on voltage is identical to separator-off at `I_app = 0`.
- Separator-on voltage decreases monotonically with `abs(I_app)`.
- Network-basis validation reaches physically plausible V-Q curves without
  hidden paper-mass normalization.
- Validation metrics explicitly state why the 10x10x10 synthetic network is not
  a quantitative geometry match to Khan 1CAL.
