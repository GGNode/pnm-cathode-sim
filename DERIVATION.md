# Derivation of the Cathode Pore Network Model

This document derives a transient pore network model (PNM) for a Li metal |
separator | NMC532 cathode half-cell under galvanostatic operation. It is
written as an implementation audit reference. No existing repository source code
was used in this derivation.

The model is a graph discretization of electrolyte salt transport, electrolyte
charge conservation, solid lithium transport, solid charge conservation, and
interfacial Butler-Volmer kinetics.

## 1. Geometry, Variables, and Sign Conventions

### 1.1 Network Sets

Use two coupled graphs plus an interface set:

- Electrolyte graph: nodes `i in E`, throats `(i,k) in T_e`.
- Solid graph: active-material/CBD nodes `m in S`, solid contacts `(m,n) in T_s`.
- Electrolyte/active-material interfaces: reactions `r in R`, connecting one
  electrolyte node `i(r)` to one active solid node `m(r)`.

Geometric quantities:

| Symbol | Unit | Meaning |
|---|---:|---|
| `V_i^e` | `m^3` | electrolyte pore control-volume volume |
| `V_m^s` | `m^3` | active solid control-volume volume |
| `A_ik^e` | `m^2` | electrolyte throat cross-sectional area |
| `L_ik^e` | `m` | electrolyte throat length |
| `A_mn^s` | `m^2` | solid contact cross-sectional area |
| `L_mn^s` | `m` | solid contact length |
| `A_r` | `m^2` | electrolyte/NMC reaction interface area |
| `L` | `m` | cathode thickness, separator at `x=0`, collector at `x=L` |

### 1.2 Unknown Fields

At each time level `t^{n+1}` the fully coupled unknown vector is

```text
x = [ c_e, phi_e, c_s, phi_s ]^T
```

with

| Symbol | Unit | Meaning |
|---|---:|---|
| `c_i^e` | `mol m^-3` | electrolyte Li salt concentration at electrolyte node `i` |
| `phi_i^e` | `V` | electrolyte potential at electrolyte node `i` |
| `c_m^s` | `mol m^-3` | solid-phase lithium concentration at active solid node `m` |
| `phi_m^s` | `V` | solid electronic potential at solid node `m` |

Constants:

| Symbol | Unit | Meaning |
|---|---:|---|
| `F` | `C mol^-1` | Faraday constant |
| `R` | `J mol^-1 K^-1` | gas constant |
| `T` | `K` | temperature |
| `D_e` | `m^2 s^-1` | electrolyte salt diffusion coefficient |
| `D_s` | `m^2 s^-1` | solid lithium diffusion coefficient |
| `kappa` | `S m^-1` | electrolyte ionic conductivity |
| `sigma` | `S m^-1` | solid electronic conductivity |
| `c_s,max` | `mol m^-3` | maximum lithium concentration in NMC |
| `U(c_s,c_e)` | `V` | equilibrium potential of NMC vs Li/Li+ |
| `alpha_a, alpha_c` | `1` | anodic and cathodic charge-transfer coefficients |

### 1.3 Reaction Sign Convention

Use the standard anodic convention for the intercalation reaction

```text
Li_s  <->  Li+_e + e-_s
```

Define:

- `q_r` [`mol m^-2 s^-1`] positive for anodic deintercalation/oxidation.
- `i_r = F q_r` [`A m^-2`] positive for anodic Faradaic current.
- During cathode discharge/lithiation, the cathode reaction is reduction, so
  `q_r < 0` and `i_r < 0`.

With this convention:

- Electrolyte Li is produced by `q_r > 0` and consumed by `q_r < 0`.
- Solid Li is consumed by `q_r > 0` and produced by `q_r < 0`.
- The reaction source appears with opposite signs in electrolyte and solid
  charge balances, so total interfacial charge is conserved exactly.

### 1.4 Current and Voltage Sign Convention

Let `I_app` [`A`] be the imposed conventional current entering the cathode
solid at the current collector. Equivalently, `I_app/A_cell` is the current
density applied at `x=L`.

- `I_app > 0`: anodic cathode operation, deintercalation/charge.
- `I_app < 0`: cathodic cathode operation, intercalation/discharge.

Many battery papers report discharge current as a positive C-rate magnitude.
For this document, a reported discharge magnitude `I_dis > 0` corresponds to
`I_app = -I_dis`.

## 2. Continuum Equations from First Principles

### 2.1 Electrolyte Nernst-Planck Transport

For species `j` with charge number `z_j`, dilute-solution Nernst-Planck flux is

```text
N_j = -D_j grad c_j - z_j u_j F c_j grad phi_e + c_j v
```

where `N_j` is molar flux [`mol m^-2 s^-1`], `u_j` is mobility
[`mol s kg^-1`], and `v` is volume-averaged electrolyte velocity. With no
convection and the Einstein relation `u_j = D_j/(R T)`,

```text
N_j = -D_j grad c_j - z_j (D_j F/(R T)) c_j grad phi_e.
```

For a binary electrolyte with `Li+` and an anion, electroneutrality gives
`c_+ = c_- = c_e`. The current density is

```text
i_e = F (N_+ - N_-).
```

For an ideal binary electrolyte this can be written

```text
i_e = -kappa(c_e) grad phi_e - kappa_D(c_e) grad ln c_e,
```

with

```text
kappa = (F^2/(R T)) (D_+ + D_-) c_e
kappa_D = (F/(R T)) (D_+ - D_-) R T c_e
```

up to the chosen dilute-solution convention. If diffusion-potential effects are
neglected, or if the electrolyte is approximated by a supporting electrolyte or
equal ion diffusivities, the charge law reduces to Ohm's law:

```text
i_e = -kappa(c_e) grad phi_e.
```

The simplified model in the project prompt uses this Ohmic form and a
Nernst-Einstein conductivity such as

```text
kappa(c_e) = F^2 D_e c_e/(R T).
```

This expression has units

```text
(C^2 mol^-2)(m^2 s^-1)(mol m^-3)/(J mol^-1)
= C^2 s kg^-1 m^-3 = A V^-1 m^-1 = S m^-1.
```

Electrolyte salt conservation is

```text
partial c_e/partial t = -div N_s + S_e,
```

where `N_s` is the salt flux. In the simplest Fickian salt model,

```text
N_s = -D_e,eff grad c_e.
```

The interfacial source depends on how the electrolyte model treats
transference:

```text
S_e = beta_e a_s i_F/F,
```

where `a_s` is interfacial area per volume [`m^-1`] and `i_F` is anodic
Faradaic current density [`A m^-2`]. For a full binary electrolyte model,
`beta_e = 1 - t_+^0`. For the simplified Li+ balance in the prompt,
`beta_e = 1`. The audit must verify which convention the implementation claims.

Thus the continuum salt equation used by the simplified PNM is

```text
partial c_e/partial t = div(D_e,eff grad c_e) + beta_e a_s i_F/F.
```

For discharge of the cathode, `i_F < 0`, so electrolyte Li is consumed.

### 2.2 Electrolyte Charge Conservation

Electroneutral electrolyte has negligible charge accumulation:

```text
div i_e = beta_q a_s i_F.
```

For the simplified single-current model, `beta_q = 1`. With `i_F > 0`,
positive ionic charge is generated at the interface and must be conducted away
through the electrolyte. With Ohm's law,

```text
div(-kappa_eff grad phi_e) = a_s i_F.
```

If the concentration-potential term is retained, replace the current by

```text
i_e = -kappa_eff grad phi_e - kappa_D,eff grad ln c_e.
```

### 2.3 Solid Lithium Diffusion from Fick's Law

Solid lithium in NMC is neutral intercalated lithium. Its flux is Fickian:

```text
N_solid = -D_s(c_s,T) grad c_s.
```

Mass conservation in active material is

```text
partial c_s/partial t = div(D_s grad c_s) - a_s i_F/F.
```

For `i_F > 0` deintercalation, solid lithium is consumed. For discharge,
`i_F < 0`, solid lithium increases.

If `D_s` depends on concentration, the continuum operator is

```text
div(D_s(c_s,T) grad c_s),
```

not `D_s grad^2 c_s` unless `D_s` is locally constant.

CBD domains conduct electrons but do not store intercalated lithium unless a
specific active storage model is added. Therefore `c_s` equations should be
assembled only on active NMC storage volumes, while `phi_s` may be assembled on
both NMC and CBD electronic nodes.

### 2.4 Solid Electronic Ohm's Law

Solid conventional current density is

```text
i_s = -sigma_eff grad phi_s.
```

Charge conservation with Faradaic transfer gives

```text
div i_s = -a_s i_F,
```

or

```text
div(-sigma_eff grad phi_s) = -a_s i_F.
```

Adding electrolyte and solid charge equations gives

```text
div(i_e + i_s) = 0,
```

so reaction only transfers charge between phases; it does not create net charge.

### 2.5 Butler-Volmer Kinetics

For

```text
Li_s <-> Li+_e + e-_s
```

the interfacial overpotential is

```text
eta_r = phi_m^s - phi_i^e - U(c_m^s, c_i^e).
```

Here `U` is the equilibrium potential of the active material vs Li/Li+ under
the same electrolyte reference. The anodic Butler-Volmer law is

```text
i_r = i0_r [ exp(alpha_a F eta_r/(R T))
             - exp(-alpha_c F eta_r/(R T)) ].
```

Therefore:

- `eta_r > 0` gives `i_r > 0`: deintercalation/oxidation.
- `eta_r < 0` gives `i_r < 0`: intercalation/reduction.
- `eta_r = 0` gives `i_r = 0`: local equilibrium.

A common exchange-current model is

```text
i0_r = F k0
       (c_i^e/c_e,ref)^gamma_e
       (c_s,max - c_m^s)^gamma_v
       (c_m^s)^gamma_s,
```

with exponents usually related to `alpha_a` and `alpha_c`. A common Li-ion
choice is

```text
gamma_e = alpha_a,
gamma_v = alpha_a,
gamma_s = alpha_c.
```

The implementation must keep units consistent. If concentrations are used in
`mol m^-3`, then `k0` must carry the complementary units needed to make
`i0` an `A m^-2`.

If no empirical OCV curve is available, an ideal Nernst form is

```text
U = U0 + (R T/F) ln[ ((c_s,max - c_s)/c_s) (c_e/c_e,ref) ].
```

For NMC532, an empirical `U(x)` with `x = c_s/c_s,max` is usually preferable.
If the empirical curve is already measured vs Li/Li+ at reference electrolyte
concentration, then electrolyte concentration should not be double-counted
unless the paper explicitly includes a Nernst correction.

## 3. Pore Network Discretization

The finite-volume control volume is the pore or solid node. Throats provide
two-point flux approximations. For every edge, use harmonic or series
averaging when material properties differ across adjacent half-throats.

### 3.1 Edge Conductances

Electrolyte diffusive conductance:

```text
K_ik^e = D_ik,eff^e A_ik^e/L_ik^e        [m^3 s^-1]
```

Electrolyte ionic conductance:

```text
G_ik^e = kappa_ik,eff A_ik^e/L_ik^e      [S = A V^-1]
```

Solid diffusive conductance:

```text
K_mn^s = D_mn^s A_mn^s/L_mn^s            [m^3 s^-1]
```

Solid electronic conductance:

```text
G_mn^s = sigma_mn,eff A_mn^s/L_mn^s      [S]
```

For Bruggeman-type corrections, for example,

```text
D_e,eff = D_e epsilon^b
kappa_eff = kappa epsilon^b
```

with `b = 1.5` in the prompt. The audit should confirm whether the correction
is applied to pore-scale throats, continuum separator regions, or both. Applying
both pore geometry and Bruggeman correction can double-count tortuosity.

### 3.2 Electrolyte Concentration Residual

For electrolyte node `i`, define the diffusive molar flow into `i` from
neighbor `k`:

```text
M_ki^e = K_ik^e (c_k^e - c_i^e)          [mol s^-1]
```

because `K [m^3 s^-1] * Delta c [mol m^-3] = mol s^-1`.

Let `R_i` be the set of reaction interfaces attached to electrolyte node `i`.
Backward Euler gives

```text
F_{c_e,i} =
  V_i^e (c_i^{e,n+1} - c_i^{e,n})/dt
  - sum_{k in N_e(i)} K_ik^e (c_k^{e,n+1} - c_i^{e,n+1})
  - sum_{r in R_i} beta_e (i_r^{n+1}/F) A_r
  = 0.
```

Sign check:

- If isolated and `i_r > 0`, then `dc_e/dt = beta_e i_r A_r/(F V_i^e) > 0`.
- If isolated and `i_r < 0`, then electrolyte concentration decreases.

### 3.3 Electrolyte Potential Residual

Define conductive current entering electrolyte node `i` from neighbor `k`:

```text
I_ki^e = G_ik^e (phi_k^e - phi_i^e)      [A]
```

The quasi-static charge balance is

```text
F_{phi_e,i} =
  sum_{k in N_e(i)} G_ik^e (phi_k^e - phi_i^e)
  + sum_{r in R_i} i_r A_r
  = 0.
```

If concentration-potential terms are retained, add to each edge current

```text
I_ki^{e,conc} = H_ik^e (ln c_k^e - ln c_i^e),
```

where `H_ik^e` has units of amperes and follows from the chosen
Nernst-Planck/concentrated-solution model.

Sign check:

- For `i_r > 0`, electrolyte charge is produced at the interface, so the
  conductive current entering from neighboring electrolyte nodes must be
  negative; current leaves the node through the electrolyte graph.
- The reaction term cancels the opposite term in the solid equation.

### 3.4 Solid Concentration Residual

For active solid node `m`, define diffusive molar flow into `m` from active
neighbor `n`:

```text
M_nm^s = K_mn^s (c_n^s - c_m^s).
```

Let `R_m` be the set of reaction interfaces attached to active solid node `m`.
The backward Euler residual is

```text
F_{c_s,m} =
  V_m^s (c_m^{s,n+1} - c_m^{s,n})/dt
  - sum_{n in N_s(m)} K_mn^s (c_n^{s,n+1} - c_m^{s,n+1})
  + sum_{r in R_m} (i_r^{n+1}/F) A_r
  = 0.
```

Sign check:

- If isolated and `i_r > 0`, then `dc_s/dt = -i_r A_r/(F V_m^s) < 0`.
- If isolated and `i_r < 0`, then solid lithium concentration increases.

For nonlinear `D_s(c_s,T)`, `K_mn^s` must be evaluated at `t^{n+1}` in a fully
implicit Newton solve. A stable edge choice is harmonic averaging of
`D_s(c_m)` and `D_s(c_n)`.

### 3.5 Solid Potential Residual

Define conventional electronic current entering solid node `m` from neighbor
`n`:

```text
I_nm^s = G_mn^s (phi_n^s - phi_m^s).
```

The quasi-static charge residual is

```text
F_{phi_s,m} =
  sum_{n in N_s(m)} G_mn^s (phi_n^s - phi_m^s)
  - sum_{r in R_m} i_r A_r
  + F_{BC,m}^s
  = 0.
```

`F_{BC,m}^s` contains applied collector current contributions. With this sign
convention, the reaction terms in `F_{phi_e}` and `F_{phi_s}` are exact
opposites.

### 3.6 Generic Coupled Interface Contribution

For one interface `r` between electrolyte node `i` and solid node `m`, the
reaction current contributes:

```text
F_{c_e,i}     += - beta_e (i_r/F) A_r
F_{phi_e,i}   += + i_r A_r
F_{c_s,m}     += + (i_r/F) A_r
F_{phi_s,m}   += - i_r A_r
```

This four-line block is the central coupling of the model.

## 4. Boundary Conditions and Voltage Reference

### 4.1 Separator/Electrolyte Boundary at `x=0`

The prompt describes a separator connected to a reservoir with
`c_e = 1200 mol m^-3` and `phi_e = 0 V`. This is a Dirichlet boundary for the
cathode electrolyte graph:

```text
c_i^e = c_e,0
phi_i^e = 0
```

for electrolyte nodes on the separator face, or an equivalent ghost-node/Robin
condition if separator resistance is explicitly modeled.

If the separator is modeled as a finite 1D region of thickness `L_sep` and
porosity `epsilon_sep`, then:

```text
D_sep,eff = D_e epsilon_sep^b
kappa_sep,eff = kappa epsilon_sep^b
```

Boundary at Li metal/reservoir:

```text
c_e = c_e,0
phi_e = 0
```

Interface to cathode: continuity of salt flux and ionic current.

If the separator is not modeled explicitly, its concentration and ohmic drops
are omitted. For high-rate validation, that omission can matter.

### 4.2 Collector Boundary at `x=L`

At the cathode current collector:

Electrolyte:

```text
N_s . n = 0
i_e . n = 0
```

There is no electrolyte beyond the current collector, so no salt or ionic
current leaves through that boundary.

Solid:

```text
integral_{collector} i_s . n dA = I_app.
```

In graph form, distribute `I_app` over collector-connected solid nodes by area
weights `w_m` satisfying `sum w_m = 1`:

```text
F_{BC,m}^s = I_app w_m
```

if the residual is written as "incoming currents plus sources equals zero" and
`I_app` is current entering the cathode solid from the external circuit. The
implementation must verify this sign by checking that a discharge input
`I_app < 0` produces negative Faradaic currents and a voltage below OCV.

Alternative formulation: pin `phi_s` at the collector and solve the resulting
total current, then use an outer scalar Newton/secant loop to adjust collector
potential until `I_total = I_app`. This is often more robust for polarization
curves, but the fully coupled Neumann formulation is natural for galvanostatic
time stepping.

### 4.3 Solid Boundary at Separator Face

No electronic current crosses into the separator:

```text
i_s . n = 0
```

unless an electronically conducting short is intentionally modeled.

For solid lithium, non-reacting external boundaries are no-flux:

```text
N_solid . n = 0.
```

Lithium enters or leaves active material only through electrolyte/NMC
interfaces in this model.

### 4.4 Potential Gauge and Reference Frame

Only potential differences are physical. The coupled equations are invariant to
adding a constant to all `phi_e` and `phi_s` unless one potential is pinned. A
well-posed half-cell cathode model must set one reference.

The prompt's reference is:

```text
phi_e at separator/Li reference = 0 V.
```

For a Li metal counter/reference electrode at equilibrium,

```text
phi_s,Li - phi_e,sep = U_Li/Li+ = 0 V
```

by definition of the Li/Li+ scale. Therefore, with `phi_e,sep = 0`,

```text
phi_s,Li = 0.
```

The cathode NMC OCV `U` is also measured vs this Li/Li+ reference, so
equilibrium at a uniform cathode state is

```text
phi_s,cathode - phi_e,cathode = U(c_s,c_e).
```

### 4.5 Correct Cell Voltage Formula

The measurable half-cell voltage is positive-current-collector potential minus
Li metal potential:

```text
V_cell = phi_s,collector - phi_s,Li.
```

Using the reference above:

```text
V_cell = phi_s,collector.
```

More generally, if the electrolyte separator potential is not zero,

```text
V_cell = phi_s,collector
         - [ phi_e,sep + U_Li/Li+ + eta_Li ].
```

For an ideal Li metal reference, `U_Li/Li+ = 0` and `eta_Li = 0`, so

```text
V_cell = phi_s,collector - phi_e,sep.
```

If `phi_e,sep` is pinned to zero, this again reduces to `phi_s,collector`.

Important audit point: `V_cell = phi_s,collector - phi_e,separator` is correct
only if the Li metal electrode is an equilibrium Li/Li+ reference and
`phi_e,separator` is the electrolyte potential at that reference. It is not a
generic full-cell voltage formula.

## 5. Newton-Raphson Formulation

### 5.1 Residual Vector

At each time step, solve

```text
F(x^{n+1}; x^n, dt, I_app) = 0
```

with block structure

```text
F =
[
  F_c_e
  F_phi_e
  F_c_s
  F_phi_s
].
```

Backward Euler makes both concentration equations implicit. Potential equations
are algebraic quasi-steady constraints. The system is a DAE-like nonlinear
algebraic solve at each time step.

Newton iteration `ell`:

```text
J(x_ell) delta x_ell = -F(x_ell)
x_{ell+1} = x_ell + lambda delta x_ell
```

where `0 < lambda <= 1` is a damping/line-search parameter.

### 5.2 Butler-Volmer Derivatives

For interface `r`, define

```text
f = F/(R T)
eta = phi_s - phi_e - U(c_s,c_e)
E_a = exp(alpha_a f eta)
E_c = exp(-alpha_c f eta)
i = i0 (E_a - E_c)
B = d i/d eta = i0 f (alpha_a E_a + alpha_c E_c).
```

Direct potential derivatives:

```text
d i/d phi_s = +B
d i/d phi_e = -B
```

Concentration derivatives:

```text
d i/d c_s =
  (d i0/d c_s)(E_a - E_c)
  - B (dU/d c_s)

d i/d c_e =
  (d i0/d c_e)(E_a - E_c)
  - B (dU/d c_e).
```

For

```text
i0 = F k0 c_e^gamma_e (c_s,max - c_s)^gamma_v c_s^gamma_s,
```

the log derivatives are

```text
d i0/d c_e = i0 gamma_e/c_e
d i0/d c_s = i0 [ gamma_s/c_s - gamma_v/(c_s,max - c_s) ].
```

These derivatives are singular at `c_e = 0`, `c_s = 0`, and
`c_s = c_s,max`; practical implementations must bound concentrations away from
these endpoints.

### 5.3 Jacobian Contributions from One Interface

Let `p` be any local variable among `c_e_i`, `phi_e_i`, `c_s_m`, `phi_s_m`.
The interface contributions to the Jacobian are:

```text
dF_{c_e,i}/dp   += - beta_e (A_r/F) d i_r/dp
dF_{phi_e,i}/dp += + A_r d i_r/dp
dF_{c_s,m}/dp   += + (A_r/F) d i_r/dp
dF_{phi_s,m}/dp += - A_r d i_r/dp
```

This is where `phi_e` and `phi_s` are coupled. The derivatives with respect to
`phi_e` and `phi_s` are equal and opposite because the reaction depends on
`phi_s - phi_e`.

### 5.4 Jacobian Contributions from Linear Edges

For electrolyte concentration residual

```text
F_i = V_i(c_i-c_i^n)/dt - sum_k K_ik(c_k-c_i) - sources,
```

with constant `K_ik`:

```text
dF_i/dc_i += V_i/dt + sum_k K_ik
dF_i/dc_k += -K_ik.
```

For electrolyte potential residual

```text
F_i = sum_k G_ik(phi_k-phi_i) + sources,
```

with constant `G_ik`:

```text
dF_i/dphi_i += -sum_k G_ik
dF_i/dphi_k += +G_ik.
```

For solid concentration:

```text
dF_m/dc_m += V_m/dt + sum_n K_mn
dF_m/dc_n += -K_mn.
```

For solid potential:

```text
dF_m/dphi_m += -sum_n G_mn
dF_m/dphi_n += +G_mn.
```

If `K` or `G` depends on concentration, additional terms are required. For an
edge residual contribution `-K(c_i,c_k)(c_k-c_i)`,

```text
d/dc_i = - (dK/dc_i)(c_k-c_i) + K
d/dc_k = - (dK/dc_k)(c_k-c_i) - K.
```

For a potential edge `G(c_i,c_k)(phi_k-phi_i)`,

```text
d/dphi_i = -G
d/dphi_k = +G
d/dc_i = (dG/dc_i)(phi_k-phi_i)
d/dc_k = (dG/dc_k)(phi_k-phi_i).
```

### 5.5 Dirichlet Boundary Rows

For a pinned unknown `y = y_B`, replace the corresponding residual row by

```text
F_y = y - y_B = 0
```

and the Jacobian row by

```text
dF_y/dy = 1
```

with all other entries in that row zero. This avoids mixing conservation rows
with boundary rows.

### 5.6 Galvanostatic Constraint Options

Two correct strategies exist.

Option A: Neumann current boundary in `F_phi_s`.

```text
sum_m F_{BC,m}^s = I_app.
```

Then `phi_s,collector` is solved as an output, and
`V_cell = average_or_terminal(phi_s,collector) - phi_s,Li`.

Option B: Unknown collector voltage plus scalar current constraint.

Pin collector solid nodes to a common unknown `V_col`, solve reaction and
transport equations, and add

```text
F_I = sum_{r in R} i_r A_r - I_app = 0
```

with the appropriate sign and any double-layer/storage terms if included. This
adds one scalar unknown and one scalar residual. It is useful when enforcing a
perfectly equipotential current collector.

In either option, at steady galvanostatic operation:

```text
sum_r i_r A_r = I_app
```

for a cathode with no capacitive side reactions, using the anodic convention.

## 6. Numerical Pitfalls and Remedies

### 6.1 Exponential Overflow in Butler-Volmer

At `T = 298 K`, `R T/F approximately 25.7 mV`. Exponents contain
`alpha F eta/(R T)`. With `alpha=0.5`, `eta = 1 V` gives exponent about
`19.5`, and larger excursions can overflow or create unusable Jacobian scales.

Remedies:

- Clip exponent arguments to a safe range, for example `[-80, 80]`, while
  preserving derivative consistency.
- Use `sinh` form for symmetric BV:
  `i = 2 i0 sinh(F eta/(2 R T))`, with stable `sinh` evaluation.
- Use Tafel asymptotes only deliberately; switching formulas can make the
  Jacobian nonsmooth.
- Apply Newton damping or line search when `||F||` increases.

### 6.2 Concentration Bounds

The model is physically meaningful only for

```text
c_e > 0
0 < c_s < c_s,max.
```

Failures near bounds:

- `ln(c)` in Nernst terms becomes undefined.
- `i0` derivatives diverge.
- Empirical OCV fits often diverge or become invalid outside calibrated `x`.
- Negative concentration can make `kappa(c)` negative.

Remedies:

- Use bound-preserving line search.
- Reject Newton steps that leave admissible concentration intervals.
- Use transformed unknowns such as `log c_e` and `logit(c_s/c_s,max)` for hard
  cases, though this complicates residual scaling.
- Keep time steps small enough that physical depletion is resolved.

### 6.3 Singular Matrices from Potential Gauge Freedom

If no electrolyte potential is pinned, adding a constant to all potentials
leaves the equations unchanged. The Jacobian has a null vector. Pin exactly one
reference potential, normally `phi_e = 0` at the separator/Li reference.

Do not also pin an inconsistent solid potential under galvanostatic operation
unless the current is solved as an output or an extra constraint is added.

### 6.4 Disconnected Pores or Solids

Disconnected electrolyte clusters not connected to the separator Dirichlet
condition have floating `phi_e` and may trap concentration. Disconnected solid
clusters not connected to the collector have floating `phi_s`. If they also
contain active interfaces, BV equations can drive unphysical local reactions
unless electronic/ionic pathways exist.

Required preprocessing:

- Keep only electrolyte components connected to the separator.
- Keep only electronic solid components connected to the collector for active
  reaction, or mark isolated active material as electrochemically inactive.
- Verify every reactive interface has both ionic and electronic connectivity.

### 6.5 Ill Conditioning from Conductivity Contrast

CBD conductivity (`~760 S m^-1`) and NMC conductivity (`~0.01 S m^-1`) can
differ by many orders of magnitude. Electrolyte and solid concentration blocks
also have different physical units from potential blocks.

Remedies:

- Use sparse matrix scaling or nondimensional residuals.
- Scale concentration residuals by `V c_ref/dt` or `I_ref/F`.
- Scale charge residuals by `I_ref`.
- Prefer robust sparse direct solvers for small/medium networks and
  preconditioned Krylov methods for large networks.

### 6.6 Time Step Selection

Backward Euler is unconditionally stable for linear diffusion, but nonlinear
accuracy and Newton convergence still impose limits.

Useful scales:

```text
tau_e ~ l_pore^2/D_e,eff
tau_s ~ l_particle^2/D_s
tau_rxn ~ F c_s,max V_active/(|I_app|)
```

Practical strategy:

- Start with small `dt` after current onset.
- Increase `dt` only when Newton converges in few iterations.
- Cut `dt` after line-search failure, concentration-bound violation, or
  excessive voltage change.
- Resolve cutoff voltage crossing by interpolation or step rejection.

### 6.7 Mass and Charge Conservation Audits

At every converged time step check:

Electrolyte lithium:

```text
Delta n_e =
sum_i V_i^e (c_i^{e,n+1} - c_i^{e,n})
```

should equal integrated boundary salt flux plus

```text
dt sum_r beta_e i_r A_r/F.
```

Solid lithium:

```text
Delta n_s =
sum_m V_m^s (c_m^{s,n+1} - c_m^{s,n})
= -dt sum_r i_r A_r/F
```

up to solid diffusive boundary fluxes, normally zero.

Charge:

```text
sum_r i_r A_r = I_app
```

for a galvanostatic cathode with no double-layer capacitance and no side
reaction.

### 6.8 Sign Convention Bugs

The most common implementation errors are:

- Using a BV convention where positive current means reduction while the charge
  equations assume positive current means oxidation.
- Reporting discharge C-rate as positive but passing it into equations as
  `I_app > 0`.
- Computing `eta = phi_e - phi_s - U` instead of `phi_s - phi_e - U`.
- Updating both electrolyte and solid concentrations with the same reaction
  sign.
- Computing `V_cell` from a local interfacial overpotential instead of the
  terminal solid potential relative to Li/Li+.

A single-pore equilibrium and small-current test should catch these errors.

## 7. Known-Limit Validation

### 7.1 Zero-Current Equilibrium (`0C`)

Set `I_app = 0`. With uniform initial concentrations and connected phases, the
expected solution is

```text
c_e = c_e,0 everywhere
c_s = c_s,0 everywhere
phi_e = 0 everywhere, if separator electrolyte is pinned to 0
eta_r = 0 at every interface
i_r = 0 at every interface
phi_s = U(c_s,0, c_e,0) everywhere in the connected solid
V_cell = U(c_s,0, c_e,0).
```

If empirical `U(x)` is used without electrolyte correction, then

```text
V_cell = U(x_0).
```

If initial `c_s` is spatially nonuniform, true equilibrium in a connected,
conductive, reactive solid requires a uniform electrochemical potential. In the
simplified model this usually relaxes toward uniform `U`, not necessarily
uniform `c_s` if a nonideal thermodynamic model is used. A code initialized
with nonuniform `c_s` should not show zero reaction everywhere unless
`phi_s - phi_e = U(c_s)` at every interface, which is impossible with uniform
connected potentials unless `U(c_s)` is uniform.

### 7.2 Small-Current / Low-C-Rate Limit

For small current, linearize BV:

```text
i_r approximately i0_r (F/(R T))(alpha_a + alpha_c) eta_r.
```

For `alpha_a + alpha_c = 1`,

```text
eta_r approximately (R T/(F i0_r)) i_r.
```

The cell voltage under discharge (`I_app < 0`) should be

```text
V_cell approximately U(x_avg) - |I_app| R_total
```

where `R_total` includes reaction, electrolyte, solid, separator, and contact
resistances. The voltage error from OCV should scale linearly with current.

### 7.3 High-C-Rate Electrolyte Depletion Limit

For cathode discharge, `i_r < 0`, electrolyte Li is consumed in the cathode.
In a 1D slab approximation with separator at `x=0`, collector at `x=L`,
uniform reaction, no-flux at the collector, and fixed `c(0)=c0`, define the
positive consumption current magnitude per area `I_dis = -I_app > 0`.

The volumetric Li consumption rate is approximately

```text
s = beta_e I_dis/(F L)       [mol m^-3 s^-1].
```

At quasi-steady state,

```text
0 = D_e,eff d^2 c/dx^2 - s
dc/dx(L) = 0
c(0) = c0.
```

The solution is

```text
c(x) = c0 + (s/D_e,eff)(x^2/2 - L x).
```

The minimum concentration occurs at the collector:

```text
c(L) = c0 - s L^2/(2 D_e,eff)
     = c0 - beta_e I_dis L/(2 F D_e,eff).
```

An approximate diffusion-limited discharge current density is therefore

```text
I_lim approximately 2 F D_e,eff c0/(beta_e L).
```

As `I_dis` approaches this scale, `c_e` near the collector approaches zero,
`i0` collapses, the required negative overpotential diverges, and accessible
capacity drops sharply.

### 7.4 High-C-Rate Ohmic Limit

Ignoring concentration gradients and reaction nonuniformity, distributed current
in a slab gives approximate half-cell polarization

```text
Delta V_ohm approximately I_dis
  [ L/(2 kappa_eff A_cell) + L/(2 sigma_eff A_cell) + R_sep + R_contact ].
```

For discharge, terminal voltage is below OCV:

```text
V_cell approximately U - Delta V_ohm - |eta_ct| - Delta V_conc.
```

The factors of `1/2` arise because ionic current is largest near the separator
and zero at the collector, while solid current is largest near the collector
and zero near the separator under uniformly distributed reaction. Strong
reaction localization invalidates the `1/2` estimate but not the sign.

### 7.5 Solid Diffusion Limit

At high discharge rate, lithium enters particle surfaces faster than it
diffuses inward. For a characteristic particle radius `R_p`,

```text
tau_s ~ R_p^2/D_s.
```

If the discharge time `tau_dis` is much smaller than `tau_s`, surface
concentration rises toward `c_s,max` while the particle core remains less
lithiated. Consequences:

- `U(c_s,surf)` shifts toward the lithiated end of the OCV curve.
- `i0` decreases as vacancy concentration `c_s,max - c_s,surf` vanishes.
- Cutoff occurs before full active-material utilization.

A single-particle validation should reproduce the known short-time diffusion
scaling. For constant flux into a semi-infinite solid, surface concentration
change scales as

```text
Delta c_s,surf ~ (2 j/F?) sqrt(t/(pi D_s))
```

with the exact prefactor depending on whether `j` is expressed as molar flux or
current density. If molar flux into the solid is `J_in = -i_F/F > 0`, then

```text
Delta c_s,surf = 2 J_in sqrt(t/(pi D_s)).
```

## 8. Ideal Module Architecture

The code should make sign conventions and units explicit. A useful architecture
is:

```text
src/
  constants.py
  network/
    geometry.py
    connectivity.py
    preprocessing.py
  physics/
    electrolyte.py
    solid.py
    reaction.py
    ocv.py
    transport.py
    boundary.py
  solver/
    state.py
    residual.py
    jacobian.py
    newton.py
    timestepping.py
  validation/
    analytic_limits.py
  post/
    voltage.py
    conservation.py
    visualization.py
```

Recommended responsibilities:

- `constants.py`: `F`, `R`, default `T`, unit helpers.
- `network.geometry`: node volumes, throat areas/lengths, interface areas.
- `network.connectivity`: graph incidence arrays, component detection.
- `network.preprocessing`: remove inactive disconnected components and build
  boundary node sets.
- `physics.transport`: conductance construction with harmonic averaging and
  Bruggeman corrections.
- `physics.electrolyte`: electrolyte concentration and charge edge residuals.
- `physics.solid`: solid concentration and charge edge residuals.
- `physics.reaction`: BV current, stable exponentials, analytical derivatives.
- `physics.ocv`: NMC532 OCV and derivatives `dU/dc_s`, optional `dU/dc_e`.
- `physics.boundary`: separator Dirichlet rows, collector Neumann/current rows,
  optional separator 1D impedance.
- `solver.state`: typed packing/unpacking of `[c_e, phi_e, c_s, phi_s]`.
- `solver.residual`: assemble `F(x)`.
- `solver.jacobian`: assemble sparse analytical `J(x)`.
- `solver.newton`: damping, bounds, convergence criteria, linear solves.
- `solver.timestepping`: adaptive `dt`, cutoff detection, output history.
- `post.voltage`: terminal voltage definitions and collector averaging.
- `post.conservation`: lithium and charge conservation audits.

Architecture rule: the reaction module should expose one sign convention and
all residual assemblers should use it. Do not duplicate BV formulas in multiple
places.

## 9. Test Strategy

### 9.1 Unit Tests

1. Butler-Volmer zero:
   - At `eta=0`, `i=0`.
   - `di/deta = i0 F(alpha_a+alpha_c)/(R T)`.

2. Butler-Volmer signs:
   - Positive `eta` gives positive anodic current.
   - Negative `eta` gives negative cathodic current.

3. BV Jacobian:
   - Compare analytical derivatives with finite differences for `phi_s`,
     `phi_e`, `c_s`, and `c_e`.

4. OCV derivatives:
   - Verify `dU/dc_s` and `dU/dc_e` against finite differences.
   - Verify behavior near concentration bounds is clipped or rejected.

5. Edge flux signs:
   - If `c_k > c_i`, diffusive flow into `i` is positive.
   - If `phi_k > phi_i`, positive current entering `i` is positive.

6. Dirichlet rows:
   - Boundary rows must be identity rows and enforce exact pinned values.

### 9.2 Small Network Tests

1. Single reactive interface, no transport:
   - With `I_app=0`, solve `eta=0` and recover `phi_s - phi_e = U`.
   - With imposed small current, compare `eta` to linearized BV.

2. Two-node diffusion:
   - Compare backward Euler update to analytical matrix solution.
   - Verify total mass conservation with no boundary flux.

3. Two-phase charge conservation:
   - One electrolyte node and one solid node connected by BV plus one current
     boundary should satisfy `i_r A = I_app`.

4. Disconnected component:
   - Solver should reject or deactivate a reactive solid island not connected
     to the collector.

### 9.3 Continuum-Limit Tests

1. 1D electrolyte slab:
   - Uniform reaction and fixed separator concentration should reproduce
     `c(x) = c0 + (s/D)(x^2/2 - Lx)`.

2. 1D Ohmic slab:
   - Distributed reaction should reproduce the `L/(2 kappa A)` electrolyte
     resistance scaling and `L/(2 sigma A)` solid scaling.

3. Low-rate voltage:
   - As `I_app -> 0`, `V_cell -> U(x)` and polarization is linear in current.

4. High-rate depletion:
   - Increasing discharge current should move the minimum `c_e` toward the
     collector and approach the analytical limiting-current estimate.

5. Global conservation:
   - At every time step, integrated solid lithium change must equal
     `-dt sum_r i_r A_r/F` within solver tolerance.
   - Total reaction current must match applied galvanostatic current.

## 10. Implementation Audit Checklist

Use this checklist when reviewing an implementation:

- Does the code document whether positive current is anodic or cathodic?
- Is discharge C-rate converted to `I_app < 0` under the anodic convention?
- Is `eta = phi_s - phi_e - U`?
- Are electrolyte and solid reaction source signs opposite?
- Does `V_cell` reference the Li/Li+ electrode correctly?
- Is exactly one potential gauge pinned?
- Are disconnected ionic/electronic components removed or deactivated?
- Are BV exponentials stabilized?
- Are concentration bounds enforced?
- Are analytical Jacobian signs consistent with residual signs?
- Are `kappa(c)` and `i0(c)` derivatives included if they are evaluated
  implicitly?
- Are units SI throughout: meters, seconds, mol/m^3, volts, amperes?
- Is Bruggeman/tortuosity correction applied once, not double-counted?
- Does `0C` recover OCV with zero reaction current?
- Does small-current polarization scale linearly with current?
- Does high-rate discharge deplete electrolyte toward the collector?
