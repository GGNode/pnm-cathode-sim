# PDF Extraction: Khan et al. 2021, JES 168 070534

Source read directly: `data/khan2021_pnm_lib_cathode.pdf`, 13 pages.

## Equations

### Eq. 2.1

```text
Li+ + e- ⇌ Li(s)
```

### Eq. 2.2

```text
v_i ((c_i^{t+Δt} - c_i^t) / Δt)
= Σ_{j=1}^{k_elec} [G_ij^d + max(-G_ij^m(φ_j^t - φ_i^t), 0)] c_i^t
  - Σ_{j=1}^{k_elec} [G_ij^d + max(G_ij^m(φ_i^t - φ_j^t), 0)] c_j^t
  + a j_{n,rp}^t
```

### Eq. 2.3

```text
G_ij^d = (1/g_i^d + 1/g_ij^d + 1/g_j^d)^-1
```

### Eq. 2.4

```text
G_ij^m = (1/g_i^m + 1/g_ij^m + 1/g_j^m)^-1
```

### Eq. 2.5

```text
g_i^d = A_i D_Li+ / l_i
```

### Eq. 2.6

```text
g_i^m = zF/(RT) g_i^d
```

### Eq. 2.7

```text
i_{i,elec} = -zF (Σ_{j=1}^{k_elec} G_ij^d(c_i^t - c_j^t))
             - Σ_{j=1}^{k_elec} K_ij^elec(φ_i^t - φ_j^t)
```

### Eq. 2.8

```text
K_ij^elec = (1/k_i^elec + 1/k_ij^elec + 1/k_j^elec)^-1
```

### Eq. 2.9

```text
k_i^elec = zF g_i^m c_i^t
```

### Eq. 2.10

```text
v_{i,AM} ((c_{i,Li}^{t+Δt} - c_{i,Li}^t) / Δt)
= Σ_{j=1}^{k_AM} G_{i,j}^{sd}(c_{i,Li}^t - c_{j,Li}^t)
```

### Eq. 2.11

```text
G_ij^sd = (1/g_i^sd + 1/g_ij^sd + 1/g_j^sd)^-1
```

### Eq. 2.12

```text
g_ij^sd = A_i D_Li / l_i
```

### Eq. 2.13

```text
i_{i,AM} = -Σ_{j=1}^{k_AM} G_ij^AM(φ_{i,AM}^t - φ_{j,AM}^t)
```

### Eq. 2.14

```text
0 = -Σ_{j=1}^{k_CBD} G_ij^CBD(φ_{i,CBD}^t - φ_{j,CBD}^t)
```

### Eq. 2.15

```text
g_ij^e = A_i σ_AM/CBD / l_i
```

### Eq. 2.16

```text
a j_{n,i}^t = ∇i_elec / zF
= i_i^{0,t} ( exp(((1 - α_c)F)/(RT) η_{i,c}^t)
              - exp((-α_c F)/(RT) η_{i,c}^t) )
```

### Eq. 2.17

```text
i_i^0 = Σ_{j=1}^{k_elec} (a_ij F k (c_{j,elec}^t)^(1-α_c)
          (c_{i,max}^t - c_{i,AM}^t)^(1-α_c) (c_{i,AM}^t)^(α_c))
```

### Eq. 2.18

```text
η_{i,c}^t = Σ_{j=1}^{k_elec} φ_{i,AM}^t - φ_{j,elec}^t - U_{i,cathode}^t
```

### Eq. 2.19

```text
U_{i,cathode} = 5744.862289 soc_i^9 - 5520.41099 soc_i^8
  + 95714.29862 soc_i^7 - 147364.5514 soc_i^6
  + 142718.3782 soc_i^5 - 90095.81521 soc_i^4
  + 37061.41195 soc_i^3 - 9578.599274 soc_i^2
  + 1409.309503 soc_i - 85.31153081
  - 0.0003 exp(7.657(soc_i^115))
```

### Eq. 2.20

```text
∇i_elec + ∇i_AM = 0
```

### Eq. 2.21

```text
∂c_{2,s}/∂t = -D_{Li,s}^eff ∇^2 c_{2,s}
```

### Eq. 2.22

```text
∂φ_{2,s}/∂x = -i_2/κ_s^eff + (RT/F)(1 - t_+) ∂ln c_{2,s}/∂x
```

### Eq. 2.23

```text
∂c_{2,s}/∂x = I(1 - t_+) / (F D_{Li,s})
```

### Eq. 2.24

```text
∇i_elec = i_foil^0 ( exp(((1 - α_a)F)/(RT) (φ_Li - φ_{2,s}))
                    - exp((-α_a F)/(RT) (φ_Li - φ_{2,s})) )
```

## Table I: Properties of XCT Images of NMC532 Cathodes

| Property | Units | NMC532 (1CAL) | NMC532 (3CAL) |
|---|---:|---:|---:|
| Coating Thickness | μm | 129 | 88 |
| Volume Fraction | % | — | — |
| Electrolyte phase | — | 36.8 | 36.6 |
| Active Material (NMC532) | — | 49.28 | 49.74 |
| Carbon and Binder Domain (CBD) | — | 13.92 | 13.66 |
| Density | mg/cm² | 29.78 | 20.40 |
| Experimental C-rate | mAh/g_NMC | 178 | 182 |
| Image size | voxels | 320 × 250 × 250 | 221 × 250 × 250 |
| Voxel Size | nm | 398 | 398 |
| Extracted Network | (Nodes, Bonds) | 4637, 31427 | 3510, 23126 |
| Electrolyte phase | — | 1238, 1943 | 982, 1459 |
| Active Material (NMC532) | — | 1722, 3736 | 1295, 2760 |
| Carbon and Binder Domain (CBD) | — | 1680, 3375 | 1233, 2310 |
| Interconnections | Bonds | 22373 | 16597 |
| NMC532-CBD | — | 7620 | 5355 |
| Electrolyte-CBD | — | 6602 | 4979 |
| Electrolyte-NMC532 | — | 8151 | 6263 |

## Table II: Model Parameters and Formula Details

Although Phase 1 asked for Table I, the requested key checks also require Table II.

| Parameter | Symbol | Units | Value |
|---|---|---:|---:|
| Separator Thickness | L_s | μm | 25 |
| Separator Porosity | ε_sep | % | 39 |
| Electronic conductivity of NMC532 | σ_AM | S/m | 0.01 |
| Electronic conductivity of CBD | σ_CBD | S/m | 760 |
| Ionic conductivity of the electrolyte in the separator | K_s | S/m | (1) |
| Initial Li+ concentration in the electrolyte | c_2,0 | mol/m³ | 1200 |
| Maximal Li concentration in NMC532 | c_1,max | mol/m³ | 48900 |
| Diffusion coefficient of Li+ in the electrolyte | D_Li+ | m²/s | `10^(-4.43 - 54.0/(T - 229 - 5c_2) - 0.22c_2)` |
| Li diffusion coefficient in AM | D_Li | m²/s | (2) |
| Faraday constant | F | C/mol | 96485 |
| Cathode rate constant | k | mol/m².s | 1e-10 |
| Bruggeman exponent coefficient for the separator | p_sep | — | 1.5 |
| Ideal gas constant | R | J/mol.K | 8.314 |
| Reference temperature | T | K | 303 |
| Li transference number in separator | t_+ | — | 0.363 |
| Cut off Voltage | V | V | 3.0 |
| Charge transfer coefficient of NMC | α_c | — | 0.5 |
| Charge transfer coefficient of Li foil | α_a | — | 0.5 |
| Exchange current density of Li foil | i_foil^0 | A/m² | 19 |

Table II footnote (1):

```text
10^(-2319soc^10 + 6642soc^9 - 5269soc^8 - 3319soc^7
    + 10038soc^6 - 9806soc^5 + 5817soc^4 - 2286soc^3
    + 575.3soc^2 - 83.16soc - 9.292)
```

Table II footnote (2):

```text
c_2(-10.5 + 0.0740T - 6.96e^-5 T^2 + 0.668c_2
   - 0.0178c_2T + 2.80e^-5 c_2T^2 + 0.494c_2^2
   - 8.86e^-4 c_2^2T)^2
```

Note: as printed, Table II assigns `(1)` to electrolyte ionic conductivity and `(2)` to AM Li diffusion, but the mathematical forms appear dimensionally/physically swapped: footnote `(1)` is a plausible solid diffusivity in `m²/s`, and footnote `(2)` is the common electrolyte conductivity expression in `S/m`.

## Methods Workflow

1. Use X-ray tomography images of two NMC532 cathodes from the NREL open-source electrode library.
2. Initial segmented cathode images contain labels 0 and 1: label 0 is electrolyte plus CBD, and label 1 is active material NMC532.
3. Add the CBD phase with image morphology functions; the paper states this method binds active-material contact points and is more suitable than other physics-based carbon binder models.
4. Extract a three-phase pore network from tomograms using an adapted watershed-based network extraction method for multiphase images.
5. Preserve interconnections between all phase pairs needed for electrolyte, NMC532, CBD, and interphase transport/reaction.
6. Model a Li foil anode, separator, and porous cathode half-cell. The cathode network contains electrolyte, NMC532 active material, and CBD.
7. During discharge, lithium ions leave Li foil, pass through the separator and cathode electrolyte, and intercalate into NMC532; electrons travel through the external circuit to the cathode current collector.
8. Do not include active-material expansion during lithium insertion.
9. In the cathode network, solve electrolyte concentration by discrete Nernst-Planck transport, electrolyte potential/current by the paper’s current-density form, NMC532 lithium transport by Fickian diffusion, NMC532/CBD electronic potential by Ohm’s law, and electrolyte/NMC reaction by Butler-Volmer kinetics.
10. Couple the separator and lithium foil regions to the pore-network cathode by 1D separator equations, not by a full separator pore network.
11. For galvanostatic discharge, set initial and boundary variables, parameters, and constants; solve concentration distribution in the electrolyte phase; solve charge conservation in electrolyte; solve lithium intercalation at the NMC532/electrolyte interface; solve lithium transport in NMC532; solve potential distribution in NMC532 and CBD; compare applied and calculated cathode current at the cathode current collector.
12. If `i_applied - i_calculated` is not below tolerance, guess a new cathode-current-collector potential and repeat the coupled solve.
13. When converged, accept the solution at the current time step, update variables, advance to the next time step, and repeat until cell voltage reaches the cutoff voltage.
14. The computational-performance section states two nested loops per time step: inner loop for the Nernst-Planck equation in electrolyte phase and outer loop for applied-vs-calculated current density convergence. Tolerances: inner `1e-8`, outer `1e-3`. Fixed time step: `2 s`.

## Boundary Conditions

Cathode region:

- No-flux boundary conditions for current and lithium-ion transport in the electrolyte phase at `x = L_s + L_c`.
- No-flux boundary conditions for lithium in the active material phase at the active material/separator interface and active material/cathode current collector interface.
- No electronic-flux boundary condition in the solid phase at the separator/solid phase (`CBD + AM`) interface at `x = L_s`.
- Continuity of mass, mass flux, charge, and charge flux at the separator/electrolyte phase interface at `x = L_s`.

Separator region:

- No-flux current boundary condition at the anode/separator interface at `x = 0`.
- Lithium solid-phase potential `φ_Li` set arbitrarily to zero at the anode/separator interface.
- Lithium-ion flux calculated by the boundary condition in Eq. 2.23.

Cell voltage:

- Defined as the difference in potential between the anode and cathode current collector.
- The separator and lithium foil regions are included by 1D modeling equations to supply separator/cathode boundary conditions and cell-potential estimate during galvanostatic discharge.

## Figure Descriptions

- Figure 1: For 1CAL and 3CAL, shows three-phase renderings, extracted pore networks, phase-pair throat interconnections, electrolyte/NMC size classes, electrolyte pore-size distribution versus separator distance, NMC particle-size distribution versus separator distance, sphericity versus volume for electrolyte pores and NMC particles, and phase size-distribution histograms.
- Figure 2: Half-cell PNM schematic with Li foil anode, separator, cathode, phase labels, cathode current collector, and unit cell/overlay views showing electrolyte, NMC532, CBD, single-phase connections, and two-phase interconnections.
- Figure 3: Galvanostatic-discharge workflow: initialize variables/constants; solve electrolyte concentration; solve electrolyte charge conservation; solve NMC532/electrolyte intercalation; solve NMC532 transport; solve NMC532/CBD potential; check applied-vs-calculated current; guess new cathode current collector potential if needed; accept current time-step solution; check cutoff; advance time.
- Figure 4: Voltage-capacity curves comparing PNM model and experimental data for 1CAL and 3CAL at 0.2C, 0.5C, and 1C, plus predicted 3C curves and a 1CAL-vs-3CAL comparison.
- Figure 5: Structural network metrics: coordination number for pores and particles by neighboring phase, slice-by-slice volume fraction along axes, electrolyte/NMC interfacial area versus separator distance for 1CAL/3CAL, and local tortuosity distributions for all phases.
- Figure 6: 3C spatial distributions at 75% SoL for 1CAL and 3CAL: electrolyte Li+ concentration, electrolyte potential, NMC state of lithiation, and solid-phase potential drop in NMC+CBD.
- Figure 7: Per-pore profiles for 1CAL and 3CAL at 1C and 3C: electrolyte Li+ concentration versus separator distance, electrolyte potential drop versus separator distance, and NMC state of lithiation, with markers categorized by pore/particle diameter classes.
- Figure 8: Deviation of current density from initial current density at 50%, 75%, and 100% SoL for 1CAL and 3CAL; marker color indicates active-material particle size.

## Specific Requested Details

- NMC532 OCV: the paper provides an explicit polynomial fit in Eq. 2.19; it is not only data points.
- Solid diffusivity: Table II prints AM diffusivity as footnote `(2)`, but footnote `(1)` is the physically plausible `D_s(soc)` expression. The printed footnote `(2)` is the electrolyte conductivity form.
- Electrolyte Li+ diffusivity: Table II gives `D_Li+ = 10^(-4.43 - 54.0/(T - 229 - 5c_2) - 0.22c_2)` with units `m²/s`.
- Exchange current density: Eq. 2.17 gives an interface-area-weighted sum over neighboring electrolyte pores: `i_i^0 = Σ(a_ij F k c_e^(1-α_c)(c_max - c_AM)^(1-α_c)c_AM^α_c)`.
- CBD conductivity: CBD is a separate conductive phase, no reaction occurs in CBD, and its electronic conductivity is `760 S/m`.
- Separator model: the separator is represented by 1D equations and boundary conditions coupled to the cathode pore network; it is not extracted as a pore network.
- Cutoff voltage: Table II and the workflow text state `3.0 V`.
- Transference number: Table II gives separator `t_+ = 0.363`, and Eqs. 2.22-2.23 use `(1 - t_+)`.
