# 阴极孔网络模型推导

本文推导一个用于 Li metal | separator | NMC532 cathode 半电池恒流工况的瞬态 pore network model (PNM)。文档作为实现审查参考，技术术语、变量名与代码标识符保留英文。

模型把电解质盐传输、电解质电荷守恒、固相锂传输、固相电荷守恒以及界面 Butler-Volmer kinetics 离散到图结构上。

## 1. 几何、变量与符号约定

### 1.1 Network Sets

模型使用两个耦合图和一个界面集合：

- Electrolyte graph：节点 `i in E`，喉道 `(i,k) in T_e`。
- Solid graph：active-material/CBD 节点 `m in S`，固相接触 `(m,n) in T_s`。
- Electrolyte/active-material interfaces：反应 `r in R`，连接一个 electrolyte 节点 `i(r)` 和一个 active solid 节点 `m(r)`。

几何量：

| Symbol | Unit | 含义 |
|---|---:|---|
| `V_i^e` | `m^3` | electrolyte pore control-volume volume |
| `V_m^s` | `m^3` | active solid control-volume volume |
| `A_ik^e` | `m^2` | electrolyte throat cross-sectional area |
| `L_ik^e` | `m` | electrolyte throat length |
| `A_mn^s` | `m^2` | solid contact cross-sectional area |
| `L_mn^s` | `m` | solid contact length |
| `A_r` | `m^2` | electrolyte/NMC reaction interface area |
| `L` | `m` | cathode thickness，separator 在 `x=0`，collector 在 `x=L` |

### 1.2 Unknown Fields

每个时间层 `t^{n+1}` 的 fully coupled unknown vector 为：

```text
x = [ c_e, phi_e, c_s, phi_s ]^T
```

| Symbol | Unit | 含义 |
|---|---:|---|
| `c_i^e` | `mol m^-3` | electrolyte 节点 `i` 的 Li salt concentration |
| `phi_i^e` | `V` | electrolyte 节点 `i` 的 potential |
| `c_m^s` | `mol m^-3` | active solid 节点 `m` 的 solid-phase lithium concentration |
| `phi_m^s` | `V` | solid 节点 `m` 的 electronic potential |

常数：

| Symbol | Unit | 含义 |
|---|---:|---|
| `F` | `C mol^-1` | Faraday constant |
| `R` | `J mol^-1 K^-1` | gas constant |
| `T` | `K` | temperature |
| `D_e` | `m^2 s^-1` | electrolyte salt diffusion coefficient |
| `D_s` | `m^2 s^-1` | solid lithium diffusion coefficient |
| `kappa` | `S m^-1` | electrolyte ionic conductivity |
| `sigma` | `S m^-1` | solid electronic conductivity |
| `c_s,max` | `mol m^-3` | NMC 最大锂浓度 |
| `U(c_s,c_e)` | `V` | NMC 相对 Li/Li+ 的 equilibrium potential |
| `alpha_a, alpha_c` | `1` | anodic 与 cathodic charge-transfer coefficients |

### 1.3 Reaction Sign Convention

采用标准 anodic convention：

```text
Li_s  <->  Li+_e + e-_s
```

- `q_r` [`mol m^-2 s^-1`]：正值表示 anodic deintercalation/oxidation。
- `i_r = F q_r` [`A m^-2`]：正值表示 anodic Faradaic current。
- cathode discharge/lithiation 时发生 reduction，因此 `q_r < 0` 且 `i_r < 0`。

因此，`q_r > 0` 时 electrolyte Li 生成、solid Li 消耗；`q_r < 0` 时 electrolyte Li 消耗、solid Li 增加。界面反应源项在 electrolyte 与 solid charge balances 中符号相反，总界面电荷严格守恒。

### 1.4 Current and Voltage Sign Convention

令 `I_app` [`A`] 为从 current collector 进入 cathode solid 的 conventional current，`I_app/A_cell` 为施加在 `x=L` 的 current density。

- `I_app > 0`：cathode anodic operation，deintercalation/charge。
- `I_app < 0`：cathode cathodic operation，intercalation/discharge。

许多电池论文把 discharge current 报告为正的 C-rate magnitude。本文中 `I_dis > 0` 对应 `I_app = -I_dis`。

## 2. Continuum Equations from First Principles

### 2.1 Electrolyte Nernst-Planck Transport

对带电数 `z_j` 的 species `j`，dilute-solution Nernst-Planck flux 为：

```text
N_j = -D_j grad c_j - z_j u_j F c_j grad phi_e + c_j v
```

其中 `N_j` 是 molar flux [`mol m^-2 s^-1`]，`u_j` 是 mobility。忽略 convection 并使用 Einstein relation `u_j = D_j/(R T)`，得到：

```text
N_j = -D_j grad c_j - z_j (D_j F/(R T)) c_j grad phi_e.
```

二元 electrolyte 中 electroneutrality 给出 `c_+ = c_- = c_e`，电流密度为：

```text
i_e = F (N_+ - N_-).
```

理想二元 electrolyte 可写为：

```text
i_e = -kappa(c_e) grad phi_e - kappa_D(c_e) grad ln c_e.
```

若忽略 diffusion-potential effects，或近似为 supporting electrolyte / equal ion diffusivities，电荷定律退化为 Ohm's law：

```text
i_e = -kappa(c_e) grad phi_e.
```

简化模型使用该 Ohmic form，并可采用 Nernst-Einstein conductivity：

```text
kappa(c_e) = F^2 D_e c_e/(R T).
```

其单位为 `S m^-1`：

```text
(C^2 mol^-2)(m^2 s^-1)(mol m^-3)/(J mol^-1)
= C^2 s kg^-1 m^-3 = A V^-1 m^-1 = S m^-1.
```

Electrolyte salt conservation：

```text
partial c_e/partial t = -div N_s + S_e.
```

最简单 Fickian salt model 中：

```text
N_s = -D_e,eff grad c_e.
```

界面源项取决于 electrolyte model 的 transference 处理：

```text
S_e = beta_e a_s i_F/F.
```

`a_s` 是 interfacial area per volume [`m^-1`]，`i_F` 是 anodic Faradaic current density [`A m^-2`]。full binary electrolyte model 常取 `beta_e = 1 - t_+^0`；简化 Li+ balance 常取 `beta_e = 1`。因此简化 PNM 的盐方程为：

```text
partial c_e/partial t = div(D_e,eff grad c_e) + beta_e a_s i_F/F.
```

cathode discharge 时 `i_F < 0`，electrolyte Li 被消耗。

### 2.2 Electrolyte Charge Conservation

Electroneutral electrolyte 的 charge accumulation 可忽略：

```text
div i_e = beta_q a_s i_F.
```

简化 single-current model 中 `beta_q = 1`。使用 Ohm's law：

```text
div(-kappa_eff grad phi_e) = a_s i_F.
```

若保留 concentration-potential term，则 edge current 需要包含：

```text
i_e = -kappa_eff grad phi_e - kappa_D,eff grad ln c_e.
```

### 2.3 Solid Lithium Diffusion from Fick's Law

NMC 中的 solid lithium 是 neutral intercalated lithium，其通量为：

```text
N_solid = -D_s(c_s,T) grad c_s.
```

active material 质量守恒：

```text
partial c_s/partial t = div(D_s grad c_s) - a_s i_F/F.
```

`i_F > 0` 时 deintercalation 消耗 solid lithium；discharge 时 `i_F < 0`，solid lithium 增加。若 `D_s` 依赖浓度，continuum operator 是：

```text
div(D_s(c_s,T) grad c_s)
```

而不是 `D_s grad^2 c_s`，除非 `D_s` 局部常数。CBD 传导电子但不储存 intercalated lithium，因此 `c_s` 方程只应在 active NMC storage volumes 上组装，`phi_s` 可在 NMC 与 CBD electronic nodes 上组装。

### 2.4 Solid Electronic Ohm's Law

Solid conventional current density：

```text
i_s = -sigma_eff grad phi_s.
```

Faradaic transfer 下的 charge conservation：

```text
div i_s = -a_s i_F
```

即：

```text
div(-sigma_eff grad phi_s) = -a_s i_F.
```

把 electrolyte 与 solid charge equations 相加：

```text
div(i_e + i_s) = 0.
```

界面反应只在相之间转移电荷，不产生净电荷。

### 2.5 Butler-Volmer Kinetics

对反应：

```text
Li_s <-> Li+_e + e-_s
```

界面 overpotential 为：

```text
eta_r = phi_m^s - phi_i^e - U(c_m^s, c_i^e).
```

`U` 是相同 electrolyte reference 下 active material 相对 Li/Li+ 的 equilibrium potential。anodic Butler-Volmer law 为：

```text
i_r = i0_r [ exp(alpha_a F eta_r/(R T))
             - exp(-alpha_c F eta_r/(R T)) ].
```

因此 `eta_r > 0` 得到 `i_r > 0`（deintercalation/oxidation），`eta_r < 0` 得到 `i_r < 0`（intercalation/reduction）。

常用 exchange-current model：

```text
i0_r = F k0 c_e^gamma_e (c_s,max - c_s)^gamma_v c_s^gamma_s.
```

常见 Li-ion 取值为：

```text
gamma_e = alpha_a
gamma_v = alpha_a
gamma_s = alpha_c
```

实现必须保持单位一致。若没有 empirical OCV curve，可用 ideal Nernst form：

```text
U = U0 + (R T/F) ln[ ((c_s,max - c_s)/c_s) (c_e/c_e,ref) ].
```

对 NMC532，通常优先使用 empirical `U(x)`，其中 `x = c_s/c_s,max`。如果 empirical curve 已相对 Li/Li+ 且在 reference electrolyte concentration 下测得，除非论文明确加入 Nernst correction，否则不要重复计入 electrolyte concentration。

## 3. Pore Network Discretization

有限体积 control volume 是 pore 或 solid node。Throats 提供 two-point flux approximations。若 edge 两侧 half-throat 的材料属性不同，应使用 harmonic 或 series averaging。

### 3.1 Edge Conductances

Electrolyte diffusive conductance：

```text
K_ik^e = D_ik,eff^e A_ik^e/L_ik^e        [m^3 s^-1]
```

Electrolyte ionic conductance：

```text
G_ik^e = kappa_ik,eff A_ik^e/L_ik^e      [S = A V^-1]
```

Solid diffusive conductance：

```text
K_mn^s = D_mn^s A_mn^s/L_mn^s            [m^3 s^-1]
```

Solid electronic conductance：

```text
G_mn^s = sigma_mn,eff A_mn^s/L_mn^s      [S]
```

Bruggeman-type corrections 例如：

```text
D_e,eff = D_e epsilon^b
kappa_eff = kappa epsilon^b
```

其中 prompt 中 `b = 1.5`。审查时应确认该修正用于 pore-scale throats、continuum separator regions 还是两者。若同时用 pore geometry 与 Bruggeman correction，可能 double-count tortuosity。

### 3.2 Electrolyte Concentration Residual

对 electrolyte node `i`，邻居 `k` 流入 `i` 的 diffusive molar flow：

```text
M_ki^e = K_ik^e (c_k^e - c_i^e)          [mol s^-1]
```

因为 `K [m^3 s^-1] * Delta c [mol m^-3] = mol s^-1`。令 `R_i` 是连接到 electrolyte node `i` 的 reaction interfaces 集合。Backward Euler residual：

```text
F_{c_e,i} =
V_i^e (c_i^{e,n+1} - c_i^{e,n})/dt
- sum_k K_ik^e (c_k^{e,n+1} - c_i^{e,n+1})
- sum_{r in R_i} beta_e (i_r/F) A_r = 0.
```

符号检查：discharge 时 `i_r < 0`，最后一项为正的消耗贡献，推动 `c_e` 下降。

### 3.3 Electrolyte Potential Residual

邻居 `k` 流入 electrolyte node `i` 的 conductive current：

```text
I_ki^e = G_ik^e (phi_k^e - phi_i^e)      [A]
```

quasi-static charge balance：

```text
F_{phi_e,i} =
sum_k G_ik^e (phi_k^e - phi_i^e)
+ sum_{r in R_i} i_r A_r = 0.
```

若保留 concentration-potential terms，应给每条 edge current 添加：

```text
I_ki^{e,conc} = H_ik^e (ln c_k^e - ln c_i^e).
```

`H_ik^e` 的单位为 ampere，取决于所选 Nernst-Planck / concentrated-solution model。

### 3.4 Solid Concentration Residual

对 active solid node `m`，邻居 `n` 流入 `m` 的 diffusive molar flow：

```text
M_nm^s = K_mn^s (c_n^s - c_m^s).
```

令 `R_m` 是连接到 active solid node `m` 的 reaction interfaces 集合。Backward Euler residual：

```text
F_{c_s,m} =
V_m^s (c_m^{s,n+1} - c_m^{s,n})/dt
- sum_n K_mn^s (c_n^{s,n+1} - c_m^{s,n+1})
+ sum_{r in R_m} (i_r/F) A_r = 0.
```

discharge 时 `i_r < 0`，该反应项使 `c_s` 增加。非线性 `D_s(c_s,T)` 下，`K_mn^s` 应在 `t^{n+1}` 评估，fully implicit Newton solve 可使用 harmonic averaging。

### 3.5 Solid Potential Residual

邻居 `n` 流入 solid node `m` 的 conventional electronic current：

```text
I_nm^s = G_mn^s (phi_n^s - phi_m^s).
```

quasi-static charge residual：

```text
F_{phi_s,m} =
sum_n G_mn^s (phi_n^s - phi_m^s)
- sum_{r in R_m} i_r A_r
+ F_{BC,m}^s = 0.
```

使用同一 anodic convention 时，`F_{phi_e}` 与 `F_{phi_s}` 的 reaction terms 完全相反。

### 3.6 Generic Coupled Interface Contribution

对 electrolyte node `i` 与 solid node `m` 之间的一个 interface `r`：

```text
F_{c_e,i}     += - beta_e (i_r/F) A_r
F_{phi_e,i}   += + i_r A_r
F_{c_s,m}     += + (i_r/F) A_r
F_{phi_s,m}   += - i_r A_r
```

这四行是模型的核心耦合块。

## 4. Boundary Conditions and Voltage Reference

### 4.1 Separator/Electrolyte Boundary at `x=0`

prompt 中 separator 连接到 reservoir。对 cathode electrolyte graph，可在 separator face 施加：

```text
c_i^e = c_e,0
phi_i^e = 0
```

若显式建模 separator resistance，则可使用等价 ghost-node / Robin condition。若 separator 是厚度 `L_sep`、孔隙率 `epsilon_sep` 的 1D 区域：

```text
D_sep,eff = D_e epsilon_sep^b
kappa_sep,eff = kappa epsilon_sep^b
```

Li metal/reservoir 边界：

```text
c_e = c_e,0
phi_e = 0
```

与 cathode 的 interface 需要 salt flux 与 ionic current 连续。若不显式建模 separator，则忽略其浓度降与 ohmic drop；高倍率验证中这一省略可能显著。

### 4.2 Collector Boundary at `x=L`

cathode current collector 处：

Electrolyte：

```text
N_s . n = 0
i_e . n = 0
```

collector 之外没有 electrolyte，因此无 salt 或 ionic current 通过该边界。

Solid：

```text
integral_{collector} i_s . n dA = I_app.
```

图形式中，可按面积权重 `w_m` 把 `I_app` 分配到 collector-connected solid nodes，并满足 `sum w_m = 1`：

```text
F_{BC,m}^s = I_app w_m
```

符号必须通过检查 discharge input 是否降低 cell voltage 来验证。另一种更稳健的 polarization curve 方法是 pin `phi_s` at collector，求总电流，再用外层 scalar Newton/secant loop 调整 collector potential，直到 `I_total = I_app`。

### 4.3 Solid Boundary at Separator Face

没有 electronic current 穿过 separator：

```text
i_s . n = 0
```

除非有意建模电子短路。Solid lithium 在非反应外边界为 no-flux：

```text
N_solid . n = 0.
```

锂只通过 electrolyte/NMC interfaces 进入或离开 active material。

### 4.4 Potential Gauge and Reference Frame

只有 potential differences 具有物理意义。若没有 pin 一个 potential，所有 `phi_e` 与 `phi_s` 同加常数不会改变方程，系统奇异。半电池模型必须设置参考点。

prompt 的参考为：

```text
phi_e at separator/Li reference = 0 V.
```

对平衡 Li metal counter/reference electrode：

```text
phi_s,Li - phi_e,sep = U_Li/Li+ = 0 V
```

因此 `phi_e,sep = 0` 时 `phi_s,Li = 0`。NMC OCV `U` 也相对该 Li/Li+ reference 测得，所以均匀 cathode equilibrium 为：

```text
phi_s,cathode - phi_e,cathode = U(c_s,c_e).
```

### 4.5 Correct Cell Voltage Formula

可测 half-cell voltage 是 positive current collector potential 减 Li metal potential：

```text
V_cell = phi_s,collector - phi_s,Li.
```

在上述参考下：

```text
V_cell = phi_s,collector.
```

更一般地，若 electrolyte separator potential 不为零：

```text
V_cell = phi_s,collector - phi_e,sep - U_Li/Li+ - eta_Li.
```

ideal Li metal reference 中 `U_Li/Li+ = 0` 且 `eta_Li = 0`，因此：

```text
V_cell = phi_s,collector - phi_e,sep.
```

审查重点：`V_cell = phi_s,collector - phi_e,separator` 只在 Li metal electrode 是 equilibrium Li/Li+ reference 且 electrolyte reference frame 一致时成立；它不是通用 full-cell voltage formula。

## 5. Newton-Raphson Formulation

### 5.1 Residual Vector

每个时间步求解：

```text
F(x^{n+1}; x^n, dt, I_app) = 0
```

块结构：

```text
F =
[ F_{c_e}     ]
[ F_{phi_e}   ]
[ F_{c_s}     ]
[ F_{phi_s}   ].
```

Backward Euler 使两个 concentration equations 隐式；potential equations 是 algebraic quasi-steady constraints。每步是 DAE-like nonlinear algebraic solve。

Newton iteration `ell`：

```text
J(x_ell) delta x_ell = -F(x_ell)
x_{ell+1} = x_ell + lambda delta x_ell
```

其中 `0 < lambda <= 1` 是 damping / line-search parameter。

### 5.2 Butler-Volmer Derivatives

对 interface `r`：

```text
f = F/(R T)
eta = phi_s - phi_e - U(c_s,c_e)
E_a = exp(alpha_a f eta)
E_c = exp(-alpha_c f eta)
i = i0 (E_a - E_c)
B = d i/d eta = i0 f (alpha_a E_a + alpha_c E_c).
```

直接 potential derivatives：

```text
d i/d phi_s = +B
d i/d phi_e = -B
```

concentration derivatives：

```text
d i/d c_s =
(d i0/d c_s)(E_a - E_c) - B (dU/dc_s)

d i/d c_e =
(d i0/d c_e)(E_a - E_c) - B (dU/dc_e)
```

若：

```text
i0 = F k0 c_e^gamma_e (c_s,max - c_s)^gamma_v c_s^gamma_s
```

log derivatives 为：

```text
d i0/d c_e = i0 gamma_e/c_e
d i0/d c_s = i0 [ gamma_s/c_s - gamma_v/(c_s,max - c_s) ].
```

这些导数在 `c_e = 0`、`c_s = 0`、`c_s = c_s,max` 处奇异，实际实现必须 clip 或 regularize。

### 5.3 Jacobian Contributions from One Interface

令 `p` 是 `c_e_i`、`phi_e_i`、`c_s_m`、`phi_s_m` 中任一 local variable，interface contribution：

```text
dF_{c_e,i}/dp   += - beta_e (A_r/F) d i_r/dp
dF_{phi_e,i}/dp += + A_r d i_r/dp
dF_{c_s,m}/dp   += + (A_r/F) d i_r/dp
dF_{phi_s,m}/dp += - A_r d i_r/dp
```

这里完成 `phi_e` 与 `phi_s` 的耦合。对浓度的导数包含 exchange-current 与 OCV derivative。

### 5.4 Jacobian Contributions from Linear Edges

Electrolyte concentration residual：

```text
F_i = V_i(c_i-c_i^n)/dt - sum_k K_ik(c_k-c_i) - sources.
```

若 `K_ik` 常数：

```text
dF_i/dc_i += V_i/dt + sum_k K_ik
dF_i/dc_k += -K_ik.
```

Electrolyte potential residual：

```text
F_i = sum_k G_ik(phi_k-phi_i) + sources.
```

若 `G_ik` 常数：

```text
dF_i/dphi_i += -sum_k G_ik
dF_i/dphi_k += +G_ik.
```

Solid concentration 与 solid potential 采用相同结构：

```text
dF_m/dc_m += V_m/dt + sum_n K_mn
dF_m/dc_n += -K_mn
dF_m/dphi_m += -sum_n G_mn
dF_m/dphi_n += +G_mn.
```

若 `K` 或 `G` 依赖浓度，还需链式法则项。

### 5.5 Dirichlet Boundary Rows

对 pinned unknown `y = y_B`，用边界 residual 替换守恒行：

```text
F_y = y - y_B = 0
dF_y/dy = 1
```

该行其他 Jacobian entries 为零，避免把 conservation rows 与 boundary rows 混合。

### 5.6 Galvanostatic Constraint Options

有两类正确策略：

Option A：在 `F_phi_s` 中施加 Neumann current boundary：

```text
sum_m F_{BC,m}^s = I_app.
```

此时 `phi_s,collector` 是求解输出。

Option B：使用 unknown collector voltage 加 scalar current constraint。把 collector solid nodes pin 到共同 unknown `V_col`，再加入：

```text
F_I = sum_{r in R} i_r A_r - I_app = 0
```

符号需与约定一致；若包含 double-layer/storage terms 也应加入。该方法适合强制 equipotential current collector。

稳态恒流且无 capacitive side reactions 时：

```text
sum_r i_r A_r = I_app.
```

## 6. Numerical Pitfalls and Remedies

### 6.1 Exponential Overflow in Butler-Volmer

`T = 298 K` 时 `R T/F approximately 25.7 mV`，指数项随 overpotential 快速增长。处理方法：

- clip exponent arguments，例如限制在 `[-500, 500]`。
- 对大 overpotential 使用 asymptotic forms。
- Newton line search 中限制 potential step。
- 对 `i0` 设置 concentration floors。

### 6.2 Concentration Bounds

物理有效范围：

```text
c_e > 0
0 < c_s < c_s,max
```

边界附近会出现 `i0` derivative singularity、OCV slope 爆炸或扩散系数不可信。处理方法：

- 对状态变量设置 floor/ceiling。
- step acceptance 检查浓度是否越界。
- 需要时缩小 `dt` 并重试。
- 报告 cutoff/depletion，而不是让求解器继续进入非物理区。

### 6.3 Singular Matrices from Potential Gauge Freedom

若没有 pin electrolyte potential，所有 potentials 加常数不改变方程，Jacobian 有 null vector。应只 pin 一个 reference potential，通常是 separator/Li reference 处 `phi_e = 0`。恒流工况下不要再 pin 不一致的 solid potential，除非把电流作为输出或加入额外约束。

### 6.4 Disconnected Pores or Solids

未连接到 separator Dirichlet condition 的 electrolyte clusters 有 floating `phi_e`；未连接到 collector 的 solid clusters 有 floating `phi_s`。若这些 clusters 还包含 active interfaces，BV 可能产生无物理路径支撑的局部反应。

必要预处理：

- 标记与 separator 连通的 electrolyte component。
- 标记与 collector 连通的 solid component。
- 只在两侧 transport paths 都存在的 interfaces 上允许反应。
- 对 disconnected floating components 施加固定值或从求解系统中移除。

### 6.5 Ill Conditioning from Conductivity Contrast

CBD conductivity (`~760 S m^-1`) 与 NMC conductivity (`~0.01 S m^-1`) 相差多个数量级，concentration blocks 与 potential blocks 的物理单位也不同。建议：

- 使用 sparse solvers 和合理 scaling。
- 对变量/残差做 nondimensionalization。
- 避免把孤立节点留在线性系统中。
- 在 Newton 中使用 damping。

### 6.6 Time Step Selection

Backward Euler 对线性扩散无条件稳定，但非线性精度和 Newton convergence 仍限制 `dt`。有用尺度：

```text
tau_e ~ l_pore^2/D_e,eff
tau_s ~ l_particle^2/D_s
tau_rxn ~ F c_s,max V_active/(|I_app|)
```

实用策略：

- 初始 `dt` 取上述尺度的保守比例。
- Newton 收敛快则增大 `dt`。
- Newton 失败或电压跳变过大则 rollback 并减小 `dt`。
- cutoff voltage 附近限制 voltage jump。

### 6.7 Mass and Charge Conservation Audits

每个收敛时间步检查：

Electrolyte lithium：

```text
Delta n_e =
sum_i V_i^e (c_i^{e,n+1} - c_i^{e,n})
```

应等于 integrated boundary salt flux 加：

```text
dt sum_r beta_e i_r A_r/F.
```

Solid lithium：

```text
Delta n_s =
sum_m V_m^s (c_m^{s,n+1} - c_m^{s,n})
```

应等于：

```text
-dt sum_r i_r A_r/F
```

再加通常为零的 solid diffusive boundary fluxes。

Charge：

```text
sum_r i_r A_r = I_app
```

适用于无 double-layer capacitance 与 side reaction 的 galvanostatic cathode。

### 6.8 Sign Convention Bugs

常见错误：

- discharge 时把 `I_app` 作为正数传入 anodic convention。
- 在 electrolyte 与 solid concentration source 中使用相同符号。
- 把 `eta = phi_e - phi_s - U` 写反。
- 使用 `V_cell = phi_e - phi_s`。
- 在 Li reference 与 separator potential 之间重复减压降。

single-pore equilibrium 与 small-current tests 应捕获这些错误。

## 7. Known-Limit Validation

### 7.1 Zero-Current Equilibrium (`0C`)

设 `I_app = 0`。均匀初始浓度且相连通时，期望：

```text
c_e = c_e,0 everywhere
c_s = c_s,0 everywhere
phi_e = 0 everywhere, if separator electrolyte is pinned to 0
eta_r = 0 at every interface
i_r = 0 at every interface
phi_s = U(c_s,0, c_e,0) everywhere in the connected solid
V_cell = U(c_s,0, c_e,0).
```

若 empirical `U(x)` 不含 electrolyte correction：

```text
V_cell = U(x_0).
```

若初始 `c_s` 空间不均匀，真实平衡需要 uniform electrochemical potential；简化模型通常向 uniform `U` 松弛，而不一定是 uniform `c_s`。

### 7.2 Small-Current / Low-C-Rate Limit

小电流下线性化 BV：

```text
i_r approximately i0_r (F/(R T))(alpha_a + alpha_c) eta_r.
```

若 `alpha_a + alpha_c = 1`：

```text
eta_r approximately (R T/(F i0_r)) i_r.
```

discharge (`I_app < 0`) 下：

```text
V_cell approximately U(x_avg) - |I_app| R_total.
```

`R_total` 包含 reaction、electrolyte、solid、separator 与 contact resistances。相对 OCV 的电压误差应随电流线性缩放。

### 7.3 High-C-Rate Electrolyte Depletion Limit

cathode discharge 时 `i_r < 0`，electrolyte Li 在 cathode 中消耗。1D slab 近似下，separator 在 `x=0`，collector 在 `x=L`，均匀反应，collector no-flux，固定 `c(0)=c0`，令 `I_dis = -I_app > 0`：

```text
s = beta_e I_dis/(F L)       [mol m^-3 s^-1].
```

quasi-steady state：

```text
D_e,eff d^2 c/dx^2 = s
dc/dx(L) = 0
c(0) = c0.
```

解为：

```text
c(x) = c0 + (s/D_e,eff)(x^2/2 - L x).
```

最低浓度在 collector：

```text
c(L) = c0 - s L^2/(2 D_e,eff)
```

近似 diffusion-limited discharge current density：

```text
I_lim approximately 2 F D_e,eff c0/(beta_e L).
```

当 `I_dis` 接近该尺度时，collector 附近 `c_e` 接近零，容量快速下降。

### 7.4 High-C-Rate Ohmic Limit

忽略浓度梯度和反应非均匀性，slab 中 distributed current 的 half-cell polarization 近似为：

```text
Delta V_ohm approximately I_dis
    ( L/(2 kappa_eff) + L/(2 sigma_eff) )
```

discharge 时 terminal voltage 低于 OCV：

```text
V_cell approximately U - Delta V_ohm - |eta_ct| - Delta V_conc.
```

`1/2` 因子来自 ionic current 在 separator 附近最大而 collector 处为零、solid current 在 collector 附近最大而 separator 附近为零。强反应定位会破坏该估计的数值大小，但不改变符号。

### 7.5 Solid Diffusion Limit

高倍率 discharge 时，lithium 进入 particle surface 的速度快于向内扩散。特征时间：

```text
tau_s ~ R_p^2/D_s.
```

若 discharge time `tau_dis` 远小于 `tau_s`，surface concentration 向 `c_s,max` 上升，而 core 仍较低锂化。后果：

- surface OCV 降低，voltage cutoff 提前。
- capacity 低于均匀粒子模型。
- 较小颗粒因 `R_p^2` 尺度更快均匀化。

单粒子验证应再现实心扩散的短时 scaling。constant flux into a semi-infinite solid 时：

```text
Delta c_s,surf = 2 J_in sqrt(t/(pi D_s))
```

其中 `J_in = -i_F/F > 0` 是进入 solid 的 molar flux。

## 8. Ideal Module Architecture

代码应显式表达 sign conventions 与 units。建议架构：

```text
src/
  pnmcathode/
    network/
      generator.py
    physics/
      electrolyte.py
      solid.py
      reaction.py
      ocv.py
      separator.py
    solver/
      steady.py
      transient.py
      single_pore.py
    post/
      analysis.py
      visualization.py
```

推荐职责：

- `network.generator`：创建 network topology、phase labels、geometric volumes/areas/lengths。
- `physics.electrolyte`：提供 `D_e(c,T)`、`kappa(c,T)` 和有效传输系数。
- `physics.solid`：提供 `D_s(c,T)`、particle discretization 与 diffusion RHS。
- `physics.reaction`：集中实现 Butler-Volmer、exchange current、constants 与 sign convention。
- `physics.ocv`：提供 NMC532 `U(x)` 与 derivative。
- `physics.separator`：提供 collapsed 1D separator boundary model。
- `solver.steady`：组装并求解 quasi-static `phi_e` / `phi_s`。
- `solver.transient`：时间推进 `c_e` / `c_s`，调用 steady solver。
- `post`：仅做分析与可视化，不改变仿真状态。

架构规则：reaction module 应只暴露一种 sign convention，所有 residual assemblers 都使用它。不要在多个位置复制 BV formulas。

## 9. Test Strategy

### 9.1 Unit Tests

- `butler_volmer(i0, 0) = 0`。
- `eta > 0` 时 `i > 0`，`eta < 0` 时 `i < 0`。
- `exchange_current_density` 对合法浓度为正，边界附近有限。
- `nmc532_ocv(x)` 在校准范围内单调/有界。
- `electrolyte_diffusion_coefficient` 与 `electrolyte_ionic_conductivity` 单位量级正确。
- `solid_diffusion_rhs` 质量守恒。
- `separator_boundary` 在 disabled 时返回零损失。

### 9.2 Small Network Tests

- 0C 均匀状态下无反应、无浓度变化。
- small-current discharge 的 voltage drop 与 `I_app` 近似线性。
- 反向 charge/discharge 的 sign 对称性。
- disconnected electrolyte / solid components 不导致 singular matrix。
- collector current balance 满足 `sum_r i_r A_r = I_app`。

### 9.3 Continuum-Limit Tests

- 1D slab 中 electrolyte depletion profile 接近解析二次曲线。
- high-rate ohmic polarization 符号正确且量级接近估计。
- single-particle diffusion limit 的 surface concentration 服从 `sqrt(t)` scaling。

## 10. Implementation Audit Checklist

审查实现时使用以下清单：

- [ ] `I_app` 的 discharge sign 与 anodic convention 一致。
- [ ] `eta = phi_s - phi_e - U`。
- [ ] electrolyte 与 solid reaction source 符号相反。
- [ ] `V_cell` 使用 Li/Li+ reference frame 正确计算。
- [ ] 至少一个且只有一个必要的 potential gauge 被 pin。
- [ ] disconnected phases 被移除、固定或禁用反应。
- [ ] `c_e` 与 `c_s` 有合理 bounds / floors。
- [ ] BV exponent 被 clip 或稳定化。
- [ ] `D_s(c_s,T)` 与 `kappa(c_e,T)` 的单位正确。
- [ ] Bruggeman correction 未与 pore geometry 重复计入。
- [ ] collector current balance 与 reaction current balance 通过测试。
- [ ] transient step 在失败时 rollback 并减小 `dt`。
- [ ] post-processing 不改变 solver state。
