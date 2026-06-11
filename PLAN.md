# P2D+PNM 耦合实现计划

本文档给出在 `pnmcathode` 项目中实现 P2D（Pseudo-two-dimensional, Doyle-Fuller-Newman/DFN）与现有 PNM（Pore Network Model）耦合的数学推导、建模边界、耦合路线和分阶段实现计划。本文是规划文档，不要求立即修改现有代码。

## 0. 符号约定与建模范围

当前项目主要面向 `Li metal | separator | NMC532 cathode` 半电池，阴极内部使用三相 PNM：electrolyte、NMC active material、CBD。P2D 模型通常用于 full cell：negative electrode、separator、positive electrode。为兼容当前项目，本文采用如下层次：

- 第一目标：实现 `separator | porous cathode` 的 half-cell DFN，边界端接 Li/Li+ reference 或简化 Li foil BV。
- 第二目标：保留 full-cell DFN 的通用结构，使 negative electrode 后续可作为同类 `ElectrodeRegion` 加入。
- 第三目标：在 cathode porous electrode 中用 PNM 替代或修正 Bruggeman 均匀化闭合关系。

本文沿用现有 `DERIVATION.md` 的 anodic convention：

- 界面反应为 `Li_s <-> Li+_e + e-_s`。
- \(i_F\) 为 Faradaic current density，单位 `A m^-2`，正值表示 anodic deintercalation/oxidation。
- 阴极放电/嵌锂时 \(i_F < 0\)。
- 项目中的外加电流密度 \(I_\mathrm{app}\) 单位 `A m^-2`，正值表示 cathode anodic operation，负值表示 cathode discharge。

### 0.1 全局变量表

| 符号 | 单位 | 含义 |
|---|---:|---|
| \(x\) | `m` | 电极厚度方向坐标 |
| \(r\) | `m` | 活性颗粒径向坐标 |
| \(t\) | `s` | 时间 |
| \(L_n,L_s,L_p\) | `m` | 负极、隔膜、正极厚度 |
| \(R_{p,k}\) | `m` | 第 \(k\) 个电极区域的代表性颗粒半径，\(k \in \{n,p\}\) |
| \(c_e(x,t)\) | `mol m^-3` | electrolyte salt concentration |
| \(c_{s,k}(x,r,t)\) | `mol m^-3` | 第 \(k\) 个固相颗粒内锂浓度 |
| \(c_{s,k}^{surf}(x,t)\) | `mol m^-3` | 固相颗粒表面锂浓度 \(c_{s,k}(x,R_{p,k},t)\) |
| \(\bar c_{s,k}(x,t)\) | `mol m^-3` | 颗粒体积平均固相锂浓度 |
| \(\phi_e(x,t)\) | `V` | electrolyte potential |
| \(\phi_{s,k}(x,t)\) | `V` | solid-phase electronic potential |
| \(i_e(x,t)\) | `A m^-2` | electrolyte ionic current density，沿 \(x\) 方向 |
| \(i_{s,k}(x,t)\) | `A m^-2` | solid electronic current density，沿 \(x\) 方向 |
| \(i_{F,k}(x,t)\) | `A m^-2` | electrolyte/solid 界面 Faradaic current density，按界面面积计 |
| \(a_{s,k}\) | `m^2 m^-3 = m^-1` | specific interfacial area |
| \(\epsilon_{e,k}\) | `1` | electrolyte volume fraction |
| \(\epsilon_{s,k}\) | `1` | active solid volume fraction |
| \(D_e\), \(D_e^\mathrm{eff}\) | `m^2 s^-1` | electrolyte diffusivity 与有效 diffusivity |
| \(D_{s,k}\) | `m^2 s^-1` | solid diffusivity |
| \(\kappa\), \(\kappa^\mathrm{eff}\) | `S m^-1` | electrolyte ionic conductivity 与有效 ionic conductivity |
| \(\sigma_k\), \(\sigma_k^\mathrm{eff}\) | `S m^-1` | solid electronic conductivity 与有效 electronic conductivity |
| \(t_+^0\) | `1` | Li+ transference number relative to solvent velocity |
| \(f_\pm\) | `1` | mean molar activity coefficient |
| \(U_k\) | `V` | open-circuit potential vs Li/Li+ |
| \(F\) | `C mol^-1` | Faraday constant |
| \(R\) | `J mol^-1 K^-1` | gas constant |
| \(T\) | `K` | temperature |

## 1. 问题定义

### 1.1 P2D 模型建模什么

P2D/DFN 是沿电池厚度方向的一维宏观 porous electrode model，但在每个宏观位置 \(x\) 内嵌一个球形活性颗粒径向扩散问题，因此称为 pseudo-two-dimensional。其核心是：

- 宏观维度 \(x\)：描述 electrolyte concentration、electrolyte potential、solid potential 沿电极厚度和隔膜的变化。
- 微观颗粒维度 \(r\)：描述活性颗粒内部 solid lithium diffusion。
- 界面反应：用 Butler-Volmer kinetics 将 \(c_e\)、\(c_s^{surf}\)、\(\phi_e\)、\(\phi_s\) 耦合。
- 多孔结构闭合：通过 \(\epsilon\)、\(a_s\)、Bruggeman/tortuosity 关系给出 \(D_e^\mathrm{eff}\)、\(\kappa^\mathrm{eff}\)、\(\sigma^\mathrm{eff}\)。

P2D 解决的问题是电池厚度方向的宏观极化、电解液耗竭、固相浓度梯度、倍率性能、容量利用率和电压曲线预测。它是工业中最常用的 physics-based cell model。

### 1.2 PNM 模型建模什么

现有 PNM 将阴极显微结构离散为图：

- electrolyte graph：孔隙节点和 electrolyte throats，承载 \(c_e\)、\(\phi_e\)。
- solid graph：NMC/CBD 节点和 solid throats，承载 \(c_s\)、\(\phi_s\)。
- electrolyte/NMC interfaces：反应界面，承载 Butler-Volmer 反应。

PNM 解决的问题是微观结构导致的局部 transport bottleneck、percolation、反应不均匀性、CBD/NMC 电子通路、电解质孔隙连通性、局部过电位分布和局部利用率。

### 1.3 为什么要耦合

P2D 的局限：

- Bruggeman 均匀化无法表达真实 XCT/合成网络的连通性、死孔、瓶颈喉道、局部界面面积分布。
- \(a_s\)、\(\tau\)、\(D^\mathrm{eff}\)、\(\kappa^\mathrm{eff}\)、\(\sigma^\mathrm{eff}\) 通常被当作标量，无法捕捉微结构沿厚度方向或局部区域的差异。
- 对高倍率、厚电极、低孔隙率、CBD 断连等情形，经验有效参数可能失效。

PNM 的局限：

- 现有模型主要是一个代表性阴极网络，没有完整宏观厚度方向的 DFN cell-level coupling。
- 隔膜当前是 folded 1D boundary，不是与阴极连续求解的完整 P2D separator domain。
- 如果只求一个 network，难以表达电池尺度上 \(x\) 方向的 electrolyte depletion 和 current redistribution。
- 计算成本高，直接将整个电极厚度的所有微结构完全解析为 PNM 不适合参数扫描。

耦合后的目标能力：

- 用 PNM 从显微结构预测 P2D 所需有效参数，减少 Bruggeman 经验假设。
- 在 P2D 宏观网格中嵌入 PNM 子问题，表达局部微结构异质性对反应分布和极化的反馈。
- 支持 half-cell 和最终 full-cell 的电压曲线、浓度场、过电位分解、局部热点和结构敏感性分析。
- 在保留 PyBaMM/Newman 级别 DFN 可验证性的同时，加入项目已有 PNM 的微结构分辨能力。

## 2. P2D 模型公式推导

### 2.1 区域与体积分数

定义 full-cell 厚度域为 \(x \in [0,L]\)，其中 \(L=L_n+L_s+L_p\)。负极为 \(\Omega_n=[0,L_n]\)，隔膜为 \(\Omega_s=[L_n,L_n+L_s]\)，正极为 \(\Omega_p=[L_n+L_s,L]\)。当前 half-cell 可只保留 \(\Omega_s\) 与 \(\Omega_p\)，并在 \(\Omega_s\) 左端接 Li metal reference。

在电极区域 \(k \in \{n,p\}\) 内，\(\epsilon_{e,k}\)、\(\epsilon_{s,k}\)、\(\epsilon_{f,k}\) 分别为 electrolyte、active solid、filler/CBD 体积分数，并满足 \(\epsilon_{e,k}+\epsilon_{s,k}+\epsilon_{f,k}=1\)。隔膜区域只有 electrolyte-filled porous phase，使用 \(\epsilon_{e,s}\)。

常用球形颗粒 specific area 为 \(a_{s,k}=3\epsilon_{s,k}/R_{p,k}\)。如果用 PNM 提供界面面积，则应替换为 \(a_{s,k}=A_{\mathrm{e/NMC},k}/V_{\mathrm{RVE},k}\)。

| 方程 | 变量定义与单位 |
|---|---|
| \(a_{s,k}=3\epsilon_{s,k}/R_{p,k}\) | \(a_{s,k}\) [`m^-1`]；\(\epsilon_{s,k}\) [`1`]；\(R_{p,k}\) [`m`] |
| \(a_{s,k}=A_{\mathrm{e/NMC},k}/V_{\mathrm{RVE},k}\) | \(A_{\mathrm{e/NMC},k}\) [`m^2`]；\(V_{\mathrm{RVE},k}\) [`m^3`] |

### 2.2 固相扩散：Fick 定律到球坐标 P2D 颗粒方程

固相中的 intercalated lithium 近似为中性 species，其 Fick 通量为 \(N_{s,k}=-D_{s,k}(c_s,T)\nabla c_{s,k}\)，单位 `mol m^-2 s^-1`。质量守恒为 \(\partial c_{s,k}/\partial t=-\nabla\cdot N_{s,k}\)，因此得到 \(\partial c_{s,k}/\partial t=\nabla\cdot(D_{s,k}\nabla c_{s,k})\)。

P2D 假设每个宏观位置 \(x\) 内的活性颗粒为球形，颗粒内部只存在径向梯度，得到固相扩散控制方程：\(\partial c_{s,k}/\partial t=(1/r^2)\partial/\partial r(r^2D_{s,k}\partial c_{s,k}/\partial r)\)，适用于 \(x\in\Omega_k\)、\(0<r<R_{p,k}\)。

球心对称边界为 \(\partial c_{s,k}/\partial r|_{r=0}=0\)。颗粒表面通量由 Faradaic reaction 决定：\(-D_{s,k}\partial c_{s,k}/\partial r|_{r=R_{p,k}}=i_{F,k}/F\)。在本文符号下，\(i_{F,k}>0\) 表示锂从 solid 离开进入 electrolyte，因此 outward solid flux 为正。

颗粒体积平均浓度为 \(\bar c_{s,k}=3R_{p,k}^{-3}\int_0^{R_{p,k}}c_{s,k}r^2dr\)。由表面积分可得平均浓度方程：\(\epsilon_{s,k}\partial \bar c_{s,k}/\partial t= -a_{s,k}i_{F,k}/F\)。

| 方程 | 变量定义与单位 |
|---|---|
| \(N_{s,k}=-D_{s,k}\nabla c_{s,k}\) | \(N_{s,k}\) [`mol m^-2 s^-1`]；\(D_{s,k}\) [`m^2 s^-1`]；\(\nabla c_{s,k}\) [`mol m^-4`] |
| \(\partial c_{s,k}/\partial t=\nabla\cdot(D_{s,k}\nabla c_{s,k})\) | \(\partial c_s/\partial t\) [`mol m^-3 s^-1`]；右端 [`mol m^-3 s^-1`] |
| \(\partial c_{s,k}/\partial t=(1/r^2)\partial/\partial r(r^2D_{s,k}\partial c_{s,k}/\partial r)\) | \(r\) [`m`]；\(c_{s,k}\) [`mol m^-3`] |
| \(\partial c_{s,k}/\partial r|_{r=0}=0\) | symmetry boundary condition；梯度单位 [`mol m^-4`] |
| \(-D_{s,k}\partial c_{s,k}/\partial r|_{R_{p,k}}=i_{F,k}/F\) | 左端与右端均 [`mol m^-2 s^-1`]；\(i_F\) [`A m^-2`]；\(F\) [`C mol^-1`] |
| \(\bar c_{s,k}=3R_{p,k}^{-3}\int_0^{R_{p,k}}c_{s,k}r^2dr\) | \(\bar c_s\) [`mol m^-3`] |
| \(\epsilon_{s,k}\partial \bar c_{s,k}/\partial t= -a_{s,k}i_{F,k}/F\) | 左右两端均 [`mol m^-3 s^-1`] |

### 2.3 液相传输：Nernst-Planck 到 concentrated-solution DFN

对 species \(j\)，稀溶液 Nernst-Planck 通量为 \(N_j=-D_j\nabla c_j-z_ju_jFc_j\nabla\phi_e+c_jv\)。忽略 convection 时 \(v=0\)，使用 Einstein relation \(u_j=D_j/(RT)\)，得到 \(N_j=-D_j\nabla c_j-z_jD_jFc_j\nabla\phi_e/(RT)\)。

Li-ion electrolyte 在 DFN 中通常采用 concentrated-solution theory。宏观 salt conservation 写为 \(\partial(\epsilon_e c_e)/\partial t = \partial/\partial x(D_e^\mathrm{eff}\partial c_e/\partial x)+(1-t_+^0)a_s i_F/F\)。该方程表示 electrolyte salt 因扩散重分布，并因界面反应与迁移数不等于 1 产生源项。

有效扩散系数在传统 P2D 中用 Bruggeman 闭合：\(D_e^\mathrm{eff}=D_e(c_e,T)\epsilon_e^{b_D}\)。PNM 耦合后，\(D_e^\mathrm{eff}\) 可由网络上扩散 cell problem 给出。

electrolyte ionic current 在 concentrated-solution theory 中为 \(i_e=-\kappa^\mathrm{eff}\partial\phi_e/\partial x +(2RT\kappa^\mathrm{eff}/F)(1-t_+^0)\chi(c_e)\partial\ln c_e/\partial x\)，其中 \(\chi(c_e)=1+\partial\ln f_\pm/\partial\ln c_e\)。若采用理想溶液，\(\chi=1\)。若采用当前 PNM 的简化 Ohmic form，则忽略第二项，得到 \(i_e=-\kappa^\mathrm{eff}\partial\phi_e/\partial x\)。

电解质电荷守恒忽略 double-layer accumulation，给出 \(\partial i_e/\partial x=a_s i_F\)。正的 anodic \(i_F\) 将正离子电流源注入 electrolyte。

| 方程 | 变量定义与单位 |
|---|---|
| \(N_j=-D_j\nabla c_j-z_ju_jFc_j\nabla\phi_e+c_jv\) | \(N_j\) [`mol m^-2 s^-1`]；\(D_j\) [`m^2 s^-1`]；\(z_j\) [`1`]；\(u_j\) [`mol s kg^-1`]；\(v\) [`m s^-1`] |
| \(u_j=D_j/(RT)\) | \(R\) [`J mol^-1 K^-1`]；\(T\) [`K`] |
| \(\partial(\epsilon_e c_e)/\partial t = \partial/\partial x(D_e^\mathrm{eff}\partial c_e/\partial x)+(1-t_+^0)a_s i_F/F\) | \(\epsilon_e\) [`1`]；\(D_e^\mathrm{eff}\) [`m^2 s^-1`]；\(t_+^0\) [`1`]；\(a_s\) [`m^-1`]；all terms [`mol m^-3 s^-1`] |
| \(D_e^\mathrm{eff}=D_e(c_e,T)\epsilon_e^{b_D}\) | \(b_D\) [`1`]；\(D_e^\mathrm{eff}\) [`m^2 s^-1`] |
| \(i_e=-\kappa^\mathrm{eff}\partial\phi_e/\partial x +(2RT\kappa^\mathrm{eff}/F)(1-t_+^0)\chi\partial\ln c_e/\partial x\) | \(i_e\) [`A m^-2`]；\(\kappa^\mathrm{eff}\) [`S m^-1`]；\(\partial\phi_e/\partial x\) [`V m^-1`]；\(\chi\) [`1`] |
| \(\chi(c_e)=1+\partial\ln f_\pm/\partial\ln c_e\) | \(f_\pm\) [`1`]；\(\chi\) [`1`] |
| \(i_e=-\kappa^\mathrm{eff}\partial\phi_e/\partial x\) | simplified Ohmic electrolyte current |
| \(\partial i_e/\partial x=a_s i_F\) | \(\partial i_e/\partial x\) [`A m^-3`]；\(a_s i_F\) [`A m^-3`] |

### 2.4 固相电荷守恒

固相电子电流服从 Ohm 定律：\(i_{s,k}=-\sigma_k^\mathrm{eff}\partial\phi_{s,k}/\partial x\)。由于界面反应将电子从 solid phase 转移到 reaction interface，固相电荷守恒为 \(\partial i_{s,k}/\partial x=-a_{s,k}i_{F,k}\)。

将 electrolyte 与 solid charge equations 相加，得到 \(\partial(i_e+i_s)/\partial x=0\)，因此总电流在同一截面守恒。

传统 Bruggeman 闭合为 \(\sigma_k^\mathrm{eff}=\sigma_k\epsilon_{s,k}^{b_\sigma}\)。在含 CBD 的正极中，\(\sigma^\mathrm{eff}\) 不应只由 NMC 体积分数决定，而应由 NMC/CBD solid graph 的电子 percolation 和 contact conductance 决定。

| 方程 | 变量定义与单位 |
|---|---|
| \(i_{s,k}=-\sigma_k^\mathrm{eff}\partial\phi_{s,k}/\partial x\) | \(i_s\) [`A m^-2`]；\(\sigma^\mathrm{eff}\) [`S m^-1`]；\(\partial\phi_s/\partial x\) [`V m^-1`] |
| \(\partial i_{s,k}/\partial x=-a_{s,k}i_{F,k}\) | \(\partial i_s/\partial x\) [`A m^-3`]；\(-a_s i_F\) [`A m^-3`] |
| \(\partial(i_e+i_s)/\partial x=0\) | total current conservation |
| \(\sigma_k^\mathrm{eff}=\sigma_k\epsilon_{s,k}^{b_\sigma}\) | \(\sigma_k\) [`S m^-1`]；\(b_\sigma\) [`1`] |

### 2.5 Butler-Volmer 反应动力学

界面 overpotential 定义为 \(\eta_k=\phi_{s,k}-\phi_e-U_k(c_{s,k}^{surf},c_e,T)\)。对反应 `Li_s <-> Li+_e + e-_s`，anodic Butler-Volmer law 为 \(i_{F,k}=i_{0,k}\{\exp(\alpha_{a,k}F\eta_k/(RT))-\exp(-\alpha_{c,k}F\eta_k/(RT))\}\)。

交换电流密度可写为 \(i_{0,k}=Fk_{0,k}c_e^{\gamma_e}(c_{s,k}^{max}-c_{s,k}^{surf})^{\gamma_v}(c_{s,k}^{surf})^{\gamma_s}\)。常用选择为 \(\gamma_e=\alpha_a\)、\(\gamma_v=\alpha_a\)、\(\gamma_s=\alpha_c\)。当使用 NMC empirical OCV \(U(x)\) 且该曲线已经相对 Li/Li+ reference 测得时，不应重复加入 electrolyte Nernst correction，除非材料模型明确要求。

线性化时需要导数。令 \(f=F/(RT)\)、\(E_a=\exp(\alpha_a f\eta)\)、\(E_c=\exp(-\alpha_c f\eta)\)，则 \(B=\partial i_F/\partial\eta=i_0f(\alpha_aE_a+\alpha_cE_c)\)。因此 \(\partial i_F/\partial\phi_s=B\)、\(\partial i_F/\partial\phi_e=-B\)。浓度导数为 \(\partial i_F/\partial c_s=(\partial i_0/\partial c_s)(E_a-E_c)-B\partial U/\partial c_s\)，\(\partial i_F/\partial c_e=(\partial i_0/\partial c_e)(E_a-E_c)-B\partial U/\partial c_e\)。

| 方程 | 变量定义与单位 |
|---|---|
| \(\eta_k=\phi_{s,k}-\phi_e-U_k(c_{s,k}^{surf},c_e,T)\) | \(\eta\) [`V`]；\(\phi_s,\phi_e,U\) [`V`] |
| \(i_{F,k}=i_{0,k}[\exp(\alpha_{a,k}F\eta_k/(RT))-\exp(-\alpha_{c,k}F\eta_k/(RT))]\) | \(i_F,i_0\) [`A m^-2`]；\(\alpha_a,\alpha_c\) [`1`] |
| \(i_{0,k}=Fk_{0,k}c_e^{\gamma_e}(c_{s,k}^{max}-c_{s,k}^{surf})^{\gamma_v}(c_{s,k}^{surf})^{\gamma_s}\) | \(k_0\) units depend on exponents；for \(\alpha_a=\alpha_c=0.5\), \(k_0\) [`m^2.5 mol^-0.5 s^-1`] |
| \(B=\partial i_F/\partial\eta=i_0f(\alpha_aE_a+\alpha_cE_c)\) | \(B\) [`A m^-2 V^-1`]；\(f\) [`V^-1`] |
| \(\partial i_F/\partial\phi_s=B\), \(\partial i_F/\partial\phi_e=-B\) | potential derivatives [`A m^-2 V^-1`] |
| \(\partial i_F/\partial c_s=(\partial i_0/\partial c_s)(E_a-E_c)-B\partial U/\partial c_s\) | concentration derivative [`A m mol^-1 s^-?`] as residual Jacobian coefficient; practically [`A m^-2` per `mol m^-3`] |
| \(\partial i_F/\partial c_e=(\partial i_0/\partial c_e)(E_a-E_c)-B\partial U/\partial c_e\) | same units as above |

### 2.6 P2D 边界条件

#### 2.6.1 Full-cell current collectors

在 negative current collector 与 positive current collector 处，electrolyte 不穿过金属集流体，因此 \(N_e\cdot n=0\) 且 \(i_e\cdot n=0\)。solid phase 承载全部外电流，边界条件为 \(i_s\cdot n=I_\mathrm{app}\)，符号应与项目约定统一。由于电势只有差值有意义，需固定一个 gauge，例如 \(\phi_s(0,t)=0\) 或 \(\phi_e\) 在 Li reference 处为 0。

| 方程 | 变量定义与单位 |
|---|---|
| \(N_e\cdot n=0\) | no electrolyte salt flux；\(N_e\) [`mol m^-2 s^-1`] |
| \(i_e\cdot n=0\) | no ionic current through metal；\(i_e\) [`A m^-2`] |
| \(i_s\cdot n=I_\mathrm{app}\) | \(I_\mathrm{app}\) [`A m^-2`] |

#### 2.6.2 Electrode/separator interfaces

在 electrode/separator 界面，electrolyte concentration 与 potential 连续，并且 salt flux 与 ionic current 连续：\(c_e^- = c_e^+\)、\(\phi_e^-=\phi_e^+\)、\(N_e^-\cdot n=N_e^+\cdot n\)、\(i_e^-\cdot n=i_e^+\cdot n\)。separator 中没有 solid phase，因此 solid current 在电极靠 separator 的边界为 \(i_s\cdot n=0\)。

| 方程 | 变量定义与单位 |
|---|---|
| \(c_e^- = c_e^+\) | electrolyte concentration continuity [`mol m^-3`] |
| \(\phi_e^-=\phi_e^+\) | electrolyte potential continuity [`V`] |
| \(N_e^-\cdot n=N_e^+\cdot n\) | salt flux continuity [`mol m^-2 s^-1`] |
| \(i_e^-\cdot n=i_e^+\cdot n\) | ionic current continuity [`A m^-2`] |
| \(i_s\cdot n=0\) | no electronic current through separator |

#### 2.6.3 Half-cell Li metal boundary

对当前项目更直接的 half-cell 边界，可在 separator/Li metal 端设置 \(\phi_e=0\) 作为 Li/Li+ reference，并可设置 \(c_e=c_{e,ref}\) 或使用 1D separator transport 给出通量边界。若包含 Li foil kinetics，则使用 \(i_\mathrm{Li}=i_{0,\mathrm{Li}}[\exp(\alpha_{a,\mathrm{Li}}F\eta_\mathrm{Li}/RT)-\exp(-\alpha_{c,\mathrm{Li}}F\eta_\mathrm{Li}/RT)]\)，其中 \(\eta_\mathrm{Li}=\phi_{s,\mathrm{Li}}-\phi_e-U_{\mathrm{Li/Li+}}\)，且 \(U_{\mathrm{Li/Li+}}=0\)。

| 方程 | 变量定义与单位 |
|---|---|
| \(\phi_e(0,t)=0\) | electrolyte reference potential [`V`] |
| \(c_e(0,t)=c_{e,ref}\) | reservoir concentration [`mol m^-3`] |
| \(i_\mathrm{Li}=i_{0,\mathrm{Li}}[\exp(\alpha_{a,\mathrm{Li}}F\eta_\mathrm{Li}/RT)-\exp(-\alpha_{c,\mathrm{Li}}F\eta_\mathrm{Li}/RT)]\) | \(i_\mathrm{Li},i_0\) [`A m^-2`] |
| \(\eta_\mathrm{Li}=\phi_{s,\mathrm{Li}}-\phi_e-U_{\mathrm{Li/Li+}}\) | \(\eta_\mathrm{Li}\) [`V`]；\(U_{\mathrm{Li/Li+}}=0\) [`V`] |

#### 2.6.4 Cell voltage

full-cell voltage 为 \(V_\mathrm{cell}=\phi_{s,p}(L,t)-\phi_{s,n}(0,t)\)。half-cell 且 Li metal 为 reference 时，若 \(\phi_{s,\mathrm{Li}}=0\)，则 \(V_\mathrm{cell}=\phi_{s,p}(L,t)\)。如果 separator/Li foil overpotential 被折叠到 \(\phi_e\) 边界，则需要明确记录 \(V_\mathrm{cell}=\phi_{s,p}(L,t)-\phi_{s,\mathrm{Li}}(0,t)\)，避免将 electrolyte gauge 当作真实电压。

| 方程 | 变量定义与单位 |
|---|---|
| \(V_\mathrm{cell}=\phi_{s,p}(L,t)-\phi_{s,n}(0,t)\) | \(V_\mathrm{cell}\) [`V`] |
| \(V_\mathrm{cell}=\phi_{s,p}(L,t)\) | half-cell Li reference 特例 |

## 3. PNM 模型公式与 P2D 对应关系

### 3.1 网络上的质量守恒

对 electrolyte node \(i\)，控制体体积为 \(V_i^e\)，邻居为 \(k\)，反应界面集合为 \(R_i\)。Backward Euler residual 为 \(F_{c_e,i}=V_i^e(c_i^{e,n+1}-c_i^{e,n})/\Delta t-\sum_kK_{ik}^e(c_k^{e,n+1}-c_i^{e,n+1})-\sum_{r\in R_i}\beta_e(i_r/F)A_r=0\)。

对 NMC solid node \(m\)，控制体体积为 \(V_m^s\)，邻居为 \(n\)，反应界面集合为 \(R_m\)。Residual 为 \(F_{c_s,m}=V_m^s(c_m^{s,n+1}-c_m^{s,n})/\Delta t-\sum_nK_{mn}^s(c_n^{s,n+1}-c_m^{s,n+1})+\sum_{r\in R_m}(i_r/F)A_r=0\)。

| 方程 | 变量定义与单位 |
|---|---|
| \(F_{c_e,i}=V_i^e(c_i^{e,n+1}-c_i^{e,n})/\Delta t-\sum_kK_{ik}^e(c_k^{e,n+1}-c_i^{e,n+1})-\sum_{r\in R_i}\beta_e(i_r/F)A_r=0\) | \(V_i^e\) [`m^3`]；\(K_{ik}^e\) [`m^3 s^-1`]；\(A_r\) [`m^2`]；each term [`mol s^-1`] |
| \(F_{c_s,m}=V_m^s(c_m^{s,n+1}-c_m^{s,n})/\Delta t-\sum_nK_{mn}^s(c_n^{s,n+1}-c_m^{s,n+1})+\sum_{r\in R_m}(i_r/F)A_r=0\) | \(V_m^s\) [`m^3`]；\(K_{mn}^s\) [`m^3 s^-1`]；each term [`mol s^-1`] |

### 3.2 喉道通量模型

PNM 的喉道通量通常由 two-point flux approximation 给出。对 electrolyte diffusion，\(M_{ik}^e=K_{ik}^e(c_i^e-c_k^e)\)，其中 \(K_{ik}^e=D_{ik}^{e,\mathrm{eff}}A_{ik}^e/L_{ik}^e\)。对 electrolyte migration/Ohmic conduction，\(I_{ik}^e=G_{ik}^e(\phi_i^e-\phi_k^e)\)，其中 \(G_{ik}^e=\kappa_{ik}^{e,\mathrm{eff}}A_{ik}^e/L_{ik}^e\)。若保留 Nernst-Planck migration 的浓度项，可增加 \(I_{ik}^{e,conc}=H_{ik}^e(\ln c_i^e-\ln c_k^e)\)。

对 solid diffusion，\(M_{mn}^s=K_{mn}^s(c_m^s-c_n^s)\)，其中 \(K_{mn}^s=D_{mn}^sA_{mn}^s/L_{mn}^s\)。对 solid electronic conduction，\(I_{mn}^s=G_{mn}^s(\phi_m^s-\phi_n^s)\)，其中 \(G_{mn}^s=\sigma_{mn}^\mathrm{eff}A_{mn}^s/L_{mn}^s\)。异质材料接触应优先使用 series/harmonic conductance。

Hagen-Poiseuille 形式可用于 pressure-driven electrolyte advection：\(Q_{ik}=g_{ik}^h(p_i-p_k)\)，圆柱喉道 \(g_{ik}^h=\pi r_{ik}^4/(8\mu L_{ik})\)。当前电池模型通常忽略 convection，因此 HP conductance 更适合作为未来 electrolyte flow 或 wetting 扩展。

| 方程 | 变量定义与单位 |
|---|---|
| \(M_{ik}^e=K_{ik}^e(c_i^e-c_k^e)\) | \(M^e\) [`mol s^-1`]；\(K^e\) [`m^3 s^-1`] |
| \(K_{ik}^e=D_{ik}^{e,\mathrm{eff}}A_{ik}^e/L_{ik}^e\) | \(D^e\) [`m^2 s^-1`]；\(A\) [`m^2`]；\(L\) [`m`] |
| \(I_{ik}^e=G_{ik}^e(\phi_i^e-\phi_k^e)\) | \(I^e\) [`A`]；\(G^e\) [`S`] |
| \(G_{ik}^e=\kappa_{ik}^{e,\mathrm{eff}}A_{ik}^e/L_{ik}^e\) | \(\kappa\) [`S m^-1`] |
| \(I_{ik}^{e,conc}=H_{ik}^e(\ln c_i^e-\ln c_k^e)\) | \(H^e\) [`A`] |
| \(M_{mn}^s=K_{mn}^s(c_m^s-c_n^s)\) | \(M^s\) [`mol s^-1`] |
| \(K_{mn}^s=D_{mn}^sA_{mn}^s/L_{mn}^s\) | \(D^s\) [`m^2 s^-1`] |
| \(I_{mn}^s=G_{mn}^s(\phi_m^s-\phi_n^s)\) | \(I^s\) [`A`] |
| \(G_{mn}^s=\sigma_{mn}^\mathrm{eff}A_{mn}^s/L_{mn}^s\) | \(\sigma\) [`S m^-1`] |
| \(Q_{ik}=g_{ik}^h(p_i-p_k)\), \(g_{ik}^h=\pi r_{ik}^4/(8\mu L_{ik})\) | \(Q\) [`m^3 s^-1`]；\(p\) [`Pa`]；\(\mu\) [`Pa s`] |

### 3.3 节点反应项

每个 electrolyte/NMC interface \(r\) 连接 electrolyte node \(i(r)\) 与 NMC node \(m(r)\)，面积为 \(A_r\)。overpotential 为 \(\eta_r=\phi_{m(r)}^s-\phi_{i(r)}^e-U(c_{m(r)}^s,c_{i(r)}^e)\)。反应电流 \(i_r\) 使用与 P2D 相同的 Butler-Volmer law。

单个 interface 对四个 residual 的贡献为 \(F_{c_e,i}+=-\beta_e(i_r/F)A_r\)、\(F_{\phi_e,i}+=+i_rA_r\)、\(F_{c_s,m}+=(i_r/F)A_r\)、\(F_{\phi_s,m}+=-i_rA_r\)。这保证质量源项符号与电荷源项符号一致，并保证 electrolyte 与 solid 之间的界面电荷守恒。

| 方程 | 变量定义与单位 |
|---|---|
| \(\eta_r=\phi_m^s-\phi_i^e-U(c_m^s,c_i^e)\) | \(\eta_r\) [`V`] |
| \(F_{c_e,i}+=-\beta_e(i_r/F)A_r\) | contribution [`mol s^-1`] |
| \(F_{\phi_e,i}+=+i_rA_r\) | contribution [`A`] |
| \(F_{c_s,m}+=(i_r/F)A_r\) | contribution [`mol s^-1`] |
| \(F_{\phi_s,m}+=-i_rA_r\) | contribution [`A`] |

### 3.4 与 P2D 的对应关系

P2D 是体积平均形式，PNM 是显微图结构有限体积形式。两者的主要对应如下：

| P2D 量 | PNM 估计或对应 | 单位 | 说明 |
|---|---|---:|---|
| \(\epsilon_e\) | \(\sum_iV_i^e/V_\mathrm{RVE}\) | `1` | electrolyte volume fraction |
| \(\epsilon_s\) | \(\sum_mV_m^\mathrm{NMC}/V_\mathrm{RVE}\) | `1` | active material volume fraction |
| \(a_s\) | \(\sum_rA_r/V_\mathrm{RVE}\) | `m^-1` | e/NMC interfacial area density |
| \(D_e^\mathrm{eff}\) | diffusion cell problem 上的 effective diffusivity | `m^2 s^-1` | 替代 \(D_e\epsilon^b\) |
| \(\kappa^\mathrm{eff}\) | ionic conduction cell problem 上的 effective conductivity | `S m^-1` | 替代 \(\kappa\epsilon^b\) |
| \(\sigma^\mathrm{eff}\) | solid graph electronic conduction cell problem | `S m^-1` | 应包含 NMC/CBD percolation |
| \(i_F(x)\) | \(\sum_{r\in \mathrm{RVE}}i_rA_r/\sum_rA_r\) | `A m^-2` | area-averaged Faradaic current |
| \(a_si_F\) | \(\sum_ri_rA_r/V_\mathrm{RVE}\) | `A m^-3` | volumetric reaction current |
| \(c_s^{surf}\) | PNM surface node concentration 或 embedded particle surface concentration | `mol m^-3` | 需要明确 node 是 bulk 还是 surface |

关键一致性要求：对任意 RVE，PNM 的体积积分守恒应还原 P2D 的源项，即 \(\int_\mathrm{RVE}\partial(\epsilon_ec_e)/\partial t\,dV=\sum_r\beta_e i_rA_r/F+\) boundary fluxes，\(\int_\mathrm{RVE}\partial(\epsilon_s\bar c_s)/\partial t\,dV=-\sum_ri_rA_r/F\)。

## 4. 耦合策略

### 4.1 方案 A：PNM 替代 P2D 的微观子模型

定义：P2D 的每个 cathode 宏观控制体 \(x_j\) 内嵌一个 PNM 子网络。P2D 负责沿厚度方向的宏观 electrolyte/solid 传输，PNM 负责该位置 RVE 内的 electrolyte/NMC/CBD 微观场、反应分布和局部有效响应。

数学形式：

- P2D 给第 \(j\) 个 PNM RVE 施加宏观状态 \(\bar c_{e,j}\)、\(\bar\phi_{e,j}\)、\(\bar\phi_{s,j}\)、\(\bar c_{s,j}\) 或宏观梯度 \(\nabla c_e\)、\(\nabla\phi_e\)、\(\nabla\phi_s\)。
- PNM 求解微观 residual \(F_\mathrm{PNM}(y_j;\bar u_j,\nabla\bar u_j)=0\)。
- PNM 返回 closure：\(\langle N_e\rangle_j\)、\(\langle i_e\rangle_j\)、\(\langle i_s\rangle_j\)、\(\langle a_si_F\rangle_j\)、局部 utilization。
- P2D 使用这些 closure 替代 \(D^\mathrm{eff}\nabla c\)、\(\kappa^\mathrm{eff}\nabla\phi\)、\(\sigma^\mathrm{eff}\nabla\phi\)、\(a_si_F\)。

优点：

- 微观结构最完整，能捕捉局部瓶颈、断连、反应热点。
- 可以研究显微结构随 \(x\) 变化的 graded electrode。
- 不依赖 Bruggeman 标量假设。

缺点：

- 计算成本最高，宏观每个 cell 都有一个 nonlinear PNM 子问题。
- Jacobian 复杂；若 monolithic Newton，需要 PNM closure 的敏感度。
- 初期验证困难，难以直接与标准 DFN 文献逐项对齐。

适用阶段：长期研究目标，不建议作为第一版实现。

### 4.2 方案 B：PNM 提供有效参数给 P2D

定义：先用 PNM 离线或在线 cell problem 计算 \(D_e^\mathrm{eff}\)、\(\kappa^\mathrm{eff}\)、\(\sigma^\mathrm{eff}\)、\(a_s\)、percolation factor、reaction accessibility 等参数，然后把这些参数注入 P2D，替代 Bruggeman 关系。

典型 cell problems：

- Diffusion：在 PNM electrolyte graph 两端施加 \(c_L-c_R\)，无反应，求总 molar flux \(Q_c\)，由 \(D_e^\mathrm{eff}=-Q_cL/(A\Delta c)\)。
- Ionic conduction：施加 \(\Delta\phi_e\)，无反应，求总 current \(I_e\)，由 \(\kappa^\mathrm{eff}=-I_eL/(A\Delta\phi)\)。
- Electronic conduction：在 solid graph 上施加 \(\Delta\phi_s\)，求 \(I_s\)，由 \(\sigma^\mathrm{eff}=-I_sL/(A\Delta\phi)\)。
- Interface area：直接计算 \(a_s=\sum_rA_r/V_\mathrm{RVE}\)。

优点：

- 与现有 P2D 方程兼容，最容易验证。
- 计算成本低，可以缓存 effective properties。
- 可以逐步替换 Bruggeman，明确量化 PNM 带来的差异。

缺点：

- 仍是均匀化模型，不能在 transient solve 中表达局部微观非均匀反应。
- 有效参数可能依赖 concentration、SOC、方向、局部 saturation，需要设计参数表或张量。

适用阶段：推荐第一版耦合方案。

### 4.3 方案 C：分层求解

定义：P2D 求解宏观场，PNM 在选定宏观位置或选定时间步求解微观场；两者通过界面条件迭代，但不一定 monolithic。可以理解为 operator splitting / heterogeneous multiscale method。

流程：

1. P2D 使用当前 closure 预测 \(\bar u^{m+1}\)。
2. 将 \(\bar u^{m+1}\) 或其局部梯度传递给 PNM RVE。
3. PNM 求解微观场，返回更新后的 closure \(\mathcal C^{m+1}\)。
4. 若 \(\|\mathcal C^{m+1}-\mathcal C^m\|\) 或 P2D residual 未收敛，继续外层 Picard/Newton 迭代。

优点：

- 比方案 A 更易实现，可按需只在少量 \(x\) 位置或时间点运行 PNM。
- 可以捕捉宏观状态对微观 closure 的非线性反馈。
- 适合 benchmark 方案 B 的误差。

缺点：

- 收敛性依赖 relaxation 与 closure 平滑性。
- 若只做 weak coupling，可能损失严格守恒或时间精度。

适用阶段：在方案 B 验证后实现，作为方案 A 的过渡。

### 4.4 推荐路线

推荐路线为：先实现完整独立 P2D，再采用方案 B 作为第一版 P2D+PNM 耦合，然后扩展到方案 C，最后保留方案 A 作为高保真研究路径。

理由：

- 数学验证顺序清晰：标准 DFN 可与 Newman/PyBaMM 对比，PNM effective properties 可与解析均匀介质和 Bruggeman 对比。
- 与现有代码风险最小：现有 PNM 已有网络、电导、BV、separator 与 transient solver，可复用为 effective property estimator。
- 先建立 P2D state/residual/solver 的可靠基线，避免把 DFN 错误和 PNM closure 错误混在一起。
- 方案 B 的产物（effective-property API、RVE metadata、方向性张量）也是方案 C/A 的基础。

## 5. 实现计划

### Phase 0：纯 P2D 模型独立实现（1D DFN）

目标：

- 实现不依赖 PNM 的 1D DFN/half-cell P2D solver。
- 支持 cathode + separator half-cell，预留 negative electrode full-cell 接口。
- 使用 finite volume in \(x\)、finite volume in particle \(r\)、Backward Euler 或 BDF1 time stepping、Newton nonlinear solve。
- 先支持 Bruggeman effective properties，作为基准。

新增模块建议：

- `src/pnmcathode/p2d/__init__.py`
- `src/pnmcathode/p2d/domain.py`
- `src/pnmcathode/p2d/state.py`
- `src/pnmcathode/p2d/materials.py`
- `src/pnmcathode/p2d/residual.py`
- `src/pnmcathode/p2d/solver.py`
- `src/pnmcathode/p2d/results.py`

关键类型与函数签名：

```python
from dataclasses import dataclass
from typing import Callable, Literal, Protocol
import numpy as np
from scipy import sparse

Array = np.ndarray

@dataclass(frozen=True)
class P2DRegion:
    name: Literal["negative", "separator", "positive"]
    x_left: float
    x_right: float
    n_cells: int
    epsilon_e: float
    epsilon_s: float = 0.0
    particle_radius: float | None = None
    specific_area: float | None = None
    bruggeman_e: float = 1.5
    bruggeman_s: float = 1.5

@dataclass(frozen=True)
class P2DMaterial:
    c_s_max: float
    diffusivity_s: Callable[[Array, float], Array]
    diffusivity_e: Callable[[Array, float], Array]
    conductivity_e: Callable[[Array, float], Array]
    conductivity_s: Callable[[Array, float], Array]
    ocv: Callable[[Array, float], Array]
    docv_dc_s: Callable[[Array, float], Array] | None = None

@dataclass
class P2DState:
    c_e: Array
    phi_e: Array
    phi_s: dict[str, Array]
    c_s: dict[str, Array]
    time: float = 0.0

@dataclass(frozen=True)
class P2DResidualContext:
    regions: tuple[P2DRegion, ...]
    material: dict[str, P2DMaterial]
    temperature: float
    current_density: float
    dt: float

def assemble_p2d_residual(
    state_new: P2DState,
    state_old: P2DState,
    context: P2DResidualContext,
) -> Array:
    ...

def assemble_p2d_jacobian(
    state_new: P2DState,
    state_old: P2DState,
    context: P2DResidualContext,
) -> sparse.csr_matrix:
    ...

class P2DSolver:
    def step(self, state: P2DState, dt: float, current_density: float) -> P2DState:
        ...

    def run(self, state0: P2DState, protocol: "P2DProtocol") -> "P2DResult":
        ...
```

数学实现重点：

- \(x\) 方向 finite volume：所有通量在 cell faces 上计算，区域界面使用 harmonic/series averaging。
- \(r\) 方向 finite volume：球壳体积 \(V_\ell=(4\pi/3)(r_{\ell+1}^3-r_\ell^3)\)，face area \(A_{\ell+1/2}=4\pi r_{\ell+1/2}^2\)。
- electrolyte concentration、solid particle concentration、\(\phi_e\)、\(\phi_s\) 全部进入同一个 nonlinear residual，避免 operator splitting 造成守恒偏差。
- potential gauge 必须固定，例如 half-cell 设 \(\phi_e\) at Li reference 为 0。
- 所有 residual 按物理量 scaling，避免 `mol/s` 和 `A` 混合导致 Newton 条件数过差。

验证方法：

- 固相球扩散：恒定表面通量下，与短时解析解或高精度径向解对比。
- 电解质扩散：无反应、固定边界浓度下，与 1D Fick 解析解对比。
- 电荷守恒：检查每个时间步 \(\max_x|\partial(i_e+i_s)/\partial x|\) 与总电流误差。
- 平衡态：\(I_\mathrm{app}=0\) 时 \(i_F=0\)，\(\phi_s-\phi_e=U\)，浓度不变。
- 小电流极限：电压接近 OCV，overpotential 与电流近似线性。

预计工作量：5-8 个工作日。

### Phase 1：P2D 验证（Newman/PyBaMM/文献基准）

目标：

- 将 Phase 0 的 P2D 与可信 DFN 实现或文献结果对比。
- 建立自动化 benchmark，锁定符号约定、单位、边界条件和电压定义。

新增模块建议：

- `tests/test_p2d_solid_particle.py`
- `tests/test_p2d_residual_conservation.py`
- `tests/test_p2d_half_cell.py`
- `examples/p2d/compare_pybamm.py`
- `examples/p2d/newman_baseline.py`

关键函数签名：

```python
@dataclass(frozen=True)
class P2DBenchmarkCase:
    name: str
    c_rate: float
    temperature: float
    cutoff_voltage: float
    expected_capacity_Ah_m2: float | None = None
    expected_voltage_path: str | None = None

def run_p2d_benchmark(case: P2DBenchmarkCase) -> "P2DResult":
    ...

def compare_voltage_curve(
    candidate_time: Array,
    candidate_voltage: Array,
    reference_time: Array,
    reference_voltage: Array,
) -> dict[str, float]:
    ...
```

验证方法：

- 与 PyBaMM DFN：同材料参数、同初始 SOC、同 C-rate，比较 \(V(t)\)、\(c_e(x,t)\)、\(c_s^{surf}(x,t)\)。
- 与 Newman/Doyle-Fuller 文献：低倍率和中倍率 discharge curves 的形状与极化分解。
- 与 Marquis et al. 参数集：使用公开 DFN benchmark 参数，避免只围绕 Khan PNM 参数调试。
- 数值收敛：将 \(N_x\)、\(N_r\)、\(\Delta t\) 加倍，确认电压和容量收敛阶。
- 守恒：全局 lithium inventory 误差应随 Newton tolerance 和 time discretization 收敛。

通过标准：

- OCV 和 \(I=0\) case 误差在 `1e-8` 到 `1e-6` 量级。
- 标准 DFN voltage curve 与参考曲线 RMS error 小于 `5-10 mV`，前提是参数完全一致。
- 总电流守恒相对误差小于 `1e-8` 或由线性求解器容差解释。

预计工作量：4-7 个工作日。

### Phase 2：耦合接口设计

目标：

- 定义 P2D 与 PNM 的最小稳定接口。
- 支持方案 B 的 effective-property injection，同时不阻碍方案 C/A。
- 统一 geometry normalization：RVE 体积、截面积、厚度方向、投影面积、interface area。

新增模块建议：

- `src/pnmcathode/coupling/__init__.py`
- `src/pnmcathode/coupling/effective.py`
- `src/pnmcathode/coupling/closures.py`
- `src/pnmcathode/coupling/rve.py`
- `tests/test_coupling_effective.py`

关键类型与函数签名：

```python
from typing import Protocol
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class RVEGeometry:
    volume: float
    area: float
    length: float
    direction: int

@dataclass(frozen=True)
class EffectiveTransport:
    epsilon_e: float
    epsilon_s: float
    specific_area: float
    diffusivity_e: float | np.ndarray
    conductivity_e: float | np.ndarray
    conductivity_s: float | np.ndarray
    percolation_e: bool
    percolation_s: bool

class EffectivePropertyProvider(Protocol):
    def evaluate(
        self,
        c_e: float,
        c_s: float,
        temperature: float,
        direction: int = 0,
    ) -> EffectiveTransport:
        ...

def estimate_effective_diffusivity(
    net: object,
    geometry: RVEGeometry,
    diffusivity_bulk: float,
    phase: str = "electrolyte",
) -> float:
    ...

def estimate_effective_conductivity(
    net: object,
    geometry: RVEGeometry,
    conductivity_bulk: float,
    phase: str,
) -> float:
    ...

def estimate_specific_area(net: object, volume: float) -> float:
    ...
```

设计要求：

- Effective properties 必须携带 direction，因为真实网络可能各向异性。
- 当 network 不 percolate 时，effective conductivity/diffusivity 应为 0，并返回 `percolation_* = False`，P2D solver 应给出清晰错误或退化处理。
- PNM 估计器不得将 pore-scale 几何 tortuosity 与 Bruggeman 再重复相乘。
- 所有 effective parameters 应可序列化，便于 benchmark 缓存。

验证方法：

- 均匀规则立方网络：effective property 应接近 bulk property 乘以几何孔隙率/连通修正。
- 人工断连网络：effective property 为 0 或接近 0。
- 方向性网络：沿不同 direction 的 effective property 可区分。
- \(a_s\) 直接求和与手工构造网络面积一致。

预计工作量：3-5 个工作日。

### Phase 3：耦合实现

目标：

- 第一子阶段实现方案 B：PNM-derived effective properties 注入 P2D。
- 第二子阶段实现方案 C：P2D-PNM 分层迭代 closure。
- 保持 P2D 独立 solver 可用，PNM coupling 是可选能力。

#### Phase 3B：PNM 提供有效参数给 P2D

新增或修改模块：

- 新增 `src/pnmcathode/coupling/p2d_adapter.py`
- 修改 `src/pnmcathode/p2d/domain.py`，允许 `EffectivePropertyProvider`
- 修改 `src/pnmcathode/p2d/residual.py`，effective coefficients 通过 provider 获取

关键函数签名：

```python
@dataclass(frozen=True)
class PNMEffectiveProvider:
    net: object
    geometry: RVEGeometry
    cache: dict[tuple[float, float, float, int], EffectiveTransport] | None = None

    def evaluate(
        self,
        c_e: float,
        c_s: float,
        temperature: float,
        direction: int = 0,
    ) -> EffectiveTransport:
        ...

def build_p2d_region_from_pnm(
    name: str,
    x_left: float,
    x_right: float,
    n_cells: int,
    provider: EffectivePropertyProvider,
    reference_c_e: float,
    reference_c_s: float,
    temperature: float,
) -> P2DRegion:
    ...
```

验证方法：

- 当 provider 返回 Bruggeman 等效参数时，coupled P2D 与 pure P2D 结果逐步一致。
- 当 provider 返回 PNM effective parameters 时，电压差、容量差和极化项方向符合物理预期：较低 \(\kappa^\mathrm{eff}\) 或 \(D_e^\mathrm{eff}\) 导致更大浓差/欧姆极化。
- 与现有 PNM 单网络 discharge 在相同体积分数和材料下做趋势对比，而不是要求逐点相同。

预计工作量：4-6 个工作日。

#### Phase 3C：分层 P2D-PNM 求解

新增模块：

- `src/pnmcathode/coupling/hierarchical.py`
- `src/pnmcathode/coupling/micro_problem.py`
- `tests/test_coupling_hierarchical.py`

关键函数签名：

```python
@dataclass(frozen=True)
class MacroCellState:
    c_e: float
    c_s_avg: float
    phi_e: float
    phi_s: float
    grad_c_e: float
    grad_phi_e: float
    grad_phi_s: float
    temperature: float

@dataclass(frozen=True)
class MicroClosure:
    electrolyte_flux: float
    ionic_current: float
    electronic_current: float
    volumetric_reaction_current: float
    average_faraday_current: float
    utilization: float

class PNMMicroProblem(Protocol):
    def solve(self, macro: MacroCellState, previous: object | None = None) -> MicroClosure:
        ...

def hierarchical_p2d_pnm_step(
    state: P2DState,
    dt: float,
    current_density: float,
    micro_problems: list[PNMMicroProblem],
    max_outer_iter: int = 20,
    relaxation: float = 0.5,
) -> P2DState:
    ...
```

验证方法：

- 如果 micro closure 固定为 linear effective transport，方案 C 应退化为方案 B。
- 外层迭代 residual 单调下降；不下降时启用 relaxation。
- 对单个 macro cell + 单个 PNM RVE，与直接 PNM transient 的体积平均结果对比。
- 检查每个 macro cell 的 micro closure 满足 \(\langle i_e\rangle+\langle i_s\rangle=I_\mathrm{total}\)。

预计工作量：8-12 个工作日。

### Phase 4：验证与 benchmark

目标：

- 建立从 pure DFN 到 PNM-corrected DFN 到 hierarchical P2D+PNM 的完整 benchmark。
- 量化 PNM coupling 相对 Bruggeman P2D 的收益、误差与成本。

新增内容：

- `examples/coupling/compare_bruggeman_vs_pnm_effective.py`
- `examples/coupling/hierarchical_single_rve.py`
- `docs/P2D_PNM_VALIDATION.md`
- benchmark 数据缓存目录，例如 `data/benchmarks/`

benchmark 维度：

- C-rate：`0.2C`, `0.5C`, `1C`, `3C`, `5C`。
- 电极厚度：Khan 1CAL/3CAL 参数、合成厚/薄电极。
- 孔隙率：低、中、高三档。
- CBD conductivity：断连、临界、充分连通。
- transport closure：Bruggeman、PNM scalar、PNM tensor、hierarchical。

核心指标：

- \(V(t)\)、\(V(Q)\)、容量到 cutoff。
- electrolyte concentration minimum 和 gradient。
- solid surface/bulk concentration gradient。
- reaction current distribution：\(\mathrm{std}(i_F)\)、Gini coefficient、active area utilization。
- 欧姆极化、浓差极化、charge-transfer overpotential 分解。
- 运行时间、Newton iterations、outer coupling iterations。

验证方法：

- Pure P2D 与 PyBaMM/Newman benchmark 对比。
- PNM effective transport 与解析网络/有限体积均匀介质对比。
- P2D+PNM 低异质性极限应收敛到 Bruggeman/均匀 P2D。
- 高异质性网络应表现出更强局部反应不均匀和更早局部 electrolyte depletion。
- 与 Khan 2021 的趋势对比：1CAL/3CAL 厚度差异、高倍率容量衰减趋势、NMC532 参数下的电压范围。

预计工作量：6-10 个工作日。

## 6. 数值与工程风险

### 6.1 符号约定风险

P2D 文献常把 discharge current 定义为正，而当前项目中 cathode discharge 对应 \(I_\mathrm{app}<0\)。实现时必须在 API 层明确：

- 用户输入 `c_rate > 0` 表示 discharge magnitude。
- solver 内部转换为 \(I_\mathrm{app}<0\)。
- \(i_F<0\) 表示 cathode lithiation。
- voltage 在 discharge 中应随时间下降。

建议测试：`test_discharge_sign_convention()` 检查 1C discharge 第一步电压低于 OCV 且 \(\bar c_s\) 增加。

### 6.2 OCV 与 electrolyte reference 风险

NMC OCV curve 通常相对 Li/Li+ reference。若再加入 \(c_e\) Nernst correction，会重复计入 electrolyte chemical potential。实现中应让 `ocv` 函数声明 reference：

```python
@dataclass(frozen=True)
class OCVModel:
    reference: Literal["li_metal_empirical", "thermodynamic_with_electrolyte"]
    value: Callable[[Array, Array, float], Array]
```

### 6.3 DAE 刚性与 Jacobian

DFN 是 stiff nonlinear DAE-like system。建议 Phase 0 先实现 analytic sparse Jacobian；若初期使用 finite-difference Jacobian，只能用于小规模验证，不适合作为目标实现。

### 6.4 Effective property 双重计数

如果 PNM 已通过真实 throat geometry 计算 tortuosity，就不能再乘 \(\epsilon^b\)。P2D region 应区分：

- `closure="bruggeman"`：使用 \(D\epsilon^b\)。
- `closure="pnm_effective"`：直接使用 \(D^\mathrm{eff}_\mathrm{PNM}\)。

### 6.5 PNM node \(c_s\) 与 P2D particle \(c_s(r)\) 的含义差异

现有 PNM solid node 的 \(c_s\) 更接近 active-material control volume 平均浓度，而 DFN Butler-Volmer 需要 surface concentration \(c_s^{surf}\)。耦合时有三种处理：

- 方案 B：P2D 仍保留 spherical particle diffusion，PNM 只给 \(a_s\)、transport coefficients。
- 方案 C：PNM microproblem 内每个 NMC node 可挂一个 spherical particle submodel。
- 方案 A：PNM 替代宏观子模型时，必须明确 node concentration 是 surface、bulk average 还是 radial subgrid state。

推荐 Phase 3B 采用第一种，避免物理含义混淆。

## 7. 建议交付物清单

| Phase | 交付物 | 完成判据 |
|---|---|---|
| Phase 0 | Pure P2D solver | 通过解析极限、守恒、平衡态测试 |
| Phase 1 | P2D benchmark suite | 与 PyBaMM/Newman voltage curve 误差达标 |
| Phase 2 | Coupling API + effective estimators | PNM effective properties 可重复、可缓存、单位正确 |
| Phase 3B | PNM-effective P2D | Bruggeman provider 退化一致，PNM provider 产生合理差异 |
| Phase 3C | Hierarchical solver | 固定 closure 退化到 3B，动态 closure 外层收敛 |
| Phase 4 | Validation report | 完整倍率、厚度、孔隙率 benchmark 与文档 |

## 8. 参考文献

1. M. Doyle, T. F. Fuller, and J. Newman, “Modeling of Galvanostatic Charge and Discharge of the Lithium/Polymer/Insertion Cell,” *Journal of The Electrochemical Society*, 140, 1526-1533, 1993.
2. T. F. Fuller, M. Doyle, and J. Newman, “Simulation and Optimization of the Dual Lithium Ion Insertion Cell,” *Journal of The Electrochemical Society*, 141, 1-10, 1994.
3. J. Newman and K. E. Thomas-Alyea, *Electrochemical Systems*, 3rd ed., Wiley-Interscience, 2004.
4. M. Doyle, J. Newman, A. S. Gozdz, C. N. Schmutz, and J. M. Tarascon, “Comparison of Modeling Predictions with Experimental Data from Plastic Lithium Ion Cells,” *Journal of The Electrochemical Society*, 143, 1890-1903, 1996.
5. S. G. Marquis, V. Sulzer, R. Timms, C. P. Please, and S. J. Chapman, “An asymptotic derivation of a single particle model with electrolyte,” *Journal of The Electrochemical Society*, 166, A3693-A3706, 2019.
6. V. Sulzer, S. G. Marquis, R. Timms, M. Robinson, and S. J. Chapman, “Python Battery Mathematical Modelling (PyBaMM),” *Journal of Open Research Software*, 9, 14, 2021.
7. Z. A. Khan, A. Elkamel, and J. T. Gostick, “Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes,” *Journal of The Electrochemical Society*, 168, 070534, 2021.
8. M. Ebner, D.-W. Chung, R. E. García, and V. Wood, “Tortuosity Anisotropy in Lithium-Ion Battery Electrodes,” *Advanced Energy Materials*, 4, 1301278, 2014.
9. D. W. Chung, M. Ebner, D. R. Ely, V. Wood, and R. E. García, “Validity of the Bruggeman relation for porous electrodes,” *Modelling and Simulation in Materials Science and Engineering*, 21, 074009, 2013.
10. J. Gostick et al., “OpenPNM: A Pore Network Modeling Package,” *Computing in Science & Engineering*, 18, 60-74, 2016.

