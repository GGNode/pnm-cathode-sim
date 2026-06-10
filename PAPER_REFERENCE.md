# Khan et al. 2021 — 完整论文数据参考

**来源**: Khan et al., "Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes", J. Electrochem. Soc. 168 (2021) 070534
**PDF**: `data/khan2021_pnm_lib_cathode.pdf`
**提取方式**: 视觉读取 (mimo-v2.5) + pymupdf 文本提取，交叉验证

---

## Table I: XCT 图像属性

| Property | NMC532 (1CAL) | NMC532 (3CAL) |
|----------|--------------|--------------|
| Coating Thickness | 129 μm | 88 μm |
| Electrolyte phase fraction | 36.8% | 36.6% |
| Active Material (NMC532) fraction | 49.28% | 49.74% |
| Carbon and Binder Domain (CBD) | 13.92% | 13.66% |
| Density | 29.78 mg/cm² | 20.40 mg/cm² |
| Experimental capacity | 178 mAh/g_NMC | 182 mAh/g_NMC |
| Image size (voxels) | 320×250×250 | 221×250×250 |
| Voxel Size | 398 nm | 398 nm |
| Extracted Network (Nodes, Bonds) | 4637, 31427 | 3510, 23126 |
| Electrolyte nodes/bonds | 1238, 1943 | 982, 1459 |
| AM nodes/bonds | 1722, 3736 | 1295, 2760 |
| CBD nodes/bonds | 1680, 3375 | 1233, 2310 |
| Electrolyte-NMC532 bonds | 8151 | 6263 |
| NMC532-CBD bonds | 7620 | 5355 |
| Electrolyte-CBD bonds | 6602 | 4979 |

## Table II: 模型参数

| Parameter | Symbol | Units | Value |
|-----------|--------|-------|-------|
| Separator Thickness | L_s | μm | 25 |
| Separator Porosity | ε_sep | % | 39 |
| Electronic conductivity of NMC532 | σ_AM | S/m | 0.01 |
| Electronic conductivity of CBD | σ_CBD | S/m | 760 |
| Ionic conductivity of electrolyte in separator | K_s | S/m | Eq.(1) |
| Initial Li+ concentration | c_2,0 | mol/m³ | 1200 |
| Max Li concentration in NMC532 | c_1,max | mol/m³ | 48900 |
| Diffusion coefficient of Li+ in electrolyte | D_Li+ | m²/s | Eq.(1) |
| Li diffusion coefficient in AM | D_Li | m²/s | Eq.(2) |
| Faraday constant | F | C/mol | 96485 |
| Cathode rate constant | k | mol/m²·s | **1e-10** |
| Bruggeman exponent (separator) | p_sep | — | 1.5 |
| Ideal gas constant | R | J/mol·K | 8.314 |
| Reference temperature | T | K | **303** |
| Li transference number | t_+ | — | 0.363 |
| Cut off Voltage | V | V | 3.0 |
| Charge transfer coefficient (NMC) | α_c | — | 0.5 |
| Charge transfer coefficient (Li foil) | α_a | — | 0.5 |
| Exchange current density of Li foil | i_foil^0 | A/m² | 19 |

### Footnote (1): 电解质离子电导率 K_s

```
K_s = 10^(-4.43 - 54.0/(T - 229 - 5*c_2) - 0.22*c_2)
```

### Footnote (1): 电解质 Li+ 扩散系数 D_Li+

```
D_Li+ = 10^(-4.43 - 54.0/(T - 229 - 5*c_2) - 0.22*c_2)
```

### Footnote (2): 固相 Li 扩散系数 D_Li

```
D_Li = 10^(-2319*soc^10 + 6642*soc^9 - 5269*soc^8 - 3319*soc^7
         + 10038*soc^6 - 9806*soc^5 + 5817*soc^4 - 2286*soc^3
         + 575.3*soc^2 - 83.16*soc - 9.292)
```

### Footnote (2): 电解质电导率 (完整公式)

```
κ_eff = c_2 * (-10.5 + 0.0740*T - 6.96e-5*T^2 + 0.668*c_2
         - 0.0178*c_2*T + 2.80e-5*c_2*T^2 + 0.494*c_2^2
         - 8.86e-4*c_2^2*T)^2
```

---

## 核心方程

### Eq. 2.2: 电解质浓度 (Nernst-Planck 离散形式)

```
v_i * (c_i^(t+Δt) - c_i^t) / Δt
  = Σ_{j=1}^{k_elec} [G_ij^d + max(-G_ij^m(φ_j^t - φ_i^t), 0)] * c_i^t
  - Σ_{j=1}^{k_elec} [G_ij^d + max(G_ij^m(φ_i^t - φ_j^t), 0)] * c_j^t
  + a * j_{n,rp}^t
```

### Eq. 2.3-2.4: 扩散/迁移电导

```
G_ij^d = (1/g_i^d + 1/g_ij^d + 1/g_j^d)^(-1)
G_ij^m = (1/g_i^m + 1/g_ij^m + 1/g_j^m)^(-1)
```

### Eq. 2.5-2.6: 单位电导

```
g_i^d = A_i * D_Li+ / l_i
g_i^m = (zF/RT) * g_i^d
```

### Eq. 2.7: 电解质电位

```
i_{i,elec} = -zF * (Σ_{j=1}^{k_elec} G_ij^d * (c_i^t - c_j^t))
             - Σ_{j=1}^{k_elec} K_ij^elec * (φ_i^t - φ_j^t)
```

### Eq. 2.8-2.9: 离子电导

```
K_ij^elec = (1/k_i^elec + 1/k_ij^elec + 1/k_j^elec)^(-1)
k_i^elec = zF * g_i^m * c_i^t
```

### Eq. 2.10-2.12: 固相扩散

```
v_{i,AM} * (c_{i,Li}^(t+Δt) - c_{i,Li}^t) / Δt
  = Σ_{j=1}^{k_AM} G_{ij}^{sd} * (c_{i,Li}^t - c_{j,Li}^t)

G_ij^sd = (1/g_i^sd + 1/g_ij^sd + 1/g_j^sd)^(-1)
g_i^sd = A_i * D_Li / l_i
```

### Eq. 2.13-2.15: 固相电位

```
i_{i,AM} = -Σ_{j=1}^{k_AM} G_{ij}^AM * (φ_{i,AM}^t - φ_{j,AM}^t)
0 = -Σ_{j=1}^{k_CBD} G_{ij}^CBD * (φ_{i,CBD}^t - φ_{j,CBD}^t)
g_ij^e = A_i * σ_{AM/CBD} / l_i
```

### Eq. 2.16: Butler-Volmer

```
a*j_{n,i}^t = i_i^{0,t} * (exp(((1-α_c)*F)/(RT) * η_{i,c}^t)
                           - exp((-α_c*F)/(RT) * η_{i,c}^t))
```

### Eq. 2.17: 交换电流密度

```
i_i^0 = Σ_{j=1}^{k_elec} (a_ij * F * k * (c_{j,elec}^t)^(1-α_c)
         * (c_{i,AM,max}^t - c_{i,AM}^t)^(1-α_c)
         * (c_{i,AM}^t)^(α_c))
```

### Eq. 2.18: 过电位

```
η_{i,c}^t = Σ_{j=1}^{k_elec} φ_{i,AM}^t - φ_{j,elec}^t - U_{i,cathode}^t
```

### Eq. 2.19: OCV (Khan 多项式)

```
U_{i,cathode} = 5744.862289*soc^9 - 5520.41099*soc^8
  + 95714.29862*soc^7 - 147364.5514*soc^6
  + 142718.3782*soc^5 - 90095.81521*soc^4
  + 37061.41195*soc^3 - 9578.599274*soc^2
  + 1409.309503*soc - 85.31153081
  - 0.0003 * exp(7.657 * (soc^115))
```

### Eq. 2.20: 电流守恒

```
∇i_eleC + ∇i_AM = 0
```

### Eq. 2.21-2.23: 分离器方程

```
∂c_{2,s}/∂t = -D_{Li,s}^eff * ∇²c_{2,s}
∂φ_{2,s}/∂x = -i_2/κ_s^eff + (RT/F)*(1-t_+) * ∂ln(c_{2,s})/∂x
∂c_{2,s}/∂x = I*(1-t_+) / (F*D_{Li,s})
```

### Eq. 2.24: Li 箔阳极 BV

```
∇i_eleC = i_foil^0 * (exp(((1-α_a)*F)/(RT) * (φ_Li - φ_{2,s}))
                      - exp((-α_a*F)/(RT) * (φ_Li - φ_{2,s})))
```

---

## Figure 4: V-Q 放电曲线 (视觉读取)

### 1CAL

| C-rate | 容量 (mAh/g) |
|--------|-------------|
| 0.2C | ~195 |
| 0.5C | ~185 |
| 1C | ~175 |
| 3C | ~150 |

### 3CAL

| C-rate | 容量 (mAh/g) |
|--------|-------------|
| 0.2C | ~190 |
| 0.5C | ~185 |
| 1C | ~175 |
| 3C | ~170 |

- X 轴: Capacity (mAh/g_NMC), 0-200
- Y 轴: Cell Voltage (V), 3.0-4.4
- 实线 = PNM 模型, 虚线 = 实验数据
- 初始电压: ~4.0V (Khan OCV at x≈0.35)
- 3CAL 在高倍率下容量衰减更慢 (更薄, 传输路径更短)

## Figure 6: 3D 空间分布 (3C, 75% SoL)

| 子图 | 电极 | 物理量 | 范围 |
|------|------|--------|------|
| (a) | 1CAL | Li⁺ 浓度 | 0-1550 mol/m³ |
| (b) | 1CAL | 电解质电位 | -0.41 ~ -0.09 V |
| (c) | 1CAL | SoL | 0.41-0.99 (41%-99%) |
| (d) | 1CAL | 固相电位 | 0.3446-0.3467 V |
| (e) | 3CAL | Li⁺ 浓度 | 0-1550 mol/m³ |
| (f) | 3CAL | 电解质电位 | -0.1 ~ -0.07 V |
| (g) | 3CAL | SoL | 0.31-0.94 (31%-94%) |
| (h) | 3CAL | 固相电位 | 0.3420-0.3427 V |

## Figure 7: 逐孔 Profile (视觉读取)

### 1CAL 厚度: ~120-125 μm
### 3CAL 厚度: ~80 μm

| 条件 | 物理量 | 隔膜端 | 集流体端 | 梯度 |
|------|--------|--------|---------|------|
| 1CAL, 1C | c_e | ~1400 | ~700 mol/m³ | 中等 |
| 1CAL, 3C | c_e | ~1400 | ~100 mol/m³ | **强** |
| 1CAL, 1C | φ_e | -0.052 | -0.068 V | 小 |
| 1CAL, 3C | φ_e | -0.010 | -0.045 V | 大 |
| 3CAL, 1C | c_e | ~1300 | ~1000 mol/m³ | 小 |
| 3CAL, 3C | c_e | ~1400 | ~400 mol/m³ | 中等 |
| 3CAL, 1C | φ_e | -0.032 | -0.038 V | 极小 |
| 3CAL, 3C | φ_e | -0.070 | -0.100 V | 中等 |

### 孔径分类
- 1CAL: 0.9-4.5 μm (小), 4.5-8.0 μm (中), 8.0-12 μm (大)
- 3CAL: 0.1-3.35 μm (小), 3.35-6.65 μm (中), 6.65-10 μm (大)

### SoL 分布
- 1CAL, 1C: SoL ~0.8-0.9, 分布较均匀
- 1CAL, 3C: SoL 从隔膜端 ~0.9 降至集流体端 ~0.45, **强梯度**
- 3CAL: SoL 分布更均匀 (更短传输路径)

---

## 求解工作流 (Page 7)

1. 指定参数、边界条件、初始值
2. 1D 分离器扩散模型 + PNM 阴极模型耦合
3. 迭代求解: 猜测集流体电位 → 解电位分布 → 计算电流 → 检查收敛 (tolerance 1e-3)
4. 更新场变量, 进入下一时间步
5. 重复直到电池电压 < 3.0V

### 关键设计选择
- 分离器用 1D 模型 (非 PNM), 因为 c_e 在分离器/阴极界面随时间变化
- CBD 相人工添加到 XCT 图像中
- 阳极 Li foil 用 Eq. 2.24 BV 表达式
- 无热效应
- 无体积膨胀

---

## 计算性能 (Page 11)

- 硬件: Lenovo ThinkStation P720, 8-core Xeon Silver 4110 @ 2.1GHz, 256GB RAM
- 0.2C: ~8 小时
- 3C: ~5 小时
- 1C: ~6 小时 (约 6x 慢于实时)
- 对比 DNS 模拟: 数天到数周, 本方法快 ~115x

---

## 与我们代码的关键差异

| 项目 | 论文 | 我们 | 状态 |
|------|------|------|------|
| k (反应速率) | **1e-10** | 5e-10 | ❌ 差 5x |
| OCV | Khan 多项式 | 文献查找表 | ❌ 不同曲线 |
| 电解质孔隙率 | **36.8%** (1CAL) | 50% | ❌ 偏高 |
| D_e | 浓度依赖 Eq.(1) | 常数 2e-10 | ❌ 简化 |
| κ | 浓度依赖 Eq.(2) | 常数 | ❌ 缺失 |
| D_s | SoC 依赖 Eq.(2) | 不同公式 | ⚠️ |
| σ_AM / σ_CBD | 0.01 / 760 S/m | 未区分 | ❌ 缺失 |
| T | 303 K | 298.15 K | ❌ 差 5K |
| 网络 | XCT 4637 节点 | 合成 1000 节点 | ⚠️ 简化 |
| 分离器 | 1D 模型 | 固定 c_e=1200 | ⚠️ 简化 |
| CBD | 人工添加 | 有 (10%) | ✅ |
