# PNM-LIB-Cathode: 复现论文

## 目标

严格复现 Khan, Elkamel, Gostick (2021) 的孔网络模型：
> "Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes"
> J. Electrochem. Soc. 168(7), 070534. DOI: 10.1149/1945-7111/ac120c

**核心要求：** 严格推导 → 设计审核 → 分阶段实现 → 验证

---

## 一、论文核心贡献（严格推导）

### 1.1 模型架构

这是一个**瞬态孔网络模型**，模拟锂离子电池半电池（Li 负极 | 隔膜 | NMC532 正极）在恒流放电下的行为。

**关键创新：** 首次将 PNM 用于 LIB 放电的瞬态模拟。此前 PNM 仅用于稳态或半稳态。

### 1.2 网络提取（3-phase XCT → PNM）

从 X 射线断层扫描图像提取三相孔网络：

1. **相分割：** 体素标签 = {0: 电解液+CBD, 1: NMC532}
2. **CBD 合成：** 在颗粒接触点通过形态学操作添加碳粘结剂域
3. **分水岭提取：** 用分水岭算法将孔空间分割为孔体(pore bodies)和孔喉(pore throats)
4. **输出：** 节点(node) = 孔体，键(bond) = 孔喉，包含几何属性（体积、面积、形状因子）

**网络规模（论文数据）：**
- 1CAL: 4637 nodes, 31427 bonds
- 3CAL: 3510 nodes, 23126 bonds

### 1.3 控制方程（严格推导）

模型在孔网络上求解以下耦合瞬态方程组：

#### A. 电解液相 — Li⁺ 浓度（Nernst-Planck 方程）

在每个孔节点上，Li⁺ 浓度 cₑ 的瞬态守恒：

```
V_pore · ∂cₑ/∂t = Σ_bonds [D_eff · A_bond / L_bond · (cₑ_neighbor - cₑ_current)] + Σ_interfaces [r_rxn · A_intf]
```

其中：
- V_pore: 孔体积 [m³]
- D_eff: 有效扩散系数 [m²/s]，考虑迂曲度
- A_bond: 孔喉截面积 [m²]
- L_bond: 孔喉长度 [m]
- r_rxn: 界面反应速率 [mol/(m²·s)]，由 Butler-Volmer 给出
- A_intf: 电解液/活性材料界面面积 [m²]

**有效扩散系数：**
```
D_eff = D_electrolyte · ε^1.5   （Bruggeman 关系）
```

#### B. 电解液相 — 电位（Ohm 定律）

在每个孔节点上，电解液电位 φₑ 的守恒：

```
Σ_bonds [κ_eff · A_bond / L_bond · (φₑ_neighbor - φₑ_current)] + Σ_interfaces [I_rxn · A_intf] = 0
```

（准稳态假设：电位弛豫远快于浓度变化）

其中：
- κ_eff: 电解液有效电导率 [S/m]
- I_rxn: 界面反应电流密度 [A/m²]

**电解液电导率（依赖于浓度）：**
```
κ(cₑ) = cₑ · F² · D_electrolyte / (R · T)   （简化 Nernst-Einstein 关系）
```

#### C. 固相（NMC532）— Li 浓度（Fick 定律）

在每个 NMC 节点上，固相 Li 浓度 cₛ 的瞬态扩散：

```
∂cₛ/∂t = D_s · ∇²cₛ
```

在网络上的离散形式：
```
V_NMC · ∂cₛ/∂t = Σ_bonds_solid [D_s · A_bond / L_bond · (cₛ_neighbor - cₛ_current)] - r_rxn · A_intf
```

其中：
- D_s: NMC532 中 Li 的固相扩散系数 [m²/s]
- D_s = f(cₛ, T) — 浓度和温度的非线性函数（见下）

**NMC532 的 Li 扩散系数（论文给出的经验公式）：**
```
D_s = cₛ² · (-10.5 + 0.0740·T - 6.96e-5·T² + 0.668·cₛ² - 0.0178·cₛ²·T + 2.80e-5·cₛ²·T² + 0.494 - 8.86e-4·T)²
```

注意：这个公式看起来很复杂，需要仔细核实。论文中 cₛ 和 T 有特定单位。

#### D. 固相（NMC532 + CBD）— 电子电位（Ohm 定律）

```
Σ_bonds_solid [σ_eff · A_bond / L_bond · (φₛ_neighbor - φₛ_current)] - I_rxn · A_intf = 0
```

（准稳态假设）

其中：
- σ_eff: 有效电子电导率 [S/m]
- NMC532: σ = 0.01 S/m
- CBD: σ = 760 S/m

#### E. 界面反应 — Butler-Volmer 动力学

在电解液/NMC 界面上：

```
I_rxn = i₀ · [exp(αₐ·F·η/(R·T)) - exp(-αc·F·η/(R·T))]
```

其中：
- i₀: 交换电流密度 [A/m²]
- αₐ = αc = 0.5: 传递系数（对称 Butler-Volmer）
- η: 过电位 [V]
- F: 法拉第常数 [C/mol]
- R: 气体常数 [J/(mol·K)]
- T: 温度 [K]

**过电位：**
```
η = φₛ - φₑ - U_eq(cₛ_surface, cₑ)
```

**平衡电位（Nernst 方程）：**
```
U_eq = U₀ + (R·T)/(F) · ln((cₛ_max - cₛ_surface)/cₛ_surface · cₑ/cₑ_ref)
```

或使用 NMC532 的经验开路电位曲线 U_eq(x)，其中 x = cₛ/cₛ_max 是锂化程度。

**交换电流密度：**
```
i₀ = k₀ · F · (cₑ)^αₐ · (cₛ_max - cₛ_surface)^αₐ · (cₛ_surface)^αc
```

### 1.4 边界条件

#### 隔膜（1D 模型）
- 厚度: 25 µm
- 孔隙率: 39%
- 边界: 电解液浓度 cₑ = 1200 mol/m³（初始值，与储液器相连）
- 电位: φₑ = 0 V（参考电位）

#### Li 负极（简化为恒定边界）
- 电位: φₛ = 0 V（短路到负极）
- Li 浓度: 恒定（无限储锂源）

#### 集流体侧
- 电子电流: I_applied = C-rate × 容量标称值 / 面积
- 无离子通量

### 1.5 恒流放电控制

```
总施加电流 I_total = C_rate × Q_nominal / t_nominal
```

在每个时间步：
1. 假设过电位分布
2. 计算 Butler-Volmer 反应电流
3. 求解电解液浓度、电位、固相浓度和电位
4. 迭代直到收敛
5. 更新端电压 V_cell = φₛ(collector) - φₑ(separator)
6. 如果 V_cell < 3.0 V（截止电压），停止

### 1.6 关键模型参数

| 参数 | 值 | 单位 | 来源 |
|---|---|---|---|
| 隔膜厚度 | 25e-6 | m | 论文 |
| 隔膜孔隙率 | 0.39 | - | 论文 |
| 初始电解液 Li⁺ 浓度 | 1200 | mol/m³ | 论文 |
| NMC532 最大 Li 浓度 | 48900 | mol/m³ | 论文 |
| NMC532 电子电导率 | 0.01 | S/m | 论文 |
| CBD 电子电导率 | 760 | S/m | 论文 |
| 电解液 Li⁺ 扩散系数 | 2.0e-10 (估) | m²/s | 文献 |
| 温度 | 298.15 | K | 标准 |
| 截止电压 | 3.0 | V | 论文 |
| αₐ, αc | 0.5, 0.5 | - | 论文 |

---

## 二、复现设计方案

### 2.1 技术栈

**核心依赖：**
- Python 3.10+
- OpenPNM (Open Pore Network Modeling) — 孔网络建模框架
- PoreSpy — 图像分析和网络提取
- scipy — 稀疏矩阵求解器
- numpy, matplotlib — 数值计算和可视化

**备选：** 如果 OpenPNM 不满足需求（如瞬态求解器），需自定义求解器。

### 2.2 模块设计

```
pnm-lib-cathode/
├── README.md                    # 项目说明
├── requirements.txt             # 依赖
├── PROMPT.md                    # 本设计文档
├── src/
│   ├── __init__.py
│   ├── network/                 # 孔网络相关
│   │   ├── __init__.py
│   │   ├── generator.py         # 生成/加载孔网络
│   │   ├── properties.py        # 计算网络几何属性
│   │   └── synthetic.py         # 合成网络（无 XCT 数据时）
│   ├── physics/                 # 物理模型
│   │   ├── __init__.py
│   │   ├── electrolyte.py       # 电解液相（Nernst-Planck + Ohm）
│   │   ├── solid.py             # 固相（Fick + Ohm）
│   │   ├── reaction.py          # Butler-Volmer 动力学
│   │   ├── boundary.py          # 边界条件（隔膜、负极、集流体）
│   │   └── ocv.py               # NMC532 开路电位曲线
│   ├── solver/                  # 求解器
│   │   ├── __init__.py
│   │   ├── transient.py         # 瞬态求解器（时间步进）
│   │   ├── steady.py            # 稳态求解器（电位分布）
│   │   └── newton.py            # Newton-Raphson 非线性迭代
│   └── post/                    # 后处理
│       ├── __init__.py
│       ├── visualization.py     # 3D 网络可视化
│       └── analysis.py          # 性能分析（V-Q 曲线等）
├── configs/
│   ├── nmc532_1cal.yaml         # 1CAL 电极参数
│   └── nmc532_3cal.yaml         # 3CAL 电极参数
├── tests/
│   ├── test_butler_volmer.py    # BV 动力学单元测试
│   ├── test_diffusion.py       # 扩散方程单元测试
│   ├── test_single_pore.py     # 单孔验证
│   └── test_voltage_curve.py   # V-Q 曲线集成测试
├── scripts/
│   ├── run_discharge.py         # 运行放电模拟
│   ├── generate_network.py      # 生成合成孔网络
│   └── plot_results.py          # 绘制结果
└── data/
    └── README.md                # 数据说明（NREL XCT 数据链接）
```

### 2.3 实现阶段（分阶段，每阶段独立验证）

#### Phase 0: 项目骨架 + 依赖
- 初始化项目结构
- 安装 OpenPNM, PoreSpy, scipy, numpy, matplotlib
- 验证 OpenPNM 基本功能

#### Phase 1: 单孔验证（核心物理验证）
- 实现 Butler-Volmer 动力学模块
- 实现 NMC532 OCV 曲线
- 实现固相扩散（单球形颗粒）
- **验证：** 单颗粒放电曲线与解析解对比（Cottrell 方程极限）

#### Phase 2: 网络生成
- 使用 PoreSpy 生成随机球堆积网络
- 或使用 OpenPNM 内置的随机网络生成器
- 分配三相属性（电解液、NMC、CBD）
- 计算几何属性（体积、面积、配位数）

#### Phase 3: 稳态求解器
- 实现电解液电位分布（Ohm 定律）
- 实现固相电位分布
- 实现 Butler-Volmer 耦合
- **验证：** 0.1C 低倍率下 V-Q 曲线（接近平衡）

#### Phase 4: 瞬态求解器
- 添加电解液浓度时间导数项
- 添加固相 Li 浓度时间导数项
- 实现时间步进（隐式或半隐式格式）
- **验证：** 0.2C, 0.5C, 1C 放电曲线

#### Phase 5: 后处理 + 分析
- 实现空间分布可视化（浓度、电位、SoL）
- 实现结构分析（孔径分布、配位数、迂曲度）
- 生成论文中 Figure 2-8 的等价图

### 2.4 求解策略

**时间离散：** 后向 Euler（隐式），保证稳定性

**空间离散：** 孔网络本身即为离散化

**非线性求解：** Newton-Raphson，每步求解稀疏线性系统

**耦合策略：** 强耦合（所有方程在同一个 Newton 迭代中联立求解）

**线性求解器：** scipy.sparse.linalg.spsolve（直接法）或 GMRES（迭代法）

### 2.5 验证标准

| 验证项 | 标准 |
|---|---|
| BV 动力学 | 单孔 i-η 曲线与解析公式一致 |
| 扩散方程 | 单颗粒 Cottrell 行为 |
| 低倍率 V-Q | 与热力学 OCV 偏差 < 5% |
| 1C V-Q | 与论文 Figure 5 一致（定性） |
| 3C 空间分布 | 电解液浓度梯度从隔膜到集流体递增 |
| 1CAL vs 3CAL | 3CAL 容量更高，SoL 更均匀 |

### 2.6 无 XCT 数据时的处理

论文使用的 NREL XCT 数据可能不再公开。备选方案：

1. **合成网络：** 用 OpenPNM 生成随机球堆积 + PoreSpy 提取网络
2. **参数化：** 调整合成网络的孔隙率、粒径、配位数以匹配论文 Table I 的统计量
3. **验证重点：** 物理方程的正确实现，而非精确匹配特定 XCT 数据

---

## 三、关键注意事项

### 3.1 单位一致性
- 所有长度单位: m（微米用 1e-6 m）
- 浓度: mol/m³
- 电位: V
- 电流密度: A/m²
- 扩散系数: m²/s
- 电导率: S/m

### 3.2 Butler-Volmer 数值稳定性
- 过电位 η 很大时，指数项会溢出 → 需要 clipping 或对数变换
- 建议: 当 |η| > 0.5V 时使用 Tafel 近似

### 3.3 NMC532 扩散系数公式
论文给出的 D_s 公式非常复杂，需要：
1. 确认公式中 cₛ 和 T 的单位
2. 验证 D_s 在合理范围内（~1e-14 到 ~1e-10 m²/s）
3. 绘制 D_s vs cₛ 曲线检查单调性

### 3.4 OCV 曲线
NMC532 的 U_eq(x) 曲线需要从文献获取（论文未给出具体公式）。常用经验公式：
```
U_eq(x) = a₀ + a₁·x + a₂·x² + a₃·x³ + a₄·x⁴ + a₅·x⁵
```
需要从文献查找 NMC532 的拟合系数。

### 3.5 开始实现前
- 先写 tests/ 里的单元测试（TDD）
- Phase 1 的单孔验证是最重要的：如果 BV + 扩散不对，后面全错
- 每完成一个 Phase 运行对应测试，通过后才进入下一 Phase

---

## 四、参考文献

1. Khan ZA, Elkamel A, Gostick JT. "Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes." J. Electrochem. Soc. 168(7), 070534 (2021). DOI: 10.1149/1945-7111/ac120c
2. OpenPNM 文档: https://openpnm.org/
3. PoreSpy 文档: https://porespy.org/
4. PMEAL lab publications: https://pmeal.com/publications.html

---

## 五、执行指令

请按 Phase 0 → 1 → 2 → 3 → 4 → 5 顺序严格实现。

每个 Phase 的要求：
1. 实现代码
2. 编写对应测试
3. 运行测试并确认通过
4. Git commit
5. 进入下一 Phase

**不要跳过任何 Phase。不要跳过任何测试。物理正确性优先于代码优雅。**
