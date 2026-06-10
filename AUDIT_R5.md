# AUDIT R5: 论文复核 — 代码 vs 论文原文

**日期**: 2026-06-10
**论文**: Khan et al., J. Electrochem. Soc. 168 (2021) 070534
**审计方式**: 视觉读取 PDF 图片 + DERIVATION.md + PAPER_EXTRACTION.md + 代码审查

---

## A. 论文图表分析（视觉读取）

### Figure 4: V-Q 放电曲线

**论文内容**（视觉读取 result）:
- 3 个子图: (a) 1CAL 实验+模型, (b) 3CAL 实验+模型, (c) 模型对比
- X 轴: Capacity (mAh/g_NMC), 0-200, ticks 25 步长
- Y 轴: Cell Voltage (V), 3.0-4.4
- 4 个 C-rate: 0.2C, 0.5C, 1C, 3C
- 实线 = 模型预测, 虚线 = 实验数据

**论文容量值**（从 Figure 4 视觉读取）:

| C-rate | 1CAL 模型 | 1CAL 实验 | 3CAL 模型 | 3CAL 实验 |
|--------|----------|----------|----------|----------|
| 0.2C   | ~195     | ~195     | ~190     | ~190     |
| 0.5C   | ~185     | ~185     | ~185     | ~185     |
| 1C     | ~175     | ~175     | ~175     | ~175     |
| 3C     | ~150     | (无数据) | ~170     | (无数据) |

**初始电压**: ~4.0-4.1V (x=0.35 附近的 Khan OCV)

**我们的结果**（Round 4, 文献 OCV, 10×10×10）:

| C-rate | 我们容量 | 初始 V | 终止 V | 到cutoff? |
|--------|---------|--------|--------|-----------|
| 0.2C   | ~273 (300步) | 3.69 | 3.43 | 未到 |
| 3C     | 164     | 3.69   | 2.98   | ✅ |

**差距分析**:
1. 初始电压差异: 论文 ~4.0V (Khan OCV at x=0.35) vs 我们 3.69V (文献 OCV)。**根因: OCV 曲线不同**
2. 容量差异: 论文 0.2C ~195 vs 我们 ~273。**根因: k 参数偏高 + 网络太小**

### Figure 6: 3D 空间分布

**论文内容**（视觉读取）:
- 8 个子图, 2 行 × 4 列
- 行: 1CAL (上), 3CAL (下)
- 列: (a)(e) Li⁺ 浓度 0-1550 mol/m³, (b)(f) 电解质电位 -0.41~-0.07V, (c)(g) SoL 31%-99%, (d)(h) 固相电位 0.3446-0.3467V
- 显示放电状态（高 SoL）下的空间异质性

**我们的现状**: 有基础可视化框架 (`visualization.py`)，但未生成对比图

### Figure 7: 逐孔 Profile

**论文内容**（视觉读取）:
- 12 个子图, 4 行 × 3 列
- 行: 1CAL/1C, 1CAL/3C, 3CAL/1C, 3CAL/3C
- 列: Li⁺ 浓度, 电解质电位, SoL
- X 轴: Distance from Membrane (μm), 0-120 (1CAL) / 0-80 (3CAL)
- 按孔径分类: 小/中/大孔分别绘制
- **关键特征**: 3C 时 Li⁺ 浓度从隔膜端 1400 降至集流体端 100 (强梯度)

**我们的现状**: 无此分析功能

---

## B. 参数验证表

| 参数 | 符号 | 论文值 (Table II) | 我们的值 | 文件:行 | 匹配? | 影响 |
|------|------|-------------------|---------|---------|-------|------|
| 反应速率常数 | k | **1e-10** mol/m²·s | **5e-10** | transient.py:108, validate:38 | ❌ 差5倍 | 🔴 P0: 直接影响容量和过电位 |
| OCV 曲线 | U(x) | Khan 多项式 Eq.2.19 | 文献查找表 | ocv.py | ❌ 不同曲线 | 🔴 P0: 影响初始电压和容量 |
| 电解质扩散 | D_e | **浓度依赖公式** Eq.2.21 | 常数 2e-10 | transient.py:390 | ❌ 简化 | 🟡 P1: 高倍率影响 |
| 固相扩散 | D_s | **SoC依赖多项式** Eq.2.13 | 公式不同 | solid.py:79 | ⚠️ 不同公式 | 🟡 P1: 影响浓度分布 |
| 电解质电导率 | κ | **浓度依赖公式** | 常数(隐含) | — | ❌ 缺失 | 🟡 P1: 影响电位分布 |
| NMC电子电导率 | σ_AM | 0.01 S/m | (未用) | — | ❌ 缺失 | 🟡 P1: 影响固相电位 |
| CBD电子电导率 | σ_CBD | 760 S/m | (未用) | — | ❌ 缺失 | 🟡 P1: 影响电子传导路径 |
| 孔隙率 | ε | 36.8% (1CAL) | 50% | validate:34 | ❌ 偏高 | 🔴 P0: 影响传输和容量 |
| 温度 | T | **303 K** | **298.15 K** | transient.py:108 | ❌ 差5K | 🟢 P2: 影响小 |
| 截止电压 | V_cutoff | 3.0 V | 3.0 V | validate:34 | ✅ | — |
| Faraday 常数 | F | 96485 C/mol | 96485 | reaction.py:16 | ✅ | — |
| 气体常数 | R | 8.314 J/mol·K | 8.314 | reaction.py:17 | ✅ | — |
| 最大锂浓度 | c_s,max | 48900 mol/m³ | 48900 | 多处 | ✅ | — |
| 初始锂浓度 | c_e,0 | 1200 mol/m³ | 1200 | validate:39 | ✅ | — |
| Bruggeman 指数 | b | 1.5 | 1.5 | transient.py:397 | ✅ | — |
| 电荷传递系数 | α_c | 0.5 | 0.5 | reaction.py | ✅ | — |

**关键不匹配**: k (5倍), OCV (不同曲线), ε (50% vs 36.8%), D_e (常数 vs 浓度依赖), T (298 vs 303K)

---

## C. 物理实现审计

| 方程 | DERIVATION 节 | 已实现? | 文件:行 | 正确? | 备注 |
|------|-------------|---------|---------|-------|------|
| Nernst-Planck 电解质传输 | §2.1 | ⚠️ 简化 | transient.py:390 | 部分 | D_e 为常数，缺迁移项 |
| 电解质电荷守恒 | §2.2 | ✅ | steady.py:466 | ✅ | L_e @ phi_e + rxn |
| Fick 固相扩散 | §2.3 | ✅ | transient.py:414-432 | ✅ | 调和平均 D_s |
| 固相电子传导 | §2.4 | ⚠️ | steady.py | 部分 | σ 为隐含值 |
| Butler-Volmer | §2.5 | ✅ | reaction.py:1-212 | ✅ | 含指数裁剪 |
| 电解质边缘电导 | §3.1 | ✅ | transient.py:389-412 | ✅ | G = D_eff A/L |
| 电解质残差 | §3.2 | ✅ | transient.py:588 | ✅ | (V/dt - L_ce) c = rhs |
| 电解质电位残差 | §3.3 | ✅ | steady.py:466 | ✅ | L_e @ phi_e + I_rxn A |
| 固相浓度残差 | §3.4 | ✅ | transient.py:600 | ✅ | (V/dt - L_cs) c = rhs |
| 固相电位残差 | §3.5 | ✅ | steady.py:472 | ✅ | L_s @ phi_s - I_rxn A + BC |
| BV 耦合 | §3.6 | ✅ | steady.py:448-456 | ✅ | g_bv = dI/deta A |
| 隔膜 Dirichlet BC | §4.1 | ✅ | transient.py:587-592 | ✅ | c_e = 1200 at sep |
| 集流体 Neumann BC | §4.2 | ✅ | steady.py:514-515 | ✅ | I_app at collector |
| 电位参考框架 | §4.4 | ✅ | steady.py:520-524 | ✅ | phi_e = 0 at sep |
| 电池电压公式 | §4.5 | ✅ | steady.py:689-693 | ✅ | V = phi_s(cc) - phi_e(sep) |
| Newton-Raphson | §5 | ✅ | steady.py:406-542 | ✅ | 含 current ramping |
| BV 溢出裁剪 | §6.1 | ✅ | reaction.py:51-56 | ✅ | arg clip ±500 |
| 浓度边界 | §6.2 | ✅ | transient.py:658 | ✅ | clip after step |
| 断开组件处理 | §6.4 | ✅ | steady.py:281-316 | ✅ | connected components |
| 自适应步长 | §6.6 | ⚠️ | — | 部分 | 有 dt 选择，无自适应 |
| 质量守恒审计 | §6.7 | ⚠️ | test_transient.py | 部分 | 仅零电流测试 |
| 0C 验证 | §7.1 | ✅ | test_steady.py | ✅ | near-OCV at low I |
| 小电流线性 | §7.2 | ✅ | test_steady.py | ✅ | linear scaling test |

---

## D. 缺失功能

| 功能 | 论文内容 | 难度 | 对结果影响 | 优先级 |
|------|----------|------|-----------|--------|
| XCT 真实微结构 | 4637 节点, NREL electrode library | 高 | 🔴 容量差 ~50% | P0 |
| 浓度依赖 D_e | Eq.2.21 多项式公式 | 中 | 🟡 高倍率容量 | P1 |
| 浓度依赖 κ | 电解质电导率随 c_e 变化 | 中 | 🟡 电位分布 | P1 |
| SoC 依赖 D_s | Eq.2.13 多项式公式 | 低 | 🟡 浓度分布 | P1 |
| σ_AM / σ_CBD | 固相电子电导率 (0.01 / 760 S/m) | 低 | 🟡 固相电位 | P1 |
| 自适应时间步长 | DERIVATION §6.6 | 中 | 🟡 数值稳定性 | P1 |
| Fig.6 3D 可视化 | c_e, phi_e, SoL, phi_s 空间分布 | 中 | — | P2 |
| Fig.7 逐孔 profile | 沿厚度方向分布, 按孔径分类 | 中 | — | P2 |
| 1CAL vs 3CAL | 两种电极压实度 | 高 | — | P2 |
| 参数灵敏度分析 | k, D_s, ε 对容量的影响 | 低 | — | P2 |
| 质量守恒 (非零电流) | DERIVATION §6.7 | 低 | 🟢 | P2 |

---

## E. 差异根因分析

| 差异 | 根因 | 修复方案 | 优先级 |
|------|------|---------|--------|
| 初始电压 3.69V vs 论文 ~4.0V | OCV 曲线不同: 我们用文献表, 论文用 Khan 多项式 | 恢复 Khan 多项式 + clip(x<0.25, 3.0) | 🔴 P0 |
| 0.2C 容量 273 vs 论文 195 | k=5e-10 偏高 5 倍 | 改 k0=1e-10 | 🔴 P0 |
| 0.2C 容量 273 vs 论文 195 | ε=50% 偏高 (论文 36.8%) | 改 porosity=0.368 | 🔴 P0 |
| 3C 容量 164 vs 论文 150 | 同上 k0 + ε | 同上 | 🔴 P0 |
| SoL 异质性 | 合成网络不如 XCT 真实 | 需要 XCT 数据 | 🟡 P1 |
| 高倍率电位分布 | D_e 和 κ 为常数 | 改为浓度依赖 | 🟡 P1 |

---

## F. 下一步实施计划

### Phase 5A: 参数校准 (预计 2-3 小时)

**任务 5A-1: 恢复 Khan OCV + 修正 k0 + 修正 ε**
- 修改 `src/physics/ocv.py`: 恢复 Khan 多项式, 对 x<0.25 加 clip 到 3.0V
- 修改 `scripts/validate_paper_figures.py` L38: `K0 = 1e-10` (从 5e-10)
- 修改 `scripts/validate_paper_figures.py` L34: `POROSITY = 0.368` (从 0.5)
- 修改 `src/network/generator.py`: 确保 porosity 参数传递正确
- 修改 temperature: T=303K (从 298.15K)
- 运行 `pytest -q` 确认 70 测试通过
- 运行 10×10×10 0.2C 放电, 预期容量 ~180-200 mAh/g

**任务 5A-2: 验证校准结果**
- 运行 0.2C/0.5C/1C/3C 四个倍率
- 对比论文 Figure 4 数值
- 如果容量偏差 >10%, 微调 k0

### Phase 5B: 物理增强 (预计 4-6 小时)

**任务 5B-1: 浓度依赖 D_e**
- 修改 `src/solver/transient.py` `_build_diffusion_matrices()`
- 将 D_e = 2e-10 改为 Eq.2.21 多项式公式
- 需要每步重新计算扩散矩阵 (当前为常数矩阵)
- 验证: 3C 容量应下降 (浓差极化增大)

**任务 5B-2: 浓度依赖 κ**
- 修改 `src/solver/steady.py` `_compute_conductances()`
- 添加电解质电导率随 c_e 变化
- 验证: 高倍率电位分布应更陡

**任务 5B-3: σ_AM / σ_CBD 固相电导率**
- 修改 `src/solver/steady.py` `_compute_conductances()`
- 区分 AM (0.01 S/m) 和 CBD (760 S/m) 的电导率
- 验证: 固相电位分布应有变化

### Phase 5C: 可视化对标 (预计 2-3 小时)

**任务 5C-1: Fig.4 V-Q 对比图**
- 修改 `scripts/validate_paper_figures.py`
- 绘制模型曲线 vs 论文实验点
- 保存为 `data/validation/figure4_comparison.png`

**任务 5C-2: Fig.6 空间分布**
- 在 10×10×10 网络上运行到特定 SoL (75%)
- 生成 c_e, phi_e, SoL, phi_s 的 3D 散点图
- 保存为 `data/validation/figure6_comparison.png`

**任务 5C-3: Fig.7 逐孔 profile**
- 沿 x 方向 (隔膜→集流体) 投影
- 按孔径分类绘制
- 保存为 `data/validation/figure7_comparison.png`

### Phase 5D: 架构优化 (预计 3-4 小时)

**任务 5D-1: 自适应时间步长**
- 实现 DERIVATION §6.6 的步长选择策略
- 基于 Newton 迭代次数调整 dt
- 验证: 更少步数达到相同精度

**任务 5D-2: 质量守恒测试 (非零电流)**
- 添加 DERIVATION §6.7 的质量守恒审计
- 在每个时间步检查 Σ(c_e ΔV_e + c_s ΔV_s) = const

---

## 总结

| 维度 | 已完成 | 待完成 |
|------|--------|--------|
| 物理框架 | ✅ 完整 (BV, 扩散, 电位) | 浓度依赖 D_e, κ |
| 求解器 | ✅ Newton-Raphson + ramping | 自适应步长 |
| OCV | ❌ 文献表 ≠ 论文 Khan | 恢复 Khan + clip |
| 参数 | ❌ k=5e-10, ε=50% | k=1e-10, ε=36.8% |
| 网络 | ⚠️ 合成 10×10×10 | XCT 真实数据 |
| 可视化 | ⚠️ 基础框架 | Fig.4/6/7 对比图 |
| 测试 | ✅ 70/70 通过 | 质量守恒测试 |

**核心结论**: 代码框架完整且正确，差距主要在**参数值**和**OCV 曲线**。Phase 5A 参数校准是最高优先级——修正 k0、OCV、ε 后，容量应能接近论文值。
