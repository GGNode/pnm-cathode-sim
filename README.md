# PNM-LIB-Cathode

**锂离子电池阴极孔隙网络放电模型**

复现论文: Khan ZA, Elkamel A, Gostick JT. *"Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes."* J. Electrochem. Soc. 168(7), 070534 (2021).

## 项目概述

本项目实现了一个瞬态孔隙网络模型 (PNM), 用于模拟 NMC532 锂离子电池阴极的恒流放电行为。该模型将多孔电极离散化为孔隙节点 (存储物质浓度和电位) 与喉道键 (传导通量) 的图结构, 并耦合以下物理过程:

- **电解质盐传输** — 液相中 Li⁺ 的扩散和迁移
- **固相锂扩散** — 球形 NMC532 颗粒中的 Fick 扩散
- **Butler-Volmer 动力学** — 电解质/NMC 界面的电化学反应
- **电荷守恒** — 电解质和固相中的欧姆传导
- **隔膜边界模型** — 折叠的一维隔膜模型, 含欧姆/浓度损失

## 论文引用

> Khan ZA, Elkamel A, Gostick JT. "Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes." *J. Electrochem. Soc.* 168(7), 070534 (2021).

## 目录结构

```
pnm-lib-cathode/
├── src/                          # 核心库
│   ├── __init__.py
│   ├── network/
│   │   └── generator.py          # 立方孔网络 + 三相标记
│   ├── physics/
│   │   ├── electrolyte.py        # D_e(c_e,T) 和 kappa(c_e,T) 物性相关性
│   │   ├── ocv.py                # NMC532 OCV 多项式 (Khan Eq. 2.19)
│   │   ├── reaction.py           # Butler-Volmer 动力学与交换电流密度
│   │   ├── separator.py          # 折叠的一维隔膜边界模型
│   │   └── solid.py              # 球形固相扩散 (有限体积壳层法)
│   ├── solver/
│   │   ├── steady.py             # 耦合 phi_e/phi_s Newton-Raphson 求解器
│   │   ├── transient.py          # 自适应后向欧拉放电求解器
│   │   └── single_pore.py        # 单孔放电仿真 (验证用)
│   └── post/
│       ├── analysis.py           # 孔径分布、锂化度、放电曲线
│       └── visualization.py      # Matplotlib 绘图工具
├── scripts/
│   ├── run_discharge.py          # 单次放电 CLI 入口
│   ├── validate_paper_figures.py # 论文图表对比 (Figure 4, 6, 7)
│   ├── parallel_validation.py    # 多 C-rate 并行验证
│   └── plot_results.py           # 快速 V-Q 绘图工具
├── tests/                        # pytest 测试套件 (88 个测试)
├── data/
│   ├── paper_figures/            # 论文裁剪图 PNG
│   └── validation/               # 预计算验证结果
├── DERIVATION.md                 # 完整数学推导
├── PAPER_REFERENCE.md            # 论文数据提取
└── AUDIT_FINAL.md                # 最终审计报告
```

## 安装指南

```bash
# 创建并激活虚拟环境
python -m venv .venv
source .venv/bin/activate

# 安装依赖
pip install -r requirements.txt
```

**依赖:** numpy, scipy, matplotlib, openpnm (>=4.0), pytest

## 快速开始

在 5x5x5 网络上运行 0.2C 放电, 仅需 5 行代码:

```python
from src.network.generator import create_cathode_network
from src.solver.transient import TransientSolver

net = create_cathode_network(shape=[5, 5, 5], spacing=1e-5, porosity=0.35, seed=42)
solver = TransientSolver(net, T=298.15, k0=5e-10)
result = solver.run_discharge(C_rate=0.2, cutoff_voltage=3.0)

print(f"Capacity: {result['capacity_Ah_m2'][-1]:.3f} Ah/m²")
print(f"Final voltage: {result['voltage'][-1]:.3f} V")
```

或通过命令行运行:

```bash
python scripts/run_discharge.py --crate 0.2 --shape 5x5x5 --cutoff 3.0
```

## 使用指南

### 单次放电 (`scripts/run_discharge.py`)

```bash
python scripts/run_discharge.py --crate 1.0 --shape 10x10x10 --cutoff 2.5 --output data/1c.npz
```

主要参数: `--crate`, `--cutoff`, `--shape`, `--spacing`, `--porosity`, `--k0`, `--temperature`, `--output`

### 论文验证 (`scripts/validate_paper_figures.py`)

在 10x10x10 网络上运行 0.2C/0.5C/1C/1C/3C 放电, 生成与论文 Figure 4、6、7 的对比图:

```bash
python scripts/validate_paper_figures.py
```

### 并行验证 (`scripts/parallel_validation.py`)

在 13x13x13 网络上并行运行 4 个 C-rate x 4 种场景 (论文/网络电流基准, 隔膜开/关):

```bash
python scripts/parallel_validation.py
```

## 物理模型概述

### 控制方程

| 方程 | 描述 | 参考 |
|------|------|------|
| `∂c_e/∂t = ∇·(D_e,eff ∇c_e) + S_e` | 电解质盐守恒 | §3.2 |
| `∇·(κ_eff ∇φ_e) = a_s·i_F` | 电解质电荷守恒 | §2.2, §3.3 |
| `∂c_s/∂t = ∇·(D_s ∇c_s)` | 固相锂扩散 (球坐标) | §2.3, §3.4 |
| `∇·(σ_eff ∇φ_s) = -a_s·i_F` | 固相电荷守恒 | §2.4, §3.5 |
| `i_r = i0·[exp(α_a f η) - exp(-α_c f η)]` | Butler-Volmer 动力学 | §2.5 |
| `η = φ_s - φ_e - U_eq(c_s)` | 过电位 | §2.5 |

### 关键材料参数 (Khan et al. Table II)

| 参数 | 符号 | 值 | 单位 |
|------|------|----|------|
| NMC532 最大锂浓度 | c_s,max | 48,900 | mol/m³ |
| NMC 电导率 | σ_nmc | 0.01 | S/m |
| CBD 电导率 | σ_cbd | 760 | S/m |
| 反应速率常数 | k0 | 5×10⁻¹⁰ | m²·⁵/(mol⁰·⁵·s) |
| 交换电流密度 | i0 | ~10⁻³ | A/m² |
| NMC 中锂扩散系数 | D_s | ~10⁻¹⁵ | m²/s |
| 电解质扩散系数 | D_e | ~10⁻¹⁰ | m²/s |

### 离散化方案

- **网络**: 立方晶格 + 三相随机标记 (电解质/NMC/CBD)
- **固相扩散**: 球形颗粒的 N 壳层有限体积法
- **电位求解**: 耦合 Newton-Raphson + 回溯线搜索
- **时间推进**: 后向欧拉 + 自适应步长 (增长/收缩/拒绝)

## 配置参数参考

所有可调参数 (含默认值):

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `shape` | [5,5,5] | 网络维度 [nx, ny, nz] |
| `spacing` | 1e-5 m | 孔间距 |
| `porosity` | 0.35 | 电解质体积分数 |
| `cbd_fraction` | 0.10 | CBD 体积分数 |
| `T` | 298.15 K | 温度 |
| `k0` | 5e-10 | BV 反应速率常数 |
| `c_e_init` | 1200 mol/m³ | 初始电解质浓度 |
| `c_s_init` | 24450 mol/m³ | 初始固相浓度 (x≈0.5) |
| `cutoff_voltage` | 2.5 V | 放电截止电压 |
| `dt_min` | 1e-6 s | 最小自适应时间步 |
| `dt_max` | 300 s | 最大自适应时间步 |

## 测试套件

```bash
# 运行全部测试
python -m pytest tests/ -q

# 详细输出
python -m pytest tests/ -v

# 运行指定模块
python -m pytest tests/test_transient.py -v
```

**测试覆盖 (88 个测试):**
- Butler-Volmer 动力学与交换电流密度
- OCV 多项式及其导数
- 电解质传输物性相关性
- 固相扩散 (球坐标有限体积法)
- 网络生成与渗透性检查
- 稳态 Newton-Raphson 求解器
- 瞬态求解器自适应步进
- 单孔放电验证
- 隔膜边界模型
- 后处理分析工具

## 已知限制

1. **合成网络 vs XCT**: 论文使用 XCT 重建的阴极微结构 (1CAL: 4637 节点, 31427 喉道)。本模型使用合成立方网络和随机相标记。几何结构差异显著, 导致:
   - 不同的质量负载和活性面积
   - 不同的相连通性和渗透率
   - 定量的容量/极化不匹配 (尤其在高 C-rate 下)

2. **质量负载差距**: 合成网络的面积比容量低于 XCT 网络, 原因是体积-面积比和活性位点连通性不同。

3. **隔膜模型**: 折叠的一维准稳态模型 — 适用于毫伏级损失, 但非完整的二维/三维隔膜仿真。

4. **等温假设**: 温度均匀恒定 (无热耦合)。

5. **无退化模型**: 未建模 SEI 生长、颗粒开裂和容量衰减。

## 引用

如果使用本代码, 请引用:

> Khan ZA, Elkamel A, Gostick JT. "Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes." *J. Electrochem. Soc.* 168(7), 070534 (2021).

## 参考文献

- [DERIVATION.md](DERIVATION.md) — 完整数学推导
- [PAPER_REFERENCE.md](PAPER_REFERENCE.md) — 论文数据提取
- [AUDIT_FINAL.md](AUDIT_FINAL.md) — 最终审计报告
