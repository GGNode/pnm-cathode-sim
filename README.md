# pnmcathode

**锂离子电池阴极孔网络放电仿真工具包**

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

复现论文: Khan ZA, Elkamel A, Gostick JT. *"Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes."* J. Electrochem. Soc. 168(7), 070534 (2021).

---

## 项目概述

`pnmcathode` 是一个基于孔网络模型 (PNM) 的锂电池阴极放电仿真工具包。它将多孔电极离散化为孔隙节点与喉道键的图结构，耦合以下物理过程：

- **电解质盐传输** — 液相中 Li⁺ 的扩散和迁移
- **固相锂扩散** — 球形活性材料颗粒中的 Fick 扩散（浓度依赖扩散系数）
- **Butler-Volmer 动力学** — 电解质/活性材料界面的电化学反应
- **电荷守恒** — 电解质和固相中的欧姆传导（耦合 Newton-Raphson 求解）
- **隔膜边界模型** — 一维折叠隔膜模型，含欧姆损失、浓度过电位、Li 箔 BV

## 快速开始

```bash
git clone https://github.com/GGNode/pnm-cathode-sim.git
cd pnm-cathode-sim
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

```python
from pnmcathode import Cathode, DischargeProtocol, Simulation
from pnmcathode.materials.presets import nmc532_khan2021, electrolyte_khan2021

# 创建阴极 + 运行 1C 放电
cathode = Cathode.cubic(
    active_material=nmc532_khan2021(),
    electrolyte=electrolyte_khan2021(),
)
result = Simulation(cathode, DischargeProtocol(c_rate=1.0)).run()

print(f"容量: {result.final_capacity:.4f} Ah/m²")
print(f"终止电压: {result.voltage[-1]:.3f} V")
```

## 包结构

```
pnmcathode/
├── config.py              # 配置 dataclass（CathodeGeometry, ActiveMaterial, Electrolyte 等）
├── cathode.py             # Cathode 类（Cathode.cubic() 工厂方法）
├── simulation.py          # Simulation 类（一站式仿真运行）
├── results.py             # DischargeResult（结果容器，支持 npz 序列化）
├── analysis.py            # 分析工具
├── plotting.py            # 绘图工具
├── materials/
│   └── presets.py         # 预设材料参数（Khan 2021 NMC532）
├── network/
│   └── generator.py       # 立方孔网络生成器 + 三相标记 + 渗流检查
├── physics/
│   ├── electrolyte.py     # D_e(c_e,T) 和 κ(c_e,T) 物性相关性
│   ├── solid.py           # D_s(c_s,T) 固相扩散系数
│   ├── reaction.py        # Butler-Volmer 动力学 + 交换电流密度
│   ├── ocv.py             # NMC532 开路电位（Khan 多项式）
│   └── separator.py       # 折叠一维隔膜边界模型
├── solver/
│   ├── steady.py          # 稳态电位求解器（耦合 Newton-Raphson）
│   ├── transient.py       # 瞬态放电求解器（自适应时间步进）
│   └── single_pore.py     # 单孔参考模型
└── post/
    ├── analysis.py        # 后处理分析
    └── visualization.py   # 可视化工具
```

## 使用指南

### 基本放电

```python
from pnmcathode import Cathode, DischargeProtocol, Simulation
from pnmcathode.materials.presets import nmc532_khan2021, electrolyte_khan2021

cathode = Cathode.cubic(
    active_material=nmc532_khan2021(),
    electrolyte=electrolyte_khan2021(),
)
result = Simulation(cathode, DischargeProtocol(c_rate=1.0, cutoff_voltage=2.5)).run()
```

### 自定义阴极几何

```python
from pnmcathode.config import CathodeGeometry

geometry = CathodeGeometry(
    shape=(13, 13, 13),       # 网格尺寸
    spacing=1e-5,             # 节点间距 [m]
    porosity=0.368,           # 孔隙率
    cbd_fraction=0.1392,      # CBD 体积分数
    seed=42,                  # 随机种子（可复现）
)
cathode = Cathode.cubic(geometry=geometry, ...)
```

### 自定义材料

```python
from pnmcathode.config import ActiveMaterial, Electrolyte

my_material = ActiveMaterial(
    name="LFP",
    cs_max=22806.0,           # 最大锂浓度 [mol/m³]
    sigma=0.1,                # 电子导电率 [S/m]
    diffusivity=my_D_s,       # D_s(c_s, T) 函数
    ocv=my_ocv_curve,         # OCV(soc) 函数
)
```

### 多倍率扫描

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
for c_rate in [0.2, 0.5, 1.0, 3.0]:
    result = Simulation(cathode, DischargeProtocol(c_rate=c_rate)).run()
    ax.plot(result.capacity_Ah_m2, result.voltage, label=f"{c_rate}C")
ax.set_xlabel("容量 [Ah/m²]")
ax.set_ylabel("电压 [V]")
ax.legend()
plt.show()
```

### 隔膜模型

```python
from pnmcathode.config import Separator

separator = Separator(enabled=True, thickness=25e-6, porosity=0.39)
result = Simulation(cathode, protocol, separator=separator).run()
```

### 结果保存与加载

```python
result.to_npz("discharge_1c.npz")
loaded = DischargeResult.from_npz("discharge_1c.npz")
```

## 物理模型概述

| 模型 | 方程 | 代码位置 |
|------|------|---------|
| 电解质扩散 | D_Li⁺ = 10^(-4.43 - 54/(T-229-5c) - 0.22c) | `physics/electrolyte.py` |
| 电解质电导 | κ = c·(-10.5 + 0.074T - ...)² | `physics/electrolyte.py` |
| 固相扩散 | D_s = 10^(多项式(soc)) | `physics/solid.py` |
| OCV | Khan 9阶多项式 + exp修正 | `physics/ocv.py` |
| 反应动力学 | Butler-Volmer: i_r = i0·[exp(α_a·F·η/RT) - exp(-α_c·F·η/RT)] | `physics/reaction.py` |
| 隔膜 | 1D 折叠边界：Ohmic + Nernst + Li箔 BV | `physics/separator.py` |
| 电位求解 | 耦合 Newton-Raphson（phi_e + phi_s） | `solver/steady.py` |
| 时间推进 | Backward Euler + 自适应步长 | `solver/transient.py` |

## 配置参数参考

### CathodeGeometry

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `shape` | (10,10,10) | 网格尺寸 [Nx, Ny, Nz] |
| `spacing` | 1e-5 | 节点间距 [m] |
| `porosity` | 0.35 | 孔隙率（电解质体积分数） |
| `cbd_fraction` | 0.10 | CBD 体积分数 |
| `seed` | None | 随机种子 |
| `throat_scale` | 1.0 | 喉道直径缩放因子 |

### DischargeProtocol

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `c_rate` | — | 放电倍率（必填） |
| `cutoff_voltage` | 2.5 | 截止电压 [V] |
| `initial_soc` | 0.5 | 初始荷电状态 |
| `max_steps` | 1000 | 最大时间步数 |

### SolverSettings

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `temperature` | 298.15 | 温度 [K] |
| `dt` | None | 初始时间步（None=自动估计） |
| `dt_max` | 300.0 | 最大时间步 [s] |
| `newton_tol` | 1e-8 | Newton 迭代收敛容差 |

## 测试

```bash
# 运行全部测试
python -m pytest tests/ -q

# 只运行 API 测试
python -m pytest tests/test_toolkit_api.py -v

# 运行物理测试
python -m pytest tests/test_separator.py tests/test_electrolyte.py -v
```

测试覆盖：
- Butler-Volmer 动力学符号和对称性
- 电解质物性函数（D_e, κ）
- 固相扩散系数
- OCV 曲线和导数
- 网络生成和渗流检查
- 稳态求解器收敛
- 瞬态放电物理守恒
- 隔膜边界模型
- 包 API 端到端测试

## 论文复现

论文复现脚本位于 `examples/reproduce_khan2021/`：

```bash
cd examples/reproduce_khan2021
python validate_paper_figures.py    # 生成 Figure 4/6/7 对比图
python parallel_validation.py       # 并行多倍率验证
```

## 已知限制

| 限制 | 说明 |
|------|------|
| 合成网络 | 使用随机立方网络而非 XCT 真实微结构 |
| 质量负载 | 合成网络质量负载低于论文 XCT 网络（~150 vs 298 g/m²） |
| 低倍率截止 | 合成网络在低 C-rate 下可能无法达到截止电压 |
| 隔膜模型 | 简化的一维折叠边界，非完整耦合隔膜网格 |

## 相关文档

- `docs/TOOLKIT_DESIGN.md` — 工具包架构设计
- `docs/TOOLKIT_EXECUTION.md` — 迁移执行方案
- `docs/AUDIT_FINAL.md` — 最终代码审计
- `DERIVATION.md` — 数学推导
- `PAPER_REFERENCE.md` — 论文数据提取

## 引用

```bibtex
@article{khan2021pore,
  title={Pore Network Modelling of Galvanostatic Discharge Behaviour of Lithium-Ion Battery Cathodes},
  author={Khan, Zaid A and Elkamel, Ali and Gostick, Jeffrey T},
  journal={Journal of The Electrochemical Society},
  volume={168},
  number={7},
  pages={070534},
  year={2021},
  publisher={IOP Publishing}
}
```
