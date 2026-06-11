# P2D+PNM Phase 1 验证任务书

本文档把 `PLAN.md` 中 Phase 1（P2D 验证）细化为可执行任务。Phase 1 的目标不是继续扩展模型，而是用可信参考数据确认当前 `src/pnmcathode/p2d/` 的 half-cell DFN/P2D solver 在符号、单位、边界条件、守恒和曲线输出上是物理正确的。

当前实现范围为 `Li metal reference | separator | porous NMC532 cathode`，状态变量为

\[
y = \left[c_e(x),\ \phi_e(x),\ \phi_s(x),\ c_s(x,r)\right],
\]

其中 \(c_s\) 和 \(\phi_s\) 只定义在 positive electrode cells。项目采用阳极约定：放电时

\[
I_{\mathrm{app}} < 0,\qquad i_F < 0,\qquad \bar c_s \ \text{增加}.
\]

## 1. 验证目标

Phase 1 要确认 P2D 的以下物理行为正确：

1. **零电流平衡态**
   - 当 \(I_{\mathrm{app}}=0\) 且 \(c_e,c_s\) 均匀时，过电位为零：
     \[
     \eta = \phi_s - \phi_e - U(c_s/c_{s,\max}) = 0.
     \]
   - Butler-Volmer 反应电流满足 \(i_F=0\)。
   - \(c_e\)、\(c_s\)、\(\phi_s-\phi_e\) 不随时间漂移。

2. **OCV 一致性**
   - half-cell 电压定义为
     \[
     V(t) = \phi_s(x=L,t) - \phi_e(x=0,t).
     \]
   - 在当前 gauge \(\phi_e(0)=0\) 下，`P2DSolver.voltage(state)` 应等于正极集流体侧 `phi_s[-1]`。
   - 零电流时：
     \[
     V = U(\mathrm{SOC})
     \]
     且使用的 OCV 曲线必须来自同一参数源。

3. **电流和电荷守恒**
   - 对每个内部 face：
     \[
     i_e + i_s = I_{\mathrm{app}}.
     \]
   - positive electrode 中：
     \[
     \frac{\partial i_e}{\partial x} = a_s i_F,\qquad
     \frac{\partial i_s}{\partial x} = -a_s i_F.
     \]
   - separator 中无反应，\(i_e\) 应承载全部电流。
   - 注意：当前 `compute_fluxes()` 对左边界令 `i_e[0]=0`。Phase 1 的非零电流边界测试应捕捉并驱动修正：Li metal/separator 边界处应为 ionic current 或 Li foil BV，而不是无电流边界。

4. **电解质浓度分布**
   - 放电时正极消耗电解质 Li salt source：
     \[
     \frac{\partial(\epsilon_e c_e)}{\partial t}
     = \nabla\cdot(D_e^{\mathrm{eff}}\nabla c_e)
       + (1-t_+^0)\frac{a_s i_F}{F}.
     \]
   - 因 \(i_F<0\)，高倍率时 cathode 内部 \(c_e\) 应相对初始值下降，且从 separator 端到 current collector 端出现明显梯度。
   - 与 Khan 2021 Figure 7 趋势一致：1CAL/3CAL 在 3C 时均出现 \(c_e(x)\) 梯度，厚电极更强。

5. **固相浓度和颗粒扩散**
   - 放电时颗粒平均固相浓度增加：
     \[
     \epsilon_s\frac{\partial \bar c_s}{\partial t}
     = -\frac{a_s i_F}{F}.
     \]
   - 表面浓度 \(c_s^{\mathrm{surf}}\) 和平均浓度 \(\bar c_s\) 应满足球形颗粒扩散的方向性：放电时表面先升高，随后向颗粒内部扩散。
   - 高倍率时 \(c_s^{\mathrm{surf}}-\bar c_s\) 的幅值应大于低倍率。

6. **电压曲线和倍率趋势**
   - 同一初始 SOC 下，倍率越高，初始 IR drop 越大，运行电压越低：
     \[
     V_{0.2C}(Q) > V_{0.5C}(Q) > V_{1C}(Q) > V_{3C}(Q)
     \]
     在共同容量区间内应大体成立。
   - 到达 cutoff voltage 的容量应随倍率升高而不增加；若没有达到 cutoff，则必须在结果 metadata 中标记。

## 2. 对标数据来源

Phase 1 支持两类对标数据，不能混用参数。

### 2.1 Khan et al. 2021 数据

项目已有：

- `data/khan2021_pnm_lib_cathode.pdf`
- `data/paper.pdf`
- `PAPER_REFERENCE.md`
- `data/paper_figures/`
- `data/validation/*.npz`

Khan 2021 作为当前 half-cell/cathode 验证的首选来源，因为项目的材料模型和 PNM 工作流已经围绕该论文建立。必须对齐的参数：

| 参数 | Khan 2021 值 | 项目入口 |
|---|---:|---|
| \(L_s\) | \(25\,\mu\mathrm{m}\) | `separator_khan2021()` |
| \(\epsilon_{\mathrm{sep}}\) | 0.39 | `separator_khan2021()` |
| \(c_{e,0}\) | \(1200\,\mathrm{mol\,m^{-3}}\) | `electrolyte_khan2021().c_init` |
| \(c_{s,\max}\) | \(48900\,\mathrm{mol\,m^{-3}}\) | `nmc532_khan2021().cs_max` |
| \(t_+^0\) | 0.363 | `Separator.t_plus` |
| \(T\) | \(303\,\mathrm{K}\) | `SolverSettings(temperature=303.0)` |
| \(\sigma_{\mathrm{AM}}\) | \(0.01\,\mathrm{S\,m^{-1}}\) | `nmc532_khan2021().sigma` |
| \(\sigma_{\mathrm{CBD}}\) | \(760\,\mathrm{S\,m^{-1}}\) | `cbd_khan2021()` |
| \(k\) | \(1\times10^{-10}\) | `Kinetics(k0=1e-10)` |
| cutoff | \(3.0\,\mathrm{V}\) | `P2DProtocol.cutoff_voltage=3.0` |
| 1CAL thickness | \(129\,\mu\mathrm{m}\) | `P2DRegion.x_right-x_left` |
| 3CAL thickness | \(88\,\mu\mathrm{m}\) | `P2DRegion.x_right-x_left` |
| 1CAL porosity | 0.368 | `P2DRegion.epsilon_e` |
| 3CAL porosity | 0.366 | `P2DRegion.epsilon_e` |
| 1CAL active fraction | 0.4928 | `P2DRegion.epsilon_s` |
| 3CAL active fraction | 0.4974 | `P2DRegion.epsilon_s` |

Khan 2021 的参考用途分两层：

1. **定性/趋势验证**：用 `PAPER_REFERENCE.md` 中 Figure 4/6/7 的视觉读取数据确认倍率趋势、厚度趋势和空间剖面趋势。
2. **项目内回归验证**：用 `data/validation/discharge_*.npz`、`spatial_*.npz` 作为已经抽取或生成的参考数据，避免 CI 依赖 PDF 解析。

### 2.2 PyBaMM DFN 默认参数

PyBaMM 当前不是项目依赖，也不在本环境中安装。Phase 1 应把 PyBaMM 对标做成可选生成步骤：

- `examples/p2d/compare_pybamm.py`：安装 PyBaMM 的开发环境中运行，生成 `data/validation/p2d_pybamm_*.npz`。
- `tests/test_p2d_against_reference.py`：CI 默认读取已经生成的 `.npz` 参考文件；若文件不存在则 `pytest.skip()`。

不要依赖 PyBaMM 隐式默认参数随版本变化。脚本必须把参数集名称、PyBaMM 版本和所有关键标量写入 metadata。推荐优先使用显式参数集：

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class PyBaMMReferenceSpec:
    parameter_set: str  # e.g. "Chen2020" or installed default recorded at runtime
    model_kind: Literal["DFN", "positive-half-cell"]
    c_rate: float
    initial_soc: float
    temperature: float
    output_path: Path
```

若 PyBaMM 支持 positive working electrode/half-cell DFN，则优先与当前 solver 对齐；否则只把 full-cell DFN 用作趋势基准，并在比较时明确扣除或记录 negative electrode 贡献，不能把 full-cell 电压误差作为当前 half-cell solver 的硬性失败标准。

## 3. 验证测试清单与具体断言

### 3.1 参数快照测试

新增 `tests/test_p2d_parameter_alignment.py`。

断言：

- Khan case 中 `temperature == 303.0`。
- `active.cs_max == 48900.0`。
- `electrolyte.c_init == 1200.0`。
- `separator.thickness == 25e-6`。
- `separator.t_plus == 0.363`。
- `kinetics.k0 == 1e-10`。
- 1CAL:
  \[
  L_p=129\,\mu\mathrm{m},\quad \epsilon_e=0.368,\quad \epsilon_s=0.4928.
  \]
- 3CAL:
  \[
  L_p=88\,\mu\mathrm{m},\quad \epsilon_e=0.366,\quad \epsilon_s=0.4974.
  \]
- 每次 benchmark 输出 `metadata["parameter_source"] in {"khan2021", "pybamm"}`。

### 3.2 零电流平衡态

新增或扩展 `tests/test_p2d_half_cell.py`。

断言：

- 对 SOC = 0.35, 0.5, 0.75 各跑 10 个步长：
  \[
  \max|c_e(t)-c_e(0)|/c_{e,0} \le 10^{-10}.
  \]
  \[
  \max|c_s(t)-c_s(0)|/c_{s,\max} \le 10^{-10}.
  \]
  \[
  \max|\eta| \le 10^{-8}\ \mathrm{V}.
  \]
  \[
  \max|i_F| \le 10^{-9}\ \mathrm{A\,m^{-2}}.
  \]
- 电压：
  \[
  |V-U(\mathrm{SOC})| \le 10^{-8}\ \mathrm{V}.
  \]

### 3.3 OCV 曲线一致性

扩展 `tests/test_ocv.py` 或新增 `tests/test_p2d_ocv_alignment.py`。

断言：

- `from_config()` 构造的 `P2DMaterial.ocv_model.value()` 与 `nmc532_ocv()` 在同一 SOC 网格上完全一致：
  \[
  \max |U_{\mathrm{p2d}} - U_{\mathrm{khan}}| \le 10^{-12}\ \mathrm{V}.
  \]
- OCV derivative 与中心差分一致，避开 clip 区间和极端 SOC：
  \[
  \max |U'(x)-U'_{\mathrm{fd}}(x)| \le 10^{-4}\ \mathrm{V}.
  \]
- 对 PyBaMM 参数 case，必须禁用 Khan OCV，改用 PyBaMM 导出的 positive electrode OCP 插值；否则测试失败。

### 3.4 非零电流边界条件和总电流守恒

新增 `tests/test_p2d_boundary_current.py`。

对 \(I_{\mathrm{app}}=-1,-10,-50\,\mathrm{A\,m^{-2}}\) 各跑一个收敛时间步后，用 `compute_fluxes()` 提取 faces：

```python
from __future__ import annotations

import numpy as np

from pnmcathode.p2d.residual import P2DResidualContext, compute_fluxes
from pnmcathode.p2d.state import P2DState


def total_current_error(state: P2DState, context: P2DResidualContext) -> float:
    fluxes = compute_fluxes(state, context)
    i_total = fluxes.ionic_current_faces + fluxes.solid_current_faces
    interior = slice(1, context.layout.n_x)
    scale = max(abs(context.current_density), 1.0)
    return float(np.max(np.abs(i_total[interior] - context.current_density)) / scale)
```

断言：

- 内部 faces：
  \[
  \frac{\max|i_e+i_s-I_{\mathrm{app}}|}{\max(|I_{\mathrm{app}}|,1)} \le 10^{-8}.
  \]
- positive/separator interface 的 solid current：
  \[
  |i_s| \le 10^{-10}\ \mathrm{A\,m^{-2}}.
  \]
- separator 内部：
  \[
  \max |i_e-I_{\mathrm{app}}|/\max(|I_{\mathrm{app}}|,1) \le 10^{-8}.
  \]
- cathode current collector：
  \[
  |i_s(L)-I_{\mathrm{app}}| \le 10^{-10}\ \mathrm{A\,m^{-2}}.
  \]

这组测试是 Phase 1 的硬门槛。若当前实现失败，应先修正 half-cell Li metal 边界，而不是放宽阈值。

### 3.5 质量守恒

扩展 `tests/test_p2d_conservation.py`。

对每个 accepted time step 检查全局 balance：

\[
\frac{dM_e}{dt}
= \int_{\Omega_p} (1-t_+^0)\frac{a_s i_F}{F}\,dV
  - \int_{\partial\Omega} N_e\cdot n\,dA,
\]

\[
\frac{dM_s}{dt}
= -\int_{\Omega_p}\frac{a_s i_F}{F}\,dV.
\]

断言：

- 电解质 salt balance 相对误差：
  \[
  \epsilon_{M_e} \le 10^{-8}.
  \]
- 固相 balance 绝对误差：
  \[
  |dM_s/dt - M_{s,\mathrm{rhs}}| \le 10^{-12}\ \mathrm{mol\,s^{-1}}.
  \]
- 若加入 Li foil boundary 或 reservoir flux，必须把边界源项写入 balance；不能断言当前 modeled subdomain 的 \(M_e+M_s\) 常数不变。

### 3.6 Khan 2021 电压曲线趋势

新增 `tests/test_p2d_khan_voltage.py`。

算例：

- `khan_1cal_0p2c`, `khan_1cal_0p5c`, `khan_1cal_1c`, `khan_1cal_3c`
- `khan_3cal_0p2c`, `khan_3cal_0p5c`, `khan_3cal_1c`, `khan_3cal_3c`

断言：

- 初始 OCV 在 Khan Figure 4 合理范围：
  \[
  3.8\ \mathrm{V} \le V(0) \le 4.2\ \mathrm{V}
  \]
  具体取决于 `initial_soc`，推荐默认 `initial_soc=0.35` 或与 `data/validation` metadata 对齐。
- 共同容量区间内电压随倍率升高下降。使用容量网格插值后：
  \[
  \mathrm{median}(V_{0.2C}-V_{1C}) \ge 0.02\ \mathrm{V},
  \]
  \[
  \mathrm{median}(V_{1C}-V_{3C}) \ge 0.02\ \mathrm{V}.
  \]
- 与参考曲线定量比较时：
  \[
  \mathrm{RMS}(V_{\mathrm{p2d}}-V_{\mathrm{ref}}) \le 0.01\ \mathrm{V}
  \]
  仅适用于参数、边界条件、OCV、current basis 完全一致的 PyBaMM/DFN 参考。
- 对 Khan PNM/实验图只要求趋势级阈值：
  \[
  \mathrm{RMS} \le 0.10\ \mathrm{V}
  \]
  或通过单调趋势和容量范围检查，因为 PNM 结构解析、current basis 和 CBD percolation 与均匀 P2D 不完全等价。

### 3.7 浓度分布形状

新增 `tests/test_p2d_spatial_profiles.py`。

在目标 SoL 或固定容量点提取 \(c_e(x)\)、\(\phi_e(x)\)、\(c_s^{\mathrm{surf}}(x)\)、\(\bar c_s(x)\)。

断言：

- 放电后：
  \[
  \mathrm{mean}(\bar c_s(t)) > \mathrm{mean}(\bar c_s(0)).
  \]
- 高倍率 cathode 中电解质浓度存在显著梯度：
  \[
  \max(c_e)-\min(c_e) \ge 0.05\,c_{e,0}
  \]
  对 1C；3C 可提高到 \(0.15\,c_{e,0}\)。
- 电解质浓度保持正值：
  \[
  \min(c_e) \ge 1.0\ \mathrm{mol\,m^{-3}}.
  \]
- 固相浓度保持物理范围：
  \[
  0 < c_s < c_{s,\max}.
  \]
- 高倍率下表面-平均浓度差大于低倍率：
  \[
  \max|c_s^{\mathrm{surf}}-\bar c_s|_{3C}
  >
  \max|c_s^{\mathrm{surf}}-\bar c_s|_{0.2C}.
  \]
- 若以 Khan Figure 7 为参考，1CAL 3C 的 SoL 沿 \(x\) 方向应表现为 separator 端较高、collector 端较低的趋势。用线性拟合斜率检查：
  \[
  \frac{d\,\mathrm{SoL}}{dx} < 0
  \]
  只作为趋势测试，不作为精确数值基准。

### 3.8 网格与时间步收敛

新增 `tests/test_p2d_convergence.py`，默认标记为 `pytest.mark.slow`。

使用同一算例跑三组网格：

- coarse: \(N_x=20, N_r=8, \Delta t=20\,s\)
- medium: \(N_x=40, N_r=16, \Delta t=10\,s\)
- fine: \(N_x=80, N_r=32, \Delta t=5\,s\)

断言：

- 电压曲线收敛：
  \[
  \mathrm{RMS}(V_{\mathrm{medium}}-V_{\mathrm{fine}})
  <
  0.7\,\mathrm{RMS}(V_{\mathrm{coarse}}-V_{\mathrm{medium}}).
  \]
- 最终容量差：
  \[
  |Q_{\mathrm{medium}}-Q_{\mathrm{fine}}|/Q_{\mathrm{fine}} \le 0.02.
  \]

当前 Jacobian 仍为 finite-difference fallback，大网格可能较慢。slow 测试可在本地运行，CI 只跑 coarse/medium smoke。

## 4. 实现细节

### 4.1 Benchmark case 数据结构

新增 `src/pnmcathode/p2d/benchmark.py` 或 `examples/p2d/benchmark.py`：

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

Array = np.ndarray


@dataclass(frozen=True)
class P2DBenchmarkCase:
    name: str
    parameter_source: Literal["khan2021", "pybamm"]
    electrode_kind: Literal["1CAL", "3CAL", "pybamm-default"]
    c_rate: float
    initial_soc: float
    temperature: float
    cutoff_voltage: float
    n_positive_cells: int
    n_separator_cells: int
    n_particle_shells: int
    dt_initial: float
    t_final: float
    expected_voltage_path: Path | None = None
```

### 4.2 从现有 P2D solver 跑出结果

使用当前 API：

```python
from __future__ import annotations

from pnmcathode.config import Kinetics, SolverSettings
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)
from pnmcathode.p2d.domain import P2DRegion
from pnmcathode.p2d.materials import from_config
from pnmcathode.p2d.results import P2DResult
from pnmcathode.p2d.solver import P2DProtocol, P2DSolver


def make_khan_solver(case: P2DBenchmarkCase) -> P2DSolver:
    active = nmc532_khan2021()
    electrolyte = electrolyte_khan2021()
    separator = separator_khan2021(enabled=True)
    settings = SolverSettings(
        temperature=case.temperature,
        newton_tol=1e-8,
        dt_max=max(case.dt_initial, 1.0),
    )
    kinetics = Kinetics(k0=1e-10, alpha_a=0.5, alpha_c=0.5)

    if case.electrode_kind == "1CAL":
        thickness = 129e-6
        epsilon_e = 0.368
        epsilon_s = 0.4928
    elif case.electrode_kind == "3CAL":
        thickness = 88e-6
        epsilon_e = 0.366
        epsilon_s = 0.4974
    else:
        raise ValueError(f"unsupported Khan electrode: {case.electrode_kind}")

    positive = P2DRegion(
        name="positive",
        x_left=separator.thickness,
        x_right=separator.thickness + thickness,
        n_cells=case.n_positive_cells,
        epsilon_e=epsilon_e,
        epsilon_s=epsilon_s,
        particle_radius=5e-6,
        sigma_s=active.sigma,
        bruggeman_e=1.5,
        bruggeman_s=1.5,
    )
    params = from_config(
        active=active,
        electrolyte=electrolyte,
        kinetics=kinetics,
        separator=separator,
        settings=settings,
        positive=positive,
        particle_shells=case.n_particle_shells,
    )
    return P2DSolver(params=params, settings=settings)


def run_p2d_benchmark(case: P2DBenchmarkCase) -> P2DResult:
    solver = make_khan_solver(case)
    state0 = solver.initial_state(soc0=case.initial_soc)
    current_density = c_rate_to_current_density(
        c_rate=case.c_rate,
        epsilon_s=solver.params.positive.epsilon_s,
        thickness=solver.params.positive.length,
        cs_max=solver.params.material.active.cs_max,
    )
    protocol = P2DProtocol(
        current_density=current_density,
        t_final=case.t_final,
        dt_initial=case.dt_initial,
        cutoff_voltage=case.cutoff_voltage,
        save_every=1,
    )
    return solver.run(state0, protocol)
```

电流密度换算建议先使用 active material theoretical areal capacity：

\[
I_{1C} =
\frac{F\,\epsilon_s L_p c_{s,\max}\,\Delta x_{\mathrm{SOC}}}{3600}.
\]

若以 Khan `data/validation` 为参考，则必须使用其中 metadata 的 `I_1C_A_m2` 或 `I_app`，避免 current basis 不一致。

```python
from __future__ import annotations

from pnmcathode.physics.reaction import F


def c_rate_to_current_density(
    c_rate: float,
    epsilon_s: float,
    thickness: float,
    cs_max: float,
    soc_window: float = 1.0,
) -> float:
    i_1c = F * epsilon_s * thickness * cs_max * soc_window / 3600.0
    return -c_rate * i_1c
```

### 4.3 提取电压曲线

当前 `P2DResult` 已提供：

```python
from __future__ import annotations

import numpy as np

from pnmcathode.p2d.results import P2DResult

Array = np.ndarray


def extract_voltage_curve(result: P2DResult) -> tuple[Array, Array, Array]:
    time_s = result.time
    voltage_v = result.voltage
    capacity_ah_m2 = result.capacity_Ah_m2
    return time_s, capacity_ah_m2, voltage_v
```

比较函数：

```python
from __future__ import annotations

import numpy as np

Array = np.ndarray


def compare_voltage_curve(
    candidate_capacity: Array,
    candidate_voltage: Array,
    reference_capacity: Array,
    reference_voltage: Array,
) -> dict[str, float]:
    q_min = max(float(candidate_capacity.min()), float(reference_capacity.min()))
    q_max = min(float(candidate_capacity.max()), float(reference_capacity.max()))
    q_grid = np.linspace(q_min, q_max, 200)
    v_candidate = np.interp(q_grid, candidate_capacity, candidate_voltage)
    v_reference = np.interp(q_grid, reference_capacity, reference_voltage)
    err = v_candidate - v_reference
    return {
        "rms_v": float(np.sqrt(np.mean(err**2))),
        "max_abs_v": float(np.max(np.abs(err))),
        "bias_v": float(np.mean(err)),
        "q_min_Ah_m2": q_min,
        "q_max_Ah_m2": q_max,
    }
```

### 4.4 提取浓度和电位分布

从 `P2DSnapshot.state` 直接读取：

```python
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pnmcathode.p2d.domain import MacroMesh, ParticleMesh
from pnmcathode.p2d.results import P2DSnapshot

Array = np.ndarray


@dataclass(frozen=True)
class P2DProfile:
    x_macro: Array
    x_positive: Array
    c_e: Array
    phi_e: Array
    phi_s: Array
    c_s_surface: Array
    c_s_average: Array
    sol_surface: Array
    sol_average: Array


def extract_profile(
    snapshot: P2DSnapshot,
    macro: MacroMesh,
    particle: ParticleMesh,
    cs_max: float,
) -> P2DProfile:
    state = snapshot.state
    pos = macro.positive_cells
    c_s_surface = state.surface_concentration()
    c_s_average = state.average_solid_concentration(particle)
    return P2DProfile(
        x_macro=macro.x_centers.copy(),
        x_positive=macro.x_centers[pos].copy(),
        c_e=state.c_e.copy(),
        phi_e=state.phi_e.copy(),
        phi_s=state.phi_s.copy(),
        c_s_surface=c_s_surface,
        c_s_average=c_s_average,
        sol_surface=c_s_surface / cs_max,
        sol_average=c_s_average / cs_max,
    )
```

### 4.5 保存可复现参考数据

`P2DResult.to_npz()` 当前只保存 time/voltage/capacity。Phase 1 需要更完整的 benchmark writer：

```python
from __future__ import annotations

from pathlib import Path

import numpy as np

from pnmcathode.p2d.results import P2DResult
from pnmcathode.p2d.solver import P2DSolver


def save_benchmark_npz(path: Path, result: P2DResult, solver: P2DSolver) -> None:
    final = result.final_state()
    np.savez(
        path,
        time=result.time,
        voltage=result.voltage,
        capacity_Ah_m2=result.capacity_Ah_m2,
        x=solver.macro.x_centers,
        x_positive=solver.macro.x_centers[solver.macro.positive_cells],
        c_e_final=final.c_e,
        phi_e_final=final.phi_e,
        phi_s_final=final.phi_s,
        c_s_surface_final=final.surface_concentration(),
        c_s_average_final=final.average_solid_concentration(solver.particle),
        current_density=result.snapshots[0].current_density,
    )
```

## 5. 边界条件验证

Phase 1 必须把边界条件作为单独测试，而不是只通过电压曲线间接判断。

### 5.1 零电流平衡态

初始条件：

\[
c_e(x,0)=c_{e,0},\quad c_s(x,r,0)=x_0 c_{s,\max},\quad
\phi_e(x,0)=0,\quad \phi_s(x,0)=U(x_0).
\]

断言：

- `equilibrium_error()["eta_max"] <= 1e-8`
- `equilibrium_error()["i_f_max"] <= 1e-9`
- `ce_variation <= 1e-12`
- `cs_variation <= 1e-12`

### 5.2 OCV 一致性

对 SOC 网格 `np.linspace(0.3, 0.9, 13)`：

\[
|V_{\mathrm{solver}} - U(\mathrm{SOC})| \le 10^{-8}\ \mathrm{V}.
\]

如果使用 Khan OCV clip 区间外的 SOC，OCV derivative 不作为硬断言。

### 5.3 Li metal/separator 左边界

当前 Phase 0 文档要求 half-cell 左端为 Li/Li+ reference。Phase 1 应明确选择一种边界：

1. **reference + galvanostatic ionic current**：
   \[
   \phi_e(0)=0,\qquad i_e(0)=I_{\mathrm{app}}.
   \]
2. **Li foil BV**：
   \[
   i_{\mathrm{Li}} = i^0_{\mathrm{foil}}
   \left[
   \exp\left(\frac{\alpha_aF\eta_{\mathrm{Li}}}{RT}\right)
   -
   \exp\left(-\frac{\alpha_cF\eta_{\mathrm{Li}}}{RT}\right)
   \right].
   \]

Phase 1 推荐先实现方案 1 以便与 P2D/Newman 对齐；Li foil BV 可作为 Khan PNM 复现增强项。

### 5.4 Separator/positive interface

断言：

- \(c_e\) 和 \(\phi_e\) 在 interface 两侧无跳变异常；有限体积中检查左右 cell 的 face flux 连续。
- separator 侧没有 solid phase：
  \[
  i_s(x=L_s)=0.
  \]
- ionic current 在 interface 连续：
  \[
  |i_e^- - i_e^+|/\max(|I_{\mathrm{app}}|,1) \le 10^{-8}.
  \]

### 5.5 Current collector 右边界

断言：

\[
i_s(L_s+L_p)=I_{\mathrm{app}},\qquad i_e(L_s+L_p)=0.
\]

对放电 \(I_{\mathrm{app}}<0\)，`solid_current_faces[-1]` 应为负值。

## 6. 参数对齐

### 6.1 Khan 2021 参数对齐规则

新增 helper：

```python
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class P2DParameterSnapshot:
    source: str
    temperature: float
    separator_thickness: float
    separator_porosity: float
    positive_thickness: float
    positive_porosity: float
    positive_active_fraction: float
    particle_radius: float
    c_e0: float
    c_s_max: float
    k0: float
    t_plus: float
    sigma_s: float
    bruggeman_e: float
    bruggeman_s: float
    ocv_reference: str
```

每个 benchmark 必须保存该 snapshot。测试断言 snapshot 与 case 名称一致。

Khan 对齐时必须：

- 使用 `nmc532_khan2021()`、`electrolyte_khan2021()`、`separator_khan2021()`。
- 显式使用 `Kinetics(k0=1e-10)`，不要沿用 `config.Kinetics` 默认 `5e-10`。
- 显式使用 `SolverSettings(temperature=303.0)`。
- `OCVModel.reference == "li_metal_empirical"`。
- \(D_e(c_e,T)\)、\(\kappa(c_e,T)\)、\(D_s(c_s,T)\) 均来自项目的 Khan correlations。
- 比较 1CAL/3CAL 时使用 Table I 的厚度、孔隙率和活性相体积分数。

### 6.2 PyBaMM 参数对齐规则

PyBaMM 对齐时必须：

- 从 PyBaMM 导出完整参数，不手工混入 Khan 参数。
- 记录 PyBaMM 版本：
  ```python
  metadata["pybamm_version"] = pybamm.__version__
  ```
- 记录 parameter set 名称。
- 把 OCP、diffusivity、conductivity、porosity、thickness、particle radius、initial concentrations、temperature 全部导出到 JSON/NPZ metadata。
- 若当前 half-cell solver 无法表达 PyBaMM full-cell negative electrode，则只比较 positive-side profiles 或趋势，不用 full-cell terminal voltage 作硬阈值。

### 6.3 参数源隔离断言

新增测试：

```python
from __future__ import annotations


def assert_no_mixed_parameter_sources(snapshot: P2DParameterSnapshot) -> None:
    if snapshot.source == "khan2021":
        assert snapshot.temperature == 303.0
        assert snapshot.c_s_max == 48900.0
        assert snapshot.k0 == 1e-10
        assert snapshot.ocv_reference == "li_metal_empirical"
    elif snapshot.source == "pybamm":
        assert snapshot.ocv_reference.startswith("pybamm:")
    else:
        raise AssertionError(f"unknown parameter source: {snapshot.source}")
```

## 7. Phase 1 完成判据

Phase 1 完成时应满足：

- `pytest tests/test_p2d_*.py` 通过。
- 新增 benchmark tests 中，零电流、OCV、质量守恒、电流守恒全部通过硬阈值。
- Khan 2021 1CAL/3CAL 的倍率趋势和空间分布趋势通过。
- 至少生成一个外部 DFN 参考数据集：
  - PyBaMM 可用时：`data/validation/p2d_pybamm_*.npz`
  - PyBaMM 不可用时：记录 skip，并保留 Khan-based reference tests。
- 对完全一致参数的 DFN 参考：
  \[
  \mathrm{RMS}(V_{\mathrm{p2d}}-V_{\mathrm{ref}}) \le 5\text{--}10\,\mathrm{mV}.
  \]
- 对 Khan PNM/实验参考：不要求 5--10 mV 硬误差，只要求趋势、范围、守恒和边界物理一致。

## 8. 建议新增文件

| 文件 | 作用 |
|---|---|
| `src/pnmcathode/p2d/benchmark.py` | benchmark case、参数快照、curve/profile extraction |
| `examples/p2d/compare_pybamm.py` | 可选 PyBaMM 参考生成 |
| `examples/p2d/run_khan_p2d.py` | Khan 1CAL/3CAL P2D benchmark |
| `tests/test_p2d_parameter_alignment.py` | Khan/PyBaMM 参数源一致性 |
| `tests/test_p2d_boundary_current.py` | half-cell 边界电流和总电流守恒 |
| `tests/test_p2d_half_cell.py` | 零电流、OCV、符号约定 |
| `tests/test_p2d_khan_voltage.py` | Khan voltage trend benchmark |
| `tests/test_p2d_spatial_profiles.py` | \(c_e,\phi_e,c_s\) 分布形状 |
| `tests/test_p2d_convergence.py` | 网格和时间步收敛，slow |

