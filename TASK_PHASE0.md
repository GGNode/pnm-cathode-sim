# P2D+PNM Phase 0 可执行实现策略

本文档把 `PLAN.md` 中 Phase 0（纯 P2D 模型独立实现）细化为可执行任务。Phase 0 的目标是先实现不依赖 PNM 的 `separator | porous cathode` half-cell 1D DFN/P2D solver，并为后续 full-cell 与 PNM coupling 保留接口。

## 0. Phase 0 范围与完成定义

Phase 0 只实现纯 P2D，不引入 PNM effective property provider。默认模型为 Li metal reference | separator | porous NMC cathode：

- 宏观 \(x\) 方向：finite volume，覆盖 separator 与 positive electrode。
- 颗粒 \(r\) 方向：每个 positive electrode macro cell 内一个 spherical particle finite volume grid。
- 时间推进：Backward Euler（BDF1）。
- 非线性求解：monolithic Newton-Raphson，状态包含 \(c_e,\phi_e,\phi_s,c_s(r)\)。
- 传输闭合：Bruggeman effective properties。
- 验证目标：守恒、平衡态、解析极限、低倍率 OCV 极限、与 PyBaMM/文献趋势对比。

完成判据：

- `pytest tests/test_p2d_*.py` 通过。
- \(I_\mathrm{app}=0\) 时浓度不变、\(i_F\approx 0\)、\(\phi_s-\phi_e=U\)。
- 每个时间步 total current conservation 相对误差 \(\le 10^{-8}\) 或由线性求解容差解释。
- 电解质扩散与固相球扩散解析极限误差随网格加密下降。
- 低倍率 discharge voltage 接近 OCV，且 cathode discharge 内部约定为 \(I_\mathrm{app}<0\)、\(i_F<0\)、\(\bar c_s\) 增加。

## 1. 模块设计

### 1.1 新建文件总览

| 文件 | 职责 | 主要依赖 |
|---|---|---|
| `src/pnmcathode/p2d/__init__.py` | P2D public API 与 lazy-safe exports | `domain`, `state`, `materials`, `solver`, `results` |
| `src/pnmcathode/p2d/domain.py` | 区域定义、宏观网格 `MacroMesh`、颗粒网格 `ParticleMesh`、face 几何与区域索引 | `dataclasses`, `numpy` |
| `src/pnmcathode/p2d/state.py` | 状态对象、状态向量布局、pack/unpack、初始状态构造 | `numpy`, `domain` |
| `src/pnmcathode/p2d/materials.py` | P2D 参数容器、Bruggeman closure、从现有 `config`/presets 构造材料参数 | `pnmcathode.config`, `pnmcathode.physics` |
| `src/pnmcathode/p2d/kinetics.py` | vectorized Butler-Volmer、exchange current、导数、OCV 导数适配 | `pnmcathode.physics.reaction`, `pnmcathode.physics.ocv` |
| `src/pnmcathode/p2d/residual.py` | monolithic residual 组装、通量计算、边界条件、diagnostic fluxes | `scipy.sparse`, `domain`, `state`, `materials`, `kinetics` |
| `src/pnmcathode/p2d/jacobian.py` | analytic sparse Jacobian 组装；初期允许 finite-difference fallback 仅用于测试 | `scipy.sparse`, `residual`, `state` |
| `src/pnmcathode/p2d/solver.py` | Newton solve、time step、run protocol、自适应时间步基础逻辑 | `scipy.sparse.linalg`, `residual`, `jacobian`, `results` |
| `src/pnmcathode/p2d/results.py` | `P2DResult`、snapshot、voltage/capacity/diagnostics 输出 | `dataclasses`, `numpy`, `state` |
| `src/pnmcathode/p2d/diagnostics.py` | 守恒检查、current balance、lithium inventory、OCV equilibrium residual | `numpy`, `residual`, `results` |
| `tests/test_p2d_domain.py` | 网格几何、区域索引、face averaging 单元测试 | `pytest`, `numpy` |
| `tests/test_p2d_state.py` | state layout、pack/unpack、初始条件测试 | `pytest`, `numpy` |
| `tests/test_p2d_particle.py` | 球扩散 finite volume 守恒与解析极限测试 | `pytest`, `numpy` |
| `tests/test_p2d_residual.py` | residual 维度、符号、平衡态、边界条件测试 | `pytest`, `numpy` |
| `tests/test_p2d_conservation.py` | 质量/电荷守恒、total current conservation 测试 | `pytest`, `numpy` |
| `tests/test_p2d_solver.py` | Newton 收敛、零电流、低倍率、极端倍率 smoke tests | `pytest`, `numpy` |
| `examples/p2d/run_half_cell.py` | 最小 half-cell discharge 示例 | `pnmcathode.p2d` |
| `examples/p2d/compare_pybamm.py` | Phase 1 可继续扩展的 PyBaMM 对比脚本骨架 | optional `pybamm` |

现有文件需要修改：

- `src/pnmcathode/__init__.py`：把 `P2DSolver`, `P2DResult`, `P2DRegion`, `P2DParameters` 加入 `_LAZY_EXPORTS`。
- `pyproject.toml`：如当前未声明 `scipy` 测试依赖，需要确保 P2D solver 使用的 `scipy.sparse` 已在依赖中。

### 1.2 `domain.py`

职责：

- 定义 Phase 0 的区域、宏观 finite volume mesh 与颗粒 finite volume mesh。
- 明确 separator cells、positive electrode cells、electrode/separator interface face。
- 提供 harmonic/series averaging 所需的 face 几何。

关键签名：

```python
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal
import numpy as np

Array = np.ndarray
RegionName = Literal["negative", "separator", "positive"]

@dataclass(frozen=True)
class P2DRegion:
    name: RegionName
    x_left: float
    x_right: float
    n_cells: int
    epsilon_e: float
    epsilon_s: float = 0.0
    particle_radius: float | None = None
    specific_area: float | None = None
    bruggeman_e: float = 1.5
    bruggeman_s: float = 1.5
    sigma_s: float | None = None

    @property
    def length(self) -> float: ...

    def area_density(self) -> float: ...

@dataclass(frozen=True)
class MacroMesh:
    regions: tuple[P2DRegion, ...]
    x_faces: Array
    x_centers: Array
    dx: Array
    volumes: Array
    face_areas: Array
    region_id: Array
    separator_cells: Array
    positive_cells: Array
    electrode_cells: Array

    @classmethod
    def from_regions(
        cls,
        regions: tuple[P2DRegion, ...],
        area: float = 1.0,
    ) -> "MacroMesh": ...

    def region_for_cell(self, cell: int) -> P2DRegion: ...
    def left_face(self, cell: int) -> int: ...
    def right_face(self, cell: int) -> int: ...

@dataclass(frozen=True)
class ParticleMesh:
    radius: float
    n_shells: int
    r_faces: Array
    r_centers: Array
    shell_volumes: Array
    face_areas: Array

    @classmethod
    def spherical(cls, radius: float, n_shells: int) -> "ParticleMesh": ...
```

依赖关系：

- `state.py` 使用 `MacroMesh` 与 `ParticleMesh` 决定 state layout。
- `residual.py` 使用 cell volumes、face areas、region masks 组装 FVM residual。
- `materials.py` 使用 `P2DRegion` 的 \(\epsilon_e,\epsilon_s,a_s,b\) 计算 effective properties。

### 1.3 `state.py`

职责：

- 组织 monolithic unknown vector。
- 将物理状态与 Newton 向量互转。
- 保证 separator 没有 \(\phi_s\) 与 \(c_s\) unknown。

推荐状态向量顺序：

\[
y = [c_e(0:N_x),\ \phi_e(0:N_x),\ \phi_s(0:N_p),\ c_s(0:N_p,0:N_r)]
\]

其中 \(N_x\) 为 separator + positive 总 cell 数，\(N_p\) 为 positive electrode cell 数，\(N_r\) 为颗粒 shell 数。`phi_s` 和 `c_s` 只在 positive cells 上定义。

关键签名：

```python
from dataclasses import dataclass
import numpy as np

from pnmcathode.p2d.domain import MacroMesh, ParticleMesh

Array = np.ndarray

@dataclass(frozen=True)
class StateLayout:
    n_x: int
    n_pos: int
    n_r: int
    c_e: slice
    phi_e: slice
    phi_s: slice
    c_s: slice
    positive_to_macro: Array
    macro_to_positive: Array

    @property
    def size(self) -> int: ...

    @classmethod
    def from_meshes(cls, macro: MacroMesh, particle: ParticleMesh) -> "StateLayout": ...

@dataclass
class P2DState:
    c_e: Array
    phi_e: Array
    phi_s: Array
    c_s: Array
    time: float = 0.0

    def copy(self) -> "P2DState": ...
    def surface_concentration(self) -> Array: ...
    def average_solid_concentration(self, particle: ParticleMesh) -> Array: ...

def pack_state(state: P2DState, layout: StateLayout) -> Array: ...

def unpack_state(y: Array, layout: StateLayout) -> P2DState: ...

def make_initial_state(
    macro: MacroMesh,
    particle: ParticleMesh,
    c_e0: float,
    soc0: float,
    c_s_max: float,
    phi_e0: float = 0.0,
) -> P2DState: ...
```

依赖关系：

- `solver.py` 对 Newton vector 只使用 `pack_state/unpack_state`。
- `diagnostics.py` 通过 `average_solid_concentration()` 做 lithium inventory。

### 1.4 `materials.py`

职责：

- 将现有 `pnmcathode.config.ActiveMaterial`、`Electrolyte`、`Kinetics`、`Separator` 适配到 P2D 参数。
- 实现 Bruggeman closure。
- 明确 OCV reference，避免重复加入 electrolyte Nernst correction。

关键签名：

```python
from dataclasses import dataclass
from typing import Callable, Literal
import numpy as np

from pnmcathode.config import ActiveMaterial, Electrolyte, Kinetics, Separator, SolverSettings
from pnmcathode.p2d.domain import P2DRegion

Array = np.ndarray

@dataclass(frozen=True)
class OCVModel:
    reference: Literal["li_metal_empirical", "thermodynamic_with_electrolyte"]
    value: Callable[[Array], Array]
    derivative: Callable[[Array], Array] | None = None

@dataclass(frozen=True)
class P2DMaterial:
    active: ActiveMaterial
    electrolyte: Electrolyte
    kinetics: Kinetics
    separator: Separator
    temperature: float
    ocv_model: OCVModel

@dataclass(frozen=True)
class P2DParameters:
    material: P2DMaterial
    positive: P2DRegion
    separator_region: P2DRegion
    particle_shells: int = 20
    area: float = 1.0
    concentration_floor: float = 1e-6
    potential_gauge: Literal["li_phi_e_left"] = "li_phi_e_left"

def from_config(
    active: ActiveMaterial,
    electrolyte: Electrolyte,
    kinetics: Kinetics,
    separator: Separator,
    settings: SolverSettings,
    positive: P2DRegion,
    particle_shells: int = 20,
    area: float = 1.0,
) -> P2DParameters: ...

def bruggeman_diffusivity_e(
    c_e: Array,
    region: P2DRegion,
    material: P2DMaterial,
) -> Array: ...

def bruggeman_conductivity_e(
    c_e: Array,
    region: P2DRegion,
    material: P2DMaterial,
) -> Array: ...

def bruggeman_conductivity_s(
    region: P2DRegion,
    material: P2DMaterial,
) -> float: ...
```

依赖关系：

- 直接复用 `pnmcathode.physics.electrolyte.electrolyte_diffusion_coefficient`、`electrolyte_ionic_conductivity`。
- 直接复用 `pnmcathode.physics.solid.nmc532_diffusion_coefficient`。
- 直接复用 `pnmcathode.physics.ocv.nmc532_ocv` 与 `ocv_derivative`。
- 直接复用 `pnmcathode.config` 的参数校验，不复制参数 dataclass。

### 1.5 `kinetics.py`

职责：

- 为 P2D residual 提供 vectorized Butler-Volmer 与 analytic derivatives。
- 复用现有标量 `butler_volmer()` 与 `exchange_current_density()` 的符号约定，但避免 residual 内逐点 Python loop。

关键签名：

```python
from dataclasses import dataclass
import numpy as np

from pnmcathode.config import Kinetics

Array = np.ndarray

@dataclass(frozen=True)
class ReactionRates:
    eta: Array
    i0: Array
    i_f: Array
    di_deta: Array
    di_dc_e: Array
    di_dc_surf: Array

def exchange_current_density_vec(
    c_e: Array,
    c_s_surf: Array,
    c_s_max: float,
    kinetics: Kinetics,
    floor: float = 1e-12,
) -> Array: ...

def butler_volmer_with_derivatives(
    c_e: Array,
    c_s_surf: Array,
    phi_e: Array,
    phi_s: Array,
    c_s_max: float,
    ocv: Array,
    docv_dc_s: Array,
    kinetics: Kinetics,
    temperature: float,
) -> ReactionRates: ...
```

依赖关系：

- `residual.py` 使用 `ReactionRates.i_f` 与导数。
- `jacobian.py` 使用 `di_deta`、`di_dc_e`、`di_dc_surf`。

### 1.6 `residual.py`

职责：

- 组装完整 nonlinear residual。
- 所有 conserved equations 使用 finite volume flux difference。
- 在同一 residual 内耦合 electrolyte concentration、electrolyte potential、solid potential、solid particle diffusion。
- 实现 half-cell 边界条件与 potential gauge。

关键签名：

```python
from dataclasses import dataclass
import numpy as np
from scipy import sparse

from pnmcathode.p2d.domain import MacroMesh, ParticleMesh
from pnmcathode.p2d.materials import P2DParameters
from pnmcathode.p2d.state import P2DState, StateLayout

Array = np.ndarray

@dataclass(frozen=True)
class P2DResidualContext:
    macro: MacroMesh
    particle: ParticleMesh
    layout: StateLayout
    params: P2DParameters
    dt: float
    current_density: float

@dataclass(frozen=True)
class P2DFluxes:
    salt_flux_faces: Array
    ionic_current_faces: Array
    solid_current_faces: Array
    faradaic_current: Array
    volumetric_reaction: Array

def compute_reaction(
    state: P2DState,
    context: P2DResidualContext,
) -> Array: ...

def compute_fluxes(
    state: P2DState,
    context: P2DResidualContext,
) -> P2DFluxes: ...

def assemble_p2d_residual(
    y_new: Array,
    state_old: P2DState,
    context: P2DResidualContext,
) -> Array: ...

def residual_norm(
    residual: Array,
    context: P2DResidualContext,
) -> float: ...
```

核心 residual 单位：

- \(c_e\) residual：`mol m^-3`，用 cell volume 或 porosity scaling 后保持量纲一致。
- \(\phi_e,\phi_s\) residual：`A m^-3` 或按 cell volume 后 `A`，同一文件内必须统一。
- \(c_s\) residual：`mol m^-3`。

建议：实现 residual 时先采用 per-volume form，所有方程都除以 macro cell volume 或 particle shell volume，再在 `residual_norm()` 中按变量尺度归一化。

### 1.7 `jacobian.py`

职责：

- 组装 sparse analytic Jacobian。
- 初期可提供 finite-difference Jacobian 用于单元测试校验，但目标 solver 默认 analytic Jacobian。

关键签名：

```python
import numpy as np
from scipy import sparse

from pnmcathode.p2d.residual import P2DResidualContext
from pnmcathode.p2d.state import P2DState

Array = np.ndarray

def assemble_p2d_jacobian(
    y_new: Array,
    state_old: P2DState,
    context: P2DResidualContext,
) -> sparse.csr_matrix: ...

def finite_difference_jacobian(
    y_new: Array,
    state_old: P2DState,
    context: P2DResidualContext,
    rel_step: float = 1e-7,
) -> sparse.csr_matrix: ...

def check_jacobian(
    y_new: Array,
    state_old: P2DState,
    context: P2DResidualContext,
    rel_step: float = 1e-7,
) -> dict[str, float]: ...
```

Jacobian 组装策略：

- 用 `scipy.sparse.lil_matrix` 或 COO triplets 收集局部 stencil。
- \(x\) 方向通量只耦合相邻 macro cells。
- \(r\) 方向固相扩散只耦合相邻 shells。
- BV 源项在同一 positive macro cell 内耦合 \(c_e,\phi_e,\phi_s,c_s^\mathrm{surf}\)。
- 若 Phase 0 先忽略物性导数 \(dD/dc\)、\(d\kappa/dc\)，必须在文档和测试中标记为 quasi-Newton approximation；推荐在实现第二步补齐。

### 1.8 `solver.py`

职责：

- Newton-Raphson nonlinear solve。
- Backward Euler time step。
- 恒流 protocol run。
- 处理 step rejection 与基础自适应 dt。

关键签名：

```python
from dataclasses import dataclass
from typing import Callable
import numpy as np

from pnmcathode.config import SolverSettings
from pnmcathode.p2d.materials import P2DParameters
from pnmcathode.p2d.results import P2DResult
from pnmcathode.p2d.state import P2DState

Array = np.ndarray

@dataclass(frozen=True)
class P2DProtocol:
    current_density: float
    t_final: float
    dt_initial: float
    cutoff_voltage: float | None = None
    save_every: int = 1

@dataclass(frozen=True)
class NewtonReport:
    converged: bool
    iterations: int
    residual_norm: float
    step_norm: float

class P2DSolver:
    def __init__(
        self,
        params: P2DParameters,
        settings: SolverSettings | None = None,
        use_analytic_jacobian: bool = True,
    ) -> None: ...

    def initial_state(self, soc0: float, c_e0: float | None = None) -> P2DState: ...

    def step(
        self,
        state: P2DState,
        dt: float,
        current_density: float,
    ) -> tuple[P2DState, NewtonReport]: ...

    def run(
        self,
        state0: P2DState,
        protocol: P2DProtocol,
        progress: Callable[[P2DState], None] | None = None,
    ) -> P2DResult: ...

    def voltage(self, state: P2DState) -> float: ...
```

Newton 默认参数：

- `newton_tol = 1e-8`，使用 scaled residual norm。
- `newton_max_iter = 30` 用于 P2D；现有 `SolverSettings.newton_max_iter=200` 可保留但 P2D 默认应更严格暴露。
- damping/line search：Armijo 或 residual-decrease backtracking，初始 damping \(1.0\)，最小 \(1/16\)。
- 变量 bounds：\(c_e>c_\mathrm{floor}\)，\(0<c_s<c_{s,\max}\)。Newton trial 超界时先 damping，不直接 clip accepted state。

### 1.9 `results.py`

职责：

- 保存 transient snapshots。
- 提供 voltage、capacity、average concentrations、diagnostics。

关键签名：

```python
from dataclasses import dataclass, field
import numpy as np

from pnmcathode.p2d.state import P2DState

Array = np.ndarray

@dataclass
class P2DSnapshot:
    time: float
    state: P2DState
    voltage: float
    current_density: float
    diagnostics: dict[str, float] = field(default_factory=dict)

@dataclass
class P2DResult:
    snapshots: list[P2DSnapshot]
    metadata: dict[str, float | str]

    @property
    def time(self) -> Array: ...

    @property
    def voltage(self) -> Array: ...

    @property
    def capacity_Ah_m2(self) -> Array: ...

    def final_state(self) -> P2DState: ...
    def to_npz(self, path: str) -> None: ...
```

### 1.10 `diagnostics.py`

职责：

- Phase 0 验证与 solver runtime diagnostics。

关键签名：

```python
import numpy as np

from pnmcathode.p2d.residual import P2DFluxes, P2DResidualContext
from pnmcathode.p2d.state import P2DState

def total_lithium_inventory(
    state: P2DState,
    context: P2DResidualContext,
) -> float: ...

def current_conservation_error(
    fluxes: P2DFluxes,
    context: P2DResidualContext,
) -> float: ...

def mass_conservation_error(
    state_new: P2DState,
    state_old: P2DState,
    fluxes: P2DFluxes,
    context: P2DResidualContext,
) -> float: ...

def equilibrium_error(
    state: P2DState,
    context: P2DResidualContext,
) -> dict[str, float]: ...
```

## 2. 公式实现清单

本节只列 Phase 0 需要数值实现的 `PLAN.md` §2 方程。PNM §3 方程不在 Phase 0 实现范围内。

| 方程编号/名称 | 方程 | 离散化方法 | 实现位置 | 验证方法 |
|---|---|---|---|---|
| 2.1 specific area，球形颗粒 | \(a_{s,k}=3\epsilon_{s,k}/R_{p,k}\) | 代数计算 | `domain.P2DRegion.area_density()` | 单元测试：给定 \(\epsilon_s,R_p\) 与手算值一致 |
| 2.1 specific area，外部给定 | \(a_{s,k}=A_{\mathrm{e/NMC},k}/V_{\mathrm{RVE},k}\) | Phase 0 不接 PNM；保留 `specific_area` override | `domain.P2DRegion.area_density()` | override 优先级测试 |
| 2.2 solid Fick flux | \(N_{s,k}=-D_{s,k}\nabla c_{s,k}\) | 颗粒 \(r\) 方向 finite volume，两点通量，face \(D_s\) harmonic/center average | `residual._solid_particle_fluxes()` | 固相总锂变化等于表面通量；与 `physics.solid.solid_diffusion_rhs` 小算例对比 |
| 2.2 spherical particle diffusion | \(\partial c_s/\partial t=\frac{1}{r^2}\frac{\partial}{\partial r}(r^2D_s\frac{\partial c_s}{\partial r})\) | 球壳 finite volume；Backward Euler | `residual._assemble_solid_diffusion_residual()` | 恒定 \(D_s\)、恒定表面通量解析/高精度参考；网格收敛 |
| 2.2 center symmetry | \(\left.\partial c_s/\partial r\right|_{r=0}=0\) | \(r=0\) face area 为 0，无通量 | `domain.ParticleMesh.spherical()`, `residual._solid_particle_fluxes()` | 中心 face flux 恒为 0 |
| 2.2 surface reaction flux | \(-D_s\left.\partial c_s/\partial r\right|_{R_p}=i_F/F\) | 外边界通量作为最外 shell source；按 PLAN 符号 \(i_F>0\) 脱锂，固相质量项为 \(-i_F/F\) | `residual._assemble_solid_diffusion_residual()` | \(\epsilon_s d\bar c_s/dt=-a_si_F/F\) 守恒检查 |
| 2.2 average concentration | \(\bar c_s=3R_p^{-3}\int_0^{R_p}c_sr^2dr\) | shell volume weighted average | `state.P2DState.average_solid_concentration()` | 均匀 \(c_s\) 返回同值；线性/二次 profile 与数值积分对比 |
| 2.2 average solid source | \(\epsilon_s\partial\bar c_s/\partial t=-a_si_F/F\) | 不作为独立 unknown；由颗粒 FVM 积分验证 | `diagnostics.mass_conservation_error()` | 每个 positive cell 颗粒积分与 BV source 相等 |
| 2.3 electrolyte salt conservation | \(\partial(\epsilon_ec_e)/\partial t=\partial_x(D_e^\mathrm{eff}\partial_xc_e)+(1-t_+^0)a_si_F/F\) | \(x\) 方向 finite volume；Backward Euler；face \(D_e^\mathrm{eff}\) harmonic average | `residual._assemble_electrolyte_concentration_residual()` | 无反应 Fick 解析解；全局 salt inventory；界面 flux 连续 |
| 2.3 Bruggeman diffusivity | \(D_e^\mathrm{eff}=D_e(c_e,T)\epsilon_e^{b_D}\) | cell property；face harmonic average | `materials.bruggeman_diffusivity_e()` | 与手算值；\(\epsilon_e=1\) 退化到 bulk \(D_e\) |
| 2.3 electrolyte current，concentrated solution | \(i_e=-\kappa^\mathrm{eff}\partial_x\phi_e+\frac{2RT\kappa^\mathrm{eff}}{F}(1-t_+^0)\chi\partial_x\ln c_e\) | \(x\) 方向 finite volume face current；\(\ln c_e\) 使用 bounded concentration；Phase 0 默认 \(\chi=1\) | `residual._electrolyte_current_faces()` | \(c_e\) 均匀时退化 Ohmic；浓差项符号测试；PyBaMM benchmark |
| 2.3 thermodynamic factor | \(\chi(c_e)=1+\partial\ln f_\pm/\partial\ln c_e\) | Phase 0 默认理想溶液 \(\chi=1\)，预留 callable | `materials.P2DMaterial` 后续字段 | 单元测试默认返回 1 |
| 2.3 simplified Ohmic current | \(i_e=-\kappa^\mathrm{eff}\partial_x\phi_e\) | concentrated current 的开关特例 | `residual._electrolyte_current_faces()` | 关闭浓差项时与解析 Ohmic slab 对比 |
| 2.3 electrolyte charge conservation | \(\partial_x i_e=a_si_F\) | finite volume current divergence | `residual._assemble_electrolyte_potential_residual()` | \(\partial_x(i_e+i_s)=0\)；separator 中 \(a_s=0\) 时 \(i_e\) 常数 |
| 2.4 solid Ohm law | \(i_s=-\sigma^\mathrm{eff}\partial_x\phi_s\) | positive electrode \(x\) 方向 finite volume | `residual._solid_current_faces()` | 无反应 Ohmic slab 解析解 |
| 2.4 solid charge conservation | \(\partial_x i_s=-a_si_F\) | finite volume current divergence | `residual._assemble_solid_potential_residual()` | 与 electrolyte charge 相加后 total current constant |
| 2.4 total current conservation | \(\partial_x(i_e+i_s)=0\) | diagnostic，不作为独立 residual | `diagnostics.current_conservation_error()` | 任意 converged step 上误差随 Newton tolerance 下降 |
| 2.4 Bruggeman solid conductivity | \(\sigma^\mathrm{eff}=\sigma\epsilon_s^{b_\sigma}\) | cell property；face harmonic average | `materials.bruggeman_conductivity_s()` | \(\epsilon_s=1\) 退化到 bulk \(\sigma\) |
| 2.5 overpotential | \(\eta=\phi_s-\phi_e-U(c_s^\mathrm{surf},c_e,T)\) | positive cell local algebraic coupling | `kinetics.butler_volmer_with_derivatives()` | 平衡态 \(\eta=0\Rightarrow i_F=0\) |
| 2.5 Butler-Volmer | \(i_F=i_0[\exp(\alpha_aF\eta/RT)-\exp(-\alpha_cF\eta/RT)]\) | algebraic nonlinear source，指数参数 clip | `kinetics.butler_volmer_with_derivatives()` | 与 `physics.reaction.butler_volmer()` 对比；小 \(\eta\) 线性化 |
| 2.5 exchange current | \(i_0=Fk_0c_e^{\gamma_e}(c_s^\max-c_s^\mathrm{surf})^{\gamma_v}(c_s^\mathrm{surf})^{\gamma_s}\) | algebraic；bounded concentrations | `kinetics.exchange_current_density_vec()` | 与 `physics.reaction.exchange_current_density()` 对比 |
| 2.5 BV potential derivatives | \(B=\partial i_F/\partial\eta=i_0f(\alpha_aE_a+\alpha_cE_c)\) | analytic Jacobian coefficient | `kinetics.butler_volmer_with_derivatives()`, `jacobian.assemble_p2d_jacobian()` | finite-difference Jacobian check |
| 2.5 BV concentration derivatives | \(\partial i_F/\partial c_s=(\partial i_0/\partial c_s)(E_a-E_c)-B\partial U/\partial c_s\)，\(\partial i_F/\partial c_e=(\partial i_0/\partial c_e)(E_a-E_c)-B\partial U/\partial c_e\) | analytic Jacobian coefficient；empirical OCV 默认 \(\partial U/\partial c_e=0\) | `kinetics.butler_volmer_with_derivatives()`, `jacobian.assemble_p2d_jacobian()` | finite-difference Jacobian check；OCV derivative 与 `physics.ocv.ocv_derivative()` 对比 |
| 2.6.1 no electrolyte flux at metal | \(N_e\cdot n=0\) | boundary face flux = 0；half-cell right current collector 适用 | `residual._salt_flux_faces()` | 边界 flux 单元测试 |
| 2.6.1 no ionic current at current collector | \(i_e\cdot n=0\) | positive current collector face \(i_e=0\) | `residual._electrolyte_current_faces()` | boundary current 单元测试 |
| 2.6.1 applied solid current | \(i_s\cdot n=I_\mathrm{app}\) | positive current collector Neumann current face；项目约定 cathode discharge \(I_\mathrm{app}<0\) | `residual._solid_current_faces()` | total current 等于输入 \(I_\mathrm{app}\) |
| 2.6.1 gauge | 固定 \(\phi_e(0,t)=0\) 或 \(\phi_s(0,t)=0\) | half-cell 使用 left Li reference \(\phi_e(0)=0\) Dirichlet residual 替换 | `residual._apply_half_cell_boundary_conditions()` | Jacobian 非奇异；整体电势平移不再是 null mode |
| 2.6.2 electrode/separator continuity | \(c_e^- = c_e^+\)，\(\phi_e^-=\phi_e^+\)，\(N_e^-=N_e^+\)，\(i_e^-=i_e^+\) | 单一连续 macro mesh；界面 face 使用 harmonic/series effective coefficient | `domain.MacroMesh.from_regions()`, `residual._face_coefficients()` | interface flux continuity；阶跃 \(\epsilon\) 解析 steady flux |
| 2.6.2 no solid current into separator | \(i_s\cdot n=0\) | positive left boundary solid current face = 0 | `residual._solid_current_faces()` | separator interface solid current 为 0 |
| 2.6.3 Li metal reference potential | \(\phi_e(0,t)=0\) | left boundary Dirichlet gauge | `residual._apply_half_cell_boundary_conditions()` | \(I=0\) 平衡态测试 |
| 2.6.3 Li reservoir concentration | \(c_e(0,t)=c_{e,\mathrm{ref}}\) | Phase 0 建议做可选 Dirichlet；默认 no-flux 或 reservoir 由参数选择 | `materials.P2DParameters`, `residual._salt_flux_faces()` | Dirichlet diffusion 解析解；no-flux 全局守恒 |
| 2.6.3 Li foil BV | \(i_\mathrm{Li}=i_{0,\mathrm{Li}}[\exp(\alpha_aF\eta_\mathrm{Li}/RT)-\exp(-\alpha_cF\eta_\mathrm{Li}/RT)]\) | Phase 0 可延后；保留 `Separator.include_li_foil_bv` 接口，不作为首版必需 | `residual._apply_half_cell_boundary_conditions()` 后续 | 开关关闭时不影响 baseline；开启后与标量 BV 对比 |
| 2.6.4 half-cell voltage | \(V_\mathrm{cell}=\phi_{s,p}(L,t)\) | 从 positive current collector solid potential 读取；Li reference 为 0 | `solver.P2DSolver.voltage()` | \(I=0\) 时等于平均 OCV；discharge 初期低于 OCV |

## 3. 数值方案选择

### 3.1 空间离散：选择 finite volume

Phase 0 应采用 finite volume，而不是 finite difference。

理由：

- DFN 的核心约束是质量与电荷守恒。finite volume 直接对 control volume 写 flux balance，能让 \(\epsilon_e c_e\)、\(\epsilon_s c_s\)、\(i_e+i_s\) 的守恒误差成为可测量的 residual。
- separator/positive electrode 界面存在 \(\epsilon\)、\(D^\mathrm{eff}\)、\(\kappa^\mathrm{eff}\) 阶跃。finite volume 可在 face 上使用 harmonic/series averaging，自然保证 flux continuity。
- 颗粒球坐标方程含 \(r^2\) 几何因子。球壳 finite volume 用 shell volume 与 face area 可避免 \(r=0\) 奇异项。
- 后续 PNM 本质也是 graph finite volume；Phase 0 使用 finite volume 能降低 Phase 3 coupling 的概念差异。

实现要求：

- \(x\) 方向所有 transport term 都写成 face flux divergence。
- \(r\) 方向所有 solid diffusion term 都写成 shell face flux divergence。
- 不在 residual 中直接写二阶中心差分形式，除非它可证明等价于 uniform finite volume。

### 3.2 时间积分：首版选择 Backward Euler / BDF1

首版选择 Backward Euler（BDF1），暂不使用 Crank-Nicolson。

理由：

- DFN 是 stiff nonlinear DAE-like system，Backward Euler L-stable，对高倍率、低扩散系数、强 BV 非线性更稳健。
- Crank-Nicolson 二阶但 A-stable 非 L-stable，在 sharp transients 或 cutoff 附近容易振荡，并且对 nonlinear solve 更敏感。
- BDF1 与 Phase 0 守恒验证最直接；后续可扩展 BDF2，但不能在首版增加历史状态复杂度。

残差形式：

\[
R_c(y^{n+1}) = M\frac{c^{n+1}-c^n}{\Delta t} - \mathcal{F}(y^{n+1}) = 0
\]

其中所有通量、反应源项、物性系数默认在 \(t^{n+1}\) 评价。

### 3.3 非线性求解：Newton-Raphson

采用 monolithic Newton-Raphson：

\[
J(y^m)\Delta y^m=-R(y^m),\qquad y^{m+1}=y^m+\lambda\Delta y^m
\]

参数建议：

- scaled residual tolerance：`1e-8`。
- absolute fallback tolerance：`1e-10`。
- max iterations：`30`。
- damping：从 \(\lambda=1\) 开始，若 residual 不下降或浓度越界，则依次尝试 \(1/2,1/4,1/8,1/16\)。
- initial guess：上一步 converged state；第一步用 equilibrium-consistent potential 初始化。

Jacobian 组装：

- 目标实现为 analytic sparse Jacobian。
- \(x\) 方向 diffusion/conduction stencil：三对角或 block 三对角。
- \(r\) 方向 particle diffusion stencil：每个 macro positive cell 内三对角。
- BV coupling：每个 positive macro cell 内局部 dense coupling，连接 \(c_e[j]\)、\(\phi_e[j]\)、\(\phi_s[j]\)、\(c_s[j,N_r-1]\)。
- 物性导数：
  - 首个可运行版本允许忽略 \(dD_e/dc_e\)、\(d\kappa/dc_e\)、\(dD_s/dc_s\)，作为 quasi-Newton。
  - 正式 Phase 0 完成前应补齐或至少用 finite-difference Jacobian test 证明 quasi-Newton 收敛范围。
- `finite_difference_jacobian()` 只用于小网格测试，不作为默认求解路径。

### 3.4 线性求解：首版使用 `scipy.sparse` 直接法

首版使用 `scipy.sparse.linalg.spsolve` 或 `splu`。

理由：

- Phase 0 网格规模预计 \(N_x\sim 20-100\)、\(N_r\sim 10-30\)，unknown 数量通常小于几千，直接法足够。
- 直接法减少 Krylov preconditioner 带来的调试变量，适合先锁定符号、单位、边界条件。
- 后续大规模或 coupled PNM 才需要 GMRES/BiCGSTAB + block preconditioner。

实现建议：

- Jacobian 使用 CSR/CSC；传入 `spsolve` 前转 CSC。
- 捕获 singular matrix 或 factorization failure，输出 gauge、边界条件、zero conductivity/percolation 提示。
- Newton report 记录 linear solve residual \(\|J\Delta y+R\|\)。

## 4. 数据结构设计

### 4.1 `MacroMesh`

`MacroMesh` 是 \(x\) 方向有限体积控制体集合，几何面积默认 \(A=1\ \mathrm{m^2}\)，因此 current density 与 total current 数值一致，后续可设置真实 area。

必须存储：

- `x_faces: Array`，shape `(n_x + 1,)`。
- `x_centers: Array`，shape `(n_x,)`。
- `dx: Array`，shape `(n_x,)`。
- `volumes: Array`，shape `(n_x,)`，\(V_j=A\Delta x_j\)。
- `face_areas: Array`，shape `(n_x + 1,)`。
- `region_id: Array`，每个 cell 所属区域。
- `separator_cells`, `positive_cells`, `electrode_cells` masks/indices。

区域顺序：

1. `separator`: \(x\in[0,L_s]\)。
2. `positive`: \(x\in[L_s,L_s+L_p]\)。

预留 full-cell 时可在左侧加入 `negative`，但 Phase 0 不要求实现 negative electrode residual。

### 4.2 `ParticleMesh`

`ParticleMesh` 是每个 positive macro cell 共享的 spherical particle radial mesh。

必须存储：

- `r_faces: Array`，shape `(n_r + 1,)`，从 0 到 \(R_p\)。
- `r_centers: Array`，shell volume centroid 或 midpoint。
- `shell_volumes: Array`，\(\frac{4\pi}{3}(r_{l+1}^3-r_l^3)\)。
- `face_areas: Array`，\(4\pi r_l^2\)，其中 \(r_0=0\) 的 face area 为 0。

如果后续不同 macro cells 有不同 \(R_p\)，可把 `ParticleMesh` 扩展为 region/cell map；Phase 0 使用 uniform particle mesh。

### 4.3 状态向量组织

物理状态对象：

```python
P2DState(
    c_e=np.ndarray,      # shape (n_x,), mol m^-3
    phi_e=np.ndarray,    # shape (n_x,), V
    phi_s=np.ndarray,    # shape (n_pos,), V
    c_s=np.ndarray,      # shape (n_pos, n_r), mol m^-3
    time=float,
)
```

Newton 向量 layout：

```text
[0 : n_x)                         c_e
[n_x : 2 n_x)                     phi_e
[2 n_x : 2 n_x + n_pos)           phi_s
[2 n_x + n_pos : end)             c_s flattened in C order
```

关键映射：

- `positive_to_macro[p] -> j`：positive local index 到 macro cell index。
- `macro_to_positive[j] -> p or -1`：macro cell 到 positive local index。
- BV source 用 positive local index `p` 访问 `phi_s[p]` 与 `c_s[p, -1]`，用 macro index `j` 访问 `c_e[j]` 与 `phi_e[j]`。

### 4.4 参数容器

Phase 0 不复制现有材料 dataclass，而是组合现有配置：

- `ActiveMaterial`：`cs_max`, `sigma`, `diffusivity`, `ocv`, `ocv_derivative`。
- `Electrolyte`：`diffusivity`, `conductivity`, `c_init`, `transference_number`, `bruggeman`。
- `Kinetics`：`k0`, `alpha_a`, `alpha_c`。
- `Separator`：`thickness`, `porosity`, `bruggeman`, `t_plus`, `c_ref`, Li foil 开关。
- `SolverSettings`：`temperature`, `dt`, `newton_tol`, `newton_max_iter`。

P2D 新增参数只描述 DFN 几何和数值布局：

- `P2DRegion`：separator 与 positive electrode 的厚度、cell 数、体积分数、Bruggeman 指数、specific area。
- `particle_shells`：径向 shell 数。
- `area`：几何截面积。
- `potential_gauge`：Phase 0 固定为 `li_phi_e_left`。
- `concentration_floor`：用于 \(\ln c_e\)、BV 与 Newton trial 检查。

## 5. 验证计划

### 5.1 守恒检查

质量守恒：

- 对 electrolyte：

\[
\sum_j V_j\epsilon_{e,j}\frac{c_{e,j}^{n+1}-c_{e,j}^{n}}{\Delta t}
= \sum_{\partial\Omega} A_f N_{e,f} + \sum_{j\in p}V_j(1-t_+^0)\frac{a_{s,j}i_{F,j}}{F}
\]

- 对 solid particle：

\[
V_j\epsilon_{s,j}\frac{\bar c_{s,j}^{n+1}-\bar c_{s,j}^{n}}{\Delta t}
= -V_j\frac{a_{s,j}i_{F,j}}{F}
\]

测试：

- `tests/test_p2d_conservation.py::test_solid_particle_mass_balance`
- `tests/test_p2d_conservation.py::test_electrolyte_salt_balance`
- no-flux electrolyte boundary 下 total lithium inventory 误差随 Newton tolerance 下降。

电荷守恒：

\[
\max_f |i_{e,f}+i_{s,f}-I_\mathrm{app}| \le \mathrm{tol}
\]

测试：

- `tests/test_p2d_conservation.py::test_total_current_conservation`
- separator 中 \(i_s=0\)，因此 \(i_e=I_\mathrm{app}\)；positive current collector 处 \(i_e=0\)，因此 \(i_s=I_\mathrm{app}\)。

### 5.2 平衡态检查

条件：

- \(I_\mathrm{app}=0\)。
- \(c_e=c_{e,0}\) uniform。
- \(c_s=c_{s,0}\) uniform。
- \(\phi_e=0\)。
- \(\phi_s=U(c_{s,0}/c_s^\max)\)。

期望：

- \(\eta=0\)。
- \(i_F=0\)。
- 所有 concentration residual 为 0。
- voltage 等于 OCV。

测试：

- `tests/test_p2d_residual.py::test_open_circuit_equilibrium_residual`
- `tests/test_p2d_solver.py::test_zero_current_state_does_not_drift`

### 5.3 解析解或半解析对比

电解质扩散：

- 关闭反应项，设置固定边界浓度或初始 sinusoidal perturbation。
- 对比 1D Fick 解析解：

\[
c(x,t)-c_0 = A\sin(\pi x/L)\exp(-D^\mathrm{eff}\pi^2t/L^2)
\]

固相球扩散：

- 恒定表面通量、恒定 \(D_s\)，对比短时近似或高精度 `solve_ivp`/细网格参考。
- 检查体积平均浓度满足：

\[
\frac{d\bar c_s}{dt}=-\frac{3i_F}{FR_p}
\]

Ohmic conduction：

- 关闭 BV 或设置均匀 source，验证 slab 中 \(\phi\) 线性或二次 profile。

### 5.4 PyBaMM 或文献对比

Phase 0 结束时至少准备一个可运行对比路径：

- 使用公开 DFN 参数集或项目 Khan 2021 参数。
- 与 PyBaMM DFN half-cell/full-cell 近似配置比较：
  - \(V(t)\)
  - \(c_e(x,t)\)
  - \(c_s^\mathrm{surf}(x,t)\)
  - current distribution
- 若 PyBaMM 未安装，测试标记为 optional，脚本保留在 `examples/p2d/compare_pybamm.py`。

误差目标：

- 参数完全一致时 voltage RMS error 目标 `5-10 mV`。
- Phase 0 若只做 half-cell 且边界条件不同，先要求趋势一致，并在结果 metadata 中记录差异。

### 5.5 边界情况

必须测试：

- 零电流：`current_density=0.0`。
- 极小倍率：例如 `abs(I_app)` 为 1C 的 `1e-4`，电压接近 OCV。
- 常规倍率：`0.2C`, `1C` smoke tests。
- 高倍率：`3C` 或 `5C`，允许更小 dt，但不能出现负浓度或 Newton silent failure。
- 极小 separator thickness 或极小 porosity 应给出清晰 `ValueError` 或 solver failure report。
- `sigma_s=0`、`kappa=0` 等不可导通参数应在参数校验或 factorization failure 中明确报错。

## 6. 与现有代码的接口

### 6.1 复用 `pnmcathode.physics`

必须复用：

- `pnmcathode.physics.electrolyte.electrolyte_diffusion_coefficient(c_e, T)`：
  - 用于 \(D_e(c_e,T)\)。
  - 在 `materials.bruggeman_diffusivity_e()` 中乘 \(\epsilon_e^{b_D}\)。
- `pnmcathode.physics.electrolyte.electrolyte_ionic_conductivity(c_e, T)`：
  - 用于 \(\kappa(c_e,T)\)。
  - 在 `materials.bruggeman_conductivity_e()` 中乘 \(\epsilon_e^{b_\kappa}\)。
- `pnmcathode.physics.solid.nmc532_diffusion_coefficient(c_s, T)`：
  - 用于 \(D_s(c_s,T)\)。
  - Phase 0 residual 可先用 shell/cell average \(D_s\)，后续补 face nonlinear averaging。
- `pnmcathode.physics.reaction.butler_volmer()` 与 `exchange_current_density()`：
  - P2D 的 `kinetics.py` 与现有标量函数保持同一符号约定。
  - 单元测试对比 vectorized 结果与现有函数。
- `pnmcathode.physics.reaction.F` 与 `R`：
  - 不重复定义 Faraday constant 与 gas constant。
- `pnmcathode.physics.ocv.nmc532_ocv()` 与 `ocv_derivative()`：
  - 用于 empirical OCV vs Li/Li+。
  - 注意 `ocv_derivative()` 当前对 SOC 求导，P2D Jacobian 需要：

\[
\frac{\partial U}{\partial c_s}=\frac{1}{c_s^\max}\frac{dU}{d\mathrm{SOC}}
\]

### 6.2 复用 `pnmcathode.config`

P2D API 示例：

```python
from pnmcathode.config import Kinetics, SolverSettings
from pnmcathode.materials.presets import (
    electrolyte_khan2021,
    nmc532_khan2021,
    separator_khan2021,
)
from pnmcathode.p2d import P2DRegion, P2DSolver, from_config

active = nmc532_khan2021()
electrolyte = electrolyte_khan2021()
separator = separator_khan2021(enabled=True)
kinetics = Kinetics(k0=5e-10)
settings = SolverSettings(temperature=303.0, newton_tol=1e-8)

positive = P2DRegion(
    name="positive",
    x_left=separator.thickness,
    x_right=separator.thickness + 75e-6,
    n_cells=40,
    epsilon_e=0.35,
    epsilon_s=0.55,
    particle_radius=5e-6,
    bruggeman_e=electrolyte.bruggeman,
    bruggeman_s=1.5,
    sigma_s=active.sigma,
)

params = from_config(
    active=active,
    electrolyte=electrolyte,
    kinetics=kinetics,
    separator=separator,
    settings=settings,
    positive=positive,
    particle_shells=20,
)
solver = P2DSolver(params, settings)
state0 = solver.initial_state(soc0=0.5)
```

设计原则：

- `config` 继续作为 material/solver 参数真源。
- `p2d.materials.P2DParameters` 只组合与补充 P2D 几何/离散参数。
- 不在 P2D 模块中硬编码 Khan 2021 参数；示例与测试可使用 `materials.presets`。

### 6.3 注册到 `pnmcathode.__init__.py`

修改 `_LAZY_EXPORTS`：

```python
_LAZY_EXPORTS = {
    # existing exports ...
    "P2DRegion": ("pnmcathode.p2d", "P2DRegion"),
    "P2DParameters": ("pnmcathode.p2d", "P2DParameters"),
    "P2DProtocol": ("pnmcathode.p2d", "P2DProtocol"),
    "P2DResult": ("pnmcathode.p2d", "P2DResult"),
    "P2DSolver": ("pnmcathode.p2d", "P2DSolver"),
}
```

`src/pnmcathode/p2d/__init__.py`：

```python
from pnmcathode.p2d.domain import MacroMesh, P2DRegion, ParticleMesh
from pnmcathode.p2d.materials import OCVModel, P2DMaterial, P2DParameters, from_config
from pnmcathode.p2d.results import P2DResult, P2DSnapshot
from pnmcathode.p2d.solver import P2DProtocol, P2DSolver
from pnmcathode.p2d.state import P2DState

__all__ = [
    "MacroMesh",
    "OCVModel",
    "P2DMaterial",
    "P2DParameters",
    "P2DProtocol",
    "P2DRegion",
    "P2DResult",
    "P2DSnapshot",
    "P2DSolver",
    "P2DState",
    "ParticleMesh",
    "from_config",
]
```

### 6.4 与现有 solver 的边界

现有 `pnmcathode.solver.transient.TransientSolver` 是 PNM transient solver，Phase 0 的 `pnmcathode.p2d.solver.P2DSolver` 不应依赖 OpenPNM，也不应修改 PNM solver 行为。

接口边界：

- P2D 与 PNM 在 Phase 0 完全分离。
- 共享材料函数与 `config` dataclass。
- 共享符号约定：anodic convention，cathode discharge 为 \(I_\mathrm{app}<0\)、\(i_F<0\)。
- 后续 Phase 2/3 通过 `EffectivePropertyProvider` 才让 PNM effective properties 注入 P2D。

## 7. 推荐实现顺序

1. 新建 `domain.py`、`state.py`、`materials.py`，完成网格、layout、参数适配测试。
2. 新建 `kinetics.py`，完成 BV vectorization 与导数 finite-difference 校验。
3. 实现 `residual.py` 的纯 residual，不先做 time loop；用 equilibrium residual 测试锁定符号。
4. 实现 `jacobian.py` 的 finite-difference fallback，再逐步替换为 analytic sparse Jacobian。
5. 实现 `solver.py` 的单步 Newton，先跑 \(I=0\)，再跑极小电流。
6. 实现 `diagnostics.py` 与 conservation tests。
7. 实现 `results.py`、example、顶层 export。
8. 增加 low-rate 与 high-rate smoke tests，记录收敛迭代数与 voltage behavior。

## 8. 风险清单与硬性约束

- 符号约定必须在测试中固定：cathode discharge 输入若用 C-rate magnitude，solver 内部必须转换为 \(I_\mathrm{app}<0\)。
- empirical NMC OCV 已相对 Li/Li+，Phase 0 默认不再加入 electrolyte Nernst correction 到 \(U\)。
- potential gauge 必须显式 residual 替换，否则 Jacobian 奇异。
- BV 指数必须 clip，但 clip 区间不应掩盖非物理解；Newton report 要记录最大 \(|\eta|\)。
- 浓度 trial 不能直接静默 clip 为 accepted state；应通过 damping 保持物理 bounds。
- PNM effective property 与 Bruggeman 不能混用；Phase 0 只有 `closure="bruggeman"`。
- 所有 tests 应使用小网格，避免把数值稳定性问题隐藏在长运行 benchmark 中。
