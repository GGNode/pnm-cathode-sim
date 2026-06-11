"""
P2D Newton-Raphson 求解器 (Newton Solver)
=========================================

职责:
- Newton-Raphson 非线性求解。
- Backward Euler 时间步进。
- 恒流 protocol run。
- 自适应时间步与步拒绝。

Newton 默认参数:
- newton_tol = 1e-8 (scaled residual norm)
- newton_max_iter = 30
- damping: 从 1.0 开始，依次尝试 1/2, 1/4, 1/8, 1/16
- variable bounds: c_e > c_floor, 0 < c_s < c_s_max

对应 TASK_PHASE0.md §1.8, §3.3
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.sparse.linalg import spsolve

from pnmcathode.config import SolverSettings
from pnmcathode.p2d.domain import MacroMesh, ParticleMesh
from pnmcathode.p2d.jacobian import assemble_p2d_jacobian
from pnmcathode.p2d.materials import P2DParameters
from pnmcathode.p2d.residual import P2DResidualContext, assemble_p2d_residual, residual_norm
from pnmcathode.p2d.results import P2DResult, P2DSnapshot
from pnmcathode.p2d.state import (
    P2DState,
    StateLayout,
    make_initial_state,
    pack_state,
    unpack_state,
)

Array = np.ndarray


@dataclass(frozen=True)
class P2DProtocol:
    """恒流放电协议。

    Attributes
    ----------
    current_density : float
        应用电流密度 [A/m²] (阳极约定, 放电 < 0)。
    t_final : float
        终止时间 [s]。
    dt_initial : float
        初始时间步 [s]。
    cutoff_voltage : float or None
        截止电压 [V]。
    save_every : int
        每隔多少步保存快照。
    """

    current_density: float
    t_final: float
    dt_initial: float
    cutoff_voltage: float | None = None
    save_every: int = 1


@dataclass(frozen=True)
class NewtonReport:
    """Newton 迭代报告。

    Attributes
    ----------
    converged : bool
        是否收敛。
    iterations : int
        迭代次数。
    residual_norm : float
        最终残差范数。
    step_norm : float
        最终步长范数。
    """

    converged: bool
    iterations: int
    residual_norm: float
    step_norm: float


class P2DSolver:
    """P2D half-cell 求解器。

    使用 monolithic Newton-Raphson 求解 c_e, phi_e, phi_s, c_s 的耦合系统。
    时间推进使用 Backward Euler (BDF1)。

    对应 TASK_PHASE0.md §1.8
    """

    def __init__(
        self,
        params: P2DParameters,
        settings: SolverSettings | None = None,
        use_analytic_jacobian: bool = True,
    ) -> None:
        """初始化 P2D 求解器。

        Parameters
        ----------
        params : P2DParameters
            模型参数。
        settings : SolverSettings or None
            求解器设置。None 时使用默认值。
        use_analytic_jacobian : bool
            是否使用解析 Jacobian (Phase 0 暂时使用有限差分)。
        """
        self.params = params
        self.settings = settings or SolverSettings()
        self.use_analytic_jacobian = use_analytic_jacobian

        # 构建网格
        self._macro = MacroMesh.from_regions(
            (params.separator_region, params.positive),
            area=params.area,
        )
        self._particle = ParticleMesh.spherical(
            params.positive.particle_radius,
            params.particle_shells,
        )
        self._layout = StateLayout.from_meshes(self._macro, self._particle)

        # Newton 参数
        self._newton_tol = self.settings.newton_tol
        self._newton_max_iter = min(self.settings.newton_max_iter, 30)

    @property
    def macro(self) -> MacroMesh:
        """宏观网格。"""
        return self._macro

    @property
    def particle(self) -> ParticleMesh:
        """颗粒网格。"""
        return self._particle

    @property
    def layout(self) -> StateLayout:
        """状态布局。"""
        return self._layout

    def initial_state(self, soc0: float, c_e0: float | None = None) -> P2DState:
        """构建初始状态。

        Parameters
        ----------
        soc0 : float
            初始 SOC (0~1)。
        c_e0 : float or None
            初始电解质浓度 [mol/m³]。None 时使用 electrolyte.c_init。

        Returns
        -------
        P2DState
        """
        if c_e0 is None:
            c_e0 = self.params.material.electrolyte.c_init
        return make_initial_state(
            macro=self._macro,
            particle=self._particle,
            c_e0=c_e0,
            soc0=soc0,
            c_s_max=self.params.material.active.cs_max,
        )

    def step(
        self,
        state: P2DState,
        dt: float,
        current_density: float,
    ) -> tuple[P2DState, NewtonReport]:
        """推进一个时间步。

        使用 Newton-Raphson 求解:
            R(y^{n+1}) = 0
        其中 y^{n+1} 是 Backward Euler 残差。

        Parameters
        ----------
        state : P2DState
            当前状态。
        dt : float
            时间步长 [s]。
        current_density : float
            电流密度 [A/m²] (阳极约定)。

        Returns
        -------
        new_state : P2DState
            新状态。
        report : NewtonReport
            Newton 迭代报告。
        """
        context = P2DResidualContext(
            macro=self._macro,
            particle=self._particle,
            layout=self._layout,
            params=self.params,
            dt=dt,
            current_density=current_density,
        )

        y = pack_state(state, self._layout)
        cs_max = self.params.material.active.cs_max
        c_floor = self.params.concentration_floor

        converged = False
        residual_norm_val = float("inf")
        step_norm_val = float("inf")

        for iteration in range(self._newton_max_iter):
            R = assemble_p2d_residual(y, state, context)
            residual_norm_val = residual_norm(R, context)

            if residual_norm_val < self._newton_tol:
                converged = True
                break

            # 组装 Jacobian
            J = assemble_p2d_jacobian(y, state, context)

            # 线性求解
            try:
                J_csc = J.tocsc()
                dy = spsolve(J_csc, -R)
            except Exception:
                break

            if not np.all(np.isfinite(dy)):
                break

            # Damping with bounds checking
            damping = 1.0
            for _ in range(5):
                y_trial = y + damping * dy
                state_trial = unpack_state(y_trial, self._layout)

                # 检查浓度 bounds
                c_e_ok = np.all(state_trial.c_e > c_floor)
                c_s_ok = (
                    np.all(state_trial.c_s > 0)
                    and np.all(state_trial.c_s < cs_max)
                )

                if c_e_ok and c_s_ok:
                    break
                damping *= 0.5

            step_norm_val = float(np.max(np.abs(damping * dy)))
            y = y + damping * dy

        new_state = unpack_state(y, self._layout)
        new_state.time = state.time + dt

        report = NewtonReport(
            converged=converged,
            iterations=iteration + 1 if converged else self._newton_max_iter,
            residual_norm=residual_norm_val,
            step_norm=step_norm_val,
        )
        return new_state, report

    def run(
        self,
        state0: P2DState,
        protocol: P2DProtocol,
        progress: Callable[[P2DState], None] | None = None,
    ) -> P2DResult:
        """运行恒流放电仿真。

        Parameters
        ----------
        state0 : P2DState
            初始状态。
        protocol : P2DProtocol
            放电协议。
        progress : callable or None
            进度回调。

        Returns
        -------
        P2DResult
        """
        I_app = protocol.current_density
        dt = protocol.dt_initial
        t_final = protocol.t_final
        save_every = protocol.save_every

        state = state0.copy()
        snapshots = []

        # 初始快照
        V0 = self.voltage(state)
        snapshots.append(P2DSnapshot(
            time=0.0,
            state=state.copy(),
            voltage=V0,
            current_density=I_app,
        ))

        step_count = 0
        t = 0.0
        dt_min = self.settings.dt_min
        dt_max = self.settings.dt_max
        growth = 1.25
        shrink = 0.5

        while t < t_final:
            # 自适应时间步
            dt_eff = min(dt, t_final - t)
            if dt_eff <= 0:
                break

            # 尝试步进
            new_state, report = self.step(state, dt_eff, I_app)

            if report.converged:
                state = new_state
                t += dt_eff
                step_count += 1

                # 电压检查
                V = self.voltage(state)
                if protocol.cutoff_voltage is not None and V < protocol.cutoff_voltage:
                    break

                # 保存快照
                if step_count % save_every == 0:
                    snapshots.append(P2DSnapshot(
                        time=t,
                        state=state.copy(),
                        voltage=V,
                        current_density=I_app,
                    ))

                # 自适应 dt
                if report.iterations <= 5:
                    dt = min(dt * growth, dt_max)
                elif report.iterations >= 15:
                    dt = max(dt * shrink, dt_min)

                if progress is not None:
                    progress(state)
            else:
                # 步拒绝: 缩小 dt
                if dt <= dt_min:
                    break
                dt = max(dt * shrink, dt_min)

        # 最终快照
        V_final = self.voltage(state)
        if len(snapshots) == 0 or snapshots[-1].time < t:
            snapshots.append(P2DSnapshot(
                time=t,
                state=state.copy(),
                voltage=V_final,
                current_density=I_app,
            ))

        return P2DResult(
            snapshots=snapshots,
            metadata={
                "current_density": I_app,
                "t_final": t,
                "steps": step_count,
                "final_voltage": V_final,
            },
        )

    def voltage(self, state: P2DState) -> float:
        """计算电池电压。

        Half-cell: V = phi_s(last positive cell) - phi_e(0)
        由于 phi_e(0) = 0 (gauge), V = phi_s[-1]。

        Parameters
        ----------
        state : P2DState
            当前状态。

        Returns
        -------
        V : float
            电池电压 [V]。
        """
        return float(state.phi_s[-1])
