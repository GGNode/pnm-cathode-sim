"""求解器模块: 瞬态放电、稳态电位、单孔仿真。

子模块:
- steady: 稳态电位求解器 (耦合 Newton-Raphson)
- transient: 瞬态放电求解器 (backward Euler + 自适应步进)
- single_pore: 单孔放电仿真器 (单粒子模型)
"""

from pnmcathode.solver.single_pore import SinglePoreDischarge
from pnmcathode.solver.steady import SteadyStateSolver
from pnmcathode.solver.transient import TransientSolver

SteadySolver = SteadyStateSolver
SinglePoreSolver = SinglePoreDischarge

__all__ = [
    "SinglePoreDischarge",
    "SinglePoreSolver",
    "SteadySolver",
    "SteadyStateSolver",
    "TransientSolver",
]
