"""Compatibility wrappers for legacy ``src.solver`` imports."""

from pnmcathode.solver.single_pore import SinglePoreDischarge
from pnmcathode.solver.steady import SteadyStateSolver
from pnmcathode.solver.transient import TransientSolver

__all__ = ["SinglePoreDischarge", "SteadyStateSolver", "TransientSolver"]
