"""P2D (Doyle-Fuller-Newman) 半电池求解器。

Phase 0: separator | porous cathode half-cell 1D DFN solver。
不依赖 PNM，使用 Bruggeman 有效物性。
"""

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
