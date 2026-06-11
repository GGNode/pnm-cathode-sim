"""可复用的阴极孔网络仿真工具包。"""

from __future__ import annotations

from importlib import import_module
from typing import Any

__version__ = "0.1.0"

_LAZY_EXPORTS = {
    "ActiveMaterial": ("pnmcathode.config", "ActiveMaterial"),
    "Cathode": ("pnmcathode.cathode", "Cathode"),
    "CathodeGeometry": ("pnmcathode.config", "CathodeGeometry"),
    "ConductiveAdditive": ("pnmcathode.config", "ConductiveAdditive"),
    "DischargeProtocol": ("pnmcathode.config", "DischargeProtocol"),
    "DischargeResult": ("pnmcathode.results", "DischargeResult"),
    "Electrolyte": ("pnmcathode.config", "Electrolyte"),
    "Kinetics": ("pnmcathode.config", "Kinetics"),
    "Separator": ("pnmcathode.config", "Separator"),
    "Simulation": ("pnmcathode.simulation", "Simulation"),
    "SolverSettings": ("pnmcathode.config", "SolverSettings"),
}


def __getattr__(name: str) -> Any:
    """按需加载顶层公共对象，避免主包导入时触发重依赖。"""
    try:
        module_name, attribute_name = _LAZY_EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module 'pnmcathode' has no attribute {name!r}") from exc
    module = import_module(module_name)
    value = getattr(module, attribute_name)
    globals()[name] = value
    return value


__all__ = [
    "__version__",
    *_LAZY_EXPORTS,
]
