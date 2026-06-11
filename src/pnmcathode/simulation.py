"""高层仿真流程编排。"""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

from pnmcathode.cathode import Cathode
from pnmcathode.config import DischargeProtocol, Kinetics, Separator, SolverSettings
from pnmcathode.physics.separator import SeparatorParams
from pnmcathode.results import DischargeResult
from pnmcathode.solver.transient import TransientSolver


def _separator_params(separator: Separator) -> SeparatorParams:
    return SeparatorParams(
        enabled=separator.enabled,
        thickness=separator.thickness,
        porosity=separator.porosity,
        bruggeman=separator.bruggeman,
        t_plus=separator.t_plus,
        c_ref=separator.c_ref,
        c_floor=separator.c_floor,
        include_concentration_overpotential=separator.include_concentration_overpotential,
        include_li_foil_bv=separator.include_li_foil_bv,
        i0_foil=separator.i0_foil,
        alpha_foil=separator.alpha_foil,
    )


class Simulation:
    """高层恒流放电仿真对象。"""

    def __init__(
        self,
        cathode: Cathode,
        protocol: DischargeProtocol,
        kinetics: Kinetics | None = None,
        separator: Separator | None = None,
        settings: SolverSettings | None = None,
    ):
        self.cathode = cathode
        self.protocol = protocol
        self.kinetics = kinetics or Kinetics()
        self.separator = separator or Separator(enabled=False)
        self.settings = settings or SolverSettings()

    def _metadata(self) -> dict[str, Any]:
        return {
            "geometry": asdict(self.cathode.geometry),
            "active_material": self.cathode.active_material.name,
            "electrolyte": self.cathode.electrolyte.name,
            "conductive_additive": self.cathode.conductive_additive.name,
            "protocol": asdict(self.protocol),
            "kinetics": asdict(self.kinetics),
            "separator": asdict(self.separator),
            "settings": asdict(self.settings),
        }

    def run(self) -> DischargeResult:
        """运行放电仿真并返回类型化结果。"""

        net = self.cathode.ensure_network()
        solver = TransientSolver(
            net,
            T=self.settings.temperature,
            k0=self.kinetics.k0,
            geometric_area=self.cathode.geometry.geometric_area,
            separator=_separator_params(self.separator),
        )
        solver.set_concentration(
            c_e=self.cathode.electrolyte.c_init,
            c_s=self.protocol.initial_soc * self.cathode.active_material.cs_max,
        )
        raw = solver.run_discharge(
            C_rate=self.protocol.c_rate,
            cutoff_voltage=self.protocol.cutoff_voltage,
            dt=self.settings.dt,
            dt_min=self.settings.dt_min,
            dt_max=self.settings.dt_max,
            max_steps=self.protocol.max_steps,
            voltage_jump_limit=self.settings.voltage_jump_limit,
        )
        return DischargeResult.from_solver_output(raw, metadata=self._metadata())
