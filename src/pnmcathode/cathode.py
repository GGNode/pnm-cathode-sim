"""Cathode 对象与网络工厂方法。"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pnmcathode.config import (
    ActiveMaterial,
    CathodeGeometry,
    ConductiveAdditive,
    Electrolyte,
)
from pnmcathode.network.generator import create_cathode_network


@dataclass
class Cathode:
    """组合阴极几何、材料模型与 OpenPNM 网络。"""

    geometry: CathodeGeometry
    active_material: ActiveMaterial
    electrolyte: Electrolyte
    conductive_additive: ConductiveAdditive = field(default_factory=ConductiveAdditive)
    network: Any | None = None

    @classmethod
    def cubic(
        cls,
        geometry: CathodeGeometry | None = None,
        active_material: ActiveMaterial | None = None,
        electrolyte: Electrolyte | None = None,
        conductive_additive: ConductiveAdditive | None = None,
    ) -> "Cathode":
        """创建合成立方阴极网络。"""

        from pnmcathode.materials.presets import (
            cbd_khan2021,
            electrolyte_khan2021,
            nmc532_khan2021,
        )

        cathode = cls(
            geometry=geometry or CathodeGeometry(),
            active_material=active_material or nmc532_khan2021(),
            electrolyte=electrolyte or electrolyte_khan2021(),
            conductive_additive=conductive_additive or cbd_khan2021(),
        )
        cathode.build_network()
        return cathode

    def build_network(self) -> Any:
        """按配置几何生成并绑定网络。"""

        net = create_cathode_network(
            shape=list(self.geometry.shape),
            spacing=self.geometry.spacing,
            porosity=self.geometry.porosity,
            cbd_fraction=self.geometry.cbd_fraction,
            seed=self.geometry.seed,
            throat_scale=self.geometry.throat_scale,
        )
        self._apply_material_properties(net)
        self.network = net
        return net

    def ensure_network(self) -> Any:
        """返回已绑定网络；必要时先生成网络。"""

        if self.network is None:
            return self.build_network()
        return self.network

    def _apply_material_properties(self, net: Any) -> None:
        nmc_mask = np.asarray(net["pore.nmc"], dtype=bool)
        cbd_mask = np.asarray(net["pore.cbd"], dtype=bool)
        net["pore.cs_max"] = np.where(nmc_mask, self.active_material.cs_max, 0.0)
        net["pore.sigma_nmc"] = np.where(nmc_mask, self.active_material.sigma, 0.0)
        net["pore.sigma_cbd"] = np.where(cbd_mask, self.conductive_additive.sigma, 0.0)
