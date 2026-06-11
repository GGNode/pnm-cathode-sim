"""物理模块: 电解质传输、固相扩散、反应动力学、OCV、隔膜边界。

子模块:
- electrolyte: 电解质扩散系数和离子电导率
- solid: NMC532 固相锂扩散 (球坐标 Fick)
- reaction: Butler-Volmer 反应动力学
- ocv: NMC532 开路电压曲线
- separator: 一维隔膜边界模型
"""

from pnmcathode.physics.electrolyte import (
    electrolyte_diffusion_coefficient,
    electrolyte_ionic_conductivity,
)
from pnmcathode.physics.ocv import nmc532_ocv, ocv_derivative
from pnmcathode.physics.reaction import F, R, butler_volmer, exchange_current_density
from pnmcathode.physics.separator import (
    SeparatorParams,
    SeparatorState,
    separator_boundary,
)
from pnmcathode.physics.solid import (
    discretize_spherical_particle,
    nmc532_diffusion_coefficient,
    solid_diffusion_rhs,
)

__all__ = [
    "F",
    "R",
    "SeparatorParams",
    "SeparatorState",
    "butler_volmer",
    "discretize_spherical_particle",
    "electrolyte_diffusion_coefficient",
    "electrolyte_ionic_conductivity",
    "exchange_current_density",
    "nmc532_diffusion_coefficient",
    "nmc532_ocv",
    "ocv_derivative",
    "separator_boundary",
    "solid_diffusion_rhs",
]
