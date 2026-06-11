"""
P2D 参数容器与 Bruggeman 有效物性 (Materials & Effective Properties)
====================================================================

职责:
- 将现有 pnmcathode.config 参数适配到 P2D。
- 实现 Bruggeman 有效传输物性闭合。
- 管理 OCV reference，避免重复加入电解质 Nernst 修正。

对应 TASK_PHASE0.md §1.4, §6.1
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np

from pnmcathode.config import ActiveMaterial, Electrolyte, Kinetics, Separator, SolverSettings
from pnmcathode.p2d.domain import P2DRegion

Array = np.ndarray


@dataclass(frozen=True)
class OCVModel:
    """OCV 模型。

    Attributes
    ----------
    reference : str
        OCV 参考类型。
    value : callable
        U(soc) → [V]，soc = c_s / c_s_max。
    derivative : callable or None
        dU/dsoc(soc) → [V]。
    """

    reference: Literal["li_metal_empirical", "thermodynamic_with_electrolyte"]
    value: Callable[[Array], Array]
    derivative: Callable[[Array], Array] | None = None


@dataclass(frozen=True)
class P2DMaterial:
    """P2D 材料参数集合。

    组合现有 config dataclass，不复制参数。

    Attributes
    ----------
    active : ActiveMaterial
        活性固相材料。
    electrolyte : Electrolyte
        电解质传输模型。
    kinetics : Kinetics
        Butler-Volmer 动力学参数。
    separator : Separator
        隔膜参数。
    temperature : float
        温度 [K]。
    ocv_model : OCVModel
        OCV 模型。
    """

    active: ActiveMaterial
    electrolyte: Electrolyte
    kinetics: Kinetics
    separator: Separator
    temperature: float
    ocv_model: OCVModel


@dataclass(frozen=True)
class P2DParameters:
    """P2D 模型完整参数。

    Attributes
    ----------
    material : P2DMaterial
        材料参数。
    positive : P2DRegion
        positive electrode 区域。
    separator_region : P2DRegion
        separator 区域。
    particle_shells : int
        颗粒径向壳层数。
    area : float
        几何截面积 [m²]。
    concentration_floor : float
        浓度下限 [mol/m³]。
    potential_gauge : str
        电势规范类型。
    """

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
    n_separator_cells: int | None = None,
) -> P2DParameters:
    """从现有 config 对象构建 P2D 参数。

    Parameters
    ----------
    active : ActiveMaterial
        活性材料参数。
    electrolyte : Electrolyte
        电解质参数。
    kinetics : Kinetics
        动力学参数。
    separator : Separator
        隔膜参数。
    settings : SolverSettings
        求解器设置 (提取温度)。
    positive : P2DRegion
        positive electrode 区域定义。
    particle_shells : int
        颗粒壳层数。
    area : float
        截面积 [m²]。

    Returns
    -------
    P2DParameters
    """
    from pnmcathode.physics.ocv import nmc532_ocv, ocv_derivative

    # 构建 OCV 模型
    def _ocv_value(soc: Array) -> Array:
        soc = np.asarray(soc, dtype=float)
        return np.asarray(nmc532_ocv(soc), dtype=float)

    def _ocv_derivative(soc: Array) -> Array:
        soc = np.asarray(soc, dtype=float)
        return np.asarray(ocv_derivative(soc), dtype=float)

    ocv_model = OCVModel(
        reference="li_metal_empirical",
        value=_ocv_value,
        derivative=_ocv_derivative,
    )

    material = P2DMaterial(
        active=active,
        electrolyte=electrolyte,
        kinetics=kinetics,
        separator=separator,
        temperature=settings.temperature,
        ocv_model=ocv_model,
    )

    # 构建 separator 区域
    sep_n_cells = (
        n_separator_cells
        if n_separator_cells is not None
        else max(1, positive.n_cells // 4)
    )
    separator_region = P2DRegion(
        name="separator",
        x_left=0.0,
        x_right=separator.thickness,
        n_cells=sep_n_cells,
        epsilon_e=separator.porosity,
        bruggeman_e=separator.bruggeman,
    )

    return P2DParameters(
        material=material,
        positive=positive,
        separator_region=separator_region,
        particle_shells=particle_shells,
        area=area,
    )


def bruggeman_diffusivity_e(
    c_e: Array,
    region: P2DRegion,
    material: P2DMaterial,
) -> Array:
    """计算 Bruggeman 有效电解质扩散系数。

    D_e^eff = D_e(c_e, T) * epsilon_e^b_D

    Parameters
    ----------
    c_e : ndarray
        电解质浓度 [mol/m³]。
    region : P2DRegion
        区域定义。
    material : P2DMaterial
        材料参数。

    Returns
    -------
    D_e_eff : ndarray
        有效扩散系数 [m²/s]。
    """
    D_e = np.asarray(
        material.electrolyte.diffusivity(c_e, material.temperature), dtype=float
    )
    return D_e * region.epsilon_e ** region.bruggeman_e


def bruggeman_conductivity_e(
    c_e: Array,
    region: P2DRegion,
    material: P2DMaterial,
) -> Array:
    """计算 Bruggeman 有效电解质电导率。

    kappa^eff = kappa(c_e, T) * epsilon_e^b_kappa

    Parameters
    ----------
    c_e : ndarray
        电解质浓度 [mol/m³]。
    region : P2DRegion
        区域定义。
    material : P2DMaterial
        材料参数。

    Returns
    -------
    kappa_eff : ndarray
        有效电导率 [S/m]。
    """
    kappa = np.asarray(
        material.electrolyte.conductivity(c_e, material.temperature), dtype=float
    )
    return kappa * region.epsilon_e ** region.bruggeman_e


def bruggeman_conductivity_s(
    region: P2DRegion,
    material: P2DMaterial,
) -> float:
    """计算 Bruggeman 有效固相电导率。

    sigma^eff = sigma * epsilon_s^b_sigma

    Parameters
    ----------
    region : P2DRegion
        区域定义。
    material : P2DMaterial
        材料参数。

    Returns
    -------
    sigma_eff : float
        有效固相电导率 [S/m]。
    """
    if region.sigma_s is None:
        return 0.0
    return region.sigma_s * region.epsilon_s ** region.bruggeman_s
