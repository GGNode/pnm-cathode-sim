"""
向量化 Butler-Volmer 动力学与解析导数 (Vectorized Kinetics)
============================================================

职责:
- 为 P2D residual 提供向量化 Butler-Volmer 反应电流。
- 计算 exchange current density 及其浓度导数。
- 计算 i_F 对 eta, c_e, c_s_surf 的解析导数 (Jacobian 系数)。
- 复用 pnmcathode.physics.reaction 的符号约定。

对应 TASK_PHASE0.md §1.5, 公式 2.5
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pnmcathode.config import Kinetics
from pnmcathode.physics.reaction import F, R

Array = np.ndarray


@dataclass(frozen=True)
class ExchangeCurrentResult:
    """交换电流密度及其导数。

    Attributes
    ----------
    i0 : ndarray
        交换电流密度 [A/m²]。
    di0_dc_surf : ndarray
        di0/d(c_s_surf) [A/m² / (mol/m³)]。
    di0_dc_e : ndarray
        di0/d(c_e) [A/m² / (mol/m³)]。
    """

    i0: Array
    di0_dc_surf: Array
    di0_dc_e: Array


@dataclass(frozen=True)
class ReactionRates:
    """Butler-Volmer 反应速率及其全部导数。

    Attributes
    ----------
    eta : ndarray
        过电位 [V]。
    i0 : ndarray
        交换电流密度 [A/m²]。
    i_f : ndarray
        Faradaic 电流密度 [A/m²]。
    di_deta : ndarray
        di_F/d(eta) [A/(m²·V)]。
    di_dc_e : ndarray
        di_F/d(c_e) [A/m² / (mol/m³)]。
    di_dc_surf : ndarray
        di_F/d(c_s_surf) [A/m² / (mol/m³)]。
    """

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
) -> ExchangeCurrentResult:
    """向量化计算交换电流密度及其浓度导数。

    i0 = F * k0 * c_e^gamma_e * (c_s_max - c_s_surf)^gamma_v * c_s_surf^gamma_s

    其中 gamma_e = alpha_a, gamma_v = alpha_a, gamma_s = alpha_c。

    对应 TASK_PHASE0.md 公式 2.5

    Parameters
    ----------
    c_e : ndarray
        电解质浓度 [mol/m³]。
    c_s_surf : ndarray
        固相表面浓度 [mol/m³]。
    c_s_max : float
        固相最大浓度 [mol/m³]。
    kinetics : Kinetics
        动力学参数。
    floor : float
        浓度下限，防止奇异点。

    Returns
    -------
    ExchangeCurrentResult
    """
    c_e = np.asarray(c_e, dtype=float)
    c_s_surf = np.asarray(c_s_surf, dtype=float)

    ce = np.maximum(c_e, floor)
    cs = np.clip(c_s_surf, floor, c_s_max - floor)
    vacancies = c_s_max - cs

    k0 = kinetics.k0
    alpha_a = kinetics.alpha_a
    alpha_c = kinetics.alpha_c

    # i0 = F * k0 * ce^αa * vacancies^αa * cs^αc
    i0 = F * k0 * (ce ** alpha_a) * (vacancies ** alpha_a) * (cs ** alpha_c)

    # di0/d(c_s_surf) = i0 * (alpha_a * (-1/vacancies) + alpha_c / cs)
    #   = i0 * (-alpha_a/(c_s_max - cs) + alpha_c/cs)
    di0_dc_surf = i0 * (-alpha_a / vacancies + alpha_c / cs)

    # di0/d(c_e) = i0 * alpha_a / c_e
    di0_dc_e = i0 * alpha_a / ce

    return ExchangeCurrentResult(
        i0=i0,
        di0_dc_surf=di0_dc_surf,
        di0_dc_e=di0_dc_e,
    )


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
) -> ReactionRates:
    """计算 Butler-Volmer 反应电流及其全部解析导数。

    eta = phi_s - phi_e - U(c_s_surf)
    i_F = i0 * [exp(alpha_a * f * eta) - exp(-alpha_c * f * eta)]

    其中 f = F / (R * T)。

    导数:
        B = di_F/d(eta) = i0 * f * (alpha_a * exp_a + alpha_c * exp_c)
        di_F/d(c_e) = (di0/d(c_e)) * (exp_a - exp_c)
        di_F/d(c_s_surf) = (di0/d(c_s_surf)) * (exp_a - exp_c) - B * dU/d(c_s)

    Parameters
    ----------
    c_e : ndarray
        电解质浓度 [mol/m³]。
    c_s_surf : ndarray
        固相表面浓度 [mol/m³]。
    phi_e : ndarray
        电解质电位 [V]。
    phi_s : ndarray
        固相电位 [V]。
    c_s_max : float
        固相最大浓度 [mol/m³]。
    ocv : ndarray
        OCV 值 U(c_s_surf) [V]。
    docv_dc_s : ndarray
        dU/d(c_s_surf) [V·m³/mol]。
    kinetics : Kinetics
        动力学参数。
    temperature : float
        温度 [K]。

    Returns
    -------
    ReactionRates
    """
    c_e = np.asarray(c_e, dtype=float)
    c_s_surf = np.asarray(c_s_surf, dtype=float)
    phi_e = np.asarray(phi_e, dtype=float)
    phi_s = np.asarray(phi_s, dtype=float)
    ocv = np.asarray(ocv, dtype=float)
    docv_dc_s = np.asarray(docv_dc_s, dtype=float)

    alpha_a = kinetics.alpha_a
    alpha_c = kinetics.alpha_c

    # 过电位
    eta = phi_s - phi_e - ocv

    # 交换电流密度及导数
    exc = exchange_current_density_vec(c_e, c_s_surf, c_s_max, kinetics)
    i0 = exc.i0
    di0_dc_surf = exc.di0_dc_surf
    di0_dc_e = exc.di0_dc_e

    # f = F/(R*T) [1/V]
    f = F / (R * temperature)

    # 指数参数 (clip 防溢出)
    arg_a = np.clip(alpha_a * f * eta, -500.0, 500.0)
    arg_c = np.clip(-alpha_c * f * eta, -500.0, 500.0)
    exp_a = np.exp(arg_a)
    exp_c = np.exp(arg_c)

    # Faradaic 电流
    i_f = i0 * (exp_a - exp_c)

    # B = di_F/d(eta) = i0 * f * (alpha_a * exp_a + alpha_c * exp_c)
    B = i0 * f * (alpha_a * exp_a + alpha_c * exp_c)

    # di_F/d(c_e) = (di0/d(c_e)) * (exp_a - exp_c)
    #   注: U 不依赖于 c_e (empirical OCV vs Li/Li+)
    di_dc_e = di0_dc_e * (exp_a - exp_c)

    # di_F/d(c_s_surf) = (di0/d(c_s_surf)) * (exp_a - exp_c) - B * dU/d(c_s)
    di_dc_surf = di0_dc_surf * (exp_a - exp_c) - B * docv_dc_s

    return ReactionRates(
        eta=eta,
        i0=i0,
        i_f=i_f,
        di_deta=B,
        di_dc_e=di_dc_e,
        di_dc_surf=di_dc_surf,
    )
