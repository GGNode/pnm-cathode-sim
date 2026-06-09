"""Single-pore discharge simulator.

A single NMC532 spherical particle discharging at constant C-rate
in a well-mixed electrolyte (no concentration gradients in liquid phase).

This is the Phase 1 verification test case:
  - Butler-Volmer kinetics at the surface
  - Solid-phase diffusion inside the particle
  - NMC532 OCV curve
  - Constant-current discharge until cutoff voltage

The particle is treated as a single pore node with:
  - No electrolyte transport limitations (c_e = const)
  - No electronic resistance (phi_s uniform)
  - Only solid diffusion + BV kinetics
"""

import numpy as np

from src.physics.ocv import nmc532_ocv
from src.physics.reaction import butler_volmer, exchange_current_density, F, R
from src.physics.solid import (
    nmc532_diffusion_coefficient,
    discretize_spherical_particle,
    solid_diffusion_rhs,
)


class SinglePoreDischarge:
    """Simulate galvanostatic discharge of a single NMC532 spherical particle."""

    def __init__(
        self,
        R_p: float = 5e-6,
        N_shell: int = 10,
        c_s_init: float = 0.5,
        c_e: float = 1200.0,
        T: float = 298.15,
        C_rate: float = 1.0,
        k0: float = 5e-10,
        cs_max: float = 48900.0,
    ):
        """
        Parameters
        ----------
        R_p : float
            Particle radius [m].
        N_shell : int
            Number of diffusion shells.
        c_s_init : float
            Initial lithiation degree x = c_s / c_s_max [0,1].
        c_e : float
            Electrolyte Li+ concentration [mol/m³].
        T : float
            Temperature [K].
        C_rate : float
            Discharge C-rate (1C = full discharge in 1 hour).
        k0 : float
            Reaction rate constant.
        cs_max : float
            Maximum Li concentration in NMC532 [mol/m³].
        """
        self.R_p = R_p
        self.N_shell = N_shell
        self.c_e = c_e
        self.T = T
        self.C_rate = C_rate
        self.k0 = k0
        self.cs_max = cs_max

        # Discretize particle
        self.dr, self.volumes = discretize_spherical_particle(R_p, N_shell)
        self.total_volume = sum(self.volumes)

        # Initial concentration (uniform)
        self.c_s_init = c_s_init * cs_max
        self.surface_area = 4.0 * np.pi * R_p**2

        # Nominal capacity: Q = n * F * V_particle * cs_max / 3600 [Ah]
        # For a single particle, this is tiny. We express current in A/m² of surface.
        # Nominal current density for 1C: I_1C = F * cs_max * R_p / 3 [A/m²]
        # (approximate: full extraction in 1 hour)
        self.I_1C = F * cs_max * R_p / 3.0 / 3600.0  # A/m² of surface

    def run(
        self,
        cutoff_V: float = 3.0,
        dt: float = 5.0,
        t_max: float = 40000.0,
    ) -> dict:
        """Run constant-current discharge.

        Parameters
        ----------
        cutoff_V : float
            Cutoff voltage [V].
        dt : float
            Time step [s].
        t_max : float
            Maximum simulation time [s].

        Returns
        -------
        result : dict
            Keys: time, voltage, capacity, c_s_surface, c_s_bulk.
        """
        N = self.N_shell
        c_s = np.ones(N) * self.c_s_init

        # Applied current density (negative for discharge: cathodic)
        I_app = -self.C_rate * self.I_1C  # A/m² (negative = reduction = lithiation)

        times = [0.0]
        voltages = []
        capacities = []
        c_s_surfaces = []
        c_s_bulks = []

        t = 0.0
        Q_delivered = 0.0  # Ah/m² delivered

        # Compute initial voltage
        cs_surface = c_s[-1]
        x_surface = cs_surface / self.cs_max
        U_eq = nmc532_ocv(x_surface)
        # At equilibrium with no current, V = U_eq
        voltages.append(U_eq)
        c_s_surfaces.append(cs_surface)
        c_s_bulks.append(np.mean(c_s))
        capacities.append(0.0)

        while t < t_max:
            # --- Solve for surface conditions ---
            # At the surface, the BV current must equal the applied current:
            #   I_rxn = I_app (galvanostatic constraint)
            #   I_rxn = i0 * [exp(a_a*F*eta/(RT)) - exp(-a_c*F*eta/(RT))]
            #   eta = phi_s - phi_e - U_eq
            # For a single pore: phi_e = 0 (reference), phi_s = eta + U_eq
            # We solve for eta such that BV(i0, eta) = I_app

            cs_surface = c_s[-1]
            x_surface = cs_surface / self.cs_max
            U_eq = nmc532_ocv(x_surface)

            i0 = exchange_current_density(
                k0=self.k0, ce=self.c_e, cs=cs_surface, cs_max=self.cs_max,
            )

            # Solve BV for eta: i0 * [exp(a*F*eta/RT) - exp(-a*F*eta/RT)] = I_app
            # For symmetric alpha=0.5: i0 * 2*sinh(0.5*F*eta/RT) = I_app
            # eta = (RT/(0.5*F)) * arcsinh(I_app / (2*i0))
            # But use Newton for generality
            eta = self._solve_eta_bisection(i0, I_app)

            # Cell voltage: V = phi_s - phi_e(separator)
            # For single pore: phi_e = 0, phi_s = eta + U_eq
            V_cell = eta + U_eq  # phi_s

            # Check cutoff
            if V_cell < cutoff_V:
                break

            # --- Compute surface flux from BV ---
            I_rxn = butler_volmer(i0=i0, eta=eta, T=self.T)
            # In a half-cell, the NMC cathode is lithiated during discharge:
            #   NMC + Li+ + e- → Li-NMC  (reduction, I_rxn < 0)
            # Li flux INTO particle = -I_rxn / F  (positive when I_rxn < 0)
            flux_surface = -I_rxn / F

            # --- Update solid concentration ---
            # Semi-implicit: use current c_s for RHS, advance with forward Euler
            D_s = nmc532_diffusion_coefficient(np.mean(c_s), self.T)
            rhs = solid_diffusion_rhs(c_s, self.R_p, D_s, flux_surface, N)
            c_s = c_s + dt * rhs

            # Clamp to physical range
            c_s = np.clip(c_s, 1.0, self.cs_max - 1.0)

            # --- Accumulate capacity ---
            Q_delivered += abs(I_rxn) * self.surface_area * dt / 3600.0  # Ah

            t += dt
            times.append(t)
            voltages.append(V_cell)
            c_s_surfaces.append(c_s[-1])
            c_s_bulks.append(np.mean(c_s))
            capacities.append(Q_delivered)

        return {
            "time": np.array(times),
            "voltage": np.array(voltages),
            "capacity": np.array(capacities),
            "c_s_surface": np.array(c_s_surfaces),
            "c_s_bulk": np.array(c_s_bulks),
        }

    def _solve_eta_bisection(
        self, i0: float, I_target: float, tol: float = 1e-8,
    ) -> float:
        """Solve BV equation for eta using bisection.

        Find eta such that butler_volmer(i0, eta) = I_target.
        """
        # For cathodic (I_target < 0), eta < 0
        # For anodic (I_target > 0), eta > 0
        # BV is monotonically increasing in eta

        # Set search bounds
        eta_low = -0.5 if I_target < 0 else -0.01
        eta_high = 0.01 if I_target < 0 else 0.5

        # Expand bounds if needed
        for _ in range(20):
            f_low = butler_volmer(i0, eta_low, self.T) - I_target
            f_high = butler_volmer(i0, eta_high, self.T) - I_target
            if f_low * f_high < 0:
                break
            eta_low *= 2
            eta_high *= 2

        # Bisection
        for _ in range(100):
            eta_mid = 0.5 * (eta_low + eta_high)
            f_mid = butler_volmer(i0, eta_mid, self.T) - I_target
            if abs(f_mid) < tol:
                return eta_mid
            if f_mid < 0:
                eta_low = eta_mid
            else:
                eta_high = eta_mid

        return 0.5 * (eta_low + eta_high)
