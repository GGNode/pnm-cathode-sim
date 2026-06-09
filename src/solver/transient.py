"""Transient discharge solver for pore network.

Combines the steady-state potential solver with time-dependent
concentration evolution for electrolyte (c_e) and solid (c_s) phases.

Uses backward Euler for concentration time stepping and the steady-state
solver for potential distribution at each time step.

References:
    Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op
from scipy import sparse
from scipy.sparse.linalg import spsolve

from src.physics.ocv import nmc532_ocv
from src.physics.reaction import butler_volmer, exchange_current_density, F, R
from src.physics.solid import nmc532_diffusion_coefficient
from src.solver.steady import SteadyStateSolver


class TransientSolver:
    """Transient discharge solver coupling potentials and concentrations."""

    def __init__(self, net: op.network.Cubic, T: float = 298.15, k0: float = 5e-10):
        self.net = net
        self.T = T
        self.k0 = k0
        self.Np = net.Np
        self.Nt = net.Nt
        self.conns = net["throat.conns"]
        self.coords = net["pore.coords"]

        self.c_e = np.ones(self.Np) * 1200.0
        self.c_s = np.ones(self.Np) * 24450.0

        self.e_mask = net["pore.electrolyte"]
        self.nmc_mask = net["pore.nmc"]
        self.cbd_mask = net["pore.cbd"]
        self.solid_mask = self.nmc_mask | self.cbd_mask

        self.e_indices = np.where(self.e_mask)[0]
        self.s_indices = np.where(self.solid_mask)[0]
        self.nmc_indices = np.where(self.nmc_mask)[0]
        self.n_e = len(self.e_indices)
        self.n_s = len(self.s_indices)
        self.n_nmc = len(self.nmc_indices)
        self.e_map = {g: l for l, g in enumerate(self.e_indices)}
        self.s_map = {g: l for l, g in enumerate(self.s_indices)}
        self.nmc_map = {g: l for l, g in enumerate(self.nmc_indices)}

        self._build_diffusion_matrices()

        # Steady-state solver for potentials
        self._steady = SteadyStateSolver(net, T=T, k0=k0)
        self._filter_reactive_interfaces()

    def set_concentration(self, c_e: float = 1200.0, c_s: float = 24450.0):
        self.c_e[:] = c_e
        self.c_s[:] = c_s
        self._steady.set_concentration(c_e, c_s)

    @property
    def collector_area(self) -> float:
        """Active current-collector area used to convert current density to current."""
        return self._steady._cc_area

    def discharge_capacity_coulombs(self, cs_max: float = 48900.0) -> float:
        """Remaining cathode lithiation capacity from current state [C]."""
        vacancies = np.maximum(cs_max - self.c_s[self.nmc_indices], 0.0)
        return float(F * np.sum(vacancies * self._vol_s))

    def current_density_for_c_rate(self, C_rate: float, cs_max: float = 48900.0) -> float:
        """Convert positive discharge C-rate to anodic-convention current density."""
        capacity_c = self.discharge_capacity_coulombs(cs_max=cs_max)
        current_a = C_rate * capacity_c / 3600.0
        return -current_a / self.collector_area

    def _build_diffusion_matrices(self):
        """Build diffusion Laplacian matrices for c_e and c_s."""
        conns = self.conns
        D_e = 2.0e-10  # electrolyte diffusion [m²/s]

        # --- Electrolyte diffusion Laplacian (n_e × n_e) ---
        # D_eff * A / L for each e-e throat
        A = self.net["throat.area"]
        L = self.net["throat.length"]
        De_eff = D_e * 0.39**1.5  # Bruggeman with porosity=0.39

        er, ec, ev = [], [], []
        self._ee_throats = []  # (local_i, local_j, conductance) for e-e throats
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            if self.e_mask[g1] and self.e_mask[g2]:
                l1, l2 = self.e_map[g1], self.e_map[g2]
                g_val = De_eff * A[t] / L[t]
                er += [l1, l2, l1, l2]; ec += [l2, l1, l1, l2]
                ev += [g_val, g_val, -g_val, -g_val]
                self._ee_throats.append((l1, l2, g_val))
        self.L_ce = sparse.csr_matrix((ev, (er, ec)), shape=(self.n_e, self.n_e))

        # --- Solid diffusion Laplacian (n_nmc × n_nmc) ---
        sr, sc, sv = [], [], []
        self._ss_throats = []
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            if self.nmc_mask[g1] and self.nmc_mask[g2]:
                l1, l2 = self.nmc_map[g1], self.nmc_map[g2]
                D_s1 = nmc532_diffusion_coefficient(self.c_s[g1], self.T)
                D_s2 = nmc532_diffusion_coefficient(self.c_s[g2], self.T)
                D_s_avg = 2.0 * D_s1 * D_s2 / (D_s1 + D_s2 + 1e-30)
                g_val = D_s_avg * A[t] / L[t]
                sr += [l1, l2, l1, l2]; sc += [l2, l1, l1, l2]
                sv += [g_val, g_val, -g_val, -g_val]
                self._ss_throats.append((l1, l2, g_val))
        self.L_cs = sparse.csr_matrix((sv, (sr, sc)), shape=(self.n_nmc, self.n_nmc))

        # Pore volumes
        self._vol_e = self.net["pore.volume"][self.e_indices]
        self._vol_s = self.net["pore.volume"][self.nmc_indices]

        # Interface data: (e_local, nmc_local, A_intf, e_global, nmc_global)
        self._interfaces = []
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            if self.e_mask[g1] and self.nmc_mask[g2]:
                self._interfaces.append((
                    self.e_map[g1], self.nmc_map[g2], A[t], g1, g2,
                ))
            elif self.nmc_mask[g1] and self.e_mask[g2]:
                self._interfaces.append((
                    self.e_map[g2], self.nmc_map[g1], A[t], g2, g1,
                ))

    def _filter_reactive_interfaces(self):
        reactive_pairs = {
            (e_g, nmc_g) for e_g, nmc_g, _area in self._steady.reactive_interfaces
            if nmc_g in self.nmc_map
        }
        self._interfaces = [
            interface for interface in self._interfaces
            if (interface[3], interface[4]) in reactive_pairs
        ]

    def step(self, dt: float, I_app: float = 0.0) -> dict:
        """Advance one time step.

        Parameters
        ----------
        dt : float
            Time step [s].
        I_app : float
            Applied current density [A/m²].

        Returns
        -------
        result : dict
            Keys: phi_e, phi_s, c_e, c_s, voltage, I_rxn.
        """
        # 1. Solve potentials at current concentrations
        self._steady.set_concentration(self.c_e, self.c_s)
        pot = self._steady.solve(I_app=I_app)

        phi_e_local = np.array([pot["phi_e"][g] for g in self.e_indices])
        phi_s_local = np.array([pot["phi_s"][g] for g in self.nmc_indices])

        # 2. Compute reaction rates at interfaces
        # For each e/NMC interface: I_rxn via Butler-Volmer
        # c_e change: +I_rxn * A_intf / F  (anodic production)
        # c_s change: -I_rxn * A_intf / F  (anodic deintercalation)
        dce_dt = np.zeros(self.n_e)
        dcs_dt = np.zeros(self.n_nmc)

        for le, ls, A_intf, e_g, nmc_g in self._interfaces:
            cs = self.c_s[nmc_g]
            ce = self.c_e[e_g]
            U_eq = nmc532_ocv(cs / 48900.0)
            i0 = exchange_current_density(self.k0, ce, cs, 48900.0)

            # Get potentials (may be NaN if solver failed)
            pe = phi_e_local[le] if le < len(phi_e_local) else 0.0
            ps = phi_s_local[ls] if ls < len(phi_s_local) else U_eq
            if np.isnan(pe):
                pe = 0.0
            if np.isnan(ps):
                ps = U_eq

            eta = ps - pe - U_eq
            I_rxn = butler_volmer(i0, eta, self.T)

            # Reaction source/sink, anodic convention:
            # I_rxn > 0: Li leaves solid and enters electrolyte.
            # I_rxn < 0: discharge/lithiation consumes electrolyte Li.
            flux = I_rxn * A_intf / F  # mol/s
            dce_dt[le] += flux / self._vol_e[le] if self._vol_e[le] > 0 else 0.0
            dcs_dt[ls] -= flux / self._vol_s[ls] if self._vol_s[ls] > 0 else 0.0

        # 3. Diffusion state vectors
        c_e_local = self.c_e[self.e_indices]
        c_s_local = self.c_s[self.nmc_indices]

        # 4. Time stepping: backward Euler for diffusion, explicit reaction.
        # V (c_new - c_old) / dt = L c_new + V source
        M_e = sparse.diags(self._vol_e / dt) - self.L_ce
        rhs_e = (self._vol_e / dt) * c_e_local + self._vol_e * dce_dt
        M_s = sparse.diags(self._vol_s / dt) - self.L_cs
        rhs_s = (self._vol_s / dt) * c_s_local + self._vol_s * dcs_dt
        c_e_new = spsolve(M_e.tocsr(), rhs_e) if self.n_e else c_e_local
        c_s_new = spsolve(M_s.tocsr(), rhs_s) if self.n_nmc else c_s_local

        # Clamp to physical range
        c_e_new = np.clip(c_e_new, 1.0, 10000.0)
        c_s_new = np.clip(c_s_new, 1.0, 48899.0)

        # Update global arrays
        self.c_e[self.e_indices] = c_e_new
        self.c_s[self.nmc_indices] = c_s_new

        return {
            "phi_e": pot["phi_e"],
            "phi_s": pot["phi_s"],
            "c_e": self.c_e.copy(),
            "c_s": self.c_s.copy(),
            "voltage": pot["voltage"],
            "I_rxn": pot["I_rxn"],
            "iterations": pot["iterations"],
            "converged": pot["converged"],
        }
