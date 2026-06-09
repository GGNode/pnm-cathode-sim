"""Steady-state potential solver for pore network.

Solves the coupled electrolyte (φ_e) and solid (φ_s) potential distributions
using Ohm's law on the pore network with Butler-Volmer coupling.

References:
    Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

import numpy as np
import openpnm as op
from scipy import sparse
from scipy.sparse.linalg import spsolve

from src.physics.ocv import nmc532_ocv
from src.physics.reaction import butler_volmer, exchange_current_density, F, R


class SteadyStateSolver:
    """Steady-state potential solver on a pore network."""

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
        self.n_e = len(self.e_indices)
        self.n_s = len(self.s_indices)
        self.e_map = {g: l for l, g in enumerate(self.e_indices)}
        self.s_map = {g: l for l, g in enumerate(self.s_indices)}

        self._compute_conductances()
        self._build_matrices()

    def set_concentration(self, c_e: float = 1200.0, c_s: float = 24450.0):
        self.c_e[:] = c_e
        self.c_s[:] = c_s

    def _compute_conductances(self):
        net = self.net
        A = net["throat.area"]
        L = net["throat.length"]
        D_e = 2.0e-10
        kappa = 1200.0 * F**2 * D_e / (R * self.T)
        sigma_nmc, sigma_cbd = 0.01, 760.0

        p1, p2 = self.conns[:, 0], self.conns[:, 1]
        self.G_e = np.where(self.e_mask[p1] & self.e_mask[p2], kappa * A / L, 0.0)
        s1, s2 = self.solid_mask[p1], self.solid_mask[p2]
        sig1 = np.where(self.nmc_mask[p1], sigma_nmc, np.where(self.cbd_mask[p1], sigma_cbd, 0.0))
        sig2 = np.where(self.nmc_mask[p2], sigma_nmc, np.where(self.cbd_mask[p2], sigma_cbd, 0.0))
        sig_avg = 2.0 * sig1 * sig2 / (sig1 + sig2 + 1e-30)
        self.G_s = np.where(s1 & s2, sig_avg * A / L, 0.0)

        self.interfaces = []
        for t in range(self.Nt):
            i, j = int(self.conns[t, 0]), int(self.conns[t, 1])
            if self.e_mask[i] and self.nmc_mask[j]:
                self.interfaces.append((i, j, A[t]))
            elif self.nmc_mask[i] and self.e_mask[j]:
                self.interfaces.append((j, i, A[t]))

    def _build_matrices(self):
        n_e, n_s, conns = self.n_e, self.n_s, self.conns

        # Electrolyte Laplacian
        er, ec, ev = [], [], []
        for t in range(self.Nt):
            if self.G_e[t] == 0:
                continue
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            l1, l2 = self.e_map[g1], self.e_map[g2]
            er += [l1, l2, l1, l2]; ec += [l2, l1, l1, l2]
            ev += [self.G_e[t], self.G_e[t], -self.G_e[t], -self.G_e[t]]
        self.L_e = sparse.csr_matrix((ev, (er, ec)), shape=(n_e, n_e))

        # Solid Laplacian
        sr, sc, sv = [], [], []
        for t in range(self.Nt):
            if self.G_s[t] == 0:
                continue
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            l1, l2 = self.s_map[g1], self.s_map[g2]
            sr += [l1, l2, l1, l2]; sc += [l2, l1, l1, l2]
            sv += [self.G_s[t], self.G_s[t], -self.G_s[t], -self.G_s[t]]
        self.L_s = sparse.csr_matrix((sv, (sr, sc)), shape=(n_s, n_s))

        # Interface vectors (pre-indexed to local)
        self._ie_e = np.array([self.e_map[ie[0]] for ie in self.interfaces])
        self._ie_s = np.array([self.s_map[ie[1]] for ie in self.interfaces])
        self._ie_A = np.array([ie[2] for ie in self.interfaces])
        self._ie_e_g = np.array([ie[0] for ie in self.interfaces])
        self._ie_s_g = np.array([ie[1] for ie in self.interfaces])

        # Boundary sets
        x = self.coords[:, 0]
        x_min, x_max = x.min(), x.max()
        span = x_max - x_min if x_max > x_min else 1.0
        self.sep_e = np.array([self.e_map[g] for g in self.e_indices if (x[g] - x_min) < span * 0.15])
        self.sep_s = np.array([self.s_map[g] for g in self.s_indices if (x[g] - x_min) < span * 0.15])
        self.cc_s = np.array([self.s_map[g] for g in self.s_indices if (x_max - x[g]) < span * 0.15])
        if len(self.sep_s) == 0 and n_s > 0:
            self.sep_s = np.array([int(np.argmin(x[self.s_indices] - x_min))])
        if len(self.cc_s) == 0 and n_s > 0:
            self.cc_s = np.array([int(np.argmax(x[self.s_indices]))])

        # Collector area: sum of throat cross-sections at collector boundary
        cc_g = set(self.s_indices[self.cc_s])
        self._cc_area = 0.0
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            if g1 in cc_g or g2 in cc_g:
                self._cc_area += self.net["throat.area"][t]
        if self._cc_area == 0:
            self._cc_area = 1.0

    def solve(self, I_app: float = 0.0, tol: float = 1e-8, max_iter: int = 100) -> dict:
        """Solve with successive-over-relaxation on the linearised BV system.

        Each iteration:
          1. Compute dI/dη at current potentials.
          2. Build coupled operator including BV Jacobian on diagonal.
          3. Solve the linear system.
        """
        n_e, n_s = self.n_e, self.n_s
        n_intf = len(self.interfaces)

        # Dirichlet values
        phi_e_bc = np.zeros(n_e)
        phi_s_bc = np.array([
            nmc532_ocv(self.c_s[self.s_indices[i]] / 48900.0) for i in range(n_s)
        ])

        phi_e = phi_e_bc.copy()
        phi_s = phi_s_bc.copy()

        for iteration in range(max_iter):
            # --- Compute BV quantities at interfaces ---
            if n_intf > 0:
                phi_e_intf = phi_e[self._ie_e]
                phi_s_intf = phi_s[self._ie_s]
                cs = self.c_s[self._ie_s_g]
                ce = self.c_e[self._ie_e_g]
                x_j = cs / 48900.0
                U_eq = nmc532_ocv(x_j)
                i0 = np.array([
                    exchange_current_density(self.k0, ce[k], cs[k], 48900.0)
                    for k in range(n_intf)
                ])
                eta = phi_s_intf - phi_e_intf - U_eq
                I_rxn = np.array([butler_volmer(i0[k], eta[k], self.T) for k in range(n_intf)])

                f_val = F / (R * self.T)
                arg_a = np.clip(0.5 * f_val * eta, -500, 500)
                arg_c = np.clip(-0.5 * f_val * eta, -500, 500)
                dI_deta = i0 * (0.5 * f_val * np.exp(arg_a) + 0.5 * f_val * np.exp(arg_c))
                g_bv = dI_deta * self._ie_A  # conductance-like BV coupling
            else:
                I_rxn = np.array([])
                g_bv = np.array([])

            # --- Build coupled system: M @ [phi_e; phi_s] = rhs ---
            # Electrolyte block: (L_e - diag(g_bv_e)) @ phi_e  +  off-diag
            # The BV term at interface i: I_rxn * A = g_bv * (phi_s_j - phi_e_i - U_eq)
            # In the electrolyte equation: +I_rxn * A (source)
            # → diagonal contribution: -g_bv (from d/dphi_e_i)
            # → off-diagonal to phi_s: +g_bv

            # Build g_bv_e[i] = sum of g_bv for all interfaces at e-pore i
            g_bv_e = np.zeros(n_e)
            if n_intf > 0:
                np.add.at(g_bv_e, self._ie_e, g_bv)

            # Modified electrolyte Laplacian: L_e - diag(g_bv_e)
            M_e = self.L_e.tolil()
            for i in range(n_e):
                M_e[i, i] -= g_bv_e[i]
            M_e = M_e.tocsr()

            # RHS for electrolyte: -sum(g_bv * (phi_s_j - U_eq)) at each e-pore
            rhs_e = np.zeros(n_e)
            if n_intf > 0:
                bv_source = g_bv * (phi_s[self._ie_s] - U_eq)
                np.add.at(rhs_e, self._ie_e, -bv_source)

            # Solid block: (L_s + diag(g_bv_s)) @ phi_s  -  off-diag
            # BV term at interface i: -I_rxn * A (sink in solid eq)
            # → diagonal contribution: +g_bv (from d/dphi_s_j)
            # → off-diagonal to phi_e: -g_bv

            g_bv_s = np.zeros(n_s)
            if n_intf > 0:
                np.add.at(g_bv_s, self._ie_s, g_bv)

            # Solid: (L_s - g_bv) @ phi_s = -g_bv*(phi_e + U_eq) - I_app
            M_s = self.L_s.tolil()
            for i in range(n_s):
                M_s[i, i] -= g_bv_s[i]
            M_s = M_s.tocsr()

            rhs_s = np.zeros(n_s)
            if n_intf > 0:
                bv_sink = g_bv * (phi_e[self._ie_e] + U_eq)
                np.add.at(rhs_s, self._ie_s, -bv_sink)

            # Applied current: positive I_app → discharge (current into solid)
            # I_app is current density [A/m²]; multiply by collector area
            # to get current [A] per pore.
            rhs_s[self.cc_s] -= I_app * self._cc_area / max(len(self.cc_s), 1)

            # Apply Dirichlet BCs
            M_e = M_e.tolil()
            for i in self.sep_e:
                M_e[i, :] = 0; M_e[i, i] = 1.0
                rhs_e[i] = phi_e_bc[i]
            # Pin disconnected e-pores
            for i in range(n_e):
                if M_e[i, :].nnz == 0:
                    M_e[i, i] = 1.0; rhs_e[i] = 0.0
            M_e = M_e.tocsr()

            M_s = M_s.tolil()
            for i in self.sep_s:
                M_s[i, :] = 0; M_s[i, i] = 1.0
                rhs_s[i] = phi_s_bc[i]
            for i in range(n_s):
                if M_s[i, :].nnz == 0:
                    M_s[i, i] = 1.0; rhs_s[i] = phi_s_bc[i]
            M_s = M_s.tocsr()

            # Solve
            try:
                phi_e_new = spsolve(M_e, rhs_e)
            except Exception:
                phi_e_new = phi_e.copy()
            try:
                phi_s_new = spsolve(M_s, rhs_s)
            except Exception:
                phi_s_new = phi_s.copy()

            # Under-relaxation for stability
            alpha = 0.5
            de = np.max(np.abs(phi_e_new - phi_e))
            ds = np.max(np.abs(phi_s_new - phi_s))
            phi_e = (1 - alpha) * phi_e + alpha * phi_e_new
            phi_s = (1 - alpha) * phi_s + alpha * phi_s_new

            if max(de, ds) < tol:
                break

        # Build full-network output
        phi_e_full = np.full(self.Np, np.nan)
        phi_s_full = np.full(self.Np, np.nan)
        for g in range(self.Np):
            if self.e_mask[g]:
                phi_e_full[g] = phi_e[self.e_map[g]]
            if self.solid_mask[g]:
                phi_s_full[g] = phi_s[self.s_map[g]]

        # Voltage
        if len(self.cc_s) > 0 and len(self.sep_e) > 0:
            V_cell = float(np.mean(phi_s[self.cc_s])) - float(np.mean(phi_e[self.sep_e]))
        else:
            V_cell = 0.0

        # I_rxn at all throats
        I_rxn_all = np.zeros(self.Nt)
        for t in range(self.Nt):
            g1, g2 = int(self.conns[t, 0]), int(self.conns[t, 1])
            if self.e_mask[g1] and self.nmc_mask[g2]:
                le, lnmc = self.e_map[g1], self.s_map[g2]
                cs = self.c_s[g2]
                U_eq = nmc532_ocv(cs / 48900.0)
                i0 = exchange_current_density(self.k0, self.c_e[g1], cs, 48900.0)
                eta = phi_s[lnmc] - phi_e[le] - U_eq
                I_rxn_all[t] = butler_volmer(i0, eta, self.T)
            elif self.nmc_mask[g1] and self.e_mask[g2]:
                le, lnmc = self.e_map[g2], self.s_map[g1]
                cs = self.c_s[g1]
                U_eq = nmc532_ocv(cs / 48900.0)
                i0 = exchange_current_density(self.k0, self.c_e[g2], cs, 48900.0)
                eta = phi_s[lnmc] - phi_e[le] - U_eq
                I_rxn_all[t] = butler_volmer(i0, eta, self.T)

        return {
            "phi_e": phi_e_full,
            "phi_s": phi_s_full,
            "voltage": V_cell,
            "I_rxn": I_rxn_all,
        }
