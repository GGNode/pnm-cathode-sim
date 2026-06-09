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
        if C_rate < 0:
            raise ValueError("C_rate must be non-negative; use positive values for discharge")
        capacity_c = self.discharge_capacity_coulombs(cs_max=cs_max)
        current_a = C_rate * capacity_c / 3600.0
        return -current_a / self.collector_area

    def characteristic_dt(
        self,
        I_app: float = 0.0,
        C_rate: float | None = None,
        safety: float = 0.05,
        dt_min: float = 1e-6,
        dt_max: float = 300.0,
    ) -> float:
        """Estimate an initial time step from diffusion and reaction timescales."""
        x = self.coords[:, 0]
        length = max(float(x.max() - x.min()), 1e-12)
        D_e = 2.0e-10 * 0.39**1.5
        tau_e = length**2 / D_e

        if self.n_nmc:
            ds_vals = np.array([
                nmc532_diffusion_coefficient(self.c_s[g], self.T)
                for g in self.nmc_indices
            ])
            D_s = max(float(np.nanmedian(ds_vals)), 1e-30)
            center_dist = np.linalg.norm(
                self.coords[self.conns[:, 1]] - self.coords[self.conns[:, 0]],
                axis=1,
            )
            spacing = max(float(np.min(center_dist)), 1e-12)
            tau_s = spacing**2 / D_s
        else:
            tau_s = np.inf

        if C_rate is not None and C_rate > 0:
            tau_rxn = 3600.0 / C_rate
        elif abs(I_app) > 0:
            capacity_c = self.discharge_capacity_coulombs()
            current_a = abs(I_app) * self.collector_area
            tau_rxn = capacity_c / current_a if current_a > 0 else np.inf
        else:
            tau_rxn = np.inf

        times = [tau for tau in (tau_e, tau_s, tau_rxn) if np.isfinite(tau) and tau > 0]
        dt = safety * min(times) if times else 1.0
        return float(np.clip(dt, dt_min, dt_max))

    def _snapshot_state(self) -> tuple[np.ndarray, np.ndarray]:
        return self.c_e.copy(), self.c_s.copy()

    def _restore_state(self, snapshot: tuple[np.ndarray, np.ndarray]) -> None:
        self.c_e[:] = snapshot[0]
        self.c_s[:] = snapshot[1]

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

        # Interface data: (throat, e_local, nmc_local, A_intf, e_global, nmc_global)
        self._interfaces = []
        for t in range(self.Nt):
            g1, g2 = int(conns[t, 0]), int(conns[t, 1])
            if self.e_mask[g1] and self.nmc_mask[g2]:
                self._interfaces.append((
                    t, self.e_map[g1], self.nmc_map[g2], A[t], g1, g2,
                ))
            elif self.nmc_mask[g1] and self.e_mask[g2]:
                self._interfaces.append((
                    t, self.e_map[g2], self.nmc_map[g1], A[t], g2, g1,
                ))

    def _filter_reactive_interfaces(self):
        reactive_pairs = {
            (e_g, nmc_g) for e_g, nmc_g, _area in self._steady.reactive_interfaces
            if nmc_g in self.nmc_map
        }
        self._interfaces = [
            interface for interface in self._interfaces
            if (interface[4], interface[5]) in reactive_pairs
        ]

    def step(self, dt: float, I_app: float = 0.0, C_rate: float | None = None) -> dict:
        """Advance one time step.

        Parameters
        ----------
        dt : float
            Time step [s].
        I_app : float
            Applied current density [A/m²].
        C_rate : float or None
            Optional positive discharge C-rate. If provided, it is converted to
            anodic-convention current density from the current electrode state.

        Returns
        -------
        result : dict
            Keys: phi_e, phi_s, c_e, c_s, voltage, I_rxn.
        """
        if C_rate is not None:
            I_app = self.current_density_for_c_rate(C_rate)

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
        eta_all = np.zeros(self.Nt)
        reaction_mol_rate = 0.0

        for t, le, ls, A_intf, e_g, nmc_g in self._interfaces:
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
            eta_all[t] = eta

            # Reaction source/sink, anodic convention:
            # I_rxn > 0: Li leaves solid and enters electrolyte.
            # I_rxn < 0: discharge/lithiation consumes electrolyte Li.
            flux = I_rxn * A_intf / F  # mol/s
            reaction_mol_rate += flux
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
            "eta": eta_all,
            "reaction_mol_rate": reaction_mol_rate,
            "I_app": I_app,
            "dt": dt,
            "iterations": pot["iterations"],
            "converged": pot["converged"],
        }

    def adaptive_step(
        self,
        dt: float,
        I_app: float,
        previous_voltage: float | None = None,
        voltage_jump_limit: float = 0.05,
        growth_factor: float = 1.25,
        shrink_factor: float = 0.5,
        fast_iterations: int = 50,
        dt_min: float = 1e-6,
        dt_max: float = 300.0,
        max_retries: int = 24,
    ) -> tuple[dict, float, float]:
        """Advance with adaptive time-step rejection and growth.

        Returns
        -------
        result, dt_used, dt_next
        """
        trial_dt = float(np.clip(dt, dt_min, dt_max))
        last_result = None

        for _ in range(max_retries):
            snapshot = self._snapshot_state()
            result = self.step(trial_dt, I_app=I_app)
            last_result = result

            voltage = result["voltage"]
            jump = 0.0 if previous_voltage is None else abs(voltage - previous_voltage)
            reject = (
                (not result["converged"])
                or (not np.isfinite(voltage))
                or jump > voltage_jump_limit
            )

            if reject and trial_dt > dt_min:
                self._restore_state(snapshot)
                trial_dt = max(trial_dt * shrink_factor, dt_min)
                continue

            dt_next = trial_dt
            if result["converged"] and result["iterations"] <= fast_iterations:
                dt_next = min(trial_dt * growth_factor, dt_max)
            elif jump > 0.5 * voltage_jump_limit:
                dt_next = max(trial_dt * shrink_factor, dt_min)

            result["dt"] = trial_dt
            result["dt_next"] = dt_next
            result["voltage_jump"] = jump
            return result, trial_dt, dt_next

        voltage = None if last_result is None else last_result.get("voltage")
        raise RuntimeError(
            f"adaptive step failed to converge after {max_retries} retries; "
            f"last voltage={voltage}"
        )

    def run_discharge(
        self,
        C_rate: float,
        cutoff_voltage: float = 2.5,
        dt: float | None = None,
        dt_min: float = 1e-6,
        dt_max: float = 300.0,
        max_steps: int = 200,
        voltage_jump_limit: float = 0.05,
    ) -> dict:
        """Run a positive-C-rate discharge with adaptive stepping."""
        if C_rate <= 0:
            raise ValueError("C_rate must be positive for discharge")

        I_app = self.current_density_for_c_rate(C_rate)
        step_dt = dt if dt is not None else self.characteristic_dt(
            I_app=I_app, C_rate=C_rate, dt_min=dt_min, dt_max=dt_max,
        )

        self._steady.set_concentration(self.c_e, self.c_s)
        initial = self._steady.solve(I_app=0.0)
        previous_voltage = initial["voltage"]
        times = [0.0]
        capacities = [0.0]
        voltages = [previous_voltage]
        dts = [0.0]
        iterations = [initial["iterations"]]

        done = previous_voltage <= cutoff_voltage
        for _ in range(max_steps):
            if done:
                break

            result, dt_used, step_dt = self.adaptive_step(
                step_dt,
                I_app=I_app,
                previous_voltage=previous_voltage,
                voltage_jump_limit=voltage_jump_limit,
                dt_min=dt_min,
                dt_max=dt_max,
            )

            t_new = times[-1] + dt_used
            q_new = abs(I_app) * t_new
            v_new = result["voltage"]

            if previous_voltage > cutoff_voltage >= v_new:
                frac = (previous_voltage - cutoff_voltage) / max(previous_voltage - v_new, 1e-30)
                times.append(times[-1] + frac * dt_used)
                capacities.append(capacities[-1] + frac * (q_new - capacities[-1]))
                voltages.append(cutoff_voltage)
                dts.append(frac * dt_used)
                iterations.append(result["iterations"])
                done = True
                break

            times.append(t_new)
            capacities.append(q_new)
            voltages.append(v_new)
            dts.append(dt_used)
            iterations.append(result["iterations"])
            previous_voltage = v_new

        return {
            "time": np.asarray(times),
            "capacity": np.asarray(capacities),
            "capacity_Ah_m2": np.asarray(capacities) / 3600.0,
            "voltage": np.asarray(voltages),
            "dt": np.asarray(dts),
            "iterations": np.asarray(iterations),
            "I_app": I_app,
            "C_rate": C_rate,
            "cutoff_voltage": cutoff_voltage,
            "reached_cutoff": done,
        }
