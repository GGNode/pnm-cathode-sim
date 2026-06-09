"""Phase 1: Single-pore integration test.

A single NMC532 spherical particle discharging at constant current.
Verifies:
  1. Voltage starts near OCV and decreases monotonically
  2. Discharge terminates near NMC532 minimum OCV (~3.53V)
  3. Higher C-rate gives lower capacity (diffusion limitation)
  4. At very low C-rate, particle fully lithiates
"""

import numpy as np
import pytest

from src.solver.single_pore import SinglePoreDischarge


class TestSinglePoreDischarge:
    """Integration tests for single spherical particle discharge."""

    def test_voltage_decreases(self):
        """Cell voltage should decrease monotonically during discharge."""
        sim = SinglePoreDischarge(
            R_p=5e-6, N_shell=10, c_s_init=0.5,
            c_e=1200.0, T=298.15, C_rate=0.1,
        )
        result = sim.run(cutoff_V=3.5, dt=10.0, t_max=50000.0)
        voltages = result["voltage"]
        # Check monotonic decrease (allow small numerical noise)
        for i in range(1, len(voltages)):
            assert voltages[i] <= voltages[i - 1] + 0.001, (
                f"Voltage increased: V[{i-1}]={voltages[i-1]:.4f} -> V[{i}]={voltages[i]:.4f}"
            )

    def test_discharge_reaches_near_minimum_ocv(self):
        """Simulation should stop near NMC532 minimum OCV (~3.53V).

        NMC532 OCV at full lithiation (x=1) is ~3.53V. With BV overpotential,
        the cell voltage reaches slightly below this. Use 3.5V as cutoff.
        """
        sim = SinglePoreDischarge(
            R_p=5e-6, N_shell=10, c_s_init=0.5,
            c_e=1200.0, T=298.15, C_rate=0.5,
        )
        result = sim.run(cutoff_V=3.5, dt=5.0, t_max=50000.0)
        final_V = result["voltage"][-1]
        # Should be near or below 3.5V (NMC532 OCV minimum + overpotential)
        assert final_V <= 3.58, f"Final voltage {final_V:.4f}V too high"
        # Particle should be nearly fully lithiated
        final_x = result["c_s_surface"][-1] / 48900.0
        assert final_x > 0.95, f"Final x={final_x:.3f}, expected >0.95"

    def test_higher_crate_lower_capacity(self):
        """Higher C-rate should deliver less capacity (diffusion limitation)."""
        params = dict(R_p=5e-6, N_shell=10, c_s_init=0.5, c_e=1200.0, T=298.15)

        r_low = SinglePoreDischarge(**params, C_rate=0.1).run(
            cutoff_V=3.5, dt=10.0, t_max=200000.0,
        )
        r_high = SinglePoreDischarge(**params, C_rate=1.0).run(
            cutoff_V=3.5, dt=1.0, t_max=50000.0,
        )

        # Compare total delivered charge (capacity)
        cap_low = r_low["capacity"][-1]
        cap_high = r_high["capacity"][-1]
        assert cap_low > cap_high, (
            f"Low C-rate capacity {cap_low:.2e} should exceed high C-rate {cap_high:.2e}"
        )

    def test_result_structure(self):
        """Result dict should contain expected keys."""
        sim = SinglePoreDischarge(
            R_p=5e-6, N_shell=10, c_s_init=0.5,
            c_e=1200.0, T=298.15, C_rate=0.5,
        )
        result = sim.run(cutoff_V=3.5, dt=5.0, t_max=50000.0)
        for key in ["time", "voltage", "capacity", "c_s_surface", "c_s_bulk"]:
            assert key in result, f"Missing key: {key}"
        assert len(result["time"]) == len(result["voltage"])
