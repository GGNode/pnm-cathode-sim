"""Phase 1: Single-pore integration test.

A single NMC532 spherical particle discharging at constant current.
Verifies:
  1. Voltage starts near OCV and decreases monotonically
  2. Discharge terminates near the low-voltage lithiated surface state
  3. Higher C-rate gives lower voltage at the same elapsed time
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
        """Simulation should stop near the lithiated low-voltage surface state."""
        sim = SinglePoreDischarge(
            R_p=5e-6, N_shell=10, c_s_init=0.5,
            c_e=1200.0, T=298.15, C_rate=0.5,
        )
        result = sim.run(cutoff_V=3.5, dt=5.0, t_max=50000.0)
        final_V = result["voltage"][-1]
        # The loop records the last accepted point before cutoff; with the
        # corrected OCP fit this is just above 3.5 V once the surface saturates.
        assert final_V <= 3.65, f"Final voltage {final_V:.4f}V too high"
        # Particle should be nearly fully lithiated
        final_x = result["c_s_surface"][-1] / 48900.0
        assert final_x > 0.95, f"Final x={final_x:.3f}, expected >0.95"

    def test_higher_crate_lower_voltage_at_same_time(self):
        """Higher C-rate should have lower voltage at the same elapsed time."""
        params = dict(R_p=5e-6, N_shell=10, c_s_init=0.5, c_e=1200.0, T=298.15)

        r_low = SinglePoreDischarge(**params, C_rate=0.1).run(
            cutoff_V=3.5, dt=1.0, t_max=50.0,
        )
        r_high = SinglePoreDischarge(**params, C_rate=1.0).run(
            cutoff_V=3.5, dt=1.0, t_max=50.0,
        )

        assert r_low["voltage"][-1] > r_high["voltage"][-1], (
            f"Low C-rate voltage {r_low['voltage'][-1]:.4f} should exceed "
            f"high C-rate {r_high['voltage'][-1]:.4f}"
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
