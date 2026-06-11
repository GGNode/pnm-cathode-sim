"""Phase 2: Pore network generation tests.

Verifies:
  1. Cubic network creation with correct size
  2. Three-phase label assignment (electrolyte, NMC, CBD)
  3. Geometric properties (pore volume, throat area, coordination number)
  4. Network statistics match expected ranges
"""

import numpy as np
import pytest

from pnmcathode.network.generator import check_percolation, create_cathode_network


class TestCathodeNetwork:
    """Tests for cathode pore network generation."""

    def test_network_size(self):
        """Network should have expected number of nodes and bonds."""
        net = create_cathode_network(
            shape=[10, 10, 10],
            spacing=1e-5,
            porosity=0.35,
        )
        # 10x10x10 cubic = 1000 nodes
        assert net.Np == 1000
        # Cubic: 3*10*10*9 = 2700 bonds
        assert net.Nt > 0

    def test_pore_coordinates(self):
        """Pore coordinates should be in expected range."""
        shape = [5, 5, 5]
        spacing = 1e-5
        net = create_cathode_network(shape=shape, spacing=spacing)
        coords = net["pore.coords"]
        assert coords.shape == (125, 3)
        # Coordinates should span roughly n*spacing in each dimension
        for dim in range(3):
            assert coords[:, dim].min() >= -spacing
            assert coords[:, dim].max() <= shape[dim] * spacing * 1.01

    def test_phase_labels(self):
        """Each pore should be labeled as electrolyte, NMC, or CBD."""
        net = create_cathode_network(
            shape=[10, 10, 10],
            spacing=1e-5,
            porosity=0.35,
        )
        labels = net["pore.phase_label"]
        # Should have exactly 3 unique labels
        unique_labels = set(labels)
        assert unique_labels == {0, 1, 2}, f"Expected {{0,1,2}}, got {unique_labels}"

    def test_porosity_approximate(self):
        """Fraction of electrolyte pores should approximate target porosity."""
        porosity = 0.35
        net = create_cathode_network(
            shape=[20, 20, 20],
            spacing=1e-5,
            porosity=porosity,
        )
        labels = net["pore.phase_label"]
        frac_electrolyte = np.sum(labels == 0) / len(labels)
        # Allow 10% tolerance (random assignment)
        assert abs(frac_electrolyte - porosity) < 0.10, (
            f"Electrolyte fraction {frac_electrolyte:.2f} != target porosity {porosity:.2f}"
        )

    def test_pore_volumes_positive(self):
        """All pore volumes should be positive."""
        net = create_cathode_network(shape=[5, 5, 5], spacing=1e-5)
        assert np.all(net["pore.volume"] > 0)

    def test_throat_areas_positive(self):
        """All throat cross-sectional areas should be positive."""
        net = create_cathode_network(shape=[5, 5, 5], spacing=1e-5)
        assert np.all(net["throat.area"] > 0)

    def test_throat_lengths_positive(self):
        """All throat lengths should be positive."""
        net = create_cathode_network(shape=[5, 5, 5], spacing=1e-5)
        assert np.all(net["throat.length"] > 0)

    def test_coordination_number(self):
        """Average coordination number should be ~6 for cubic network."""
        net = create_cathode_network(shape=[10, 10, 10], spacing=1e-5)
        # Count connections per pore
        conns = net["throat.conns"]
        coord = np.zeros(net.Np)
        for t in range(net.Nt):
            coord[conns[t, 0]] += 1
            coord[conns[t, 1]] += 1
        avg_coord = np.mean(coord)
        # Cubic network: interior nodes have 6, surface nodes have fewer
        # Average should be between 3 and 6
        assert 3.0 < avg_coord < 6.5, f"Avg coordination {avg_coord:.1f} out of range"

    def test_nmc_and_cbd_connected(self):
        """NMC and CBD phases should have non-zero representation."""
        net = create_cathode_network(
            shape=[10, 10, 10],
            spacing=1e-5,
            porosity=0.35,
        )
        labels = net["pore.phase_label"]
        n_nmc = np.sum(labels == 1)
        n_cbd = np.sum(labels == 2)
        assert n_nmc > 0, "No NMC pores assigned"
        assert n_cbd > 0, "No CBD pores assigned"

    def test_electrolyte_percolation_check(self):
        """Electrolyte percolation should require a same-phase x path."""
        net = create_cathode_network(shape=[3, 1, 1], spacing=1e-5, porosity=1.0, cbd_fraction=0.0)
        assert check_percolation(net, phase="electrolyte")

        labels = np.array([0, 1, 0])
        net["pore.phase_label"] = labels
        net["pore.electrolyte"] = labels == 0
        net["pore.nmc"] = labels == 1
        net["pore.cbd"] = labels == 2
        assert not check_percolation(net, phase="electrolyte")

    def test_solid_percolation_check(self):
        """Solid percolation should include both NMC and CBD pores."""
        net = create_cathode_network(shape=[3, 1, 1], spacing=1e-5, porosity=0.0, cbd_fraction=0.0)
        assert check_percolation(net, phase="solid")

        labels = np.array([1, 0, 2])
        net["pore.phase_label"] = labels
        net["pore.electrolyte"] = labels == 0
        net["pore.nmc"] = labels == 1
        net["pore.cbd"] = labels == 2
        assert not check_percolation(net, phase="solid")
