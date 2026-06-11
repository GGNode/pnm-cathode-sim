"""Phase 0: Verify project setup, dependencies, and OpenPNM basics."""

import numpy as np
import pytest


class TestImports:
    """Verify all required packages are importable."""

    def test_numpy(self):
        import numpy
        assert numpy.__version__

    def test_scipy(self):
        import scipy
        assert scipy.__version__

    def test_matplotlib(self):
        import matplotlib
        assert matplotlib.__version__

    def test_openpnm(self):
        import openpnm
        assert openpnm.__version__

    @pytest.mark.skip(reason="PoreSpy has scikit-image compat issue on Python 3.14")
    def test_porespy(self):
        import porespy
        assert porespy.__version__

    def test_yaml(self):
        import yaml
        assert yaml.__version__

    def test_pnmcathode_import(self):
        import pnmcathode
        assert pnmcathode.__version__


class TestOpenPNMBasics:
    """Verify OpenPNM can create and manipulate pore networks."""

    def test_create_cubic_network(self):
        """OpenPNM should create a simple cubic pore network."""
        import openpnm as op

        net = op.network.Cubic(shape=[5, 5, 5], spacing=1e-5)
        # 5x5x5 cubic = 125 nodes, 300 throats
        assert net.Np == 125
        assert net.Nt == 300

    def test_network_geometry(self):
        """Network nodes should have coordinates and throat connectivity."""
        import openpnm as op

        net = op.network.Cubic(shape=[3, 3, 3], spacing=1e-5)
        coords = net["pore.coords"]
        assert coords.shape == (27, 3)
        conns = net["throat.conns"]
        assert conns.shape[1] == 2

    def test_assign_custom_properties(self):
        """Should be able to assign custom pore/throat properties."""
        import openpnm as op

        net = op.network.Cubic(shape=[3, 3, 3], spacing=1e-5)
        # Assign pore volumes
        vol = np.ones(net.Np) * 1e-15  # 1 µm³
        net["pore.volume"] = vol
        np.testing.assert_array_equal(net["pore.volume"], vol)

    def test_add_geometry_models(self):
        """Should be able to add geometry models to network."""
        import openpnm as op

        net = op.network.Cubic(shape=[3, 3, 3], spacing=1e-5)
        geo = op.models.geometry

        net.add_model(
            propname="pore.diameter",
            model=geo.pore_size.random,
            seed=42,
        )
        net.add_model(
            propname="throat.diameter",
            model=geo.throat_size.random,
            seed=42,
        )
        net.add_model(
            propname="pore.volume",
            model=geo.pore_volume.sphere,
            pore_diameter="pore.diameter",
        )
        assert "pore.diameter" in net.keys()
        assert "throat.diameter" in net.keys()
        assert "pore.volume" in net.keys()
        assert net["pore.diameter"].shape == (net.Np,)
        assert net["throat.diameter"].shape == (net.Nt,)
        assert net["pore.volume"].shape == (net.Np,)


class TestPhysicalConstants:
    """Verify constants are correctly defined."""

    def test_faraday(self, constants):
        assert abs(constants["F"] - 96485.3329) < 0.01

    def test_gas_constant(self, constants):
        assert abs(constants["R"] - 8.314462) < 1e-5

    def test_temperature(self, constants):
        assert abs(constants["T"] - 298.15) < 0.01
