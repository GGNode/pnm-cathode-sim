"""Shared test fixtures for PNM-LIB-Cathode."""

import numpy as np
import pytest


# --- Physical Constants ---
F = 96485.3329  # Faraday constant [C/mol]
R = 8.314462    # Gas constant [J/(mol·K)]
T = 298.15      # Temperature [K]


@pytest.fixture
def constants():
    """Return physical constants as a dict."""
    return {"F": F, "R": R, "T": T}
