"""PNM-LIB-Cathode 共享测试夹具。"""

import numpy as np
import pytest


# --- 物理常数 ---
F = 96485.3329  # 法拉第常数 [C/mol]
R = 8.314462    # 气体常数 [J/(mol·K)]
T = 298.15      # 温度 [K]


@pytest.fixture
def constants():
    """以字典形式返回物理常数。"""
    return {"F": F, "R": R, "T": T}
