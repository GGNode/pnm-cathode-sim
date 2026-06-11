"""子包公共导出测试。"""


def test_network_exports():
    """network 子包应导出网络构建与渗流检查函数。"""
    from pnmcathode.network import check_percolation, create_cathode_network

    assert callable(create_cathode_network)
    assert callable(check_percolation)


def test_solver_exports():
    """solver 子包应导出主要求解器类。"""
    from pnmcathode.solver import SinglePoreSolver, SteadySolver, TransientSolver

    assert TransientSolver.__name__ == "TransientSolver"
    assert SteadySolver.__name__ == "SteadyStateSolver"
    assert SinglePoreSolver.__name__ == "SinglePoreDischarge"


def test_physics_exports():
    """physics 子包应导出常用物理函数和隔膜数据类。"""
    from pnmcathode.physics import (
        F,
        SeparatorParams,
        butler_volmer,
        electrolyte_diffusion_coefficient,
        nmc532_diffusion_coefficient,
        nmc532_ocv,
        separator_boundary,
    )

    assert F > 0.0
    assert callable(butler_volmer)
    assert callable(electrolyte_diffusion_coefficient)
    assert callable(nmc532_diffusion_coefficient)
    assert callable(nmc532_ocv)
    assert callable(separator_boundary)
    assert SeparatorParams(enabled=False).enabled is False


def test_post_and_material_exports():
    """post 和 materials 子包应导出分析函数与材料预设。"""
    from pnmcathode.materials import nmc532_khan2021
    from pnmcathode.post import (
        coordination_number,
        plot_discharge_curve,
        pore_size_distribution,
        state_of_lithiation,
    )

    assert callable(nmc532_khan2021)
    assert callable(coordination_number)
    assert callable(plot_discharge_curve)
    assert callable(pore_size_distribution)
    assert callable(state_of_lithiation)
