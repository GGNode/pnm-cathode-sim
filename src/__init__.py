"""PNM-LIB-Cathode: 锂离子电池阴极孔隙网络放电模型。

孔隙网络模型 (PNM) 将多孔电极离散化为节点 (孔隙体) 和键 (喉道) 的图结构,
耦合电解质传输、固相扩散和 Butler-Volmer 界面反应动力学,
仿真 NMC532 阴极的恒流放电行为。

参考文献: Khan et al. (2021), J. Electrochem. Soc. 168(7), 070534.
"""

__version__ = "0.1.0"
