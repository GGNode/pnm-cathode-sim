# 数据

## 论文参考

- `khan2021_pnm_lib_cathode.pdf` — Khan et al. (2021) 原始论文
- `paper.pdf` — 论文副本

## 论文图表

`paper_figures/` — 从论文中裁剪的图表 PNG, 用于验证对比:

- `figure4_crop_200dpi.png` — 0.2C、0.5C、1C、3C 的 V-Q 放电曲线
- `figure6_crop_200dpi.png` — 1C 和 3C (75% SoL) 的空间锂化度分布
- `figure7_crop_200dpi.png` — 电解质浓度分布
- `page_08_200dpi.png` — `page_11_200dpi.png` — 完整页面渲染

## 验证结果

`validation/` — 预计算的验证输出:

- `discharge_*.npz` — 放电仿真结果 (0.2C、0.5C、1C、3C)
- `spatial_*_75sol.npz` — 75% 锂化度时的空间快照
- `figure*_comparison.png` — 生成的对比图
- `validation_metrics.json` — 定量验证指标
- `parallel_results.json` — 并行多场景验证结果

## 备注

- 如果原始 Khan et al. (2021) NREL XCT 数据可用, 可通过 `scripts/generate_network.py` 生成几何匹配的网络。
- 当前验证使用合成立方网络; 与论文图表的定量对比仅为定性验证。
