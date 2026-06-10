# 数据目录

## 文件说明

| 文件 | 说明 |
|------|------|
| `khan2021_pnm_lib_cathode.pdf` | Khan et al. 2021 原始论文 PDF |
| `paper_figures/` | 从 PDF 提取的论文图表（用于对比验证） |
| `validation/` | 验证结果数据（.npz, .json, .png） |

## 验证结果

| 文件 | 说明 |
|------|------|
| `parallel_results.json` | 并行验证结果（16 场景：4 C-rate × 2 电流基准 × 2 隔膜状态） |
| `discharge_*.npz` | 各倍率放电数据（时间、电压、容量） |
| `spatial_*.npz` | 空间分布数据（浓度、电位、SoL） |
| `figure*.png` | 论文 Figure 4/6/7 对比图 |
| `validation_metrics.json` | 验证指标汇总 |

## 数据格式

### .npz 文件

使用 `numpy.load()` 加载：

```python
import numpy as np
data = np.load("data/validation/discharge_1c.npz")
print(data.files)  # ['time', 'voltage', 'capacity', ...]
```

### 论文参考数据

详见 `PAPER_REFERENCE.md`（根目录），包含：
- Table I: XCT 图像属性
- Table II: 模型参数
- Figure 4/6/7: 定量数据提取
