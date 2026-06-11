# Toolkit Execution Status

本文档记录 `pnmcathode` 工具包当前状态与维护验证流程。旧的迁移设计与审查材料已归档到 `docs/archive/`。

## 当前状态

- pip 包边界仅包含 `pnmcathode` 及其子包。
- 旧的 `src.network`、`src.physics`、`src.solver`、`src.post` 不再进入分发包。
- 测试已迁移到 `pnmcathode.*` import 路径。
- 子包 `network`、`physics`、`solver`、`post`、`materials` 均提供公共导出。
- 命令行入口为 `pnmcathode-discharge = "pnmcathode.cli:main"`。
- 生成数据和验证图像已加入 `.gitignore`，保留应跟踪的论文 PDF 与 `data/README.md`。

## 维护约束

- 代码标识符保持英文。
- 文档、注释和 docstring 使用中文；技术术语可保留英文。
- 数值核心修改后必须运行对应测试。
- 包边界修改后必须运行 packaging smoke test。
- 不要重新引入 `src.*` 测试 import。

## 常用验证命令

```bash
python -m pytest tests/test_packaging.py -q
python -m pytest tests/test_subpackage_imports.py -q
python -m pytest tests/test_toolkit_api.py -q
python -m pytest tests -v
```

检查旧 import：

```bash
rg 'from src\.|import src\.' tests
```

检查包发现范围：

```bash
python -m pip install -e .
python - <<'PY'
import pnmcathode
import pnmcathode.network
import pnmcathode.physics
import pnmcathode.solver
import pnmcathode.post
print(pnmcathode.__version__)
PY
```

## CLI 示例

```bash
pnmcathode-discharge --c-rate 0.2 --max-steps 2 --dt 0.01 --output data/example_discharge.npz
```

该命令会运行一次短放电仿真并写出压缩 NPZ 结果。
