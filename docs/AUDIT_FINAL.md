# 最终审计报告 — Khan 2021 PNM 阴极复现

## 验证

- 请求的命令 `python -m pytest tests/ -v` 在收集阶段因默认解释器缺少 `numpy` 而失败。
- 使用项目虚拟环境的等效运行通过: `MPLCONFIGDIR="$PWD/data/.matplotlib" .venv/bin/python -m pytest tests/ -v`
  - 结果: `88 passed, 1 skipped`。

## A. 可能导致错误物理结果的剩余代码缺陷

- 未发现会导致报告的隔膜关闭放电结果出错的核心稳态/瞬态求解器缺陷。
- Butler-Volmer 符号约定内部一致:
  - `eta = phi_s - phi_e - U_eq`
  - 放电/嵌锂时 `I_rxn < 0`
  - 放电过程中电解质源项减小, 固相源项增大。
- 耦合 Newton 残差/Jacobian 符号与阳极电流约定一致。
- 浓度依赖的 `D_e(c_e)`、`kappa(c_e)` 和 `D_s(c_s)` 已接入相关矩阵。
- 在 `scripts/parallel_validation.py` 中发现了真正的问题 (非核心求解器):
  - 隔膜开启场景在构建后修改 `solver._steady.separator`, 而非将 `separator=SeparatorParams(enabled=True)` 传入 `TransientSolver`。
  - 这使得稳态求解中启用了隔膜电压边界条件, 但 `TransientSolver.separator.enabled` 仍为 `False`, 因此瞬态隔膜浓度边界条件未被施加。
  - 对当前结果的影响较小, 因为隔膜浓度/欧姆损失仅几毫伏, 但隔膜开启验证数据并非完全自洽。
- `scripts/parallel_validation.py` 中的小问题:
  - 返回的 `I_1C` 在放电后从剩余容量重新计算, 因此可能报告的是最终剩余电流基准, 而非运行时使用的初始 1C 基准。

## B. 隔膜模型集成

- 当隔膜参数通过构造函数传入时, 核心集成是正确的。
- 稳态求解器:
  - 根据 `I_app` 计算隔膜状态。
  - 放电时施加阴极侧电解质 Dirichlet 电位 `phi_e_sep < 0`。
  - 当隔膜启用时, 电池电压报告为 `phi_s(collector)`, 因为 Li/Li+ 参考电极位于隔膜外部, 电位为 0 V。
- 瞬态求解器:
  - 正常构建时与 `SteadyStateSolver` 共享隔膜对象。
  - 将隔膜阴极侧浓度作为电解质边界浓度施加。
- 测试覆盖: 零电流恒等式、隔膜开启电压降、隔膜浓度单调耗尽、禁用/默认向后兼容性。
- 注意: `parallel_validation.py` 应在隔膜开启运行时构建 `TransientSolver(..., separator=SeparatorParams(enabled=True))`。

## C. 3C = 30.6 mAh/g 的物理合理性

- 对于当前合成网络几何结构, 该值在物理上是合理的。
- 它不是 Khan 1CAL 的定量复现。
- 3C 提前截止与以下因素一致:
  - 合成网络几何结构而非 XCT 重建的 1CAL 拓扑,
  - 较低/改变的质量负载和活性界面连通性,
  - 高电流下的大局部极化,
  - 许多活性位点在阴极完全利用之前就达到了传输限制。
- 隔膜损失太小, 无法弥补差距:
  - 保存的并行结果表明 3C 下仅有约 `5.3 mV` 欧姆损失和 `2.0 mV` 浓度降。
  - 隔膜开启的 3C 容量仅从约 `30.6` 变为 `29.7 mAh/g`。
- 低 C 率下的时间受限容量不应被解释为真实的截止容量; 它们是电压未达到 3.0 V 时的运行时长伪影。

## D. 必须修改的内容 / 就绪状态

- 对于可发表/可展示的方法开发框架, 没有必须修改的核心求解器缺陷。
- 在将隔膜开启验证作为最终结果之前, 需修复 `scripts/parallel_validation.py` 中的隔膜构建和 `I_1C` 报告。
- 在声称定量复现 Khan Figure 4 之前, 几何/质量负载必须匹配论文 XCT 网络; 当前结果应定性地描述为方法验证。

## E. 单一最具影响力的改进

- 替换或校准合成立方网络, 使其与几何匹配:
  - 匹配 Khan 1CAL 质量负载,
  - 匹配电极厚度/投影面积,
  - 匹配活性界面面积和相连通性,
  - 最好使用 XCT 推导的拓扑或经校准的代理模型。
- 这将主导于进一步的隔膜调优之上, 因为隔膜仅贡献毫伏级变化, 而几何结构控制了大容量和高倍率极化的不匹配。
