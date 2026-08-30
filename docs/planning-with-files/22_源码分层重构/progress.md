# JuBat 源码物理分层重构进度

## Session: 2026-08-05

### Phase 0：检查点、分支与基线

- **Status:** complete
- 已完整读取 `planning-with-files` 技能。
- 已确认规划文件保存约定为 `docs/planning-with-files/<中文任务名>/`。
- 已确认当前分支、HEAD 和大规模脏工作树状态。
- 已创建本任务的 `task_plan.md`、`findings.md`、`progress.md`。
- 仓库外检查点已生成 tracked binary patch、Git 状态、身份文件及包含 43 个未跟踪文件的 zip；zip SHA-256 为 `488C0B3E0B53662E2504974B763D1F471F4E7735E43DC0D6414847D804E7D040`。
- 首次创建重构分支被 `.git/refs` 沙箱写权限阻止，待受控权限重试。
- 已通过受控权限创建并切换到 `codex/src-physics-modularization`。
- 迁移前 22/22 个现有测试文件通过；用户随后将强制旧示例回归范围明确为 7 个指定脚本。
- 已运行全部 7 个指定示例：`testexample` 与 `thermal_example` 通过，其余 5 个因同一既有温度索引缺陷失败。

## Test Results

| Test | Result | Status |
|---|---|---|
| 当前 Git 状态预检 | 分支 `czm-refactor`，HEAD `df8aab3`，存在大量既有改动 | PASS |
| 现有测试健康度 | 22 个 `test/*.jl` 独立 Julia 进程 | PASS（22/22） |
| `change_current.jl` preflight | `StandardVariables` 缺少 `temperature` 键 | FAIL（既有缺陷） |
| `change_model.jl` preflight | 同一 `temperature` 键缺陷 | FAIL（既有缺陷） |
| `change_time_step.jl` preflight | 同一 `temperature` 键缺陷 | FAIL（既有缺陷） |
| `mechanical_example.jl` preflight | 同一缺陷，尚未执行到旧绝对 CSV 路径 | FAIL（既有缺陷） |
| `minimal_example.jl` preflight | 同一 `temperature` 键缺陷 | FAIL（既有缺陷） |
| `testexample.jl` preflight | 冻结科学指标与 PNG SHA-256 精确一致 | PASS |
| `thermal_example.jl` preflight | exit 0，V/T PDF 已生成 | PASS |

### Phase 1：标准 Julia package 与架构测试

- **Status:** in_progress
- 已添加 `test/test_example_preconditions.jl`，先固定无热模型温度分配与 mechanical 数据路径要求。
- 首轮修复后 4 个此前失败示例通过；`change_model` 在 SPMe 分支暴露第二处温度索引缺陷，已扩展测试。
- `Variables.jl` 统一分配 1 行温度结果，不向无热模型状态向量添加 DOF；SPMe 两个变量入口在无热 DOF 时使用环境初温。
- 已将该缺陷固化为重构约束：固定结果维度禁止使用 `haskey(case.index, ...)` 推导；`variables["temperature"]` 必须直接赋值为 `1 × num`，无热模型禁止访问或伪造温度状态索引。
- 已审计整个 `src/` 并将 SPM、SPMe、P2D 的 4 处温度索引存在性判断收敛到 `representative_temperature`；普通工作区和网格映射的成员检查不受该约束影响。
- `mechanical_example.jl` 改用仓库相对数据路径并删除不存在的 `JuBat.Citation()` 调用。
- 聚焦测试扩展为 20/20 通过；五个原失败示例全部 exit 0。
- 全源码温度状态访问统一后，用户指定的 7 个示例全部 exit 0。
- `testexample` 仍为 1682 个热单元、1763 个节点、19 步、4.0367→3.9438 V、0.0833 Ah、298.15–299.00 K，PNG 88,744 字节且 SHA-256 为 `3538fe6ab336f9852e90566b17edbd2cd6c2c14b93c8eeff5aed01f7037df9d5`。
- 九个子模块创建批次已按用户要求完整回退；保留其之前已完成的 package 文件与温度访问修复。
- 已审计 `src/` 中 `?: nothing`/可选缓存赋值：删除 ThermalDistributed 的三处冗余 geometry 回退和 CallModel 的伪 `I_e_cache` 探测；保留 CZM 等真实可选参数分支。
- 进一步将 `MeshGeometry.boundary_edges` 收紧为必需的 `BoundaryEdgeCache`，并删除 CallModel/ThermalDistributed 对固定 `czm_element_map` 字段的重复存在性判断。
- 缓存不变量测试 11/11、热边界原位测试 18/18 通过；用户指定的 7 个示例全部通过。最后一次 `testexample` 复验仍为 1682/1763/19，PNG 88,744 字节且 SHA-256 精确一致。

## Files Created

- `docs/planning-with-files/22_源码分层重构/task_plan.md`
- `docs/planning-with-files/22_源码分层重构/findings.md`
- `docs/planning-with-files/22_源码分层重构/progress.md`
