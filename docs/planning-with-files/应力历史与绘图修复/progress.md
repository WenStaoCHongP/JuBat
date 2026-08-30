# 应力历史与绘图修复：进度记录

## 会话：2026-08-30

### Phase 1：调用链与基线核实

- **状态：** completed
- 已确认本轮授权范围和 `alphaT=0` 保持决定。
- 已建立独立规划记录，后续同步更新总索引。
- 已定位时间错位根因和 `outputType` 条件下不能只复制上一列的边界。
- 已盘点散点图源码、两处生成物及当前基线/规则引用。

### Phase 2：窄范围实现

- **状态：** completed
- 计划以独立 `latest_macro_stress` 表示最近一次实际力学恢复结果；CZM 更新刷新它，每个输出列显式写入它。
- 已完成 Mechanical/Solve 时序和历史写入修改，并保持 `export_macro_stress` 入口兼容。
- 已新增 `czm_update_interval=2` 的小型端到端回归：检查快照时间属于结果时间轴、非更新列严格等于最近更新列且非零。
- 已删除示例中的额外散点图代码，并把当前绘图门改为三张 PNG。

### Phase 3：验证与交付

- **状态：** completed
- 层分辨应力专项与新增时间/间隔回归共 20/20 通过。
- CZM 分层应变输入分发测试 20,194/20,194 通过。
- 60 s 快速门和完整绘图例均 exit 0；当前输出严格为三张云图。
- 47 个源码清单条目逐一匹配，标准 TSV 聚合哈希与清单文件哈希均复算一致。
- 总规划索引已同步为 36 个任务目录、121 个文件。

## 测试结果

| 检查 | 预期 | 实际 | 状态 |
|---|---|---|---|
| `test/test_layer_resolved_stress.jl` | 原专项 + 间隔/时间回归通过 | 20/20 通过，其中新增回归 8/8 | pass |
| 散点代码静态搜索 | `couple_example.jl` 零 `scatter/radial_plot` | 零匹配 | pass |
| `SP/PCC/NCC.alphaT` | 保持 0 | 本轮未修改，现值仍为 0 | pass |
| `git diff --check`（相关文件） | 本轮新增行干净 | 被工作树原有 `Solve.jl` 行尾空白阻塞 | blocked-by-existing |
| `example/testexample.jl` 60 s 快速门 | 非力学指标不变，修复相关力学指标可解释变化 | exit 0；网格/19 步/电热/零损伤不变；分离与应力小幅变化 | pass-authorized-change |
| `example/couple_example.jl` | exit 0，只产三张最终云图 | exit 0，PNG 数量 3，径向散点不存在 | pass |
| 三张 PNG 视觉检查 | 温度热网格、应力力学子网格、无散点 | 三张均符合 | pass |
| 源码清单 | 47 文件哈希同步 | 聚合 `cbc96a2d...bab3a`，清单文件 `3e3fc600...e73ee` | pass |
| `test/test_czm_strain_inputs.jl` | 分层材料输入映射不回归 | 20,194/20,194 通过 | pass |
| 基线 TOML/清单复核 | 可解析且每个文件哈希匹配 | TOML 可解析；47/47 匹配；聚合复算一致 | pass |

## 错误日志

| 时间 | 错误 | 尝试 | 处理 |
|---|---|---:|---|
| 2026-08-30 | `git diff --check` 报 `Solve.jl` 大量行尾空白 | 1 | 属于本轮前已存在的混合 diff；不整文件格式化，改用本轮行范围检查并在交付中披露 |
| 2026-08-30 | TOML 校验脚本误按不存在的 `[baseline]` 表读取 | 1 | 按顶层 `baseline_id` 等键重新解析，文件有效 |
