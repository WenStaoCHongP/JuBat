# CZM 最大分离位移核查进度

## Session: 2026-09-01

### Phase 1：数据流与现状复现

- **Status:** completed
- **Started:** 2026-09-01
- Actions taken:
  - 确认用户要求核查最大分离并运行单元级脚本。
  - 选择系统化调试：数据流追踪、受控单元验证、端到端独立复算。
  - 记录当前示例值及边界修改背景。
  - 完成全局符号索引，初步区分当前步 `CZMResult.separation_n` 与历史峰值 `DamageState.δ_max_n`。
  - 追踪局部标架、2 点线积分、损伤状态更新和最终尺度还原；识别出混合模式下 `δ_max_n` 受 `δ_max_eff` 门控的潜在一般性语义风险。
  - 确认 `Solve` 在力学更新前生成峰值摘要，更新后未刷新摘要；且未把返回的 `CZMResult` 写回逐单元 CZM 变量，待实跑确认输出影响。
- Files created/modified:
  - `docs/planning-with-files/43_CZM最大分离位移核查/task_plan.md`
  - `docs/planning-with-files/43_CZM最大分离位移核查/findings.md`
  - `docs/planning-with-files/43_CZM最大分离位移核查/progress.md`

### Phase 2：单元级与手算交叉验证

- **Status:** completed
- Actions taken:
  - 运行现有 cohesive 网格/法向测试、尺度重设计测试和分界面双线性本构测试，全部通过。
  - 新增并运行诊断产物 `output/czm_separation_audit/unit_probe.jl`：三个不同位置/界面类型的非均匀分离与独立手算、装配值和物理换算一致，刚体平移为零。
  - 用非比例混合模态序列复现 `DamageState.δ_max_n` 受 `δ_max_eff` 门控而遗漏法向峰值的语义问题。
- Files created/modified:
  - `output/czm_separation_audit/unit_probe.jl`

### Phase 3：端到端定位与敏感性核查

- **Status:** completed
- Actions taken:
  - 按用户更正中止旧边界配置下的 3600 s 探针。
  - 修改 `identify_bc_nodes_czm`，使 `fix_inner=false` 也共用分层端点逻辑；直接计数确认净增 23 点。
  - 现有 `true` 模式测试 11/11、因子分解缓存测试 17/17、弧长/几何测试 32/32 通过。
  - 按用户要求将 `testexample.jl` 与探针改回 60 s，中止尚未完成的 3600 s 运行。
  - 完成 60 s 快照、节点手算、helper 与结果字典对齐：定位到最后时间层的 NE–NCC 单元 5127。
  - 确认真实峰值 `1.101106737768e-12 m`；公开结果少报为 `8.730984248720e-13 m`，根因为峰值摘要滞后一步，同时逐单元分离历史未写回。
- Files created/modified:
  - 暂无。

### Phase 4：结论与交付

- **Status:** completed
- Actions taken:
  - 汇总边界逻辑修改、单元级正确性、结果输出滞后缺陷与预应力/屈服语义。
  - 保留分离输出链的已知缺陷，未经用户授权不修复。
- Files created/modified:
  - 暂无。

## Test Results

| Test | Input | Expected | Actual | Status |
|---|---|---|---|---|
| `test_create_czm_mesh.jl` | 实际 Jellyroll CZM 网格 | 外向定向与纯张开分离正确 | 80,747/80,747 pass | pass |
| `test_czm_scale_redesign.jl` | 尺度锚点与装配换算 | `Λ=L/δ_czm` 与物理值一致 | 28/28 pass | pass |
| `test_bilinear_per_interface.jl` | PE–PCC/NE–NCC 本构 | 张开、混合模式、批更新通过 | 71/71 pass | pass |
| `unit_probe.jl` | 非均匀端点跳跃、刚体平移、混合历史 | 独立手算与尺度一致；识别历史语义 | 分离/尺度 pass；历史语义缺口复现 | mixed |
| `test_czm_layer_endpoint_bc.jl` | `fix_inner=true` | 复用现有端点集合测试 | 11/11 pass | pass |
| `fix_inner=false` 直接计数 | nθ=80 | 原外圈 +23 | 80 → 103 | pass |
| `test_czm_factorization_cache.jl` | basic + `fix_inner=false` | 缓存与直接解一致 | 17/17 pass | pass |
| `test_czm_arc_geo.jl` | 弧长/自由芯部路径 | 受新边界影响的路径仍通过 | 32/32 pass | pass |
| `end_to_end_probe.jl` | 60 s，`fix_inner=false` +23 | 快照/手算/helper 一致，辨识输出时序 | 单元值一致；滞后误差 0 | pass with output defect |
| `example/testexample.jl` | 正式 60 s 快速基线 | 退出码、网格/步数、文字指标 | exit 0，19 步，打印 `8.7310e-13 m` | pass; printed max is stale |

## Error Log

| Timestamp | Error | Attempt | Resolution |
|---|---|---:|---|
| 2026-09-01 | 暂无 | 1 | — |

## 5-Question Reboot Check

| Question | Answer |
|---|---|
| Where am I? | Phase 1：数据流与现状复现 |
| Where am I going? | 单元级独立复算、端到端最大值定位、结论交付 |
| What's the goal? | 判断最大法向分离位移是否计算正确 |
| What have I learned? | 当前值随新增边界改变，但尚不能据此判断正确性 |
| What have I done? | 建立调查记录并确定证据链 |
