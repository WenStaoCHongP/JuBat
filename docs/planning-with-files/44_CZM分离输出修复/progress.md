# CZM 分离输出修复进度

## Session: 2026-09-01

### Phase 1：输出契约 RED 测试

- **Status:** completed
- Actions taken:
  - 用户确认保留原法向键、新增有效峰值键，并按模式打印。
  - 追踪现有状态、变量初始化、求解后更新、后处理与示例消费者。
  - 新增 `test/test_czm_separation_output.jl`，受控结果预期有效分离 `[5,3]`，首次运行得到 2 个预期失败。0 error。
- Files created/modified:
  - `docs/planning-with-files/44_CZM分离输出修复/task_plan.md`
  - `docs/planning-with-files/44_CZM分离输出修复/findings.md`
  - `docs/planning-with-files/44_CZM分离输出修复/progress.md`
  - `test/test_czm_separation_output.jl`

### Phase 2：最小生产修复

- **Status:** completed
- Actions taken:
  - 新增有效分离的变量历史、逐单元场和物理单位还原。
  - 在 CZM 收敛后、历史写入前刷新当前场及 `δ_max_n/δ_max_eff`。
  - 增加模式到最大分离结果键的严格映射，更新两个正式示例的打印。
  - GREEN 复跑全部 11 个断言通过。

### Phase 3：验证与交付

- **Status:** completed
- Actions taken:
  - 完成尺度、本构、basic 缓存回归，全部通过。
  - 完成 60 s `testexample.jl` mix 文字基线：最大有效分离 `1.1207e-12 m`。
  - 完成 60 s `couple_example.jl` model1 求解与绘图：最大法向分离 `1.1011e-12 m`。
  - 更新诊断探针为修复后契约，全尺寸比对法向/有效值的结果历史、快照、状态和节点手算，全部一致。
  - 完成前复跑新增契约测试 11/11，并执行作用域 `git diff --check`，退出码 0（仅 Git 行尾转换提示）。

## Test Results

| Test | Expected | Actual | Status |
|---|---|---|---|
| `test_czm_separation_output.jl` RED | 缺失有效分离键导致失败 | 2 fail, 0 error | expected fail |
| `test_czm_separation_output.jl` GREEN | 有效值、模式映射、本步无滞后 | 11/11 pass | pass |
| `test_czm_scale_redesign.jl` | 尺度与装配换算 | 28/28 pass | pass |
| `test_czm_factorization_cache.jl` | basic 求解与缓存 | 17/17 pass | pass |
| `test_bilinear_per_interface.jl` | 法向/混合本构 | 71/71 pass | pass |
| `example/testexample.jl` | 60 s mix 打印有效峰值 | exit 0, 19 步, `1.1207e-12 m` | pass |
| `example/couple_example.jl` | 60 s model1 打印法向峰值与绘图 | exit 0, `1.1011e-12 m`, 3 PNG | pass |
| `end_to_end_probe.jl` | nθ=80 快照/历史/状态/手算 | 法向误差 0，有效误差 `1.26e-29 m` | pass |

## Error Log

| Timestamp | Error | Attempt | Resolution |
|---|---|---:|---|
| 2026-09-01 | 暂无 | 1 | — |
