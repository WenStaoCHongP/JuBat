# 源码分层重构迁移前示例基线

日期：2026-08-05

环境：Julia 1.11.2、1 thread、`GKSwstype=100`。

| 示例 | 迁移前结果 | 关键证据 |
|---|---|---|
| `change_current.jl` | FAIL | 无热模型下 `StandardVariables` 读取不存在的温度索引 |
| `change_model.jl` | FAIL | 同上 |
| `change_time_step.jl` | FAIL | 同上 |
| `mechanical_example.jl` | FAIL | 同上；尚未执行到旧绝对 CSV 路径 |
| `minimal_example.jl` | FAIL | 同上 |
| `testexample.jl` | PASS | PNG SHA-256 `3538fe6ab336f9852e90566b17edbd2cd6c2c14b93c8eeff5aed01f7037df9d5` |
| `thermal_example.jl` | PASS | V PDF `18e776385f0180aaa92ada4623c42ac67a549b829501293457799767a2f933ec`；T PDF `9479317e14eb36c017fdba4740a946f400a5dd582df5f2cda9b7072697b1f46e` |

## 前置缺陷修复后的成功基线

| 示例 | 状态 | 产物 SHA-256 |
|---|---|---|
| `change_current.jl` | PASS | `drive_UDDS_V.pdf`: `ebe0f821dc927c2a3f1181706d3397a5b295b16e60aabd46a237c67b28cd59f8` |
| `change_model.jl` | PASS | `change_model.pdf`: `2871735dc0c1a54de0fed2e8e7cc417fe60fc6ef710fcd051df453d453709355` |
| `change_time_step.jl` | PASS | `change_dt.pdf`: `561deafb0c5ffe97a9966c18fa8f73c4f0d42f587c939a9a8579ad17793863d0` |
| `mechanical_example.jl` | PASS | `mechanical_example_SPM.pdf`: `cd405f80af859d14f45c9f158e239835332d5d9e8c17230cf15d135f6027c07e` |
| `minimal_example.jl` | PASS | `minimal_example.pdf`: `75839d02d71c3eefbfc2d537b6eb782a22799986ce66cd85d3fd98ae5674064c` |
| `testexample.jl` | PASS | `testexample_results.png`: `3538fe6ab336f9852e90566b17edbd2cd6c2c14b93c8eeff5aed01f7037df9d5` |
| `thermal_example.jl` | PASS | V: `18e776385f0180aaa92ada4623c42ac67a549b829501293457799767a2f933ec`; T: `9479317e14eb36c017fdba4740a946f400a5dd582df5f2cda9b7072697b1f46e` |

完整输出保存在同目录的七份 `.log` 文件中。以上成功结果作为后续模块迁移回归基线。
