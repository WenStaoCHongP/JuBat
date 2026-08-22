# example/testexample.jl 代码简化基线

- **Baseline ID**: `testexample-20260815T011730-0600`
- **状态**: PASS（exit code 0）
- **入口**: `example/testexample.jl`
- **命令**: `D:/Julia-1.11.2/bin/julia.exe --startup-file=no example/testexample.jl`
- **环境**: Julia 1.11.2，Plots 1.40.9，1 thread，`GKSwstype=100`
- **Git HEAD**: `df8aab3e278243fa9efe0759b666f4e386628276`
- **脚本 SHA-256**: `8c19e11619444bac6954f01845eed676d281ec22b7aa40e2a88b7e367f06b9e9`
- **预检源码聚合 SHA-256**: `1cf823c5c77e23f85a8353b0af0e34a3ecdb2d6c21866f2334266a287bd867de`（兼容旧预检：字面量 `` `t`` 分隔）
- **标准 TSV 源码聚合 SHA-256**: `77cc5df31b1d1ce21c6329052546db8b12bb9f2cc3b6e74df7e1d76e36b9a0df`（TAB 分隔、LF 拼接）
- **说明**: 基线建立时工作树已有用户修改；因此以脚本哈希和 46 个 Julia 文件的内容清单为准，而不是仅以 HEAD 为准。

本基线取代 `testexample-20260806T031217-0600`。当前力学周向划分直接继承热网格实际角段，CZM 与 1682 个父热单元严格对应；同时源码清单纳入截至本次已审阅的严格循环状态契约。`testexample` 的 19 次 CZM 更新全部实际执行并收敛。

## 冻结结果

| 指标 | 基线值 |
|---|---:|
| thermal elements | 1682 |
| thermal nodes | 1763 |
| result time steps | 19 |
| initial voltage | 4.0367 V |
| final voltage | 3.9438 V |
| voltage drop | 0.0929 V |
| final capacity | 0.0833 Ah |
| minimum temperature | 298.15 K |
| maximum temperature | 299.00 K |
| final CZM D_max | 0.0000% |
| final CZM D_mean | 0.0000% |
| maximum normal separation | 1.2572e-13 m |
| fractured elements | 0 |
| CZM converged updates | 19 / 19 |
| result PNG SHA-256 | `b31ffb494c8e83701caad1b1181be1373b50236019a580b72c795c891ebd05ac` |

## 后续比较规则

每个代码简化批次完成后，必须用同一 Julia 版本、线程数、入口和参数重新运行。

- 必须 exit code 0。
- 网格规模、时间步数和上表全部科学结果必须与基线在脚本打印精度下完全一致。
- `output/testexample/testexample_results.png`（2026-08-21 起 testexample.jl 按 AGENTS.md §9.9 输出到该子目录）的 SHA-256 应一致；如图片哈希不同，即使打印指标相同也需检查曲线数据或绘图库环境变化。
- wall-clock、模块耗时、耗时占比不作严格相等要求，只记录趋势。
- 任一科学指标不一致时，该简化批次不得继续，先定位差异或回滚。
- 本次未生成 `output/simple_coupling_debug.log`，因此该文件不属于基线。

## 文件

- `metrics.toml`: 机器可读的关键指标与比较策略。
- `source_manifest.tsv`: 运行时所有 Julia 源文件和入口脚本的 SHA-256。
- `preflight.log`: 运行前身份和旧输出状态。
- `run.log`: 完整控制台输出。

## 基线 v2（2026-08-22 重冻结，用户宏观参数修正）

- 触发：用户修正 `src/parameters/Jellyroll.jl` 宏观力学参数（PCC.E→70 GPa/0.33、NCC.E→110 GPa/0.34、PE/NE.E_coat→1 GPa、SP.E→750 MPa/0.35；提交 `1a74411`）。电化学/热全部指标与 v1 一致；仅两项力学输出移动：`maximum normal separation 1.2557e-14 → 1.2572e-13 m`（刚性箔下微小变形放大一个量级，仍为 ~零损伤工况）与 PNG SHA（上表已更新为 v2 值）。
- v1 SHA：`4ba6207c…e932`（历史比对记录见 `Simplify/baseline.md` 批次表）。
