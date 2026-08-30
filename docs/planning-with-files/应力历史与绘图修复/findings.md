# 应力历史与绘图修复：发现与决策

## 用户决策

- 修复应力历史错位和零应力问题。
- SP/PCC/NCC 热膨胀系数暂时保持 0，后续另行添加。
- 删除额外径向散点图，只保留三张最终场云图。

## 已知问题来源

- `Variable_update!` 先写第 `v` 列，`t` 随后推进，CZM 与 `export_macro_stress` 最后运行，导致新 CZM 解写回旧时间标签。
- 应力历史矩阵全零预分配，但只在 `czm_update_interval` 命中时写入；非更新列因此是“未计算零”，不是物理解。
- `couple_example.jl` 仍创建并保存 `final_hoop_stress_radial_profile.png`。

## 目标语义

- 每次 CZM 更新必须使用当前 `CallModel` 的 SOC 与同一时刻的活动温度场，并在时间推进前完成。
- 每个被记录的输出列都写入 `latest_macro_stress`：更新步刷新它，非更新步保持最近一次实际求解结果。
- 初始列由当前 `u_prev` 和初始载荷显式恢复，不依赖数组初始化零。

## 调用链事实

- `Solve` 当前在 `Variable_update!` 记录第 `v` 列后更新 `y_old`、推进 `t`，最后才执行 CZM 与应力收割；快照使用推进后时间，应力却写回旧列。
- 当前循环已在求解后把 `T_nodes_carry` 绑定到 `y_c` 的新温度，而 SOC 仍来自本轮 `CallModel(y_old, t)`；修复需在覆盖温度前保留当前时间层温度，并在推进 `t` 前执行 CZM。
- 仅把非更新列复制上一列不足以覆盖 `outputType != auto`：CZM 可能在没有新历史列时更新。需要维护独立的 `latest_macro_stress`，记录输出列时再写入。
- 初始 `CallModel` 后已具备 SOC、温度和 `variables_hist`；可创建 `CzmLayout` 后显式恢复初始应力列。

## 散点图与基线范围

- 额外散点图代码集中在 `example/couple_example.jl` 的径向坐标、`scatter`、保存路径和提示文本。
- 当前生成物在 `output/couple_example/` 与遗留 `output/testexample/` 各有一份径向剖面 PNG。
- 当前规则/基线引用包括 `AGENTS.md`、`Simplify/baseline/testexample/README.md`、`metrics.toml` 和 `preflight.log`；历史进度记录保留其当时事实，不改写。

## 已实施设计

- `Mechanical.jl` 将“计算当前应力”和“写历史列”拆开：`compute_macro_stress` 返回四分量，`write_macro_stress!` 负责显式写列，原 `export_macro_stress` API 保持可用。
- `Solve` 在初始列显式计算应力；每次实际 CZM 更新刷新 `latest_macro_stress`，每次输出都写入它。因此非更新列是最近一次真实力学解的保持值。
- CZM 更新块已移动到 `t += dt` 之前，并使用本次 `CallModel` 返回的当前温度副本及当前 SOC；快照时间、应力载荷和历史列时间一致。
- `couple_example.jl` 已删除径向坐标、分组、`scatter`、第四张保存路径和提示；AGENTS 当前门禁改为三张 PNG。
- `SP/PCC/NCC.alphaT` 未修改。
- 当前值明确核对为 `SP.alphaT=PCC.alphaT=NCC.alphaT=0`；PE/NE 既有涂层系数不属于本轮新增或修改范围。

## 快速基线差异归因

- 60 s 快速门 exit 0，热网格 1682/1763、19 步、电压、容量、温度、零损伤和断裂数均与 v7 一致。
- 同时间层修复后，最大法向分离由 `9.4407e-13` 调整为 `9.6486e-13 m`；环向应力由 `-1.3952~3.8570` 调整为 `-1.3954~3.8603 MPa`；切向剪应力由 `-0.37392~0.39515` 调整为 `-0.37425~0.39548 MPa`。变化仅限修复直接影响的力学输出。

## 三张云图验证

- `couple_example.jl` exit 0，输出目录中 PNG 数量严格为 3，径向散点文件不存在。
- 温度图保持原 SHA-256 `540fe42f...f978c`；环向应力新 SHA-256 为 `9726a7c8...050d`；切向剪应力新 SHA-256 为 `b3b43d7c...5445`。
- 视觉检查：温度仍为热网格面云图；环向和切向剪应力仍呈细密力学层条纹，均为 Q4 面填色而非散点，坐标/单位/色标正常。
- 当前三张图完整哈希：温度 `540fe42f1039cb4446046a6da96f1ecf18d162b643bc40a40a099ca2b93f978c`，环向 `9726a7c84cf59a932a2cf45346770c6c046a6a81a7ca36f951714339e558050d`，切向剪应力 `b3b43d7cee67ed2a6fca33e0fcffee899629b84998852ba05aeb08fdaa105445`。

## 最终一致性检查

- `metrics.toml` 可由 Julia `TOML.parsefile` 正常读取，当前基线 ID 为 `testexample-20260830T172856+0800`。
- `source_manifest.tsv` 的 47/47 个文件哈希匹配；标准 TSV 聚合 SHA-256 复算为 `cbc96a2d665d498a458a51b558e87c6568014e73e1c9cf688f89f6aff34bab3a`，清单文件 SHA-256 为 `3e3fc6009159c614f30645b78a6912572e0a09e17782974f210521b8618e73ee`。
- `test/test_czm_strain_inputs.jl` 20,194/20,194 通过，确认分层应变输入映射未因时序修改回归。
