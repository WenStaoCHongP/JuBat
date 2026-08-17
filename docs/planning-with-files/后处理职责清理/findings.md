# 后处理职责清理发现

## Requirements

- 审计 `PostProcessing.jl`、`CsvExport.jl`、`CycleData.jl`、`CzmPostProcess.jl`。
- 重命名下划线前缀函数并更新全部调用者。
- 检查职责交叉与功能重复。
- 将 `PostProcessing.jl` 第 117 行后的函数分配到更合适的文件。
- 保持数值与输出契约不变。

## Research Findings

- 实际源码规模：`PostProcessing.jl` 350 行、`CsvExport.jl` 620 行、`CycleData.jl` 504 行、`CzmPostProcess.jl` 117 行；`CsvExport.md` 仍记录 637 行，已过期。
- `PostProcessing.jl` L1-L117 只有通用 `PostProcessing(case, variables, v)` 结果单位还原；L119-L350 全部是循环阶段/循环汇总、SOH、终止判据、打印和绘图逻辑，与通用结果还原职责明显交叉。
- `PostProcessing.jl` L119-L309 有 9 个下划线前缀 helper；L310-L350 的 `plot_cycling_results` 已是公开名称。
- `CsvExport.jl` 有 15 个下划线前缀函数/方法，覆盖配置验证、采样策略、受保护写入、七类文件 writer、安全索引、循环结果查找和单元面积计算。
- `CycleData.jl` 与 `CzmPostProcess.jl` 当前没有下划线前缀函数。
- `CzmPostProcess.jl` 已专注 CZM 损伤统计、判据、重置/累积和 variables 映射；没有明显 CSV 或循环汇总逻辑。
- `PostProcessing.jl` L119-L309 的 9 个 helper 只被自身或 `CycleSolver.jl` 调用；`plot_cycling_results` 是循环结果专属公开绘图入口。因此 L119-L350 可整体迁入新的 `CyclePostProcess.jl`，而不是污染 `CycleData.jl` 或 `CzmPostProcess.jl`。
- `CsvExport.jl` 的下划线函数除 `write_csv_guarded!` 的测试调用外均为文件内调用；可一次性改名并同步测试，不需要兼容 façade。
- `CycleData.jl` 的 `export_cycle_data_to_csv` 属于 CSV 序列化职责，和专职的 `CsvExport.jl` 存在文件职责交叉；其数据容器和在线求解采集仍应留在 `CycleData.jl`。
- 全仓未发现第二个同名或同算法的 Q4 Shoelace 面积 helper；`CsvExport.compute_csv_element_areas` 属于 CSV 自包含几何还原，不与 `MeshGeometry` 的缓存构建实现直接重复，本轮保留算法并仅重命名。
- `CycleData.export_cycle_data_to_csv` 与 `CsvExport.export_cycling_csv` 都写 CSV，但输入边界不同：前者接收单循环 `CycleExportData` 并写原始快照/网格六文件，后者接收多循环 `CyclingResult` 并写七类分析表。应统一文件归属而不合并两个公开入口。
- 实施后文件规模：`PostProcessing.jl` 117 行、`CyclePostProcess.jl` 235 行、`CycleData.jl` 410 行、`CsvExport.jl` 718 行、`CzmPostProcess.jl` 117 行。
- 指定后处理范围内已无下划线前缀函数定义，源码/测试/示例中也无旧函数名调用残留。
- 现有定向测试只有 `test_csv_export_guard.jl` 和 `test_cycle_data_import_removed.jl` 覆盖 CSV guard 与单循环导出；需要新增职责架构测试覆盖文件归属和旧名消失。
- `PostProcessing.md` 仍描述迁移前的 351 行/11 函数；`CsvExport.md` 仍描述旧 637 行和下划线名称；`CycleData.md` 仍把 CSV 序列化列为自身职责，三者均需重建。
- `CycleSolver.md` 有 4 处段落仍引用旧下划线名称，并错误建议 helper 可放 CsvExport/PostProcessing；应改为指向 `CyclePostProcess.jl`。
- `JuBat.jl` 新增 include 后为 94 行；源码入口索引中的 include/export 行号整体需要同步。

## Technical Decisions

| Decision | Rationale |
|---|---|
| 以实际调用图和数据类型归属决定目标文件 | 避免仅按函数名或行号移动造成加载依赖倒置 |
| 暂不把循环汇总 helper 移入 `CycleData.jl` | `CycleData.jl` 当前职责是在线快照采集，不应混入 `CyclingResult` 汇总或 CSV 序列化 |
| 新增 `CyclePostProcess.jl` 承接 `PostProcessing.jl` L119-L350 | 循环汇总、SOH、终止、打印和绘图是独立职责；放回 `CycleSolver.jl` 会令求解器进一步膨胀 |
| 将 `export_cycle_data_to_csv` 从 `CycleData.jl` 移至 `CsvExport.jl` | 统一所有 CSV 序列化入口，消除 CycleData 与 CsvExport 的文件职责交叉 |
| 下划线函数不保留旧名别名 | 所有调用点均在仓库内可控；清除 façade 才能真正落实命名约束 |

## Rename Map

### CyclePostProcess

- `_phase_termination_symbol` → `phase_termination_symbol`
- `_state_concentration_variance` → `state_concentration_variance`
- `_postprocess_phase_result` → `postprocess_phase_result`
- `_postprocess_cycle_result!` → `postprocess_cycle_result!`
- `_append_cycle_result!` → `append_cycle_result!`
- `_update_soh_and_capacity!` → `update_soh_and_capacity!`
- `_print_cycle_summary` → `print_cycle_summary`
- `_check_cycle_termination` → `check_cycle_termination`
- `_print_cycling_summary` → `print_cycling_summary`

### CsvExport

- `_validate_csv_opt` → `validate_csv_options`
- `_should_output_step` → `should_export_step`
- `_compute_last_snap_indices` → `compute_last_snapshot_indices`
- `_should_output_snapshot` → `should_export_snapshot`
- `_write_csv_guarded!` → `write_csv_guarded!`
- `_write_cycle_summary` → `write_cycle_summary_csv`
- `_write_element_currents` → `write_element_currents_csv`
- `_write_node_temperature` → `write_node_temperature_csv`
- `_write_cohesive_damage` → `write_cohesive_damage_csv`
- `_write_node_displacement` → `write_node_displacement_csv`
- `_write_driving_force` → `write_cohesive_driving_force_csv`
- `_write_czm_diagnostics` → `write_czm_diagnostics_csv`
- `_safe_get` → `safe_csv_matrix_value`
- `_find_solve_result` → `find_cycle_solve_result`
- `_compute_element_areas` → `compute_csv_element_areas`

## Issues Encountered

| Issue | Resolution |
|---|---|
| 暂无 | — |

## Validation Findings

- 职责架构测试、CSV 防护测试和单循环 CSV 公共路径测试全部通过，共 51 项断言。
- `example/testexample.jl` 保持 1682 个热单元、1763 个节点和 19 个结果步；电压、容量、温度及 CZM 指标与冻结基线一致。
- 输出 PNG 为 88,744 字节，SHA-256 为 `3538fe6ab336f9852e90566b17edbd2cd6c2c14b93c8eeff5aed01f7037df9d5`。
- 指定五个后处理文件无下划线前缀函数定义，活动源码、测试和示例无旧名称调用；`export_cycle_data_to_csv` 仅在 `CsvExport.jl` 定义一次。

## Resources

- `src/PostProcessing.jl`
- `src/CsvExport.jl`
- `src/CycleData.jl`
- `src/CzmPostProcess.jl`
- `md/源码函数索引/PostProcessing.md`
- `md/源码函数索引/CsvExport.md`
- `md/源码函数索引/CycleData.md`
- `md/源码函数索引/CzmPostProcess.md`
