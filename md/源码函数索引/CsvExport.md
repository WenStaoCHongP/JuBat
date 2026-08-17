# CsvExport.jl

- **源文件**: `src/CsvExport.jl`
- **行数**: 718 行
- **函数/struct 计数**: 1 个 struct、2 个便捷构造器、18 个函数/方法定义
- **职责**: JuBat CSV 序列化单一归属——导出多循环分析结果与单循环原始快照；采样、越界保护和单文件失败隔离均在本文件实现。
- **相关技术文档**: `md/10_参数传递与模块架构.md`、`md/06_内聚力模型_CZM.md`

## 数据结构

### `struct CsvExportOptions` — L20-L25

配置多循环 CSV 采样方式：`mode`、`save_every`、`full_output_cycles`、`skip_files`。

### 便捷构造器 — L27-L28

提供默认完整输出与指定模式两种构造方式。

## 多循环 CSV 入口与 helper

### `validate_csv_options(csv_opt)` — L30-L37

验证模式和自定义采样步长。

### `should_export_step(...)` — L40-L46

判断循环时间步是否需要输出。

### `compute_last_snapshot_indices(czm_snapshots)` — L50-L65

计算每个 cycle/phase 分组最后一条 CZM 快照索引。

### `should_export_snapshot(...)` — L68-L75

按采样配置判断快照是否输出。

### `write_csv_guarded!(write_fn, filename, files_written, files_skipped)` — L81-L91

隔离单个 CSV writer 的异常并维护成功/跳过列表。

### `export_cycling_csv(result, case, czm_mesh; ...)` — L100-L192

多循环分析 CSV 主入口，最多写出七类文件。

### Writer 函数

| 函数 | 行号 | 输出 |
|---|---:|---|
| `write_cycle_summary_csv` | L198-L225 | `cycle_summary.csv` |
| `write_element_currents_csv` | L231-L298 | `element_currents.csv` |
| `write_node_temperature_csv` | L304-L351 | `node_temperature.csv` |
| `write_cohesive_damage_csv` | L357-L405 | `cohesive_damage.csv` |
| `write_node_displacement_csv` | L411-L445 | `node_displacement.csv` |
| `write_cohesive_driving_force_csv` | L452-L551 | `cohesive_driving_force.csv` |
| `write_czm_diagnostics_csv` | L557-L573 | `czm_solver_diagnostics.csv` |

### 共享 CSV helper

- `safe_csv_matrix_value` — L579-L585：二维数值矩阵安全访问及非匹配类型 fallback。
- `find_cycle_solve_result` — L588-L600：按循环和阶段查找 `solve_result`。
- `compute_csv_element_areas` — L606-L620：用 Shoelace 公式计算 Q4 归一化面积。

## 单循环原始数据入口

### `export_cycle_data_to_csv(export_data, output_dir; prefix="cycle")` — L639-L718

接收 `CycleExportData` 并写出：

1. `{prefix}_timesteps.csv`
2. `{prefix}_T_nodes.csv`
3. `{prefix}_soc_n.csv`
4. `{prefix}_soc_p.csv`
5. `{prefix}_mesh_nodes.csv`
6. `{prefix}_mesh_elements.csv`

公开函数名、CSV 列名、数值格式和返回路径 tuple 均保持不变。

## 省略项

无。全部 struct、构造器和函数/方法均有条目。

### [DEBUG]

- 多循环 writer 的 Skipping/Written 与汇总输出是用户可见状态。
- 单循环导出在六个文件写出后分别打印保存路径。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|---|---|---|
| L210 | SOH 越界时写 `NaN` | 下游统计必须处理 NaN |
| L484-L486 | driving-force 参数仍为零占位 | 对应 CSV 会提前跳过 |
| L583/L585 | 安全矩阵访问越界或类型不匹配返回 `NaN` | 上游字段错误可能延迟到 CSV 检查才暴露 |
| L393 | `D >= 0.95` 判为 fractured | 与其他损伤阈值语义不同 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|---|---|---|
| L109-L183 | 七类文件具有相似的 skip/guard/write 调度结构 | 可进一步表驱动，但需保留每类前置条件 |
| L280-L288 | 九个字段连续安全访问 | 可抽取 NamedTuple 快照 |
| L386-L404 | cohesive 快照有多组长度守卫 | 可抽取安全向量读取 helper |
| L493-L556 | driving-force 输出包含多层数据存在性与索引检查 | 参数接通后再单独简化 |
