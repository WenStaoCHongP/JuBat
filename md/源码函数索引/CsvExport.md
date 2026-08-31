# CsvExport.jl

- **源文件**: `src/CsvExport.jl`
- **行数**: 750 行
- **函数/struct 计数**: 22 个顶层定义（1 struct + 21 function，含便捷构造器）
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
| `write_cycle_summary_csv` | L198-L221 | `cycle_summary.csv` |
| `write_element_currents_csv` | L227-L282 | `element_currents.csv` |
| `write_node_temperature_csv` | L289-L331 | `node_temperature.csv` |
| `write_cohesive_damage_csv` | L338-L400 | `cohesive_damage.csv` |
| `write_node_displacement_csv` | L406-L442 | `node_displacement.csv` |
| `write_cohesive_driving_force_csv` | L450-L547 | `cohesive_driving_force.csv` |
| `write_czm_diagnostics_csv` | L555-L570 | `czm_solver_diagnostics.csv` |

`cohesive_damage.csv` 的角度列 `theta_cum_deg`（2026-08-29 由 `theta_deg` 改名）为半径反解**累计角** `(r−a)/b·180/π`（同 `edge_boundary`，`a=cell.Rin`、`b=cell.layer/2π`），多匝单调、无 atan ±π 接缝；同一物理方位在不同匝上相差整数个 360°。

### 共享 CSV helper

- `require_csv_solve_result` — L576-L581：阶段结果存在性校验，缺失即报错（不静默跳过）。
- `require_csv_vector` — L583-L590：按 key 取一维向量，长度/类型不符报错。
- `require_csv_matrix` — L592-L601：按 key 取二维矩阵，形状不符报错。
- `require_csv_length` — L603-L609：数组长度守卫。
- `existing_cycle_phases` — L611-L617：列出已有阶段的 phase 结果。
- `find_cycle_solve_result` — L620-L635：按循环和阶段查找 `solve_result`。
- `compute_csv_element_areas` — L638-L652：用 Shoelace 公式计算 Q4 归一化面积。

## 单循环原始数据入口

### `export_cycle_data_to_csv(export_data, output_dir; prefix="cycle")` — L671-L750

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
| L394 | `D >= 0.95` 判为 fractured | 与其他损伤阈值语义不同 |

（原表中 SOH 越界写 NaN、driving-force 零占位、安全矩阵访问 NaN fallback 三行对应的代码已不存在——现为 `require_*` 家族显式报错，2026-08-29 核实删除。）

### [COMPLEX-CHECK]

| 位置 | 内容 | 简化建议 |
|---|---|---|
| `export_cycling_csv`（L100-L192） | 七类文件具有相似的 skip/guard/write 调度结构 | 可进一步表驱动，但需保留每类前置条件 |
| `write_element_currents_csv`（L227-L282） | 多字段连续安全访问 | 可抽取 NamedTuple 快照 |
| `write_cohesive_damage_csv`（L368-L377） | cohesive 快照有多组长度守卫 | 可抽取安全向量读取 helper |
| `write_cohesive_driving_force_csv`（L450-L547） | driving-force 输出包含多层数据存在性与索引检查 | 参数接通后再单独简化 |
