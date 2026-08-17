# CycleData.jl

- **源文件**: `src/CycleData.jl`
- **行数**: 410 行
- **函数/struct 计数**: 2 个 struct + 2 个函数 = 4 条目
- **职责**: 循环数据在线采集——定义时间步/单循环快照容器，并执行带数据采集的阶段和单循环求解；不承担 CSV 序列化。
- **相关技术文档**: `md/09_分流求解器.md`、`md/10_参数传递与模块架构.md`

## 数据结构

### `struct TimeStepData` — L7-L18

单时间步快照，保存时间、阶段、电压、电流、节点温度、温度统计和正负极 SOC。

### `struct CycleExportData` — L23-L30

单循环导出数据容器，保存时间步序列、节点坐标、单元连接以及网格规模。

## 函数清单

### `solve_phase_with_export(case, phase_type, t_max, I_current, V_limit, initial_state; ...)` — L32-L258

执行带在线快照采集的单阶段求解，返回 `(PhaseResult, Vector{TimeStepData})`。

- 支持普通模型与逐单元 SPMe 状态。
- 使用现有时间离散和自适应步长。
- 记录电压、温度、SOC 与容量。
- 按电压截止或阶段时长终止。

### `solve_cycle_with_export(case, cycle_opt; verbose, export_interval)` — L268-L410

按放电、静置、充电、静置顺序聚合阶段快照，返回 `CycleExportData`。

## 职责边界

- `export_cycle_data_to_csv` 已迁移至 `CsvExport.jl`；公开函数名和输出六文件格式保持不变。
- `CycleData.jl` 不读取外部 CSV，也不负责多循环 `CyclingResult` 汇总。

## 省略项

无。

### [DEBUG]

阶段进度打印属于用户可见运行状态，不是调试代码。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|---|---|---|
| L245-L247 | `PhaseResult` 的 CZM 损伤字段固定为 0 | 此在线采集入口不替代完整 CZM 循环耦合求解 |
| L294 | 初始状态中的 `"V" => 3.7` | 典型电压初值尚未参数化 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|---|---|---|
| L56-L82 | multi-SPMe/普通模型初态和温度场分支组装 | 可抽取状态初始化 helper |
| L165-L199 | SOC、温度滤波、面积均温与快照构造集中在时间循环 | 可抽取 `TimeStepData` 构造函数 |
| L300-L391 | 四阶段调用结构重复 | 可在行为夹具完善后抽取阶段执行 helper |
