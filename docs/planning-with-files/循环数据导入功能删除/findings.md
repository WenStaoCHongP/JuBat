# CycleData 外部循环数据导入求解功能删除调研

## Scope Findings

- `CycleData.jl` 定义 2 个数据结构和 4 个函数：`solve_phase_with_export`、`solve_cycle_with_export`、`export_cycle_data_to_csv`、`load_cycle_data_from_csv`。
- 真正读取外部循环 CSV 的源码入口只有 `load_cycle_data_from_csv`；它返回松散 `Dict`，并非 `CycleExportData`，无法作为对称序列化 API 使用。
- 唯一实际调用点是 `example/循环验证/czm_from_precomputed_example.jl`，该示例加载预计算温度/SOC 后独立驱动 CZM。
- `solve_phase_with_export` 另被 `example/电化学-热耦合验证/不同倍率温度曲线.jl` 用于在线生成温度曲线，因此不能仅因名称带 export 就归为“导入外部数据求解”。
- `JuBat.jl` 顶层单独导出 `load_cycle_data_from_csv`；删除实现时必须同步移除该 export。
- `md/源码函数索引/CycleData.md` 仍完整记录 loader，且 `md/源码函数索引/JuBat.md` 也列出该顶层导出，需要同步更新。
- `example/循环验证/czm_from_precomputed_example.jl` 的全部职责都是读取 CSV 后跳过电化学-热模型直接驱动 CZM；删除 loader 后该文件没有可保留的可运行主路径，因此整体删除。
- `example/循环验证/export_cycle_data_example.jl` 仍可用于外部分析/归档，但其末尾原本引导用户运行已删除的预计算 CZM 示例；改为明确 CSV 不再被 JuBat 读取来驱动求解。
- `docs/文件对比/*.md` 是历史分支对比快照，不作为当前源码索引；其中旧 API 记录保留历史语义，不纳入当前索引一致性门禁。

## Deletion Decision

- 删除：`load_cycle_data_from_csv`、对应顶层 export、预计算 CSV 驱动 CZM 示例及当前源码索引条目。
- 保留：`TimeStepData`、`CycleExportData`、`solve_phase_with_export`、`solve_cycle_with_export`、`export_cycle_data_to_csv`。
