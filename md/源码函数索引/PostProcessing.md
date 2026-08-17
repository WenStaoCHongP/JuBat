# PostProcessing.jl

- **源文件**: `src/PostProcessing.jl`
- **行数**: 117 行
- **函数/struct 计数**: 1 个函数（无独立 struct）
- **职责**: 通用求解结果提取与去归一化还原；不包含循环汇总、CSV 写出或 CZM 损伤管理。
- **相关技术文档**: `md/10_参数传递与模块架构.md`、`md/13_耦合验证方案.md`

## 数据结构

本文件无独立 struct。返回值为以稳定字符串键组织的结果 `Dict`。

## 函数清单

### `PostProcessing(case, variables, v) -> Dict` — L1-L117

从无量纲中央 `variables` 提取前 `v` 个时间步，并恢复物理单位。

- 通用量（L3-L12）：时间、电压、电流、温度、颗粒应力和位移。
- 电化学量（L13-L47）：按 SPM、SPMe、P2D/sP2D 恢复浓度、电位、电流密度等。
- 热量（L49-L98）：恢复 lumped 或 distributed2D 温度、分层热源和单元结果。
- CZM 量（L99-L115）：恢复损伤、位移、牵引与分离。
- 结果字符串键、单位和数组切片语义保持不变。

## 职责边界

- 循环阶段/周期汇总、SOH、终止与绘图：`CyclePostProcess.jl`
- 在线循环快照采集：`CycleData.jl`
- CSV 序列化：`CsvExport.jl`
- CZM 损伤统计与 variables 映射：`CzmPostProcess.jl`

## 省略项

无。

### [DEBUG]

无。

### [PLACEHOLDER]

无。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|---|---|---|
| L13-L47 | 单函数内按四类电化学模型恢复不同字段 | 后续可按模型抽取字段恢复 helper，但必须固定所有结果键 |
| L85-L90 | 九个无量纲热单元键使用字符串数组循环直传 | 后续可提升为共享常量，避免与变量定义漂移 |
