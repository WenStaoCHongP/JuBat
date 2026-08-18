# StateAccess.jl

- **源文件**: `src/StateAccess.jl`
- **行数**: 18 行
- **函数/struct 计数**: 1 个独立函数；0 个 struct 定义
- **职责**: 电化学模型统一取温入口——由热学配置（`opt.thermalmodel`）而非字典键的偶然存在性决定温度来源。
- **相关技术文档**: `md/04_电化学模型_SPMe.md` §6.3、`md/10_参数传递与模块架构.md`
- **新增背景**: b4c0cde 模块化重构新增文件（2026-08-17）。

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `representative_temperature(case, state; supplied=nothing)` — L8-L18

返回电化学模型应当使用的温度：

- `supplied !== nothing` → 直接返回（`SPMe_element` 经此注入单元温度 `T_e`，实现逐单元温度反馈）
- `opt.thermalmodel == "none"` → 返回 `case.param.cell.T0`（恒温）
- `"lumped"` / `"distributed2D"` → 返回 `only(state[case.index["temperature"]])`
- 其他值 → 抛 `ArgumentError`（无对应的热-电化学温度状态规则）

**调用点**: `src/SPM.jl:59`、`src/SPMe.jl:128,231`、`src/P2D.jl:216`。

## 跨文件依赖

- `SetCase.jl`：`Case`（`case.opt`、`case.param`、`case.index`）

## 省略项

无。全部 function 均有独立条目。

### [DEBUG]

无。

### [PLACEHOLDER]

无。

### [COMPLEX-CHECK]

无。
