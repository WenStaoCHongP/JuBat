# parameters/Enertech.jl

## 文件状态: 修改 (M)

## main分支
- 行数: 174
- 内容: Enertech 电池参数集
- 结构:
  - PE (正极) 参数: 完整电化学 + 力学
  - NE (负极) 参数: 完整电化学 + 力学
  - EL (电解液) 参数: De, kappa, tplus
  - SP (隔膜) 参数
  - cell (电池级) 参数
  - NCC/PCC (集流体) 参数
  - tab (极耳) 参数
  - binder (粘结剂) 参数
  - scale (归一化尺度)
  - `param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)` - **9个参数组**

## Parameters_Design分支
- 行数: 175 (+1)

## 变更详情

### 新增内容
1. **新增 `cohesive = Cohesive()`** (1行)
   - 创建默认空的内聚力模型参数

2. **修改 `param_dim` 构造** (末尾行)
   ```julia
   # main:
   param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)
   # Parameters_Design:
   param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale, cohesive)
   ```
   - 新增第 10 个参数组 `cohesive`

### 未修改
- 所有物理参数值完全不变
- 两个分支均无文件末尾换行符

## 依赖关系

### 依赖
- 依赖 `SetParams.jl` 中定义的所有结构体
- 通过 `ChooseCell()` 中的 `include` 加载

### 被依赖
- 被 `ChooseCell("Enertech")` 调用
- 产出的 `param_dim` 被整个框架使用

## 耦合分析

**直接耦合到 multi-SPMe+distributed2D+CZM**: 间接（通过 `Cohesive` 结构）

与 LGM50.jl 和 Northrop.jl 完全相同的适配性变更模式。Enertech 不是 Jellyroll 电池，不包含 Jellyroll 特有参数，CZM 参数为零值默认。

变更性质: 纯适配性变更（兼容 `Params` 结构新增字段），不涉及物理参数修改。
