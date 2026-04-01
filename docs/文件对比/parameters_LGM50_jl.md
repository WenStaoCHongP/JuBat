# parameters/LGM50.jl

## 文件状态: 修改 (M)

## main分支
- 行数: 133
- 内容: LG M50 电池参数集
- 结构:
  - PE (正极) 参数: 电化学 + 力学 (E, nu, Omega)
  - NE (负极) 参数: 电化学 + 力学 (E, nu, Omega)
  - EL (电解液) 参数: De, kappa, tplus
  - SP (隔膜) 参数
  - cell (电池级) 参数
  - NCC/PCC (集流体) 参数
  - tab (极耳) 参数
  - binder (粘结剂) 参数
  - scale (归一化尺度)
  - `param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)` - **9个参数组**

## Parameters_Design分支
- 行数: 136 (+3)

## 变更详情

### 新增内容
1. **新增 `cohesive = Cohesive()`** (1行)
   - 创建默认空的内聚力模型参数（所有字段使用默认零值）

2. **修改 `param_dim` 构造** (1行)
   ```julia
   # main:
   param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)
   # Parameters_Design:
   param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale, cohesive)
   ```
   - 新增第 10 个参数组 `cohesive`

3. **新增空行** (2处)
   - PE 和 NE 参数块之间各添加了一个空行（格式化）

### 未修改
- 所有物理参数值（电化学、力学、热学等）完全不变
- 函数/闭包定义不变
- 文件末尾无换行符（两个分支均如此）

## 依赖关系

### 依赖
- 依赖 `SetParams.jl` 中定义的所有结构体（`Electrode`, `Cell`, `Scale`, `Cohesive`, `Params` 等）
- 通过 `ChooseCell()` 中的 `include` 加载

### 被依赖
- 被 `ChooseCell("LG M50")` 调用
- 产出的 `param_dim` 被整个框架使用

## 耦合分析

**直接耦合到 multi-SPMe+distributed2D+CZM**: 间接（通过 `Cohesive` 结构）

此参数文件的变更非常小，仅新增了默认空的 `Cohesive` 参数。LG M50 不是 Jellyroll 电池，因此不包含 Jellyroll 特有的几何/热学参数。CZM 参数全部为零值默认，意味着 LG M50 参数集默认不启用 CZM。

变更性质: 适配性变更（兼容 `Params` 结构新增的 `cohesive` 字段），不涉及物理参数修改。
