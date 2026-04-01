# Thermal.jl

## 文件状态
修改 (modified)

## main分支
- 行数: 49
- 主要函数列表:
  - `ThermalLumped(case, variables)` -- 集总热模型，计算质量矩阵和净热源

## Parameters_Design分支
- 行数: 80 (+31 行, +63%)
- 主要函数列表:
  - `ThermalLumped(case, variables)` -- 扩展版，新增 SPMe 专用欧姆热计算和内部产热输出

## 变更详情

### 新增函数
无新增函数。

### 修改函数

#### `ThermalLumped(case::Case, variables)`

**变更内容：**

1. **SPMe 从 SPM 分支拆出**：
   - main: `if case.opt.model == "SPM" || case.opt.model == "SPMe"` 合并处理
   - HEAD: SPM 和 SPMe 各自独立分支，SPMe 分支包含更完整的欧姆热计算

2. **新增 SPMe 专用欧姆热计算** (行 15-44)：
   - 引入电解液浓度的高斯点值：`ce_n_gs`, `ce_p_gs`, `ce_sp_gs`
   - 引入区域网格：`mesh_ne`, `mesh_pe`, `mesh_sp`
   - 计算温度相关的有效电解液电导率：`kappa_ne_gs`, `kappa_pe_gs`, `kappa_sp_gs`
   - 通过 `IntV` 积分求平均电导率：`kappa_ne_av`, `kappa_pe_av`, `kappa_sp_av`
   - 计算固体相欧姆电势降 `dphi_S`：
     ```julia
     dphi_S = I_app / 3 * (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig)
     ```
   - 计算电解液相电阻 `R_EL`：
     ```julia
     R_EL = param.NE.thickness / (3 * kappa_ne_av) + param.SP.thickness / kappa_sp_av + param.PE.thickness / (3 * kappa_pe_av)
     ```
   - 计算电解液电势降 `dphi_e`（含浓差电势项）：
     ```julia
     dphi_e = 2 * T * (1 - tplus) * (csp_av - csn_av) / ce0 - I_app * R_EL - dphi_S
     ```
   - 计算欧姆热 `Q_ohm = -I_app * (dphi_S + dphi_e)`

3. **新增内部产热输出** (行 75-76)：
   - 计算总内部产热 `Q_in = Q_rxn + Q_ohm + Q_rev`
   - 写入变量字典 `variables["thermal lumped internal heat"] = [Q_in]`
   - 这使得内部产热可以在外部被读取（之前只有净热源 = 内部产热 - 散热）

4. **返回值微调** (行 78)：
   - main: `return MT, Q_rxn + Q_ohm + Q_rev - q`（内联计算）
   - HEAD: `return MT, Q_in - q`（使用预计算的 `Q_in`，语义相同）

### 删除函数
无删除函数。

## 依赖关系

### 该文件依赖哪些其他文件
- `src/Option.jl` -- 通过 `case.opt.model` 区分模型类型
- `src/SetCase.jl` -- `Case` 类型，`case.param` 热物性参数
- `src/SetMesh.jl` -- 网格对象 (`case.mesh["negative electrode"]` 等)
- `src/Assemble.jl` -- `IntV` 积分函数（SPMe 分支新增依赖）
- `src/ElectrolyteDiffusion.jl` (间接) -- 电解液浓度变量

### 哪些文件依赖该文件
- `src/JuBat.jl` -- `include("Thermal.jl")`, 导出 `ThermalLumped`
- `src/Solve.jl` -- 在集总热模型模式下调用 `ThermalLumped(case, variables)`

### 新增的外部依赖
无新增外部包依赖。

## 耦合分析

### 与 multi-SPMe + distributed2D + CZM 耦合的关系
- Thermal.jl 实现的是**集总热模型** (lumped)，与 distributed2D 是平行的、互斥的关系。
- 在 multi-SPMe 架构中，如果 `thermalmodel == "lumped"`，则使用此文件；如果 `thermalmodel == "distributed2D"`，则使用 ThermalDistributed.jl。
- SPMe 分支的欧姆热完善化（电解液相欧姆热）提高了电-热耦合的精度，但不涉及逐单元架构。

### 哪些变更是耦合相关的
- SPMe 欧姆热计算完善（电解液相电导率、浓差电势） -- 提高电-热耦合精度
- `thermal lumped internal heat` 输出 -- 为外部耦合分析（如与 distributed2D 结果对比）提供数据
- 这两项变更虽然不直接参与 distributed2D 路径，但服务于整体耦合架构的验证和对比

### 哪些变更是独立的
- 代码重构（SPMe 从 SPM 拆出独立分支）属于代码组织改进
- 返回值使用 `Q_in` 变量替代内联表达式属于代码可读性改进
