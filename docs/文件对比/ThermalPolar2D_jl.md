# ThermalPolar2D.jl

## 文件状态
新增 (new file)

## main分支
- 文件不存在于 main 分支

## Parameters_Design分支
- 行数: 120
- 主要函数列表:
  1. `ThermalPolar2D_Ring(case, variables, mesh_data)` -- 极坐标 FVM 环形热模型组装

## 变更详情

### 新增函数

#### 1. `ThermalPolar2D_Ring(case::Case, variables::Dict{String,Any}, mesh_data)` (行 1-120)
- **功能**：使用有限体积法 (FVM) 在极坐标 (r, theta) 下组装环形热模型的质量矩阵、刚度矩阵和载荷向量
- **函数签名**：
  ```julia
  function ThermalPolar2D_Ring(case::Case, variables::Dict{String,Any}, mesh_data)
  # 返回: MT (稀疏对角), KT (稀疏), F (向量)
  ```
- **实现要点**：
  1. **网格结构**：使用 `mesh_data.r` (径向节点) 和 `mesh_data.ntheta` (周向节点数) 定义极坐标网格
  2. **热物性归一化**：
     - 体积热容 `(rho_c)* = C* / V*` = `param.cell.heat_Q / param.cell.volume`
     - 各向异性导热率 `k_r` (径向) 和 `k_t` (切向) 已在 `param.cell` 中归一化
  3. **FVM 离散**：
     - 径向：使用半节点间距 `r_imh`, `r_iph` 确定控制体边界
     - 周向：周期性边界条件（theta 方向首尾相连）
     - 控制体面积 `vol = 0.5 * (r_iph^2 - r_imh^2) * dtheta`
  4. **径向传导**：
     - 内邻居传导率 `a_rm = k_r * r_imh * dtheta / dr_im`
     - 外邻居传导率 `a_rp = k_r * r_iph * dtheta / dr_ip`
  5. **切向传导**：
     - 传导率 `a_t = k_t * area_theta / (r_i * dtheta)`，周期性处理
  6. **外边界对流**：
     - 在最外层径向节点施加 Newton 对流 BC：`Bi * area_bc * (T - T_amb)`
  7. **热源映射**：将 `heat_source_fields` 从单元平均值映射到节点
  8. **矩阵格式**：返回稀疏矩阵 `MT` (对角), `KT` (稀疏 COO -> CSR)

### 修改函数
不适用（新文件）。

### 删除函数
不适用（新文件）。

## 依赖关系

### 该文件依赖哪些其他文件
- `src/Option.jl` -- 通过 `case.opt` 获取模型配置（隐式，通过 Case 类型）
- `src/SetCase.jl` -- `Case` 类型定义，`case.param`, `case.mesh`
- 外部包：`SparseArrays` -- `sparse()`, `spdiagm()` 用于稀疏矩阵构造

### 哪些文件依赖该文件
- `src/JuBat.jl` -- 导出 `ThermalPolar2D_Ring`
- `src/Solve.jl` -- 在 `thermalmodel == "ring2D_polar"` 条件下调用 `ThermalPolar2D_Ring(case, vars, mesh_data)`

### 新增的外部依赖
- `SparseArrays` (Julia 标准库，可能已在 JuBat.jl 中引入)

## 耦合分析

### 与 multi-SPMe + distributed2D + CZM 耦合的关系
- **替代热求解器**：`ThermalPolar2D_Ring` 提供了与 `ThermalDistributed2D_Ring` (FEM) 不同的数值方法 (FVM) 来求解同一个环形热模型。
- 通过 `variables["heat_source_fields"]` 接收来自 SPMe 的热源（电-热耦合输入）。
- 计算得到的温度场通过 `MT * dT/dt + KT * T = F` 传递给主求解器，再反馈给 SPMe 模型（热-电耦合输出）。
- 不直接涉及 CZM 耦合（边界条件仅包含对流，不含界面热阻）。

### 哪些变更是耦合相关的
- 热源映射逻辑（行 ~40-50）：接收 SPMe 产热 -- 电-热耦合输入
- 对流边界条件（行 ~85-92）：外边界冷却 -- 热边界耦合
- 全部函数均为耦合服务

### 哪些变更是独立的
无独立变更（全部为耦合服务）。

### 与 ThermalDistributed.jl 的关系
- `ThermalDistributed2D_Ring` (FEM) 和 `ThermalPolar2D_Ring` (FVM) 是同一物理模型的两种数值实现
- FVM 版本更适合极坐标几何，天然保持守恒性
- FEM 版本可以复用已有的矩阵组装工具 (`Assemble`)
- 两者共用 `ThermalRing2D_BC` 的边界条件概念，但实现不同（FVM 将对流直接嵌入矩阵组装）
