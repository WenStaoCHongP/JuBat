# P2D.jl

## 文件状态: 修改 (M)

## main分支
- 行数: 292
- 函数列表:
  - `P2D(case, yt, t; jacobi)` - P2D 模型主函数
  - `P2D_mass_BC(case, variables)` - P2D 质量边界条件
  - `P2D_charge_BC(case, variables)` - P2D 电荷边界条件
  - `P2D_potentials(case, yt, t, K_pot, variables)` - P2D 电势迭代求解
  - `P2D_variables(case, yt, t)` - P2D 变量计算

## Parameters_Design分支
- 行数: 287 (-5, -2%)

## 变更详情

### 修改函数

#### `P2D_potentials()`
**静水应力计算简化**:
- 删除了 `stress_rn_surf_gs = 0` 和 `stress_rp_surf_gs = 0`（硬编码零径向应力）
- 删除了完整的静水应力计算:
  ```julia
  # main 分支 (已删除):
  hydrostatic_stress_n_gs = (1/3) * (2 * stress_theta_n_surf_gs .+ stress_rn_surf_gs)
  hydrostatic_stress_p_gs = (1/3) * (2 * stress_theta_p_surf_gs .+ stress_rp_surf_gs)
  ```
- 将过电势 eta 计算从静水应力改为直接使用切向应力:
  ```julia
  # main 分支:
  eta_n_gs_rel = ... - hydrostatic_stress_n_gs .* case.param.NE.Omega
  # Parameters_Design 分支:
  eta_n_gs_rel = ... - (2/3) * stress_theta_n_surf_gs * case.param.NE.Omega
  ```

**分析**: 在 main 分支中，`stress_rn_surf_gs = 0`，因此:
- `hydrostatic_stress = (1/3) * (2 * stress_theta + 0) = (2/3) * stress_theta`

所以 Parameters_Design 分支的 `(2/3) * stress_theta_n_surf_gs * Omega` 与 main 分支的 `hydrostatic_stress_n_gs * Omega` 在数学上**完全等价**。变更目的是简化代码，去掉不必要的中间变量。

**其他变更**:
- 移除了一些尾随空格
- 去掉了文件末尾的空行

#### `P2D_variables()`
- 移除了一行末尾空格（第 258 行 `eta_n_gs` 行）

### 删除的代码（5行）
```julia
-    stress_rn_surf_gs = 0
-    stress_rp_surf_gs = 0
-    hydrostatic_stress_n_gs = (1/3) * (2 * stress_theta_n_surf_gs .+ stress_rn_surf_gs)
-    hydrostatic_stress_p_gs = (1/3) * (2 * stress_theta_p_surf_gs .+ stress_rp_surf_gs)
```

## 依赖关系

### 无变更
依赖关系与 main 分支完全一致。

## 耦合分析

**直接耦合到 multi-SPMe+distributed2D+CZM**: 否（无实质逻辑变更）

P2D 模型的变更是纯粹的代码简化（消除冗余中间变量），数学行为完全不变。P2D 模型在当前 Jellyroll 框架中主要作为参考/验证模型使用，multi-SPMe 架构基于 SPMe 而非 P2D。

此文件的变更独立于 multi-SPMe + distributed2D + CZM 的开发。
