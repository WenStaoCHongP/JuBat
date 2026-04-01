# Mechanical.jl

## 文件状态: 修改

## main分支
- 行数: 149
- 主要函数列表:
  - `Mechanical(case, variables)` — 颗粒级扩散应力计算（SPM/SPMe/P2D三种模型）
  - `Calstressdisp(electrode, mesh, cs, T)` — 球形颗粒扩散应力计算（径向中心应力、表面切向应力、表面位移、应力耦合扩散系数）

## Parameters_Design分支
- 行数: 354
- 主要函数列表:
  - `Mechanicaloutput(case, variables)` — 重命名的原 `Mechanical` 函数
  - `Calstressdisp(electrode, mesh, cs, T)` — 球形颗粒扩散应力计算（小幅修改）
  - `thermal_diffusion_stress_2D(case, variables)` — **新增**：宏观2D有限元热-扩散应力计算

## 变更详情

### 新增函数

#### `thermal_diffusion_stress_2D(case, variables)` (L165-L354)
- **功能**：计算宏观层面的热和扩散应力，使用2D平面应力有限元方法
- **输入**：case（含热网格）、variables（含温度场、SOC分布）
- **输出**：更新variables字典，添加6个应力/位移场：
  - `"diffusion stress xx"`, `"diffusion stress yy"`, `"diffusion stress xy"` — 各应力分量
  - `"diffusion stress vonMises"` — Von Mises等效应力
  - `"displacement x"`, `"displacement y"` — 位移场
  - `"thermal stress vonMises"` — 热应力分量
  - `"diffusion stress vonMises only"` — 扩散应力分量

- **算法**：
  1. 提取温度场和SOC分布（ne个单元级别）
  2. 计算有效材料参数（E_eff, nu_eff, alpha_eff, beta_n, beta_p），厚度加权平均
  3. 组装力学刚度矩阵K_mech（使用全局Assemble函数）
  4. 计算热-扩散载荷向量F_mech：epsilon_0 = alpha*deltaT + beta_n*deltaSoc_n + beta_p*deltaSoc_p
  5. 施加边界条件（内外壁fixed_xy，罚函数法）
  6. 求解位移场 U_M = K_mech \ F_mech
  7. 恢复应力场（通过q4_center_gradients计算单元应变，减去初始应变得到弹性应变）
  8. 分离热应力和扩散应力贡献（按应变比例分配）
  9. 结果转换为有量纲（乘以L_ref）

### 修改函数

#### `Mechanical` -> `Mechanicaloutput` (L1)
- **变更**：函数名从 `Mechanical` 改为 `Mechanicaloutput`
- **原因**：原函数名与模块文件名冲突，且新函数名更准确反映其功能（输出应力相关的变量更新）

#### `Mechanicaloutput` 中静水应力计算简化 (L11-L19, L74-L108)
- **变更**：删除 `stress_rn_surf = 0`、`stress_rp_surf = 0` 及显式 `hydrostatic_stress` 计算
- **原代码**：
  ```
  hydrostatic_stress_n = (1/3) * (2 * stress_theta_n_surf .+ stress_rn_surf)  # stress_rn_surf = 0
  eta_p_new = eta_p - hydrostatic_stress_p * param.PE.Omega
  ```
- **新代码**：
  ```
  eta_p_new = eta_p - (2/3) * stress_theta_p_surf * param.PE.Omega
  ```
- **分析**：由于 `stress_rn_surf = 0`，`hydrostatic_stress = (1/3) * 2 * stress_theta_surf = (2/3) * stress_theta_surf`，新代码是数学等价简化，但更清晰地表达了"球形颗粒表面径向应力为零时，静水应力等于2/3切向应力"这一物理事实。

#### `Calstressdisp` 中 `pi` 规范化 (L141)
- **变更**：`pi` -> `π`（Unicode pi符号）
- **原代码**：`cs_av = (3 /(4 * pi * (rs .^ 3))) * ...`
- **新代码**：`cs_av = (3 /(4 * π * (rs ^ 3))) * ...`
- **附加变更**：`rs .^ 3` 改为 `rs ^ 3`（rs是标量，无需广播），更准确

### 删除函数

无删除。所有原main分支的函数均保留（仅改名）。

## 依赖关系

### main分支依赖
- `src/SetParams.jl` — Electrode参数
- `src/SetMesh.jl` — Mesh、PickElement
- `src/Variables.jl` — Case类型

### Parameters_Design分支新增依赖
- `src/Jellyrollmodel.jl` — `identify_boundary_nodes`（用于thermal_diffusion_stress_2D的边界条件）
- `src/Assemble.jl`（或等效） — `Assemble`、`Assemble1D`、`element_nodal_mean`、`q4_center_gradients`全局函数

### 调用该文件的文件
- `src/SPMe.jl` — 调用 `Mechanicaloutput(case, variables)`（L4, L47）
- `src/SPM.jl` — 调用 `Mechanicaloutput(case, variables)`（L4）
- `src/P2D.jl` — 调用 `Mechanicaloutput(case, variables)`（L4, L193）
- `src/JuBat.jl` — `include("mechanical.jl")`（L28），`export thermal_diffusion_stress_2D`（L52）

## 耦合分析

### main分支
- 仅提供颗粒级扩散应力计算，通过修改过电位（eta）和电压（V_cell）影响电化学动力学
- 力学-电化学单向耦合：应力计算结果反馈到Butler-Volmer方程

### Parameters_Design分支新增
- **宏观2D应力场**：`thermal_diffusion_stress_2D` 是全新的宏观力学分析功能
  - 与distributed2D热模型紧密耦合：使用相同的热网格和温度场
  - 与multi-SPMe耦合：使用各单元SOC分布计算扩散应变
  - 提供Von Mises等效应力场，可用于评估电池结构完整性
  - 分离热应力和扩散应力的贡献，便于分析各物理场的力学影响

- **与CZM的关系**：`thermal_diffusion_stress_2D` 的应力场可以作为CZM载荷的参考值。当前CZM使用独立的热-化学载荷计算（`assemble_thermal_chemical_load`），但两者使用相同的物理原理（epsilon_0 = alpha*deltaT + beta*deltaSoc）。

- **设计意图**：`thermal_diffusion_stress_2D` 作为独立后处理工具函数使用（被export），不在主求解循环中调用。它提供应力场的可视化输出，而CZM求解器（czm.jl）内部有自己的应力计算流程。
