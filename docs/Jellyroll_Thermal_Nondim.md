# Jellyroll 热归一化理论（2D 分布式热模型）

本文整理并统一 Jellyroll 电池在本仓库中的热学无量纲化方案（Scheme B），给出 SI → 无量纲的参考尺度、控制方程、边界条件与热源缩放，并标注与代码实现的对应关系，便于验证与扩展。

## 1. 模型范围与耦合

- 热模型：二维平面热传导（电芯横截面，Q4 单元）
  - SI 方程： (ρ c) ∂T/∂t = ∇·(K ∇T) + q    （q 为体积内热源，W/m³）
  - 外圈对流、内孔默认绝热；可选极耳边界处理
- 电化学耦合：SPMe 供给反应热、可逆热与欧姆热；并联系统决定各热单元电流分配 I_e
- 维度匹配：1D 电化学 → 2D 热 通过“分层平均 + 几何映射（层权重）”闭合

主要代码位置：
- 无量纲标尺：`src/SetParams.jl` 中 `ChooseCell`（SI）与 `NormaliseParam`（电化学 nd）
- 热装配与边界：`src/ThermalDistributed.jl`（`ThermalDistributed2D`、`ThermalDistributed2D_BC`、`heatQ_Source`）
- 示例：`example/spme_thermal2d_example.jl`

## 2. 参考尺度（Scheme B）

由 `ChooseCell` 计算并存入 `param_dim.scale`：

- L_th：特征热长度（优先用 Rout），`scale.L_th`
- k_ref：参考导热率（优先取电极导热率之一），`scale.k_th`
- (ρc)_ref：参考体积热容（ρ·c），`scale.rho_c_th`
- T_ref：温度基准（通常 298 K），`scale.T_ref`
- q_ref：参考体热源 k_ref·T_ref/L_th²，`scale.q_th`
- t_th：热扩散时间 (ρc)_ref·L_th²/k_ref，`scale.t_th`
- Bi：对流 Biot 数 h·L_th/k_ref，`scale.h_th`

注：电化学仍采用其自身的 t0=3600 s、ϕ_ref=RT_ref/F 等标尺；本文仅针对热场无量纲化。

## 3. 无量纲变量与算子

- 温度：T* = T / T_ref
- 时间：t* = t / t_th
- 坐标：x* = x / L_th，dΩ* = dΩ / L_th²，dL* = dL / L_th
- 材料：
  - K* = K / k_ref（各向异性时径向/切向先聚合再旋转）
  - (ρc)* = (ρc) / (ρc)_ref
- 内热源：q* = q / q_ref

无量纲控制方程：

  (ρc)* ∂T*/∂t* = ∇*·( K* ∇*T* ) + q*

离散装配（`ThermalDistributed2D`）：
- MT = ∬ (ρc)* NᵀN dΩ*
- KT = ∬ Bᵀ K* B dΩ*
- FT = ∬ q* N dΩ*

实现要点：
- 质量项与右端的面积权重统一使用 dΩ* = dΩ/L_th²（代码：用 `wJ./(L_th^2)` 缩放）
- K* 可按“层权重 f_k(e)”在局部径/切向上聚合（串联调和平均/并联算术平均），再按单元极角旋转成全局分量（Kxx/Kxy/Kyy）

## 4. 边界条件（对流）

SI：−n·(K ∇T) = h (T − T_amb)

无量纲：− n·(K* ∇*T*) = Bi (T* − T_amb*)

离散边界项（`ThermalDistributed2D_BC`）：
- 刚度贡献：KT += ∫_{Γ_out*} Bi NᵀN dL*
- 载荷贡献：FT += ∫_{Γ_out*} Bi T_amb* N dL*

实现细节：
- 自动识别外边界（出现一次的边，且半径接近 r_max）
- 线积分统一用 dL* = dL / L_th
- 内边界默认绝热（可扩展对流）
- 极耳（可选）：通过网格定位节点并对这些节点施加惩罚法近似 Dirichlet，使 T 节点等于随时间变化的极耳温度 T_tab(t)

## 5. 内热源 q 的缩放与映射

电化学侧（SPMe）给出层平均热源：
- 反应热：Q_rxn = a_s j η
- 可逆热：Q_rev = a_s j T (∂U/∂T)
- 欧姆热（简化一致流近似）：
  - 固相：q_ohm,s ≈ I_app² / (3 σ_eff)
  - 液相：q_ohm,e ≈ I_app² / (3 κ_eff)（隔膜为 I_app² / κ_sp）

层合成到单元：
- 若单元层权重为 f_NE, f_SP, f_PE, f_PCC, f_NCC，则
  q_e ≈ (f_NE+f_SP+f_PE)·Q_ele + f_PCC·Q_PCC + f_NCC·Q_NCC  （均为体热源，W/m³）

单位路径：
- 维度化路径：`variables["heat_source_fields"] = q_e`，并置 `heat_source_units_code=1.0`
- 无量纲路径：`variables["heat_source_fields"] = q_e / q_ref`，并置 `heat_source_units_code=0.0`

装配时（`ThermalDistributed2D`）：如标记为 SI，则自动除以 q_ref 转为 q* 后进入 FT。

## 6. 时间步与两套时标的关系

- 电化学步长使用 t0=3600 s 的无量纲 dt
- 热步长使用 t_th 的无量纲 dt_th
- 两者换算（示例 `spme_thermal2d_example.jl`）：

  dt_th = dt · (t0 / t_th)

这样热方程矩阵方程写成：

  [(MT/dt_th) + KT] T^{n+1} = (MT/dt_th) T^n + FT + 边界项

## 7. 快速校核清单（量纲一致性）

- 面积、长度：所有积分权重中 dΩ* 与 dL* 的缩放是否一致使用 L_th
- 热源：`heat_source_units_code` 与 `q_ref` 的除/乘是否在唯一位置完成，不重复
- K*、(ρc)*：无量纲化是否均采用 k_ref 与 (ρc)_ref
- Bi：边界项中是否统一使用 `scale.h_th = h·L_th/k_ref`
- 时间：是否使用 dt_th（非电化学 dt）进入热步矩阵

## 8. 与代码的映射

- 参考尺度：`src/SetParams.jl` 中 `ChooseCell` 末尾设置
  - `scale.L_th, k_th, rho_c_th, q_th, t_th, h_th`
- 热装配：`src/ThermalDistributed.jl/ThermalDistributed2D`
  - MT：`ρc_weights .= ρc_e[ele] .* (wJ./(L_th^2))`
  - KT：各向异性旋转或各向同性回退，均是 K* 进入装配
  - FT：若 SI，先 `/ q_ref` 得 q* 再用 `(wJ./(L_th^2))` 装配
- 边界：`ThermalDistributed2D_BC`
  - Bi = `scale.h_th`，`T_amb* = T_amb / T_ref`，线权重用 `1/L_th`
- 热源生成：`heatQ_Source`
  - 计算 Q_rxn/Q_rev/Q_ohm、层合成 q_e，由 `opt.units_thermal` 决定 SI/nd 路径

## 9. 小结

本方案将 Jellyroll 的 2D 热问题严格以 Scheme B 无量纲化：使用 L_th、k_ref、(ρc)_ref、T_ref，使控制方程、边界与热源在数值装配中具有一致的尺度。与电化学（SPMe）耦合时，通过“层平均热源 + 几何层权重”映射到 2D 单元，使 1D-2D 维度闭合；时间刻度通过 `dt_th = dt · (t0/t_th)` 对齐两套时标。配合 `heat_source_units_code` 的单一转换点，避免单位重复换算，便于在 SI/nd 两路径间切换与验证。
