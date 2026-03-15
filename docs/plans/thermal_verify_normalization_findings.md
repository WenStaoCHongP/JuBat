# 热验证脚本归一化审查发现

## 审查范围

- `example/热模块验证/thermal_verify.jl`
- `example/热模块验证/thermal_error_source_analysis.jl`
- `example/热模块验证/thermal_equivalent_lumped_compare.jl`

## 审查标准

基于 `md/01_参数定义与归一化.md` 中定义的统一归一化方案：

### 正确的归一化模式
1. **热源归一化**：`q_nd = q_phys / scale.q`，其中 `scale.q = P_ref / L³`
2. **时间输入**：`opt.time` 和 `opt.dt` 使用物理单位（秒）
3. **时间处理**：Solve.jl 自动归一化（当 `model == "thermal"` 时）
4. **update_fn 时间参数**：接收的是无量纲时间 `t_nd`
5. **结果时间**：输出的 `result.time` 是物理时间（已还原）
6. **边界条件**：使用物理参数 `param_dim.cell.h`、`param_dim.cell.T_amb`
7. **几何积分**：网格坐标已归一化，`wJ` 是无量纲值

### 错误的归一化模式（需要识别）
1. 使用旧的热尺度（如 `ρ_th`, `c_th`, `k_th`）
2. 手动进行时间归一化（应由 Solve.jl 自动处理）
3. 在 update_fn 中将时间参数当作物理时间使用
4. 对边界条件参数进行不必要的归一化
5. 对几何积分权重 `wJ` 除以 `L²`（已经是无量纲）

## 发现记录

### thermal_verify.jl

#### ✓ 正确使用的部分
1. **热源归一化**（Line 334, 398, 502）：
   ```julia
   q_ref = scale.q  # 统一能量尺度热源参考 (P_ref / L^3)
   ```
   ✓ 使用正确的 `scale.q`

2. **时间输入**（Line 328-329, 391-392, 495-496）：
   ```julia
   opt.time = [0.0, 3600]
   opt.dt = [1.0, 10]
   ```
   ✓ 使用物理时间单位

3. **热源计算**（Line 344, 408, 512）：
   ```julia
   q0 = 2.0e5
   q0_nd = q0 / q_ref
   ```
   ✓ 正确归一化

#### ⚠️ 需要验证的部分
1. **update_fn 时间参数**（Line 355-361, 419-425, 523-529）：
   ```julia
   update_fn = (t, vars) -> begin
       if model == "ring2D_polar"
           vars["heat_source_fields"] = fill(q0_nd, ne)
       else
           vars["heat_source_fields"] = compute_q_elem(mesh, q_func, t) ./ q_ref
       end
   end
   ```
   - `q_func` 定义为 `q_func = (r, theta, t) -> q0`（不依赖时间）
   - ✓ 时间参数 `t` 未被使用，因此无影响
   - ✓ 热源归一化正确（除以 `q_ref`）

2. **几何积分**（Line 453-456）：
   ```julia
   Rin_nd = Rin / scale.L  # Normalized for mask calculation
   Rout_nd = Rout / scale.L
   ```
   - ✓ 仅用于绘图遮罩，不影响计算

3. **精确解对比**（Line 372-376, 440-444, 463-464, 547-548）：
   ```julia
   r_nodes_dim = hypot.(mesh.node[:, 1], mesh.node[:, 2]) .* scale.L
   T_exact = analytical_solution_ring.(r_nodes_dim, q0, k_r, h, Rin, Rout, T_f)
   ```
   - ✓ 网格坐标乘以 `scale.L` 还原为物理坐标
   - ✓ 使用物理参数计算精确解

#### 结论
**thermal_verify.jl 已正确使用新归一化方案，无需修改。**

---

### thermal_error_source_analysis.jl

#### ✓ 正确使用的部分
1. **热源归一化**（Line 155-156）：
   ```julia
   q_ref = case.param_dim.scale.q  # P_ref / L^3
   q_phys_hist = q_nd_hist .* q_ref
   ```
   ✓ 使用正确的 `scale.q` 还原物理热源

2. **时间输入**（Line 104-105）：
   ```julia
   opt.time = [0.0, 60.0]
   opt.dt = [0.5, 10.0]
   ```
   ✓ 使用物理时间单位

3. **体积积分**（Line 144-149）：
   ```julia
   A_elem = zeros(Float64, size(mesh.element, 1))
   @inbounds for g in 1:ngs
       e = mesh.gs.ele[g]
       A_elem[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
   end
   ```
   ✓ 直接使用 `wJ = weight * detJ`（已归一化）

4. **内热源功率**（Line 185）：
   ```julia
   P_internal[k] = sum(q_phys .* A_elem) * H
   ```
   - `q_phys` 是物理热源密度 [W/m³]
   - `A_elem` 是无量纲面积（需要乘以 `L²` 才是物理面积）
   - `H` 是物理高度 [m]
   - ⚠️ **潜在问题**：`A_elem` 是无量纲的，应该乘以 `L²`

#### ⚠️ 需要验证的部分
1. **面积计算的量纲**（Line 144-150）：
   ```julia
   A_elem[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
   ```
   - 根据文档（md/01_参数定义与归一化.md, section 3.4），网格坐标已归一化
   - `detJ` 是无量纲雅可比行列式
   - `A_elem` 是无量纲面积
   - **问题**：Line 185 中 `P_internal[k] = sum(q_phys .* A_elem) * H`
     - `q_phys` [W/m³] × `A_elem` [无量纲] × `H` [m] = ?
     - 应该是：`q_phys` [W/m³] × `A_elem * L²` [m²] × `H` [m] = [W]
   - ⚠️ **需要修正**：`P_internal[k] = sum(q_phys .* A_elem) * H * scale.L^2`

2. **边界散热功率**（Line 193-201）：
   ```julia
   for (a, b) in outer_edges
       xa, ya = mesh.node[a, 1], mesh.node[a, 2]
       xb, yb = mesh.node[b, 1], mesh.node[b, 2]
       L = hypot(xb - xa, yb - ya)
       T_edge = 0.5 * (T_nodes[a] + T_nodes[b])
       p_out += h * L * H * (T_edge - Tamb)
   end
   ```
   - `mesh.node` 是无量纲坐标
   - `L = hypot(xb - xa, yb - ya)` 是无量纲长度
   - `h` 是物理传热系数 [W/(m²·K)]
   - `H` 是物理高度 [m]
   - **问题**：`h` [W/(m²·K)] × `L` [无量纲] × `H` [m] × `ΔT` [K] = ?
   - 应该是：`h` [W/(m²·K)] × `L * scale.L` [m] × `H` [m] × `ΔT` [K] = [W]
   - ⚠️ **需要修正**：`p_out += h * L * scale.L * H * (T_edge - Tamb)`

3. **表面散热功率**（Line 205）：
   ```julia
   P_boundary_surface[k] = 2.0 * h * sum((T_elem .- Tamb) .* A_elem)
   ```
   - `h` [W/(m²·K)] × `A_elem` [无量纲] × `ΔT` [K] = ?
   - 应该是：`h` [W/(m²·K)] × `A_elem * scale.L²` [m²] × `ΔT` [K] = [W]
   - ⚠️ **需要修正**：`P_boundary_surface[k] = 2.0 * h * sum((T_elem .- Tamb) .* A_elem) * scale.L^2`

#### 结论
**thermal_error_source_analysis.jl 存在几何尺度问题，需要修正。**

**需要修正的地方**：
1. Line 185: `P_internal[k] = sum(q_phys .* A_elem) * H * scale.L^2`
2. Line 199: `p_out += h * L * scale.L * H * (T_edge - Tamb)`
3. Line 205: `P_boundary_surface[k] = 2.0 * h * sum((T_elem .- Tamb) .* A_elem) * scale.L^2`

---

### thermal_equivalent_lumped_compare.jl

#### ✓ 正确使用的部分
1. **热源归一化**（Line 102-103）：
   ```julia
   q_ref = case.param_dim.scale.q  # 统一能量尺度热源参考 (P_ref / L^3)
   q_phys_hist = q_nd_hist .* q_ref
   ```
   ✓ 使用正确的 `scale.q` 还原物理热源

2. **时间输入**（Line 63-64）：
   ```julia
   opt.time = [0.0, 60.0]
   opt.dt = [0.5, 10.0]
   ```
   ✓ 使用物理时间单位

3. **体积积分**（Line 95-99）：
   ```julia
   A_elem = zeros(Float64, size(mesh.element, 1))
   @inbounds for g in 1:length(mesh.gs.detJ)
       e = mesh.gs.ele[g]
       A_elem[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
   end
   ```
   ✓ 直接使用 `wJ = weight * detJ`（已归一化）

#### ⚠️ 需要验证的部分
1. **内热源功率**（Line 133）：
   ```julia
   P_internal[k] = sum(q_phys_hist[:, k] .* A_elem) * H
   ```
   - 与 thermal_error_source_analysis.jl 相同的问题
   - ⚠️ **需要修正**：`P_internal[k] = sum(q_phys_hist[:, k] .* A_elem) * H * scale.L^2`

2. **分项热源功率**（Line 134-136）：
   ```julia
   P_rxn[k] = sum((q_rxn_ne[:, k] .+ q_rxn_pe[:, k]) .* A_elem) * H
   P_reversible[k] = sum((q_rev_ne[:, k] .+ q_rev_pe[:, k]) .* A_elem) * H
   P_ohmic[k] = sum((q_ohm_s_ne[:, k] .+ ...) .* A_elem) * H
   ```
   - 与 Line 133 相同的问题
   - ⚠️ **需要修正**：所有这些行都需要乘以 `scale.L^2`

3. **边界散热功率**（Line 137）：
   ```julia
   P_boundary_eq[k] = h * A_cool * (T_vol[k] - Tamb)
   ```
   - `A_cool = case.param_dim.cell.cooling_surface` 是物理面积 [m²]
   - ✓ 这里使用的是物理面积，正确

#### 结论
**thermal_equivalent_lumped_compare.jl 存在几何尺度问题，需要修正。**

**需要修正的地方**：
1. Line 133: `P_internal[k] = sum(q_phys_hist[:, k] .* A_elem) * H * scale.L^2`
2. Line 134: `P_rxn[k] = sum((q_rxn_ne[:, k] .+ q_rxn_pe[:, k]) .* A_elem) * H * scale.L^2`
3. Line 135: `P_reversible[k] = sum((q_rev_ne[:, k] .+ q_rev_pe[:, k]) .* A_elem) * H * scale.L^2`
4. Line 136: `P_ohmic[k] = sum((q_ohm_s_ne[:, k] .+ ...) .* A_elem) * H * scale.L^2`

---

## 问题分类

### 类型 1: 几何尺度问题（高优先级）
**影响文件**：
- thermal_error_source_analysis.jl
- thermal_equivalent_lumped_compare.jl

**问题描述**：
网格积分得到的 `A_elem` 是无量纲面积，在计算物理功率时需要乘以 `scale.L^2` 转换为物理面积。

**根本原因**：
根据文档（md/01_参数定义与归一化.md, section 3.3-3.4），网格坐标已归一化为 `x* = x / L`，因此：
- 无量纲面积元：`dA* = dx* dy* = (dx/L)(dy/L) = dA/L²`
- 物理面积：`dA = dA* × L²`

**修复策略**：
在所有涉及功率计算的地方，将 `A_elem` 乘以 `scale.L^2`。

### 类型 2: 边界积分尺度问题（高优先级）
**影响文件**：
- thermal_error_source_analysis.jl

**问题描述**：
边界边长 `L = hypot(xb - xa, yb - ya)` 是无量纲的，在计算边界散热功率时需要乘以 `scale.L` 转换为物理长度。

**根本原因**：
网格节点坐标是无量纲的，边长也是无量纲的。

**修复策略**：
在边界散热功率计算中，将边长 `L` 乘以 `scale.L`。

### 类型 3: 无问题（已正确）
**影响文件**：
- thermal_verify.jl

**说明**：
该文件已正确使用新归一化方案，无需修改。

---

## 修复优先级

### 高优先级（影响结果正确性）
1. ✅ thermal_error_source_analysis.jl - 几何尺度修正
2. ✅ thermal_equivalent_lumped_compare.jl - 几何尺度修正

### 低优先级（代码改进）
1. 添加注释说明归一化约定
2. 统一代码风格

---

## 下一步行动

1. 实施 thermal_error_source_analysis.jl 的修复
2. 实施 thermal_equivalent_lumped_compare.jl 的修复
3. 运行修复后的脚本验证结果
4. 更新文档注释
