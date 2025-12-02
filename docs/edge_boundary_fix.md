# edge_boundary 终点角度判定修复

## 问题描述

`edge_boundary` 函数在判定边界节点时，其默认的θ范围计算逻辑与 `jellyroll_collector_seed_mesh` 网格生成时的终点判定逻辑存在差异，导致外边界节点识别数量不足。

## 根本原因

### 修改前的逻辑

**网格生成** (`jellyroll_collector_seed_mesh`)：
```julia
s_in  = 0.0
s_out = p.t_repeat
θ0 = max(0.0, (Rin - a - s_in) / b)      # 内螺旋起点约束
θ1 = min((Rout - a - s_out) / b, (Rout - a) / b)  # 外螺旋终点约束
```
- 内外螺旋共享同一个θ范围 `[θ0, θ1]`
- θ0 由内螺旋到达 Rin 的角度决定
- θ1 由外螺旋到达 Rout 和内螺旋到达 Rout 的角度中较小者决定

**边界判定** (`edge_boundary` 修改前)：
```julia
θ_end = (p.Rout - p.a - p.t_repeat) / bval
θ_min, θ_max = θ_end-2π, θ_end
```
- 使用固定的 `θ_max - 2π` 作为起点
- 没有考虑网格实际的起点约束 `max(0.0, (Rin - a - s_in) / b)`
- 没有使用 `min()` 函数来处理多个约束条件

### 问题影响

1. **起点不一致**：edge_boundary 使用 `θ_end - 2π` 作为起点，而网格实际起点是 `max(0.0, (Rin - a) / b)`
2. **终点判定不完整**：edge_boundary 只使用了外螺旋约束，没有考虑 `min(..., (Rout - a) / b)` 的第二个约束
3. **导致节点遗漏**：当节点的 θ_cum 在网格范围内但在 edge_boundary 旧范围外时，会被错误地排除

## 修复方案

### 修改后的逻辑

```julia
# 与 jellyroll_collector_seed_mesh 保持完全一致
if theta_range === nothing
    # 网格生成时内外螺旋共享同一个θ范围 [θ0, θ1]：
    s_in = 0.0
    s_out = p.t_repeat
    θ_start = max(0.0, (p.Rin - p.a - s_in) / bval)
    θ_end = min((p.Rout - p.a - s_out) / bval, (p.Rout - p.a) / bval)
    θ_min, θ_max = θ_start, θ_end 
else
    θ_min, θ_max = theta_range
end
```

### 关键改进

1. **起点一致**：θ_start = `max(0.0, (Rin - a) / b)`，与网格生成的 θ0 完全一致
2. **终点一致**：θ_end = `min((Rout - a - t_repeat) / b, (Rout - a) / b)`，与网格生成的 θ1 完全一致
3. **内外螺旋统一**：无论 `which=:inner` 还是 `which=:outer`，都使用相同的θ范围（与网格生成逻辑一致）

## 验证逻辑

### 数学推导

假设：
- a = Rin（螺旋参数）
- b = t_repeat / (2π)（螺距）
- Rout - Rin = N * t_repeat（N圈螺旋）

则：
```
θ0 = max(0.0, (Rin - Rin - 0) / b) = max(0.0, 0) = 0

θ1 = min((Rout - Rin - t_repeat) / b, (Rout - Rin) / b)
   = min((N*t_repeat - t_repeat) / (t_repeat/2π), N*t_repeat / (t_repeat/2π))
   = min((N-1)*2π, N*2π)
   = (N-1)*2π
```

所以网格覆盖 θ ∈ [0, 2π(N-1)]，即从第0圈到第N-1圈。

### 节点验证

对于网格上的任意节点：
- 坐标：(x, y)
- 半径：r = √(x² + y²)
- 反推角度：θ_cum = (r - a - offset) / b
  - 内螺旋：offset = 0
  - 外螺旋：offset = t_repeat

修改后的 edge_boundary 保证：
- 当 θ0 ≤ θ_cum ≤ θ1 时，节点被识别为边界节点
- 这与网格生成的θ范围完全一致，不会遗漏任何节点

## 兼容性说明

### 对现有代码的影响

1. **使用默认参数的调用**：行为改变
   ```julia
   # 修改前：识别 θ ∈ [θ_end-2π, θ_end] 范围的节点
   # 修改后：识别 θ ∈ [θ0, θ1] 范围的节点（与网格一致）
   is_outer = edge_boundary(mesh, i, param_dim; which=:outer)
   ```

2. **显式指定 theta_range 的调用**：行为不变
   ```julia
   # 用户显式指定θ范围，不受默认值修改影响
   is_outer = edge_boundary(mesh, i, param_dim; which=:outer, theta_range=(2π*(N-1), 2π*N))
   ```

### 建议

- 对于需要识别特定圈层的场景，应显式传入 `theta_range` 参数
- 对于需要识别网格所有边界节点的场景，使用默认参数即可

## 测试建议

1. **基本验证**：
   ```julia
   mesh = jellyroll_collector_seed_mesh(param_dim; nθ=160)
   
   # 检查内螺旋节点数
   inner_count = count(i -> edge_boundary(mesh, i, param_dim; which=:inner), 1:mesh.nlen)
   
   # 检查外螺旋节点数
   outer_count = count(i -> edge_boundary(mesh, i, param_dim; which=:outer), 1:mesh.nlen)
   
   # 预期：inner_count ≈ nθ+1 = 161，outer_count ≈ nθ+1 = 161
   ```

2. **拓扑一致性验证**：
   ```julia
   # 网格拓扑：节点1-161为内圈，节点162-322为外圈
   inner_nodes_topo = 1:(nθ+1)
   outer_nodes_topo = (nθ+2):(2*(nθ+1))
   
   # 边界识别
   inner_nodes_edge = findall(i -> edge_boundary(mesh, i, param_dim; which=:inner), 1:mesh.nlen)
   outer_nodes_edge = findall(i -> edge_boundary(mesh, i, param_dim; which=:outer), 1:mesh.nlen)
   
   # 验证匹配度
   @assert inner_nodes_topo ⊆ inner_nodes_edge  # 拓扑内圈应全部被识别
   @assert outer_nodes_topo ⊆ outer_nodes_edge  # 拓扑外圈应全部被识别
   ```

3. **距离容差验证**：
   ```julia
   # 检查识别到的边界节点到理论螺旋线的距离
   for i in outer_nodes_edge
       x, y = mesh.node[i, :]
       r = hypot(x, y)
       θ_cum = (r - p.a - p.t_repeat) / p.b
       r_theo = p.a + p.b * θ_cum + p.t_repeat
       dist = abs(r - r_theo)
       @assert dist < 1e-4  # 应小于默认容差
   end
   ```

## 修改文件

- `src/Jellyrollmodel.jl`：修改 `edge_boundary` 函数的默认θ范围计算逻辑（第301-316行）

## 修改日期

2025-12-02
