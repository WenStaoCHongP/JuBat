# Jellyrollmodel.jl 代码精简建议

## 文件概况
- **总行数**: 694 行
- **函数数量**: 13 个主要函数
- **核心功能**: 果冻卷几何、网格生成、边界识别、极耳处理

---

## 1. 重复逻辑（高优先级）

### 🔴 问题1：`material_at` 中的重复代码（~30行重复）

**位置**: Lines 113-139

**问题描述**:
```julia
# :spiral 模式
acc = 0.0
for (name, w) in p.order
    if s <= acc + w
        return name, s - acc
    end
    acc += w
end

# :rings 模式（几乎相同）
acc = 0.0
for (name, w) in p.order
    if x <= acc + w
        return name, x - acc
    end
    acc += w
end
```

**精简方案**:
```julia
"""提取为辅助函数"""
function _find_layer_in_period(offset::Float64, order)
    acc = 0.0
    for (name, w) in order
        if offset <= acc + w
            return name, offset - acc
        end
        acc += w
    end
    # 边界回退
    return order[end][1], offset - (sum(w for (_, w) in order) - order[end][2])
end

function material_at(r::Real, θ::Real, p; logic::Symbol=:spiral)
    r <= p.Rin && return :inner, 0.0
    r >= p.Rout && return :outer, 0.0
    
    offset = if logic === :spiral
        mod(r - (p.a + p.b*θ), p.t_repeat)
    elseif logic === :rings
        mod(r - p.Rin, p.t_repeat)
    else
        error("Unknown logic=$(logic). Use :spiral or :rings")
    end
    
    return _find_layer_in_period(offset, p.order)
end
```

**收益**: 
- 消除 20+ 行重复代码
- 提高可读性
- 便于维护

---

### 🔴 问题2：`edge_boundary` 中的重复检查（~15行）

**位置**: Lines 310-333, 341-345

**问题描述**:
```julia
# 重复检查 which 三次
if which === :inner
    θ_min, θ_max = 0.0, 2.0*pi
elseif which === :outer
    ...
else
    error("which must be :inner or :outer")
end

# 又重复检查
if which === :inner
    θ_cum = (r - p.a) / bval
elseif which === :outer
    θ_cum = (r - p.a - p.t_repeat) / bval
else
    error("which must be :inner or :outer")
end

# 再次重复检查
if which === :inner
    r_theo = p.a + p.b * θ_cum
else
    r_theo = p.a + p.b * θ_cum + p.t_repeat
end
```

**精简方案**:
```julia
function edge_boundary(mesh, nidx::Int, param_dim; 
                       which::Symbol=:inner, 
                       theta_range::Union{Tuple{Float64,Float64},Nothing}=nothing, 
                       tol::Float64=1e-4)
    # 验证输入并获取偏移
    offset = if which === :inner
        0.0
    elseif which === :outer
        p.t_repeat
    else
        error("which must be :inner or :outer")
    end
    
    # 获取节点坐标和参数
    x, y = mesh.node[nidx, 1], mesh.node[nidx, 2]
    p = jellyroll_spiral_params(param_dim)
    bval = max(p.b, 1e-12)
    r = hypot(x, y)
    
    # 获取 θ 范围（只需一次判断）
    if theta_range === nothing
        θ_min, θ_max = which === :inner ? 
            (0.0, 2.0*pi) : 
            (2.0*pi*(Int(p.n_wind)-1), 2.0*pi*Int(p.n_wind))
    else
        θ_min, θ_max = theta_range
    end
    
    # 计算累计角度和理论位置（使用 offset 统一处理）
    θ_cum = (r - p.a - offset) / bval
    θ_cum < θ_min || θ_cum > θ_max && return false
    
    # 验证距离
    r_theo = p.a + p.b * θ_cum + offset
    dist = hypot(x - r_theo*cos(θ_cum), y - r_theo*sin(θ_cum))
    return dist <= tol
end
```

**收益**: 
- 消除 3 次重复的条件判断
- 代码从 ~55 行减少到 ~35 行
- 逻辑更清晰

---

### 🟡 问题3：`jellyroll_tab_node_indices` 中正负极耳逻辑重复（~100行）

**位置**: Lines 635-686

**问题描述**: 正极耳和负极耳的处理逻辑几乎完全相同，只是参数不同。

**精简方案**:
```julia
"""辅助函数：通用极耳节点识别"""
function _find_tab_nodes(mesh, p, tab_angles, θ_cum_nodes, θ_cum_range, 
                         delta_theta_fn, tw, Rin, Rout; reverse_range=false)
    idx = Int[]
    θc_min, θc_max = θ_cum_range
    twoπ = 2.0*pi
    
    for θ0_orig in tab_angles
        θ0 = Float64(θ0_orig)
        # 归一化到节点覆盖范围
        while θ0 > θc_max; θ0 -= twoπ; end
        while θ0 < θc_min; θ0 += twoπ; end
        
        Δθ = delta_theta_fn(θ0, tw)
        θ_start, θ_end = reverse_range ? (θ0 - Δθ, θ0) : (θ0, θ0 + Δθ)
        
        for i in 1:length(θ_cum_nodes)
            r = hypot(mesh.node[i,1], mesh.node[i,2])
            θ_cum = θ_cum_nodes[i]
            if (Rin - 1e-8 <= r <= Rout + 1e-8) && 
               (θ_start <= θ_cum <= θ_end)
                push!(idx, i)
            end
        end
    end
    return unique(idx)
end

function jellyroll_tab_node_indices(mesh, param_dim)
    @assert mesh.dimension == 2 "jellyroll_tab_node_indices 仅适用于 2D 网格"
    p = jellyroll_spiral_params(param_dim)
    tab = param_dim.tab
    
    # 预计算所有节点的累计角度
    nn = size(mesh.node, 1)
    θ_cum_in = [(hypot(mesh.node[i,1], mesh.node[i,2]) - p.a) / p.b for i in 1:nn]
    θ_cum_out = [(hypot(mesh.node[i,1], mesh.node[i,2]) - p.a - p.t_repeat) / p.b for i in 1:nn]
    
    tw = get(tab, :width, 0.0)
    delta_theta_fn = (θ, w) -> _delta_theta_from_width(p.a, p.b, θ, w)
    
    # 正极耳（内螺旋）
    pos_idx = _find_tab_nodes(mesh, p, tab.theta_pos, θ_cum_in, 
                               (minimum(θ_cum_in), maximum(θ_cum_in)),
                               delta_theta_fn, tw, p.Rin, p.Rout)
    
    # 负极耳（外螺旋）
    neg_idx = _find_tab_nodes(mesh, p, tab.theta_neg, θ_cum_out,
                               (minimum(θ_cum_out), maximum(θ_cum_out)),
                               delta_theta_fn, tw, p.Rin, p.Rout; reverse_range=true)
    
    return pos_idx, neg_idx
end

"""从弧长计算角度增量（提取为独立函数）"""
function _delta_theta_from_width(a, b, θ0, width)
    width <= 0.0 || b <= 0.0 && return 0.0
    u0 = a + b * θ0
    F(u) = (u * sqrt(u^2 + b^2) + b^2 * asinh(u / b)) / (2.0 * b)
    # ... 二分查找逻辑 ...
end
```

**收益**:
- 消除 ~50 行重复代码
- 函数从 140 行减少到 ~70 行
- 更易测试和维护

---

## 2. 冗余计算（中优先级）

### 🟡 问题4：`jellyroll_spiral_params` 中的重复计算

**位置**: Lines 56-97

**问题描述**:
```julia
# widths 和 fracs 重复定义层次结构
widths = (PCC = ..., PE = ..., SP = ..., NE = ..., NCC = ...)
fracs = (PCC = widths.PCC/t_repeat, ...)
order = [(:PCC, widths.PCC), (:PE, widths.PE), ...]
```

**精简方案**:
```julia
function jellyroll_spiral_params(param_dim)
    cell = param_dim.cell
    Rin, Rout = cell.Rin, cell.Rout
    
    # 统一定义层结构（避免重复）
    layers = [:PCC, :PE, :SP, :NE, :NCC]
    widths = NamedTuple{Tuple(layers)}(
        getfield(getfield(param_dim, layer), :thickness) for layer in layers
    )
    lambdas = NamedTuple{Tuple(layers)}(
        getfield(getfield(param_dim, layer), :lambda) for layer in layers
    )
    
    t_repeat = sum(widths)
    fracs = map(w -> w/t_repeat, widths)
    order = collect(zip(layers, widths))
    boundaries = cumsum([0.0; collect(widths)])
    
    # 有效热导率（向量化）
    λ_r_eff = 1.0 / sum(fracs ./ lambdas)  # 调和平均
    λ_t_eff = sum(fracs .* lambdas)        # 算术平均
    
    return (; 
        a = Rin, 
        b = t_repeat/(2π),
        n_wind = Int(floor((Rout - Rin)/t_repeat)),
        t_repeat, fracs, names=layers, widths, order, boundaries,
        λ_r_eff, λ_t_eff, Rin, Rout
    )
end
```

**收益**:
- 消除重复定义
- 代码从 ~50 行减少到 ~30 行
- 易于扩展新层

---

### 🟡 问题5：Q4 形函数重复定义

**位置**: Lines 374-379

**问题描述**: 
`jellyroll_element_layer_weights` 中每次调用都重新定义形函数。

**精简方案**:
```julia
# 在文件顶部定义为常量
"""Q4 等参元形函数"""
const N_Q4 = (xi, eta) -> 0.25 .* [
    (1 - xi)*(1 - eta),
    (1 + xi)*(1 - eta),
    (1 + xi)*(1 + eta),
    (1 - xi)*(1 + eta)
]

function jellyroll_element_layer_weights(mesh, param_dim; 
                                         nsamples_per_dim::Int=4, 
                                         logic::Symbol=:spiral)
    # 直接使用常量 N_Q4
    # ...
end
```

**收益**:
- 避免每次调用重新定义
- 代码更清晰

---

## 3. 可删除内容（低优先级）

### 🟢 问题6：未使用的导入

**位置**: Line 8

```julia
using Plots  # ❌ 文件中未使用
```

**精简方案**: 删除该行

**收益**: 减少依赖，加快加载速度

---

### 🟢 问题7：未使用的坐标变换函数

**位置**: Lines 20-32

**问题描述**: 
`cart2pol` 和 `pol2cart` 功能简单，且Julia有类似标准函数。

**建议**: 
- 保留 `cart2pol`（内部多处使用）
- 删除 `pol2cart`（未使用，且 `pol2cart(r,θ) = (r*cos(θ), r*sin(θ))` 太简单）

**收益**: 减少 8 行代码

---

## 4. 代码组织优化（低优先级）

### 🟢 问题8：长函数可拆分

**需要拆分的函数**:
1. `jellyroll_tab_node_indices` (140行) → 拆分为 3-4 个子函数
2. `jellyroll_collector_seed_mesh` (90行) → 拆分为网格生成和权重计算

**建议结构**:
```julia
# 主函数
function jellyroll_collector_seed_mesh(param_dim; kwargs...)
    nodes = _generate_spiral_nodes(param_dim, kwargs)
    elements = _create_band_elements(nodes)
    mesh = _assemble_mesh(nodes, elements, kwargs)
    _cache_layer_weights!(mesh, param_dim)
    return mesh
end

# 子函数
function _generate_spiral_nodes(p, nθ, phase) ... end
function _create_band_elements(nodes) ... end
function _assemble_mesh(nodes, elements, gsorder) ... end
function _cache_layer_weights!(mesh, p) ... end
```

**收益**:
- 提高可读性
- 易于测试
- 便于复用

---

### 🟢 问题9：可合并的判断逻辑

**位置**: Lines 168-189 (`jellyroll_Q4_mesh`)

**当前代码**:
```julia
if crop_mode === :center
    # ... 处理 :center
elseif crop_mode === :inscribed
    # ... 处理 :inscribed
else
    error(...)
end
```

**优化方案**:
```julia
crop_fn = if crop_mode === :center
    (base, e) -> begin
        center = jellyroll_element_centers(base)[e,:]
        r = hypot(center...)
        Rin <= r <= Rout
    end
elseif crop_mode === :inscribed
    (base, e) -> begin
        nodes = base.element[e, :]
        rs = hypot.(base.node[nodes, 1], base.node[nodes, 2])
        minimum(rs) >= Rin && maximum(rs) <= Rout
    end
else
    error("Unknown crop_mode=$(crop_mode)")
end

keep = [e for e in 1:ne if crop_fn(base, e)]
```

**收益**: 更函数式，易于扩展新模式

---

## 5. 性能优化（可选）

### 🟢 问题10：向量化机会

**位置**: Lines 197-207

**当前代码**:
```julia
function jellyroll_element_centers(mesh)
    centers = zeros(Float64, ne, 2)
    for e in 1:ne
        nodes = mesh.element[e, :]
        xy = mesh.node[nodes, :]
        centers[e, 1] = mean(xy[:,1])
        centers[e, 2] = mean(xy[:,2])
    end
    return centers
end
```

**优化方案**:
```julia
function jellyroll_element_centers(mesh)
    ne = size(mesh.element, 1)
    # 向量化计算（Q4 固定4个节点）
    return [mean(mesh.node[mesh.element[e, :], d]) 
            for e in 1:ne, d in 1:2]
end
```

**收益**: Julia 向量化通常更快

---

## 精简总结

### 预期成果

| 项目 | 精简前 | 精简后 | 减少 |
|------|--------|--------|------|
| 总行数 | 694 | ~520 | 25% |
| 主函数数 | 13 | 10 | 3个 |
| 重复代码 | ~120行 | 0 | 100% |

### 优先级路线图

#### 🔴 高优先级（立即实施）
1. 精简 `material_at` 重复逻辑
2. 精简 `edge_boundary` 重复检查
3. 精简 `jellyroll_tab_node_indices` 重复代码

**预估工作量**: 1-2 小时  
**预估减少**: 100+ 行

#### 🟡 中优先级（本周）
4. 优化 `jellyroll_spiral_params` 计算
5. 提取 Q4 形函数为常量

**预估工作量**: 30分钟  
**预估减少**: 30+ 行

#### 🟢 低优先级（可选）
6. 删除未使用导入和函数
7. 拆分长函数
8. 向量化优化

**预估工作量**: 2-3 小时  
**预估减少**: 50+ 行，提升可读性

---

## 实施建议

### 第一阶段：消除重复（1-2小时）
1. 提取 `_find_layer_in_period` 辅助函数
2. 重构 `edge_boundary` 使用统一 offset
3. 提取 `_find_tab_nodes` 通用函数

### 第二阶段：优化计算（30分钟）
4. 向量化 `jellyroll_spiral_params`
5. 提取 Q4 形函数常量

### 第三阶段：代码组织（可选，2-3小时）
6. 删除未使用代码
7. 拆分长函数
8. 添加单元测试

---

## 预期收益

### 代码质量
- ✅ 消除 ~120 行重复代码（17%）
- ✅ 提高可维护性（统一逻辑）
- ✅ 改善可测试性（函数更小）

### 性能
- ✅ 轻微性能提升（减少重复计算）
- ✅ 更快的加载时间（删除未使用导入）

### 可读性
- ✅ 函数更短更专注
- ✅ 逻辑更清晰
- ✅ 更易理解和修改

---

**编写日期**: 2025-11-19  
**文件版本**: Jellyrollmodel.jl (694行)  
**建议类型**: 代码精简与优化
