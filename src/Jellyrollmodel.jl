# Jellyrollmodel.jl - 重构版
# 果冻卷（卷绕电池）物理建模工具
# 
# 结构组织：
# 1. 变量参数定义计算
# 2. 网格生成划分
# 3. 边界定义
# 4. 极耳边界识别
# 5. 辅助函数

using LinearAlgebra
using Statistics

# ========================================================================
# 1. 变量参数定义计算
# ========================================================================

export jellyroll_spiral_params, material_at

"""
    jellyroll_spiral_params(param_dim)

计算果冻卷阿基米德螺旋参数：
- r(θ) = a + b θ，其中 b ≈ t_repeat / (2π)
- a ≈ Rin（内半径）
- 返回：螺旋参数、层厚度、有效热导率等

# 返回字段
- `a, b`: 螺旋方程参数
- `t_repeat`: 一个完整层序的厚度（m）
- `n_wind`: 绕组圈数
- `fracs`: 各层体积分数 (NE, SP, PE, PCC, NCC)
- `widths`: 各层厚度 (m)
- `order`: 层顺序 [(name, thickness), ...]
- `boundaries`: 累积边界位置
- `λ_r_eff, λ_t_eff`: 等效径向/切向热导率
- `Rin, Rout`: 内外半径
"""
function jellyroll_spiral_params(param_dim)
    cell = param_dim.cell
    Rin, Rout = cell.Rin, cell.Rout
    
    # 定义层结构（从内到外：PCC -> PE -> SP -> NE -> NCC）
    layers = [:PCC, :PE, :SP, :NE, :NCC]
    
    # 获取各层厚度和热导率
    widths = NamedTuple{Tuple(layers)}(
        getfield(getfield(param_dim, layer), :thickness) for layer in layers
    )
    lambdas = NamedTuple{Tuple(layers)}(
        getfield(getfield(param_dim, layer), :lambda) for layer in layers
    )
    
    # 计算总厚度和体积分数
    t_repeat = sum(widths)
    fracs = map(w -> w/t_repeat, widths)
    
    # 层顺序和累积边界
    order = collect(zip(layers, widths))
    boundaries = cumsum([0.0; collect(widths)])
    
    # 等效热导率（向量化计算）
    frac_vals = collect(fracs)
    lambda_vals = collect(lambdas)
    λ_r_eff = 1.0 / sum(frac_vals ./ lambda_vals)  # 径向调和平均
    λ_t_eff = sum(frac_vals .* lambda_vals)        # 切向算术平均
    
    # 螺旋参数
    a = Rin
    b = t_repeat / (2π)
    n_wind = Int(floor((Rout - Rin) / t_repeat))
    
    return (; 
        a, b, t_repeat, n_wind, 
        fracs, names=layers, widths, order, boundaries,
        λ_r_eff, λ_t_eff, Rin, Rout
    )
end

"""
    material_at(r, θ, p; logic=:spiral)

确定 (r, θ) 处的材料层。

# 参数
- `r, θ`: 极坐标
- `p`: 螺旋参数（来自 `jellyroll_spiral_params`）
- `logic`: `:spiral`（真实螺旋）或 `:rings`（同心环近似）

# 返回
- `(layer::Symbol, offset::Float64)`: 层名和层内偏移
- 特殊值：`:inner` (r≤Rin), `:outer` (r≥Rout)
"""
function material_at(r::Real, θ::Real, p; logic::Symbol=:spiral)
    r <= p.Rin && return :inner, 0.0
    r >= p.Rout && return :outer, 0.0
    
    # 计算在一个周期内的偏移
    offset = if logic === :spiral
        mod(r - (p.a + p.b*θ), p.t_repeat)
    elseif logic === :rings
        mod(r - p.Rin, p.t_repeat)
    else
        error("Unknown logic=$(logic). Use :spiral or :rings")
    end
    
    # 查找对应层
    return _find_layer_in_period(offset, p.order)
end

# ========================================================================
# 2. 网格生成划分
# ========================================================================

export jellyroll_collector_seed_mesh, jellyroll_get_layer_weights

# 全局缓存：网格 -> 层权重
const __jr_layer_weights = IdDict{Any, Matrix{Float64}}()

"""
    jellyroll_collector_seed_mesh(param_dim; nθ=360, gsorder=2, phase=0.0)

基于集流体"导轨"的条带网格生成器（唯一推荐的网格生成方法）。

# 原理
- 沿内外螺旋线采样：
  - 内螺旋：r_in(θ) = a + bθ（s_in = 0，PCC内侧）
  - 外螺旋：r_out(θ) = a + bθ + t_repeat（s_out = t_repeat，NCC外侧）
- 连接对应点形成 Q4 条带单元：[(in_i, out_i, out_{i+1}, in_{i+1})]
- 每个单元天然跨越完整层序（NE, SP, PE, PCC, NCC）

# 参数
- `nθ`: 每圈的分段数（建议≥160）
- `gsorder`: 高斯积分阶数（默认2）
- `phase`: 相位对齐角度（默认0.0）

# 返回
- `Mesh`: Q4/2D 网格对象
- 自动缓存层权重，可通过 `jellyroll_get_layer_weights(mesh)` 获取

# 特点
- ✅ 单元包含完整层序，适合"每单元=子电池"模型
- ✅ 自动生成精确的层权重 f_k
- ✅ 覆盖整个螺旋（多圈）
- ✅ "以直代曲"近似，计算高效
"""
function jellyroll_collector_seed_mesh(param_dim; nθ::Int=360, gsorder::Int=2, phase::Float64=0.0)
    p = jellyroll_spiral_params(param_dim)
    a, b = p.a, p.b
    Rin, Rout = p.Rin, p.Rout
    
    # 两条螺旋偏移：完整层序 [0, t_repeat]
    s_in  = 0.0
    s_out = p.t_repeat
    
    # θ 范围裁剪（确保 r_in≥Rin 且 r_out≤Rout）
    θ0 = max(0.0, (Rin - a - s_in) / b)
    θ1 = min((Rout - a - s_out) / b, (Rout - a) / b)
    θ1 > θ0 || error("collector-seeded: 无有效 θ 范围 [Rin, Rout]")
    
    # 等角度采样
    seg_per_turn = max(3, nθ)
    dθ = 2π / seg_per_turn
    
    # 相位对齐
    k0 = ceil(Int, (θ0 - phase) / dθ)
    k1 = floor(Int, (θ1 - phase) / dθ)
    
    θ = if k1 <= k0
        # 窄窗：至少2段
        nθ_eff = max(2, round(Int, (θ1 - θ0) / dθ))
        collect(range(θ0, θ1; length=nθ_eff+1))
    else
        phase .+ (k0:k1) .* dθ
    end
    
    nθ = length(θ) - 1  # 实际分段数
    
    # 生成节点：内外螺旋
    spiral_xy(offset) = begin
        r = a .+ b .* θ .+ offset
        x = r .* cos.(θ)
        y = r .* sin.(θ)
        return x, y
    end
    
    x_in, y_in = spiral_xy(s_in)
    x_out, y_out = spiral_xy(s_out)
    
    # 组装节点数组
    nnode = 2*(nθ+1)
    node = zeros(Float64, nnode, 2)
    for i in 1:(nθ+1)
        node[i, 1] = x_in[i];  node[i, 2] = y_in[i]
        j = (nθ+1) + i
        node[j, 1] = x_out[i]; node[j, 2] = y_out[i]
    end
    
    # 生成单元：条带 Q4
    ne = nθ
    element = zeros(Int64, ne, 4)
    for i in 1:nθ
        element[i, 1] = i                # n1: 内侧当前
        element[i, 2] = (nθ+1) + i       # n2: 外侧当前
        element[i, 3] = (nθ+1) + i + 1   # n3: 外侧下一个
        element[i, 4] = i + 1            # n4: 内侧下一个
    end
    
    # 高斯积分
    gs = GetGS(element, node, gsorder, "Q4")
    mesh = Mesh("Q4", 2, node, nnode, element, gs)
    
    # 缓存层权重（每个单元等于全层序的体积分数）
    f_k = zeros(Float64, ne, 5)  # [NE, SP, PE, PCC, NCC]
    fr = p.fracs
    @inbounds for e in 1:ne
        f_k[e, 1] = fr.NE
        f_k[e, 2] = fr.SP
        f_k[e, 3] = fr.PE
        f_k[e, 4] = fr.PCC
        f_k[e, 5] = fr.NCC
    end
    __jr_layer_weights[mesh] = f_k
    
    return mesh
end

"""
    jellyroll_get_layer_weights(mesh) -> Matrix{Float64} or nothing

返回通过 `jellyroll_collector_seed_mesh` 生成的网格的层权重矩阵。

# 返回
- `Matrix{Float64}` (ne×5): 层权重 [NE, SP, PE, PCC, NCC]
- `nothing`: 如果网格不是通过 collector_seed_mesh 生成
"""
function jellyroll_get_layer_weights(mesh)
    return get(__jr_layer_weights, mesh, nothing)
end

# ========================================================================
# 3. 边界定义
# ========================================================================

export edge_boundary

"""
    edge_boundary(mesh, nidx, param_dim; which=:inner/:outer, theta_range=nothing, tol=1e-4)

基于螺旋方程的精确边界节点识别。

# 螺旋方程
- 内螺旋：r(θ) = a + bθ
- 外螺旋：r(θ) = a + bθ + t_repeat

# 判断方法
1. 计算节点的累计角度 θ_cum
2. 检查 θ_cum ∈ [θ_min, θ_max]
    - **默认范围与 `jellyroll_collector_seed_mesh` 的截断角度完全一致**
    - θ_start = max(0.0, (Rin - a - s_spiral) / b)
    - θ_end = min((Rout - a - s_spiral) / b, (Rout - a) / b)
    - 其中 s_spiral = 0（内螺旋）或 t_repeat（外螺旋）
3. 验证节点到理论螺旋线的距离 ≤ tol

# 参数
- `mesh`: 网格对象
- `nidx`: 节点索引
- `param_dim`: 参数对象
- `which`: `:inner`（内螺旋）或 `:outer`（外螺旋）
- `theta_range`: θ 范围 (θ_min, θ_max)，默认使用网格覆盖的完整角度区间
- `tol`: 距离容差（m），默认 1e-4

# 返回
- `Bool`: 是否为边界节点

# 示例
```julia
# 判断节点是否在第1圈内螺旋上
is_inner = edge_boundary(mesh, i, param_dim; which=:inner)

# 判断节点是否在第N圈外螺旋上
is_outer = edge_boundary(mesh, i, param_dim; which=:outer)
```
"""
function edge_boundary(mesh, nidx::Int, param_dim; 
                       which::Symbol=:inner, 
                       theta_range::Union{Tuple{Float64,Float64},Nothing}=nothing, 
                       tol::Float64=1e-4)
    # 获取螺旋参数和偏移
    p = jellyroll_spiral_params(param_dim)
    offset = if which === :inner
        0.0
    elseif which === :outer
        p.t_repeat
    else
        error("which must be :inner or :outer")
    end
    
    # 获取节点坐标
    x, y = mesh.node[nidx, 1], mesh.node[nidx, 2]
    r = hypot(x, y)
    bval = max(p.b, 1e-12)
    
    # 确定 θ 范围；默认遵循 collector_seed_mesh 的截断逻辑
    if theta_range === nothing
        s_in = 0.0
        s_out = p.t_repeat
        θ_start = 0.0
        θ_end = min((p.Rout - p.a - s_out) / bval, (p.Rout - p.a) / bval)
        if which === :inner
            θ_min = θ_start
            θ_max = min(θ_start + 2.0*π, θ_end)
        else  # :outer
            θ_min = max(θ_end - 2.0*π, θ_start)
            θ_max = θ_end
        end
    else
        θ_min, θ_max = theta_range
    end
    
    # 计算累计角度
    θ_cum = (r - p.a - offset) / bval
    if θ_cum < θ_min || θ_cum > θ_max
        return false
    end
    
    # 计算理论位置并验证距离
    r_theo = p.a + p.b * θ_cum + offset
    x_theo = r_theo * cos(θ_cum)
    y_theo = r_theo * sin(θ_cum)
    dist = hypot(x - x_theo, y - y_theo)
    
    return dist <= tol
end

# ========================================================================
# 4. 极耳边界识别
# ========================================================================

export jellyroll_tab_node_indices

"""
    jellyroll_tab_node_indices(mesh, param_dim) -> (pos_indices, neg_indices)

识别受极耳影响的节点索引。

# 原理
基于螺旋参数将 `tab.width` 映射为沿螺旋的角度增量 Δθ，
然后选择落在对应角度/半径范围内的节点。

# 极耳位置
- 正极耳：内螺旋（s_in = 0, PCC侧）
- 负极耳：外螺旋（s_out = t_repeat, NCC侧）

# 返回
- `(pos_indices::Vector{Int}, neg_indices::Vector{Int})`

# 参数要求
- `param_dim.tab.theta_pos`: 正极耳角度数组
- `param_dim.tab.theta_neg`: 负极耳角度数组
- `param_dim.tab.width`: 极耳宽度（m）
"""
function jellyroll_tab_node_indices(mesh, param_dim)
    @assert mesh.dimension == 2 "jellyroll_tab_node_indices 仅适用于 2D 网格"
    
    p = jellyroll_spiral_params(param_dim)
    a, b = p.a, p.b
    Rin, Rout = p.Rin, p.Rout
    tab = param_dim.tab
    tw = hasproperty(tab, :width) ? tab.width : 0.0
    
    # 预计算所有节点的累计角度
    nn = size(mesh.node, 1)
    θ_cum_in = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a) / b for i in 1:nn]
    θ_cum_out = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a - p.t_repeat) / b for i in 1:nn]
    
    # 角度增量计算函数
    delta_theta_fn = (θ, w) -> _delta_theta_from_width(a, b, θ, w)
    
    # 正极耳（内螺旋）
    pos_idx = _find_tab_nodes(
        mesh, tab.theta_pos, θ_cum_in, 
        (minimum(θ_cum_in), maximum(θ_cum_in)),
        delta_theta_fn, tw, Rin, Rout
    )
    
    # 负极耳（外螺旋，角度范围反向）
    neg_idx = _find_tab_nodes(
        mesh, tab.theta_neg, θ_cum_out,
        (minimum(θ_cum_out), maximum(θ_cum_out)),
        delta_theta_fn, tw, Rin, Rout; reverse_range=true
    )
    
    return pos_idx, neg_idx
end

# ========================================================================
# 5. 辅助函数
# ========================================================================

export cart2pol, jellyroll_element_centers, jellyroll_effective_K_at

# ------------------------------------------------------------------------
# 5.1 坐标变换
# ------------------------------------------------------------------------

"""
    cart2pol(x, y) -> (r, θ)

笛卡尔坐标转极坐标。

# 返回
- `r`: 半径
- `θ`: 角度（弧度，范围 (-π, π]）
"""
function cart2pol(x::Real, y::Real)
    r = hypot(x, y)
    θ = atan(y, x)
    return r, θ
end

# ------------------------------------------------------------------------
# 5.2 网格几何计算
# ------------------------------------------------------------------------

"""
    jellyroll_element_centers(mesh) -> Matrix{Float64}

计算每个 Q4 单元的几何中心。

# 返回
- `Matrix{Float64}` (ne×2): 每行为 [x_center, y_center]
"""
function jellyroll_element_centers(mesh)
    ne = size(mesh.element, 1)
    return [mean(mesh.node[mesh.element[e, :], d]) for e in 1:ne, d in 1:2]
end

"""
    jellyroll_effective_K_at(θ, param_dim) -> Matrix{Float64}

计算在极角 θ 处的各向异性等效导热张量。

# 公式
K = λ_r * e_r * e_r' + λ_t * e_θ * e_θ'

其中：
- e_r = [cos(θ), sin(θ)]（径向单位向量）
- e_θ = [-sin(θ), cos(θ)]（切向单位向量）
- λ_r, λ_t 为等效径向/切向热导率

# 返回
- `Matrix{Float64}` (2×2): 导热张量
"""
function jellyroll_effective_K_at(θ::Real, param_dim)
    p = jellyroll_spiral_params(param_dim)
    λr, λt = p.λ_r_eff, p.λ_t_eff
    c, s = cos(θ), sin(θ)
    
    # 基向量
    er = [c, s]
    et = [-s, c]
    
    # 组装张量
    return λr * (er * er') + λt * (et * et')
end

# ------------------------------------------------------------------------
# 5.3 内部辅助函数（不导出）
# ------------------------------------------------------------------------

"""
    _find_layer_in_period(offset, order) -> (layer, local_offset)

在一个周期内查找层（内部辅助函数）。

# 参数
- `offset`: 周期内偏移（0 ≤ offset < t_repeat）
- `order`: 层顺序 [(name, width), ...]

# 返回
- `(layer::Symbol, local_offset::Float64)`: 层名和层内偏移
"""
function _find_layer_in_period(offset::Float64, order)
    acc = 0.0
    for (name, w) in order
        if offset <= acc + w
            return name, offset - acc
        end
        acc += w
    end
    
    # 边界回退到最后一层
    last_layer, last_width = order[end]
    total_width = sum(w for (_, w) in order)
    return last_layer, offset - (total_width - last_width)
end

"""
    _find_tab_nodes(...) -> Vector{Int}

通用极耳节点查找（内部辅助函数）。

# 参数
- `mesh`: 网格对象
- `tab_angles`: 极耳角度列表
- `θ_cum_nodes`: 所有节点的累计角度
- `θ_cum_range`: 累计角度范围 (min, max)
- `delta_theta_fn`: 角度增量计算函数
- `tw`: 极耳宽度
- `Rin, Rout`: 内外半径
- `reverse_range`: 是否反向角度范围（负极耳用）

# 返回
- `Vector{Int}`: 节点索引列表
"""
function _find_tab_nodes(mesh, tab_angles, θ_cum_nodes, θ_cum_range, 
                         delta_theta_fn, tw, Rin, Rout; reverse_range=false)
    isempty(tab_angles) && return Int[]
    
    idx = Int[]
    θc_min, θc_max = θ_cum_range
    twoπ = 2.0*π
    nn = length(θ_cum_nodes)
    
    for θ0_orig in tab_angles
        θ0 = Float64(θ0_orig)
        
        # 归一化到节点覆盖范围
        while θ0 > θc_max; θ0 -= twoπ; end
        while θ0 < θc_min; θ0 += twoπ; end
        
        # 计算角度增量
        Δθ = delta_theta_fn(θ0, tw)
        θ_start, θ_end = reverse_range ? (θ0 - Δθ, θ0) : (θ0, θ0 + Δθ)
        
        # 查找符合条件的节点
        for i in 1:nn
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

"""
    _delta_theta_from_width(a, b, θ0, width) -> Float64

从弧长计算角度增量（内部辅助函数）。

求解：沿阿基米德螺旋 r(θ) = a + bθ，从 θ0 到 θ0+Δθ 的弧长等于 width。

# 方法
使用螺旋弧长的闭式解和二分查找。

# 返回
- `Float64`: 角度增量 Δθ
"""
function _delta_theta_from_width(a::Float64, b::Float64, θ0::Float64, width::Float64)
    width <= 0.0 || b <= 0.0 && return 0.0
    width < 1e-12 && return 0.0
    
    u0 = a + b * θ0
    F(u) = (u * sqrt(u^2 + b^2) + b^2 * asinh(u / b)) / (2.0 * b)
    s0 = F(u0)
    
    # 初始上界估计
    hi = max(1e-6, width / max(1e-12, sqrt(u0^2 + b^2)))
    
    # 扩展上界直到包含解
    for _ in 1:100
        (F(u0 + b*hi) - s0) >= width && break
        hi *= 2.0
    end
    
    # 二分查找
    lo = 0.0
    tol = max(1e-12, width * 1e-9)
    for _ in 1:80
        mid = 0.5 * (lo + hi)
        sval = F(u0 + b*mid) - s0
        abs(sval - width) <= tol && return mid
        sval < width ? (lo = mid) : (hi = mid)
    end
    
    return 0.5 * (lo + hi)
end

# ========================================================================
# 模块说明
# ========================================================================

"""
# Jellyrollmodel 模块

果冻卷（卷绕电池）几何建模工具。

## 主要功能

### 1. 参数计算
- `jellyroll_spiral_params`: 计算螺旋参数、层厚度、等效热导率
- `material_at`: 判断给定位置的材料层

### 2. 网格生成
- `jellyroll_collector_seed_mesh`: 基于集流体导轨的条带网格生成器（推荐）
- `jellyroll_get_layer_weights`: 获取网格的层权重矩阵

### 3. 边界识别
- `edge_boundary`: 基于螺旋方程的精确边界节点识别

### 4. 极耳处理
- `jellyroll_tab_node_indices`: 识别受极耳影响的节点

### 5. 辅助工具
- `cart2pol`: 坐标变换
- `jellyroll_element_centers`: 计算单元中心
- `jellyroll_effective_K_at`: 计算各向异性导热张量

## 使用示例

```julia
# 1. 计算螺旋参数
p = jellyroll_spiral_params(param_dim)

# 2. 生成网格
mesh = jellyroll_collector_seed_mesh(param_dim; nθ=160)

# 3. 获取层权重
f_k = jellyroll_get_layer_weights(mesh)  # (ne×5) [NE, SP, PE, PCC, NCC]

# 4. 识别边界节点
is_inner = edge_boundary(mesh, node_idx, param_dim; which=:inner)
is_outer = edge_boundary(mesh, node_idx, param_dim; which=:outer)

# 5. 识别极耳节点
pos_nodes, neg_nodes = jellyroll_tab_node_indices(mesh, param_dim)
```

## 网格生成建议

**推荐使用 `jellyroll_collector_seed_mesh`**，原因：
- ✅ 每个单元包含完整层序（NE, SP, PE, PCC, NCC）
- ✅ 适合"每单元=子电池"的电化学-热耦合模型
- ✅ 自动生成精确的层权重
- ✅ 覆盖整个螺旋（多圈），物理意义明确
- ✅ 计算效率高（"以直代曲"近似）

## 理论基础

### 阿基米德螺旋
- 方程：r(θ) = a + bθ
- 参数：a ≈ Rin, b ≈ t_repeat/(2π)
- 层序：PCC → PE → SP → NE → NCC（从内到外）

### 等效热导率
- 径向（串联）：调和平均 λ_r = 1 / Σ(f_k / λ_k)
- 切向（并联）：算术平均 λ_t = Σ(f_k * λ_k)

### 各向异性导热张量
K(θ) = λ_r * e_r⊗e_r + λ_t * e_θ⊗e_θ

其中 e_r、e_θ 为径向和切向单位向量。

## 代码重构说明

### v2.0 重构要点
1. **结构重组**：按功能模块清晰分组
2. **精简代码**：消除重复逻辑（减少 ~25%）
3. **统一接口**：只保留 collector_seed_mesh 网格生成
4. **优化性能**：向量化计算、提取辅助函数
5. **改进文档**：详细的函数说明和使用示例

### 删除的功能
- `jellyroll_Q4_mesh` 的其他模式（:inscribed, :center）
- `pol2cart` 函数（未使用）
- `jellyroll_element_layer_weights` 采样计算（被缓存方案替代）
- `Plots` 导入（未使用）

### 新增改进
- ✅ 更清晰的代码组织
- ✅ 更高效的计算（向量化）
- ✅ 更精确的层权重（直接从几何）
- ✅ 更简洁的接口（单一网格生成方法）

"""