# Jellyrollmodel.jl
# 果冻卷（卷绕电池）物理建模工具
# 
# 结构组织：
# 1. 几何参数计算
# 2. 网格生成
# 3. 单元属性计算（面积、层权重）
# 4. 边界识别
# 5. 极耳边界识别
# 6. 辅助函数

using LinearAlgebra
using Statistics
using Dates

# ========================================================================
# 1. 几何参数计算
# ========================================================================

export jellyroll_spiral_params

"""
    jellyroll_spiral_params(param_dim)

计算果冻卷阿基米德螺旋参数。

# 层序结构（从内到外）
PE → PCC → PE → SP → NE → NCC → NE → SP → (重复)

# 完整周期厚度
t_repeat = t_PCC + 2*(t_PE + t_SP + t_NE) + t_NCC

# 螺旋方程
r(θ) = a + b*θ，其中 b = t_repeat / (2π)

# 返回字段
- `a, b`: 螺旋方程参数
- `t_repeat`: 完整周期厚度 (m)
- `n_wind`: 绕组圈数
- `Rin, Rout`: 内外半径
- `λ_r, λ_t`: 等效径向/切向热导率（直接引用 cell.lambda_r/t）
"""
function jellyroll_spiral_params(param_dim)
    cell = param_dim.cell
    Rin, Rout = cell.Rin, cell.Rout
    
    # 获取各层厚度
    t_PE = param_dim.PE.thickness
    t_NE = param_dim.NE.thickness
    t_SP = param_dim.SP.thickness
    t_PCC = param_dim.PCC.thickness
    t_NCC = param_dim.NCC.thickness
    
    # 完整周期厚度：PCC + 2*(PE + SP + NE) + NCC
    t_repeat = t_PCC + 2*(t_PE + t_SP + t_NE) + t_NCC
    
    # 螺旋参数
    a = Rin
    b = t_repeat / (2π)
    n_wind = Int(floor((Rout - Rin) / t_repeat))
    
    # 等效热导率（直接使用 Jellyroll.jl 中已计算的值）
    λ_r = cell.lambda_r
    λ_t = cell.lambda_t
    
    return (; a, b, t_repeat, n_wind, Rin, Rout, λ_r, λ_t)
end

# ========================================================================
# 2. 网格生成
# ========================================================================

export jellyroll_collector_seed_mesh

# 全局缓存
const __jr_layer_weights = IdDict{Any, Matrix{Float64}}()
const __jr_element_areas = IdDict{Any, Vector{Float64}}()

"""
    jellyroll_collector_seed_mesh(param_dim; nθ=360, gsorder=2, phase=0.0)

基于集流体"导轨"的条带网格生成器。

# 原理
沿内外螺旋线采样，连接对应点形成 Q4 条带单元。
每个单元天然跨越完整层序。

# 参数
- `nθ`: 每圈的分段数（建议≥160）
- `gsorder`: 高斯积分阶数（默认2）
- `phase`: 相位对齐角度（默认0.0）

# 返回
- `Mesh`: Q4/2D 网格对象

# 自动缓存
- 层权重：通过 `jellyroll_get_layer_weights(mesh)` 获取
- 单元面积：通过 `jellyroll_get_element_areas(mesh)` 获取
"""
function jellyroll_collector_seed_mesh(param_dim; nθ::Int=360, gsorder::Int=2, phase::Float64=0.0)
    p = jellyroll_spiral_params(param_dim)
    a, b = p.a, p.b
    Rin, Rout = p.Rin, p.Rout
    
    # 两条螺旋偏移：完整层序 [0, t_repeat]
    s_in  = 0.0
    s_out = p.t_repeat
    
    # θ 范围裁剪
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
        nθ_eff = max(2, round(Int, (θ1 - θ0) / dθ))
        collect(range(θ0, θ1; length=nθ_eff+1))
    else
        phase .+ (k0:k1) .* dθ
    end
    
    nθ_actual = length(θ) - 1
    
    # 生成节点
    spiral_xy(offset) = begin
        r = a .+ b .* θ .+ offset
        x = r .* cos.(θ)
        y = r .* sin.(θ)
        return x, y
    end
    
    x_in, y_in = spiral_xy(s_in)
    x_out, y_out = spiral_xy(s_out)
    
    # 组装节点数组
    nnode = 2*(nθ_actual+1)
    node = zeros(Float64, nnode, 2)
    for i in 1:(nθ_actual+1)
        node[i, 1] = x_in[i];  node[i, 2] = y_in[i]
        j = (nθ_actual+1) + i
        node[j, 1] = x_out[i]; node[j, 2] = y_out[i]
    end
    
    # 生成单元
    ne = nθ_actual
    element = zeros(Int64, ne, 4)
    for i in 1:nθ_actual
        element[i, 1] = i
        element[i, 2] = (nθ_actual+1) + i
        element[i, 3] = (nθ_actual+1) + i + 1
        element[i, 4] = i + 1
    end
    
    # 高斯积分
    gs = GetGS(element, node, gsorder, "Q4")
    mesh = Mesh("Q4", 2, node, nnode, element, gs)
    
    # 计算并缓存单元面积
    areas = _compute_element_areas(mesh)
    __jr_element_areas[mesh] = areas
    
    # 计算并缓存层权重（基于网格几何的面积权重）
    # 层序（从内到外）: PE → PCC → PE → SP → NE → NCC → NE → SP
    # 输出格式: [NE, SP, PE, PCC, NCC]
    f_k = _compute_layer_area_weights(mesh, param_dim)
    __jr_layer_weights[mesh] = f_k
    
    return mesh
end

# ========================================================================
# 3. 单元属性获取
# ========================================================================

export jellyroll_get_layer_weights, jellyroll_get_element_areas, jellyroll_element_properties

"""
    jellyroll_get_layer_weights(mesh) -> Matrix{Float64}

返回层权重矩阵 (ne×5): [NE, SP, PE, PCC, NCC]

# 说明
层权重基于网格几何计算（面积加权）。
对于螺旋结构，不同径向位置的单元各层面积比例略有不同：
- 靠近中心（小半径）：内层（PE/PCC）面积占比略大
- 靠近外侧（大半径）：外层（NE/NCC）面积占比略大

# 层序（从内到外）
PE → PCC → PE → SP → NE → NCC → NE → SP
"""
function jellyroll_get_layer_weights(mesh)
    return get(__jr_layer_weights, mesh, nothing)
end

"""
    jellyroll_get_element_areas(mesh) -> Vector{Float64}

返回各单元面积向量 (ne,)
"""
function jellyroll_get_element_areas(mesh)
    return get(__jr_element_areas, mesh, nothing)
end

"""
    jellyroll_element_properties(mesh, param_dim) -> (areas, layer_weights)

统一接口：获取单元面积和层权重。
如果缓存不存在，会重新计算。

# 层序（从内到外）
PE → PCC → PE → SP → NE → NCC → NE → SP

# 返回
- `areas`: Vector{Float64}(ne) - 各单元面积
- `layer_weights`: Matrix{Float64}(ne, 5) - 层权重 [NE, SP, PE, PCC, NCC]

# 说明
层权重基于网格几何计算（面积加权），不同单元的权重会略有不同。
"""
function jellyroll_element_properties(mesh, param_dim)
    # 获取或计算单元面积
    areas = jellyroll_get_element_areas(mesh)
    if areas === nothing
        areas = _compute_element_areas(mesh)
        __jr_element_areas[mesh] = areas
    end
    
    # 获取或计算层权重（基于网格几何的面积权重）
    fks = jellyroll_get_layer_weights(mesh)
    if fks === nothing
        fks = _compute_layer_area_weights(mesh, param_dim)
        __jr_layer_weights[mesh] = fks
    end
    
    return areas, fks
end

"""计算单元面积（内部函数）"""
function _compute_element_areas(mesh)
    ne = size(mesh.element, 1)
    A = zeros(Float64, ne)
    ngs = length(mesh.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        A[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    return A
end

# ========================================================================
# 4. 边界识别
# ========================================================================

export edge_boundary

"""
    edge_boundary(mesh, nidx, param_dim; which=:inner/:outer, theta_range=nothing, tol=1e-4)

基于螺旋方程的精确边界节点识别。

# 参数
- `mesh`: 网格对象
- `nidx`: 节点索引
- `param_dim`: 参数对象
- `which`: `:inner`（内螺旋）或 `:outer`（外螺旋）
- `theta_range`: θ 范围，默认使用网格覆盖的完整角度区间
- `tol`: 距离容差（m）

# 返回
- `Bool`: 是否为边界节点
"""
function edge_boundary(mesh, nidx::Int, param_dim; which::Symbol=:inner, theta_range::Union{Tuple{Float64,Float64},Nothing}=nothing, tol::Float64=1e-4)
    p = jellyroll_spiral_params(param_dim)
    offset = which === :inner ? 0.0 : (which === :outer ? p.t_repeat : error("which must be :inner or :outer"))
    
    x, y = mesh.node[nidx, 1], mesh.node[nidx, 2]
    r = hypot(x, y)
    bval = max(p.b, 1e-12)
    
    # 确定 θ 范围
    if theta_range === nothing
        s_in, s_out = 0.0, p.t_repeat
        θ_start = 0.0
        θ_end = min((p.Rout - p.a - s_out) / bval, (p.Rout - p.a) / bval)
        if which === :inner
            θ_min, θ_max = θ_start, min(θ_start + 2.0*π, θ_end)
        else
            θ_min, θ_max = max(θ_end - 2.0*π, θ_start), θ_end
        end
    else
        θ_min, θ_max = theta_range
    end
    
    # 计算累计角度
    θ_cum = (r - p.a - offset) / bval
    (θ_cum < θ_min || θ_cum > θ_max) && return false
    
    # 验证距离
    r_theo = p.a + p.b * θ_cum + offset
    x_theo = r_theo * cos(θ_cum)
    y_theo = r_theo * sin(θ_cum)
    dist = hypot(x - x_theo, y - y_theo)
    
    return dist <= tol
end

# ========================================================================
# 5. 极耳边界识别
# ========================================================================

export jellyroll_tab_node_indices

"""
    jellyroll_tab_node_indices(mesh, param_dim) -> (pos_indices, neg_indices)

识别受极耳影响的节点索引。

# 返回
- `(pos_indices::Vector{Int}, neg_indices::Vector{Int})`
"""
function jellyroll_tab_node_indices(mesh, param_dim)
    @assert mesh.dimension == 2 "jellyroll_tab_node_indices 仅适用于 2D 网格"
    
    p = jellyroll_spiral_params(param_dim)
    a, b = p.a, p.b
    Rin, Rout = p.Rin, p.Rout
    tw = param_dim.tab.width
    
    nn = size(mesh.node, 1)
    θ_cum_in = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a) / b for i in 1:nn]
    θ_cum_out = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a - p.t_repeat) / b for i in 1:nn]
    
    delta_theta_fn = (θ, w) -> _delta_theta_from_width(a, b, θ, w)
    
    pos_idx = _find_tab_nodes(mesh, param_dim.tab.theta_pos, θ_cum_in, (minimum(θ_cum_in), maximum(θ_cum_in)), delta_theta_fn, tw, Rin, Rout)
    neg_idx = _find_tab_nodes(mesh, param_dim.tab.theta_neg, θ_cum_out, (minimum(θ_cum_out), maximum(θ_cum_out)), delta_theta_fn, tw, Rin, Rout; reverse_range=true)
    
    return pos_idx, neg_idx
end

# ========================================================================
# 6. 辅助函数
# ========================================================================

export cart2pol, jellyroll_element_centers

"""
    cart2pol(x, y) -> (r, θ)

笛卡尔坐标转极坐标。
"""
function cart2pol(x::Real, y::Real)
    return hypot(x, y), atan(y, x)
end

"""
    jellyroll_element_centers(mesh) -> Matrix{Float64}

计算每个 Q4 单元的几何中心 (ne×2)。
"""
function jellyroll_element_centers(mesh)
    ne = size(mesh.element, 1)
    return [mean(mesh.node[mesh.element[e, :], d]) for e in 1:ne, d in 1:2]
end

# ------------------------------------------------------------------------
# 内部辅助函数
# ------------------------------------------------------------------------

"""
    _compute_layer_area_weights(mesh, param_dim) -> Matrix{Float64}

基于网格几何计算每个单元的层面积权重。

# 层序（从内螺旋到外螺旋）
PE → PCC → PE → SP → NE → NCC → NE → SP

# 返回
Matrix{Float64}(ne, 5): 各单元的层权重 [NE, SP, PE, PCC, NCC]

# 原理
对于每个Q4单元：
1. 从节点坐标计算内边半径 r_in 和外边半径 r_out
2. 计算单元跨越的角度 Δθ
3. 根据层序，计算每层的面积（使用扇形面积公式）
4. 归一化得到权重

面积计算：A_layer = 0.5 * (r_outer² - r_inner²) * Δθ
"""
function _compute_layer_area_weights(mesh, param_dim)
    ne = size(mesh.element, 1)
    
    # 获取各层厚度
    t_PE = param_dim.PE.thickness
    t_NE = param_dim.NE.thickness
    t_SP = param_dim.SP.thickness
    t_PCC = param_dim.PCC.thickness
    t_NCC = param_dim.NCC.thickness
    
    # 层序（从内到外）及其厚度
    # PE → PCC → PE → SP → NE → NCC → NE → SP
    # 对应材料类型: :PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP
    layer_sequence = [
        (:PE,  t_PE),   # 层1
        (:PCC, t_PCC),  # 层2
        (:PE,  t_PE),   # 层3
        (:SP,  t_SP),   # 层4
        (:NE,  t_NE),   # 层5
        (:NCC, t_NCC),  # 层6
        (:NE,  t_NE),   # 层7
        (:SP,  t_SP),   # 层8
    ]
    
    # 输出权重矩阵：[NE, SP, PE, PCC, NCC]
    f_k = zeros(Float64, ne, 5)
    
    @inbounds for e in 1:ne
        # 获取单元4个节点的坐标
        # 节点顺序: 1-内侧起点, 2-外侧起点, 3-外侧终点, 4-内侧终点
        n1, n2, n3, n4 = mesh.element[e, :]
        x1, y1 = mesh.node[n1, 1], mesh.node[n1, 2]
        x2, y2 = mesh.node[n2, 1], mesh.node[n2, 2]
        x3, y3 = mesh.node[n3, 1], mesh.node[n3, 2]
        x4, y4 = mesh.node[n4, 1], mesh.node[n4, 2]
        
        # 计算各节点的极坐标
        r1, θ1 = hypot(x1, y1), atan(y1, x1)
        r2, θ2 = hypot(x2, y2), atan(y2, x2)
        r3, θ3 = hypot(x3, y3), atan(y3, x3)
        r4, θ4 = hypot(x4, y4), atan(y4, x4)
        
        # 单元内边半径（内螺旋边：节点1-4的平均）
        r_in = 0.5 * (r1 + r4)
        
        # 计算单元跨越的角度Δθ
        # 处理角度跨越 ±π 的情况
        dθ_14 = _angle_diff(θ1, θ4)  # 内边角度差
        dθ_23 = _angle_diff(θ2, θ3)  # 外边角度差
        Δθ = 0.5 * (abs(dθ_14) + abs(dθ_23))
        
        # 确保Δθ为正且合理
        Δθ = max(Δθ, 1e-10)
        
        # 从内边开始，依次计算各层的面积
        r_current = r_in
        A_NE, A_SP, A_PE, A_PCC, A_NCC = 0.0, 0.0, 0.0, 0.0, 0.0
        
        for (mat_type, t_layer) in layer_sequence
            r_inner = r_current
            r_outer = r_current + t_layer
            
            # 扇形面积: A = 0.5 * (r_outer² - r_inner²) * Δθ
            A_layer = 0.5 * (r_outer^2 - r_inner^2) * Δθ
            
            # 累加到对应材料
            if mat_type == :NE
                A_NE += A_layer
            elseif mat_type == :SP
                A_SP += A_layer
            elseif mat_type == :PE
                A_PE += A_layer
            elseif mat_type == :PCC
                A_PCC += A_layer
            elseif mat_type == :NCC
                A_NCC += A_layer
            end
            
            r_current = r_outer
        end
        
        # 总面积
        A_total = A_NE + A_SP + A_PE + A_PCC + A_NCC
        
        # 归一化得到权重
        if A_total > 0
            f_k[e, 1] = A_NE / A_total   # NE
            f_k[e, 2] = A_SP / A_total   # SP
            f_k[e, 3] = A_PE / A_total   # PE
            f_k[e, 4] = A_PCC / A_total  # PCC
            f_k[e, 5] = A_NCC / A_total  # NCC
        else
            # 回退到厚度权重
            t_total = 2*t_NE + 2*t_SP + 2*t_PE + t_PCC + t_NCC
            f_k[e, 1] = 2*t_NE / t_total
            f_k[e, 2] = 2*t_SP / t_total
            f_k[e, 3] = 2*t_PE / t_total
            f_k[e, 4] = t_PCC / t_total
            f_k[e, 5] = t_NCC / t_total
        end
    end
    
    return f_k
end

"""计算两个角度之间的差值（处理周期性）"""
function _angle_diff(θ1::Float64, θ2::Float64)
    d = θ2 - θ1
    # 将差值归一化到 [-π, π]
    while d > π
        d -= 2π
    end
    while d < -π
        d += 2π
    end
    return d
end

"""极耳节点查找（内部函数）"""
function _find_tab_nodes(mesh, tab_angles, θ_cum_nodes, θ_cum_range, delta_theta_fn, tw, Rin, Rout; reverse_range=false)
    isempty(tab_angles) && return Int[]
    
    idx = Int[]
    θc_min, θc_max = θ_cum_range
    nn = length(θ_cum_nodes)
    
    for θ0_orig in tab_angles
        θ0 = Float64(θ0_orig)
        
        while θ0 > θc_max; θ0 -= 2.0*π; end
        while θ0 < θc_min; θ0 += 2.0*π; end
        
        Δθ = delta_theta_fn(θ0, tw)
        θ_start, θ_end = reverse_range ? (θ0 - Δθ, θ0) : (θ0, θ0 + Δθ)
        
        for i in 1:nn
            r = hypot(mesh.node[i,1], mesh.node[i,2])
            θ_cum = θ_cum_nodes[i]
            if (Rin - 1e-8 <= r <= Rout + 1e-8) && (θ_start <= θ_cum <= θ_end)
                push!(idx, i)
            end
        end
    end
    
    return unique(idx)
end

"""从弧长计算角度增量（内部函数）"""
function _delta_theta_from_width(a::Float64, b::Float64, θ0::Float64, width::Float64)
    (width <= 0.0 || b <= 0.0) && return 0.0
    width < 1e-12 && return 0.0
    
    u0 = a + b * θ0
    F(u) = (u * sqrt(u^2 + b^2) + b^2 * asinh(u / b)) / (2.0 * b)
    s0 = F(u0)
    
    hi = max(1e-6, width / max(1e-12, sqrt(u0^2 + b^2)))
    for _ in 1:100
        (F(u0 + b*hi) - s0) >= width && break
        hi *= 2.0
    end
    
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
# 调试日志工具
# ========================================================================

"""写入调试日志"""
function _debug_log(opt, msg::String)
    opt.debug_coupling || return
    try
        open(opt.debug_log_path, "a") do f
            println(f, "[$(Dates.now())] $msg")
        end
    catch
        # 忽略写入错误
    end
end
