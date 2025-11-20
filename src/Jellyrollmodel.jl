# Jellyrollmodel.jl
# Physical modeling helpers for a jellyroll (wound) cylindrical cell
# - Archimedean spiral description of layers
# - Coordinate transforms (polar/cartesian)
# - Material layout across a winding period

using LinearAlgebra
using Plots
using Statistics

# ---------------------------
# Coordinate transforms
# ---------------------------
export cart2pol, pol2cart

"""
    cart2pol(x, y)
Return (r, θ) with θ in radians in (-π, π].
"""
function cart2pol(x::Real, y::Real)
    r = hypot(x, y)
    θ = atan(y, x)
    return r, θ
end

"""
    pol2cart(r, θ)
Return (x, y) from radius r and angle θ (radians).
"""
function pol2cart(r::Real, θ::Real)
    return r*cos(θ), r*sin(θ)
end

# ---------------------------
# Archimedean spiral and winding
# ---------------------------
export jellyroll_spiral_params, material_at
export jellyroll_Q4_mesh, jellyroll_element_centers, jellyroll_effective_K_at
export jellyroll_element_layer_weights
export jellyroll_collector_seed_mesh, jellyroll_get_layer_weights

"""
    jellyroll_spiral_params(param_dim::Params)
Compute Archimedean spiral parameters for the jellyroll:
- r(θ) = a + b θ, where b ≈ t_repeat / (2π) for a thin winding
- a ≈ Rin
Also returns per-winding layer thickness fractions and effective conductivities.
"""
function jellyroll_spiral_params(param_dim)
    cell = param_dim.cell
    Rin = getfield(cell, :Rin)
    Rout = getfield(cell, :Rout)

    # Layer thicknesses per repeating winding (m) from parameters
    # Order within a repeat: PCC -> PE -> SP -> NE -> NCC (inner to outer)
    widths = (
        PCC = getfield(param_dim.PCC, :thickness),
        PE  = getfield(param_dim.PE,  :thickness),
        SP  = getfield(param_dim.SP,  :thickness),
        NE  = getfield(param_dim.NE,  :thickness),
        NCC = getfield(param_dim.NCC, :thickness),
    )
    t_repeat = widths.PCC + widths.PE + widths.SP + widths.NE + widths.NCC

    a = Rin
    b = t_repeat/(2*pi)
    n_wind = Int(floor((Rout - Rin)/t_repeat))

    # Fractions and names in one repeat (radial order from inner to outer within one period)
    # As requested: PCC -> PE -> SP -> NE -> NCC
    fracs = (
        PCC = widths.PCC/t_repeat,
        PE  = widths.PE /t_repeat,
        SP  = widths.SP /t_repeat,
        NE  = widths.NE /t_repeat,
        NCC = widths.NCC/t_repeat,
    )
    names = (:PCC, :PE, :SP, :NE, :NCC)
    # Order and cumulative boundaries (offsets within one repeat)
    order = [(:PCC, widths.PCC), (:PE, widths.PE), (:SP, widths.SP), (:NE, widths.NE), (:NCC, widths.NCC)]
    boundaries = Float64[0.0]
    begin
        acc = 0.0
        for (_, w) in order
            acc += w
            push!(boundaries, acc)
        end
    end

    λ_an = param_dim.NE.lambda
    λ_ca = param_dim.PE.lambda
    λ_sep = param_dim.SP.lambda
    λ_Al = param_dim.PCC.lambda
    λ_Cu = param_dim.NCC.lambda
    λ_r_eff = 1 / (fracs.NE/λ_an + fracs.SP/λ_sep + fracs.PE/λ_ca + fracs.PCC/λ_Al + fracs.NCC/λ_Cu)
    λ_t_eff = fracs.NE*λ_an + fracs.SP*λ_sep + fracs.PE*λ_ca + fracs.PCC*λ_Al + fracs.NCC*λ_Cu

    return (; a, b, t_repeat, n_wind, fracs, names, widths, order, boundaries, λ_r_eff, λ_t_eff, Rin, Rout)
end

"""
    material_at(r, θ, p; logic=:spiral)
统一材料判定函数，返回 (层名::Symbol, 层内局部偏移::Float64)。
- logic = :spiral：按阿基米德螺旋 r = a + bθ 进行判定（真实几何）。
- logic = :rings ：按同心环近似进行判定（与 θ 无关）。
特殊返回：当 r ≤ Rin 返回 (:inner, 0.0)，当 r ≥ Rout 返回 (:outer, 0.0)。
"""
function material_at(r::Real, θ::Real, p; logic::Symbol=:spiral)
    Rin, Rout = p.Rin, p.Rout
    r <= Rin && return :inner, 0.0
    r >= Rout && return :outer, 0.0

    if logic === :spiral
        # 相对一个周期的径向偏移：δ = r - (a + bθ)，再折返到 [0, t_repeat)
        δ = r - (p.a + p.b*θ)
        s = mod(δ, p.t_repeat)
        acc = 0.0
        for (name, w) in p.order
            if s <= acc + w
                return name, s - acc
            end
            acc += w
        end
        return :NCC, s - (p.t_repeat - p.widths.NCC) # 数值边界回退
    elseif logic === :rings
        # 同心环：径向相对位移从 Rin 计，忽略 θ
        ρ = r - Rin
        x = ρ % p.t_repeat
        acc = 0.0
        for (name, w) in p.order
            if x <= acc + w
                return name, x - acc
            end
            acc += w
        end
        return :NCC, x - (p.t_repeat - p.widths.NCC)
    else
        error("Unknown logic=$(logic). Use :spiral or :rings")
    end
end

# ---------------------------
# Q4 mesh for physical modeling (straight-edge approx.)
# ---------------------------
"""
    jellyroll_Q4_mesh(param_dim; nx=100, ny=100, gsorder=2, crop_to_annulus=true, crop_mode=:inscribed)
基于直角 Q4 单元的俯视网格：在 [-Rout,Rout]×[-Rout,Rout] 规则划分；
若 crop_to_annulus=true，则仅保留环域 [Rin,Rout] 内的单元。
- crop_mode=:inscribed（默认）：仅保留四个节点全部位于环域内的单元（保证“完整小电池”假设）。
- crop_mode=:center：按单元中心点位于环域内裁剪（可能保留跨边界的单元）。
返回 Mesh（见 SetMesh.jl 定义）。
"""
function jellyroll_Q4_mesh(param_dim; nx::Int=100, ny::Int=100, gsorder::Int=2, crop_to_annulus::Bool=true, crop_mode::Symbol=:inscribed)
    # 新增模式：:collector-seeded → 生成条带单元，每个单元跨越完整层序
    if crop_mode === :collector_seeded
        # 复用 nx 作为沿 θ 的离散段数；ny 被忽略
        return jellyroll_collector_seed_mesh(param_dim; nθ=nx, gsorder=gsorder)
    end
    p = jellyroll_spiral_params(param_dim)
    Rout, Rin = p.Rout, p.Rin
    domain = [-Rout, Rout, -Rout, Rout]
    base = SetMesh(domain, [nx, ny], "Q4", gsorder)
    if !crop_to_annulus
        return base
    end
    ne = size(base.element,1)
    keep = Int[]
    if crop_mode === :center
        centers = jellyroll_element_centers(base)
        for i in 1:ne
            x = centers[i,1]; y = centers[i,2]
            r = hypot(x, y)
            if (r >= Rin) && (r <= Rout)
                push!(keep, i)
            end
        end
    elseif crop_mode === :inscribed
        # 四节点全部在环域内才保留，避免“半格子”导致热单元不完整
        for i in 1:ne
            nodes = base.element[i, :]
            xy = base.node[nodes, :]
            rs = hypot.(xy[:,1], xy[:,2])
            if minimum(rs) >= Rin && maximum(rs) <= Rout
                push!(keep, i)
            end
        end
    else
        error("Unknown crop_mode=$(crop_mode). Use :inscribed or :center")
    end
    return PickElement(base, collect(keep))
end

"""
    jellyroll_element_centers(mesh::Mesh)
返回每个 Q4 单元的几何中心 (x,y)（四节点坐标的平均）。
"""
function jellyroll_element_centers(mesh)
    ne = size(mesh.element,1)
    centers = zeros(Float64, ne, 2)
    for e in 1:ne
        nodes = mesh.element[e, :]
        xy = mesh.node[nodes, :]
        centers[e, 1] = mean(xy[:,1])
        centers[e, 2] = mean(xy[:,2])
    end
    return centers
end

"""
    jellyroll_effective_K_at(θ, param_dim)
返回在极角 θ 处的各向异性等效导热张量 K(2×2)，
K = λ_r e_r e_r' + λ_t e_θ e_θ'，其中 e_r=[cosθ,sinθ], e_θ=[-sinθ,cosθ]。
当前使用 jellyroll 的径/切向等效热导率 (λ_r_eff, λ_t_eff)。
"""
function jellyroll_effective_K_at(θ::Real, param_dim)
    p = jellyroll_spiral_params(param_dim)
    λr = p.λ_r_eff
    λt = p.λ_t_eff
    c, s = cos(θ), sin(θ)
    # basis vectors
    er = [c, s]
    et = [-s, c]
    # assemble tensor
    K = λr .* (er*er') + λt .* (et*et')
    return K
end

# ---------------------------
# Plot helpers
# ---------------------------
export tab_positions
export jellyroll_tab_node_indices
export edge_boundary

"""
        tab_positions(param_dim)
Return tab anchors on the outer radius for positive (PCC) and negative (NCC) tabs.
Returns a named tuple:
    (pos = [(x,y,θ,side), ...], neg = [(x,y,θ,side), ...])
Angles read from `param_dim.tab.theta_pos/theta_neg` or auto-spaced if empty.
"""
function tab_positions(param_dim)
    p = jellyroll_spiral_params(param_dim)
    a = p.a; b = p.b; Rout = p.Rout
        tab = param_dim.tab
        # if theta arrays are empty -> interpret as no tabs (empty list)
        if isempty(tab.theta_pos)
            θp = Float64[]
        else
            θp = length(tab.theta_pos) == tab.n_pos ? tab.theta_pos : collect(range(0, 2*pi, length=tab.n_pos+1))[1:end-1]
        end
        if isempty(tab.theta_neg)
            θn = Float64[]
        else
            θn = length(tab.theta_neg) == tab.n_neg ? tab.theta_neg : collect(range(pi/tab.n_neg, 2*pi + pi/tab.n_neg, length=tab.n_neg+1))[1:end-1]
        end
    # positive tabs attach to the inner spiral (s_in = 0)
    s_in = 0.0
    pos = [( (a + b*θ + s_in)*cos(θ), (a + b*θ + s_in)*sin(θ), θ, tab.side_pos) for θ in θp]
    # negative tabs attach to the outer spiral (s_out = t_repeat)
    s_out = p.t_repeat
    neg = [( (a + b*θ + s_out)*cos(θ), (a + b*θ + s_out)*sin(θ), θ, tab.side_neg) for θ in θn]
        return (; pos, neg)
end


"""
    edge_boundary(kind::Symbol, args...; kwargs...)

统一边界工具：当 kind==:inner 或 :outer 时，返回与原有函数等价的结果。

此函数集中实现了内/外边界的共同逻辑，保留原有函数作为小包装以保证向后兼容。
"""
function edge_boundary(kind::Symbol, args...; kwargs...)
    if kind === :r
        theta = args[1]; param_dim = args[2]
        p = jellyroll_spiral_params(param_dim)
        # 第二个参数决定 inner/outer，由 kwargs[:which] 指定为 :inner 或 :outer
        which = get(kwargs, :which, :inner)
        if which === :inner
            return p.a + p.b * theta
        elseif which === :outer
            return p.a + p.b * theta + p.t_repeat
        else
            error("edge_boundary(:r) which must be :inner or :outer")
        end
    elseif kind === :theta_range
        param_dim = args[1]
        which = get(kwargs, :which, :inner)
        if which === :inner
            return 0.0, 2.0*pi
        elseif which === :outer
            p = jellyroll_spiral_params(param_dim)
            N = max(1, Int(p.n_wind))
            return 2.0*pi*(N-1), 2.0*pi*N
        else
            error("edge_boundary(:theta_range) which must be :inner or :outer")
        end
    elseif kind === :is_on
        x = args[1]; y = args[2]; param_dim = args[3]
        which = get(kwargs, :which, :inner)
        p = jellyroll_spiral_params(param_dim)
        
        # 自适应容差：对于小螺距，使用更小的长度容差以避免角度容差过大
        # tol_default = min(1e-4, 0.05 * abs(p.b)) 确保 eps_theta < 0.05 弧度（约3度）
        bval = p.b == 0.0 ? 1e-12 : p.b
        tol_default = min(1e-4, 0.05 * abs(bval))
        tol = get(kwargs, :tol, tol_default)
        
        r, φ = cart2pol(x, y)
        # cumulative theta for inner/outer spirals
        # θ_cum_in = (r - a - s_in)/b with s_in=0
        # θ_cum_out = (r - a - s_out)/b with s_out = t_repeat
        θ_cum_in = (r - p.a) / bval
        θ_cum_out = (r - p.a - p.t_repeat) / bval
        eps_theta = tol / max(abs(bval), 1e-12)
        if which === :inner
            θ_start, θ_end = 0.0, 2.0*pi
            return (θ_cum_in >= θ_start - eps_theta) && (θ_cum_in <= θ_end + eps_theta)
        elseif which === :outer
            N = max(1, Int(p.n_wind))
            θ_start, θ_end = 2.0*pi*(N-1), 2.0*pi*N
            return (θ_cum_out >= θ_start - eps_theta) && (θ_cum_out <= θ_end + eps_theta)
        else
            error("edge_boundary(:is_on) which must be :inner or :outer")
        end
    elseif kind === :node_on
        mesh = args[1]; nidx = args[2]; param_dim = args[3]
        tol = get(kwargs, :tol, 1e-4)
        which = get(kwargs, :which, :inner)
        x = mesh.node[nidx,1]; y = mesh.node[nidx,2]
        return edge_boundary(:is_on, x, y, param_dim; tol=tol, which=which)
    else
        error("Unsupported edge_boundary(kind=$(kind)). Use :r, :theta_range, :is_on, or :node_on")
    end
end

"""
    jellyroll_element_layer_weights(mesh::Mesh, param_dim; nsamples_per_dim=4, logic::Symbol=:spiral)

为 Q4 热网格计算每个单元的层权重 f_k（面积分数），顺序为:
    [NE, SP, PE, PCC, NCC]
采用单元内规则采样（nsamples_per_dim×nsamples_per_dim 个点），按螺旋/同心环判层后统计比例。
逻辑:
- logic=:spiral 使用真实螺旋边界 material_at(r,θ, p; logic=:spiral)
- logic=:rings  使用同心环近似 material_at(r,θ, p; logic=:rings)
返回 Matrix{Float64} (ne×5)，各行和≤1（剔除了 inner/outer 区域），若总计数为0则全零。
"""
function jellyroll_element_layer_weights(mesh, param_dim; nsamples_per_dim::Int=4, logic::Symbol=:spiral)
    @assert mesh.type == "Q4" && mesh.dimension == 2 "layer_weights 仅支持 Q4/2D 网格"
    p = jellyroll_spiral_params(param_dim)
    ne = size(mesh.element, 1)
    f = zeros(Float64, ne, 5)  # [NE, SP, PE, PCC, NCC]
    # 局部规则采样点（均匀）s,t ∈ [-1,1]
    ns = max(2, nsamples_per_dim)
    grid = collect(range(-1.0, 1.0; length=ns))
    # Q4 局部形函数
    N_Q4 = function (xi, eta)
        return 0.25 .* [(1 - xi)*(1 - eta),
                        (1 + xi)*(1 - eta),
                        (1 + xi)*(1 + eta),
                        (1 - xi)*(1 + eta)]
    end
    for e in 1:ne
        nodes = mesh.element[e, :]
        xy = mesh.node[nodes, :]
        cnt_NE = 0; cnt_SP = 0; cnt_PE = 0; cnt_PCC = 0; cnt_NCC = 0; cnt_tot = 0
        @inbounds for eta in grid, xi in grid
            N = N_Q4(xi, eta)
            x = N[1]*xy[1,1] + N[2]*xy[2,1] + N[3]*xy[3,1] + N[4]*xy[4,1]
            y = N[1]*xy[1,2] + N[2]*xy[2,2] + N[3]*xy[3,2] + N[4]*xy[4,2]
            r, th = cart2pol(x, y)
            layer, _ = material_at(r, th, p; logic=logic)
            if layer === :NE
                cnt_NE += 1; cnt_tot += 1
            elseif layer === :SP
                cnt_SP += 1; cnt_tot += 1
            elseif layer === :PE
                cnt_PE += 1; cnt_tot += 1
            elseif layer === :PCC
                cnt_PCC += 1; cnt_tot += 1
            elseif layer === :NCC
                cnt_NCC += 1; cnt_tot += 1
            else
                # 忽略 :inner/:outer
            end
        end
        if cnt_tot > 0
            f[e,1] = cnt_NE / cnt_tot
            f[e,2] = cnt_SP / cnt_tot
            f[e,3] = cnt_PE / cnt_tot
            f[e,4] = cnt_PCC / cnt_tot
            f[e,5] = cnt_NCC / cnt_tot
        else
            f[e,1:5] .= 0.0
        end
    end
    return f
end

# ------------------------------------------------------------------
# Collector-seeded band mesh: each quad spans full layer sequence
# ------------------------------------------------------------------

const __jr_layer_weights = IdDict{Any, Matrix{Float64}}()

"""
    jellyroll_collector_seed_mesh(param_dim; nθ=360, gsorder=2)

基于集流体“导轨”的条带网格：
- 构造两条螺旋折线：r_in(θ) = a + bθ + s_in（取 s_in = 0, 对应一周期的起始边，含 PCC 内侧）
- r_out(θ) = a + bθ + s_out（取 s_out = t_repeat，跨越完整层序到下一周期起点前）
- 在 θ∈[θ0,θ1] 上等距采样，连接对应点形成 Q4 条带单元 [(in_i, out_i, out_{i+1}, in_{i+1})]
- 每个单元天然跨越完整层序，因此元素层权重 f_k 恒等于一次层序的体积分数（[NE,SP,PE,PCC,NCC]）。

注意：此网格为“以直代曲”的条带近似，覆盖从 r_in 到 r_out 的有效区域；θ 范围自动裁剪，保证 r_in≥Rin 且 r_out≤Rout。
返回：Mesh（Q4/2D）。可用 `jellyroll_get_layer_weights(mesh)` 获取与该网格关联的 f_k (ne×5)。
"""
function jellyroll_collector_seed_mesh(param_dim; nθ::Int=360, gsorder::Int=2, phase::Float64=0.0)
    p = jellyroll_spiral_params(param_dim)
    a, b = p.a, p.b
    Rin, Rout = p.Rin, p.Rout
    # 两条偏移：完整层序 [0, t_repeat]
    s_in  = 0.0
    s_out = p.t_repeat
    # θ 范围裁剪，确保 r_in>=Rin 且 r_out<=Rout
    θ0 = max(0.0, (Rin - a - s_in)/b)
    θ1 = min((Rout - a - s_out)/b, (Rout - a)/b)
    θ1 > θ0 || error("collector-seeded: no valid theta window within [Rin,Rout]")
    # 等角度布点：将 nθ 直接解释为“每圈分段数”（与 plot.py 的 seg_per_turn 一致）
    # 这样不会因估算圈数过大而把每圈分段数压得过小，避免出现三扇区畸形
    seg_per_turn = max(3, nθ)
    dθ = 2*pi / seg_per_turn
    # 相位对齐：θ = phase + k*dθ，k 取使 θ ∈ [θ0, θ1] 的整数区间
    k0 = ceil(Int, (θ0 - phase) / dθ)
    k1 = floor(Int, (θ1 - phase) / dθ)
    if k1 <= k0
        # 极端窄窗：按步长 dθ 估算应有的段数，至少 2 段
        nθ_eff = max(2, round(Int, (θ1 - θ0)/dθ))
        θ = collect(range(θ0, θ1; length=nθ_eff+1))
    else
        θ = phase .+ (k0:k1) .* dθ
    end
    nθ = length(θ) - 1  # 实际分段数
    # 组装节点：先 in 曲线，后 out 曲线
    function spiral_xy(offset)
        r = a .+ b .* θ .+ offset
        x = r .* cos.(θ)
        y = r .* sin.(θ)
        return x, y
    end
    x_in, y_in = spiral_xy(s_in)
    x_out, y_out = spiral_xy(s_out)
    nnode = 2*(nθ+1)
    node = zeros(Float64, nnode, 2)
    for i in 1:(nθ+1)
        node[i,1] = x_in[i];  node[i,2] = y_in[i]
        j = (nθ+1) + i
        node[j,1] = x_out[i]; node[j,2] = y_out[i]
    end
    # 元素：[(in_i, out_i, out_{i+1}, in_{i+1})], i=1..nθ
    ne = nθ
    element = zeros(Int64, ne, 4)
    for i in 1:nθ
        n1 = i
        n2 = (nθ+1) + i
        n3 = (nθ+1) + i + 1
        n4 = i + 1
        element[i,1] = n1
        element[i,2] = n2
        element[i,3] = n3
        element[i,4] = n4
    end
    gs = GetGS(element, node, gsorder, "Q4")
    mesh = Mesh("Q4", 2, node, size(node,1), element, gs)

    # 为该网格建立并缓存元素层权重（每个元素等同于全层序的体积分数）
    fr = p.fracs
    f_k = zeros(Float64, ne, 5)  # [NE, SP, PE, PCC, NCC]
    @inbounds for e in 1:ne
        f_k[e,1] = fr.NE
        f_k[e,2] = fr.SP
        f_k[e,3] = fr.PE
        f_k[e,4] = fr.PCC
        f_k[e,5] = fr.NCC
    end
    __jr_layer_weights[mesh] = f_k
    return mesh
end

"""
    jellyroll_get_layer_weights(mesh) -> Matrix{Float64} or nothing

返回通过 `jellyroll_collector_seed_mesh` 生成并绑定到该 mesh 的元素层权重 f_k（ne×5），
若不存在（如普通裁剪网格），返回 nothing。
"""
function jellyroll_get_layer_weights(mesh)
    return get(__jr_layer_weights, mesh, nothing)
end

"""
    get_element_layer_weights(mesh, param_dim; nsamples_per_dim=4, logic=:spiral)

Fetch per-element layer weight fractions for a jellyroll thermal mesh.

Preference order:
1. If the mesh already has cached weights from `jellyroll_collector_seed_mesh`, reuse them.
2. Otherwise, compute weights via `jellyroll_element_layer_weights` (requires Q4/2D).

The returned matrix has size (ne, 5) following the ordering `[NE, SP, PE, PCC, NCC]`.
"""
function get_element_layer_weights(mesh, param_dim; nsamples_per_dim::Int=4, logic::Symbol=:spiral)
    # honour cached weights first to maintain consistency with collector-seeded meshes
    cached = jellyroll_get_layer_weights(mesh)
    if cached !== nothing
        return cached
    end

    if !(mesh isa Mesh)
        error("get_element_layer_weights expects a Mesh-like object")
    end
    if mesh.type == "Q4" && mesh.dimension == 2
        return jellyroll_element_layer_weights(mesh, param_dim; nsamples_per_dim=nsamples_per_dim, logic=logic)
    end
    error("get_element_layer_weights currently supports 2D Q4 meshes only")
end

"""
    jellyroll_tab_node_indices(mesh::Mesh, param_dim)

返回受极耳影响的节点索引（pos_indices, neg_indices）。
基于螺旋参数将 tab.width 映射为沿螺旋的角度增量 Δθ，然后在网格节点上选择落在对应角度/半径段附近的节点。

返回：(pos_indices::Vector{Int}, neg_indices::Vector{Int})
"""
function jellyroll_tab_node_indices(mesh, param_dim)
    # 仅处理二维 Q4 网格
    @assert mesh.dimension == 2 "jellyroll_tab_node_indices 仅适用于 2D 网格"
    p = jellyroll_spiral_params(param_dim)
    a, b = p.a, p.b
    Rin, Rout = p.Rin, p.Rout
    tab = param_dim.tab

    # numeric inverse: 给定 arc-length s (沿螺旋线，线性于 θ 近似)，求对应的 Δθ，使得沿螺旋弧长近似为 tab.width
    # 使用简单的 bisection 在 θ 空间内查找 Δθ，使得 |r(θ+Δθ)-r(θ)| ≈ tab.width (近似)
    function delta_theta_from_width(θ0, width)
        # compute Δθ such that arc length along the Archimedean spiral
        # between θ0 and θ0+Δθ equals `width`.
        # Spiral: r(θ) = a + b θ, dr/dθ = b. Arc length integrand: sqrt(r^2 + (dr/dθ)^2)
        # We use closed-form primitive in u = a + b θ: F(u) = (u*sqrt(u^2+b^2) + b^2*asinh(u/b)) / (2b)
        # arc length S(Δθ) = F(u0 + b Δθ) - F(u0). Solve S(Δθ) = width by bisection.
        if width <= 0.0 || b <= 0.0
            return 0.0
        end
        u0 = a + b * θ0
        F(u) = (u * sqrt(u^2 + b^2) + b^2 * asinh(u / b)) / (2.0 * b)
        s0 = F(u0)
        # quick check: tiny-width approximation
        if width < 1e-12
            return 0.0
        end

        # initial upper bound guess (use small-angle approx as starting point)
        hi = max(1e-6, width / max(1e-12, sqrt(u0^2 + b^2)))
        # expand hi until it brackets the solution or reach a generous limit
        max_expand = 100
        for _ in 1:max_expand
            uhi = u0 + b * hi
            if (F(uhi) - s0) >= width
                break
            end
            hi *= 2.0
        end

        lo = 0.0
        # bisection to solve F(u0+b*Δθ)-F(u0) = width
        tol = max(1e-12, width * 1e-9)
        for _ in 1:80
            mid = 0.5 * (lo + hi)
            umid = u0 + b * mid
            sval = F(umid) - s0
            if abs(sval - width) <= tol
                return mid
            end
            if sval < width
                lo = mid
            else
                hi = mid
            end
        end
        return 0.5 * (lo + hi)
    end

    # 将节点转换为极坐标并按累积螺旋角 θ_cum = (r - a)/b 分组
    nn = size(mesh.node, 1)
    pos_idx = Int[]
    neg_idx = Int[]

    # use s_in = 0 (inner) and s_out = t_repeat (outer) when computing cumulative spiral angle
    tw = hasproperty(tab, :width) ? tab.width : 0.0
    s_in = 0.0
    s_out = p.t_repeat

    # precompute node θ_cum ranges for inner/outer spirals
    θ_cum_in_nodes = zeros(Float64, nn)
    θ_cum_out_nodes = zeros(Float64, nn)
    for i in 1:nn
        x = mesh.node[i,1]; y = mesh.node[i,2]
        r, _ = cart2pol(x, y)
        θ_cum_in_nodes[i] = (r - a - s_in) / b
        θ_cum_out_nodes[i] = (r - a - s_out) / b
    end
    θc_in_min = minimum(θ_cum_in_nodes)
    θc_in_max = maximum(θ_cum_in_nodes)
    θc_out_min = minimum(θ_cum_out_nodes)
    θc_out_max = maximum(θ_cum_out_nodes)

    # positive tabs (inner spiral): match θ_cum_in = (r - a - s_in) / b
    if !isempty(tab.theta_pos)
        for θ0_orig in tab.theta_pos
            θ0 = Float64(θ0_orig)
            # normalize θ0 into the node coverage range by adding/subtracting multiples of 2π
            twoπ = 2.0*pi
            while θ0 > θc_in_max
                θ0 -= twoπ
            end
            while θ0 < θc_in_min
                θ0 += twoπ
            end
            Δθ = delta_theta_from_width(θ0, tw)
            θ_start = θ0
            θ_end = θ0 + Δθ
            for i in 1:nn
                r = hypot(mesh.node[i,1], mesh.node[i,2])
                θ_cum_in = θ_cum_in_nodes[i]
                if (r >= Rin - 1e-8) && (r <= Rout + 1e-8)
                    if θ_cum_in >= θ_start && θ_cum_in <= θ_end
                        push!(pos_idx, i)
                    end
                end
            end
        end
    end

    # negative tabs (outer spiral): match θ_cum_out = (r - a - s_out) / b
    if !isempty(tab.theta_neg)
        for θ0_orig in tab.theta_neg
            θ0 = Float64(θ0_orig)
            twoπ = 2.0*pi
            # normalize into outer node coverage range
            while θ0 > θc_out_max
                θ0 -= twoπ
            end
            while θ0 < θc_out_min
                θ0 += twoπ
            end
            Δθ = delta_theta_from_width(θ0, tw)
            θ_start = θ0 - Δθ
            θ_end = θ0
            for i in 1:nn
                r = hypot(mesh.node[i,1], mesh.node[i,2])
                θ_cum_out = θ_cum_out_nodes[i]
                if (r >= Rin - 1e-8) && (r <= Rout + 1e-8)
                    if θ_cum_out >= θ_start && θ_cum_out <= θ_end
                        push!(neg_idx, i)
                    end
                end
            end
        end
    end

    # 去重并返回
    pos_idx = unique(pos_idx)
    neg_idx = unique(neg_idx)
    return pos_idx, neg_idx
end

