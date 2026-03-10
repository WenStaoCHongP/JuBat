function jellyroll_collector_seed_mesh(param_dim; nθ::Int=360, gsorder::Int=2, phase::Float64=0.0, tol::Float64=1e-8)
    Rin, Rout = param_dim.cell.Rin, param_dim.cell.Rout
    t_repeat = param_dim.PCC.thickness + 2 * (param_dim.PE.thickness + param_dim.SP.thickness + param_dim.NE.thickness) + param_dim.NCC.thickness
    a = Rin
    b = t_repeat / (2 * pi)

    # 两条螺旋偏移：完整层序 [0, t_repeat]
    s_in = 0.0
    s_out = t_repeat

    # theta 范围裁剪
    theta0 = max(0.0, (Rin - a - s_in) / b)
    theta1 = min((Rout - a - s_out) / b, (Rout - a) / b)
    theta1 > theta0 || error("collector-seeded: no valid theta range [Rin, Rout]")

    # 等角度采样
    seg_per_turn = max(3, nθ)
    dtheta = 2 * pi / seg_per_turn

    # 相位对齐
    k0 = ceil(Int, (theta0 - phase) / dtheta)
    k1 = floor(Int, (theta1 - phase) / dtheta)

    theta = if k1 <= k0
        n_theta_eff = max(2, round(Int, (theta1 - theta0) / dtheta))
        collect(range(theta0, theta1; length=n_theta_eff + 1))
    else
        phase .+ (k0:k1) .* dtheta
    end

    n_theta_actual = length(theta) - 1

    # 生成节点
    r_in = a .+ b .* theta .+ s_in
    r_out = a .+ b .* theta .+ s_out
    x_in = r_in .* cos.(theta)
    y_in = r_in .* sin.(theta)
    x_out = r_out .* cos.(theta)
    y_out = r_out .* sin.(theta)

    # 组装节点数组
    nnode = 2 * (n_theta_actual + 1)
    node = zeros(Float64, nnode, 2)
    for i in 1:(n_theta_actual + 1)
        node[i, 1] = x_in[i]
        node[i, 2] = y_in[i]
        j = (n_theta_actual + 1) + i
        node[j, 1] = x_out[i]
        node[j, 2] = y_out[i]
    end

    # 生成单元
    ne = n_theta_actual
    element = zeros(Int64, ne, 4)
    for i in 1:n_theta_actual
        element[i, 1] = i
        element[i, 2] = (n_theta_actual + 1) + i
        element[i, 3] = (n_theta_actual + 1) + i + 1
        element[i, 4] = i + 1
    end

    # 高斯积分
    gs = GetGS(element, node, gsorder, "Q4")
    mesh_unmerged = Mesh("Q4", 2, node, nnode, element, gs)

    # 识别界面节点对（使用未合并网格）
    inner_nodes = collect(1:(ne + 1))
    outer_nodes = collect((ne + 2):(2 * (ne + 1)))
    interface_pairs = Tuple{Int64, Int64}[]
    for n_out in outer_nodes
        x_out = mesh_unmerged.node[n_out, 1]
        y_out = mesh_unmerged.node[n_out, 2]
        for n_in in inner_nodes
            x_in = mesh_unmerged.node[n_in, 1]
            y_in = mesh_unmerged.node[n_in, 2]
            if abs(x_out - x_in) < tol && abs(y_out - y_in) < tol
                push!(interface_pairs, (n_out, n_in))
                break
            end
        end
    end
    sort!(interface_pairs, by = p -> atan(mesh_unmerged.node[p[1], 2], mesh_unmerged.node[p[1], 1]))

    # 计算热单元分层信息（使用未合并网格）
    element_layer = ones(Int64, ne)
    is_inner_layer = zeros(Bool, ne)
    for e in 1:ne
        nodes = mesh_unmerged.element[e, :]
        x_c = mean(mesh_unmerged.node[nodes, 1])
        y_c = mean(mesh_unmerged.node[nodes, 2])
        r_c = hypot(x_c, y_c)

        layer = max(1, Int(floor((r_c - Rin) / t_repeat) + 1))
        element_layer[e] = layer

        n2, n3 = mesh_unmerged.element[e, 2], mesh_unmerged.element[e, 3]
        r_outer = 0.5 * (
            hypot(mesh_unmerged.node[n2, 1], mesh_unmerged.node[n2, 2]) +
            hypot(mesh_unmerged.node[n3, 1], mesh_unmerged.node[n3, 2])
        )
        is_inner_layer[e] = r_outer < (Rout - t_repeat * 0.1)
    end

    # 计算热单元到CZM单元的映射关系（使用未合并网格）
    czm_element_map = Dict{Int64, Vector{Int64}}()
    for e in 1:ne
        czm_element_map[e] = Int64[]
    end

    n_pairs = length(interface_pairs)
    sorted_pairs = sort(interface_pairs, by = p -> atan(mesh_unmerged.node[p[1], 2], mesh_unmerged.node[p[1], 1]))
    node_to_elem = Dict{Int64, Vector{Int64}}(n => Int64[] for n in 1:mesh_unmerged.nlen)
    for e in 1:ne
        for n in mesh_unmerged.element[e, :]
            push!(node_to_elem[n], e)
        end
    end

    for czm_idx in 1:(n_pairs - 1)
        n_out_1, n_in_1 = sorted_pairs[czm_idx]
        n_out_2, n_in_2 = sorted_pairs[czm_idx + 1]
        related_elems = Int64[]
        for n in (n_out_1, n_out_2, n_in_1, n_in_2)
            append!(related_elems, node_to_elem[n])
        end
        related_elems = unique(related_elems)
        for e in related_elems
            if !(czm_idx in czm_element_map[e])
                push!(czm_element_map[e], czm_idx)
            end
        end
    end

    # 识别极耳节点（使用未合并网格）
    pos_tab_nodes, neg_tab_nodes = jellyroll_tab_node_indices(mesh_unmerged, param_dim)

    # 合并重合节点生成备选热网格
    node_orig = mesh_unmerged.node
    element_orig = mesh_unmerged.element
    nnode_orig = mesh_unmerged.nlen

    merge_map = collect(1:nnode_orig)
    merged_to = zeros(Int, nnode_orig)
    for i in 1:nnode_orig
        merged_to[i] != 0 && continue
        merged_to[i] = i
        for j in (i + 1):nnode_orig
            merged_to[j] != 0 && continue
            dx = node_orig[i, 1] - node_orig[j, 1]
            dy = node_orig[i, 2] - node_orig[j, 2]
            if dx * dx + dy * dy < tol * tol
                merged_to[j] = i
                merge_map[j] = i
            end
        end
    end

    unique_nodes = findall(i -> merged_to[i] == i, 1:nnode_orig)
    nnode_new = length(unique_nodes)
    old_to_new = zeros(Int, nnode_orig)
    for (new_idx, old_idx) in enumerate(unique_nodes)
        old_to_new[old_idx] = new_idx
    end
    for i in 1:nnode_orig
        target = merge_map[i]
        merge_map[i] = old_to_new[target]
    end

    node_new = node_orig[unique_nodes, :]
    element_new = zeros(Int64, ne, 4)
    for e in 1:ne
        for k in 1:4
            old_node = element_orig[e, k]
            element_new[e, k] = merge_map[old_node]
        end
    end

    gs_new = GetGS(element_new, node_new, gsorder, "Q4")
    thermal2D_merged = Mesh("Q4", 2, node_new, nnode_new, element_new, gs_new)

    return (
        thermal2D = mesh_unmerged,
        thermal2D_merged = thermal2D_merged,
        Jellyroll_czm = mesh_unmerged,
        merge_map = merge_map,
        interface_pairs = interface_pairs,
        czm_element_map = czm_element_map,
        element_layer = element_layer,
        is_inner_layer = is_inner_layer,
        inner_nodes = inner_nodes,
        outer_nodes = outer_nodes,
        pos_tab_nodes = pos_tab_nodes,
        neg_tab_nodes = neg_tab_nodes,
        ne = ne,
        nnode = nnode
    )
end

# ========================================================================
# jellyroll_element_properties - 单元属性计算
# ========================================================================

"""
    jellyroll_element_properties(mesh, param_dim) -> (areas, layer_weights)

计算各单元面积和各单元层面积权重。

# 层序（从内到外）
PE → PCC → PE → SP → NE → NCC → NE → SP

# 返回
- `areas`: Vector{Float64}(ne) - 各单元面积
- `layer_weights`: Matrix{Float64}(ne, 5) - 层权重 [NE, SP, PE, PCC, NCC]

# 说明
层权重基于网格几何计算（面积加权）。
对于螺旋结构，不同径向位置的单元各层面积比例略有不同：
- 靠近中心（小半径）：内层（PE/PCC）面积占比略大
- 靠近外侧（大半径）：外层（NE/NCC）面积占比略大
"""
function jellyroll_element_properties(mesh, param_dim)
    ne = size(mesh.element, 1)
    
    # ============ 计算单元面积 ============
    areas = zeros(Float64, ne)
    ngs = length(mesh.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        areas[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    
    # ============ 计算层面积权重 ============
    # 获取各层厚度
    t_PE = param_dim.PE.thickness
    t_NE = param_dim.NE.thickness
    t_SP = param_dim.SP.thickness
    t_PCC = param_dim.PCC.thickness
    t_NCC = param_dim.NCC.thickness
    
    # 层序（从内到外）及其厚度
    # PE → PCC → PE → SP → NE → NCC → NE → SP
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
    layer_weights = zeros(Float64, ne, 5)
    
    @inbounds for e in 1:ne
        # 获取单元4个节点的坐标
        # 节点顺序: 1-内侧起点, 2-外侧起点, 3-外侧终点, 4-内侧终点
        n1, n2, n3, n4 = mesh.element[e, :]
        x1, y1 = mesh.node[n1, 1], mesh.node[n1, 2]
        x2, y2 = mesh.node[n2, 1], mesh.node[n2, 2]
        x3, y3 = mesh.node[n3, 1], mesh.node[n3, 2]
        x4, y4 = mesh.node[n4, 1], mesh.node[n4, 2]
        
        # 计算各节点的极坐标
        r1, theta1 = hypot(x1, y1), atan(y1, x1)
        r2, theta2 = hypot(x2, y2), atan(y2, x2)
        r3, theta3 = hypot(x3, y3), atan(y3, x3)
        r4, theta4 = hypot(x4, y4), atan(y4, x4)
        
        # 单元内边半径（内螺旋边：节点1-4的平均）
        r_in = 0.5 * (r1 + r4)
        
        # 计算单元跨越的角度dtheta（处理角度周期性）
        dtheta_14 = theta4 - theta1
        while dtheta_14 > pi; dtheta_14 -= 2 * pi; end
        while dtheta_14 < -pi; dtheta_14 += 2 * pi; end
        
        dtheta_23 = theta3 - theta2
        while dtheta_23 > pi; dtheta_23 -= 2 * pi; end
        while dtheta_23 < -pi; dtheta_23 += 2 * pi; end
        
        dtheta = 0.5 * (abs(dtheta_14) + abs(dtheta_23))
        dtheta = max(dtheta, 1e-10)
        
        # 从内边开始，依次计算各层的面积
        r_current = r_in
        A_NE, A_SP, A_PE, A_PCC, A_NCC = 0.0, 0.0, 0.0, 0.0, 0.0
        
        for (mat_type, t_layer) in layer_sequence
            r_inner = r_current
            r_outer = r_current + t_layer
            
            # 扇形面积: A = 0.5 * (r_outer^2 - r_inner^2) * dtheta
            A_layer = 0.5 * (r_outer^2 - r_inner^2) * dtheta
            
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
            layer_weights[e, 1] = A_NE / A_total   # NE
            layer_weights[e, 2] = A_SP / A_total   # SP
            layer_weights[e, 3] = A_PE / A_total   # PE
            layer_weights[e, 4] = A_PCC / A_total  # PCC
            layer_weights[e, 5] = A_NCC / A_total  # NCC
        else
            # 回退到厚度权重
            t_total = 2*t_NE + 2*t_SP + 2*t_PE + t_PCC + t_NCC
            layer_weights[e, 1] = 2*t_NE / t_total
            layer_weights[e, 2] = 2*t_SP / t_total
            layer_weights[e, 3] = 2*t_PE / t_total
            layer_weights[e, 4] = t_PCC / t_total
            layer_weights[e, 5] = t_NCC / t_total
        end
    end
    
    return areas, layer_weights
end

# ========================================================================
# edge_boundary - 边界识别
# ========================================================================

"""
    edge_boundary(mesh, nidx, param_dim; which=:inner/:outer, theta_range=nothing, tol=1e-4)

基于螺旋方程的精确边界节点识别。

# 参数
- `mesh`: 网格对象
- `nidx`: 节点索引
- `param_dim`: 参数对象
- `which`: `:inner`（内螺旋）或 `:outer`（外螺旋）
- `theta_range`: theta 范围，默认使用网格覆盖的完整角度区间
- `tol`: 距离容差（m）

# 返回
- `Bool`: 是否为边界节点
"""
function edge_boundary(mesh, nidx::Int, param_dim; which::Symbol=:inner, theta_range::Union{Tuple{Float64,Float64},Nothing}=nothing, tol::Float64=1e-4)
    Rin, Rout = param_dim.cell.Rin, param_dim.cell.Rout
    t_repeat = param_dim.PCC.thickness + 2 * (param_dim.PE.thickness + param_dim.SP.thickness + param_dim.NE.thickness) + param_dim.NCC.thickness
    a = Rin
    b = t_repeat / (2 * pi)
    offset = which === :inner ? 0.0 : (which === :outer ? t_repeat : error("which must be :inner or :outer"))
    
    x, y = mesh.node[nidx, 1], mesh.node[nidx, 2]
    r = hypot(x, y)
    bval = max(b, 1e-12)
    
    # 确定 theta 范围
    if theta_range === nothing
        s_in, s_out = 0.0, t_repeat
        theta_start = 0.0
        theta_end = min((Rout - a - s_out) / bval, (Rout - a) / bval)
        if which === :inner
            theta_min, theta_max = theta_start, min(theta_start + 2.0 * pi, theta_end)
        else
            theta_min, theta_max = max(theta_end - 2.0 * pi, theta_start), theta_end
        end
    else
        theta_min, theta_max = theta_range
    end
    
    # 计算累计角度
    theta_cum = (r - a - offset) / bval
    (theta_cum < theta_min || theta_cum > theta_max) && return false
    
    # 验证距离
    r_theo = a + b * theta_cum + offset
    x_theo = r_theo * cos(theta_cum)
    y_theo = r_theo * sin(theta_cum)
    dist = hypot(x - x_theo, y - y_theo)
    
    return dist <= tol
end

# ========================================================================
# jellyroll_element_centers - 单元中心计算
# ========================================================================

"""
    jellyroll_element_centers(mesh) -> Matrix{Float64}

计算每个 Q4 单元的几何中心 (ne×2)。
"""
function jellyroll_element_centers(mesh)
    ne = size(mesh.element, 1)
    return [mean(mesh.node[mesh.element[e, :], d]) for e in 1:ne, d in 1:2]
end

# ========================================================================
# jellyroll_tab_node_indices - 极耳边界识别
# ========================================================================

"""
    jellyroll_tab_node_indices(mesh, param_dim) -> (pos_indices, neg_indices)

识别受极耳影响的节点索引。

# 返回
- `(pos_indices::Vector{Int}, neg_indices::Vector{Int})`
"""
function jellyroll_tab_node_indices(mesh, param_dim)
    @assert mesh.dimension == 2 "jellyroll_tab_node_indices 仅适用于 2D 网格"

    Rin, Rout = param_dim.cell.Rin, param_dim.cell.Rout
    t_repeat = param_dim.PCC.thickness + 2 * (param_dim.PE.thickness + param_dim.SP.thickness + param_dim.NE.thickness) + param_dim.NCC.thickness
    a = Rin
    b = t_repeat / (2 * pi)
    tw = param_dim.tab.width

    nn = size(mesh.node, 1)
    theta_cum_in = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a) / b for i in 1:nn]
    theta_cum_out = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a - t_repeat) / b for i in 1:nn]

    pos_idx = Int[]
    theta_min, theta_max = minimum(theta_cum_in), maximum(theta_cum_in)
    for theta0_orig in param_dim.tab.theta_pos
        theta0 = Float64(theta0_orig)
        while theta0 > theta_max; theta0 -= 2.0 * pi; end
        while theta0 < theta_min; theta0 += 2.0 * pi; end

            delta_theta = 0.0
            if tw > 1e-12 && b > 0.0
                u0 = a + b * theta0
                F(u) = (u * sqrt(u^2 + b^2) + b^2 * asinh(u / b)) / (2.0 * b)
                s0 = F(u0)

                hi = max(1e-6, tw / max(1e-12, sqrt(u0^2 + b^2)))
                for _ in 1:100
                    (F(u0 + b * hi) - s0) >= tw && break
                    hi *= 2.0
                end

                lo = 0.0
                tol_bsearch = max(1e-12, tw * 1e-9)
                for _ in 1:80
                    mid = 0.5 * (lo + hi)
                    sval = F(u0 + b * mid) - s0
                    if abs(sval - tw) <= tol_bsearch
                        delta_theta = mid
                        break
                    end
                    sval < tw ? (lo = mid) : (hi = mid)
                end
                delta_theta == 0.0 && (delta_theta = 0.5 * (lo + hi))
            end

            theta_start = theta0
            theta_end = theta0 + delta_theta
        for i in 1:nn
            r = hypot(mesh.node[i,1], mesh.node[i,2])
            theta_cum = theta_cum_in[i]
            if (Rin - 1e-8 <= r <= Rout + 1e-8) && (theta_start <= theta_cum <= theta_end)
                push!(pos_idx, i)
            end
        end
    end

    neg_idx = Int[]
    theta_min, theta_max = minimum(theta_cum_out), maximum(theta_cum_out)
    for theta0_orig in param_dim.tab.theta_neg
        theta0 = Float64(theta0_orig)
        while theta0 > theta_max; theta0 -= 2.0 * pi; end
        while theta0 < theta_min; theta0 += 2.0 * pi; end

            delta_theta = 0.0
            if tw > 1e-12 && b > 0.0
                u0 = a + b * theta0
                F(u) = (u * sqrt(u^2 + b^2) + b^2 * asinh(u / b)) / (2.0 * b)
                s0 = F(u0)

                hi = max(1e-6, tw / max(1e-12, sqrt(u0^2 + b^2)))
                for _ in 1:100
                    (F(u0 + b * hi) - s0) >= tw && break
                    hi *= 2.0
                end

                lo = 0.0
                tol_bsearch = max(1e-12, tw * 1e-9)
                for _ in 1:80
                    mid = 0.5 * (lo + hi)
                    sval = F(u0 + b * mid) - s0
                    if abs(sval - tw) <= tol_bsearch
                        delta_theta = mid
                        break
                    end
                    sval < tw ? (lo = mid) : (hi = mid)
                end
                delta_theta == 0.0 && (delta_theta = 0.5 * (lo + hi))
            end

            theta_start = theta0 - delta_theta
            theta_end = theta0
        for i in 1:nn
            r = hypot(mesh.node[i,1], mesh.node[i,2])
            theta_cum = theta_cum_out[i]
            if (Rin - 1e-8 <= r <= Rout + 1e-8) && (theta_start <= theta_cum <= theta_end)
                push!(neg_idx, i)
            end
        end
    end

    return unique(pos_idx), unique(neg_idx)
end



# ========================================================================
# setup_thermal2D_mesh - 设置热网格并保存界面信息
# ========================================================================

"""
    setup_thermal2D_mesh(case, mesh_data; use_merged=false)

设置热网格并保存层间界面信息到新的case中。

# 参数
- `case`: Case对象
- `mesh_data`: 由 `jellyroll_collector_seed_mesh` 返回的网格数据
- `use_merged`: 是否使用合并节点的网格（默认false，使用未合并节点）

# 说明
此函数返回新的case对象，保持原case不变。

# 示例
```julia
param_dim = JuBat.ChooseCell("Jellyroll")
case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
```
"""
function setup_thermal2D_mesh(case, mesh_data; use_merged::Union{Nothing,Bool}=nothing)
    case_new = deepcopy(case)

    # 选择使用哪个网格（默认: czm_enabled=false 使用合并网格）
    use_merged_flag = use_merged === nothing ? !case_new.opt.czm_enabled : use_merged
    if use_merged_flag
        mesh_th = mesh_data.thermal2D_merged
        interface_pairs = Tuple{Int64, Int64}[]
    else
        mesh_th = mesh_data.thermal2D
        interface_pairs = mesh_data.interface_pairs
    end

    case_new.mesh["thermal2D"] = mesh_th

    case_new.multi_spme_layout["interface_pairs"] = interface_pairs
    case_new.multi_spme_layout["element_layer"] = mesh_data.element_layer
    case_new.multi_spme_layout["is_inner_layer"] = mesh_data.is_inner_layer
    case_new.multi_spme_layout["czm_element_map"] = mesh_data.czm_element_map
    case_new.multi_spme_layout["inner_nodes"] = mesh_data.inner_nodes
    case_new.multi_spme_layout["outer_nodes"] = mesh_data.outer_nodes

    ne = size(mesh_th.element, 1)
    nnode = mesh_th.nlen
    n_pairs = length(interface_pairs)
    @debug "Thermal2D mesh setup" ne=ne nnode=nnode n_interface_pairs=n_pairs use_merged=use_merged_flag

    return case_new
end