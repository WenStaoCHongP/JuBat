struct JellyrollMesh
    thermal2D::Mesh
    thermal2D_merged::Mesh
    merge_map::Vector{Int}
    interface_pairs::Vector{Tuple{Int,Int}}
    czm_element_map::Dict{Int,Vector{Int}}
    element_layer::Vector{Int}
    is_inner_layer::Vector{Bool}
    inner_nodes::Vector{Int}
    outer_nodes::Vector{Int}
    pos_tab_nodes::Vector{Int}
    neg_tab_nodes::Vector{Int}
    ne::Int
    nnode::Int
    czm_submesh::Union{Nothing, CzmSubmesh}   # v3 新增
end

function jellyroll_collector_seed_mesh(param; nθ::Int=360, gsorder::Int=2, phase::Float64=0.0,
        tol::Float64=1e-8, czm_enabled::Bool=false, thin_subdiv::Int=1)
    nθ >= 3 || throw(ArgumentError("collector-seeded: nθ must be at least 3, got $nθ"))
    isfinite(phase) || throw(ArgumentError("collector-seeded: phase must be finite, got $phase"))
    isfinite(tol) && tol > 0.0 || throw(ArgumentError(
        "collector-seeded: tol must be finite and positive, got $tol"))

    # 参数已归一化，直接使用无量纲值
    a = param.cell.Rin
    b = param.cell.layer / (2 * pi)
    isfinite(a) || throw(ArgumentError("collector-seeded: Rin must be finite, got $a"))
    isfinite(param.cell.Rout) || throw(ArgumentError(
        "collector-seeded: Rout must be finite, got $(param.cell.Rout)"))
    isfinite(b) && b > 0.0 || throw(ArgumentError(
        "collector-seeded: cell layer thickness must be finite and positive, got $(param.cell.layer)"))

    # 两条螺旋偏移：完整层序 [0, param.cell.layer]
    s_in = 0.0
    s_out = param.cell.layer

    # theta 范围裁剪
    theta0 = max(0.0, (param.cell.Rin - a - s_in) / b)
    theta1 = min((param.cell.Rout - a - s_out) / b, (param.cell.Rout - a) / b)
    theta1 > theta0 || error("collector-seeded: no valid theta range [Rin, Rout]")

    # 等角度采样
    seg_per_turn = nθ
    dtheta = 2 * pi / seg_per_turn

    # 相位对齐
    k0 = ceil(Int, (theta0 - phase) / dtheta)
    k1 = floor(Int, (theta1 - phase) / dtheta)
    k1 > k0 || error(
        "collector-seeded: phase-aligned theta range contains no angular segment " *
        "(theta0=$theta0, theta1=$theta1, phase=$phase, nθ=$nθ)")

    theta = phase .+ (k0:k1) .* dtheta

    n_theta_actual = length(theta) - 1
    n_phi_nodes = n_theta_actual + 1 - seg_per_turn
    n_phi_nodes >= 2 || error(
        "collector-seeded: winding span must exceed one full turn to form a Φ interface, " *
        "got $n_theta_actual segments with $seg_per_turn segments per turn")

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

    # Φ 配对由均匀角网格的一圈索引偏移直接生成：outer(θ) ↔ inner(θ+2π)
    inner_nodes = collect(1:(ne + 1))
    outer_nodes = collect((ne + 2):(2 * (ne + 1)))
    n_theta_nodes = ne + 1
    interface_pairs = Vector{Tuple{Int64, Int64}}(undef, n_phi_nodes)
    for outer_col in 1:n_phi_nodes
        n_out = n_theta_nodes + outer_col
        n_in = outer_col + seg_per_turn
        dx = mesh_unmerged.node[n_out, 1] - mesh_unmerged.node[n_in, 1]
        dy = mesh_unmerged.node[n_out, 2] - mesh_unmerged.node[n_in, 2]
        dist = hypot(dx, dy)
        dist <= tol || error(
            "collector-seeded: Φ topology pair ($n_out, $n_in) is not coincident; " *
            "distance=$dist exceeds tol=$tol")
        interface_pairs[outer_col] = (n_out, n_in)
    end

    # 计算热单元分层信息（使用未合并网格）
    element_layer = ones(Int64, ne)
    is_inner_layer = zeros(Bool, ne)
    for e in 1:ne
        nodes = mesh_unmerged.element[e, :]
        x_c = mean(mesh_unmerged.node[nodes, 1])
        y_c = mean(mesh_unmerged.node[nodes, 2])
        r_c = hypot(x_c, y_c)

        layer = Int(floor((r_c - param.cell.Rin) / param.cell.layer) + 1)
        layer >= 1 || error(
            "collector-seeded: element $e has invalid winding layer $layer at radius $r_c")
        element_layer[e] = layer

        n2, n3 = mesh_unmerged.element[e, 2], mesh_unmerged.element[e, 3]
        r_outer = 0.5 * (hypot(mesh_unmerged.node[n2, 1], mesh_unmerged.node[n2, 2])+hypot(mesh_unmerged.node[n3, 1], mesh_unmerged.node[n3, 2]))
        is_inner_layer[e] = r_outer < (param.cell.Rout - param.cell.layer * 0.1)
    end

    # 计算热单元到CZM单元的映射关系（使用未合并网格）
    czm_element_map = Dict{Int64, Vector{Int64}}()
    for e in 1:ne
        czm_element_map[e] = Int64[]
    end

    n_pairs = length(interface_pairs)
    node_to_elem = Dict{Int64, Vector{Int64}}(n => Int64[] for n in 1:mesh_unmerged.nlen)
    for e in 1:ne
        for n in mesh_unmerged.element[e, :]
            push!(node_to_elem[n], e)
        end
    end

    for czm_idx in 1:(n_pairs - 1)
        n_out_1, n_in_1 = interface_pairs[czm_idx]
        n_out_2, n_in_2 = interface_pairs[czm_idx + 1]
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
    pos_tab_nodes, neg_tab_nodes = jellyroll_tab_node_indices(mesh_unmerged, param)

    # 合并重合节点生成备选热网格
    node_orig = mesh_unmerged.node
    element_orig = mesh_unmerged.element
    nnode_orig = mesh_unmerged.nlen

    merge_map = collect(1:nnode_orig)
    for (n_out, n_in) in interface_pairs
        merge_map[n_out] = n_in
    end

    unique_nodes = findall(i -> merge_map[i] == i, 1:nnode_orig)
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

    # 分层力学网格直接继承热网格的实际周向角节点和 Φ 配对
    czm_submesh = czm_enabled ?
        build_czm_submesh(param, mesh_unmerged, theta, interface_pairs; gsorder=gsorder, thin_subdiv=thin_subdiv) : nothing

    Jellyroll_Mesh = JellyrollMesh(mesh_unmerged, thermal2D_merged, merge_map, interface_pairs,czm_element_map, element_layer, is_inner_layer,inner_nodes, outer_nodes, pos_tab_nodes, neg_tab_nodes,ne, nnode, czm_submesh)
    return Jellyroll_Mesh
end

# ========================================================================
# jellyroll_element_properties - 单元属性计算
# ========================================================================

"""
    jellyroll_element_properties(mesh, param) -> (areas, layer_weights)

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
function jellyroll_element_properties(mesh, param)
    ne = size(mesh.element, 1)
    
    # ============ 计算单元面积 ============
    areas = zeros(Float64, ne)
    ngs = length(mesh.gs.detJ)
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        areas[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    
    # ============ 计算层面积权重 ============
    
    # 层序（从内到外）及其厚度
    # PE → PCC → PE → SP → NE → NCC → NE → SP
    layer_sequence = [
        (:PE,  param.PE.thickness),   # 层1
        (:PCC, param.PCC.thickness),  # 层2
        (:PE,  param.PE.thickness),   # 层3
        (:SP,  param.SP.thickness),   # 层4
        (:NE,  param.NE.thickness),   # 层5
        (:NCC, param.NCC.thickness),  # 层6
        (:NE,  param.NE.thickness),   # 层7
        (:SP,  param.SP.thickness),   # 层8
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
        A_total > 0 || error("jellyroll_element_properties: element $e has zero total area")

        # 归一化得到权重
        layer_weights[e, 1] = A_NE / A_total   # NE
        layer_weights[e, 2] = A_SP / A_total   # SP
        layer_weights[e, 3] = A_PE / A_total   # PE
        layer_weights[e, 4] = A_PCC / A_total  # PCC
        layer_weights[e, 5] = A_NCC / A_total  # NCC
    end
    
    return areas, layer_weights
end

# ========================================================================
# edge_boundary - 边界识别
# ========================================================================

"""
    edge_boundary(mesh, nidx, param; which=:inner/:outer, theta_range, tol=1e-4)

基于螺旋方程的精确边界节点识别。

# 参数
- `mesh`: 网格对象
- `nidx`: 节点索引
- `param`: 无量纲参数对象
- `which`: `:inner`（内螺旋）或 `:outer`（外螺旋）
- `theta_range`: theta 范围 (theta_min, theta_max)
- `tol`: 距离容差（无量纲）

# 返回
- `Bool`: 是否为边界节点
"""
function edge_boundary(mesh, nidx::Int, param; which::Symbol=:inner, theta_range::Tuple{Float64,Float64}, tol::Float64=1e-4)
    a = param.cell.Rin
    b = param.cell.layer / (2 * pi)
    offset = which === :inner ? 0.0 : (which === :outer ? param.cell.layer : error("which must be :inner or :outer"))

    x, y = mesh.node[nidx, 1], mesh.node[nidx, 2]
    r = hypot(x, y)
    b > 0 || error("edge_boundary: param.cell.layer must be positive, got $(param.cell.layer)")

    theta_min, theta_max = theta_range

    # 计算累计角度
    theta_cum = (r - a - offset) / b
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
    jellyroll_tab_node_indices(mesh, param) -> (pos_indices, neg_indices)

识别受极耳影响的节点索引。

# 返回
- `(pos_indices::Vector{Int}, neg_indices::Vector{Int})`
"""
function jellyroll_tab_node_indices(mesh, param)
    @assert mesh.dimension == 2 "jellyroll_tab_node_indices 仅适用于 2D 网格"

    a = param.cell.Rin
    b = param.cell.layer / (2 * pi)
    tw = param.tab.width
    Rin = param.cell.Rin
    Rout = param.cell.Rout

    nn = size(mesh.node, 1)
    t_repeat = param.cell.layer  # 一层完整卷绕的厚度
    theta_cum_in = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a) / b for i in 1:nn]
    theta_cum_out = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a - t_repeat) / b for i in 1:nn]

    pos_idx = Int[]
    theta_min, theta_max = minimum(theta_cum_in), maximum(theta_cum_in)
    for theta0_orig in param.tab.theta_pos
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
    for theta0_orig in param.tab.theta_neg
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
    setup_thermal2D_mesh(case, mesh_data; use_merged=nothing)

设置热网格并保存层间界面信息到新的case中。

# 参数
- `case`: Case对象（包含无量纲参数 param 和选项 opt）
- `mesh_data`: 由 `jellyroll_collector_seed_mesh` 返回的网格数据
- `use_merged`: 是否使用合并节点的网格
  - `nothing`（默认）：根据 CZM 是否启用自动决定
    - CZM 未启用 → 使用合并网格（保证径向导热路径）
    - CZM 启用 → 使用未合并网格 + 界面热阻模型
  - `true`: 强制使用合并网格
  - `false`: 强制使用未合并网格

# 说明
此函数返回新的case对象，保持原case不变。

# 示例
```julia
param_dim = JuBat.ChooseCell("Jellyroll")
param = JuBat.NormaliseParam(param_dim)
case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param; nθ=80, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
```
"""
function setup_thermal2D_mesh(case, mesh_data; use_merged::Union{Bool,Nothing}=nothing)
    case_new = deepcopy(case)

    # ============== [v2 修订 2026-07-21] 界面热阻暂禁用（spec §2.4，先验证 CZM 本构）==============
    # 原代码：
    #     # 根据 CZM 启用状态自动决定是否使用合并网格
    #     if isnothing(use_merged)
    #         use_merged = !getfield(case_new.opt, :czm_enabled)
    #         @debug "Auto-selecting thermal mesh" czm_enabled=case_new.opt.czm_enabled use_merged=use_merged
    #     end
    # 替代：强制使用合并网格（径向连续导热），CZM 启用与否不影响热网格选择
    # =========================================================================================
    if isnothing(use_merged)
        use_merged = true   # v2 修订：强制合并网格
        @debug "Auto-selecting thermal mesh (v2: forced merged)" czm_enabled=case_new.opt.czm_enabled use_merged=use_merged
    end

    if use_merged
        mesh_th = mesh_data.thermal2D_merged
        interface_pairs = Tuple{Int64, Int64}[]
    else
        mesh_th = mesh_data.thermal2D
        interface_pairs = mesh_data.interface_pairs
    end

    case_new.mesh["thermal2D"] = mesh_th

    _, layer_weights = jellyroll_element_properties(case_new.mesh["thermal2D"], case_new.param)

    # 预计算边界边缓存（网格不变量）
    is_inner, is_outer = identify_boundary_nodes(mesh_th, case_new.param)
    boundary_cache = compute_boundary_edge_cache(mesh_th, is_outer)

    case_new.geometry = MeshGeometry(
        mesh_data.element_layer,
        mesh_data.is_inner_layer,
        layer_weights,
        interface_pairs,
        mesh_data.czm_element_map,
        mesh_data.inner_nodes,
        mesh_data.outer_nodes,
        boundary_cache
    )

    ne = size(mesh_th.element, 1)
    nnode = mesh_th.nlen
    n_pairs = length(interface_pairs)
    @debug "Thermal2D mesh setup" ne=ne nnode=nnode n_interface_pairs=n_pairs use_merged=use_merged

    return case_new
end

# ========================================================================
# build_czm_submesh - v3 内部辅助：构造径向 8 层分层 Q4 子网格
# ========================================================================

"""
    build_czm_submesh(param, thermal2D, theta, thermal_phi_pairs; gsorder) -> CzmSubmesh

构造径向 8 层分层 Q4 子网格。周向节点直接继承热网格，
`thermal_elem_map` 与力学 Φ 节点对均由热网格拓扑直接生成。Φ 只表示
跨匝 outer/inner 配对，不参与 cohesive 面计数；cohesive 总数由四个真实面
乘整条螺旋分段总数 `length(theta)-1` 决定。
不导出（仅由 jellyroll_collector_seed_mesh 调用）。

见 spec §4.1.1。
"""
function build_czm_submesh(param, thermal2D, theta::AbstractVector{<:Real},
        thermal_phi_pairs::Vector{Tuple{Int,Int}}; gsorder::Int, thin_subdiv::Int=1)
    thin_subdiv >= 1 || error("build_czm_submesh: thin_subdiv must be >= 1, got $thin_subdiv")
    # 螺旋几何参数（与粗热网格一致，使用归一化值）
    a = param.cell.Rin
    s_total = param.cell.layer
    b = s_total / (2 * pi)

    # 径向 8 层厚度（按层序 PE→PCC→PE→SP→NE→NCC→NE→SP）
    layer_thicknesses = [
        param.PE.thickness, param.PCC.thickness,
        param.PE.thickness, param.SP.thickness,
        param.NE.thickness, param.NCC.thickness,
        param.NE.thickness, param.SP.thickness,
    ]
    material_sequence = [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    if thin_subdiv > 1
        # Batch 2''（D-B2''-1）：厚涂层（PE/NE）径向等厚细分（用户语义：细分承载压缩/起皱的厚层），材料继承；同材内部边界豁免
        ts = Float64[]; ms = Symbol[]
        for (t, m) in zip(layer_thicknesses, material_sequence)
            k = m in (:PE, :NE) ? thin_subdiv : 1   # 细分厚涂层（PE/NE）：本征应变压缩载体、起皱物理所在
            for _ in 1:k
                push!(ts, t / k); push!(ms, m)
            end
        end
        layer_thicknesses = ts
        material_sequence = ms
    end
    n_layers = length(layer_thicknesses)
    @assert sum(layer_thicknesses) ≈ s_total rtol=1e-6

    n_thermal = size(thermal2D.element, 1)
    n_segments = length(theta) - 1
    n_segments == n_thermal || throw(DimensionMismatch(
        "build_czm_submesh: thermal theta segments $n_segments do not match thermal elements $n_thermal"))
    n_segments > 0 || error("build_czm_submesh: thermal mesh has no angular segment")

    # 节点：(n_layers+1) 条螺旋 × (n_segments+1) 点
    n_spirals = n_layers + 1
    nnode = n_spirals * (n_segments + 1)
    node = zeros(Float64, nnode, 2)

    s_offsets = [0.0; cumsum(layer_thicknesses)]
    for layer_idx in 0:n_layers
        s_offset = s_offsets[layer_idx + 1]
        r = a .+ b .* theta .+ s_offset
        x = r .* cos.(theta)
        y = r .* sin.(theta)
        for k in 1:(n_segments + 1)
            node_idx = layer_idx * (n_segments + 1) + k
            node[node_idx, 1] = x[k]
            node[node_idx, 2] = y[k]
        end
    end

    # 单元
    ne = n_layers * n_segments
    element = zeros(Int64, ne, 4)
    material_type = Vector{Symbol}(undef, ne)
    winding_turn = Vector{Int}(undef, ne)
    thermal_elem_map = Vector{Int}(undef, ne)

    elem_idx = 0
    for layer_idx in 1:n_layers
        s_offset = s_offsets[layer_idx]
        for seg in 1:n_segments
            elem_idx += 1
            inner_base = (layer_idx - 1) * (n_segments + 1)
            outer_base = layer_idx * (n_segments + 1)
            element[elem_idx, 1] = inner_base + seg
            element[elem_idx, 2] = outer_base + seg
            element[elem_idx, 3] = outer_base + seg + 1
            element[elem_idx, 4] = inner_base + seg + 1

            material_type[elem_idx] = material_sequence[layer_idx]

            theta_center = 0.5 * (theta[seg] + theta[seg + 1])
            winding_turn[elem_idx] = floor(Int, (theta_center - theta[1]) / (2 * pi)) + 1
            thermal_elem_map[elem_idx] = seg
        end
    end

    n_theta_nodes = n_segments + 1
    phi_pairs = Tuple{Int,Int}[]
    sizehint!(phi_pairs, length(thermal_phi_pairs))
    for (thermal_outer, thermal_inner) in thermal_phi_pairs
        outer_col = thermal_outer - n_theta_nodes
        inner_col = thermal_inner
        1 <= outer_col <= n_theta_nodes || error(
            "build_czm_submesh: thermal outer Φ node $thermal_outer has invalid column $outer_col")
        1 <= inner_col <= n_theta_nodes || error(
            "build_czm_submesh: thermal inner Φ node $thermal_inner has invalid column $inner_col")
        mechanical_outer = n_layers * n_theta_nodes + outer_col
        mechanical_inner = inner_col
        push!(phi_pairs, (mechanical_outer, mechanical_inner))
    end

    gs = GetGS(element, node, gsorder, "Q4")
    mesh = Mesh("Q4", 2, node, nnode, element, gs)
    return CzmSubmesh(mesh, material_type, winding_turn, thermal_elem_map, phi_pairs)
end
