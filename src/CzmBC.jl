# CZM boundary-node identification and Dirichlet enforcement.

# ========================================================================
# 6. 边界条件
# ========================================================================

function apply_bc_czm(K::SparseMatrixCSC{Float64,Int64}, F::Vector{Float64}; bc_nodes=nothing, bc_dofs=nothing, bc_vals=nothing)
    nrow, ncol = size(K)
    nrow == ncol || throw(DimensionMismatch(
        "apply_bc_czm requires a square stiffness matrix, got size $(size(K))"
    ))
    nrow > 0 || throw(ArgumentError("apply_bc_czm requires a nonempty stiffness matrix"))
    iseven(nrow) || throw(DimensionMismatch(
        "apply_bc_czm requires two DOFs per node, got matrix size $(size(K))"
    ))
    length(F) == nrow || throw(DimensionMismatch(
        "load vector length $(length(F)) does not match stiffness size $nrow"
    ))

    node_mode = bc_nodes !== nothing
    dof_mode = bc_dofs !== nothing || bc_vals !== nothing
    node_mode == dof_mode && throw(ArgumentError(
        "provide exactly one BC specification: bc_nodes or the pair bc_dofs/bc_vals"
    ))

    if node_mode
        isempty(bc_nodes) && throw(ArgumentError("bc_nodes must not be empty"))
        nnode = nrow ÷ 2
        for (node, bc_type) in bc_nodes
            node isa Integer || throw(ArgumentError(
                "BC node index must be an integer, got $(typeof(node))"
            ))
            1 <= node <= nnode || throw(ArgumentError(
                "BC node index $node is outside 1:$nnode"
            ))
            bc_type in (:fixed_x, :fixed_y, :fixed_xy) || throw(ArgumentError(
                "unknown CZM boundary type $bc_type for node $node"
            ))
        end
    else
        bc_dofs !== nothing && bc_vals !== nothing || throw(ArgumentError(
            "bc_dofs and bc_vals must be provided together"
        ))
        length(bc_dofs) == length(bc_vals) || throw(DimensionMismatch(
            "bc_dofs length $(length(bc_dofs)) does not match bc_vals length $(length(bc_vals))"
        ))
        isempty(bc_dofs) && throw(ArgumentError("bc_dofs and bc_vals must not be empty"))
        length(unique(bc_dofs)) == length(bc_dofs) || throw(ArgumentError(
            "bc_dofs contains duplicate constraints"
        ))
        for (dof, val) in zip(bc_dofs, bc_vals)
            dof isa Integer || throw(ArgumentError(
                "BC DOF index must be an integer, got $(typeof(dof))"
            ))
            1 <= dof <= nrow || throw(ArgumentError(
                "BC DOF index $dof is outside 1:$nrow"
            ))
            val isa Real && isfinite(val) || throw(ArgumentError(
                "BC value for DOF $dof must be finite, got $val"
            ))
        end
    end

    # 相对罚（重设计 v2 §6）：跟随矩阵对角量级，避免固定罚值掩盖非法系统刚度。
    dmax = maximum(abs, diag(K))
    isfinite(dmax) && dmax > 0.0 || throw(ArgumentError(
        "stiffness matrix must have a finite positive diagonal scale, got $dmax"
    ))
    penalty = 1e6 * dmax
    isfinite(penalty) || throw(ArgumentError(
        "CZM penalty is non-finite for stiffness scale $dmax"
    ))

    K_new = copy(K)
    F_new = copy(F)
    if node_mode
        for (node, bc_type) in bc_nodes
            if bc_type == :fixed_x
                dof = 2 * node - 1
                K_new[dof, dof] += penalty
                F_new[dof] = 0.0
            elseif bc_type == :fixed_y
                dof = 2 * node
                K_new[dof, dof] += penalty
                F_new[dof] = 0.0
            elseif bc_type == :fixed_xy
                dof_x = 2 * node - 1
                dof_y = 2 * node
                K_new[dof_x, dof_x] += penalty
                K_new[dof_y, dof_y] += penalty
                F_new[dof_x] = 0.0
                F_new[dof_y] = 0.0
            end
        end
    else
        for (dof, val) in zip(bc_dofs, bc_vals)
            K_new[dof, dof] += penalty
            F_new[dof] = penalty * val
        end
    end

    return K_new, F_new
end

function identify_bc_nodes_czm(czm_mesh::CohesiveMesh, param; opt=nothing, fix_inner::Bool=true)
    nnode = czm_mesh.nnode
    bc_nodes = Dict{Int64, Symbol}()

    is_inner, is_outer = identify_boundary_nodes(czm_mesh, param, opt)
    inner_count = 0
    outer_count = 0
    for i in 1:nnode
        if fix_inner && is_inner[i]
            bc_nodes[i] = :fixed_xy
            inner_count += 1
        end
        if is_outer[i]
            bc_nodes[i] = :fixed_xy
            outer_count += 1
        end
    end

    # 分层端点拓扑（任务59 P2）：从真实 bulk_element 读取 S/E 与 a/b，
    # 不硬编码层节点号、不用坐标重合代替拓扑身份。
    # Q4 节点序为 [靠内起点, 靠外起点, 靠外终点, 靠内终点]；CZM 节点复制后
    # 必须从 bulk_element 读取，才能约束实际承载层上的节点。
    submesh = czm_mesh.czm_submesh
    submesh === nothing && error(
        "identify_bc_nodes_czm requires a layered CzmSubmesh to identify spiral endpoints")
    n_segments = maximum(submesh.thermal_elem_map)
    n_bulk = size(czm_mesh.bulk_element, 1)
    n_bulk % n_segments == 0 || throw(DimensionMismatch(
        "CZM bulk element count $n_bulk is not divisible by angular segment count $n_segments"))
    n_layers = n_bulk ÷ n_segments

    start_nodes = Set{Int64}()
    end_nodes = Set{Int64}()
    layer_materials = Vector{Symbol}(undef, n_layers)
    for layer in 1:n_layers
        first_elem = (layer - 1) * n_segments + 1
        last_elem = layer * n_segments
        layer_materials[layer] = submesh.material_type[first_elem]
        union!(start_nodes,
            (czm_mesh.bulk_element[first_elem, 1],
             czm_mesh.bulk_element[first_elem, 2]))
        union!(end_nodes,
            (czm_mesh.bulk_element[last_elem, 4],
             czm_mesh.bulk_element[last_elem, 3]))
    end

    # a=第二个 SP 层靠外起点（用户指定外圈固定属性）；b=第一层 PE 靠内终点（内圈属性）
    second_sp_layer = findall(==(:SP), layer_materials)[2]
    first_pe_layer = findfirst(==(:PE), layer_materials)
    first_pe_layer === nothing && error(
        "identify_bc_nodes_czm requires at least one PE layer")
    second_sp_first_elem = (second_sp_layer - 1) * n_segments + 1
    first_pe_last_elem = first_pe_layer * n_segments
    a = czm_mesh.bulk_element[second_sp_first_elem, 2]
    b = czm_mesh.bulk_element[first_pe_last_elem, 4]

    if fix_inner
        # 任务59 P2：固定模式 F = O ∪ I ∪ S ∪ E——a、b 补回，不再排除两处开口端
        for node in union(start_nodes, end_nodes)
            bc_nodes[node] = :fixed_xy
        end
    else
        # 任务59 P2：自由模式 F = O ∪ (E\{b})——不再叠加起点集合 S；
        # a 具外圈属性补全固定；b 具自由内圈属性保持自由；其他终点保留原规则
        bc_nodes[a] = :fixed_xy
        for node in end_nodes
            node == b && continue
            bc_nodes[node] = :fixed_xy
        end
    end

    # 统计（P2-R2）：圈属性计数含新增端点——按补全后的 O/I 去重，
    # a 补入外圈计数、b（仅固定模式实际固定时）补入内圈计数；总数为最终固定集合大小
    if !is_outer[a]
        outer_count += 1
    end
    if fix_inner && !is_inner[b]
        inner_count += 1
    end
    return bc_nodes, inner_count, outer_count, length(bc_nodes)
end
