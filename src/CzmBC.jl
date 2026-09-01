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

    # 分层力学网格在已选定的内/外螺旋边界基础上，
    # 始终叠加每个材料层的螺旋起点与终点。`fix_inner`
    # 只决定是否固定内圈，不影响分层端点约束。
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

    endpoint_nodes = Set{Int64}()
    layer_materials = Vector{Symbol}(undef, n_layers)
    for layer in 1:n_layers
        first_elem = (layer - 1) * n_segments + 1
        last_elem = layer * n_segments
        layer_materials[layer] = submesh.material_type[first_elem]
        union!(endpoint_nodes,
            (czm_mesh.bulk_element[first_elem, 1],
             czm_mesh.bulk_element[first_elem, 2]),
            (czm_mesh.bulk_element[last_elem, 4],
             czm_mesh.bulk_element[last_elem, 3]))
    end

    # 开口卷绕端的两个特定节点保持自由：第二个 SP 的靠外起点、
    # 第一层 PE 的靠内终点。其余端点与原边界集合取并集。
    second_sp_layer = findall(==(:SP), layer_materials)[2]
    first_pe_layer = findfirst(==(:PE), layer_materials)
    first_pe_layer === nothing && error(
        "identify_bc_nodes_czm requires at least one PE layer")
    second_sp_first_elem = (second_sp_layer - 1) * n_segments + 1
    first_pe_last_elem = first_pe_layer * n_segments
    delete!(endpoint_nodes, czm_mesh.bulk_element[second_sp_first_elem, 2])
    delete!(endpoint_nodes, czm_mesh.bulk_element[first_pe_last_elem, 4])

    for node in endpoint_nodes
        bc_nodes[node] = :fixed_xy
    end

    return bc_nodes, inner_count, outer_count
end
