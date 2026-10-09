# CZM boundary-node identification and Dirichlet enforcement.

# ========================================================================
# 6. 边界条件
# ========================================================================

function apply_bc_czm(K::SparseMatrixCSC{Float64,Int64}, F::Vector{Float64}; bc_nodes=nothing, bc_dofs=nothing, bc_vals=nothing)
    if bc_dofs isa MechanicalBCDOFs
        bc_nodes === nothing || throw(ArgumentError("provide only one BC specification"))
        bc_vals === nothing && throw(ArgumentError("bc_vals required"))
        return apply_mechanical_bc(K,F,bc_dofs,bc_vals)
    end
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

# Internal mechanical boundary description. Directional enforcement lives in
# MechanicalBCSolve.jl; the public Cartesian boundary helpers retain their API.
struct MechanicalBoundary
    nodes::Dict{Int64,Symbol}
    coordinates::Matrix{Float64}
    provenance::Dict{Int64,Vector{Symbol}}
    start_nodes::Set{Int64}
    end_nodes::Set{Int64}
    geometric_inner::Set{Int64}
    geometric_outer::Set{Int64}
    a::Int64
    b::Int64
    inner_count::Int64
    outer_count::Int64
end

function mechanical_endpoint_topology(submesh::CzmSubmesh, elements::Matrix{Int64})
    n_bulk, n_vertices = size(elements)
    n_vertices == 4 || throw(DimensionMismatch("mechanical endpoint topology requires Q4 elements"))
    length(submesh.material_type) == n_bulk || throw(DimensionMismatch(
        "mechanical material count does not match bulk element count $n_bulk"))
    length(submesh.thermal_elem_map) == n_bulk || throw(DimensionMismatch(
        "mechanical thermal map does not match bulk element count $n_bulk"))
    isempty(submesh.thermal_elem_map) && throw(ArgumentError("mechanical thermal map must not be empty"))
    n_segments = maximum(submesh.thermal_elem_map)
    n_segments > 0 || throw(ArgumentError("mechanical angular segment count must be positive"))
    n_bulk % n_segments == 0 || throw(DimensionMismatch(
        "mechanical bulk element count $n_bulk is not divisible by angular segment count $n_segments"))
    n_layers = n_bulk ÷ n_segments

    start_nodes = Set{Int64}()
    end_nodes = Set{Int64}()
    layer_materials = Vector{Symbol}(undef, n_layers)
    for layer in 1:n_layers
        first_elem = (layer - 1) * n_segments + 1
        last_elem = layer * n_segments
        material = submesh.material_type[first_elem]
        all(==(material), @view(submesh.material_type[first_elem:last_elem])) ||
            throw(ArgumentError("mechanical layer $layer has inconsistent material topology"))
        layer_materials[layer] = material
        union!(start_nodes, (elements[first_elem, 1], elements[first_elem, 2]))
        union!(end_nodes, (elements[last_elem, 4], elements[last_elem, 3]))
    end
    sp_layers = findall(==(:SP), layer_materials)
    length(sp_layers) >= 2 || throw(ArgumentError("mechanical endpoints require two SP layers"))
    pe_layer = findfirst(==(:PE), layer_materials)
    pe_layer === nothing && throw(ArgumentError("mechanical endpoints require a PE layer"))
    a = elements[(sp_layers[2] - 1) * n_segments + 1, 2]
    b = elements[pe_layer * n_segments, 4]
    return start_nodes, end_nodes, a, b
end

function _resolve_mechanical_bc(mesh, elements, submesh, param, czm_opt;
    endpoint_variant::Symbol, opt, bonded::Bool)
    czm_opt.outer_bc in (:fixed_xy, :radial_slide) || throw(ArgumentError(
        "unknown mechanical outer boundary $(czm_opt.outer_bc)"))
    endpoint_variant in (:keep, :omit_a, :omit_b, :omit_ab) || throw(ArgumentError(
        "unknown mechanical endpoint variant $endpoint_variant"))
    czm_opt.outer_bc === :radial_slide && czm_opt.geo_nonlinear && throw(ArgumentError(
        "radial_slide supports reference small-displacement geometry; geo_nonlinear must be false"))

    starts, ends, a, b = mechanical_endpoint_topology(submesh, elements)
    is_inner, is_outer = identify_boundary_nodes(mesh, param, opt)
    geometric_inner = Set{Int64}(findall(is_inner))
    geometric_outer = Set{Int64}(findall(is_outer))
    omit_a = endpoint_variant in (:omit_a, :omit_ab)
    omit_b = endpoint_variant in (:omit_b, :omit_ab)
    nodes = Dict{Int64,Symbol}()
    provenance = Dict{Int64,Vector{Symbol}}()

    # Filtering precedes merging: geometric support survives an endpoint variant.
    function contribute!(node, mode, source)
        if source in (:start, :end, :dual_a, :dual_b)
            (omit_a && node == a || omit_b && node == b) && return nothing
        end
        push!(get!(provenance, node, Symbol[]), source)
        if mode === :fixed_xy || !haskey(nodes, node)
            nodes[node] = mode
        end
        return nothing
    end

    outer_mode = czm_opt.outer_bc === :fixed_xy ? :fixed_xy : :radial
    for node in geometric_outer
        contribute!(node, outer_mode, :geometric_outer)
    end
    if czm_opt.fix_inner
        for node in geometric_inner
            contribute!(node, :fixed_xy, :geometric_inner)
        end
    end
    contribute!(a, outer_mode, :dual_a)
    czm_opt.fix_inner && contribute!(b, :fixed_xy, :dual_b)

    add_start = czm_opt.fix_start === nothing ? (!bonded && czm_opt.fix_inner) : czm_opt.fix_start
    add_end = czm_opt.fix_end === nothing ? !bonded : czm_opt.fix_end
    if add_start
        for node in starts
            contribute!(node, :fixed_xy, :start)
        end
    end
    if add_end
        for node in ends
            # Only the legacy CZM end selector omits b when the inner ring is free.
            czm_opt.fix_end === nothing && !czm_opt.fix_inner && node == b && continue
            contribute!(node, :fixed_xy, :end)
        end
    end

    inner_count = czm_opt.fix_inner ? length(geometric_inner) +
        Int(!omit_b && !(b in geometric_inner)) : 0
    outer_count = length(geometric_outer) + Int(!omit_a && !(a in geometric_outer))
    return MechanicalBoundary(nodes, mesh.node, provenance, starts, ends,
        geometric_inner, geometric_outer, a, b, inner_count, outer_count)
end

"""
    resolve_mechanical_bc(czm_mesh, param, czm_opt; endpoint_variant=nothing, opt=nothing)

Resolve circle and layered-endpoint contributions on the actual CZM bulk
topology. `nothing` selects the existing per-field CZM rule at each call.
`endpoint_variant` defaults to the `czm_opt.endpoint_variant` field; an
explicit keyword argument overrides it (the frozen diagnostic collector
relies on this override).
"""
function resolve_mechanical_bc(czm_mesh::CohesiveMesh, param, czm_opt::CzmOptions;
    endpoint_variant::Union{Nothing,Symbol}=nothing, opt=nothing)
    submesh = czm_mesh.czm_submesh
    submesh === nothing && throw(ArgumentError(
        "resolve_mechanical_bc requires a layered CzmSubmesh"))
    return _resolve_mechanical_bc(czm_mesh, czm_mesh.bulk_element, submesh, param, czm_opt;
        endpoint_variant=something(endpoint_variant, czm_opt.endpoint_variant), opt=opt, bonded=false)
end

"""
    resolve_mechanical_bc(submesh, mesh_bonded, param, czm_opt; endpoint_variant=nothing, opt=nothing)

Resolve the corresponding bonded topology. Legacy `nothing` selectors add no
start/end constraints; explicit selectors use this mesh's own Q4 node identities.
The endpoint variant defaults to the `czm_opt.endpoint_variant` field.
"""
function resolve_mechanical_bc(submesh::CzmSubmesh, mesh_bonded::Mesh, param, czm_opt::CzmOptions;
    endpoint_variant::Union{Nothing,Symbol}=nothing, opt=nothing)
    return _resolve_mechanical_bc(mesh_bonded, mesh_bonded.element, submesh, param, czm_opt;
        endpoint_variant=something(endpoint_variant, czm_opt.endpoint_variant), opt=opt, bonded=true)
end

# Internal directional constraint carrier; old Cartesian DOF arrays keep their
# original enforcement and arithmetic. The IDs are the full Cartesian pins;
# radial relations are represented by the orthogonal free-space projector.
struct MechanicalBCDOFs <: AbstractVector{Int64}
    ids::Vector{Int64}
    boundary::MechanicalBoundary
    projector::SparseMatrixCSC{Float64,Int64}
    gauge::Union{Nothing,Vector{Float64}}
end
Base.size(d::MechanicalBCDOFs) = size(d.ids)
Base.getindex(d::MechanicalBCDOFs, i::Int) = d.ids[i]
Base.IndexStyle(::Type{MechanicalBCDOFs}) = IndexLinear()

function mechanical_boundary_dofs(bc::MechanicalBoundary)
    n = size(bc.coordinates, 1)
    ndof = 2n
    ids = Int64[]
    rows = Int64[]; cols = Int64[]; values = Float64[]
    rigid = zeros(ndof, 3)
    rigid[1:2:end, 1] .= 1.0
    rigid[2:2:end, 2] .= 1.0
    radius = maximum(hypot.(bc.coordinates[:,1], bc.coordinates[:,2]))
    isfinite(radius) && radius > 0 || throw(ArgumentError("invalid mechanical coordinates"))
    rigid[1:2:end, 3] .= -bc.coordinates[:,2] ./ radius
    rigid[2:2:end, 3] .= bc.coordinates[:,1] ./ radius
    gram = zeros(3,3)
    for node in 1:n
        i = 2node-1; j = i+1
        kind = get(bc.nodes, node, :free)
        if kind == :fixed_xy
            append!(ids, (i,j))
            for d in (i,j)
                v = rigid[d,:]
                gram .+= v*v'
            end
        elseif kind == :radial
            x,y = bc.coordinates[node,1], bc.coordinates[node,2]
            r = hypot(x,y)
            r > 0 || throw(ArgumentError("radial support at polar origin"))
            nx,ny = x/r,y/r
            tx,ty = -ny,nx
            append!(rows,(i,i,j,j)); append!(cols,(i,j,i,j))
            append!(values,(tx*tx,tx*ty,ty*tx,ty*ty))
            v = nx*rigid[i,:] + ny*rigid[j,:]
            gram .+= v*v'
        elseif kind == :free
            append!(rows,(i,j)); append!(cols,(i,j)); append!(values,(1.0,1.0))
        else
            throw(ArgumentError("unknown mechanical constraint $kind"))
        end
    end
    P = sparse(rows,cols,values,ndof,ndof)
    eigenvalues = eigvals(Symmetric(gram))
    rank = count(v -> v > 1e-10*max(1.0,maximum(eigenvalues)), eigenvalues)
    rank >= 2 || throw(ArgumentError("mechanical support leaves unrestrained translation modes"))
    gauge = nothing
    if rank == 2
        z = rigid[:,3]
        norm(z-P*z) <= 1e-10*norm(z) || throw(ArgumentError(
            "remaining rigid mode is not pure rotation about the polar origin"))
        gauge = z/norm(z)
    end
    dofs = MechanicalBCDOFs(ids,bc,P,gauge)
    return dofs, zeros(length(ids))
end

function apply_mechanical_bc(K::SparseMatrixCSC{Float64,Int64}, F::Vector{Float64},
                             d::MechanicalBCDOFs, vals)
    size(K,1) == size(K,2) == length(F) == size(d.projector,1) ||
        throw(DimensionMismatch("directional BC matrix/load/mesh dimensions differ"))
    length(vals) == length(d) && all(iszero, vals) || throw(ArgumentError(
        "directional mechanical support requires homogeneous prescribed values"))
    scale = maximum(abs,diag(K))
    isfinite(scale) && scale > 0 || throw(ArgumentError("invalid directional stiffness scale"))
    P = d.projector
    # Orthogonal elimination in ambient coordinates. The complement has no
    # coupling to free DOFs, so this is exact rather than a finite-penalty BC.
    complement = spdiagm(0=>ones(length(F))) - P
    return P*K*P + scale*complement, P*F
end

function apply_czm_dirichlet!(u::AbstractVector{Float64}, d::MechanicalBCDOFs,
                              vals::AbstractVector{Float64})
    length(vals) == length(d) && all(iszero,vals) || throw(ArgumentError(
        "directional support is homogeneous"))
    u .= d.projector*u
    d.gauge === nothing || (u .-= d.gauge*dot(d.gauge,u))
    return u
end
function zero_czm_bc_entries!(v::AbstractVector{Float64}, d::MechanicalBCDOFs)
    v .= d.projector*v
    return v
end

function set_czm_residual!(R, u, d::AbstractVector{Int64}, vals)
    for (dof,val) in zip(d,vals)
        R[dof] = val-u[dof]
    end
    return R
end
function set_czm_residual!(R, u, d::MechanicalBCDOFs, vals)
    free_u = d.projector*u
    R .= d.projector*R + free_u-u
    # The gauge fixes only coordinate choice. Leave the physical torque in the
    # equilibrium residual: an unsupported torque must never be hidden.
    d.gauge === nothing || (R .-= d.gauge*dot(d.gauge,u))
    return R
end

function czm_linear_solve(K, R, d::AbstractVector{Int64}; ws=nothing, equilibrated=false)
    equilibrated && return solve_equilibrated(K,R)
    return ws === nothing ? K\R : solve_czm_linear_system_cached!(ws,K,R)
end
function czm_linear_solve(K, R, d::MechanicalBCDOFs; ws=nothing, equilibrated=false)
    if d.gauge === nothing
        delta = equilibrated ? solve_equilibrated(K,R) :
            (ws === nothing ? K\R : solve_czm_linear_system_cached!(ws,K,R))
    else
        z = d.gauge
        norm(K*z) <= 1e-9*max(1.0,opnorm(K,Inf)) || throw(ArgumentError(
            "rotation gauge is not a null mode of the mechanical tangent"))
        abs(dot(z,R)) <= 1e-8*max(1.0,norm(R)) || throw(ArgumentError(
            "radial support cannot balance external torque; rotation gauge is not a support"))
        # A sparse bordered equation imposes zero mean rotation without a
        # dense rank-one stiffness, or a hidden tangential pin at an endpoint.
        column = sparse(reshape(z,:,1))
        A = [K column; transpose(column) spzeros(1,1)]
        rhs = vcat(R,0.0)
        solution = solve_equilibrated(A,rhs)
        delta = solution[1:length(R)]
    end
    apply_czm_dirichlet!(delta,d,zeros(length(d)))
    return delta
end

check_czm_external_load(d::AbstractVector{Int64}, F) = nothing
function check_czm_external_load(d::MechanicalBCDOFs, F)
    d.gauge === nothing && return nothing
    abs(dot(d.gauge,F)) <= 1e-8*max(1.0,norm(F)) || throw(ArgumentError(
        "radial support cannot balance external torque; rotation gauge is not a support"))
    return nothing
end

function validate_mechanical_boundary_mode(boundary,geo_nl)
    boundary !== nothing && geo_nl && any(==(:radial),values(boundary.nodes)) &&
        throw(ArgumentError("radial_slide requires geo_nonlinear=false"))
    return nothing
end

function mechanical_boundary_signature(bc::MechanicalBoundary)
    d,_ = mechanical_boundary_dofs(bc)
    return repr((size(bc.coordinates),d.projector.colptr,d.projector.rowval,
                 d.projector.nzval,d.gauge))
end

function mechanical_boundary_diagnostics(u, raw_R, bc::MechanicalBoundary)
    d,_ = mechanical_boundary_dofs(bc)
    constrained = u-d.projector*u
    gauge_coordinate = d.gauge === nothing ? 0.0 : dot(d.gauge,u)
    free_R = d.projector*raw_R
    z = zeros(length(u))
    z[1:2:end] .= -bc.coordinates[:,2]
    z[2:2:end] .= bc.coordinates[:,1]
    return (free_residual=norm(free_R),
        constraint_violation=hypot(norm(constrained),gauge_coordinate),
        gauge_present=d.gauge !== nothing, gauge_coordinate=gauge_coordinate,
        gauge_moment=d.gauge === nothing ? 0.0 : dot(z,raw_R),
        constraint_count=2count(==(:fixed_xy),values(bc.nodes))+
                         count(==(:radial),values(bc.nodes)),
        reaction=raw_R-free_R)
end
