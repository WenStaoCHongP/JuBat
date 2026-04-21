mutable struct CohesiveElement <: AbstractCohesiveElement
    id::Int64
    nodes::Vector{Int64}           # [n1, n2, n3, n4]
    nodes_bottom::Vector{Int64}    # [n1, n2] 底面节点
    nodes_top::Vector{Int64}       # [n4, n3] 顶面节点（注意顺序与底面一致）
    length::Float64                # 单元长度
    layer_idx::Int64               # 层间界面索引
end

"""
    DamageState - 内聚力单元的损伤状态

存储每个内聚力单元（或高斯点）的损伤历史。

# 字段
- `D`: 当前损伤变量 [0, 1]，0=完好，1=完全断裂
- `δ_max_n`: 历史最大法向分离位移
- `δ_max_t`: 历史最大切向分离位移
- `δ_max_eff`: 历史最大等效分离位移
- `fractured`: 是否已完全断裂
- `accumulated_damage`: 累积损伤（用于循环加载）
"""
mutable struct DamageState <: AbstractDamageState
    D::Float64                     # 当前损伤变量
    δ_max_n::Float64              # 历史最大法向分离
    δ_max_t::Float64              # 历史最大切向分离
    δ_max_eff::Float64            # 历史最大等效分离
    fractured::Bool               # 是否断裂
    accumulated_damage::Float64   # 累积损伤（循环）
    
    # 默认构造函数
    DamageState() = new(0.0, 0.0, 0.0, 0.0, false, 0.0)
end

"""
    create_czm_mesh(thermal_mesh, param_dim; tol=1e-8)

基于热网格创建内聚力网格。

# 核心算法
通过检查节点坐标是否重合来识别层间界面：
- 外螺旋在theta位置的点 与 内螺旋在theta+2*pi位置的点 坐标重合
- 这些重合点就是相邻卷绕圈之间的界面
- 在界面处复制节点并插入内聚力单元

# 参数
- `thermal_mesh`: 热分析网格（Q4单元）
- `param_dim`: 参数对象（包含螺旋几何信息）
- `tol`: 坐标重合判断容差（默认1e-8）

# 返回
- `CohesiveMesh`: 内聚力网格对象
"""
function create_czm_mesh(thermal_mesh::Mesh, param_dim; tol::Float64=1e-8)
    @assert thermal_mesh.type == "Q4" "create_czm_mesh requires Q4 mesh"
    @assert thermal_mesh.dimension == 2 "create_czm_mesh requires 2D mesh"
    ne = size(thermal_mesh.element, 1)
    nnode = thermal_mesh.nlen
    n_theta = ne
    inner_nodes = collect(1:(n_theta + 1))
    outer_nodes = collect((n_theta + 2):(2 * (n_theta + 1)))

    # 通过坐标重合检测界面节点对
    interface_pairs = Tuple{Int64, Int64}[]
    for n_out in outer_nodes
        x_out = thermal_mesh.node[n_out, 1]
        y_out = thermal_mesh.node[n_out, 2]

        for n_in in inner_nodes
            x_in = thermal_mesh.node[n_in, 1]
            y_in = thermal_mesh.node[n_in, 2]

            if abs(x_out - x_in) < tol && abs(y_out - y_in) < tol
                push!(interface_pairs, (n_out, n_in))
                break
            end
        end
    end

    sort!(interface_pairs, by = p -> atan(thermal_mesh.node[p[1], 2], thermal_mesh.node[p[1], 1]))

    # 创建内聚力单元（直接使用原网格节点）
    cohesive_elements = CohesiveElement[]
    n_pairs = length(interface_pairs)
    for i in 1:(n_pairs - 1)
        n_out_1, n_in_1 = interface_pairs[i]
        n_out_2, n_in_2 = interface_pairs[i + 1]

        x_out_1, y_out_1 = thermal_mesh.node[n_out_1, 1], thermal_mesh.node[n_out_1, 2]
        x_out_2, y_out_2 = thermal_mesh.node[n_out_2, 1], thermal_mesh.node[n_out_2, 2]
        x_in_1, y_in_1 = thermal_mesh.node[n_in_1, 1], thermal_mesh.node[n_in_1, 2]
        x_in_2, y_in_2 = thermal_mesh.node[n_in_2, 1], thermal_mesh.node[n_in_2, 2]

        L_out = hypot(x_out_2 - x_out_1, y_out_2 - y_out_1)
        L_in = hypot(x_in_2 - x_in_1, y_in_2 - y_in_1)
        elem_length = 0.5 * (L_out + L_in)

        coh_elem = CohesiveElement(
            i,
            [n_in_1, n_in_2, n_out_2, n_out_1],
            [n_in_1, n_in_2],
            [n_out_1, n_out_2],
            elem_length,
            1
        )

        push!(cohesive_elements, coh_elem)
    end

    n_cohesive = length(cohesive_elements)
    damage_states = [DamageState() for _ in 1:n_cohesive]

    node_map = Dict{Int64, Vector{Int64}}()
    for i in 1:nnode
        node_map[i] = [i]
    end
    for (n_out, n_in) in interface_pairs
        push!(node_map[n_out], n_in)
    end

    czm_mesh = CohesiveMesh()
    czm_mesh.bulk_mesh = thermal_mesh
    czm_mesh.node = copy(thermal_mesh.node)
    czm_mesh.nnode = nnode
    czm_mesh.bulk_element = copy(thermal_mesh.element)
    czm_mesh.cohesive_elements = cohesive_elements
    czm_mesh.n_cohesive = n_cohesive
    czm_mesh.n_layers = 2
    czm_mesh.node_map = node_map
    czm_mesh.interface_nodes = [interface_pairs]
    czm_mesh.damage_states = damage_states

    return czm_mesh
end



# ========================================================================
# 2. 系统组装
# ========================================================================


"""
    assemble_czm_system(czm_mesh, u, cohesive_params; damage_states=nothing)

组装内聚力单元的全局刚度矩阵和内力向量。

# 返回
- `K_coh`: 内聚力刚度矩阵 (ndof × ndof)
- `f_int_coh`: 内聚力内力向量 (ndof,)
- `separations`: 每个单元的分离位移
- `tractions`: 每个单元的牵引力
"""
function assemble_czm_system(czm_mesh::CohesiveMesh, u::Vector{Float64}, cohesive_params::Cohesive; damage_states=nothing, geom_cache::Union{Nothing, Vector{CohesiveElementGeom}}=nothing, ws::Union{Nothing, CZMAssemblyWorkspace}=nothing)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive
    states = damage_states === nothing ? czm_mesh.damage_states : damage_states

    # 使用或创建工作区
    if ws === nothing
        ws = CZMAssemblyWorkspace(ndof, n_coh)
    end

    # 每轮重置
    fill!(ws.f_int_coh, 0.0)

    K_coh = ws.K_coh
    # 首次调用时构建稀疏结构（基于 DOF 映射），后续只清零 nonzero 值
    if size(K_coh, 1) != ndof || nnz(K_coh) == 0
        # 构建 sparsity pattern：每个 cohesive 单元贡献 8×8 块
        I_pat = Int64[]
        J_pat = Int64[]
        sizehint!(I_pat, n_coh * 64)
        sizehint!(J_pat, n_coh * 64)
        for i in 1:n_coh
            if geom_cache !== nothing
                dofs = geom_cache[i].dofs
            else
                elem = czm_mesh.cohesive_elements[i]
                n1, n2 = elem.nodes_bottom
                n4, n3 = elem.nodes_top
                dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3, 2*n4-1, 2*n4]
            end
            for a in 1:8
                for b in 1:8
                    push!(I_pat, dofs[a])
                    push!(J_pat, dofs[b])
                end
            end
        end
        V_pat = zeros(Float64, length(I_pat))
        K_coh = sparse(I_pat, J_pat, V_pat, ndof, ndof)
        ws.K_coh = K_coh
    end

    # 清零 nonzero 值（O(nnz)，无分配）
    fill!(nonzeros(K_coh), 0.0)

    @inbounds for i in 1:n_coh
        damage_state = states[i]

        fill!(ws.K_e, 0.0)
        fill!(ws.f_int_e, 0.0)

        T_n_avg = 0.0
        T_t_avg = 0.0
        δ_n_avg = 0.0
        δ_t_avg = 0.0
        w_sum = 0.0

        if geom_cache !== nothing
            geom = geom_cache[i]
            L = geom.length
            R = geom.R
            dofs = geom.dofs
            wts = geom.gauss_wts
            pts = geom.gauss_pts
        else
            elem = czm_mesh.cohesive_elements[i]
            n1, n2 = elem.nodes_bottom
            n4, n3 = elem.nodes_top
            L = elem.length

            if L >= 1e-15
                x1, y1 = czm_mesh.node[n1, 1], czm_mesh.node[n1, 2]
                x2, y2 = czm_mesh.node[n2, 1], czm_mesh.node[n2, 2]
                dx, dy = x2 - x1, y2 - y1
                t_vec = [dx / L, dy / L]
                n_vec = [-t_vec[2], t_vec[1]]
                R = [n_vec[1] n_vec[2]; t_vec[1] t_vec[2]]
            else
                R = [0.0 1.0; 1.0 0.0]
            end
            dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3, 2*n4-1, 2*n4]
            order = czm_mesh.bulk_mesh.gs.order
            wts, pts = NCweight(order)
        end

        if L >= 1e-15
            # 提取单元位移
            ws.u_e[1] = u[dofs[1]]; ws.u_e[2] = u[dofs[2]]
            ws.u_e[3] = u[dofs[3]]; ws.u_e[4] = u[dofs[4]]
            ws.u_e[5] = u[dofs[5]]; ws.u_e[6] = u[dofs[6]]
            ws.u_e[7] = u[dofs[7]]; ws.u_e[8] = u[dofs[8]]

            for (ξ, w) in zip(pts, wts)
                N1 = 0.5 * (1.0 - ξ)
                N2 = 0.5 * (1.0 + ξ)

                fill!(ws.B_global, 0.0)
                ws.B_global[1, 1] = -N1; ws.B_global[2, 2] = -N1
                ws.B_global[1, 3] = -N2; ws.B_global[2, 4] = -N2
                ws.B_global[1, 5] = N2;  ws.B_global[2, 6] = N2
                ws.B_global[1, 7] = N1;  ws.B_global[2, 8] = N1

                # B_local = R * B_global  （mul! 无分配）
                mul!(ws.B_local, R, ws.B_global)

                # δ_local = B_local * u_e  （mul! 无分配）
                mul!(ws.δ_local, ws.B_local, ws.u_e)
                δ_n = ws.δ_local[1]
                δ_t = ws.δ_local[2]

                T_n, T_t, _, _ = bilinear_traction_state(δ_n, δ_t, damage_state, cohesive_params)
                dT_dδ = bilinear_tangent(δ_n, δ_t, damage_state, cohesive_params)

                J = L / 2.0
                wJ = w * J

                # BL_dT = B_local' * dT_dδ  （mul! 无分配）
                mul!(ws.BL_dT, transpose(ws.B_local), dT_dδ)

                # K_e += wJ * BL_dT * B_local  （mul! 无分配）
                mul!(ws.BL_dT_B, ws.BL_dT, ws.B_local)
                for a in 1:8
                    for b in 1:8
                        ws.K_e[a, b] += wJ * ws.BL_dT_B[a, b]
                    end
                end

                # f_int_e += wJ * B_local' * [T_n, T_t]  （mul! 无分配）
                ws.T_vec[1] = T_n
                ws.T_vec[2] = T_t
                mul!(ws.BLtT, transpose(ws.B_local), ws.T_vec)
                for a in 1:8
                    ws.f_int_e[a] += wJ * ws.BLtT[a]
                end

                T_n_avg += w * T_n
                T_t_avg += w * T_t
                δ_n_avg += w * δ_n
                δ_t_avg += w * δ_t
                w_sum += w
            end
        end

        if w_sum > 0.0
            T_n_avg /= w_sum
            T_t_avg /= w_sum
            δ_n_avg /= w_sum
            δ_t_avg /= w_sum
        end

        ws.separations[i] = (δ_n_avg, δ_t_avg)
        ws.tractions[i] = (T_n_avg, T_t_avg)

        # 组装到预分配稀疏矩阵（直接索引，无 sparse() 重建）
        for a in 1:8
            ws.f_int_coh[dofs[a]] += ws.f_int_e[a]
            for b in 1:8
                K_coh[dofs[a], dofs[b]] += ws.K_e[a, b]
            end
        end
    end

    return K_coh, ws.f_int_coh, ws.separations, ws.tractions
end

"""
    assemble_bulk_stiffness(czm_mesh, E_eff, ν_eff)

组装固体单元（Q4）的刚度矩阵。
使用与 mechanical.jl 中相同的方法。

# 返回
- `K_bulk`: 固体刚度矩阵 (ndof × ndof)
"""
function assemble_bulk_stiffness(czm_mesh::CohesiveMesh, E_eff::Float64, ν_eff::Float64)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    
    # 需要重新计算高斯积分点，因为节点可能已经改变
    element = czm_mesh.bulk_element
    node = czm_mesh.node
    ne = size(element, 1)
    gsorder = 2
    
    # 稀疏矩阵组装
    I_idx = Int64[]
    J_idx = Int64[]
    K_vals = Float64[]
    
    # 弹性矩阵（平面应力）
    E = E_eff
    ν = ν_eff
    D_mat = E / (1.0 - ν^2) * [1.0 ν 0.0;
                               ν 1.0 0.0;
                               0.0 0.0 (1.0-ν)/2.0]
    
    for e in 1:ne
        # 单元节点
        elem_nodes = element[e, :]
        x_e = node[elem_nodes, 1]
        y_e = node[elem_nodes, 2]
        
        # 单元刚度矩阵
        K_e = zeros(Float64, 8, 8)
        
        IntQ4(x_e, y_e; order=gsorder) do ξ, η, w, dNdx, dNdy, detJ
            # B矩阵
            B = zeros(Float64, 3, 8)
            for i in 1:4
                B[1, 2*i-1] = dNdx[i]
                B[2, 2*i] = dNdy[i]
                B[3, 2*i-1] = dNdy[i]
                B[3, 2*i] = dNdx[i]
            end

            # 积分
            K_e += w * detJ * (B' * D_mat * B)
        end
        
        # 组装到全局
        dofs = Int64[]
        for n in elem_nodes
            push!(dofs, 2*n - 1)
            push!(dofs, 2*n)
        end
        
        for a in 1:8
            for b in 1:8
                push!(I_idx, dofs[a])
                push!(J_idx, dofs[b])
                push!(K_vals, K_e[a, b])
            end
        end
    end
    
    K_bulk = sparse(I_idx, J_idx, K_vals, ndof, ndof)
    
    return K_bulk
end

function assemble_thermal_chemical_load(czm_mesh::CohesiveMesh, E_eff::Float64, ν_eff::Float64,α_eff::Float64, β_n::Float64, β_p::Float64,dT_elem::Vector{Float64}, Δsoc_n_elem::Vector{Float64}, Δsoc_p_elem::Vector{Float64})
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    element = czm_mesh.bulk_element
    node = czm_mesh.node
    ne = size(element, 1)
    
    E = E_eff
    ν = ν_eff
    
    # 计算每个单元的初始应变
    # ε_0 = α*ΔT + β_n*Δsoc_n + β_p*Δsoc_p
    epsilon_0_elem = zeros(Float64, ne)
    @inbounds for e in 1:ne
        epsilon_0_elem[e] = α_eff * dT_elem[e] + β_n * Δsoc_n_elem[e] + β_p * Δsoc_p_elem[e]
    end
    
    gsorder = 2
    
    # 载荷向量
    F_thermo_chem = zeros(Float64, ndof)
    
    for e in 1:ne
        elem_nodes = element[e, :]
        x_e = node[elem_nodes, 1]
        y_e = node[elem_nodes, 2]
        ε_0 = epsilon_0_elem[e]
        
        # 单元载荷向量
        f_e = zeros(Float64, 8)
        
        IntQ4(x_e, y_e; order=gsorder) do ξ, η, w, dNdx, dNdy, detJ
            # 载荷贡献
            # F = ∫ B^T D ε_0 dΩ
            # D ε_0 = E/(1-ν²) * [ε_0*(1+ν), ε_0*(1+ν), 0]^T
            factor = E / (1.0 - ν^2) * ε_0 * (1.0 + ν) * w * detJ

            for i in 1:4
                # F_ux = ∫ dN/dx * σ_xx dΩ = ∫ dN/dx * D11 * ε_0 dΩ
                f_e[2*i - 1] += dNdx[i] * factor
                # F_uy = ∫ dN/dy * σ_yy dΩ = ∫ dN/dy * D11 * ε_0 dΩ
                f_e[2*i] += dNdy[i] * factor
            end
        end
        
        # 组装到全局
        for i in 1:4
            n = elem_nodes[i]
            F_thermo_chem[2*n - 1] += f_e[2*i - 1]
            F_thermo_chem[2*n] += f_e[2*i]
        end
    end
    
    return F_thermo_chem
end

# ========================================================================
# CZM Assembly Cache Builder
# ========================================================================

"""
    build_czm_cache(czm_mesh, E_eff, ν_eff, param)

构建 CZM 装配缓存，包括 K_bulk、cohesive 几何、边界条件。
当 E/ν 不变时，整个缓存可在多次 Newton 迭代中复用。
"""
function build_czm_cache(czm_mesh::CohesiveMesh, E_eff::Float64, ν_eff::Float64, param)
    cache = CZMAssemblyCache()

    # 1. 缓存 K_bulk（最高 ROI：消除 ~60 次/更新 的冗余 bulk 重组）
    cache.K_bulk = assemble_bulk_stiffness(czm_mesh, E_eff, ν_eff)

    # 2. 缓存 bulk DOF 映射
    element = czm_mesh.bulk_element
    ne = size(element, 1)
    cache.bulk_dofs = Vector{Vector{Int64}}(undef, ne)
    for e in 1:ne
        elem_nodes = element[e, :]
        dofs = Vector{Int64}(undef, 8)
        for (i, n) in enumerate(elem_nodes)
            dofs[2*i - 1] = 2*n - 1
            dofs[2*i]     = 2*n
        end
        cache.bulk_dofs[e] = dofs
    end

    # 3. 缓存 cohesive 单元几何
    n_coh = czm_mesh.n_cohesive
    cache.cohesive_geom = Vector{CohesiveElementGeom}(undef, n_coh)
    for (i, elem) in enumerate(czm_mesh.cohesive_elements)
        n1, n2 = elem.nodes_bottom
        n4, n3 = elem.nodes_top
        L = elem.length

        if L >= 1e-15
            x1, y1 = czm_mesh.node[n1, 1], czm_mesh.node[n1, 2]
            x2, y2 = czm_mesh.node[n2, 1], czm_mesh.node[n2, 2]
            dx, dy = x2 - x1, y2 - y1
            t_vec = [dx / L, dy / L]
            n_vec = [-t_vec[2], t_vec[1]]
            R = [n_vec[1] n_vec[2]; t_vec[1] t_vec[2]]
        else
            t_vec = [1.0, 0.0]
            n_vec = [0.0, 1.0]
            R = [0.0 1.0; 1.0 0.0]
        end

        dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3, 2*n4-1, 2*n4]

        order = czm_mesh.bulk_mesh.gs.order
        wts, pts = NCweight(order)

        cache.cohesive_geom[i] = CohesiveElementGeom(
            L, n_vec, t_vec, R, dofs,
            [n1, n2], [n4, n3],
            wts, pts
        )
    end

    # 4. 缓存边界条件
    bc_nodes, _, _ = identify_bc_nodes_czm(czm_mesh, param)
    bc_dofs = Int64[]
    bc_vals = Float64[]
    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_xy
            push!(bc_dofs, 2 * node - 1); push!(bc_vals, 0.0)
            push!(bc_dofs, 2 * node);     push!(bc_vals, 0.0)
        elseif bc_type == :fixed_x
            push!(bc_dofs, 2 * node - 1); push!(bc_vals, 0.0)
        elseif bc_type == :fixed_y
            push!(bc_dofs, 2 * node);     push!(bc_vals, 0.0)
        end
    end
    cache.bc_dofs = bc_dofs
    cache.bc_vals = bc_vals

    # 5. 记录参数用于失效判断
    cache.E_eff = E_eff
    cache.ν_eff = ν_eff

    # 6. 创建可复用工作区（ndof × n_coh，跨时间步复用）
    ndof = 2 * czm_mesh.nnode
    cache.ws = CZMAssemblyWorkspace(ndof, n_coh)

    cache.valid = true

    return cache
end

"""
    ensure_czm_cache(case, czm_mesh, E_eff, ν_eff)

确保 `case.czm_cache` 可用且未过期。如果缓存不存在或参数不匹配则重建。
"""
function ensure_czm_cache(case::Case, czm_mesh::CohesiveMesh, E_eff::Float64, ν_eff::Float64)
    cache = case.czm_cache
    if cache === nothing || !cache.valid ||
       cache.E_eff != E_eff || cache.ν_eff != ν_eff ||
       length(cache.cohesive_geom) != czm_mesh.n_cohesive
        cache = build_czm_cache(czm_mesh, E_eff, ν_eff, case.param)
        case.czm_cache = cache
    end
    return cache
end

function assemble_coupled_system(czm_mesh::CohesiveMesh, u::Vector{Float64},E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive;F_ext::Union{Vector{Float64}, Nothing}=nothing,F_thermo_chem::Union{Vector{Float64}, Nothing}=nothing,damage_states=nothing,K_bulk_cached::Union{Nothing, SparseMatrixCSC{Float64, Int64}}=nothing,geom_cache::Union{Nothing, Vector{CohesiveElementGeom}}=nothing,ws::Union{Nothing, CZMAssemblyWorkspace}=nothing)
    ndof = 2 * czm_mesh.nnode

    # 固体刚度（使用缓存或重新计算）
    K_bulk = K_bulk_cached !== nothing ? K_bulk_cached : assemble_bulk_stiffness(czm_mesh, E_eff, ν_eff)

    # 内聚力刚度和内力（使用几何缓存和工作区）
    K_coh, f_int_coh, separations, tractions = assemble_czm_system(
        czm_mesh, u, cohesive_params; damage_states=damage_states,
        geom_cache=geom_cache, ws=ws)
    
    # 固体内力（线性弹性：f_int = K * u）
    f_int_bulk = K_bulk * u
    
    # 总刚度矩阵
    K_total = K_bulk + K_coh
    
    # 总内力 = 固体内力 + 内聚力内力
    f_int_total = f_int_bulk + f_int_coh
    
    return K_total, f_int_total, separations, tractions
end

function assemble_coupled_system_full(czm_mesh::CohesiveMesh, u::Vector{Float64},E_eff::Float64, ν_eff::Float64,α_eff::Float64, β_n::Float64, β_p::Float64,cohesive_params::Cohesive,dT_elem::Vector{Float64},Δsoc_n_elem::Vector{Float64},Δsoc_p_elem::Vector{Float64};F_ext::Union{Vector{Float64}, Nothing}=nothing,damage_states=nothing,K_bulk_cached::Union{Nothing, SparseMatrixCSC{Float64, Int64}}=nothing,geom_cache::Union{Nothing, Vector{CohesiveElementGeom}}=nothing,ws::Union{Nothing, CZMAssemblyWorkspace}=nothing)
    ndof = 2 * czm_mesh.nnode

    # 组装基本系统
    K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; damage_states=damage_states, K_bulk_cached=K_bulk_cached, geom_cache=geom_cache, ws=ws)
    
    # 热-化学载荷
    F_thermo_chem = assemble_thermal_chemical_load(czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p,dT_elem, Δsoc_n_elem, Δsoc_p_elem)
    
    # 外部载荷
    F_external = F_ext === nothing ? zeros(Float64, ndof) : F_ext
    
    # 残差 = 外力 + 热化学力 - 内力
    R = F_external + F_thermo_chem - f_int_total
    
    return K_total, R, F_thermo_chem, separations, tractions
end


# ========================================================================
# 6. 边界条件
# ========================================================================

function apply_bc_czm(K::SparseMatrixCSC{Float64,Int64}, F::Vector{Float64}; bc_nodes=nothing, bc_dofs=nothing, bc_vals=nothing)
    K_new = copy(K)
    F_new = copy(F)
    penalty = 1e12

    if bc_nodes !== nothing
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
    elseif bc_dofs !== nothing && bc_vals !== nothing
        for (dof, val) in zip(bc_dofs, bc_vals)
            K_new[dof, dof] += penalty
            F_new[dof] = penalty * val
        end
    end

    return K_new, F_new
end

function identify_bc_nodes_czm(czm_mesh::CohesiveMesh, param; opt=nothing)
    nnode = czm_mesh.nnode
    bc_nodes = Dict{Int64, Symbol}()
    
    is_inner, is_outer = identify_boundary_nodes(czm_mesh, param, opt)
    inner_count = 0
    outer_count = 0
    for i in 1:nnode
        if is_inner[i]
            bc_nodes[i] = :fixed_xy
            inner_count += 1
        end
        if is_outer[i]
            bc_nodes[i] = :fixed_xy
            outer_count += 1
        end
    end
    
    return bc_nodes, inner_count, outer_count
end
