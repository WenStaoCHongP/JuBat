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
function assemble_czm_system(czm_mesh::CohesiveMesh, u::Vector{Float64}, cohesive_params::Cohesive; damage_states=nothing)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive
    states = damage_states === nothing ? czm_mesh.damage_states : damage_states
    
    # 稀疏矩阵组装
    I_idx = Int64[]
    J_idx = Int64[]
    K_vals = Float64[]
    
    f_int_coh = zeros(Float64, ndof)
    
    # 存储每个单元的分离和牵引
    separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
    tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)
    
    for (i, elem) in enumerate(czm_mesh.cohesive_elements)
        damage_state = states[i]
        
        # 单元有4个节点，每节点2个DOF，共8个DOF
        K_e = zeros(Float64, 8, 8)
        f_int_e = zeros(Float64, 8)
        
        n1, n2 = elem.nodes_bottom
        n4, n3 = elem.nodes_top
        x1, y1 = czm_mesh.node[n1, 1], czm_mesh.node[n1, 2]
        x2, y2 = czm_mesh.node[n2, 1], czm_mesh.node[n2, 2]
        L = elem.length
        
        T_n_avg = 0.0
        T_t_avg = 0.0
        δ_n_avg = 0.0
        δ_t_avg = 0.0
        w_sum = 0.0
        
        if L >= 1e-15
            dx = x2 - x1
            dy = y2 - y1
            t_vec = [dx / L, dy / L]
            n_vec = [-t_vec[2], t_vec[1]]
            R = [n_vec[1] n_vec[2]; t_vec[1] t_vec[2]]
            
            order = czm_mesh.bulk_mesh.gs.order
            wts, pts = NCweight(order)

            for (ξ, w) in zip(pts, wts)
                N1 = 0.5 * (1.0 - ξ)
                N2 = 0.5 * (1.0 + ξ)
                
                B_global = zeros(Float64, 2, 8)
                B_global[1, 1] = -N1; B_global[2, 2] = -N1
                B_global[1, 3] = -N2; B_global[2, 4] = -N2
                B_global[1, 5] = N2; B_global[2, 6] = N2
                B_global[1, 7] = N1; B_global[2, 8] = N1
                
                B_local = R * B_global
                
                u_e = zeros(Float64, 8)
                u_e[1] = u[2 * n1 - 1]; u_e[2] = u[2 * n1]
                u_e[3] = u[2 * n2 - 1]; u_e[4] = u[2 * n2]
                u_e[5] = u[2 * n3 - 1]; u_e[6] = u[2 * n3]
                u_e[7] = u[2 * n4 - 1]; u_e[8] = u[2 * n4]
                
                δ_local = B_local * u_e
                δ_n = δ_local[1]
                δ_t = δ_local[2]
                
                T_n, T_t, _, _ = bilinear_traction_state(δ_n, δ_t, damage_state, cohesive_params)
                dT_dδ = bilinear_tangent(δ_n, δ_t, damage_state, cohesive_params)
                
                J = L / 2.0
                K_e += w * J * (B_local' * dT_dδ * B_local)
                T_local = [T_n, T_t]
                f_int_e += w * J * (B_local' * T_local)
                
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
        
        δ_n = δ_n_avg
        δ_t = δ_t_avg
        T_n = T_n_avg
        T_t = T_t_avg
        
        separations[i] = (δ_n, δ_t)
        tractions[i] = (T_n, T_t)
        
        # 获取单元DOF编号
        n1, n2 = elem.nodes_bottom
        n4, n3 = elem.nodes_top
        dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3, 2*n4-1, 2*n4]
        
        # 组装到全局
        for a in 1:8
            f_int_coh[dofs[a]] += f_int_e[a]
            for b in 1:8
                push!(I_idx, dofs[a])
                push!(J_idx, dofs[b])
                push!(K_vals, K_e[a, b])
            end
        end
    end
    
    K_coh = sparse(I_idx, J_idx, K_vals, ndof, ndof)
    
    return K_coh, f_int_coh, separations, tractions
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

function assemble_coupled_system(czm_mesh::CohesiveMesh, u::Vector{Float64},E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive;F_ext::Union{Vector{Float64}, Nothing}=nothing,F_thermo_chem::Union{Vector{Float64}, Nothing}=nothing,damage_states=nothing)
    ndof = 2 * czm_mesh.nnode
    
    # 固体刚度
    K_bulk = assemble_bulk_stiffness(czm_mesh, E_eff, ν_eff)
    
    # 内聚力刚度和内力
    K_coh, f_int_coh, separations, tractions = assemble_czm_system(
        czm_mesh, u, cohesive_params; damage_states=damage_states)
    
    # 固体内力（线性弹性：f_int = K * u）
    f_int_bulk = K_bulk * u
    
    # 总刚度矩阵
    K_total = K_bulk + K_coh
    
    # 总内力 = 固体内力 + 内聚力内力
    f_int_total = f_int_bulk + f_int_coh
    
    return K_total, f_int_total, separations, tractions
end

function assemble_coupled_system_full(czm_mesh::CohesiveMesh, u::Vector{Float64},E_eff::Float64, ν_eff::Float64,α_eff::Float64, β_n::Float64, β_p::Float64,cohesive_params::Cohesive,dT_elem::Vector{Float64},Δsoc_n_elem::Vector{Float64},Δsoc_p_elem::Vector{Float64};F_ext::Union{Vector{Float64}, Nothing}=nothing,damage_states=nothing)
    ndof = 2 * czm_mesh.nnode
    
    # 组装基本系统
    K_total, f_int_total, separations, tractions = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; damage_states=damage_states)
    
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
