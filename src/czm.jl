# CZM damage state, material lookup, and system assembly.

"""
    DamageState - 内聚力单元的损伤状态

存储每个内聚力单元（或高斯点）的损伤历史。

# 字段
- `D`: 等效损伤 D_eq [0, 1]，由当前双线性律计算
- `D_visc`: 粘性有效损伤，用于牵引力和切线（粘性正则化关闭时 D_visc = D）
- `δ_max_n`: 历史最大法向分离位移
- `δ_max_t`: 历史最大切向分离位移
- `δ_max_eff`: 历史最大等效分离位移
- `fractured`: 是否已完全断裂
- `accumulated_damage`: 累积损伤（用于循环加载）
"""
mutable struct DamageState <: AbstractDamageState
    D::Float64                     # 等效损伤 D_eq
    D_visc::Float64                # 粘性有效损伤（用于牵引和切线）
    δ_max_n::Float64              # 历史最大法向分离
    δ_max_t::Float64              # 历史最大切向分离
    δ_max_eff::Float64            # 历史最大等效分离
    fractured::Bool               # 是否断裂
    accumulated_damage::Float64   # 累积损伤（循环）

    # 默认构造函数
    DamageState() = new(0.0, 0.0, 0.0, 0.0, 0.0, false, 0.0)
end

# ========================================================================
# 2. 系统组装
# ========================================================================

"""
    moduli_of(param, mt::Symbol) -> (E, ν)

按材料类型从 param 读取体模量、泊松比，并**统一到 CZM 应力空间（σ_czm 参考）**。
供 assemble_bulk_stiffness / assemble_thermal_chemical_load 复用。

材料类型对应 CzmSubmesh.material_type 中的 Symbol：
- :PE  → 涂层模量 param.PE.E_coat / param.PE.nu_coat
- :NE  → 涂层模量 param.NE.E_coat / param.NE.nu_coat
- :SP  → param.SP.E / param.SP.nu
- :PCC → param.PCC.E / param.PCC.nu
- :NCC → param.NCC.E / param.NCC.nu

双重再缩放（重设计 v2 §3）：模量字段在 NormaliseParam 中以 scale.E_coat 归一，
而 CZM 牵引-分离律以 scale.σ_czm 归一。体刚度与内聚力刚度装配到同一残差，
必须共享应力参考，故此处乘 `scale.E_coat / scale.σ_czm` 转到 σ_czm 空间
（与 CzmInterfaceParams.E_eff 的构造一致）。

注意：α 已从此函数移除（I2-a 修复）。两个调用者均不使用 α，且
SP/PCC/NCC.alphaT 字段在 Jellyroll.jl 中未设置，silently 取 0 易踩坑。
如未来热-化学载荷需要 α，应显式新建 ``alpha_of(param, mt)`` helper。
"""
function moduli_of(param, mt::Symbol)
    s = param.scale.E_coat / param.scale.σ_czm
    mt === :PE  && return (param.PE.E_coat * s,  param.PE.nu_coat)
    mt === :NE  && return (param.NE.E_coat * s,  param.NE.nu_coat)
    mt === :SP  && return (param.SP.E * s,       param.SP.nu)
    mt === :PCC && return (param.PCC.E * s,      param.PCC.nu)
    mt === :NCC && return (param.NCC.E * s,      param.NCC.nu)
    error("moduli_of: unknown material_type $mt")
end


# 列内二分查找稀疏矩阵 (i,j) 的 nzval 下标；不存在返回 0（仅 pattern 初始化时使用）
function _nz_index(K::SparseMatrixCSC, i::Int, j::Int)
    lo, hi = K.colptr[j], K.colptr[j+1] - 1
    while lo <= hi
        mid = (lo + hi) >>> 1
        if K.rowval[mid] == i; return mid
        elseif K.rowval[mid] < i; lo = mid + 1
        else; hi = mid - 1; end
    end
    return 0
end

"""
    assemble_czm_system(czm_mesh, u, param_cache; damage_states=nothing, ...)

组装内聚力单元的全局刚度矩阵和内力向量。
按 cohesive 单元 `interface_type` 从 `param_cache.by_interface` 取 `CzmInterfaceParams`。

# 返回
- `K_coh`: 内聚力刚度矩阵 (ndof × ndof)
- `f_int_coh`: 内聚力内力向量 (ndof,)
- `separations`: 每个单元的分离位移
- `tractions`: 每个单元的牵引力
"""
function assemble_czm_system(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param_cache::CzmParamCache;
    damage_states=nothing,
    geom_cache::Union{Nothing, Vector{CohesiveElementGeom}}=nothing,
    ws::Union{Nothing, CZMAssemblyWorkspace}=nothing,
    visc_beta::Float64=1.0
)
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
        # pattern 固定，同步缓存每单元 64 项的 nzval 下标（与组装循环 a/b 顺序一一对应），
        # 之后组装用 K_coh.nzval[idx] += 直加，替代标量稀疏索引的查找开销（数值路径不变）
        nzidx = zeros(Int, n_coh, 64)
        for i in 1:n_coh
            if geom_cache !== nothing
                dofs_ = geom_cache[i].dofs
            else
                elem = czm_mesh.cohesive_elements[i]
                n1, n2 = elem.nodes_bottom
                n4, n3 = elem.nodes_top
                dofs_ = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3, 2*n4-1, 2*n4]
            end
            c = 0
            for a in 1:8, b in 1:8
                c += 1
                nzidx[i, c] = _nz_index(K_coh, dofs_[a], dofs_[b])
            end
        end
        ws.cohesive_nzidx = nzidx
    end

    # 清零 nonzero 值（O(nnz)，无分配）
    fill!(nonzeros(K_coh), 0.0)

    @inbounds for i in 1:n_coh
        damage_state = states[i]

        # 按 interface_type 从 param_cache 取本构参数（spec §7.1）
        iface = czm_mesh.cohesive_elements[i].interface_type
        params = param_cache.by_interface[iface]
        # Λ：位移空间（L 归一）→ 分离空间（δ_czm 归一）换算因子（重设计 v2 §5）。
        # 虚功一致性：δ̃ = Λ·B·ũ；内力 f = ∫BᵀT̃ dΓ 不乘 Λ；切线刚度乘一次 Λ。
        Λ = params.Λ

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
            L >= 1e-15 || error("degenerate cohesive element $i: tangential length is $L")
            R = geom.R
            dofs = geom.dofs
            wts = geom.gauss_wts
            pts = geom.gauss_pts
        else
            elem = czm_mesh.cohesive_elements[i]
            n1, n2 = elem.nodes_bottom
            n4, n3 = elem.nodes_top
            L = elem.length
            L >= 1e-15 || error("degenerate cohesive element $i: tangential length is $L")

            x1, y1 = czm_mesh.node[n1, 1], czm_mesh.node[n1, 2]
            x2, y2 = czm_mesh.node[n2, 1], czm_mesh.node[n2, 2]
            dx, dy = x2 - x1, y2 - y1
            t_vec = [dx / L, dy / L]
            n_vec = [-t_vec[2], t_vec[1]]
            R = [n_vec[1] n_vec[2]; t_vec[1] t_vec[2]]
            dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3, 2*n4-1, 2*n4]
            order = czm_mesh.bulk_mesh.gs.order
            wts, pts = NCweight(order)
        end

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
                # 换算到分离空间（δ_czm 归一），与本构参数 δ_0*/δ_c*/K* 同参考
                δ_n = Λ * ws.δ_local[1]
                δ_t = Λ * ws.δ_local[2]

                T_n, T_t, _, _ = bilinear_traction_state(δ_n, δ_t, damage_state, params; visc_beta=visc_beta)
                dT_dδ = bilinear_tangent(δ_n, δ_t, damage_state, params; visc_beta=visc_beta)

                J = L / 2.0
                wJ = w * J

                # BL_dT = B_local' * dT_dδ  （mul! 无分配）
                mul!(ws.BL_dT, transpose(ws.B_local), dT_dδ)

                # K_e += wJ * Λ * BL_dT * B_local  （mul! 无分配）
                # dT̃/dũ = (dT̃/dδ̃)·Λ·B —— 切线刚度含一次 Λ（重设计 v2 §5 式(3)）
                mul!(ws.BL_dT_B, ws.BL_dT, ws.B_local)
                wJΛ = wJ * Λ
                for a in 1:8
                    for b in 1:8
                        ws.K_e[a, b] += wJΛ * ws.BL_dT_B[a, b]
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

        if w_sum > 0.0
            T_n_avg /= w_sum
            T_t_avg /= w_sum
            δ_n_avg /= w_sum
            δ_t_avg /= w_sum
        end

        ws.separations[i] = (δ_n_avg, δ_t_avg)
        ws.tractions[i] = (T_n_avg, T_t_avg)

        # 组装到预分配稀疏矩阵（nzval 直索引：getindex/setindex 语义 = 读存储值→加→写回，
        # 与原 K_coh[dofs[a], dofs[b]] += 逐位一致；无 nzidx 时回退原写法）
        if ws.cohesive_nzidx !== nothing
            c = 0
            for a in 1:8
                ws.f_int_coh[dofs[a]] += ws.f_int_e[a]
                for b in 1:8
                    c += 1
                    K_coh.nzval[ws.cohesive_nzidx[i, c]] += ws.K_e[a, b]
                end
            end
        else
            for a in 1:8
                ws.f_int_coh[dofs[a]] += ws.f_int_e[a]
                for b in 1:8
                    K_coh[dofs[a], dofs[b]] += ws.K_e[a, b]
                end
            end
        end
    end

    return K_coh, ws.f_int_coh, ws.separations, ws.tractions
end

"""
    assemble_bulk_stiffness(czm_mesh, param_cache)

组装固体单元（Q4）的刚度矩阵。按 `czm_submesh.material_type` 分组取
体模量（PE/NE 用 E_coat，SP/PCC/NCC 用连续层 E），不再使用全栈均一模量。

# 返回
- `K_bulk`: 固体刚度矩阵 (ndof × ndof)
"""
function assemble_bulk_stiffness(czm_mesh::CohesiveMesh, param_cache::CzmParamCache)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    param = param_cache.param_ref
    submesh = czm_mesh.czm_submesh
    submesh === nothing && error(
        "assemble_bulk_stiffness: czm_submesh is nothing " *
        "(must be built via jellyroll_collector_seed_mesh with czm_enabled=true)")

    # 需要重新计算高斯积分点，因为节点可能已经改变
    element = czm_mesh.bulk_element
    node = czm_mesh.node
    ne = size(element, 1)
    gsorder = 2

    # 稀疏矩阵组装
    I_idx = Int64[]
    J_idx = Int64[]
    K_vals = Float64[]

    for e in 1:ne
        # 按材料类型查表（PE/NE 用涂层模量，SP/PCC/NCC 用连续层模量）
        E_e, ν_e = moduli_of(param, submesh.material_type[e])

        # 弹性矩阵（平面应力）
        D_mat = E_e / (1.0 - ν_e^2) * [1.0 ν_e 0.0;
                                       ν_e 1.0 0.0;
                                       0.0 0.0 (1.0-ν_e)/2.0]

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

"""
    assemble_thermal_chemical_load(czm_mesh, param_cache, α_eff, β_n, β_p, dT_elem, Δsoc_n_elem, Δsoc_p_elem)

组装热-化学载荷向量。按 `czm_submesh.material_type` 分组取 (E, ν)（PE/NE 用
涂层模量，SP/PCC/NCC 用连续层模量），保留 α_eff / β_n / β_p 为位置参数
（电化学浓度膨胀系数，跨材料统一）。
"""
function assemble_thermal_chemical_load(
    czm_mesh::CohesiveMesh,
    param_cache::CzmParamCache,
    α_eff::Float64, β_n::Float64, β_p::Float64,
    dT_elem::Vector{Float64}, Δsoc_n_elem::Vector{Float64}, Δsoc_p_elem::Vector{Float64}
)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    param = param_cache.param_ref
    submesh = czm_mesh.czm_submesh
    submesh === nothing && error(
        "assemble_thermal_chemical_load: czm_submesh is nothing " *
        "(must be built via jellyroll_collector_seed_mesh with czm_enabled=true)")
    element = czm_mesh.bulk_element
    node = czm_mesh.node
    ne = size(element, 1)

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
        # 按材料类型查表（PE/NE 用涂层模量，SP/PCC/NCC 用连续层模量）
        E_e, ν_e = moduli_of(param, submesh.material_type[e])

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
            factor = E_e / (1.0 - ν_e^2) * ε_0 * (1.0 + ν_e) * w * detJ

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
    build_czm_cache(czm_mesh, param_cache; fix_inner=true)

构建 CZM 装配缓存，包括 K_bulk、cohesive 几何、边界条件。

失效判据由 `czm_mesh_id = objectid(czm_mesh)` 与 `param_cache_id = param_cache.id`
共同决定（见 `ensure_czm_cache`），整个缓存可在多次 Newton 迭代中复用，
只要 mesh 对象与 param_cache 对象不变。

# 参数
- `czm_mesh::CohesiveMesh`: CZM 网格
- `param_cache::CzmParamCache`: per-interface 参数缓存（提供 `param_ref` 给
  `assemble_bulk_stiffness` 与 `identify_bc_nodes_czm`）
- `fix_inner::Bool=true`: 是否固定内圈节点（影响 BC 构造）

# 返回
- `CZMAssemblyCache`: 填充好的缓存，挂载到 `case.czm_cache` 上跨步复用
"""
function build_czm_cache(czm_mesh::CohesiveMesh, param_cache::CzmParamCache; fix_inner::Bool=true)
    param = param_cache.param_ref
    cache = CZMAssemblyCache()

    # 1. 缓存 K_bulk（最高 ROI：消除 ~60 次/更新 的冗余 bulk 重组）
    cache.K_bulk = assemble_bulk_stiffness(czm_mesh, param_cache)

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
        L >= 1e-15 || error("degenerate cohesive element $i: tangential length is $L")

        x1, y1 = czm_mesh.node[n1, 1], czm_mesh.node[n1, 2]
        x2, y2 = czm_mesh.node[n2, 1], czm_mesh.node[n2, 2]
        dx, dy = x2 - x1, y2 - y1
        t_vec = [dx / L, dy / L]
        n_vec = [-t_vec[2], t_vec[1]]
        R = [n_vec[1] n_vec[2]; t_vec[1] t_vec[2]]

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
    bc_nodes, _, _ = identify_bc_nodes_czm(czm_mesh, param; fix_inner=fix_inner)
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

    # 5. 失效判据标记
    cache.fix_inner = fix_inner
    cache.czm_mesh_id = objectid(czm_mesh)
    cache.param_cache_id = param_cache.id

    # 6. 创建可复用工作区（ndof × n_coh，跨时间步复用）
    ndof = 2 * czm_mesh.nnode
    cache.ws = CZMAssemblyWorkspace(ndof, n_coh)

    cache.valid = true

    return cache
end

"""
    ensure_czm_cache(case, czm_mesh, param_cache; fix_inner=true)

确保 `case.czm_cache` 可用且未过期。失效条件（任一触发即重建）：
1. `cache === nothing` 或 `!cache.valid`
2. `cache.czm_mesh_id != objectid(czm_mesh)`：mesh 对象变了
3. `cache.param_cache_id != param_cache.id`：param_cache 对象变了
4. `cache.fix_inner != fix_inner`：BC 配置切换

判据：`czm_mesh_id` 用 `objectid(czm_mesh)` 检测网格对象替换；`param_cache_id` 用
`param_cache.id`（`compute_czm_params_per_interface` 计算的内容哈希）检测参数内容变化。
内容修改后重新调用 `compute_czm_params_per_interface` 拿到新对象即可触发失效。
"""
function ensure_czm_cache(case::Case, czm_mesh::CohesiveMesh, param_cache::CzmParamCache; fix_inner::Bool=true)
    cache = case.czm_cache
    if cache === nothing || !cache.valid ||
       cache.czm_mesh_id != objectid(czm_mesh) ||
       cache.param_cache_id != param_cache.id ||
       cache.fix_inner != fix_inner
        cache = build_czm_cache(czm_mesh, param_cache; fix_inner=fix_inner)
        case.czm_cache = cache
    end
    return cache
end

"""
    assemble_coupled_system(czm_mesh, u, param_cache; F_ext=nothing, ...)

组装耦合系统（体刚度 + 内聚力）。签名按 spec v2 §7.1 改为接受 `param_cache`，
体内刚度按 `czm_submesh.material_type` 分组取模量。
"""
function assemble_coupled_system(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param_cache::CzmParamCache;
    F_ext::Union{Vector{Float64}, Nothing}=nothing,
    F_thermo_chem::Union{Vector{Float64}, Nothing}=nothing,
    damage_states=nothing,
    K_bulk_cached::Union{Nothing, SparseMatrixCSC{Float64, Int64}}=nothing,
    geom_cache::Union{Nothing, Vector{CohesiveElementGeom}}=nothing,
    ws::Union{Nothing, CZMAssemblyWorkspace}=nothing,
    visc_beta::Float64=1.0,
    assemble_K::Bool=true
)
    ndof = 2 * czm_mesh.nnode

    # 固体刚度（使用缓存或重新计算）
    K_bulk = K_bulk_cached !== nothing ? K_bulk_cached : assemble_bulk_stiffness(czm_mesh, param_cache)

    # 内聚力刚度和内力（使用几何缓存和工作区，透传 param_cache）
    K_coh, f_int_coh, separations, tractions = assemble_czm_system(
        czm_mesh, u, param_cache; damage_states=damage_states,
        geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)

    # 固体内力（预分配 mul!：稀疏 matvec 与 K_bulk * u 同一实现路径）
    if ws !== nothing
        if ws.f_int_bulk_buf === nothing || length(ws.f_int_bulk_buf) != ndof
            ws.f_int_bulk_buf = zeros(Float64, ndof)
        end
        mul!(ws.f_int_bulk_buf, K_bulk, u)
        f_int_bulk = ws.f_int_bulk_buf
    else
        f_int_bulk = K_bulk * u
    end

    # 总内力 = 固体内力 + 内聚力内力（copyto! + .+= 与向量 + 的逐元素加同序，逐位一致）
    if ws !== nothing
        if ws.f_int_total_buf === nothing || length(ws.f_int_total_buf) != ndof
            ws.f_int_total_buf = zeros(Float64, ndof)
        end
        copyto!(ws.f_int_total_buf, f_int_bulk)
        ws.f_int_total_buf .+= f_int_coh
        f_int_total = ws.f_int_total_buf
    else
        f_int_total = f_int_bulk + f_int_coh
    end

    # 总刚度矩阵（assemble_K=false 时跳过——线搜索等只消费 f_int 的调用点专用）
    if assemble_K
        if ws !== nothing && ws.K_total_buf !== nothing
            # 预分配同序加法：K_total 每存储位置恰一次 K_bulk 值 + K_coh 值（等价性实验已验证
            # 与 SparseArrays + 的合并语义逐位一致）
            nzT = nonzeros(ws.K_total_buf)
            mapK = ws.K_total_mapK
            mapC = ws.K_total_mapC
            for k in 1:length(nzT)
                nzT[k] = (mapK[k] > 0 ? K_bulk.nzval[mapK[k]] : 0.0) +
                         (mapC[k] > 0 ? K_coh.nzval[mapC[k]] : 0.0)
            end
            K_total = ws.K_total_buf
        else
            K_total = K_bulk + K_coh
            if ws !== nothing
                # 首次：缓存 K_total pattern 与 K_bulk/K_coh nzval 下标映射
                ws.K_total_buf = copy(K_total)
                mapK = zeros(Int, nnz(K_total))
                mapC = zeros(Int, nnz(K_total))
                for col in 1:size(K_total, 2)
                    for k in K_total.colptr[col]:(K_total.colptr[col+1]-1)
                        mapK[k] = _nz_index(K_bulk, K_total.rowval[k], col)
                        mapC[k] = _nz_index(K_coh, K_total.rowval[k], col)
                    end
                end
                ws.K_total_mapK = mapK
                ws.K_total_mapC = mapC
            end
        end
    else
        K_total = nothing
    end

    return K_total, f_int_total, separations, tractions
end

"""
    assemble_coupled_system_full(czm_mesh, u, param_cache, α_eff, β_n, β_p, dT_elem, ...; kwargs...)

耦合系统组装 + 热-化学载荷 + 残差计算。按 spec v2 §7.1 改为接受 `param_cache`，
透传给 `assemble_coupled_system` 与 `assemble_thermal_chemical_load`。
"""
function assemble_coupled_system_full(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param_cache::CzmParamCache,
    α_eff::Float64, β_n::Float64, β_p::Float64,
    dT_elem::Vector{Float64}, Δsoc_n_elem::Vector{Float64}, Δsoc_p_elem::Vector{Float64};
    F_ext::Union{Vector{Float64}, Nothing}=nothing,
    damage_states=nothing,
    K_bulk_cached::Union{Nothing, SparseMatrixCSC{Float64, Int64}}=nothing,
    geom_cache::Union{Nothing, Vector{CohesiveElementGeom}}=nothing,
    ws::Union{Nothing, CZMAssemblyWorkspace}=nothing,
    visc_beta::Float64=1.0
)
    ndof = 2 * czm_mesh.nnode

    # 组装基本系统（透传 param_cache）
    K_total, f_int_total, separations, tractions = assemble_coupled_system(
        czm_mesh, u, param_cache;
        damage_states=damage_states, K_bulk_cached=K_bulk_cached,
        geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)

    # 热-化学载荷（透传 param_cache）
    F_thermo_chem = assemble_thermal_chemical_load(
        czm_mesh, param_cache, α_eff, β_n, β_p,
        dT_elem, Δsoc_n_elem, Δsoc_p_elem)

    # 外部载荷
    F_external = F_ext === nothing ? zeros(Float64, ndof) : F_ext

    # 残差 = 外力 + 热化学力 - 内力
    R = F_external + F_thermo_chem - f_int_total

    return K_total, R, F_thermo_chem, separations, tractions
end
