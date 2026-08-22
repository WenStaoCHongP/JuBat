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
    gl_element_residual_tangent(x_e, y_e, u_e, D_mat, ε0, gsorder) -> (f_e, K_e)

单个 Q4 的完全 Green-Lagrange 残差/切线（Batch 2，spec §3.2，D9）。坐标与位移均在
L 归一系（生产 `czm_mesh.node` 已 x/L 归一，u 为 L 归一位移），F = I + ∇u 无量纲。

- 应变（工程剪切约定，E_vec = [E11, E22, 2E12]，直接由 E = ½(FᵀF−I) 计算）；
- 本构 S = C:(E_vec − ε₀[1,1,0])（逐层各向同性 SVK，D_mat 同线性路径；D-B2-1）；
- f = ∫ B_GLᵀ S dA；K = ∫ B_GLᵀ C B_GL dA + ∫ Gᵀ Ŝ G dA（标准初应力 K_G，
  无手工曲率项——曲率由物理坐标网格几何携带）。
"""
function gl_element_residual_tangent(x_e, y_e, u_e::Vector{Float64},
                                     D_mat::Matrix{Float64}, ε0::Float64, gsorder::Int;
                                     plastic::Union{Nothing, Tuple}=nothing,
                                     commit_to::Union{Nothing, Vector{Tuple{NTuple{3,Float64},Float64}}}=nothing,
                                     σ0::NTuple{3, Float64}=(0.0, 0.0, 0.0))
    f_e = zeros(Float64, 8)
    K_e = zeros(Float64, 8, 8)
    gp = 0
    IntQ4(x_e, y_e; order=gsorder) do ξ, η, w, dNdx, dNdy, detJ
        gp += 1
        # ∇u
        uxx = 0.0; uxy = 0.0; uyx = 0.0; uyy = 0.0
        for i in 1:4
            uxx += dNdx[i] * u_e[2*i-1]
            uxy += dNdy[i] * u_e[2*i-1]
            uyx += dNdx[i] * u_e[2*i]
            uyy += dNdy[i] * u_e[2*i]
        end
        # 完全 GL 应变（E = ½(FᵀF − I)）
        E11 = uxx + 0.5 * (uxx * uxx + uyx * uyx)
        E22 = uyy + 0.5 * (uxy * uxy + uyy * uyy)
        g12 = (uxy + uyx) + (uxx * uxy + uyx * uyy)
        # 应力与（塑性时的）算法一致切线
        local S1::Float64, S2::Float64, S3::Float64, D_tan::Matrix{Float64}
        if plastic === nothing
            S1 = D_mat[1, 1] * (E11 - ε0) + D_mat[1, 2] * (E22 - ε0)
            S2 = D_mat[1, 2] * (E11 - ε0) + D_mat[2, 2] * (E22 - ε0)
            S3 = D_mat[3, 3] * g12
            D_tan = D_mat
        else
            σ_y, H, eps_p, κ = plastic
            e_mech = [E11 - ε0 - eps_p[gp][1], E22 - ε0 - eps_p[gp][2], g12 - eps_p[gp][3]]
            σ, C_ep, Δp, Δκ = return_mapping_plane_stress(e_mech, D_mat, σ_y, H, eps_p[gp], κ[gp])
            S1, S2, S3 = σ
            D_tan = C_ep
            if commit_to !== nothing
                commit_to[gp] = ((eps_p[gp][1] + Δp[1], eps_p[gp][2] + Δp[2], eps_p[gp][3] + Δp[3]),
                                 κ[gp] + Δκ)
            end
        end
        # 卷绕预应力 σ₀（Batch 2'，D-B2'-3：零值旁路保逐位；残差与 K_G 的 Ŝ 均用总应力）
        if σ0 != (0.0, 0.0, 0.0)
            S1 += σ0[1]; S2 += σ0[2]; S3 += σ0[3]
        end
        # B_GL = ∂E_vec/∂u_e（u=0 时退化为线性 B）
        B = zeros(Float64, 3, 8)
        for i in 1:4
            dx, dy = dNdx[i], dNdy[i]
            B[1, 2*i-1] = (1.0 + uxx) * dx
            B[1, 2*i]   = uyx * dx
            B[2, 2*i-1] = uxy * dy
            B[2, 2*i]   = (1.0 + uyy) * dy
            B[3, 2*i-1] = (1.0 + uxx) * dy + uxy * dx
            B[3, 2*i]   = (1.0 + uyy) * dx + uyx * dy
        end
        # 残差 f = ∫ Bᵀ S dA
        BtS = B' * [S1, S2, S3]
        axpy!(w * detJ, BtS, f_e)
        # 材料切线 ∫ Bᵀ C B dA（塑性时为算法一致切线 C_ep）
        DB = D_tan * B
        K_e .+= (w * detJ) .* (B' * DB)
        # 标准初应力 K_G = ∫ Gᵀ Ŝ G dA
        G = zeros(Float64, 4, 8)
        for i in 1:4
            G[1, 2*i-1] = dNdx[i]
            G[2, 2*i-1] = dNdy[i]
            G[3, 2*i]   = dNdx[i]
            G[4, 2*i]   = dNdy[i]
        end
        Sh = [S1 S3 0.0 0.0; S3 S2 0.0 0.0; 0.0 0.0 S1 S3; 0.0 0.0 S3 S2]
        K_e .+= (w * detJ) .* (G' * Sh * G)
    end
    return f_e, K_e
end

"""
    assemble_bulk_residual_tangent(czm_mesh, u, param_cache, mech_state=nothing;
                                  geo_nl=false, plasticity=false, K_bulk_cached=nothing,
                                  eigenstrain=nothing)
        -> (f_int_bulk, K_tangent)

bulk 残差/切线的统一入口（spec 2026-08-20-core-collapse-mechanics-design.md §4.2）。

三个槽位，前两个已实现：

1. **线弹性**（`geo_nl=false, plasticity=false`）：`f_int = K_bulk*u`，`K_tangent = K_bulk`，
   与既有 `assemble_bulk_stiffness` 路径逐位等价。
2. **几何非线性**（`geo_nl=true`，Batch 2）：完全 Green-Lagrange 全 Lagrangian——
   `gl_element_residual_tangent` 逐单元装配，S = C:(E_GL − ε₀I)（D-B2-1，ε₀ 与
   `assemble_thermal_chemical_load` 同式）；切线含标准初应力 K_G。切线依赖 u，
   禁止传 `K_bulk_cached`。`eigenstrain` 为 NamedTuple `(α_eff, β_n, β_p, dT, Δsn, Δsp)`，
   `nothing` 表示本次调用无本征应变（ε₀≡0，运动学场景合法状态）。
3. **J2 塑性**（`plasticity=true`）：PCC/NCC 平面应力一致返回映射，Batch 3；届时经
   `mech_state` 传入 `PlasticState`。

未实现的槽位传入非默认值一律 `error`——静默走线弹性会让上层误以为已生效
（AGENTS 9.7）。

`mech_state` 按 spec §4.2 保留为尾置可选位置参数，使 Batch 3 引入塑性状态时无需改签名。
"""
function assemble_bulk_residual_tangent(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param_cache::CzmParamCache,
    mech_state=nothing;
    geo_nl::Bool=false,
    plasticity::Bool=false,
    K_bulk_cached::Union{Nothing, SparseMatrixCSC{Float64, Int64}}=nothing,
    eigenstrain=nothing,
    commit_plastic::Bool=false,
    prestress=nothing
)
    ndof = 2 * czm_mesh.nnode
    length(u) == ndof || throw(DimensionMismatch(
        "assemble_bulk_residual_tangent: u 长度为 $(length(u))，应为 $ndof " *
        "(2 × nnode=$(czm_mesh.nnode))"))

    if plasticity && !geo_nl
        error("assemble_bulk_residual_tangent: plasticity=true 需要 geo_nl=true（D-B3-1：塑性逐 GP " *
              "应力评估仅在 GL 单元路径实现，线弹性常刚度路径无消费者）。")
    end
    if !plasticity && mech_state !== nothing
        error("assemble_bulk_residual_tangent: mech_state 的消费者为 plasticity（Batch 3 起）；" *
              "未开塑性时必须传 nothing，收到 $(typeof(mech_state))。")
    end

    if prestress !== nothing
        geo_nl || error(
            "assemble_bulk_residual_tangent: prestress 需要 geo_nl=true（D-B2'-2：初应力由 GL 残差与 K_G 消费）。")
        length(prestress) == size(czm_mesh.bulk_element, 1) || throw(DimensionMismatch(
            "assemble_bulk_residual_tangent: prestress 长度 $(length(prestress)) 应为 bulk 单元数 $(size(czm_mesh.bulk_element, 1))"))
    end
    if geo_nl
        K_bulk_cached !== nothing && error(
            "assemble_bulk_residual_tangent: geo_nl=true 时切线依赖 u，" *
            "不得传 K_bulk_cached（须逐迭代重组）。")
        if plasticity
            # D-B3-1：塑性逐 GP 应力评估仅在 GL 单元路径实现
            mech_state isa Matrix{PlasticState} || error(
                "assemble_bulk_residual_tangent: plasticity=true 需要 mech_state::Matrix{PlasticState}（[ne, 4] 高斯点状态），" *
                "收到 $(typeof(mech_state))。")
        end
        α_eff = eigenstrain === nothing ? 0.0 : eigenstrain.α_eff
        β_n   = eigenstrain === nothing ? 0.0 : eigenstrain.β_n
        β_p   = eigenstrain === nothing ? 0.0 : eigenstrain.β_p
        dT_el = eigenstrain === nothing ? nothing : eigenstrain.dT
        Δsn   = eigenstrain === nothing ? nothing : eigenstrain.Δsn
        Δsp   = eigenstrain === nothing ? nothing : eigenstrain.Δsp
        param = param_cache.param_ref
        submesh = czm_mesh.czm_submesh
        element = czm_mesh.bulk_element
        node = czm_mesh.node
        ne0 = size(element, 1)
        if dT_el !== nothing
            (length(dT_el) == ne0 && length(Δsn) == ne0 && length(Δsp) == ne0) ||
                throw(DimensionMismatch(
                    "assemble_bulk_residual_tangent: eigenstrain 向量长度应为 bulk 单元数 $ne0"))
        end
        if plasticity && size(mech_state, 1) != ne0
            throw(DimensionMismatch(
                "assemble_bulk_residual_tangent: mech_state 行数 $(size(mech_state,1)) 应为 bulk 单元数 $ne0"))
        end
        I_idx = Int64[]; J_idx = Int64[]; K_vals = Float64[]
        sizehint!(I_idx, ne0 * 64); sizehint!(J_idx, ne0 * 64); sizehint!(K_vals, ne0 * 64)
        f_gl = zeros(Float64, ndof)
        for e in 1:ne0
            mt = submesh.material_type[e]
            local D_mat::Matrix{Float64}, plastic::Union{Nothing, Tuple}
            # 模量统一走 moduli_of（E 为 ÷E_coat 归一，×E_coat/σ_czm 链到 σ_czm 系；
            # 用户参数修正后 PCC.E/NCC.E 即物理箔模量，塑性与弹性路径共用）
            E_e, ν_e = moduli_of(param, mt)
            D_mat = E_e / (1.0 - ν_e^2) * [1.0 ν_e 0.0;
                                          ν_e 1.0 0.0;
                                          0.0 0.0 (1.0 - ν_e) / 2.0]
            plastic = nothing
            if plasticity && (mt === :PCC || mt === :NCC)
                σ_y, H = foil_params_of(param, mt)
                σ_y > 0.0 || error(
                    "assemble_bulk_residual_tangent: czm_j2_plasticity=true 但 $mt 的 sigma_y ≤ 0（未设置）。" *
                    "缺参即拦截，不默认、不置零（AGENTS 9.4/9.7）。")
                plastic = (σ_y, H,
                           NTuple{3, Float64}[mech_state[e, g].eps_p for g in 1:4],
                           Float64[mech_state[e, g].kappa for g in 1:4])
            end
            ε0 = 0.0
            if dT_el !== nothing
                ε0 = α_eff * dT_el[e] + β_n * Δsn[e] + β_p * Δsp[e]
            end
            commit_to = (plasticity && commit_plastic && plastic !== nothing) ?
                        Vector{Tuple{NTuple{3,Float64},Float64}}(undef, 4) : nothing
            elem_nodes = element[e, :]
            x_e = node[elem_nodes, 1]
            y_e = node[elem_nodes, 2]
            u_e = zeros(Float64, 8)
            for (k, n) in enumerate(elem_nodes)
                u_e[2*k-1] = u[2*n-1]
                u_e[2*k] = u[2*n]
            end
            f_e, K_e = gl_element_residual_tangent(x_e, y_e, u_e, D_mat, ε0, 2;
                                                   plastic=plastic, commit_to=commit_to,
                                                   σ0=prestress === nothing ? (0.0, 0.0, 0.0) : prestress[e])
            if commit_to !== nothing
                for g in 1:4
                    mech_state[e, g].eps_p = commit_to[g][1]
                    mech_state[e, g].kappa = commit_to[g][2]
                end
            end
            dofs = Int64[]
            for n in elem_nodes
                push!(dofs, 2*n - 1)
                push!(dofs, 2*n)
            end
            for a in 1:8
                f_gl[dofs[a]] += f_e[a]
                for b in 1:8
                    push!(I_idx, dofs[a])
                    push!(J_idx, dofs[b])
                    push!(K_vals, K_e[a, b])
                end
            end
        end
        return f_gl, sparse(I_idx, J_idx, K_vals, ndof, ndof)
    else
        K_tangent = K_bulk_cached !== nothing ? K_bulk_cached :
                    assemble_bulk_stiffness(czm_mesh, param_cache)
        f_int_bulk = K_tangent * u
        return f_int_bulk, K_tangent
    end
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
    geo_nl::Bool=false,
    eigenstrain=nothing,
    plasticity::Bool=false,
    mech_state=nothing,
    commit_plastic::Bool=false,
    prestress=nothing
)
    ndof = 2 * czm_mesh.nnode

    # 固体残差与切线（统一入口，spec §4.2）。geo_nl=false 为 Batch 1 线弹性槽位
    # （与 K_bulk*u 逐位等价）；geo_nl=true 走完全 GL（Batch 2），ε* 内嵌（D-B2-1）；
    # plasticity=true 时 PCC/NCC 用物理箔 J2（Batch 3，D-B3-0）；
    # prestress 叠加卷绕预应力 σ₀ 进残差与 K_G（Batch 2'，D-B2'-1/2）。
    f_int_bulk, K_bulk = assemble_bulk_residual_tangent(
        czm_mesh, u, param_cache, mech_state; K_bulk_cached=K_bulk_cached,
        geo_nl=geo_nl, eigenstrain=eigenstrain,
        plasticity=plasticity, commit_plastic=commit_plastic,
        prestress=prestress)

    # 内聚力刚度和内力（使用几何缓存和工作区，透传 param_cache）
    K_coh, f_int_coh, separations, tractions = assemble_czm_system(
        czm_mesh, u, param_cache; damage_states=damage_states,
        geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)

    # 总刚度矩阵
    K_total = K_bulk + K_coh

    # 总内力 = 固体内力 + 内聚力内力
    f_int_total = f_int_bulk + f_int_coh

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
