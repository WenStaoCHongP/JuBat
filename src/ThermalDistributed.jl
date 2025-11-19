"""
    ThermalDistributed1D(case, variables)
框架：装配一维分布式热传导方程的 M/K/F（暂用占位，后续接入具体网格/源项）
"""
function ThermalDistributed1D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
    # Mesh & params (1D thermal mesh 复用现有 1D Mesh 结构)
    @assert haskey(case.mesh, "thermal1D") "thermal1D mesh is missing in case.mesh"
    mesh = case.mesh["thermal1D"]
    nnode = mesh.nlen
    # consistent mass matrix: ∫ ρc NᵀN dS ；刚度：∫ k dNᵀ/dx dN/dx dS
    ρc = case.param.cell.rho * case.param.cell.heat_Q
    k  = case.param.PE.lambda  # 初版占位：后续按区域/等效处理
    # 组装框架（占位稀疏矩阵/向量）
    MT = spzeros(nnode, nnode)
    KT = spzeros(nnode, nnode)
    FT = zeros(Float64, nnode)
    return MT, KT, FT
end
"""
    ThermalDistributed2D(case, variables)

 果冻卷热传导装配：
  方程: (ρ c) ∂T/∂t = ∇·(k ∇T) + q
  
  **重要约定**：为与电化学部分统一，返回的 KT 已包含负号
  即：M dT/dt = KT T + F （其中 KT = -∫ k ∇N^T ∇N dΩ）
  
  极坐标: ∇·(k ∇T) = 1/r ∂/∂r(r k_rr ∂T/∂r) + 1/r^2 k_θθ ∂²T/∂θ² (忽略耦合项, k 对角)。
  实现策略：
    1. 使用参数域 Q4 元素 (r,θ)；在每个 Gauss 点将 (∂N/∂r, ∂N/∂θ) → 笛卡尔梯度形式, 但对各向同性 k 时可直接在极坐标形式简化。
    2. 此处采用 jellyroll 几何工具提供的 r、θ；假设各向异性 。
    3. 质量矩阵采用一致 (consistent) Nᵀ N ρc r w h。
    4. 内热源 q 已映射为节点值 (variables["heat source total"])；装配右端 F = ∫ N q r w h。
    5. 边界条件: 
        - r = Rout 外边界：对流 -k ∂T/∂n = h (T - T_amb)  ⇒ 贡献到刚度 (h) 与载荷 (-h T_amb)。
        - r = Rin 内侧：若存在冷却 (cell.h_inner>0) 同理；否则绝热。
    6. 返回 (MT,KT,FT) 。
"""
function ThermalDistributed2D(case::Case, variables::Dict{String,Union{Array{Float64},Float64}})
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh is missing"
    mesh = case.mesh["thermal2D"]
    @assert mesh.type == "Q4" "ThermalDistributed2D currently assumes Q4 mesh"
    nnode = mesh.nlen

    # Thermal Scheme B scales
    scale = case.param_dim.scale
    ρc_ref = scale.rho_c_th
    k_ref  = scale.k_th
    q_ref  = scale.q_th
    L_th   = scale.L_th

     # 高斯点数据与局部形函数
    Ni   = mesh.gs.Ni                    # (ngs x 4)
    dNdx = mesh.gs.dNidx[:, 1:4]        # (ngs x 4) dN/dx
    dNdy = mesh.gs.dNidx[:, 5:8]        # (ngs x 4) dN/dy
    wJ   = mesh.gs.weight .* mesh.gs.detJ  # (ngs)
    ne   = size(mesh.element, 1)
    ngs  = size(Ni, 1)

    # 将每个高斯点映射到对应的全局节点索引表 (Vi, Vj)
    Vi = mesh.element[mesh.gs.ele, :]    # (ngs x 4)
    Vj = Vi                              # 对称装配

    # 质量矩阵 MT = ∬ (ρc/ρc_ref) NᵀN dΩ*
    # 优先使用元素层权重 f_k 聚合 (ρc)_eff；无则回退到全局 cell 常数
    ρc_weights = similar(wJ)
    # 正确做法是从字典中取出层权重矩阵 fks。
    fks = haskey(variables, "thermal2D layer_weights") 
    if fks !== nothing && isa(fks, AbstractArray)  # accept Matrix{Float64}
        ne = size(mesh.element, 1)
        # 各层 (ρc)* = (ρ·c)/ρc_ref（无量纲）
        ρc_NE = (getfield(case.param_dim.NE, :rho)  * getfield(case.param_dim.NE, :heat_Q))  / ρc_ref
        ρc_SP = (getfield(case.param_dim.SP, :rho)  * getfield(case.param_dim.SP, :heat_Q))  / ρc_ref
        ρc_PE = (getfield(case.param_dim.PE, :rho)  * getfield(case.param_dim.PE, :heat_Q))  / ρc_ref
        ρc_PCC= (getfield(case.param_dim.PCC,:rho)  * getfield(case.param_dim.PCC,:heat_Q))  / ρc_ref
        ρc_NCC= (getfield(case.param_dim.NCC,:rho)  * getfield(case.param_dim.NCC,:heat_Q))  / ρc_ref
        ρc_e = zeros(Float64, ne)
        @inbounds for e in 1:ne
            f_NE, f_SP, f_PE, f_PCC, f_NCC = fks[e,1], fks[e,2], fks[e,3], fks[e,4], fks[e,5]
            ρc_e[e] = f_NE*ρc_NE + f_SP*ρc_SP + f_PE*ρc_PE + f_PCC*ρc_PCC + f_NCC*ρc_NCC
        end
        # dΩ* = dΩ / L_th^2
        ρc_weights .= ρc_e[mesh.gs.ele] .* (wJ ./ (L_th^2))
    else
        # fallback: use cell-averaged (ρc)* and apply area scaling
        ρc_cell_nd = (case.param_dim.cell.rho * case.param_dim.cell.heat_Q) / ρc_ref
        ρc_weights .= ρc_cell_nd .* (wJ ./ (L_th^2))
    end
    MT = Assemble(Vi, Vj, Ni, Ni, ρc_weights, nnode)
    
    # 刚度矩阵（各向异性/各向同性自适应）：KT = ∬ Bᵀ K B dΩ
    # 当可用 jellyroll 有效张量时采用各向异性；否则退化为各向同性 k_iso.
    ngs = size(Ni,1)
    gx = mesh.gs.x[:,1]; gy = mesh.gs.x[:,2]
    use_aniso = false
    Rin = hasproperty(case.param_dim.cell, :Rin)
    Rout = hasproperty(case.param_dim.cell, :Rout) 
    if !(Rout > Rin)
        # fallback from mesh extents
        r_all = hypot.(gx, gy)
        Rin = minimum(r_all)
        Rout = maximum(r_all)
    end
    use_aniso = (Rout > Rin) && (case.opt.thermalmodel == "distributed2D")
    # 可选开关：当用户要求强制各向同性时，覆盖上面的判断
    if haskey(variables, "force_isotropic_k")
        flag = variables["force_isotropic_k"]
        if isa(flag, Float64)
            if flag > 0.5
                use_aniso = false
            end
        end
    end

    if use_aniso
        # 各向异性：优先使用元素层权重 f_k 聚合的 λ_r/λ_t（无量纲），再按 θ 旋转至 Kxx/Kxy/Kyy
        ne = size(mesh.element, 1)
        lam_r_e = zeros(Float64, ne)
        lam_t_e = zeros(Float64, ne)
        fks = haskey(variables, "thermal2D layer_weights") ? variables["thermal2D layer_weights"] : nothing
        if fks !== nothing && isa(fks, AbstractArray) && size(fks,1) == ne && size(fks,2) >= 5
            # 材料热导率（无量纲）：λ*_k = λ_k / k_ref
            λ_NE = max(getfield(case.param_dim.NE, :lambda), 0.0) / k_ref
            λ_SP = max(getfield(case.param_dim.SP, :lambda), 0.0) / k_ref
            λ_PE = max(getfield(case.param_dim.PE, :lambda), 0.0) / k_ref
            λ_PCC = max(getfield(case.param_dim.PCC, :lambda), 0.0) / k_ref
            λ_NCC = max(getfield(case.param_dim.NCC, :lambda), 0.0) / k_ref
            ϵ = 1e-12
            @inbounds for e in 1:ne
                f_NE, f_SP, f_PE, f_PCC, f_NCC = fks[e,1], fks[e,2], fks[e,3], fks[e,4], fks[e,5]
                # 径向（串联）调和平均；切向（并联）算术平均
                denom = (f_NE/max(λ_NE,ϵ)) + (f_SP/max(λ_SP,ϵ)) + (f_PE/max(λ_PE,ϵ)) + (f_PCC/max(λ_PCC,ϵ)) + (f_NCC/max(λ_NCC,ϵ))
                lam_r_e[e] = denom > 0 ? (1.0/denom) : 0.0
                lam_t_e[e] = f_NE*λ_NE + f_SP*λ_SP + f_PE*λ_PE + f_PCC*λ_PCC + f_NCC*λ_NCC
            end
        else
            # 回退：使用全局 jellyroll 等效（并确保无量纲化）
            pgeo = jellyroll_spiral_params(case.param_dim)
            λr_nd = pgeo.λ_r_eff / k_ref
            λt_nd = pgeo.λ_t_eff / k_ref
            lam_r_e .= λr_nd
            lam_t_e .= λt_nd
        end

        # 在各高斯点按极角旋转得到 Kxx/Kxy/Kyy
        Kxx = zeros(Float64, ngs)
        Kxy = zeros(Float64, ngs)
        Kyy = zeros(Float64, ngs)
        e_of_g = mesh.gs.ele
        @inbounds for g in 1:ngs
            θ = atan(gy[g], gx[g])
            c = cos(θ); s = sin(θ)
            lr = lam_r_e[e_of_g[g]]
            lt = lam_t_e[e_of_g[g]]
            # K = λr e_r e_r' + λt e_θ e_θ' → 分量
            # e_r=[c,s], e_θ=[-s,c]
            Kxx[g] = lr*c*c + lt*s*s
            Kxy[g] = (lt - lr)*s*c
            Kyy[g] = lr*s*s + lt*c*c
        end
        # 加负号以与电化学约定统一：M dT/dt = KT T（KT 包含负号）
        cxx = -Kxx .* wJ
        cxy = -Kxy .* wJ
        cyy = -Kyy .* wJ
        KT_xx = Assemble(Vi, Vj, dNdx, dNdx, cxx, nnode)
        KT_xy = Assemble(Vi, Vj, dNdx, dNdy, cxy, nnode)
        KT_yx = Assemble(Vi, Vj, dNdy, dNdx, cxy, nnode) # Kxy=Kyx
        KT_yy = Assemble(Vi, Vj, dNdy, dNdy, cyy, nnode)
        KT = KT_xx + KT_xy + KT_yx + KT_yy
    else
        # isotropic fallback with λ_iso* = λ_iso/k_ref
        # choose mean of jellyroll effective λr/λt if available; else 1.0
        let λ_iso_nd = begin
                pgeo = jellyroll_spiral_params(case.param_dim)
                if pgeo.λ_r_eff > 0 && pgeo.λ_t_eff > 0
                    ((pgeo.λ_r_eff + pgeo.λ_t_eff)/2) / k_ref
                else
                    1.0
                end
            end
            # 加负号以与电化学约定统一：M dT/dt = KT T（KT 包含负号）
            KT_x = Assemble(Vi, Vj, dNdx, dNdx, (-λ_iso_nd .* wJ), nnode)
            KT_y = Assemble(Vi, Vj, dNdy, dNdy, (-λ_iso_nd .* wJ), nnode)
            KT = KT_x + KT_y
        end
    end

    # 载荷向量 FT = ∬ q N dΩ （q 为单元平均体热源时，扩展到高斯点）
    FT = zeros(Float64, nnode)
    if haskey(variables, "heat_source_fields")
        q_elem = variables["heat_source_fields"]  # expected SI or already converted? We convert to dimensionless below if flagged
        # read unit flag as numeric code (0 = nd, 1 = SI)
        is_SI = false
        if haskey(variables, "heat_source_units_code")
            code = variables["heat_source_units_code"]
            if isa(code, Float64)
                is_SI = code > 0.5
            elseif isa(code, Array{Float64}) && length(code) > 0
                is_SI = code[1] > 0.5
            end
        end
        if is_SI
            # q_elem in W/m^3 -> convert to dimensionless q* = q / q_ref
            q_elem = q_elem ./ q_ref
        end
        # 支持 q_elem 形状 (ne,) 或 (ne,1) 或 (ne,num)
        if isa(q_elem, Array{Float64}) && size(q_elem,1) == ne
            # 取当前时间步列或第一列
            qe = ndims(q_elem) == 1 ? q_elem : q_elem[:,1]
            q_gs = qe[mesh.gs.ele]                 # (ngs)
            # dΩ* = dΩ / L_th^2
            coeff_f = q_gs .* (wJ ./ (L_th^2))    # (ngs)
            FT .+= Assemble1D(Vi, Ni, coeff_f, nnode)
        end
    end
    return MT, KT, FT
end

"""对流/辐射边界条件：外径 r=Rout (及可选内径 r=Rin) 处理，修改 K 与 F。
在此函数内进行热边界的自动识别与分类（内侧/外侧）。
实现要点：
1) 扫描所有 Q4 单元的四条边，统计出现次数，仅出现一次的为外边界边；
2) 依据边两端点半径与 (Rin, Rout) 比较分类：:inner 或 :outer；
3) 对 :outer 边组装对流项：KT += ∫ h NᵀN dL，FT += ∫ h T_amb Nᵀ dL；:inner 默认绝热（不贡献）。
备注：采用边上 2 点高斯积分；边上形函数取线性 2 节点形式 N=[(1-s)/2, (1+s)/2]，s∈[-1,1]。
"""
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64=0.0)
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh is missing in case.mesh"
    mesh = case.mesh["thermal2D"]
    @assert mesh.type == "Q4" && mesh.dimension == 2 "Boundary BC currently implemented for Q4/2D"
    # 使用 edge_boundary(:node_on) 判定外边界节点：
    ne = size(mesh.element, 1)
    nnode = mesh.nlen
    is_outer_node = falses(nnode)
    for i in 1:nnode
        # which=:outer 对应外螺旋终圈 θ ∈ [2π(N-1), 2πN]
        is_outer_node[i] = edge_boundary(:node_on, mesh, i, case.param_dim; which=:outer)
    end

    # 3) 对外边界边组装对流项（内侧默认绝热）。扫描单元四条边，若两端均为外边界节点则施加换热。
    # 使用维度化或无量纲参数需全局一致；此处沿用 param 中的 h 与 T_amb。
    # Dimensionless convection: Bi = h * L_th / k_ref already stored as scale.h_th.
    # We need effective (h_dim / k_ref) * L_th for boundary integral -> scale.h_th / L_th? For standard form with x*:
    # Original dimensional BC: -k ∂T/∂n = h (T - T_amb)
    # Non-dim: - (k/k_ref) ∂T*/∂n* = (h L_th / k_ref) (T* - T_amb*) = Bi (T* - T_amb*)
    scale = case.param_dim.scale
    Bi = scale.h_th
    L_th = max(scale.L_th, 1e-16)
    h_coeff = Bi  # operator already uses (k/k_ref) in domain, boundary uses Bi with dL*
    T_amb = case.param_dim.cell.T_amb / scale.T_ref  # dimensionless ambient
    if Bi != 0
        s_vals = (-0.577350269189626, 0.577350269189626)
        w_vals = (1.0, 1.0)
        x = mesh.node[:,1]; y = mesh.node[:,2]
        seen = Set{Tuple{Int,Int}}()  # 去重
        for e in 1:ne
            n1, n2, n3, n4 = mesh.element[e,1], mesh.element[e,2], mesh.element[e,3], mesh.element[e,4]
            for (a,b) in ((n1,n2), (n2,n3), (n3,n4), (n4,n1))
                if !(is_outer_node[a] && is_outer_node[b])
                    continue
                end
                key = a < b ? (a,b) : (b,a)
                if key in seen
                    continue
                end
                push!(seen, key)
                L = hypot(x[b]-x[a], y[b]-y[a])
                J = L/2
                ke11 = 0.0; ke12 = 0.0; ke22 = 0.0
                fe1 = 0.0; fe2 = 0.0
                for (s,w) in zip(s_vals, w_vals)
                    N1 = 0.5*(1 - s)
                    N2 = 0.5*(1 + s)
                    wt = h_coeff * w * (J / L_th)  # dL* 缩放
                    # 加负号以与体内扩散项的约定统一（KT 包含负号）
                    ke11 += -wt * (N1*N1)
                    ke12 += -wt * (N1*N2)
                    ke22 += -wt * (N2*N2)
                    fe1  += wt * (T_amb * N1)
                    fe2  += wt * (T_amb * N2)
                end
                KT[a,a] += ke11; KT[a,b] += ke12
                KT[b,a] += ke12; KT[b,b] += ke22
                FT[a]   += fe1
                FT[b]   += fe2
            end
        end
    end

    # -----------------------------------------------------------------
    # Tab (极耳) 强制温度处理：使用 jellyroll 工具定位受极耳影响的节点并施加惩罚法锚定温度
    try
        pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
        tab_nodes = unique(vcat(pos_idx, neg_idx))
        # optional debug: print counts and sample nodes
        debug_on = hasproperty(case.opt, :debug_coupling) && case.opt.debug_coupling
        if debug_on
            @info "[thermal BC] tab node counts" pos=length(pos_idx) neg=length(neg_idx) total=length(tab_nodes)
            # print up to N sample nodes with coordinates and cumulative theta used for selection
            Nsample = hasproperty(case.opt, :debug_sample_tab_nodes) ? case.opt.debug_sample_tab_nodes : 6
            ns = collect(tab_nodes)[1:min(Nsample, length(tab_nodes))]
            if !isempty(ns)
                x = mesh.node[:,1]; y = mesh.node[:,2]
                pgeo = jellyroll_spiral_params(case.param_dim)
                a, b = pgeo.a, pgeo.b
                s_in = 0.0
                s_out = pgeo.t_repeat
                for n in ns
                    xn, yn = x[n], y[n]
                    r = hypot(xn, yn)
                    # decide if this node was chosen as pos (inner) or neg (outer)
                    ispos = n in pos_idx
                    θ_cum = ispos ? ((r - a - s_in) / b) : ((r - a - s_out) / b)
                    @info "[thermal BC] tab node sample" id=n x=xn y=yn r=r theta_cum=θ_cum is_pos=ispos
                end
            end
        end

        if !isempty(tab_nodes)
            # 读取配置：加热速率（K/s）与惩罚强度
            rate_Ks = hasproperty(case.opt, :tab_heating_rate) ? case.opt.tab_heating_rate : 0.1
            penalty = hasproperty(case.opt, :tab_penalty) ? case.opt.tab_penalty : 1e12
            scale = case.param_dim.scale
            T_amb_nd = case.param_dim.cell.T_amb / scale.T_ref
            T_tab_nd = T_amb_nd + (rate_Ks * t) / max(scale.T_ref, 1e-16)
            for n in tab_nodes
                KT[n,n] += penalty
                FT[n]   += penalty * T_tab_nd
            end
        end
    catch err
        @warn "Tab BC processing failed" err
    end

    return nothing
end

"""
    heatQ_Source(case, variables, t)

计算并映射内热源 (reaction / ohmic / reversible) 到热网格节点(或元素)平均, 返回更新字典。
策略:
1. 从电化学 1D 变量中提取必要量 。
2. 对每个电化学控制体生成体积平均热源 。
3. 映射到 2D 热网格: 由于假设单元均匀。
"""
function heatQ_Source(case::Case, variables::Dict{String,Union{Array{Float64},Float64}}, t::Float64, y_state)
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh is missing"
    param = case.param
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)

    # 面积缓存
    areas = if haskey(variables, "thermal2D element area") && isa(variables["thermal2D element area"], Array{Float64})
        variables["thermal2D element area"]
    else
        A = zeros(Float64, ne)
        ngs = length(mesh_th.gs.detJ)
        @inbounds for g in 1:ngs
            e = mesh_th.gs.ele[g]
            A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
        end
        variables["thermal2D element area"] = A
        A
    end
    A_global = sum(areas)

    # 元素平均温度（使用节点场 T_nodes）
    T_n = if haskey(variables, "T_nodes")
        variables["T_nodes"]
    else
        fill(case.param.cell.T0, mesh_th.nlen)
    end
    T_e = zeros(Float64, ne)
    @inbounds for e in 1:ne
        nds = mesh_th.element[e, :]
        T_e[e] = sum(T_n[nds]) / length(nds)
    end
    I_e = variables["thermal2D element current"]
    # 正确获取层权重矩阵（而非布尔值）；此处仅取值，不更改/检查其内容
    fks = if haskey(variables, "thermal2D layer_weights")
        variables["thermal2D layer_weights"]
    else
        # 默认假设所有单元包含所有层（权重都为1）
        ones(Float64, ne, 5)
    end
    # 小工具
    to_vec(x) = isa(x, Number) ? [Float64(x)] : (isa(x, AbstractVector) ? Vector{Float64}(x) : (isa(x, AbstractArray) ? Vector{Float64}(x[:,1]) : Float64[]))
    vec_mean(x) = (isempty(x) ? 0.0 : sum(x) / length(x))
    fixval(x) = (isfinite(x) ? x : 0.0)

    q_elem = zeros(Float64, ne)
    sample_log = 0
    # 热尺度（用于长度无量纲）
    L_th = case.param_dim.scale.L_th
    for e in 1:ne
    # I_app = I_e / I1C_total
        I_nd = I_e[e]
        T_nd_e = T_e[e]
        if case.opt.model == "SPM" || case.opt.model == "SPMe"
            I_app = variables["cell current"]
            # 取界面电流密度、过电位与 dUdT
            eta_n = variables["negative electrode overpotential"][1]
            eta_p = variables["positive electrode overpotential"][end]
            csn_surf = variables["negative particle surface lithium concentration"][1]
            csp_surf = variables["positive particle surface lithium concentration"][end]
            Q_rxn = abs(I_app * (eta_p - eta_n) ) # reaction heat is always positive
            Q_rev = abs(I_app) * T_nd_e * (param.PE.dUdT(csp_surf) - param.NE.dUdT(csn_surf)) 
            # 欧姆热（电极固相/电解液/隔膜）
            sig_n_eff = param.NE.sig .* param.NE.eps_s
            sig_p_eff = param.PE.sig .* param.PE.eps_s
            kappa_ne = param.EL.kappa(param.EL.ce0, T_nd_e) * param.NE.eps ^ param.NE.brugg
            kappa_pe = param.EL.kappa(param.EL.ce0, T_nd_e) * param.PE.eps ^ param.PE.brugg
            kappa_sp = param.EL.kappa(param.EL.ce0, T_nd_e) * param.SP.eps ^ param.SP.brugg
            t_n = param.NE.thickness / L_th
            t_p = param.PE.thickness / L_th
            t_sp = param.SP.thickness / L_th
            P_s_ne = I_nd^2 * (t_n / sig_n_eff) / 3.0
            P_s_pe = I_nd^2 * (t_p / sig_p_eff) / 3.0
            P_e_ne = I_nd^2 * (t_n / kappa_ne) / 3.0
            P_e_sp = I_nd^2 * (t_sp / kappa_sp)
            P_e_pe = I_nd^2 * (t_p / kappa_pe) / 3.0
            Q_ohm = P_e_ne / t_n + P_s_ne / t_n + P_e_pe / t_p + P_s_pe / t_p + P_e_sp / t_sp 
            Q_ele = Q_rxn + Q_rev + Q_ohm
            σ_PCC = (hasproperty(case.param, :PCC) && hasproperty(case.param.PCC, :sig)) ? max(case.param.PCC.sig, 1e-12) : 1e12
            σ_NCC = (hasproperty(case.param, :NCC) && hasproperty(case.param.NCC, :sig)) ? max(case.param.NCC.sig, 1e-12) : 1e12
            Q_PCC = I_nd^2 / (3.0 * σ_PCC)
            Q_NCC = I_nd^2 / (3.0 * σ_NCC)
            # 使用已获取的 fks 矩阵进行索引（脚本已提供 fks）
            q_elem[e] = (fks[e,1] + fks[e,2] + fks[e,3]) * Q_ele + fks[e,4] * Q_PCC + fks[e,5] * Q_NCC 
        end
        if hasproperty(case.opt, :debug_coupling) && case.opt.debug_coupling && hasproperty(case.opt, :debug_sample_elems) && sample_log < case.opt.debug_sample_elems
            @info "[thermal] heat elem=$(e)" I_nd=I_nd T_nd=T_nd_e Q_PCC=Q_PCC Q_NCC=Q_NCC q_e=q_elem[e]
            sample_log += 1
        end
    end

    # 写入热源（根据单位选项）
    if hasproperty(case.opt, :units_thermal) && case.opt.units_thermal == "SI"
        variables["heat_source_fields"] = q_elem
        variables["heat_source_units_code"] = 1.0
    else
        q_ref = case.param_dim.scale.q_th
        variables["heat_source_fields"] = q_elem ./ q_ref
        variables["heat_source_units_code"] = 0.0
    end

    # 供无 y 的重载/其他模块可能使用
    try
        variables["__last_y_state"] = y_state
    catch
    end

    if hasproperty(case.opt, :debug_coupling) && case.opt.debug_coupling
        qe = variables["heat_source_fields"]
        qstats = (min=minimum(qe), max=maximum(qe), mean=mean(qe))
        units_code = variables["heat_source_units_code"]
        @info "[thermal] heat sources" units=(units_code>0.5 ? "SI W/m^3" : "nd") q_min=qstats.min q_max=qstats.max q_mean=qstats.mean
    end

    return variables
end

"""
    energy_balance_log!(case, MT, T_prev::Vector{Float64}, T_new::Vector{Float64}, dt_th::Float64, variables)

简易能量守恒日志（无量纲）：
  Q_gen* ≈ ∑_e q*_e · A*_e,   Q_conv* ≈ ∑_{outer edges} Bi · (T* - T_amb*) · L*,
  dE*/dt* ≈ (M_lumped*/dt*)·(T_new* - T_prev*).
打印 residual* = (Q_gen* - Q_conv*) - dE*/dt* 及相对残差。
仅在 opt.debug_coupling 为 true 时输出。
"""
function energy_balance_log!(case::Case, MT, T_prev::AbstractVector{<:Real}, T_new::AbstractVector{<:Real}, dt_th::Float64, variables::Dict{String,Union{Array{Float64},Float64}})
    if !hasproperty(case.opt, :debug_coupling) || !case.opt.debug_coupling
        return
    end
    haskey(case.mesh, "thermal2D") || return
    mesh = case.mesh["thermal2D"]
    scale = case.param_dim.scale
    L_th = max(scale.L_th, 1e-16)
    q_ref = scale.q_th
    Bi = scale.h_th
    T_amb_nd = case.param_dim.cell.T_amb / scale.T_ref

    # 元素面积（若未缓存则计算一次）
    areas = if haskey(variables, "thermal2D element area") && isa(variables["thermal2D element area"], Array{Float64})
        variables["thermal2D element area"]
    else
        ne = size(mesh.element,1)
        A = zeros(Float64, ne)
        ngs = length(mesh.gs.detJ)
        for g in 1:ngs
            e = mesh.gs.ele[g]
            A[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
        end
        A
    end
    # 发热 q*
    if !haskey(variables, "heat_source_fields")
        return
    end
    q_elem = variables["heat_source_fields"]
    is_SI = false
    if haskey(variables, "heat_source_units_code")
        code = variables["heat_source_units_code"]
        is_SI = (isa(code, Float64) && code > 0.5) || (isa(code, Array{Float64}) && length(code)>0 && code[1]>0.5)
    end
    q_nd = is_SI ? (q_elem ./ q_ref) : q_elem
    # ∫ q* dΩ* ≈ Σ q*_e · (A_e / L_th^2)
    A_nd = areas ./ (L_th^2)
    Q_gen_nd = sum(q_nd .* A_nd)

    # 外边界对流功率：用 edge_boundary 判定外边界节点并遍历单元边
    ne = size(mesh.element, 1)
    nnode = mesh.nlen
    x = mesh.node[:,1]; y = mesh.node[:,2]
    is_outer_node = falses(nnode)
    for i in 1:nnode
        is_outer_node[i] = edge_boundary(:node_on, mesh, i, case.param_dim; which=:outer)
    end
    seen = Set{Tuple{Int,Int}}()
    Q_conv_nd = 0.0
    for e in 1:ne
        n1, n2, n3, n4 = mesh.element[e,1], mesh.element[e,2], mesh.element[e,3], mesh.element[e,4]
        for (a,b) in ((n1,n2), (n2,n3), (n3,n4), (n4,n1))
            if !(is_outer_node[a] && is_outer_node[b])
                continue
            end
            key = a < b ? (a,b) : (b,a)
            if key in seen
                continue
            end
            push!(seen, key)
            L_nd = hypot(x[b]-x[a], y[b]-y[a]) / L_th
            Tbar = 0.5 * (T_new[a] + T_new[b])
            Q_conv_nd += Bi * (Tbar - T_amb_nd) * L_nd
        end
    end

    # 储能变化率（采用 lumped 近似）：(1/dt_th) · Σ_i M_lumped*_i · ΔT*_i
    # 使用致密化后的行和以提高兼容性
    Mrow = dropdims(sum(Matrix(MT), dims=2), dims=2)
    # 注意：MT 构造时已包含 dΩ* = dΩ/L_th^2 的缩放，这里不应再次除以 L_th^2。
    M_lumped = Mrow
    dE_dt_nd = dot(M_lumped, (T_new .- T_prev)) / max(dt_th, 1e-16)

    residual = (Q_gen_nd - Q_conv_nd) - dE_dt_nd
    rel = residual / max(abs(Q_gen_nd), 1e-12)
    @info "[thermal] Energy balance (nd)" Q_gen_nd=Q_gen_nd Q_conv_nd=Q_conv_nd dE_dt_nd=dE_dt_nd residual_nd=residual rel_residual=rel
    # 可选：写入变量以便后续诊断
    try
        variables["thermal2D energy residual"] = residual
    catch
    end
    return
end

