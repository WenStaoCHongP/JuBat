"""
极耳边界条件模块（新版）

提供三种边界处理方式：
1. 惩罚法（penalty）：强制温度，数值不稳定（不推荐）
2. 一般表面散热（surface_convection）：模拟整体z方向冷却
3. 极耳强化散热（tab_convection）：模拟极耳区域局部强化冷却

使用方法：
在 opt 中设置 opt.tab_bc_type 选择边界类型
"""

# ========================================================================
# 方式1：惩罚法（传统方法，不推荐）
# ========================================================================

"""
    _apply_tab_bc_penalty!(KT, FT, mesh, case, t)

惩罚法强制极耳节点温度（传统方法，数值不稳定）

# 原理
通过在刚度矩阵对角线添加极大的惩罚值 penalty，强制节点温度等于指定值：
    K[i,i] += penalty
    F[i] += penalty * T_prescribed

# 缺点
- penalty 过大（默认1e12）导致矩阵病态
- 与小时间步长结合时，M - K*dt 可能出现负对角元素
- 引发数值不稳定，产生 NaN

# 参数
- `penalty`: 惩罚系数，建议 1e6-1e8（如果必须使用惩罚法）
- `tab_heating_rate`: 极耳升温速率 [K/s]，默认 0.1
"""
function _apply_tab_bc_penalty!(KT, FT, mesh, case, t)
    try
        pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
        tab_nodes = unique(vcat(pos_idx, neg_idx))
        
        isempty(tab_nodes) && return
        
        # 参数
        rate_Ks = hasproperty(case.opt, :tab_heating_rate) ? case.opt.tab_heating_rate : 0.1
        penalty = hasproperty(case.opt, :tab_penalty) ? case.opt.tab_penalty : 1e6  # 降低默认值
        
        scale = case.param_dim.scale
        T_amb_nd = case.param_dim.cell.T_amb / scale.T_ref
        T_tab_nd = T_amb_nd + (rate_Ks * t) / scale.T_ref
        
        # 惩罚法
        for n in tab_nodes
            KT[n, n] += penalty
            FT[n] += penalty * T_tab_nd
        end
        
        # 警告
        if penalty > 1e8
            @warn "惩罚值过大，可能导致数值不稳定" penalty=penalty
        end
        
    catch err
        @warn "惩罚法边界条件失败" exception=(err, catch_backtrace())
    end
end

# ========================================================================
# 方式2：一般表面散热（推荐）
# ========================================================================

"""
    _apply_surface_convection_bc!(KT, FT, mesh, case, t)

一般表面散热：模拟电池上下表面（z方向）整体对流散热

# 物理模型
电池通过上下表面与环境换热：
    q_z = h * (T - T_amb)

在2D模型中，将z方向的表面积投影到节点：
    A_z(i) = A_voronoi(i) * H

其中 H 为电池高度（z方向）。

# 弱形式贡献
刚度矩阵：K_ij += ∫ h * (H / L_th²) * N_i * N_j dΩ
载荷向量：F_i += ∫ h * T_amb * (H / L_th²) * N_i dΩ

# 无量纲化
Biot数：Bi_z = h * H / k_th
对流项：K_ij^* = ∫ Bi_z * (H / L_th) * N_i * N_j dΩ^*

# 参数
- `h_surface`: 表面对流换热系数 [W/(m²·K)]，典型值：5-50
- `cell_height`: 电池高度 [m]，默认使用 case.param_dim.cell.width
"""
function _apply_surface_convection_bc!(KT, FT, mesh, case, t)
    try
        # 获取参数
        h_surface = hasproperty(case.opt, :h_surface) ? case.opt.h_surface : 10.0  # W/(m²·K)
        H = hasproperty(case.param_dim.cell, :height) ? case.param_dim.cell.height : case.param_dim.cell.width
        
        scale = case.param_dim.scale
        k_th = scale.k_th
        L_th = scale.L_th
        T_ref = scale.T_ref
        T_amb_nd = case.param_dim.cell.T_amb / T_ref
        
        # Biot数（无量纲）
        Bi_z = h_surface * H / k_th
        
        # 尺度因子：将z方向面密度转换为2D积分
        # factor = Bi_z * (H / L_th) / L_th²
        # 简化：conv_factor = Bi_z * H / L_th³
        conv_factor = Bi_z * H / L_th^3
        
        if conv_factor < 1e-12
            return  # 对流可忽略
        end
        
        # 高斯积分数据
        ngs = length(mesh.gs.detJ)
        Ni = mesh.gs.Ni  # (ngs, nnodes_per_elem)
        wJ = mesh.gs.weight .* mesh.gs.detJ
        ele = mesh.gs.ele
        
        nn_per_elem = size(mesh.element, 2)  # Q4 = 4
        nnode = mesh.nlen
        
        # 装配对流项到所有节点
        for g in 1:ngs
            e = ele[g]
            nodes = mesh.element[e, :]
            wt = conv_factor * wJ[g]
            
            # 单元对流矩阵和载荷
            for i in 1:nn_per_elem
                ni = nodes[i]
                Ni_g = Ni[g, i]
                
                for j in 1:nn_per_elem
                    nj = nodes[j]
                    Nj_g = Ni[g, j]
                    
                    # K_ij += wt * N_i * N_j
                    KT[ni, nj] += wt * Ni_g * Nj_g
                end
                
                # F_i += wt * T_amb * N_i
                FT[ni] += wt * T_amb_nd * Ni_g
            end
        end
        
        # 调试信息
        if hasproperty(case.opt, :debug_coupling) && case.opt.debug_coupling
            @info "[surface_convection] 应用整体表面散热" h=h_surface Bi_z=Bi_z H=H
        end
        
    catch err
        @warn "表面对流边界条件失败" exception=(err, catch_backtrace())
    end
end

# ========================================================================
# 方式3：极耳强化散热（推荐）
# ========================================================================

"""
    _apply_tab_convection_bc!(KT, FT, mesh, case, t)

极耳强化散热：仅在极耳节点施加增强的z方向对流散热

# 物理模型
极耳区域通过与外壳接触或其他冷却措施，散热能力强于其他区域：
    q_z,tab = h_tab * (T - T_amb)

h_tab >> h_surface（极耳换热系数远大于一般表面）

# 实现方式
1. 识别极耳节点（通过 jellyroll_tab_node_indices）
2. 计算每个极耳节点的"影响面积"（Voronoi面积）
3. 施加增强对流边界条件

# 节点面积估算
使用单元面积平均分配：
    A_node(i) = Σ_{e ∋ i} A_element(e) / n_nodes_per_elem

# 参数
- `h_tab`: 极耳对流换热系数 [W/(m²·K)]，典型值：50-500
- `cell_height`: 电池高度 [m]
"""
function _apply_tab_convection_bc!(KT, FT, mesh, case, t)
    try
        # 识别极耳节点
        pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
        tab_nodes = unique(vcat(pos_idx, neg_idx))
        
        isempty(tab_nodes) && return
        
        # 获取参数
        h_tab = hasproperty(case.opt, :h_tab) ? case.opt.h_tab : 100.0  # W/(m²·K)
        H = hasproperty(case.param_dim.cell, :height) ? case.param_dim.cell.height : case.param_dim.cell.width
        
        scale = case.param_dim.scale
        k_th = scale.k_th
        L_th = scale.L_th
        T_ref = scale.T_ref
        T_amb_nd = case.param_dim.cell.T_amb / T_ref
        
        # Biot数（无量纲）
        Bi_z_tab = h_tab * H / k_th
        
        # 计算节点影响面积（通过单元面积分配）
        node_areas = _compute_node_areas(mesh)
        
        # 对每个极耳节点施加对流边界条件
        for n in tab_nodes
            A_node = node_areas[n]  # 节点在x-y平面的面积
            A_z = A_node * H  # z方向投影面积 = A_xy * H
            
            # 无量纲化
            A_z_nd = A_z / L_th^2
            
            # K[n,n] += h_tab * A_z / (k_th * L_th)
            #         = (h_tab * H / k_th) * (A_xy / L_th^2)
            #         = Bi_z_tab * A_z_nd
            KT[n, n] += Bi_z_tab * A_z_nd
            
            # F[n] += h_tab * T_amb * A_z / (k_th * L_th * T_ref)
            #       = Bi_z_tab * T_amb_nd * A_z_nd
            FT[n] += Bi_z_tab * T_amb_nd * A_z_nd
        end
        
        # 调试信息
        if hasproperty(case.opt, :debug_coupling) && case.opt.debug_coupling
            @info "[tab_convection] 应用极耳强化散热" h_tab=h_tab Bi_z=Bi_z_tab n_nodes=length(tab_nodes)
        end
        
    catch err
        @warn "极耳对流边界条件失败" exception=(err, catch_backtrace())
    end
end

"""
    _compute_node_areas(mesh) -> Vector{Float64}

计算每个节点的影响面积（2D平面）

# 方法
将每个单元的面积平均分配给其节点：
    A_node(i) = Σ_{e ∋ i} A_elem(e) / n_nodes_per_elem

对于Q4单元，使用高斯积分计算精确面积：
    A_elem = Σ_g weight[g] * detJ[g]
"""
function _compute_node_areas(mesh)
    nnode = mesh.nlen
    ne = size(mesh.element, 1)
    nn_per_elem = size(mesh.element, 2)
    
    node_areas = zeros(Float64, nnode)
    
    # 计算单元面积
    elem_areas = zeros(Float64, ne)
    ngs = length(mesh.gs.detJ)
    
    for g in 1:ngs
        e = mesh.gs.ele[g]
        elem_areas[e] += mesh.gs.weight[g] * mesh.gs.detJ[g]
    end
    
    # 分配给节点
    for e in 1:ne
        A_e = elem_areas[e]
        A_per_node = A_e / nn_per_elem
        
        for i in 1:nn_per_elem
            n = mesh.element[e, i]
            node_areas[n] += A_per_node
        end
    end
    
    return node_areas
end

# ========================================================================
# 统一接口：根据 opt.tab_bc_type 选择边界条件类型
# ========================================================================

"""
    _apply_tab_bc!(KT, FT, mesh, case, t)

极耳边界条件统一接口（新版）

# 使用方法
在 case.opt 中设置 tab_bc_type 选择边界类型：

```julia
opt.tab_bc_type = "surface_convection"  # 方式1：整体表面散热
opt.h_surface = 10.0

opt.tab_bc_type = "tab_convection"      # 方式2：极耳强化散热
opt.h_tab = 100.0

opt.tab_bc_type = "penalty"             # 方式3：惩罚法（不推荐）
opt.tab_penalty = 1e6
```

# 默认行为
如果未设置 tab_bc_type，则：
- 如果设置了 h_tab，使用 tab_convection
- 否则，使用 penalty（兼容旧代码）
"""
function _apply_tab_bc!(KT, FT, mesh, case, t)
    # 确定边界类型
    bc_type = if hasproperty(case.opt, :tab_bc_type)
        case.opt.tab_bc_type
    else
        # 向后兼容：如果设置了 h_tab，默认使用对流法
        hasproperty(case.opt, :h_tab) ? "tab_convection" : "penalty"
    end
    
    # 调用对应的函数
    if bc_type == "surface_convection"
        _apply_surface_convection_bc!(KT, FT, mesh, case, t)
    elseif bc_type == "tab_convection"
        _apply_tab_convection_bc!(KT, FT, mesh, case, t)
    elseif bc_type == "penalty"
        _apply_tab_bc_penalty!(KT, FT, mesh, case, t)
    else
        @warn "未知的极耳边界类型，使用默认惩罚法" bc_type=bc_type
        _apply_tab_bc_penalty!(KT, FT, mesh, case, t)
    end
end
