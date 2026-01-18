# czm.jl - 内聚力区域模型 (Cohesive Zone Model)
# 
# 功能：
# 1. 基于热网格生成分层内聚力网格
# 2. 双线性牵引力-分离本构关系
# 3. 损伤演化与加卸载准则
# 4. 与固体力学耦合的系统组装
# 5. 牛顿-拉弗森非线性求解
# 6. 损伤累积与断裂判据
#
# 内聚力单元插入位置：
# - 仅在环向的层与层之间（相邻卷绕圈的界面）
# - 不在径向的单元与单元之间
#
# 作者：JuBat Team
# 日期：2024

using LinearAlgebra
using SparseArrays
using Statistics

# ========================================================================
# 1. 数据结构定义
# ========================================================================

export CohesiveElement, CohesiveMesh, DamageState, CZMResult

"""
    CohesiveElement - 四节点零厚度内聚力单元

存储单个内聚力单元的拓扑和几何信息。

# 节点编号约定
```
      4 -------- 3   (上界面 / 外圈)
      |          |
      |  厚度=0  |   → 环向 (θ)
      |          |
      1 -------- 2   (下界面 / 内圈)
```

# 字段
- `id`: 单元编号
- `nodes`: 节点编号 [n1, n2, n3, n4]（底面2个 + 顶面2个）
- `nodes_bottom`: 底面（内圈）节点 [n1, n2]
- `nodes_top`: 顶面（外圈）节点 [n4, n3]
- `length`: 单元长度（沿环向）[m]
- `layer_idx`: 所属层间界面索引（第几圈与第几圈之间）
"""
mutable struct CohesiveElement
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
mutable struct DamageState
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
    CohesiveMesh - 内聚力网格

基于热网格扩展的分层网格，包含：
- 按圈分层的固体单元节点
- 层间界面的内聚力单元
- 节点映射关系

# 字段
- `bulk_mesh`: 原始固体网格（Q4单元）
- `node`: 扩展后的节点坐标数组 (nnode × 2)
- `nnode`: 总节点数
- `cohesive_elements`: 内聚力单元数组
- `n_cohesive`: 内聚力单元总数
- `n_layers`: 卷绕圈数
- `node_map_bottom`: 底面节点映射（原节点 → 新底面节点）
- `node_map_top`: 顶面节点映射（原节点 → 新顶面节点）
- `interface_nodes`: 每个界面的节点对列表
- `damage_states`: 每个内聚力单元的损伤状态
"""
mutable struct CohesiveMesh
    bulk_mesh::Mesh                           # 原始固体网格
    node::Matrix{Float64}                     # 扩展后节点坐标
    nnode::Int64                              # 总节点数
    bulk_element::Matrix{Int64}               # 更新后的固体单元连接关系
    cohesive_elements::Vector{CohesiveElement} # 内聚力单元
    n_cohesive::Int64                         # 内聚力单元数
    n_layers::Int64                           # 卷绕圈数
    node_map::Dict{Int64, Vector{Int64}}      # 原节点 → [分层后的节点们]
    interface_nodes::Vector{Vector{Tuple{Int64,Int64}}} # 每个界面的节点对
    damage_states::Vector{DamageState}        # 损伤状态
    
    # 内部构造函数（空初始化）
    function CohesiveMesh()
        new(Mesh("Q4", 2, zeros(0,2), 0, zeros(Int64,0,4), 
            GaussPoint(zeros(0,2), zeros(0,2), zeros(0), zeros(0), zeros(Int64,0), zeros(0,4), zeros(0,8), 2)),
            zeros(0, 2), 0, zeros(Int64, 0, 4),
            CohesiveElement[], 0, 0, Dict{Int64, Vector{Int64}}(),
            Vector{Vector{Tuple{Int64,Int64}}}(), DamageState[])
    end
end

"""
    CZMResult - 内聚力分析结果

存储求解后的结果。

# 字段
- `displacement`: 位移场 (ndof,)
- `damage`: 每个内聚力单元的损伤值
- `traction_n`: 法向牵引力
- `traction_t`: 切向牵引力
- `separation_n`: 法向分离位移
- `separation_t`: 切向分离位移
- `converged`: 是否收敛
- `iterations`: 迭代次数
- `residual_norm`: 最终残差范数
"""
mutable struct CZMResult
    displacement::Vector{Float64}
    damage::Vector{Float64}
    traction_n::Vector{Float64}
    traction_t::Vector{Float64}
    separation_n::Vector{Float64}
    separation_t::Vector{Float64}
    converged::Bool
    iterations::Int64
    residual_norm::Float64
    
    CZMResult(ndof::Int, n_coh::Int) = new(
        zeros(ndof), zeros(n_coh), zeros(n_coh), zeros(n_coh),
        zeros(n_coh), zeros(n_coh), false, 0, Inf)
end


# ========================================================================
# 2. 网格划分
# ========================================================================

export create_czm_mesh, identify_layer_interfaces

"""
    create_czm_mesh(thermal_mesh, param_dim; nθ_per_turn=nothing)

基于热网格创建内聚力网格。

# 算法
1. 识别卷绕圈数和每圈的单元
2. 在相邻圈的界面处复制节点
3. 创建四节点零厚度内聚力单元
4. 更新固体单元的节点连接关系

# 参数
- `thermal_mesh`: 热分析网格（Q4单元）
- `param_dim`: 参数对象（包含螺旋几何信息）
- `nθ_per_turn`: 每圈的单元数（可选，自动检测）

# 返回
- `CohesiveMesh`: 内聚力网格对象
"""
function create_czm_mesh(thermal_mesh::Mesh, param_dim; nθ_per_turn::Union{Int,Nothing}=nothing)
    @assert thermal_mesh.type == "Q4" "create_czm_mesh requires Q4 mesh"
    @assert thermal_mesh.dimension == 2 "create_czm_mesh requires 2D mesh"
    
    # 获取螺旋参数
    p = jellyroll_spiral_params(param_dim)
    
    # 识别层数和每层单元
    ne = size(thermal_mesh.element, 1)
    nnode_orig = thermal_mesh.nlen
    
    # 估算每圈单元数（如果未指定）
    if nθ_per_turn === nothing
        # 从螺旋参数估算圈数
        n_turns = max(1, round(Int, (p.Rout - p.Rin) / p.t_repeat))
        nθ_per_turn = max(1, div(ne, n_turns))
    end
    
    # 计算层数（圈数）
    n_layers = max(1, div(ne, nθ_per_turn))
    n_interfaces = n_layers - 1  # 层间界面数
    
    if n_interfaces < 1
        @warn "Only one layer detected, no cohesive elements will be created"
        # 返回一个空的CohesiveMesh但保留原始网格信息
        czm_mesh = CohesiveMesh()
        czm_mesh.bulk_mesh = thermal_mesh
        czm_mesh.node = copy(thermal_mesh.node)
        czm_mesh.nnode = nnode_orig
        czm_mesh.bulk_element = copy(thermal_mesh.element)
        czm_mesh.n_layers = 1
        return czm_mesh
    end
    
    # 识别每层的单元和界面节点
    layer_elements, interface_info = _identify_layers_and_interfaces(
        thermal_mesh, nθ_per_turn, n_layers)
    
    # 复制界面节点并创建映射
    node_new, node_map, bulk_element_new = _duplicate_interface_nodes(
        thermal_mesh, layer_elements, interface_info, n_interfaces)
    
    # 创建内聚力单元
    cohesive_elements = _create_cohesive_elements(
        thermal_mesh, node_new, interface_info, node_map, n_interfaces)
    
    # 初始化损伤状态
    n_cohesive = length(cohesive_elements)
    damage_states = [DamageState() for _ in 1:n_cohesive]
    
    # 构建CohesiveMesh
    czm_mesh = CohesiveMesh()
    czm_mesh.bulk_mesh = thermal_mesh
    czm_mesh.node = node_new
    czm_mesh.nnode = size(node_new, 1)
    czm_mesh.bulk_element = bulk_element_new
    czm_mesh.cohesive_elements = cohesive_elements
    czm_mesh.n_cohesive = n_cohesive
    czm_mesh.n_layers = n_layers
    czm_mesh.node_map = node_map
    czm_mesh.interface_nodes = [info.node_pairs for info in interface_info]
    czm_mesh.damage_states = damage_states
    
    @info "Created CZM mesh" n_layers=n_layers n_interfaces=n_interfaces n_cohesive=n_cohesive nnode_new=czm_mesh.nnode
    
    return czm_mesh
end

"""
层间界面信息结构
"""
struct InterfaceInfo
    layer_below::Int64           # 下层索引
    layer_above::Int64           # 上层索引
    elements_below::Vector{Int64} # 下层相邻单元
    elements_above::Vector{Int64} # 上层相邻单元
    node_pairs::Vector{Tuple{Int64,Int64}}  # (底面节点, 顶面节点) 对
end

"""
识别层和界面
"""
function _identify_layers_and_interfaces(mesh::Mesh, nθ_per_turn::Int, n_layers::Int)
    ne = size(mesh.element, 1)
    nnode_orig = mesh.nlen
    
    # 按角度位置分配单元到层
    # 假设单元按θ顺序编号，每nθ_per_turn个单元为一层
    layer_elements = Vector{Vector{Int64}}(undef, n_layers)
    for layer in 1:n_layers
        start_e = (layer - 1) * nθ_per_turn + 1
        end_e = min(layer * nθ_per_turn, ne)
        layer_elements[layer] = collect(start_e:end_e)
    end
    
    # 识别层间界面
    # 对于条带网格，界面在每层的外边界与下一层的内边界之间
    n_interfaces = n_layers - 1
    interface_info = Vector{InterfaceInfo}(undef, n_interfaces)
    
    for iface in 1:n_interfaces
        layer_below = iface
        layer_above = iface + 1
        
        # 获取相邻层的单元
        elem_below = layer_elements[layer_below]
        elem_above = layer_elements[layer_above]
        
        # 识别界面节点对
        # 在条带网格中，层i的外边界节点与层i+1的内边界节点位置相同
        node_pairs = _find_interface_node_pairs(mesh, elem_below, elem_above)
        
        interface_info[iface] = InterfaceInfo(
            layer_below, layer_above,
            elem_below, elem_above,
            node_pairs
        )
    end
    
    return layer_elements, interface_info
end

"""
查找界面节点对
"""
function _find_interface_node_pairs(mesh::Mesh, elem_below::Vector{Int64}, elem_above::Vector{Int64})
    # 获取下层的外边界节点（节点2和3，即外螺旋上的节点）
    outer_nodes_below = Set{Int64}()
    for e in elem_below
        push!(outer_nodes_below, mesh.element[e, 2])  # 外侧当前
        push!(outer_nodes_below, mesh.element[e, 3])  # 外侧下一个
    end
    
    # 获取上层的内边界节点（节点1和4，即内螺旋上的节点）
    inner_nodes_above = Set{Int64}()
    for e in elem_above
        push!(inner_nodes_above, mesh.element[e, 1])  # 内侧当前
        push!(inner_nodes_above, mesh.element[e, 4])  # 内侧下一个
    end
    
    # 找到坐标匹配的节点对
    node_pairs = Tuple{Int64, Int64}[]
    tol = 1e-10  # 坐标匹配容差
    
    for n_below in outer_nodes_below
        x_below = mesh.node[n_below, 1]
        y_below = mesh.node[n_below, 2]
        
        for n_above in inner_nodes_above
            x_above = mesh.node[n_above, 1]
            y_above = mesh.node[n_above, 2]
            
            if abs(x_below - x_above) < tol && abs(y_below - y_above) < tol
                push!(node_pairs, (n_below, n_above))
                break
            end
        end
    end
    
    # 按角度排序节点对
    if !isempty(node_pairs)
        sort!(node_pairs, by = p -> atan(mesh.node[p[1], 2], mesh.node[p[1], 1]))
    end
    
    return node_pairs
end

"""
复制界面节点并创建映射
"""
function _duplicate_interface_nodes(mesh::Mesh, layer_elements::Vector{Vector{Int64}}, 
                                    interface_info::Vector{InterfaceInfo}, n_interfaces::Int)
    nnode_orig = mesh.nlen
    ne = size(mesh.element, 1)
    
    # 初始化节点映射：原节点 → [新节点列表]
    node_map = Dict{Int64, Vector{Int64}}()
    for i in 1:nnode_orig
        node_map[i] = [i]  # 初始时每个节点映射到自己
    end
    
    # 复制节点坐标数组
    node_list = [mesh.node[i, :] for i in 1:nnode_orig]
    
    # 复制单元连接关系
    bulk_element = copy(mesh.element)
    
    # 处理每个界面
    for iface in 1:n_interfaces
        info = interface_info[iface]
        
        # 对于每个界面节点对，创建新节点（用于上层）
        for (n_below, n_above) in info.node_pairs
            # 如果 n_below 和 n_above 是同一个节点（共享），需要分离
            if n_below == n_above
                # 创建新节点（坐标相同）
                new_node_idx = length(node_list) + 1
                push!(node_list, mesh.node[n_above, :])
                
                # 更新映射
                push!(node_map[n_above], new_node_idx)
                
                # 更新上层单元的节点连接
                for e in info.elements_above
                    for k in 1:4
                        if bulk_element[e, k] == n_above
                            bulk_element[e, k] = new_node_idx
                        end
                    end
                end
            end
        end
    end
    
    # 转换为矩阵
    nnode_new = length(node_list)
    node_new = zeros(Float64, nnode_new, 2)
    for i in 1:nnode_new
        node_new[i, :] = node_list[i]
    end
    
    return node_new, node_map, bulk_element
end

"""
创建内聚力单元
"""
function _create_cohesive_elements(mesh::Mesh, node_new::Matrix{Float64},
                                   interface_info::Vector{InterfaceInfo},
                                   node_map::Dict{Int64, Vector{Int64}},
                                   n_interfaces::Int)
    cohesive_elements = CohesiveElement[]
    elem_id = 0
    
    for iface in 1:n_interfaces
        info = interface_info[iface]
        n_pairs = length(info.node_pairs)
        
        # 每两个相邻节点对形成一个内聚力单元
        for i in 1:(n_pairs - 1)
            n1_orig, n1_top_orig = info.node_pairs[i]
            n2_orig, n2_top_orig = info.node_pairs[i + 1]
            
            # 获取映射后的节点编号
            # 底面节点：使用原节点（属于下层）
            n1_bottom = n1_orig
            n2_bottom = n2_orig
            
            # 顶面节点：使用新创建的节点（属于上层）
            # 查找映射中的最后一个节点（新创建的）
            n1_top = node_map[n1_top_orig][end]
            n2_top = node_map[n2_top_orig][end]
            
            # 如果没有创建新节点（节点对不是共享的情况），使用原节点
            if n1_top == n1_top_orig && length(node_map[n1_top_orig]) == 1
                n1_top = n1_top_orig
            end
            if n2_top == n2_top_orig && length(node_map[n2_top_orig]) == 1
                n2_top = n2_top_orig
            end
            
            # 计算单元长度
            x1 = node_new[n1_bottom, 1]
            y1 = node_new[n1_bottom, 2]
            x2 = node_new[n2_bottom, 1]
            y2 = node_new[n2_bottom, 2]
            elem_length = hypot(x2 - x1, y2 - y1)
            
            elem_id += 1
            
            # 创建内聚力单元
            # 节点顺序：[底1, 底2, 顶2, 顶1] 形成逆时针
            coh_elem = CohesiveElement(
                elem_id,
                [n1_bottom, n2_bottom, n2_top, n1_top],
                [n1_bottom, n2_bottom],
                [n1_top, n2_top],
                elem_length,
                iface
            )
            
            push!(cohesive_elements, coh_elem)
        end
    end
    
    return cohesive_elements
end


# ========================================================================
# 3. 双线性本构模型
# ========================================================================

export bilinear_traction, bilinear_tangent, update_damage!

"""
    bilinear_traction(δ_n, δ_t, damage_state, cohesive_params; update=false)

计算双线性本构的牵引力。

# 参数
- `δ_n`: 法向分离位移（正为张开）
- `δ_t`: 切向分离位移
- `damage_state`: 损伤状态
- `cohesive_params`: 内聚力参数
- `update`: 是否更新损伤状态

# 返回
- `T_n`: 法向牵引力
- `T_t`: 切向牵引力
- `D`: 当前损伤变量
"""
function bilinear_traction(δ_n::Float64, δ_t::Float64, damage_state::DamageState, 
                          cohesive_params::Cohesive; update::Bool=false)
    # 提取参数
    K_n = cohesive_params.K_n
    K_t = cohesive_params.K_t
    δ_0_n = cohesive_params.δ_0_n
    δ_c_n = cohesive_params.δ_c_n
    δ_0_t = cohesive_params.δ_0_t
    δ_c_t = cohesive_params.δ_c_t
    η = cohesive_params.eta
    
    # 已断裂检查
    if damage_state.fractured
        return 0.0, 0.0, 1.0
    end
    
    # 计算等效分离位移（混合模式）
    # 使用Macaulay括号处理压缩：只有张开才贡献
    δ_n_pos = max(0.0, δ_n)  # 张开部分
    δ_eff = sqrt(δ_n_pos^2 + δ_t^2)
    
    # 混合模式临界值（BK准则）
    if δ_eff > 1e-15
        β = abs(δ_t) / δ_eff  # 模式混合比
        δ_0_eff = sqrt(δ_0_n^2 + (δ_0_t^2 - δ_0_n^2) * β^η)
        δ_c_eff = sqrt(δ_c_n^2 + (δ_c_t^2 - δ_c_n^2) * β^η)
    else
        δ_0_eff = δ_0_n
        δ_c_eff = δ_c_n
    end
    
    # 获取历史最大分离（加卸载判断）
    δ_max_hist = damage_state.δ_max_eff
    
    # 计算当前损伤变量
    D = damage_state.D
    
    if δ_eff > δ_max_hist  # 加载
        if δ_eff <= δ_0_eff
            # 弹性阶段
            D = 0.0
        elseif δ_eff >= δ_c_eff
            # 完全断裂
            D = 1.0
        else
            # 损伤软化阶段
            D = δ_c_eff * (δ_eff - δ_0_eff) / (δ_eff * (δ_c_eff - δ_0_eff))
        end
        
        # 更新历史
        if update
            damage_state.δ_max_eff = δ_eff
            damage_state.δ_max_n = max(damage_state.δ_max_n, δ_n_pos)
            damage_state.δ_max_t = max(damage_state.δ_max_t, abs(δ_t))
            damage_state.D = D
            damage_state.accumulated_damage = max(damage_state.accumulated_damage, D)
            
            if D >= 1.0 - 1e-10
                damage_state.fractured = true
            end
        end
    else
        # 卸载：D保持不变，使用历史值
        D = damage_state.D
    end
    
    # 计算牵引力
    # 法向：考虑压缩时的惩罚接触
    if δ_n >= 0
        T_n = (1.0 - D) * K_n * δ_n
    else
        # 压缩时使用纯弹性接触（不受损伤影响）
        T_n = K_n * δ_n
    end
    
    # 切向
    T_t = (1.0 - D) * K_t * δ_t
    
    return T_n, T_t, D
end

"""
    bilinear_tangent(δ_n, δ_t, damage_state, cohesive_params)

计算双线性本构的切线刚度矩阵。

# 返回
- `dT_dδ`: 2×2 切线刚度矩阵 [dTn/dδn, dTn/dδt; dTt/dδn, dTt/dδt]
"""
function bilinear_tangent(δ_n::Float64, δ_t::Float64, damage_state::DamageState,
                         cohesive_params::Cohesive)
    # 提取参数
    K_n = cohesive_params.K_n
    K_t = cohesive_params.K_t
    δ_0_n = cohesive_params.δ_0_n
    δ_c_n = cohesive_params.δ_c_n
    δ_0_t = cohesive_params.δ_0_t
    δ_c_t = cohesive_params.δ_c_t
    η = cohesive_params.eta
    
    # 初始化切线矩阵
    dT_dδ = zeros(Float64, 2, 2)
    
    # 已断裂
    if damage_state.fractured
        # 返回小刚度避免奇异
        dT_dδ[1, 1] = 1e-10 * K_n
        dT_dδ[2, 2] = 1e-10 * K_t
        return dT_dδ
    end
    
    # 计算等效分离
    δ_n_pos = max(0.0, δ_n)
    δ_eff = sqrt(δ_n_pos^2 + δ_t^2)
    
    # 混合模式参数
    if δ_eff > 1e-15
        β = abs(δ_t) / δ_eff
        δ_0_eff = sqrt(δ_0_n^2 + (δ_0_t^2 - δ_0_n^2) * β^η)
        δ_c_eff = sqrt(δ_c_n^2 + (δ_c_t^2 - δ_c_n^2) * β^η)
    else
        δ_0_eff = δ_0_n
        δ_c_eff = δ_c_n
    end
    
    D = damage_state.D
    δ_max_hist = damage_state.δ_max_eff
    
    # 判断加载/卸载状态
    is_loading = (δ_eff > δ_max_hist - 1e-15)
    
    if δ_eff <= δ_0_eff || !is_loading
        # 弹性阶段或卸载
        # dT_n/dδ_n
        if δ_n >= 0
            dT_dδ[1, 1] = (1.0 - D) * K_n
        else
            dT_dδ[1, 1] = K_n  # 压缩接触
        end
        # dT_t/dδ_t
        dT_dδ[2, 2] = (1.0 - D) * K_t
        
    elseif δ_eff >= δ_c_eff
        # 完全断裂
        dT_dδ[1, 1] = 1e-10 * K_n
        dT_dδ[2, 2] = 1e-10 * K_t
        
    else
        # 损伤软化阶段 - 需要计算完整的切线
        # D = δ_c * (δ_eff - δ_0) / (δ_eff * (δ_c - δ_0))
        # dD/dδ_eff = δ_c * δ_0 / (δ_eff^2 * (δ_c - δ_0))
        
        dD_dδeff = δ_c_eff * δ_0_eff / (δ_eff^2 * (δ_c_eff - δ_0_eff))
        
        # dδ_eff/dδ_n = δ_n_pos / δ_eff (当 δ_n > 0)
        # dδ_eff/dδ_t = δ_t / δ_eff
        
        if δ_n >= 0 && δ_eff > 1e-15
            dδeff_dδn = δ_n_pos / δ_eff
            dδeff_dδt = δ_t / δ_eff
            
            # T_n = (1-D) * K_n * δ_n
            # dT_n/dδ_n = (1-D)*K_n - K_n*δ_n*dD/dδ_eff*dδeff/dδ_n
            dT_dδ[1, 1] = (1.0 - D) * K_n - K_n * δ_n * dD_dδeff * dδeff_dδn
            dT_dδ[1, 2] = -K_n * δ_n * dD_dδeff * dδeff_dδt
            
            # T_t = (1-D) * K_t * δ_t
            dT_dδ[2, 1] = -K_t * δ_t * dD_dδeff * dδeff_dδn
            dT_dδ[2, 2] = (1.0 - D) * K_t - K_t * δ_t * dD_dδeff * dδeff_dδt
        else
            # 压缩或零分离
            dT_dδ[1, 1] = K_n
            dT_dδ[2, 2] = (1.0 - D) * K_t
        end
    end
    
    return dT_dδ
end

"""
    update_damage!(damage_states, separations, cohesive_params)

批量更新所有内聚力单元的损伤状态。

# 参数
- `damage_states`: 损伤状态数组
- `separations`: 分离位移数组 [(δ_n, δ_t), ...]
- `cohesive_params`: 内聚力参数
"""
function update_damage!(damage_states::Vector{DamageState}, 
                       separations::Vector{Tuple{Float64, Float64}},
                       cohesive_params::Cohesive)
    n = length(damage_states)
    @assert length(separations) == n "Mismatch in array lengths"
    
    for i in 1:n
        δ_n, δ_t = separations[i]
        bilinear_traction(δ_n, δ_t, damage_states[i], cohesive_params; update=true)
    end
end


# ========================================================================
# 4. 单元计算
# ========================================================================

export cohesive_element_matrices, compute_separation

"""
    compute_separation(elem, node, u)

计算内聚力单元的分离位移。

# 参数
- `elem`: 内聚力单元
- `node`: 节点坐标数组
- `u`: 位移向量 (ndof,)

# 返回
- `δ_n`: 法向分离位移（在高斯点）
- `δ_t`: 切向分离位移（在高斯点）
- `n_vec`: 法向向量
- `t_vec`: 切向向量
"""
function compute_separation(elem::CohesiveElement, node::Matrix{Float64}, u::Vector{Float64})
    # 获取节点坐标
    n1, n2 = elem.nodes_bottom
    n4, n3 = elem.nodes_top  # 注意：nodes_top = [n4, n3]，但顺序对应 [n1, n2]
    
    x1, y1 = node[n1, 1], node[n1, 2]
    x2, y2 = node[n2, 1], node[n2, 2]
    
    # 计算单元的局部坐标系
    # 切向：沿着单元长度方向
    dx = x2 - x1
    dy = y2 - y1
    L = hypot(dx, dy)
    
    if L < 1e-15
        return 0.0, 0.0, [0.0, 1.0], [1.0, 0.0]
    end
    
    t_vec = [dx / L, dy / L]  # 切向单位向量
    n_vec = [-t_vec[2], t_vec[1]]  # 法向单位向量（逆时针90度）
    
    # 获取节点位移
    # DOF编号：节点i -> [2i-1, 2i] = [u_x, u_y]
    u_bottom = zeros(Float64, 2, 2)  # [节点, 方向]
    u_top = zeros(Float64, 2, 2)
    
    u_bottom[1, 1] = u[2*n1 - 1]; u_bottom[1, 2] = u[2*n1]
    u_bottom[2, 1] = u[2*n2 - 1]; u_bottom[2, 2] = u[2*n2]
    u_top[1, 1] = u[2*n4 - 1]; u_top[1, 2] = u[2*n4]
    u_top[2, 1] = u[2*n3 - 1]; u_top[2, 2] = u[2*n3]
    
    # 在单元中心（ξ=0）计算分离
    # 形函数：N1 = 0.5, N2 = 0.5
    u_bottom_mid = 0.5 * (u_bottom[1, :] + u_bottom[2, :])
    u_top_mid = 0.5 * (u_top[1, :] + u_top[2, :])
    
    # 分离位移 = 顶面位移 - 底面位移
    Δu = u_top_mid - u_bottom_mid
    
    # 投影到局部坐标系
    δ_n = dot(Δu, n_vec)  # 法向分离（正为张开）
    δ_t = dot(Δu, t_vec)  # 切向分离
    
    return δ_n, δ_t, n_vec, t_vec
end

"""
    cohesive_element_matrices(elem, node, u, damage_state, cohesive_params)

计算单个内聚力单元的刚度矩阵和内力向量。

# 参数
- `elem`: 内聚力单元
- `node`: 节点坐标数组
- `u`: 当前位移向量
- `damage_state`: 损伤状态
- `cohesive_params`: 内聚力参数

# 返回
- `K_e`: 单元刚度矩阵 (8×8)
- `f_int_e`: 单元内力向量 (8,)
- `δ_n, δ_t`: 分离位移
- `T_n, T_t`: 牵引力
"""
function cohesive_element_matrices(elem::CohesiveElement, node::Matrix{Float64},
                                  u::Vector{Float64}, damage_state::DamageState,
                                  cohesive_params::Cohesive)
    # 单元有4个节点，每节点2个DOF，共8个DOF
    ndof_e = 8
    K_e = zeros(Float64, ndof_e, ndof_e)
    f_int_e = zeros(Float64, ndof_e)
    
    # 获取节点编号
    n1, n2 = elem.nodes_bottom
    n4, n3 = elem.nodes_top
    
    # 计算单元几何
    x1, y1 = node[n1, 1], node[n1, 2]
    x2, y2 = node[n2, 1], node[n2, 2]
    L = elem.length
    
    if L < 1e-15
        return K_e, f_int_e, 0.0, 0.0, 0.0, 0.0
    end
    
    # 局部坐标系
    dx = x2 - x1
    dy = y2 - y1
    t_vec = [dx / L, dy / L]
    n_vec = [-t_vec[2], t_vec[1]]
    
    # 旋转矩阵：将全局坐标系下的位移转换到局部坐标系
    # R = [n_x, n_y; t_x, t_y]
    R = [n_vec[1] n_vec[2]; t_vec[1] t_vec[2]]
    
    # 2点高斯积分
    gauss_pts = [-1.0/sqrt(3.0), 1.0/sqrt(3.0)]
    gauss_wts = [1.0, 1.0]
    
    # 累积牵引力用于输出
    T_n_avg = 0.0
    T_t_avg = 0.0
    δ_n_avg = 0.0
    δ_t_avg = 0.0
    
    for (ξ, w) in zip(gauss_pts, gauss_wts)
        # 形函数
        N1 = 0.5 * (1.0 - ξ)
        N2 = 0.5 * (1.0 + ξ)
        
        # B矩阵：将节点位移映射到高斯点分离
        # Δu = u_top - u_bottom
        # 底面：节点1,2；顶面：节点4,3
        # DOF顺序：[u1x, u1y, u2x, u2y, u3x, u3y, u4x, u4y]
        #          = [底1x, 底1y, 底2x, 底2y, 顶2x, 顶2y, 顶1x, 顶1y]
        
        # 分离位移：Δu_x = N1*u4x + N2*u3x - N1*u1x - N2*u2x
        #          Δu_y = N1*u4y + N2*u3y - N1*u1y - N2*u2y
        
        # B_global: [Δu_x; Δu_y] = B_global * u_e
        B_global = zeros(Float64, 2, 8)
        # 底面贡献（负号）
        B_global[1, 1] = -N1; B_global[2, 2] = -N1  # 节点1
        B_global[1, 3] = -N2; B_global[2, 4] = -N2  # 节点2
        # 顶面贡献（正号）
        B_global[1, 5] = N2; B_global[2, 6] = N2    # 节点3 (对应底面节点2)
        B_global[1, 7] = N1; B_global[2, 8] = N1    # 节点4 (对应底面节点1)
        
        # 转换到局部坐标系
        # [δ_n; δ_t] = R * [Δu_x; Δu_y] = R * B_global * u_e
        B_local = R * B_global
        
        # 计算当前高斯点的分离位移
        # 获取单元位移向量
        u_e = zeros(Float64, 8)
        u_e[1] = u[2*n1 - 1]; u_e[2] = u[2*n1]      # 底1
        u_e[3] = u[2*n2 - 1]; u_e[4] = u[2*n2]      # 底2
        u_e[5] = u[2*n3 - 1]; u_e[6] = u[2*n3]      # 顶2
        u_e[7] = u[2*n4 - 1]; u_e[8] = u[2*n4]      # 顶1
        
        δ_local = B_local * u_e
        δ_n = δ_local[1]
        δ_t = δ_local[2]
        
        # 计算牵引力和切线刚度
        T_n, T_t, D = bilinear_traction(δ_n, δ_t, damage_state, cohesive_params; update=false)
        dT_dδ = bilinear_tangent(δ_n, δ_t, damage_state, cohesive_params)
        
        # 雅可比行列式（1D）
        J = L / 2.0
        
        # 组装单元刚度矩阵
        # K_e = ∫ B^T * D_tan * B * dA
        # 对于2D内聚力单元，dA = 1 * J * dξ（单位厚度）
        K_e += w * J * (B_local' * dT_dδ * B_local)
        
        # 组装单元内力向量
        # f_int = ∫ B^T * T * dA
        T_local = [T_n, T_t]
        f_int_e += w * J * (B_local' * T_local)
        
        # 累积平均值
        T_n_avg += 0.5 * T_n
        T_t_avg += 0.5 * T_t
        δ_n_avg += 0.5 * δ_n
        δ_t_avg += 0.5 * δ_t
    end
    
    return K_e, f_int_e, δ_n_avg, δ_t_avg, T_n_avg, T_t_avg
end


# ========================================================================
# 5. 系统组装
# ========================================================================

export assemble_czm_system, assemble_coupled_system

"""
    assemble_czm_system(czm_mesh, u, cohesive_params)

组装内聚力单元的全局刚度矩阵和内力向量。

# 返回
- `K_coh`: 内聚力刚度矩阵 (ndof × ndof)
- `f_int_coh`: 内聚力内力向量 (ndof,)
- `separations`: 每个单元的分离位移
- `tractions`: 每个单元的牵引力
"""
function assemble_czm_system(czm_mesh::CohesiveMesh, u::Vector{Float64}, 
                            cohesive_params::Cohesive)
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive
    
    # 稀疏矩阵组装
    I_idx = Int64[]
    J_idx = Int64[]
    K_vals = Float64[]
    
    f_int_coh = zeros(Float64, ndof)
    
    # 存储每个单元的分离和牵引
    separations = Vector{Tuple{Float64, Float64}}(undef, n_coh)
    tractions = Vector{Tuple{Float64, Float64}}(undef, n_coh)
    
    for (i, elem) in enumerate(czm_mesh.cohesive_elements)
        damage_state = czm_mesh.damage_states[i]
        
        K_e, f_int_e, δ_n, δ_t, T_n, T_t = cohesive_element_matrices(
            elem, czm_mesh.node, u, damage_state, cohesive_params)
        
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
    
    # 高斯积分点
    gauss_pts = [-1.0/sqrt(3.0), 1.0/sqrt(3.0)]
    gauss_wts = [1.0, 1.0]
    
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
        
        for (ξ, wξ) in zip(gauss_pts, gauss_wts)
            for (η, wη) in zip(gauss_pts, gauss_wts)
                # Q4形函数导数
                dNdxi = 0.25 * [-(1-η), (1-η), (1+η), -(1+η)]
                dNdeta = 0.25 * [-(1-ξ), -(1+ξ), (1+ξ), (1-ξ)]
                
                # 雅可比矩阵
                J11 = sum(dNdxi .* x_e)
                J12 = sum(dNdxi .* y_e)
                J21 = sum(dNdeta .* x_e)
                J22 = sum(dNdeta .* y_e)
                
                detJ = J11 * J22 - J12 * J21
                
                if abs(detJ) < 1e-15
                    continue
                end
                
                # 物理导数
                invJ = [J22 -J12; -J21 J11] / detJ
                dNdx = invJ[1,1] * dNdxi + invJ[1,2] * dNdeta
                dNdy = invJ[2,1] * dNdxi + invJ[2,2] * dNdeta
                
                # B矩阵
                B = zeros(Float64, 3, 8)
                for i in 1:4
                    B[1, 2*i-1] = dNdx[i]
                    B[2, 2*i] = dNdy[i]
                    B[3, 2*i-1] = dNdy[i]
                    B[3, 2*i] = dNdx[i]
                end
                
                # 积分
                K_e += wξ * wη * detJ * (B' * D_mat * B)
            end
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
    assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params)

组装耦合的固体-内聚力系统。

# 返回
- `K_total`: 总刚度矩阵
- `f_int_total`: 总内力向量
- `separations`: 分离位移
- `tractions`: 牵引力
"""
function assemble_coupled_system(czm_mesh::CohesiveMesh, u::Vector{Float64},
                                E_eff::Float64, ν_eff::Float64, 
                                cohesive_params::Cohesive)
    # 固体刚度
    K_bulk = assemble_bulk_stiffness(czm_mesh, E_eff, ν_eff)
    
    # 内聚力刚度
    K_coh, f_int_coh, separations, tractions = assemble_czm_system(
        czm_mesh, u, cohesive_params)
    
    # 固体内力（线性弹性）
    f_int_bulk = K_bulk * u
    
    # 总系统
    K_total = K_bulk + K_coh
    f_int_total = f_int_bulk + f_int_coh
    
    return K_total, f_int_total, separations, tractions
end


# ========================================================================
# 6. 边界条件
# ========================================================================

export apply_bc_czm!

"""
    apply_bc_czm!(K, f, bc_dofs, bc_vals)

应用位移边界条件（惩罚法）。

# 参数
- `K`: 刚度矩阵（会被修改）
- `f`: 载荷向量（会被修改）
- `bc_dofs`: 约束DOF列表
- `bc_vals`: 约束值
"""
function apply_bc_czm!(K::SparseMatrixCSC{Float64,Int64}, f::Vector{Float64},
                      bc_dofs::Vector{Int64}, bc_vals::Vector{Float64})
    penalty = 1e20
    
    for (dof, val) in zip(bc_dofs, bc_vals)
        K[dof, dof] += penalty
        f[dof] = penalty * val
    end
end

"""
    identify_bc_nodes_czm(czm_mesh, param_dim)

识别需要施加边界条件的节点。

# 返回
- `bc_nodes`: Dict{Int, Symbol}，节点 → 边界类型
"""
function identify_bc_nodes_czm(czm_mesh::CohesiveMesh, param_dim)
    bc_nodes = Dict{Int64, Symbol}()
    
    # 获取螺旋参数
    p = jellyroll_spiral_params(param_dim)
    
    # 识别内外边界节点
    tol = 1e-4
    for i in 1:czm_mesh.nnode
        x, y = czm_mesh.node[i, 1], czm_mesh.node[i, 2]
        r = hypot(x, y)
        
        # 内边界：固定
        if abs(r - p.Rin) < tol
            bc_nodes[i] = :fixed_xy
        end
        
        # 外边界：可以自由或施加载荷
        # 默认不约束
    end
    
    return bc_nodes
end


# ========================================================================
# 7. 牛顿-拉弗森求解器
# ========================================================================

export newton_raphson_czm, solve_czm_step

"""
    newton_raphson_czm(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param_dim;
                       max_iter=50, tol=1e-8, u0=nothing)

牛顿-拉弗森非线性求解。

# 参数
- `czm_mesh`: 内聚力网格
- `F_ext`: 外力向量
- `E_eff`: 有效杨氏模量
- `ν_eff`: 有效泊松比
- `cohesive_params`: 内聚力参数
- `param_dim`: 参数对象
- `max_iter`: 最大迭代次数
- `tol`: 收敛容差
- `u0`: 初始位移（可选）

# 返回
- `result`: CZMResult 结果对象
"""
function newton_raphson_czm(czm_mesh::CohesiveMesh, F_ext::Vector{Float64},
                           E_eff::Float64, ν_eff::Float64,
                           cohesive_params::Cohesive, param_dim;
                           max_iter::Int=50, tol::Float64=1e-8,
                           u0::Union{Vector{Float64},Nothing}=nothing)
    
    nnode = czm_mesh.nnode
    ndof = 2 * nnode
    n_coh = czm_mesh.n_cohesive
    
    # 初始化结果
    result = CZMResult(ndof, n_coh)
    
    # 初始位移
    u = u0 === nothing ? zeros(Float64, ndof) : copy(u0)
    
    # 识别边界条件
    bc_nodes = identify_bc_nodes_czm(czm_mesh, param_dim)
    bc_dofs = Int64[]
    bc_vals = Float64[]
    
    for (node, bc_type) in bc_nodes
        if bc_type == :fixed_xy
            push!(bc_dofs, 2*node - 1)
            push!(bc_vals, 0.0)
            push!(bc_dofs, 2*node)
            push!(bc_vals, 0.0)
        elseif bc_type == :fixed_x
            push!(bc_dofs, 2*node - 1)
            push!(bc_vals, 0.0)
        elseif bc_type == :fixed_y
            push!(bc_dofs, 2*node)
            push!(bc_vals, 0.0)
        end
    end
    
    # 牛顿-拉弗森迭代
    for iter in 1:max_iter
        # 组装系统
        K_total, f_int_total, separations, tractions = assemble_coupled_system(
            czm_mesh, u, E_eff, ν_eff, cohesive_params)
        
        # 残差
        R = F_ext - f_int_total
        
        # 应用边界条件到残差
        for (dof, val) in zip(bc_dofs, bc_vals)
            R[dof] = val - u[dof]
        end
        
        # 检查收敛
        R_norm = norm(R)
        
        if iter == 1
            R_norm_0 = max(R_norm, 1e-10)
        end
        
        rel_norm = R_norm / R_norm_0
        
        if R_norm < tol || rel_norm < tol
            result.converged = true
            result.iterations = iter
            result.residual_norm = R_norm
            result.displacement = u
            
            # 更新损伤状态
            update_damage!(czm_mesh.damage_states, separations, cohesive_params)
            
            # 存储结果
            for i in 1:n_coh
                result.damage[i] = czm_mesh.damage_states[i].D
                result.separation_n[i] = separations[i][1]
                result.separation_t[i] = separations[i][2]
                result.traction_n[i] = tractions[i][1]
                result.traction_t[i] = tractions[i][2]
            end
            
            @info "Newton-Raphson converged" iterations=iter residual=R_norm
            return result
        end
        
        # 应用边界条件到刚度矩阵
        K_bc = copy(K_total)
        R_bc = copy(R)
        apply_bc_czm!(K_bc, R_bc, bc_dofs, zeros(length(bc_dofs)))
        
        # 求解增量
        Δu = K_bc \ R_bc
        
        # 更新位移
        u += Δu
        
        # 强制边界条件
        for (dof, val) in zip(bc_dofs, bc_vals)
            u[dof] = val
        end
    end
    
    # 未收敛
    @warn "Newton-Raphson did not converge" max_iter=max_iter residual=norm(R)
    result.converged = false
    result.iterations = max_iter
    result.displacement = u
    
    return result
end

"""
    solve_czm_step(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param_dim, u_prev;
                   max_iter=50, tol=1e-8)

求解单个时间步的内聚力问题。
保持损伤历史的连续性。

# 参数
- `czm_mesh`: 内聚力网格（损伤状态会被更新）
- `F_ext`: 外力向量
- `u_prev`: 上一步的位移

# 返回
- `result`: CZMResult
"""
function solve_czm_step(czm_mesh::CohesiveMesh, F_ext::Vector{Float64},
                       E_eff::Float64, ν_eff::Float64,
                       cohesive_params::Cohesive, param_dim,
                       u_prev::Vector{Float64};
                       max_iter::Int=50, tol::Float64=1e-8)
    
    # 使用上一步位移作为初始猜测
    result = newton_raphson_czm(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param_dim;
                                max_iter=max_iter, tol=tol, u0=u_prev)
    
    return result
end


# ========================================================================
# 8. 损伤累积与断裂判据
# ========================================================================

export get_damage_statistics, check_fracture_criterion, reset_damage_states!

"""
    get_damage_statistics(czm_mesh)

获取损伤统计信息。

# 返回
- NamedTuple: 包含 max_D, mean_D, n_fractured, fraction_damaged 等
"""
function get_damage_statistics(czm_mesh::CohesiveMesh)
    n = czm_mesh.n_cohesive
    
    if n == 0
        return (max_D=0.0, mean_D=0.0, min_D=0.0, 
                n_fractured=0, fraction_damaged=0.0,
                total_accumulated=0.0)
    end
    
    D_vals = [s.D for s in czm_mesh.damage_states]
    n_fractured = count(s -> s.fractured, czm_mesh.damage_states)
    n_damaged = count(d -> d > 0.01, D_vals)  # 超过1%认为有损伤
    
    accumulated = [s.accumulated_damage for s in czm_mesh.damage_states]
    
    return (
        max_D = maximum(D_vals),
        mean_D = mean(D_vals),
        min_D = minimum(D_vals),
        n_fractured = n_fractured,
        fraction_damaged = n_damaged / n,
        total_accumulated = sum(accumulated)
    )
end

"""
    check_fracture_criterion(czm_mesh; threshold=0.99)

检查是否满足宏观断裂判据。

# 判据
1. 连续断裂路径：相邻内聚力单元全部断裂
2. 整体损伤阈值：平均损伤超过阈值

# 返回
- `(is_fractured::Bool, fracture_info::NamedTuple)`
"""
function check_fracture_criterion(czm_mesh::CohesiveMesh; threshold::Float64=0.99)
    stats = get_damage_statistics(czm_mesh)
    
    # 简单判据：平均损伤超过阈值
    is_fractured_avg = stats.mean_D >= threshold
    
    # 或者：断裂单元比例超过50%
    is_fractured_count = (stats.n_fractured / max(1, czm_mesh.n_cohesive)) > 0.5
    
    is_fractured = is_fractured_avg || is_fractured_count
    
    fracture_info = (
        is_fractured = is_fractured,
        criterion = is_fractured_avg ? :average_damage : (is_fractured_count ? :fractured_count : :none),
        stats = stats
    )
    
    return is_fractured, fracture_info
end

"""
    reset_damage_states!(czm_mesh)

重置所有损伤状态（用于新分析）。
"""
function reset_damage_states!(czm_mesh::CohesiveMesh)
    for state in czm_mesh.damage_states
        state.D = 0.0
        state.δ_max_n = 0.0
        state.δ_max_t = 0.0
        state.δ_max_eff = 0.0
        state.fractured = false
        state.accumulated_damage = 0.0
    end
end

"""
    accumulate_cycle_damage!(czm_mesh, cycle_damage_increment)

累积循环损伤（用于疲劳分析）。

# 参数
- `czm_mesh`: 内聚力网格
- `cycle_damage_increment`: 每个循环的损伤增量
"""
function accumulate_cycle_damage!(czm_mesh::CohesiveMesh, cycle_damage_increment::Float64)
    for state in czm_mesh.damage_states
        if !state.fractured
            state.accumulated_damage += cycle_damage_increment
            
            # 如果累积损伤超过1，标记为断裂
            if state.accumulated_damage >= 1.0
                state.D = 1.0
                state.fractured = true
            end
        end
    end
end


# ========================================================================
# 9. 后处理与输出
# ========================================================================

export czm_output_to_variables

"""
    czm_output_to_variables(czm_mesh, result, variables)

将CZM结果写入variables字典。
"""
function czm_output_to_variables(czm_mesh::CohesiveMesh, result::CZMResult,
                                variables::Dict{String, Union{Array{Float64}, Float64}})
    # 位移
    nnode = czm_mesh.nnode
    u_x = result.displacement[1:2:end]
    u_y = result.displacement[2:2:end]
    variables["czm displacement x"] = u_x
    variables["czm displacement y"] = u_y
    
    # 损伤场
    variables["czm damage"] = result.damage
    
    # 牵引力
    variables["czm traction normal"] = result.traction_n
    variables["czm traction tangent"] = result.traction_t
    
    # 分离位移
    variables["czm separation normal"] = result.separation_n
    variables["czm separation tangent"] = result.separation_t
    
    # 统计信息
    stats = get_damage_statistics(czm_mesh)
    variables["czm max damage"] = stats.max_D
    variables["czm mean damage"] = stats.mean_D
    variables["czm fractured elements"] = Float64(stats.n_fractured)
    
    return variables
end
