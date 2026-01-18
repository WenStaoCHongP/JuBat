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

export create_czm_mesh, identify_interface_node_pairs

"""
    create_czm_mesh(thermal_mesh, param_dim; tol=1e-8)

基于热网格创建内聚力网格。

# 核心算法
通过检查节点坐标是否重合来识别层间界面：
- 外螺旋在θ位置的点 与 内螺旋在θ+2π位置的点 坐标重合
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
    nnode_orig = thermal_mesh.nlen
    
    # 识别内外螺旋节点
    # 在 jellyroll_collector_seed_mesh 中：
    # - 节点 1 到 nθ+1 在内螺旋上
    # - 节点 nθ+2 到 2*(nθ+1) 在外螺旋上
    # nθ = ne (单元数)
    nθ = ne
    inner_nodes = collect(1:(nθ+1))
    outer_nodes = collect((nθ+2):(2*(nθ+1)))
    
    # 通过坐标重合检测界面节点对
    interface_pairs = _find_coincident_node_pairs(thermal_mesh, inner_nodes, outer_nodes, tol)
    
    if isempty(interface_pairs)
        @warn "No interface nodes found (no coincident nodes). Check mesh structure."
        czm_mesh = CohesiveMesh()
        czm_mesh.bulk_mesh = thermal_mesh
        czm_mesh.node = copy(thermal_mesh.node)
        czm_mesh.nnode = nnode_orig
        czm_mesh.bulk_element = copy(thermal_mesh.element)
        czm_mesh.n_layers = 1
        return czm_mesh
    end
    
    # 按角度排序界面节点对
    sort!(interface_pairs, by = p -> atan(thermal_mesh.node[p[1], 2], thermal_mesh.node[p[1], 1]))
    
    # 将界面节点对分组（按连续性分成多个界面）
    interface_groups = _group_interface_pairs(thermal_mesh, interface_pairs)
    n_interfaces = length(interface_groups)
    
    @info "Detected interfaces" n_interfaces=n_interfaces total_pairs=length(interface_pairs)
    
    # 复制界面节点并更新单元连接
    node_new, node_map, bulk_element_new, interface_node_info = _split_interface_nodes(
        thermal_mesh, interface_groups)
    
    # 创建内聚力单元
    cohesive_elements = _create_cohesive_elements_from_pairs(
        node_new, interface_node_info)
    
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
    czm_mesh.n_layers = n_interfaces + 1
    czm_mesh.node_map = node_map
    czm_mesh.interface_nodes = [grp.node_pairs for grp in interface_node_info]
    czm_mesh.damage_states = damage_states
    
    @info "Created CZM mesh" n_layers=czm_mesh.n_layers n_interfaces=n_interfaces n_cohesive=n_cohesive nnode_new=czm_mesh.nnode
    
    return czm_mesh
end

"""
    _find_coincident_node_pairs(mesh, inner_nodes, outer_nodes, tol)

通过坐标重合检测界面节点对。

外螺旋在θ位置的点与内螺旋在θ+2π位置的点坐标重合，
这些重合点就是相邻卷绕圈之间的界面。

# 返回
- `Vector{Tuple{Int64, Int64}}`: (外螺旋节点, 内螺旋节点) 对列表
"""
function _find_coincident_node_pairs(mesh::Mesh, inner_nodes::Vector{Int64}, 
                                     outer_nodes::Vector{Int64}, tol::Float64)
    pairs = Tuple{Int64, Int64}[]
    
    # 对于每个外螺旋节点，检查是否有内螺旋节点与其坐标重合
    for n_out in outer_nodes
        x_out = mesh.node[n_out, 1]
        y_out = mesh.node[n_out, 2]
        
        for n_in in inner_nodes
            x_in = mesh.node[n_in, 1]
            y_in = mesh.node[n_in, 2]
            
            # 检查坐标是否重合
            if abs(x_out - x_in) < tol && abs(y_out - y_in) < tol
                # (外螺旋节点, 内螺旋节点)
                # 外螺旋节点属于"下层"单元的外边界
                # 内螺旋节点属于"上层"单元的内边界
                push!(pairs, (n_out, n_in))
                break  # 每个外螺旋节点最多匹配一个内螺旋节点
            end
        end
    end
    
    return pairs
end

"""
界面节点组信息
"""
struct InterfaceNodeGroup
    interface_idx::Int64                      # 界面索引
    node_pairs::Vector{Tuple{Int64, Int64}}   # (外螺旋节点, 内螺旋节点) 对
    new_node_map::Dict{Int64, Int64}          # 原内螺旋节点 → 新节点编号
end

"""
    _group_interface_pairs(mesh, pairs)

将界面节点对按连续性分组（识别不同的界面）。

由于螺旋结构，相邻卷绕圈的界面在θ方向上是连续的。
通过检查节点的θ位置来判断是否属于同一界面。
"""
function _group_interface_pairs(mesh::Mesh, pairs::Vector{Tuple{Int64, Int64}})
    if isempty(pairs)
        return InterfaceNodeGroup[]
    end
    
    # 计算每个节点对的θ位置
    θ_vals = [atan(mesh.node[p[1], 2], mesh.node[p[1], 1]) for p in pairs]
    
    # 将θ归一化到 [0, 2π]
    θ_vals = [θ < 0 ? θ + 2π : θ for θ in θ_vals]
    
    # 按θ排序
    sorted_idx = sortperm(θ_vals)
    sorted_pairs = pairs[sorted_idx]
    sorted_θ = θ_vals[sorted_idx]
    
    # 分组：检查θ的跳变来识别不同的界面
    # 如果θ跳变超过π，认为是新的界面
    groups = Vector{Vector{Tuple{Int64, Int64}}}()
    current_group = [sorted_pairs[1]]
    
    for i in 2:length(sorted_pairs)
        Δθ = sorted_θ[i] - sorted_θ[i-1]
        
        # 如果θ跳变较大，但这只是因为回到了2π的边界，不应该分组
        # 真正的界面分隔应该是节点在不同的"圈"上
        # 简单起见：如果相邻节点的θ差距很大（>π），可能是边界回绕，不分组
        # 更准确的方法是检查半径
        
        # 这里我们使用一个简单的方法：所有重合的节点对都属于同一类界面
        # 因为它们都是"外螺旋-内螺旋"的重合点
        push!(current_group, sorted_pairs[i])
    end
    
    push!(groups, current_group)
    
    # 创建InterfaceNodeGroup对象
    result = InterfaceNodeGroup[]
    for (idx, grp) in enumerate(groups)
        push!(result, InterfaceNodeGroup(idx, grp, Dict{Int64, Int64}()))
    end
    
    return result
end

"""
    _split_interface_nodes(mesh, interface_groups)

在界面处复制节点，将原本重合的节点分离。

对于每个界面节点对 (n_out, n_in)：
- n_out（外螺旋节点）保持不变，属于"下层"单元
- n_in（内螺旋节点）被复制为新节点，新节点属于"上层"单元
- 原节点位置保持不变，但连接关系更新

# 返回
- `node_new`: 扩展后的节点坐标数组
- `node_map`: 节点映射关系
- `bulk_element_new`: 更新后的单元连接
- `interface_node_info`: 更新后的界面信息（包含新节点编号）
"""
function _split_interface_nodes(mesh::Mesh, interface_groups::Vector{InterfaceNodeGroup})
    nnode_orig = mesh.nlen
    ne = size(mesh.element, 1)
    
    # 收集所有需要复制的内螺旋节点
    nodes_to_duplicate = Set{Int64}()
    for grp in interface_groups
        for (n_out, n_in) in grp.node_pairs
            push!(nodes_to_duplicate, n_in)
        end
    end
    
    # 创建新节点数组
    n_new_nodes = length(nodes_to_duplicate)
    node_new = zeros(Float64, nnode_orig + n_new_nodes, 2)
    node_new[1:nnode_orig, :] = mesh.node
    
    # 创建节点映射
    node_map = Dict{Int64, Vector{Int64}}()
    for i in 1:nnode_orig
        node_map[i] = [i]
    end
    
    # 复制节点并记录映射
    # old_to_new: 原内螺旋节点 → 新节点编号
    old_to_new = Dict{Int64, Int64}()
    new_node_idx = nnode_orig
    
    for n_in in nodes_to_duplicate
        new_node_idx += 1
        # 新节点坐标与原节点相同
        node_new[new_node_idx, :] = mesh.node[n_in, :]
        old_to_new[n_in] = new_node_idx
        push!(node_map[n_in], new_node_idx)
    end
    
    # 更新单元连接关系
    # 识别哪些单元使用了被复制的内螺旋节点，并且这些单元属于"上层"
    # "上层"单元的内边界节点应该使用新节点
    
    bulk_element_new = copy(mesh.element)
    nθ = ne
    
    # 在条带网格中，单元的内边界节点是节点1和4
    # 如果这些节点是被复制的节点，检查该单元是否是"上层"单元
    # "上层"单元的定义：其内边界节点参与了界面（与某个外螺旋节点重合）
    
    # 创建集合：哪些节点是界面上的内螺旋节点
    interface_inner_nodes = Set{Int64}()
    for grp in interface_groups
        for (n_out, n_in) in grp.node_pairs
            push!(interface_inner_nodes, n_in)
        end
    end
    
    # 更新单元的内边界节点
    for e in 1:ne
        for local_idx in [1, 4]  # 内边界节点的局部索引
            n = bulk_element_new[e, local_idx]
            if n in interface_inner_nodes
                # 这个单元的内边界节点在界面上，使用新节点
                bulk_element_new[e, local_idx] = old_to_new[n]
            end
        end
    end
    
    # 更新interface_groups中的映射
    interface_node_info = InterfaceNodeGroup[]
    for grp in interface_groups
        new_map = Dict{Int64, Int64}()
        for (n_out, n_in) in grp.node_pairs
            new_map[n_in] = old_to_new[n_in]
        end
        push!(interface_node_info, InterfaceNodeGroup(grp.interface_idx, grp.node_pairs, new_map))
    end
    
    return node_new, node_map, bulk_element_new, interface_node_info
end

"""
    _create_cohesive_elements_from_pairs(node, interface_info)

基于界面节点对创建内聚力单元。

每两个相邻的界面节点对形成一个四节点内聚力单元：
```
    n_in_new[i+1] -------- n_in_new[i]    (上层，使用新节点)
          |                    |
          |      厚度=0        |
          |                    |
    n_out[i+1] -------- n_out[i]          (下层，使用原节点)
```
"""
function _create_cohesive_elements_from_pairs(node::Matrix{Float64}, 
                                              interface_info::Vector{InterfaceNodeGroup})
    cohesive_elements = CohesiveElement[]
    elem_id = 0
    
    for grp in interface_info
        pairs = grp.node_pairs
        new_map = grp.new_node_map
        n_pairs = length(pairs)
        
        if n_pairs < 2
            continue
        end
        
        # 按θ排序节点对
        sorted_pairs = sort(pairs, by = p -> atan(node[p[1], 2], node[p[1], 1]))
        
        # 每两个相邻节点对形成一个内聚力单元
        for i in 1:(n_pairs - 1)
            n_out_1, n_in_1 = sorted_pairs[i]
            n_out_2, n_in_2 = sorted_pairs[i + 1]
            
            # 获取新节点编号（上层使用新节点）
            n_in_new_1 = new_map[n_in_1]
            n_in_new_2 = new_map[n_in_2]
            
            # 计算单元长度
            x1, y1 = node[n_out_1, 1], node[n_out_1, 2]
            x2, y2 = node[n_out_2, 1], node[n_out_2, 2]
            elem_length = hypot(x2 - x1, y2 - y1)
            
            elem_id += 1
            
            # 创建内聚力单元
            # 节点顺序：[底1, 底2, 顶2, 顶1] 
            # 底面（下层）：使用原外螺旋节点
            # 顶面（上层）：使用新内螺旋节点
            coh_elem = CohesiveElement(
                elem_id,
                [n_out_1, n_out_2, n_in_new_2, n_in_new_1],
                [n_out_1, n_out_2],           # 底面节点
                [n_in_new_1, n_in_new_2],     # 顶面节点
                elem_length,
                grp.interface_idx
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
