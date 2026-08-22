# CZM cohesive mesh topology and construction.

mutable struct CohesiveElement <: AbstractCohesiveElement
    id::Int64
    nodes::Vector{Int64}           # [n1, n2, n3, n4]
    nodes_bottom::Vector{Int64}    # [n1, n2] 底面节点
    nodes_top::Vector{Int64}       # [n4, n3] 顶面节点（顺序与底面一致）
    length::Float64                # 单元长度
    interface_type::Symbol         # 2 种本构/材料类型之一；不表示 4 个真实面的计数
    host_outer_elem::Int           # 外层 Q4 单元 id（在 czm_submesh.mesh.element 中的行号）
    host_inner_elem::Int           # 内层 Q4 单元 id
end

"""
    create_czm_mesh(czm_submesh::CzmSubmesh, thermal_mesh::Mesh, param) -> CohesiveMesh

基于细化 CZM 子网格构造内聚力网格（spec §4.4）。

# 核心算法
1. 建立 共边(2 节点) → 单元对 映射，识别每个周向分段的 4 个真实箔–涂层面；按材料组合归入 PE-PCC / NE-NCC 两种类型
2. 节点复制：对每个界面对的共边 2 节点生成副本（memoized），cohesive 单元 4 节点 = [n_a, n_b, n_b', n_a']
3. **重写外层 bulk 单元连接**：把外层单元共边位置的原节点替换为副本节点（关键：否则分离位移恒为 0）
4. 构造 cohesive_to_thermal[e_coh] = thermal_elem_map[e_outer]

# 参数
- `czm_submesh`: 细化 CZM 子网格（含 material_type / thermal_elem_map）
- `thermal_mesh`: 粗热网格（用于校验 thermal_elem_map 索引范围；本函数不复制其连接）
- `param`: 参数对象（保留签名一致性，当前未使用）

# 返回
- `CohesiveMesh`: 内聚力网格对象，bulk_mesh 指向 czm_submesh.mesh
"""
function create_czm_mesh(czm_submesh::CzmSubmesh, thermal_mesh::Mesh, param)
    # v1.5 §3.4：拓扑消费 Φ 合并网格（完美粘结）；插值按 .mesh 未合并布局构造后按 phi_keep 裁剪
    sub_mesh = czm_submesh.mesh_bonded
    ne_sub = size(sub_mesh.element, 1)
    nnode_sub = sub_mesh.nlen
    n_thermal = size(thermal_mesh.element, 1)

    # Step 1: 建立 共边 → 单元对 映射
    edge_to_elems = Dict{Tuple{Int, Int}, Vector{Int}}()
    for e in 1:ne_sub
        n1, n2, n3, n4 = sub_mesh.element[e, :]
        for edge in ((n1, n2), (n2, n3), (n3, n4), (n4, n1))
            key = (min(edge[1], edge[2]), max(edge[1], edge[2]))
            push!(get!(edge_to_elems, key, Int[]), e)
        end
    end

    # Step 2: 遍历共边，识别 PE-PCC / NE-NCC 径向界面
    interface_pairs = Tuple{Int, Int, Symbol}[]   # (e_inner, e_outer, interface_type)
    for (edge, elems) in edge_to_elems
        length(elems) == 2 || continue   # 周向相邻同层（材料相同）自动过滤
        e1, e2 = elems[1], elems[2]
        m1, m2 = czm_submesh.material_type[e1], czm_submesh.material_type[e2]
        iface = if (m1 == :PE && m2 == :PCC) || (m1 == :PCC && m2 == :PE)
            :PE_PCC
        elseif (m1 == :NE && m2 == :NCC) || (m1 == :NCC && m2 == :NE)
            :NE_NCC
        else
            nothing
        end
        if iface !== nothing
            # 判断哪个是内层（径向更小）——用单元质心判断，鲁棒于节点排列变化
            cx1 = sum(sub_mesh.node[sub_mesh.element[e1, c], 1] for c in 1:4) / 4
            cy1 = sum(sub_mesh.node[sub_mesh.element[e1, c], 2] for c in 1:4) / 4
            cx2 = sum(sub_mesh.node[sub_mesh.element[e2, c], 1] for c in 1:4) / 4
            cy2 = sum(sub_mesh.node[sub_mesh.element[e2, c], 2] for c in 1:4) / 4
            r1 = hypot(cx1, cy1)
            r2 = hypot(cx2, cy2)
            if r1 < r2
                push!(interface_pairs, (e1, e2, iface))
            else
                push!(interface_pairs, (e2, e1, iface))
            end
        end
    end

    # Step 3: 节点复制 + 重写外层 bulk 连接
    # 排序保证 cohesive 单元 id 跨运行一致（不依赖 Dict 哈希顺序）
    sort!(interface_pairs, by = first)
    node_copy = Dict{Int, Int}()

    n_cohesive = length(interface_pairs)
    max_new_nodes = 2 * n_cohesive
    extended_node = zeros(Float64, nnode_sub + max_new_nodes, 2)
    extended_node[1:nnode_sub, :] = sub_mesh.node
    new_node_count = nnode_sub

    bulk_element_new = Matrix{Int}(sub_mesh.element)

    cohesive_elements = CohesiveElement[]
    cohesive_to_thermal = Vector{Int}(undef, n_cohesive)
    sizehint!(cohesive_elements, n_cohesive)

    for (i, (e_inner, e_outer, iface)) in enumerate(interface_pairs)
        inner_nodes = sub_mesh.element[e_inner, :]
        outer_nodes = sub_mesh.element[e_outer, :]
        common_set = intersect(Set(inner_nodes), Set(outer_nodes))
        @assert length(common_set) == 2 "共边应有 2 节点，实际 $(length(common_set))"
        common = collect(common_set)

        # build_czm_submesh 按 layer-outer / segment-inner 顺序赋 node id，
        # 同一螺旋上相邻 segment 的 node id 单调递增（差为 1），等价于 θ 单调递增。
        # 用 id 排序可避免 atan(y,x) 在 ±π 分支切割处的翻转 bug。
        sort!(common)
        n_lo = common[1]   # θ 较小（id 较小）
        n_hi = common[2]   # θ 较大（id 较大）

        for n in (n_lo, n_hi)
            if !haskey(node_copy, n)
                new_node_count += 1
                extended_node[new_node_count, :] = sub_mesh.node[n, :]
                node_copy[n] = new_node_count
            end
        end
        n_lo_copy = node_copy[n_lo]
        n_hi_copy = node_copy[n_hi]

        # 重写外层 bulk 单元连接
        for col in 1:4
            if bulk_element_new[e_outer, col] == n_lo
                bulk_element_new[e_outer, col] = n_lo_copy
            elseif bulk_element_new[e_outer, col] == n_hi
                bulk_element_new[e_outer, col] = n_hi_copy
            end
        end

        # cohesive 单元几何长度
        x_lo, y_lo = sub_mesh.node[n_lo, 1], sub_mesh.node[n_lo, 2]
        x_hi, y_hi = sub_mesh.node[n_hi, 1], sub_mesh.node[n_hi, 2]
        elem_length = hypot(x_hi - x_lo, y_hi - y_lo)

        coh = CohesiveElement(
            i,
            [n_lo, n_hi, n_hi_copy, n_lo_copy],   # 逆时针
            [n_lo, n_hi],                          # nodes_bottom
            [n_lo_copy, n_hi_copy],                # nodes_top
            elem_length,
            iface,
            e_outer,
            e_inner,
        )
        push!(cohesive_elements, coh)

        # cohesive_to_thermal
        thermal_elem_of_outer = czm_submesh.thermal_elem_map[e_outer]
        @assert thermal_elem_of_outer > 0 "外层单元 $e_outer 的 thermal_elem_map 无效"
        @assert thermal_elem_of_outer <= n_thermal "外层单元 $e_outer thermal_elem_map=$thermal_elem_of_outer 超出热网格范围 $n_thermal"
        cohesive_to_thermal[i] = thermal_elem_of_outer
    end

    # 裁剪 extended_node
    extended_node = extended_node[1:new_node_count, :]

    # Step 4: 组装 CohesiveMesh
    damage_states = [DamageState() for _ in 1:n_cohesive]

    czm_mesh = CohesiveMesh()
    czm_mesh.bulk_mesh = sub_mesh
    czm_mesh.node = extended_node
    czm_mesh.nnode = new_node_count
    czm_mesh.bulk_element = bulk_element_new
    czm_mesh.cohesive_elements = cohesive_elements
    czm_mesh.n_cohesive = n_cohesive
    czm_mesh.n_layers = 2   # 遗留字段：N_type^coh=2；真实面数/重复单元=4，总离散数=4*n_segments
    czm_mesh.node_map = Dict(n => [n, c] for (n, c) in node_copy)
    czm_mesh.interface_nodes = [[]]   # 旧字段，保留兼容
    czm_mesh.damage_states = damage_states
    czm_mesh.czm_submesh = czm_submesh
    # 插值按 .mesh 未合并布局构造（索引算术定位原封不动），行数契约 size(M,1)==mesh.nlen
    # 维持不变（T_czm_nodes 无下游消费者，未合并行序仅影响该返回值）
    czm_mesh.thermal_to_czm = build_thermal_to_czm_interp(thermal_mesh, czm_submesh)
    czm_mesh.cohesive_to_thermal = cohesive_to_thermal

    # 正确性自检（spec §4.3）
    for coh in cohesive_elements
        n_lo, n_hi, n_hi_copy, n_lo_copy = coh.nodes
        @assert czm_mesh.node[n_lo, :] ≈ czm_mesh.node[n_lo_copy, :] atol=1e-12 "副本坐标不一致"
        @assert czm_mesh.node[n_hi, :] ≈ czm_mesh.node[n_hi_copy, :] atol=1e-12 "副本坐标不一致"
        @assert length(unique(coh.nodes)) == 4 "cohesive 单元 4 节点重复"
    end

    return czm_mesh
end
