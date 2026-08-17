"""
    create_unit_czm_strip(param; width=nothing, y0=1.0, gsorder=2)

平直 8 层单元条带：8×Q4 + 经 create_czm_mesh 得到 4×COH2D4。
层序自下而上 PE→PCC→PE→SP→NE→NCC→NE→SP。底边 y=y0>0，节点 id 沿 x 递增。
"""
function create_unit_czm_strip(param; width=nothing, y0::Float64=1.0, gsorder::Int=2)
    layer_materials = [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    heights = [
        param.PE.thickness, param.PCC.thickness,
        param.PE.thickness, param.SP.thickness,
        param.NE.thickness, param.NCC.thickness,
        param.NE.thickness, param.SP.thickness,
    ]
    H = sum(heights)
    W = width === nothing ? H : Float64(width)
    y_interfaces = Vector{Float64}(undef, 9)
    y_interfaces[1] = y0
    for i in 1:8
        y_interfaces[i + 1] = y_interfaces[i] + heights[i]
    end

    # 节点：行优先，每行左→右（id 沿 x 递增）
    # 行 r=0..8 对应 y_interfaces[r+1]；列 c=0..1 对应 x=0,W
    nnode = 18
    node = zeros(Float64, nnode, 2)
    for r in 0:8
        for c in 0:1
            idx = r * 2 + c + 1
            node[idx, 1] = c == 0 ? 0.0 : W
            node[idx, 2] = y_interfaces[r + 1]
        end
    end

    # Q4：下层左、下层右、上层右、上层左（逆时针，inner=bottom）
    element = zeros(Int64, 8, 4)
    for e in 1:8
        bl = (e - 1) * 2 + 1   # bottom-left
        br = bl + 1
        tl = e * 2 + 1         # top-left
        tr = tl + 1
        element[e, :] = [bl, br, tr, tl]
    end

    gs = GetGS(element, node, gsorder, "Q4")
    bulk = Mesh("Q4", 2, node, nnode, element, gs)
    material_type = copy(layer_materials)
    winding_turn = ones(Int, 8)
    thermal_elem_map = ones(Int, 8)
    submesh = CzmSubmesh(bulk, material_type, winding_turn, thermal_elem_map, Tuple{Int,Int}[])

    # 哑热网格：单 Q4 覆盖条带 bbox（供 build_thermal_to_czm_interp）
    pad = 1e-6 * max(W, H)
    tn = zeros(Float64, 4, 2)
    tn[1, :] = [-pad, y0 - pad]
    tn[2, :] = [W + pad, y0 - pad]
    tn[3, :] = [W + pad, y_interfaces[end] + pad]
    tn[4, :] = [-pad, y_interfaces[end] + pad]
    te = reshape(Int64[1, 2, 3, 4], 1, 4)
    tgs = GetGS(te, tn, gsorder, "Q4")
    dummy_thermal = Mesh("Q4", 2, tn, 4, te, tgs)

    czm_mesh = create_czm_mesh(submesh, dummy_thermal, param)

    # ---- 硬断言（方案 C）----
    @assert czm_mesh.n_cohesive == 4 "expected 4 cohesive, got $(czm_mesh.n_cohesive)"
    types = [e.interface_type for e in czm_mesh.cohesive_elements]
    @assert count(==(:PE_PCC), types) == 2 "PE_PCC count"
    @assert count(==(:NE_NCC), types) == 2 "NE_NCC count"
    for coh in czm_mesh.cohesive_elements
        n_lo, n_hi, n_hi_c, n_lo_c = coh.nodes
        @assert czm_mesh.node[n_lo, :] ≈ czm_mesh.node[n_lo_c, :] atol=1e-12
        @assert czm_mesh.node[n_hi, :] ≈ czm_mesh.node[n_hi_c, :] atol=1e-12
        @assert length(unique(coh.nodes)) == 4
        # 外层 bulk 共边应为副本
        e_out = coh.host_outer_elem
        outer_nodes = Set(czm_mesh.bulk_element[e_out, :])
        @assert n_lo_c in outer_nodes && n_hi_c in outer_nodes
        @assert !(n_lo in outer_nodes) && !(n_hi in outer_nodes)
    end

    bottom_nodes = [1, 2]
    # 顶层（层 8）在 create_czm_mesh 后：外层若被重写则取 bulk_element[8,:] 上边
    top_row_orig = [17, 18]
    e8 = czm_mesh.bulk_element[8, :]
    ys = [czm_mesh.node[n, 2] for n in e8]
    ytop = maximum(ys)
    top_nodes_after = sort([Int(n) for n in e8 if abs(czm_mesh.node[n, 2] - ytop) < 1e-14])

    # 集流体 bulk 节点（create_czm_mesh 后，含界面副本）
    pcc_nodes = sort(unique(Int.(czm_mesh.bulk_element[2, :])))
    ncc_nodes = sort(unique(Int.(czm_mesh.bulk_element[6, :])))

    meta = (
        y_interfaces = y_interfaces,
        bottom_nodes = bottom_nodes,
        top_nodes = top_row_orig,
        top_nodes_after_czm = top_nodes_after,
        pcc_nodes = pcc_nodes,
        ncc_nodes = ncc_nodes,
        cohesive_ids = collect(1:4),
        interface_types = types,
        layer_materials = layer_materials,
        width = W,
        heights = heights,
    )
    return czm_mesh, meta
end
