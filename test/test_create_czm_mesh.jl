using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "create_czm_mesh from CzmSubmesh" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm_enabled = true                   # 选择未合并热网格
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    submesh = mesh_data.czm_submesh
    @test submesh !== nothing
    czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)

    @test czm_mesh isa JuBat.CohesiveMesh
    @test czm_mesh.n_cohesive > 0

    # interface_type 取值合法
    for elem in czm_mesh.cohesive_elements
        @test elem.interface_type in (:PE_PCC, :NE_NCC)
        @test elem.host_outer_elem >= 1
        @test elem.host_inner_elem >= 1
    end

    # 2 种 interface_type；每个 8 层重复单元有 4 个真实面；整条螺旋每个周向分段各生成 4 个 cohesive 单元
    n_segments = size(submesh.mesh.element, 1) ÷ 8
    @test czm_mesh.n_cohesive == 4 * n_segments

    # czm_submesh 字段已设置
    @test czm_mesh.czm_submesh === submesh

    # thermal_to_czm 字段已设置（在 Chunk 5 中真正填充，此处先 nothing）
    @test hasfield(JuBat.CohesiveMesh, :thermal_to_czm)

    # cohesive_to_thermal 长度 = n_cohesive，值合法
    @test length(czm_mesh.cohesive_to_thermal) == czm_mesh.n_cohesive
    n_thermal = size(case.mesh["thermal2D"].element, 1)
    @test all(1 <= e <= n_thermal for e in czm_mesh.cohesive_to_thermal)

    # 节点复制 + 重写外层 bulk 正确性自检（spec §4.3）
    for coh in czm_mesh.cohesive_elements
        n_lo, n_hi, n_hi_copy, n_lo_copy = coh.nodes
        # 副本节点坐标与原节点一致
        @test czm_mesh.node[n_lo, :] ≈ czm_mesh.node[n_lo_copy, :] atol=1e-12
        @test czm_mesh.node[n_hi, :] ≈ czm_mesh.node[n_hi_copy, :] atol=1e-12
        # 4 节点不重复
        @test length(unique(coh.nodes)) == 4
    end

    # 外层 bulk 单元的共边位置必须是副本节点
    for coh in czm_mesh.cohesive_elements
        outer_nodes = czm_mesh.bulk_element[coh.host_outer_elem, :]
        n_lo, n_hi, n_hi_copy, n_lo_copy = coh.nodes
        # 副本 n_lo_copy/n_hi_copy 必须在外层单元中
        @test n_lo_copy in outer_nodes
        @test n_hi_copy in outer_nodes
        # 原节点 n_lo/n_hi 不应在外层单元中（已被副本替换）
        @test !(n_lo in outer_nodes)
        @test !(n_hi in outer_nodes)
    end
end

@testset "create_czm_mesh cohesive normal 方向" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

    # 副本节点坐标与原节点一致（top 与 bottom 几何重合），
    # 因此 cohesive 法向由拓扑（host_inner 在内、host_outer 在外）决定，而非节点几何。
    # 用相邻 bulk 单元质心方向作为参考外向，验证 host_outer 比 host_inner 更靠外。
    for coh in czm_mesh.cohesive_elements
        inner_nodes = czm_mesh.bulk_element[coh.host_inner_elem, :]
        outer_nodes = czm_mesh.bulk_element[coh.host_outer_elem, :]
        cx_inner = sum(czm_mesh.node[n, 1] for n in inner_nodes) / 4
        cy_inner = sum(czm_mesh.node[n, 2] for n in inner_nodes) / 4
        cx_outer = sum(czm_mesh.node[n, 1] for n in outer_nodes) / 4
        cy_outer = sum(czm_mesh.node[n, 2] for n in outer_nodes) / 4
        r_inner = hypot(cx_inner, cy_inner)
        r_outer = hypot(cx_outer, cy_outer)
        # host_outer 质心应径向更靠外
        @test r_outer > r_inner

        # 进一步验证：lo/hi 两节点 id 相差 1（同一螺旋上的相邻 segment 节点）
        # 这正是 C1 修复所依赖的前提——id 顺序等价于 θ 顺序。
        n_lo, n_hi, n_hi_copy, n_lo_copy = coh.nodes
        @test n_hi - n_lo == 1
    end
end
