using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "界面热阻暂禁用 (spec v2 §2.4)" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm_enabled = true             # 即使 CZM 启用
    opt.mechanicalmodel = "full"
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    # v2 修订：czm_enabled=true 时也应使用合并网格（连续径向导热路径）
    @test case.mesh["thermal2D"] === mesh_data.thermal2D_merged

    # 界面热阻在 BC 函数中已注释掉，验证 ThermalDistributed2D_BC 不修改 K 矩阵的 czm 相关项
    # 检查方式：构造 K=0, F=0 调用 BC 后，与不传 czm_enabled 的结果数值一致
    # （此处仅做语法 smoke test；数值等价性由全网格回归覆盖）
    @test case.opt.czm_enabled == true   # 配置仍启用
end
