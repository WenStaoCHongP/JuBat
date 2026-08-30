using Test
using Statistics

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# 层分辨宏观应力：耦合在线收割（export_macro_stress）与固体按需工具
# （thermal_diffusion_stress_2D，mesh_bonded 层分辨实现）的契约测试。
# 物理基准：PE.Omega = +7.88e-7（2026-08-29 由 −7.28e-7 更正）——放电时 NE 脱锂收缩
# （受拉）、PE 嵌锂膨胀（受压），涂层间出现拉压交替；刚性集流体承担约束拉力。

function layer_stress_fixture(; nθ::Int=8, czm_enabled::Bool=false, fix_inner::Bool=true)
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm.enabled = czm_enabled
    opt.czm.fix_inner = fix_inner
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=nθ, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    case.mech = JuBat.MechState(case.czm_mesh)
    return case
end

function discharge_variables(case)
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    return Dict{String, Union{Array{Float64},Float64}}(
        "T_nodes" => fill(case.param.cell.T0, mesh.nlen),
        "thermal2D element soc_n" => fill(case.param.NE.cs0 - 0.02, ne),
        "thermal2D element soc_p" => fill(case.param.PE.cs0 + 0.02, ne))
end

function interval_history_case()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Current = _ -> 5.0
    opt.Nn = 4; opt.Ns = 3; opt.Np = 4
    opt.Nrn = 4; opt.Nrp = 4
    opt.gsorder = 2
    opt.time = [0.0, 2.0]
    opt.dt = [0.5, 0.5]
    opt.dtType = "constant"
    opt.outputType = "auto"
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm.enabled = true
    opt.czm.fix_inner = false
    opt.czm.update_interval = 2
    case = JuBat.SetCase(param_dim, opt)
    md = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=8, czm_enabled=true, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, md)
    case.czm_mesh = JuBat.create_czm_mesh(md.czm_submesh, case.mesh["thermal2D"], case.param)
    case.mech = JuBat.MechState(case.czm_mesh)
    return case
end

@testset "固体工具：零扰动应力严格为零" begin
    case = layer_stress_fixture()
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    vars = Dict{String, Union{Array{Float64},Float64}}(
        "T_nodes" => fill(case.param.cell.T0, mesh.nlen),
        "thermal2D element soc_n" => fill(case.param.NE.cs0, ne),
        "thermal2D element soc_p" => fill(case.param.PE.cs0, ne))
    out = JuBat.thermal_diffusion_stress_2D(case, vars)
    @test maximum(abs, out["diffusion stress vonMises [Pa]"]) == 0.0
    @test maximum(abs, out["displacement x [m]"]) < 1e-12
    @test maximum(abs, out["displacement y [m]"]) < 1e-12
end

@testset "固体工具：放电载荷下涂层一拉一压" begin
    case = layer_stress_fixture(; fix_inner=false)
    out = JuBat.thermal_diffusion_stress_2D(case, discharge_variables(case))
    sxx = out["diffusion stress xx [Pa]"]
    mt = case.czm_mesh.czm_submesh.material_type
    ne_mean = mean(sxx[mt .== :NE])
    pe_mean = mean(sxx[mt .== :PE])
    # 涂层间拉压交替是稳定物理图案；集流体均值受端部集中主导、对离散敏感，不断言其排序
    @test ne_mean > 0 > pe_mean
    @test all(isfinite, sxx)
    @test maximum(abs, out["displacement x [m]"]) > 0
end

@testset "耦合收割：export_macro_stress 写入历史列" begin
    case = layer_stress_fixture(; czm_enabled=true, fix_inner=false)
    case.mech = JuBat.MechState(case.czm_mesh)
    vars = discharge_variables(case)
    T_nodes = vars["T_nodes"]
    JuBat.update_czm_damage!(case, vars, T_nodes)
    ne_czm = size(case.czm_mesh.bulk_element, 1)
    vh = Dict{String, Union{Array{Float64},Float64}}(
        "diffusion stress xx" => zeros(ne_czm, 2),
        "diffusion stress yy" => zeros(ne_czm, 2),
        "diffusion stress xy" => zeros(ne_czm, 2),
        "diffusion stress vonMises" => zeros(ne_czm, 2))
    JuBat.export_macro_stress(case, vars, vh, 2, T_nodes)
    sxx = vh["diffusion stress xx"][:, 2]
    mt = case.czm_mesh.czm_submesh.material_type
    @test maximum(abs, sxx) > 0
    @test mean(sxx[mt .== :NE]) > mean(sxx[mt .== :PE])
end

@testset "Solve：CZM 间隔更新保持最近有效应力且快照时间对齐" begin
    case = interval_history_case()
    snapshots = JuBat.CZMSnapshot[]
    result = JuBat.Solve(case; czm_snapshots=snapshots)
    times = result["time [s]"]
    sxx = result["diffusion stress xx [Pa]"]

    @test !isempty(snapshots)
    @test all(s -> any(t -> isapprox(t, s.time_s; atol=1e-10), times), snapshots)

    first_update_col = findfirst(t -> isapprox(t, snapshots[1].time_s; atol=1e-10), times)
    @test first_update_col !== nothing
    if first_update_col !== nothing
        held_col = first_update_col + 1
        @test held_col <= length(times)
        if held_col <= length(times)
            @test !any(s -> isapprox(s.time_s, times[held_col]; atol=1e-10), snapshots)
            @test sxx[:, held_col] == sxx[:, first_update_col]
            @test maximum(abs, sxx[:, held_col]) > 0.0
        end
    end
    @test all(isfinite, sxx)
end

@testset "门控：czm-off 时 export_macro_stress 原样返回" begin
    case = layer_stress_fixture(; czm_enabled=false)
    vars = discharge_variables(case)
    vh = Dict{String, Union{Array{Float64},Float64}}()
    out = JuBat.export_macro_stress(case, vars, vh, 1, vars["T_nodes"])
    @test out === vh
    @test !haskey(vh, "diffusion stress xx")
end

@testset "E_coat 缺参拦截（LGM50）" begin
    param_lgm = JuBat.ChooseCell("LG M50")
    case_lgm = JuBat.SetCase(param_lgm, JuBat.Option())
    vars = Dict{String, Union{Array{Float64},Float64}}()
    try
        JuBat.thermal_diffusion_stress_2D(case_lgm, vars)
        error("FAIL: 期望抛 AssertionError 但未抛")
    catch e
        @test e isa AssertionError
        @test occursin("E_coat", e.msg)
    end
end
