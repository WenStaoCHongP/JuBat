using Test

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# spec 2026-08-20-core-collapse-mechanics-design.md §5 兼容性契约：堆芯塌陷力学建模的所有新能力 opt-in，默认关。
# 本测试锁定"默认关"这一契约本身——任何把默认值改为开的提交都会在此失败。

@testset "堆芯塌陷子选项默认全关（spec §5）" begin
    opt = JuBat.Option()

    @test opt.czm_geo_nonlinear == false
    @test opt.czm_winding_prestress == false
    @test opt.czm_j2_plasticity == false
    @test opt.czm_phi_bond == false
    @test opt.czm_continuous_feedback == false

    # Batch 8 预留：当前无消费者，只锁定默认值
    @test opt.czm_friction_mu == 0.10
end

@testset "czm_enabled 是唯一主开关（spec §5）" begin
    # czm_enabled=true 且子选项全关时，子选项仍为 false——
    # 主开关不得隐式打开任何新能力
    opt = JuBat.Option()
    opt.czm_enabled = true

    @test opt.czm_enabled == true
    @test opt.czm_geo_nonlinear == false
    @test opt.czm_winding_prestress == false
    @test opt.czm_j2_plasticity == false
    @test opt.czm_phi_bond == false
    @test opt.czm_continuous_feedback == false
end

@testset "子选项可显式开启（关键字构造与字段赋值两条路径）" begin
    opt_kw = JuBat.Option(czm_geo_nonlinear=true)
    @test opt_kw.czm_geo_nonlinear == true
    @test opt_kw.czm_j2_plasticity == false

    opt_set = JuBat.Option()
    opt_set.czm_j2_plasticity = true
    @test opt_set.czm_j2_plasticity == true
    @test opt_set.czm_geo_nonlinear == false
end
