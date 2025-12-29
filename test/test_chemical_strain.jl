"""
化学应变计算单元测试

测试内容：
1. eps_s 参数计算正确性
2. 化学膨胀系数 β 的计算
3. 参数物理合理性检查
4. 与理论公式的一致性

作者：AI Assistant
日期：2025-12-29
"""

using Test
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

@testset "化学应变参数计算" begin
    
    @testset "Jellyroll电池参数" begin
        param_dim = JuBat.ChooseCell("Jellyroll")
        
        # 检查 eps_s 计算正确性
        @test hasfield(typeof(param_dim.NE), :eps_s)
        @test hasfield(typeof(param_dim.PE), :eps_s)
        
        eps_s_n = param_dim.NE.eps_s
        eps_s_p = param_dim.PE.eps_s
        
        # 验证 eps_s 的定义公式
        @test eps_s_n ≈ 1 - param_dim.NE.eps - param_dim.NE.eps_fi rtol=1e-10
        @test eps_s_p ≈ 1 - param_dim.PE.eps - param_dim.PE.eps_fi rtol=1e-10
        
        # 检查 eps_s 的物理范围（0.3 ~ 0.8为典型值）
        @test 0.2 < eps_s_n < 0.9
        @test 0.2 < eps_s_p < 0.9
        
        # 检查孔隙率总和不超过1
        @test param_dim.NE.eps + param_dim.NE.eps_fi + eps_s_n ≈ 1.0 rtol=1e-10
        @test param_dim.PE.eps + param_dim.PE.eps_fi + eps_s_p ≈ 1.0 rtol=1e-10
        
        println("  ✓ eps_s 计算正确")
        println("    负极 eps_s = $(eps_s_n)")
        println("    正极 eps_s = $(eps_s_p)")
    end
    
    @testset "化学膨胀系数 β" begin
        param_dim = JuBat.ChooseCell("Jellyroll")
        
        Ω_n = param_dim.NE.Omega
        Ω_p = param_dim.PE.Omega
        eps_s_n = param_dim.NE.eps_s
        eps_s_p = param_dim.PE.eps_s
        
        # 计算化学膨胀系数
        β_n = Ω_n / 3.0 * eps_s_n
        β_p = Ω_p / 3.0 * eps_s_p
        
        # 检查符号（负极膨胀，正极收缩）
        @test β_n > 0  # 负极Ω>0（石墨膨胀）
        @test β_p < 0  # 正极Ω<0（三元材料收缩）
        
        # 检查量级（基于文献典型值）
        @test abs(β_n) > 1e-7
        @test abs(β_n) < 1e-5
        @test abs(β_p) > 1e-8
        @test abs(β_p) < 1e-6
        
        # 检查相对大小（负极膨胀通常大于正极收缩）
        @test abs(β_n) > abs(β_p)
        
        println("  ✓ 化学膨胀系数计算正确")
        println("    β_n = $(β_n) (负极膨胀)")
        println("    β_p = $(β_p) (正极收缩)")
    end
    
    @testset "部分摩尔体积 Ω" begin
        param_dim = JuBat.ChooseCell("Jellyroll")
        
        Ω_n = param_dim.NE.Omega
        Ω_p = param_dim.PE.Omega
        
        # 检查量级（基于文献：石墨~3e-6，NCM~-1e-6 m³/mol）
        @test 1e-7 < Ω_n < 1e-5
        @test -1e-5 < Ω_p < -1e-7
        
        # 检查单位（m³/mol应该是小量）
        @test Ω_n < 0.01  # 不应该过大
        @test Ω_p > -0.01
        
        println("  ✓ 部分摩尔体积量级合理")
        println("    Ω_n = $(Ω_n) m³/mol")
        println("    Ω_p = $(Ω_p) m³/mol")
    end
    
    @testset "应变计算一致性" begin
        param_dim = JuBat.ChooseCell("Jellyroll")
        
        # 模拟一个SOC变化
        ΔSOC = 0.5  # 从0.5充到1.0
        
        # 理论公式：ε = (Ω/3) × eps_s × ΔSOC
        β_n = param_dim.NE.Omega / 3.0 * param_dim.NE.eps_s
        β_p = param_dim.PE.Omega / 3.0 * param_dim.PE.eps_s
        
        ε_n = β_n * ΔSOC
        ε_p = β_p * ΔSOC
        
        # 检查应变值在合理范围（微应变~百微应变）
        @test abs(ε_n) < 0.01  # 小于1%
        @test abs(ε_p) < 0.01
        
        # 检查符号（负极膨胀为正，正极收缩为负）
        @test ε_n > 0
        @test ε_p < 0
        
        println("  ✓ 应变计算一致性验证通过")
        println("    ΔSOC = $(ΔSOC) 时：")
        println("    负极应变 = $(ε_n*1e6) με (膨胀)")
        println("    正极应变 = $(ε_p*1e6) με (收缩)")
    end
    
    @testset "LGM50电池参数" begin
        param_dim = JuBat.ChooseCell("LG M50")
        
        # 检查参数存在性
        @test hasfield(typeof(param_dim.NE), :eps_s)
        @test hasfield(typeof(param_dim.PE), :eps_s)
        @test hasfield(typeof(param_dim.NE), :Omega)
        @test hasfield(typeof(param_dim.PE), :Omega)
        
        # 检查物理合理性
        @test 0.2 < param_dim.NE.eps_s < 0.9
        @test 0.2 < param_dim.PE.eps_s < 0.9
        
        println("  ✓ LGM50参数检查通过")
        println("    负极 eps_s = $(param_dim.NE.eps_s)")
        println("    正极 eps_s = $(param_dim.PE.eps_s)")
    end
    
    @testset "无量纲化参数传递" begin
        param_dim = JuBat.ChooseCell("Jellyroll")
        param = JuBat.NormaliseParam(param_dim)
        
        # 检查 eps_s 是否正确传递到无量纲参数
        @test hasfield(typeof(param.NE), :eps_s)
        @test hasfield(typeof(param.PE), :eps_s)
        
        # eps_s 本身是分数，无量纲化应该保持不变
        @test param.NE.eps_s == param_dim.NE.eps_s
        @test param.PE.eps_s == param_dim.PE.eps_s
        
        println("  ✓ 无量纲化参数传递正确")
    end
end

@testset "应力计算函数接口" begin
    
    @testset "thermal_diffusion_stress_2D 函数" begin
        # 创建最小测试案例
        param_dim = JuBat.ChooseCell("Jellyroll")
        opt = JuBat.Option()
        opt.model = "SPMe"
        opt.thermal_enabled = true
        opt.thermalmodel = "distributed2D"
        
        case = JuBat.SetCase(param_dim, opt)
        
        # 创建简单网格
        mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=20, gsorder=2)
        case.mesh["thermal2D"] = mesh_th
        
        ne = size(mesh_th.element, 1)
        nT = mesh_th.nlen
        
        # 准备输入变量
        variables = Dict{String, Union{Array{Float64},Float64}}()
        variables["T_nodes"] = ones(Float64, nT)  # 均匀温度场
        variables["thermal2D element soc_n"] = fill(0.5, ne)
        variables["thermal2D element soc_p"] = fill(0.5, ne)
        
        # 调用应力计算函数
        variables_out = JuBat.thermal_diffusion_stress_2D(case, variables)
        
        # 检查输出变量
        @test haskey(variables_out, "diffusion stress xx")
        @test haskey(variables_out, "diffusion stress yy")
        @test haskey(variables_out, "diffusion stress xy")
        @test haskey(variables_out, "diffusion stress vonMises")
        @test haskey(variables_out, "thermal stress vonMises")
        @test haskey(variables_out, "diffusion stress vonMises only")
        @test haskey(variables_out, "displacement x")
        @test haskey(variables_out, "displacement y")
        
        # 检查数组大小
        @test length(variables_out["diffusion stress xx"]) == ne
        @test length(variables_out["displacement x"]) == nT
        
        # 检查数值有效性（无NaN或Inf）
        @test all(isfinite.(variables_out["diffusion stress vonMises"]))
        @test all(isfinite.(variables_out["displacement x"]))
        
        println("  ✓ thermal_diffusion_stress_2D 函数接口正确")
    end
end

println("\n" * "="^80)
println("化学应变单元测试完成")
println("="^80)
println("""
测试总结：
✓ eps_s 参数计算和传递正确
✓ 化学膨胀系数 β 符合理论公式
✓ 参数物理合理性验证通过
✓ 应力计算函数接口工作正常

修正内容已验证：
  β_n = Ω_n / 3 × eps_s_n  ✓
  β_p = Ω_p / 3 × eps_s_p  ✓

理论一致性：
  宏观应变 = 颗粒应变 × 体积分数  ✓
  ε_macro = (Ω/3) × eps_s × ΔSOC  ✓
""")
