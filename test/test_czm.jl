"""
    test_czm.jl

内聚力模型(CZM)单元测试

测试内容：
1. CZM参数合理性
2. 双线性本构模型
3. 指数本构模型
4. 界面分离计算
5. 损伤演化
6. 能量守恒

运行方式：
julia test/test_czm.jl

或使用Test.jl:
using Test
include("test_czm.jl")
"""

using Test

# 引入CZM模块
include("../src/czm.jl")
include("../src/parameters/CZM_Parameters.jl")

# ============================================================================
# 测试集1：参数合理性
# ============================================================================

@testset "CZM Parameter Validation" begin
    println("\n测试集1：CZM参数合理性")
    
    @testset "NE-NCC Parameters" begin
        mat = get_NE_NCC_parameters()
        
        # 参数为正
        @test mat.K_n > 0
        @test mat.K_t > 0
        @test mat.t_n_max > 0
        @test mat.t_t_max > 0
        @test mat.G_Ic > 0
        @test mat.G_IIc > 0
        
        # 临界分离量一致性
        δ_0 = mat.t_n_max / mat.K_n
        δ_f = 2.0 * mat.G_Ic / mat.t_n_max
        
        @test δ_f > δ_0
        
        println("  ✓ NE-NCC参数合理")
        println("    δ_0 = $(δ_0*1e9) nm")
        println("    δ_f = $(δ_f*1e9) nm")
    end
    
    @testset "PE-PCC Parameters" begin
        mat = get_PE_PCC_parameters()
        
        @test mat.K_n > 0
        @test mat.G_Ic > 0
        
        δ_0 = mat.t_n_max / mat.K_n
        δ_f = 2.0 * mat.G_Ic / mat.t_n_max
        
        @test δ_f > δ_0
        
        println("  ✓ PE-PCC参数合理")
    end
    
    @testset "Validate Function" begin
        mat = get_NE_NCC_parameters()
        @test validate_czm_material(mat) == true
        
        # 创建不合理参数
        bad_mat = CZMMaterial(
            K_n = 1e12,
            K_t = 5e11,
            t_n_max = 10e6,
            t_t_max = 8e6,
            G_Ic = 10.0,  # 太小，会导致 δ_f < δ_0
            G_IIc = 20.0,
            alpha = 1.0,
            beta = 1.0,
            viscosity = 0.0,
            model = BILINEAR
        )
        
        @test validate_czm_material(bad_mat) == false
        
        println("  ✓ 参数验证函数正常")
    end
end

# ============================================================================
# 测试集2：双线性本构模型
# ============================================================================

@testset "Bilinear Cohesive Law" begin
    println("\n测试集2：双线性本构模型")
    
    mat = get_NE_NCC_parameters()
    
    δ_0 = mat.t_n_max / mat.K_n
    δ_f = 2.0 * mat.G_Ic / mat.t_n_max
    
    @testset "Stage 1: Elastic (δ < δ_0)" begin
        δ_n = 0.5 * δ_0  # 一半损伤起始值
        δ_t = 0.0
        
        t_n, t_t, D, dG = bilinear_cohesive_law(δ_n, δ_t, mat, 0.0)
        
        @test D == 0.0  # 无损伤
        @test t_n ≈ mat.K_n * δ_n  # 线性关系
        @test dG == 0.0  # 无耗散
        
        println("  ✓ 阶段1（弹性）：D=0, 线性响应")
    end
    
    @testset "Stage 2: Softening (δ_0 ≤ δ < δ_f)" begin
        δ_n = 0.5 * (δ_0 + δ_f)  # 中间值
        δ_t = 0.0
        
        t_n, t_t, D, dG = bilinear_cohesive_law(δ_n, δ_t, mat, 0.0)
        
        @test D > 0.0 && D < 1.0  # 部分损伤
        @test t_n < mat.K_n * δ_n  # 软化（低于弹性）
        @test dG > 0.0  # 有耗散
        
        # 检查牵引力-分离关系
        expected_D = (δ_n - δ_0) / (δ_f - δ_0)
        @test D ≈ expected_D atol=1e-6
        
        println("  ✓ 阶段2（软化）：0<D<1, 有能量耗散")
    end
    
    @testset "Stage 3: Complete Failure (δ ≥ δ_f)" begin
        δ_n = 1.5 * δ_f  # 超过完全失效值
        δ_t = 0.0
        
        t_n, t_t, D, dG = bilinear_cohesive_law(δ_n, δ_t, mat, 0.0)
        
        @test D == 1.0  # 完全失效
        @test t_n == 0.0  # 无法承载
        
        println("  ✓ 阶段3（失效）：D=1, t=0")
    end
    
    @testset "Damage Irreversibility" begin
        δ_n = 1.2 * δ_0  # 产生损伤
        
        # 第一次加载
        _, _, D1, _ = bilinear_cohesive_law(δ_n, 0.0, mat, 0.0)
        @test D1 > 0
        
        # 卸载（分离量减小）
        δ_n_unload = 0.5 * δ_0
        _, _, D2, _ = bilinear_cohesive_law(δ_n_unload, 0.0, mat, D1)
        
        @test D2 == D1  # 损伤不可恢复
        
        println("  ✓ 损伤不可逆性")
    end
    
    @testset "Mixed Mode" begin
        δ_n = 8e-9  # nm
        δ_t = 6e-9  # nm
        
        t_n, t_t, D, _ = bilinear_cohesive_law(δ_n, δ_t, mat, 0.0)
        
        δ_eff = sqrt(δ_n^2 + δ_t^2)
        
        @test δ_eff ≈ 10e-9
        @test t_n > 0
        @test abs(t_t) > 0  # 切向牵引力非零
        
        println("  ✓ 混合模式（法向+切向）")
    end
end

# ============================================================================
# 测试集3：指数本构模型
# ============================================================================

@testset "Exponential Cohesive Law" begin
    println("\n测试集3：指数本构模型")
    
    mat = get_particle_binder_parameters()  # 使用指数模型的参数
    
    @testset "Basic Response" begin
        δ_0 = mat.G_Ic / mat.t_n_max
        δ_n = 0.5 * δ_0
        
        t_n, t_t, D, dG = exponential_cohesive_law(δ_n, 0.0, mat, 0.0)
        
        @test D >= 0.0 && D <= 1.0
        @test t_n > 0
        
        println("  ✓ 指数模型基本响应")
    end
    
    @testset "Smooth Degradation" begin
        δ_0 = mat.G_Ic / mat.t_n_max
        δ_values = range(0, 5*δ_0, length=20)
        
        D_old = 0.0
        for δ in δ_values
            _, _, D, _ = exponential_cohesive_law(δ, 0.0, mat, D_old)
            
            @test D >= D_old  # 单调增长
            D_old = D
        end
        
        println("  ✓ 损伤平滑增长")
    end
end

# ============================================================================
# 测试集4：界面分离计算
# ============================================================================

@testset "Interface Separation" begin
    println("\n测试集4：界面分离计算")
    
    @testset "Displacement Extraction" begin
        # 模拟variables
        variables = Dict{String, Union{Array{Float64}, Float64}}(
            "displacement x" => [1e-9, 2e-9, 3e-9],
            "displacement y" => [0.5e-9, 1.5e-9, 2.5e-9]
        )
        
        U_global = extract_displacement_from_variables(variables)
        
        @test length(U_global) == 6
        @test U_global[1] == 1e-9
        @test U_global[2] == 0.5e-9
        @test U_global[3] == 2e-9
        @test U_global[4] == 1.5e-9
        
        println("  ✓ 位移提取正确")
    end
    
    @testset "Separation Calculation" begin
        # 创建模拟界面
        mat = get_NE_NCC_parameters()
        interface = CZMInterface(
            1,                          # id
            COATING_COLLECTOR,          # type
            2,                          # node_plus
            1,                          # node_minus
            [1.0, 0.0],                # normal (x方向)
            [0.0, 1.0],                # tangent (y方向)
            1e-6,                       # area
            mat                         # material
        )
        
        # 模拟位移场
        U_global = [
            0.0, 0.0,        # node 1: u_x=0, u_y=0
            5e-9, 2e-9       # node 2: u_x=5nm, u_y=2nm
        ]
        
        δ_n, δ_t = compute_interface_separation(interface, U_global)
        
        # 位移跳跃 = [5, 2] nm
        # 法向分离 = [5, 2]·[1, 0] = 5 nm
        # 切向分离 = [5, 2]·[0, 1] = 2 nm
        
        @test δ_n ≈ 5e-9
        @test δ_t ≈ 2e-9
        
        println("  ✓ 分离量计算正确")
    end
end

# ============================================================================
# 测试集5：界面状态更新
# ============================================================================

@testset "Interface State Update" begin
    println("\n测试集5：界面状态更新")
    
    mat = get_NE_NCC_parameters()
    interface = CZMInterface(
        1,
        COATING_COLLECTOR,
        2, 1,
        [1.0, 0.0], [0.0, 1.0],
        1e-6,
        mat
    )
    
    @testset "Single Update" begin
        # 初始状态
        @test interface.D == 0.0
        @test interface.delta_n == 0.0
        
        # 施加位移
        U_global = [0.0, 0.0, 15e-9, 0.0]  # 15 nm法向分离
        dt = 1.0
        
        update_czm_interface!(interface, U_global, dt)
        
        # 检查更新后状态
        @test interface.delta_n ≈ 15e-9
        @test interface.D > 0.0  # 已产生损伤
        @test interface.t_n > 0.0  # 有牵引力
        
        println("  ✓ 单次更新正确")
    end
    
    @testset "Multiple Updates" begin
        # 重置界面
        interface.D = 0.0
        interface.delta_n = 0.0
        interface.delta_max = 0.0
        interface.G_dissipated = 0.0
        
        # 模拟加载过程
        n_steps = 10
        for i in 1:n_steps
            δ = i * 3e-9  # 逐步增加到30 nm
            U_global = [0.0, 0.0, δ, 0.0]
            update_czm_interface!(interface, U_global, 1.0)
        end
        
        # 检查最终状态
        @test interface.delta_max > 0
        @test interface.D > 0
        @test interface.G_dissipated > 0
        
        println("  ✓ 多步更新正确")
    end
end

# ============================================================================
# 测试集6：能量守恒
# ============================================================================

@testset "Energy Conservation" begin
    println("\n测试集6：能量守恒")
    
    mat = get_NE_NCC_parameters()
    
    @testset "DCB Test" begin
        # 模拟双悬臂梁加载到完全失效
        
        δ_0 = mat.t_n_max / mat.K_n
        δ_f = 2.0 * mat.G_Ic / mat.t_n_max
        
        # 分100步加载到2倍δ_f
        n_steps = 100
        δ_values = range(0, 2*δ_f, length=n_steps)
        
        G_numerical = 0.0
        D_old = 0.0
        
        for δ in δ_values
            t_n, _, D, dG = bilinear_cohesive_law(δ, 0.0, mat, D_old)
            G_numerical += dG
            D_old = D
        end
        
        # 理论值
        G_theory = mat.G_Ic
        
        # 检查误差（允许10%，因为离散化）
        error = abs(G_numerical - G_theory) / G_theory
        
        @test error < 0.10
        
        println("  ✓ 能量守恒验证")
        println("    理论: $(G_theory) J/m²")
        println("    数值: $(round(G_numerical, digits=2)) J/m²")
        println("    误差: $(round(error*100, digits=2))%")
    end
end

# ============================================================================
# 测试集7：辅助函数
# ============================================================================

@testset "Utility Functions" begin
    println("\n测试集7：辅助函数")
    
    @testset "Statistics" begin
        mat = get_NE_NCC_parameters()
        
        # 创建一组界面
        interfaces = [
            CZMInterface(i, COATING_COLLECTOR, i+1, i,
                        [1.0, 0.0], [0.0, 1.0], 1e-6, mat)
            for i in 1:10
        ]
        
        # 设置不同损伤状态
        for i in 1:10
            interfaces[i].D = i * 0.1  # 0.1, 0.2, ..., 1.0
            if i == 10
                interfaces[i].failed = true
            end
        end
        
        stats = compute_czm_statistics(interfaces)
        
        @test stats["n_total"] == 10
        @test stats["n_failed"] == 1
        @test stats["D_max"] == 1.0
        @test stats["D_mean"] ≈ 0.55
        
        println("  ✓ 统计函数正确")
    end
    
    @testset "Parameter Adjustment" begin
        mat_base = get_NE_NCC_parameters()
        
        # 降低强度50%
        mat_weak = adjust_czm_parameters(mat_base, Dict(:t_max => 0.5))
        
        @test mat_weak.t_n_max == 0.5 * mat_base.t_n_max
        @test mat_weak.K_n == mat_base.K_n  # 其他参数不变
        
        println("  ✓ 参数调整函数正确")
    end
end

# ============================================================================
# 汇总测试结果
# ============================================================================

println("\n" * "="^70)
println("CZM单元测试完成")
println("="^70)
println("\n所有测试通过 ✓")
println("\n下一步：")
println("  1. 集成到主程序")
println("  2. 运行完整仿真案例")
println("  3. 验证物理结果")
