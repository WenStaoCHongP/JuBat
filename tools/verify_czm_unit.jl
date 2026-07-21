# verify_czm_unit.jl
# 单元级别验证：2个实体单元 + 1个内聚力单元

using LinearAlgebra
using SparseArrays
using Printf
using Plots
using DelimitedFiles

# 包含JuBat模块
# 假设脚本在 tools/ 目录下运行
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# ========================================================================
# 辅助函数
# ========================================================================

"""
双线性本构解析解 (从 verify_czm_analytical.jl 复制)
"""
function analytical_bilinear(δ, δ_0, δ_c, K, σ_max)
    if δ <= 0
        return 0.0, K * δ, 0.0  # 压缩：纯弹性
    elseif δ <= δ_0
        return 0.0, K * δ, 0.0  # 弹性阶段
    elseif δ >= δ_c
        return 1.0, 0.0, 1.0    # 完全断裂
    else
        # 软化阶段
        D = δ_c * (δ - δ_0) / (δ * (δ_c - δ_0))
        T = (1 - D) * K * δ
        return D, T, D
    end
end

"""
简单测试问题：
    [固定]----[单元1]----[CZM]----[单元2]----[位移u]
"""
function test_2element_tension(czm_model::String)
    println("\n" * "="^60)
    println("测试2：双实体单元-内聚力单元拉伸验证 (Mode I)")
    println("当前模型: $(czm_model)")
    println("="^60)

    # ----------------------------------------------------------------
    # 1. 定义几何与网格
    # ----------------------------------------------------------------
    # 两个1x1的方形单元
    # 节点编号：
    # 单元1: 1(0,0), 2(1,0), 3(1,1), 4(0,1)
    # 单元2: 5(1,0), 6(2,0), 7(2,1), 8(1,1)
    # 界面: 单元1的右边(2-3) 与 单元2的左边(5-8)
    # CZM连接: Bottom=[3,2], Top=[8,5] (法向指向+x)

    L_ele = 100e-6 # 100 um
    nnode = 8
    node_coords = [
        0.0 0.0;          # 1
        L_ele 0.0;        # 2
        L_ele L_ele;      # 3
        0.0 L_ele;        # 4
        L_ele 0.0;        # 5
        2*L_ele 0.0;      # 6
        2*L_ele L_ele;    # 7
        L_ele L_ele       # 8
    ]

    # 实体单元 (Q4)
    bulk_elements = [
        1 2 3 4;
        5 6 7 8
    ]

    # 手动构建 CohesiveMesh
    czm_mesh = JuBat.CohesiveMesh()
    czm_mesh.node = node_coords
    czm_mesh.nnode = nnode
    czm_mesh.bulk_element = bulk_elements
    
    # ----------------------------------------------------------------
    # 2. 定义内聚力单元
    # ----------------------------------------------------------------
    # Bottom: 3->2 (Nodes 3, 2). Vector (0, -1). Normal (1, 0)
    # Top: 8->5 (Nodes 8, 5). Vector (0, -1).
    # CohesiveElement struct: nodes=[n1, n2, n3, n4]
    # nodes_bottom=[n1, n2], nodes_top=[n4, n3]
    # n1=3, n2=2, n4=8, n3=5
    # nodes = [3, 2, 5, 8]
    
    ce = JuBat.CohesiveElement(
        1,              # id
        [3, 2, 5, 8],   # nodes
        [3, 2],         # nodes_bottom
        [8, 5],         # nodes_top
        L_ele,          # length
        :PE_PCC,        # interface_type  # TODO Chunk 3 重写
        0,              # host_outer_elem
        0               # host_inner_elem
    )
    
    czm_mesh.cohesive_elements = [ce]
    czm_mesh.n_cohesive = 1
    czm_mesh.damage_states = [JuBat.DamageState()]

    # ----------------------------------------------------------------
    # 3. 材料参数
    # ----------------------------------------------------------------
    # 实体参数
    E = 100.0e9  # 100 GPa
    ν = 0.3
    
    # 内聚力参数
    σ_max = 100.0e6  # 100 MPa
    δ_c = 1.0e-5     # 10 um
    δ_0 = 1.0e-6     # 1 um (K = 100 MPa / 1 um = 1e14)
    K_init = σ_max / δ_0
    
    cohesive_params = JuBat.Cohesive()
    cohesive_params.czm_model = czm_model
    # cohesive_params.σ_max_n = σ_max                   # TODO Chunk 2 Task 2.1
    # cohesive_params.δ_0_n = δ_0                        # TODO Chunk 2 Task 2.1
    # cohesive_params.δ_c_n = δ_c                        # TODO Chunk 2 Task 2.1
    # cohesive_params.K_n = K_init                       # TODO Chunk 2 Task 2.1

    # 设置切向参数（虽然Mode I下理论上不用，但防止数值噪声引其调用）
    # cohesive_params.τ_max_t = σ_max                   # TODO Chunk 2 Task 2.1
    # cohesive_params.δ_0_t = δ_0                        # TODO Chunk 2 Task 2.1
    # cohesive_params.δ_c_t = δ_c                        # TODO Chunk 2 Task 2.1
    # cohesive_params.K_t = K_init                       # TODO Chunk 2 Task 2.1
    cohesive_params.eta = 1.45 # BK准则参数

    println("参数:")
    println("  实体 E = $(E/1e9) GPa")
    println("  CZM σ_max = $(σ_max/1e6) MPa, δ_0 = $(δ_0*1e6) um, δ_c = $(δ_c*1e6) um")

    # ----------------------------------------------------------------
    # 4. 加载循环 (位移控制)
    # ----------------------------------------------------------------
    # 固定单元1左边 (Nodes 1, 4)
    # 拉伸单元2右边 (Nodes 6, 7)
    
    fixed_dofs = [
        1*2-1; 1*2;  # Node 1 x,y
        4*2-1; 4*2;  # Node 4 x,y
        6*2;         # Node 6 y (Allow x movement)
        7*2          # Node 7 y
    ]
    
    prescribed_dofs = [
        6*2-1;       # Node 6 x
        7*2-1        # Node 7 x
    ]
    
    free_dofs = setdiff(1:(nnode*2), [fixed_dofs; prescribed_dofs])
    
    u_max = 1.5 * δ_c
    n_steps = 50
    
    # 存储结果
    results_disp = Float64[]
    results_traction = Float64[]
    results_damage = Float64[]
    
    u_current = zeros(Float64, nnode*2)
    
    # 计算实体刚度矩阵 (线性)
    K_bulk_all = JuBat.assemble_bulk_stiffness(czm_mesh, E, ν)
    
    println("\n开始加载...")
    println("Step | u_applied (um) | Sep (um) | Traction (MPa) | Damage (%) | Iter")
    println("-"^70)
    
    for step in 1:n_steps
        u_app = u_max * step / n_steps
        
        # 更新强制位移
        u_current[prescribed_dofs] .= u_app
        
        # Newton-Raphson 求解
        converged = false
        iter = 0
        max_iter = 20
        tol = 1e-6
        
        while !converged && iter < max_iter
            iter += 1
            
            # 1. 组装 CZM 部分 (依赖当前位移 u_current)
            K_czm, f_int_czm, seps, tracts = JuBat.assemble_czm_system(czm_mesh, u_current, cohesive_params)
            
            # 2. 组装 Total K 和 F_int
            # F_int_bulk = K_bulk * u
            f_int_bulk = K_bulk_all * u_current
            
            R = - (f_int_bulk + f_int_czm) # 改：F_ext = 0 (在自由节点上)，残差 R = F_ext - F_int
            # 注意：对于强制位移节点，F_ext不为0，我们在求解时只考虑自由自由度
            
            K_total = K_bulk_all + K_czm
            
            # 3. 求解修正量 du
            # 仅针对自由自由度求解
            K_free = K_total[free_dofs, free_dofs]
            R_free = R[free_dofs]
            
            du_free = K_free \ R_free
            
            u_current[free_dofs] += du_free
            
            # 检查收敛
            norm_R = norm(R_free)
            norm_du = norm(du_free)
            
            if norm_R < tol || norm_du < 1e-10
                converged = true
            end
        end
        
        # 记录每步结果
        damage_state = czm_mesh.damage_states[1]
        δ_n, δ_t = JuBat.compute_separation(czm_mesh.cohesive_elements[1], czm_mesh.node, u_current)
        T_n, T_t, D = JuBat.bilinear_traction(δ_n, δ_t, damage_state, cohesive_params; update=true)
        
        # 保存
        push!(results_disp, δ_n)
        push!(results_traction, T_n)
        push!(results_damage, D)
        
        if step % 5 == 0 || step == 1
             @printf("%4d | %12.4f | %8.4f | %12.4f | %8.2f | %d\n", 
                step, u_app*1e6, δ_n*1e6, T_n/1e6, D*100, iter)
        end
    end
    
    # ----------------------------------------------------------------
    # 5. 结果分析与验证
    # ----------------------------------------------------------------
    println("\n验证结果:")

    # 绘制曲线
    p1 = plot(results_disp*1e6, results_traction/1e6, 
        label="Simulated", xlabel="Separation (um)", ylabel="Traction (MPa)",
        title="Traction-Separation Curve", lw=2, marker=:circle, markersize=3)
    
    # 理论曲线
    δ_theory = range(0, 1.5*δ_c, length=100)
    T_theory = [analytical_bilinear(d, δ_0, δ_c, K_init, σ_max)[2] for d in δ_theory]
    plot!(p1, δ_theory*1e6, T_theory/1e6, label="Analytical", linestyle=:dash, color=:red, lw=2)

    p2 = plot(results_disp*1e6, results_damage*100,
        label="Damage", xlabel="Separation (um)", ylabel="Damage (%)",
        title="Damage Evolution", lw=2, color=:red)
    
    p = plot(p1, p2, layout=(1,2), size=(1000, 400))
    output_fig = "czm_unit_verification_$(czm_model).png"
    savefig(p, output_fig)
    println("  已保存结果图像: $(output_fig)")
    
    # 检查峰值应力
    max_T = maximum(results_traction)
    err_max = abs(max_T - σ_max) / σ_max
    passed_max = err_max < 0.05
    @printf("  峰值应力: %.2f MPa (理论 %.2f MPa) -> 误差 %.2f%% [%s]\n", 
        max_T/1e6, σ_max/1e6, err_max*100, passed_max ? "PASS" : "FAIL")

    # 检查完全失效时的位移
    # 理论上 δ_c 时 T 降为 0
    # 找到最后一个 T > 0 的点
    idx_failure = findfirst(x -> x >= 0.99, results_damage)
    if isnothing(idx_failure)
        println("  警告: 未达到完全失效状态")
    else
        δ_fail = results_disp[idx_failure]
        @printf("  失效分离: %.2f um (理论 %.2f um)\n", δ_fail*1e6, δ_c*1e6)
    end
    
    # 检查弹性刚度
    # 取第一个点
    K_measured = results_traction[1] / results_disp[1]
    err_K = abs(K_measured - K_init) / K_init
    passed_K = err_K < 0.05
    @printf("  初始刚度: %.2e (理论 %.2e) -> 误差 %.2f%% [%s]\n", 
        K_measured, K_init, err_K*100, passed_K ? "PASS" : "FAIL")

    return passed_max && passed_K
end

# 运行测试
println("\n运行model1验证...")
passed_model1 = test_2element_tension("model1")
println("\n运行mix验证...")
passed_mix = test_2element_tension("mix")

println("\n验证结果汇总：")
println("  model1: ", passed_model1 ? "PASS" : "FAIL")
println("  mix: ", passed_mix ? "PASS" : "FAIL")
