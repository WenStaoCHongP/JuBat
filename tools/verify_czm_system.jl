# tools/verify_czm_system.jl

using LinearAlgebra
using SparseArrays
using Printf
using Plots
using DelimitedFiles
using Statistics

# 包含JuBat模块 (假设脚本在 tools/ 目录下运行)
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# 引入 CZM 相关 (如果 JuBat 没有导出这两个，则手动 includes)
# 通常 JuBat.jl 已经 include 了 czm.jl

"""
系统级验证 Stage 3: 纯热/化学工况
"""
function verify_system_level()
    println("============================================================")
    println("Stage 3: System Level Verification (Thermal/Chemical)")
    println("============================================================")

    # 1. 准备模型（果冻卷）
    # ----------------------------------------------------------------
    # 使用较少的圈数或扇区，简化问题以便观察
    param_dim = JuBat.ChooseCell("Jellyroll")
    
    # 修改参数仅用于测试：单层材料均匀化
    # 我们可以暂时把正负极活性材料参数设为一样，方便解析解计算
    # 但为了真实性，我们尽量使用 compute_effective_properties 之后的有效值
    
    opt = JuBat.Option()
    opt.gsorder = 2
    
    # 生成热网格 (2D截面)
    # nθ 控制周向划分密度
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=60, gsorder=2)
    thermal_mesh = mesh_data.Jellyroll_czm
    
    # 转换网格为 CZM 网格
    # create_czm_mesh 会自动识别重合节点并在其间建立 Cohesive Elements
    czm_mesh = JuBat.create_czm_mesh(thermal_mesh, param_dim)
    
    println("Mesh Stats:")
    println("  Total Nodes: $(czm_mesh.nnode)")
    println("  Bulk Elements: $(size(czm_mesh.bulk_element, 1))")
    println("  Cohesive Elements: $(czm_mesh.n_cohesive)")
    
    # 计算有效性质 (Effective Properties)
    # 根据 mechanical.jl 里的逻辑，这里需要 E_eff, nu_eff, alpha_eff
    # 假设均质化参数
    E_eff = 100.0e9 # 100 GPa
    ν_eff = 0.3
    α_eff = 1.0e-5  # 1e-5 /K
    β_n = 0.0 # 暂时只测纯热，化学膨胀设为0
    β_p = 0.0
    
    # CZM 参数 (硬)
    # 系统级验证如果是测自由膨胀，CZM 应该足够强，保持界面不分离
    # 或者如果界面是完美的，CZM 刚度应很大
    cohesive_params = JuBat.Cohesive()
    # cohesive_params.σ_max_n = 50.0e6                                                  # TODO Chunk 2 Task 2.1
    # cohesive_params.δ_0_n = 1.0e-9                                                    # TODO Chunk 2 Task 2.1
    # cohesive_params.δ_c_n = 1.0e-6                                                    # TODO Chunk 2 Task 2.1
    # cohesive_params.K_n = cohesive_params.σ_max_n / cohesive_params.δ_0_n # 5e16 Very Stiff  # TODO Chunk 2 Task 2.1
    # cohesive_params.K_t = cohesive_params.K_n                                         # TODO Chunk 2 Task 2.1
    
    # ----------------------------------------------------------------
    # Case 1: 纯热膨胀 (Free Thermal Expansion)
    # 预测：
    #   径向位移 u_r = α * ΔT * r
    #   切向位移 u_t = 0 (对于轴对称/中心对称)
    #   实际上螺旋结构不是完美轴对称，但近似应该符合。
    # ----------------------------------------------------------------
    println("\nCase 1: Pure Thermal Expansion (Free)")
    
    ΔT = 10.0 # 10 K
    T_nodes = zeros(Float64, czm_mesh.nnode)
    fill!(T_nodes, ΔT) # 均匀温升
    
    Δsoc_n = zeros(Float64, size(czm_mesh.bulk_element, 1))
    Δsoc_p = zeros(Float64, size(czm_mesh.bulk_element, 1))
    dT_elem = zeros(Float64, size(czm_mesh.bulk_element, 1))
    fill!(dT_elem, ΔT)
    
    # 组装刚度矩阵
    K_bulk = JuBat.assemble_bulk_stiffness(czm_mesh, E_eff, ν_eff)
    
    # 初始位移为0
    ndof = czm_mesh.nnode * 2
    u = zeros(Float64, ndof)
    
    # Newton-Raphson 求解 loop (虽然这里是线性的，但CZM是非线性的)
    for iter in 1:10
        K_czm, f_int_czm, _, _ = JuBat.assemble_czm_system(czm_mesh, u, cohesive_params)
        
        # 热载荷 (F_th)
        # 注意: assemble_thermal_chemical_load 返回的是 F_load
        # 平衡方程: K*u = F_load + F_ext
        # 残差 R = F_load + F_ext - F_int
        # 这里 F_int_bulk = K_bulk * u
        # F_load = \int B^T D epsilon_0
        F_th = JuBat.assemble_thermal_chemical_load(czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p, dT_elem, Δsoc_n, Δsoc_p)
        
        f_int_bulk = K_bulk * u
        
        # 边界条件：固定一个点 + 去除刚体旋转 (无罚函数)
        # 1. 固定参考点 A (最近接近 (0,0))
        dist_sq = [czm_mesh.node[k,1]^2 + czm_mesh.node[k,2]^2 for k in 1:czm_mesh.nnode]
        node_A = argmin(dist_sq)
        pos_A = czm_mesh.node[node_A, :]
        
        # 2. 选一个“接近 x 轴且最远”的点 B，用于约束转动
        # 这里选择 |y| 最小的若干节点中，半径最大的那个
        y_abs = [abs(czm_mesh.node[k,2]) for k in 1:czm_mesh.nnode]
        idx_sorted = sortperm(y_abs)
        candidates = idx_sorted[1:min(20, length(idx_sorted))]
        node_B = candidates[argmax(dist_sq[candidates])]
        pos_B = czm_mesh.node[node_B, :]
        
        if iter == 1
            println("BC Nodes: A=$node_A (Fixed ux,uy), B=$node_B (Fixed uy)")
        end
        
        # 残差
        R = F_th - (f_int_bulk + f_int_czm)
        K_total = K_bulk + K_czm
        
        # 施加BC: 固定点 A 的 ux, uy；点 B 的 uy (消除刚体转动)
        free_dofs = trues(ndof)
        free_dofs[2*node_A-1] = false
        free_dofs[2*node_A] = false
        free_dofs[2*node_B] = false
        
        # 求解增量
        K_free = K_total[free_dofs, free_dofs]
        R_free = R[free_dofs]
        
        du_free = K_free \ R_free
        
        du_full = zeros(Float64, ndof)
        du_full[free_dofs] = du_free
        
        u += du_full
        
        if norm(R_free) < 1e-6
            println("  Iter $iter: Converged (Res Norm = $(norm(R_free)))")
            break
        else
            println("  Iter $iter: Res Norm = $(norm(R_free))")
        end
    end
    
    # 结果验证
    # 理论解调整：相对于固定点 A
    # u_theo(P) = epsilon * (P - A)
    
    # Re-calculate node_A (as it was local to loop scope in previous implementation)
    dist_sq = [czm_mesh.node[k,1]^2 + czm_mesh.node[k,2]^2 for k in 1:czm_mesh.nnode]
    node_A = argmin(dist_sq)
    
    u_sim_vec = reshape(u, 2, :) # 2 x nnode
    pos = czm_mesh.node' # 2 x nnode
    pos_A = czm_mesh.node[node_A, :]
    
    # 相对位置
    rel_pos = pos .- pos_A # 2 x nnode
    
    # 理论位移
    u_theo_vec = α_eff * ΔT .* rel_pos
    
    # 计算误差
    # 距离 A 的距离
    dist_A = [norm(rel_pos[:, i]) for i in 1:czm_mesh.nnode]
    
    # 投影到径向 (相对于 A)
    # 对于每个点 P，定义局部径向 = (P-A) / |P-A|
    u_sim_rad = zeros(Float64, czm_mesh.nnode)
    u_theo_rad = zeros(Float64, czm_mesh.nnode)
    
    for i in 1:czm_mesh.nnode
        if dist_A[i] > 1e-6
            dir = rel_pos[:, i] / dist_A[i]
            u_sim_rad[i] = dot(u_sim_vec[:, i], dir)
            u_theo_rad[i] = dot(u_theo_vec[:, i], dir)
        else
            u_sim_rad[i] = 0.0
            u_theo_rad[i] = 0.0
        end
    end
    
    # 绘图比对
    p1 = scatter(dist_A*1000, u_sim_rad*1e6, label="Simulated", xlabel="Distance from Fixed Node (mm)", ylabel="Radial Disp (um)", title="Thermal Expansion", markersize=2, alpha=0.6)
    plot!(p1, dist_A*1000, u_theo_rad*1e6, label="Analytical", lw=2, color=:red)
    
    savefig(p1, "verify_system_thermal.png")
    println("  Saved plot: verify_system_thermal.png")
    
    # 误差计算
    err_rel = abs.(u_sim_rad - u_theo_rad) ./ (abs.(u_theo_rad) .+ 1e-9)
    # 去除固定点附近的误差
    valid_idx = dist_A .> 1e-4

    println("\nDetailed check of first 10 valid points:")
    println("dist (mm) | u_sim (um) | u_theo (um) | err (%)")
    count = 0
    for i in 1:length(dist_A)
        if valid_idx[i]
            count += 1
            @printf("%.4f | %.4f | %.4f | %.2f\n", dist_A[i]*1000, u_sim_rad[i]*1e6, u_theo_rad[i]*1e6, err_rel[i]*100)
            if count > 10
                break
            end
        end
    end
    
    avg_err = mean(err_rel[valid_idx])
    rmse = sqrt(mean((u_sim_rad[valid_idx] - u_theo_rad[valid_idx]).^2))
    linf = maximum(abs.(u_sim_rad[valid_idx] - u_theo_rad[valid_idx]))
    println("  Average Relative Error (dist > 0.1mm): $(avg_err * 100) %")
    println("  RMSE (um): $(rmse * 1e6)")
    println("  L_inf (um): $(linf * 1e6)")
    
    if avg_err < 0.05
        println("  [PASS] Thermal Expansion Verification")
    else
        println("  [FAIL] Thermal Expansion Verification")
    end

    # ----------------------------------------------------------------
    # Case 2: Mixed-mode sanity check (model1 vs mix should differ)
    # 说明：纯热膨胀是法向主导，model1 与 mix 会趋同。
    # 这里构造一个非零切向分离，直接调用本构以验证差异。
    # ----------------------------------------------------------------
    println("\nCase 2: Mixed-mode sanity check (model1 vs mix)")

    # 选择一个会触发软化且含切向分离的点
    δ_n_test = 6.0e-6   # 6 um
    δ_t_test = 4.0e-6   # 4 um

    # 准备参数（与上面的法向参数一致，切向参数与法向不同以放大差异）
    czm_params_base = JuBat.Cohesive()
    # czm_params_base.σ_max_n = 50.0e6                                              # TODO Chunk 2 Task 2.1
    # czm_params_base.δ_0_n = 1.0e-6                                                # TODO Chunk 2 Task 2.1
    # czm_params_base.δ_c_n = 1.0e-5                                                # TODO Chunk 2 Task 2.1
    # czm_params_base.K_n = czm_params_base.σ_max_n / czm_params_base.δ_0_n         # TODO Chunk 2 Task 2.1
    # czm_params_base.τ_max_t = 30.0e6                                              # TODO Chunk 2 Task 2.1
    # czm_params_base.δ_0_t = 2.0e-6                                                # TODO Chunk 2 Task 2.1
    # czm_params_base.δ_c_t = 2.0e-5                                                # TODO Chunk 2 Task 2.1
    # czm_params_base.K_t = czm_params_base.τ_max_t / czm_params_base.δ_0_t         # TODO Chunk 2 Task 2.1
    czm_params_base.eta = 1.45

    # model1
    czm_params_m1 = deepcopy(czm_params_base)
    czm_params_m1.czm_model = "model1"
    damage_m1 = JuBat.DamageState()
    Tn_m1, Tt_m1, D_m1 = JuBat.bilinear_traction(δ_n_test, δ_t_test, damage_m1, czm_params_m1; update=true)

    # mix
    czm_params_mix = deepcopy(czm_params_base)
    czm_params_mix.czm_model = "mix"
    damage_mix = JuBat.DamageState()
    Tn_mix, Tt_mix, D_mix = JuBat.bilinear_traction(δ_n_test, δ_t_test, damage_mix, czm_params_mix; update=true)

    println("  δ_n = $(δ_n_test*1e6) um, δ_t = $(δ_t_test*1e6) um")
    println("  model1: Tn=$(Tn_m1/1e6) MPa, Tt=$(Tt_m1/1e6) MPa, D=$(D_m1)")
    println("  mix:    Tn=$(Tn_mix/1e6) MPa, Tt=$(Tt_mix/1e6) MPa, D=$(D_mix)")

    if abs(D_m1 - D_mix) > 1e-6 || abs(Tt_m1 - Tt_mix) > 1e3
        println("  [PASS] model1 与 mix 存在差异（混合模式生效）")
    else
        println("  [FAIL] 未观察到差异，请检查参数或本构实现")
    end

    # ----------------------------------------------------------------
    # Case 3: Iteration method check (basic / arc_length / load_substep)
    # 使用同一热载荷，调用 solve_czm_step 验证三种迭代方式能正常运行
    # ----------------------------------------------------------------
    println("\nCase 3: Iteration method check (basic / arc_length / load_substep)")

    ndof = 2 * czm_mesh.nnode
    F_ext = zeros(Float64, ndof)
    u_prev = zeros(Float64, ndof)

    methods = ["basic", "arc_length", "load_substep"]
    results = Dict{String, JuBat.CZMResult}()
    for method in methods
        local_czm = deepcopy(czm_mesh)
        result = first(JuBat.solve_czm_step(
            local_czm, F_ext, E_eff, ν_eff, cohesive_params, param_dim, u_prev;
            α_eff=α_eff, β_n=β_n, β_p=β_p,
            dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n, Δsoc_p_elem=Δsoc_p,
            max_iter=50, tol=1e-6, n_load_steps=10, arc_length_alpha=1.0, iter_method=method
        ))

        results[method] = result
        println("  method=$(method), converged=$(result.converged), iters=$(result.iterations), residual=$(result.residual_norm)")
    end

    # 结果对比：分离位移与损伤
    ref = results["load_substep"]
    for method in methods
        if method == "load_substep"
            continue
        end
        res = results[method]
        d_sep_n = maximum(abs.(res.separation_n .- ref.separation_n))
        d_sep_t = maximum(abs.(res.separation_t .- ref.separation_t))
        d_D = maximum(abs.(res.damage .- ref.damage))
        println("  compare $(method) vs load_substep: max|Δδ_n|=$(d_sep_n), max|Δδ_t|=$(d_sep_t), max|ΔD|=$(d_D)")
    end

end

# 运行
verify_system_level()
