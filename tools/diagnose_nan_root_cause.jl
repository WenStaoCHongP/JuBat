"""
诊断工具：分析 NaN 产生的根本原因

根据错误信息：
1. 第 110 行 CallModel 成功（V = 3.966 V 正常）
2. 第 165 行初始求解步骤产生 NaN
3. NaN 出现在热自由度（thermal_nan = 6962）

关键：分析第 165 行的矩阵求解为何失败
y_c = (M_old - K_old * dt_init) \ (M_old * y0[vc] + F_old * dt_init)
其中 dt_init = 1e-8
"""

using LinearAlgebra, SparseArrays, Statistics, Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function diagnose_nan_root_cause()
    println("="^80)
    println("NaN 根因诊断：矩阵求解数值稳定性分析")
    println("="^80)
    
    # 设置参数（与用户相同）
    param_dim = JuBat.ChooseCell("Jellyroll")
    param_dim.cell.v_l = 2.5
    param_dim.cell.v_h = 4.2
    param_dim.tab.theta_pos = [0.0]
    param_dim.tab.theta_neg = [20.0 * π]
    
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.Current = x -> 5.0
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.time = [0.0, 60.0]
    opt.dt = [0.5, 10]
    
    case = JuBat.SetCase(param_dim, opt)
    
    # 创建网格
    nθ = 80
    mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=nθ, gsorder=2)
    case.mesh["thermal2D"] = mesh_th
    
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    println("\n[1] 网格信息:")
    println("  热单元数 ne = $ne")
    println("  热节点数 nT = $nT")
    
    # 初始化
    println("\n[2] 执行初始化...")
    y0 = JuBat.ModelInitialisation_MultiSPMe(case)
    layout = case.multi_spme_layout
    n_chem = layout["n_chem"]
    
    println("  状态向量维度:")
    println("    总维度 = $(length(y0))")
    println("    电化学自由度 n_chem = $n_chem")
    println("    热自由度 nT = $nT")
    println("    总和 = $(n_chem + nT) (应该等于总维度)")
    
    # 检查 y0
    nan_in_y0 = sum(.!isfinite.(y0))
    if nan_in_y0 > 0
        @error "❌ 初始状态向量已包含 NaN/Inf！" count=nan_in_y0
        return
    else
        println("    ✓ y0 不包含 NaN/Inf")
    end
    
    # 调用 CallModel（模拟第 110 行）
    println("\n[3] 调用 CallModel（模拟 Solve.jl 第 110 行）...")
    t0 = 0.0
    
    try
        M_old, K_old, F_old, variables, y_phi = JuBat.CallModel(case, y0, t0, jacobi="update")
        
        println("  ✓ CallModel 成功")
        println("  返回矩阵维度:")
        println("    M_old: $(size(M_old))")
        println("    K_old: $(size(K_old))")
        println("    F_old: $(length(F_old))")
        
        V_init = variables["cell voltage"] * case.param.scale.phi
        println("  初始电压 V = $(round(V_init, digits=4)) V")
        
        # 分析矩阵的结构
        println("\n[4] 分析矩阵结构...")
        
        # 提取化学和热部分
        M_chem = M_old[1:n_chem, 1:n_chem]
        K_chem = K_old[1:n_chem, 1:n_chem]
        
        M_therm = M_old[(n_chem+1):end, (n_chem+1):end]
        K_therm = K_old[(n_chem+1):end, (n_chem+1):end]
        
        println("  化学部分:")
        if nnz(M_chem) > 0
            m_vals = nonzeros(M_chem)
            println("    M_chem 元素范围: [$(minimum(m_vals)), $(maximum(m_vals))]")
        end
        if nnz(K_chem) > 0
            k_vals = nonzeros(K_chem)
            println("    K_chem 元素范围: [$(minimum(k_vals)), $(maximum(k_vals))]")
        end
        
        println("\n  热部分:")
        if nnz(M_therm) > 0
            mt_vals = nonzeros(M_therm)
            println("    M_therm 元素范围: [$(minimum(mt_vals)), $(maximum(mt_vals))]")
            println("    M_therm 对角元素示例（前10个）:")
            for i in 1:min(10, nT)
                println("      M_therm[$i,$i] = $(M_therm[i,i])")
            end
        end
        
        if nnz(K_therm) > 0
            kt_vals = nonzeros(K_therm)
            println("    K_therm 元素范围: [$(minimum(kt_vals)), $(maximum(kt_vals))]")
            println("    K_therm 对角元素示例（前10个）:")
            for i in 1:min(10, nT)
                println("      K_therm[$i,$i] = $(K_therm[i,i])")
            end
            
            # 检查是否有极大的对角元素（惩罚法）
            diag_kt = [K_therm[i,i] for i in 1:nT]
            max_diag = maximum(abs.(diag_kt))
            if max_diag > 1e10
                @warn "⚠️  K_therm 包含极大的对角元素（惩罚法）" max=max_diag
                penalty_nodes = findall(abs.(diag_kt) .> 1e10)
                println("    惩罚节点数量: $(length(penalty_nodes))")
                if length(penalty_nodes) > 0
                    println("    惩罚节点索引（前10个）: $(penalty_nodes[1:min(10, length(penalty_nodes))])")
                    println("    对应的对角元素值:")
                    for i in penalty_nodes[1:min(5, length(penalty_nodes))]
                        println("      K_therm[$i,$i] = $(K_therm[i,i])")
                    end
                end
            end
        end
        
        # 分析初始求解步骤（模拟第 165 行）
        println("\n[5] 模拟初始求解步骤（Solve.jl 第 165 行）...")
        dt_init = 1e-8
        
        println("  dt_init = $dt_init")
        
        # 计算左端矩阵 A = M_old - K_old * dt_init
        A = M_old - K_old * dt_init
        
        println("  计算 A = M_old - K_old * dt_init")
        
        # 检查热部分的对角元素
        A_therm = A[(n_chem+1):end, (n_chem+1):end]
        diag_A_therm = [A_therm[i,i] for i in 1:nT]
        
        println("\n  热部分 A_therm 的对角元素分析:")
        println("    范围: [$(minimum(diag_A_therm)), $(maximum(diag_A_therm))]")
        println("    均值: $(mean(diag_A_therm))")
        println("    标准差: $(std(diag_A_therm))")
        
        # 检查是否有负数或接近零的对角元素
        negative_diag = sum(diag_A_therm .< 0)
        near_zero_diag = sum(abs.(diag_A_therm) .< 1e-10)
        
        if negative_diag > 0
            @error "❌ 发现负数对角元素！" count=negative_diag
            neg_indices = findall(diag_A_therm .< 0)
            println("    负数对角元素索引（前10个）: $(neg_indices[1:min(10, length(neg_indices))])")
            for i in neg_indices[1:min(5, length(neg_indices))]
                println("      A_therm[$i,$i] = $(A_therm[i,i])")
                println("        M_therm[$i,$i] = $(M_therm[i,i])")
                println("        K_therm[$i,$i] = $(K_therm[i,i])")
                println("        K_therm[$i,$i] * dt_init = $(K_therm[i,i] * dt_init)")
            end
            
            println("\n  ❌ 根因分析：")
            println("    M_therm[i,i] - K_therm[i,i] * dt_init < 0")
            println("    原因：K_therm[i,i] 由于惩罚法过大（例如 1e12）")
            println("          K_therm[i,i] * dt_init = 1e12 * 1e-8 = 1e4")
            println("          如果 M_therm[i,i] < 1e4，则 A_therm[i,i] < 0")
            println("    结果：矩阵变为病态或奇异，求解产生 NaN")
        end
        
        if near_zero_diag > 0
            @warn "⚠️  发现接近零的对角元素" count=near_zero_diag
        end
        
        # 尝试求解（如果安全的话）
        if negative_diag == 0 && near_zero_diag < nT * 0.1
            println("\n  尝试求解...")
            vc = 1:size(M_old,1)
            b = M_old * y0[vc] + F_old * dt_init
            
            try
                y_c = A \ b
                nan_count = sum(.!isfinite.(y_c))
                
                if nan_count > 0
                    @error "❌ 求解产生 NaN！" count=nan_count
                    
                    # 分析哪些自由度变成了 NaN
                    chem_nan = sum(.!isfinite.(y_c[1:n_chem]))
                    therm_nan = sum(.!isfinite.(y_c[(n_chem+1):end]))
                    
                    println("    化学自由度 NaN: $chem_nan / $n_chem")
                    println("    热自由度 NaN: $therm_nan / $nT")
                    
                    if therm_nan > 0
                        println("\n    热自由度中 NaN 的位置（前10个）:")
                        nan_thermal_indices = findall(.!isfinite.(y_c[(n_chem+1):end]))
                        for i in nan_thermal_indices[1:min(10, length(nan_thermal_indices))]
                            println("      y_c[$(n_chem + i)] = $(y_c[n_chem + i])")
                        end
                    end
                else
                    println("    ✓ 求解成功，未产生 NaN")
                    println("    y_c 范围: [$(minimum(y_c)), $(maximum(y_c))]")
                end
            catch err
                @error "❌ 求解失败" exception=(err, catch_backtrace())
            end
        else
            println("\n  ⚠️  跳过求解（矩阵明显病态）")
        end
        
        # 给出解决方案
        println("\n" * "="^80)
        println("解决方案")
        println("="^80)
        
        if negative_diag > 0
            println("\n✅ 根因确认：惩罚法导致的矩阵病态")
            println("\n方案 1：降低惩罚值（推荐）")
            println("  在 testexample.jl 中添加：")
            println("    opt.tab_penalty = 1e6  # 或 1e8")
            println("\n  解释：")
            println("    - 默认 penalty = 1e12 太大")
            println("    - penalty * dt_init = 1e12 * 1e-8 = 1e4")
            println("    - 如果 M_therm[i,i] < 1e4，则 A[i,i] < 0")
            println("    - 降低 penalty 可以避免这个问题")
            
            println("\n方案 2：增大初始时间步长")
            println("  修改 Solve.jl 第 163 行：")
            println("    dt_init = 1e-6  # 从 1e-8 改为 1e-6 或更大")
            println("\n  解释：")
            println("    - 增大 dt_init 可以减小 K * dt 的影响")
            println("    - 但这只是权宜之计，不解决根本问题")
            
            println("\n方案 3：使用迭代求解器")
            println("  修改求解方法，使用预条件迭代求解器")
            println("  （需要修改源码，较复杂）")
        end
        
    catch err
        @error "❌ CallModel 失败" exception=(err, catch_backtrace())
        return
    end
    
    println("\n" * "="^80)
end

# 运行诊断
diagnose_nan_root_cause()
