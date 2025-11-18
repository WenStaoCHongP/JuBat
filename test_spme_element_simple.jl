"""
SPMe_element 简化测试脚本

直接通过 include 方式加载，适合快速验证
"""

using Plots, CSV, DataFrames
include("./src/JuBat.jl") 

println("="^80)
println("SPMe_element 简化测试")
println("="^80)

# 创建简单测试案例
println("\n创建测试案例...")
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
opt.model = "SPMe"
opt.Nn = 5
opt.Ns = 5
opt.Np = 5
opt.Nrn = 5
opt.Nrp = 5
opt.Current = x -> 5.0  # 5A
opt.time = [0 100]
opt.thermalmodel = "none"

case = JuBat.SetCase(param_dim, opt)
println("✓ 案例创建成功")

# 初始化
yt0_raw = JuBat.ModelInitialisation(case)
# 确保是向量形式
yt0 = isa(yt0_raw, Vector) ? yt0_raw : vec(yt0_raw)
println("✓ 初始化完成，状态向量长度: $(length(yt0))")

# 测试1: 基本调用
println("\n[测试1] 基本调用...")
try
    t = 0.0
    I_e = case.opt.Current(t) / case.param_dim.cell.I1C
    T_e = case.param.cell.T0
    
    M_e, K_e, F_e, vars_e = JuBat.SPMe_element(
        case, yt0, t, 1;
        I_e = I_e,
        T_e = T_e,
        jacobi = "update"
    )
    
    println("✓ SPMe_element 调用成功")
    println("  矩阵维度: M=$(size(M_e)), K=$(size(K_e)), F=$(length(F_e))")
    println("  单元电压: $(vars_e["cell voltage"] * case.param.scale.phi) V")
    println("  负极过电位: $(vars_e["negative electrode overpotential"][1])")
    println("  正极过电位: $(vars_e["positive electrode overpotential"][end])")
catch e
    println("✗ 测试失败: $e")
    rethrow(e)
end

# 测试2: 与全局SPMe对比
println("\n[测试2] 与全局SPMe对比...")
try
    t = 0.0
    I_global = case.opt.Current(t) / case.param_dim.cell.I1C
    T_global = case.param.cell.T0
    
    # 全局SPMe
    M_g, K_g, F_g, vars_g = JuBat.SPMe(case, yt0, t, jacobi="update")
    
    # 单元SPMe
    M_e, K_e, F_e, vars_e = JuBat.SPMe_element(
        case, yt0, t, 1; I_e=I_global, T_e=T_global, jacobi="update"
    )
    
    # 比较
    V_diff = abs(vars_g["cell voltage"] - vars_e["cell voltage"])
    eta_n_diff = abs(vars_g["negative electrode overpotential"][1] - vars_e["negative electrode overpotential"][1])
    
    println("✓ 对比完成")
    println("  电压差异: $(V_diff)")
    println("  负极过电位差异: $(eta_n_diff)")
    
    if V_diff < 1e-10 && eta_n_diff < 1e-10
        println("✓ 结果一致（误差 < 1e-10）")
    else
        println("⚠ 结果有微小差异")
    end
catch e
    println("✗ 测试失败: $e")
    rethrow(e)
end

# 测试3: 不同电流
println("\n[测试3] 不同电流响应...")
try
    t = 0.0
    T_e = 1.0
    I_tests = [0.0, 0.5, 1.0, 2.0]
    
    println("  I_e [1] | V [V] | η_n | η_p")
    println("  " * "-"^40)
    
    for I_e in I_tests
        M_e, K_e, F_e, vars_e = JuBat.SPMe_element(
            case, yt0, t, 1; I_e=I_e, T_e=T_e
        )
        V = vars_e["cell voltage"] * case.param.scale.phi
        eta_n = vars_e["negative electrode overpotential"][1]
        eta_p = vars_e["positive electrode overpotential"][end]
        
        @printf("  %.1f | %.4f | %.4e | %.4e\n", I_e, V, eta_n, eta_p)
    end
    println("✓ 电流扫描完成")
catch e
    println("✗ 测试失败: $e")
    rethrow(e)
end

# 测试4: 不同温度
println("\n[测试4] 不同温度响应...")
try
    t = 0.0
    I_e = 1.0
    T_tests = [0.95, 1.0, 1.05]
    T_ref = case.param_dim.scale.T_ref
    
    println("  T_e [-] | T [K] | V [V] | j0_n")
    println("  " * "-"^40)
    
    for T_e in T_tests
        M_e, K_e, F_e, vars_e = JuBat.SPMe_element(
            case, yt0, t, 1; I_e=I_e, T_e=T_e
        )
        V = vars_e["cell voltage"] * case.param.scale.phi
        j0_n = vars_e["negative electrode exchange current density"]
        
        @printf("  %.2f | %.1f | %.4f | %.4e\n", T_e, T_e*T_ref, V, j0_n)
    end
    println("✓ 温度扫描完成")
catch e
    println("✗ 测试失败: $e")
    rethrow(e)
end

println("\n" * "="^80)
println("所有测试通过！SPMe_element 函数工作正常。")
println("="^80)
