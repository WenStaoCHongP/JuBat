"""
诊断 M 矩阵全零问题

检查：
1. M_np, M_pp, M_el 的各个组件
2. 时间尺度系数
3. 网格初始化
"""

push!(LOAD_PATH, "./src")
using JuBat
using LinearAlgebra
using SparseArrays

println("="^80)
println("M矩阵诊断测试")
println("="^80)

# 创建测试案例
function create_test_case()
    param = LGM50_parameters()
    param_dim = LGM50_parameters_dimensional()
    
    opt = Option()
    opt.model = "SPMe"
    opt.Nn = 5
    opt.Ns = 5
    opt.Np = 5
    opt.Nrn = 5
    opt.Nrp = 5
    opt.gsorder = 2
    opt.time = [0, 3600]
    opt.Current = t -> 5.0
    opt.solveType = "Crank-Nicolson"
    opt.jacobi = "update"
    opt.thermalmodel = "none"
    opt.mechanicalmodel = "none"
    
    case = SetCase(param, param_dim, opt)
    return case
end

case = create_test_case()
yt0 = ModelInitialisation(case)

println("\n[1/4] 检查网格初始化...")
println("  网格信息:")
for (name, mesh) in case.mesh
    println("    - $name: nlen=$(mesh.nlen)")
end

println("\n[2/4] 检查时间尺度参数...")
println("  无量纲参数:")
println("    ts_n = $(case.param.scale.ts_n)")
println("    ts_p = $(case.param.scale.ts_p)")
println("    te = $(case.param.scale.te)")
println("  有量纲参数:")
println("    t0 = $(case.param_dim.scale.t0) s")

println("\n[3/4] 手动构建M矩阵组件...")

# 调用 SPMe_variables 获取浓度场
t = 0.0
variables = SPMe_variables(case, vec(yt0), t)

# 提取高斯点浓度
csn_gs = variables["negative particle concentration at gauss point"]
csp_gs = variables["positive particle concentration at gauss point"]

println("  浓度场检查:")
println("    csn_gs: min=$(minimum(csn_gs)), max=$(maximum(csn_gs)), mean=$(mean(csn_gs))")
println("    csp_gs: min=$(minimum(csp_gs)), max=$(maximum(csp_gs)), mean=$(mean(csp_gs))")

# 构建粒子扩散矩阵
mesh_np = case.mesh["negative particle"]
mesh_pp = case.mesh["positive particle"]
mesh_el = case.mesh["electrolyte"]

param = case.param

# 调用 ElectrodeDiffusion
println("\n  调用 ElectrodeDiffusion...")
M_np_raw, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, 0.0)
M_pp_raw, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, 0.0)

println("    M_np_raw: size=$(size(M_np_raw)), nnz=$(nnz(M_np_raw)), norm=$(norm(Matrix(M_np_raw)))")
println("    M_pp_raw: size=$(size(M_pp_raw)), nnz=$(nnz(M_pp_raw)), norm=$(norm(Matrix(M_pp_raw)))")

# 应用时间尺度
M_np = M_np_raw .* param.scale.ts_n / case.param_dim.scale.t0
M_pp = M_pp_raw .* param.scale.ts_p / case.param_dim.scale.t0

println("    M_np (scaled): nnz=$(nnz(M_np)), norm=$(norm(Matrix(M_np)))")
println("    M_pp (scaled): nnz=$(nnz(M_pp)), norm=$(norm(Matrix(M_pp)))")

# 电解液扩散矩阵
println("\n  调用 ElectrolyteDiffusion...")
M_el_raw, K_el = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, variables)

println("    M_el_raw: size=$(size(M_el_raw)), nnz=$(nnz(M_el_raw)), norm=$(norm(Matrix(M_el_raw)))")

M_el = M_el_raw .* param.scale.te / case.param_dim.scale.t0

println("    M_el (scaled): nnz=$(nnz(M_el)), norm=$(norm(Matrix(M_el)))")

# 装配全局M矩阵
M = blockdiag(M_np, M_pp, M_el)

println("\n[4/4] 全局M矩阵:")
println("  size: $(size(M))")
println("  nnz: $(nnz(M))")
println("  norm: $(norm(Matrix(M)))")
println("  最小元素: $(minimum(Matrix(M)))")
println("  最大元素: $(maximum(Matrix(M)))")

if nnz(M) == 0
    println("\n  ❌ M矩阵确实为全零矩阵！")
    println("  可能原因:")
    println("    1. ts_n/t0 或 ts_p/t0 或 te/t0 接近零")
    println("    2. ElectrodeDiffusion 或 ElectrolyteDiffusion 返回零矩阵")
    println("    3. 网格未正确初始化")
    
    # 检查比值
    println("\n  时间尺度比值:")
    println("    ts_n/t0 = $(param.scale.ts_n / case.param_dim.scale.t0)")
    println("    ts_p/t0 = $(param.scale.ts_p / case.param_dim.scale.t0)")
    println("    te/t0 = $(param.scale.te / case.param_dim.scale.t0)")
else
    println("\n  ✓ M矩阵非零")
end

println("\n" * "="^80)
println("诊断完成")
println("="^80)
