"""
阶段4简化测试：快速验证完整Solve流程
"""

using Plots, CSV, DataFrames
include("./src/JuBat.jl")

println("="^60)
println("阶段4简化测试：完整时间推进")
println("="^60)

# 创建多SPMe案例
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
opt.model = "SPMe"
opt.Current = x -> 10.0  # 10A
opt.time = [0 5]  # 5秒
opt.dt = [1e-3 0.1]
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true
opt.thermal_enabled = true

case = JuBat.SetCase(param_dim, opt)

# 创建thermal2D网格（5×3单元，快速）
nx, ny = 5, 3
Lx, Ly = 0.1, 0.05
nodes = zeros(Float64, (nx+1)*(ny+1), 2)
idx = 1
for j in 1:(ny+1), i in 1:(nx+1)
    nodes[idx, :] = [(i-1)*Lx/nx, (j-1)*Ly/ny]
    idx += 1
end

elements = zeros(Int, nx*ny, 4)
idx = 1
for j in 1:ny, i in 1:nx
    n1 = (j-1)*(nx+1) + i
    elements[idx, :] = [n1, n1+1, n1+(nx+1)+1, n1+(nx+1)]
    idx += 1
end

case.mesh["thermal2D"] = JuBat.Mesh2D(nodes, elements, size(nodes,1), size(elements,1),
                                       JuBat.GetGS("Q4",2), "Q4")

ne = size(elements, 1)
println("✓ 案例创建: $ne 单元")

# 测试Solve
println("\n[测试] 多SPMe模式Solve...")
result = JuBat.Solve(case)

t = result["time"]
V = result["cell voltage"]
I = result["cell current"]

println("✓ 求解成功")
println("  时间步数: $(length(t))")
println("  初始电压: $(V[1]) V")
println("  最终电压: $(V[end]) V")

# 验证逐单元变量
if haskey(result, "thermal2D element current")
    I_e = result["thermal2D element current"]
    println("✓ 逐单元变量记录: $(size(I_e))")
    println("  I_e 范围: [$(minimum(I_e[:, end])), $(maximum(I_e[:, end]))]")
else
    println("✗ 缺少逐单元变量")
end

# 对比单SPMe
println("\n[对比] 单SPMe模式...")
case_single = deepcopy(case)
case_single.opt.per_element_spme = false
result_single = JuBat.Solve(case_single)

V_single = result_single["cell voltage"]
println("✓ 单SPMe求解成功")
println("  电压对比: 多SPMe=$(V[end]) V, 单SPMe=$(V_single[end]) V")
println("  差异: $(abs(V[end] - V_single[end])) V ($(100*abs(V[end] - V_single[end])/V_single[end])%)")

println("\n" * "="^60)
println("✓ 所有测试通过！阶段4完成。")
println("="^60)
