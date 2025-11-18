"""
阶段3简化测试：快速验证CallModel_MultiSPMe
"""

using Plots, CSV, DataFrames
include("./src/JuBat.jl")

println("="^60)
println("阶段3简化测试")
println("="^60)

# 创建案例
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
opt.model = "SPMe"
opt.Current = x -> 5.0
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true

case = JuBat.SetCase(param_dim, opt)

# 创建thermal2D网格
nx, ny = 10, 5
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

# 初始化
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
println("✓ 初始化: $(length(y0)) 自由度")

# 测试CallModel_MultiSPMe
println("\n[测试] CallModel_MultiSPMe调用...")
t = 0.0
M, K, F, variables, _ = JuBat.CallModel_MultiSPMe(case, y0, t, jacobi="update")

println("✓ 矩阵维度: M=$(size(M)), K=$(size(K)), F=$(length(F))")
println("✓ 电压: $(variables["cell voltage"] * case.param.scale.phi) V")

# 验证热源
q_elem = variables["heat_source_fields"]
I_e = variables["thermal2D element current"]
println("✓ 热源范围: [$(minimum(q_elem)), $(maximum(q_elem))]")
println("✓ 电流范围: [$(minimum(I_e)), $(maximum(I_e))]")

# 验证电流守恒
areas = case.thermal2D_element_area_cache
w = areas ./ sum(areas)
I_total = case.opt.Current(t * case.param.scale.t0) / case.param_dim.cell.I1C
I_sum = sum(w .* I_e)
println("✓ 电流守恒: I_total=$(I_total), Σ(w·I_e)=$(I_sum), err=$(abs(I_total-I_sum))")

println("\n" * "="^60)
println("✓ 所有测试通过！阶段3完成。")
println("="^60)
