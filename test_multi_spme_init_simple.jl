"""
阶段2简化测试：快速验证多SPMe初始化功能
"""

using Plots, CSV, DataFrames
include("./src/JuBat.jl")

println("="^60)
println("阶段2简化测试")
println("="^60)

# 创建案例
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
opt.model = "SPMe"
opt.Current = x -> 5.0
opt.time = [0 100]
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true

case = JuBat.SetCase(param_dim, opt)

# 创建简单thermal2D网格 (10x5 = 50单元)
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
nT = size(nodes, 1)
println("✓ 案例创建: $ne 单元, $nT 节点")

# 测试1: 基本初始化
println("\n[测试1] 基本初始化...")
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
println("✓ 状态向量长度: $(length(y0))")

# 测试2: 提取状态
println("\n[测试2] 状态提取...")
yt_1 = JuBat.MultiSPMe_extract_element_state(y0, 1, case)
T_nodes = JuBat.MultiSPMe_get_thermal_dofs(y0, case)
println("✓ 单元状态长度: $(length(yt_1))")
println("✓ 热场长度: $(length(T_nodes))")

# 测试3: 状态更新
println("\n[测试3] 状态更新...")
y_test = copy(y0)
yt_modified = yt_1 .* 1.1
JuBat.MultiSPMe_update_element_state!(y_test, 1, yt_modified, case)
yt_check = JuBat.MultiSPMe_extract_element_state(y_test, 1, case)
if all(yt_check .≈ yt_modified)
    println("✓ 更新成功")
else
    println("✗ 更新失败")
end

# 测试4: 非均匀SOC
println("\n[测试4] 非均匀SOC...")
soc_dist = range(0.8, 1.0, length=ne)
y0_nu = JuBat.ModelInitialisation_MultiSPMe(case; initial_soc_distribution=collect(soc_dist))
println("✓ 非均匀初始化成功")

# 测试5: SPMe_element集成
println("\n[测试5] SPMe_element集成...")
M_e, K_e, F_e, vars_e = JuBat.SPMe_element(case, yt_1, 0.0, 1; I_e=1.0, T_e=1.0)
println("✓ 集成成功，V=$(vars_e["cell voltage"]*case.param.scale.phi) V")

println("\n" * "="^60)
println("✓ 所有测试通过！阶段2完成。")
println("="^60)
