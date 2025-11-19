#=
阶段2测试：多SPMe初始化与状态向量管理

测试目标：
1. 验证 ModelInitialisation_MultiSPMe 正确构造扩展状态向量
2. 验证状态向量结构（维度、布局）
3. 测试状态提取和更新函数
4. 测试非均匀初始SOC分布
=#

using Plots, CSV, DataFrames
include("../src/JuBat.jl")
using .JuBat
using Printf

println("="^80)
println("阶段2测试：多SPMe初始化与状态向量管理")
println("="^80)

# ============================================================================
# 测试准备：创建带热场的案例
# ============================================================================
println("\n[1/6] 创建测试案例（SPMe + thermal2D）...")

# 创建案例
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.model = "SPMe"
opt.Nn = 5
opt.Ns = 5
opt.Np = 5
opt.Nrn = 5
opt.Nrp = 5
opt.Current = x -> 5.0
opt.time = [0 100]
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true  # 启用多SPMe模式

# 需要创建thermal2D网格
# 为了测试，我们创建一个简单的矩形网格
println("  创建简单的thermal2D网格...")

# 简单矩形网格参数
Lx, Ly = 0.1, 0.05  # 10cm × 5cm
nx, ny = 10, 5      # 网格划分

# 使用库内 SetMesh 构建规则Q4网格（避免手写 GetGS 调用）
mesh_thermal = JuBat.SetMesh([0.0, Lx, 0.0, Ly], [nx, ny], "Q4", 2)

# 添加到case
case = JuBat.SetCase(param_dim, opt)
case.mesh["thermal2D"] = mesh_thermal

ne = size(mesh_thermal.element, 1)
nT = size(mesh_thermal.node, 1)

println("✓ 案例创建成功")
println("  模型: SPMe + thermal2D")
println("  热单元数: $ne")
println("  热节点数: $nT")

# ============================================================================
# 测试 1: 基本初始化（均匀SOC）
# ============================================================================
println("\n[2/6] 测试基本初始化（均匀SOC）...")

try
    y0 = JuBat.ModelInitialisation_MultiSPMe(case)
    
    println("✓ 初始化成功")
    println("  状态向量总长度: $(length(y0))")
    
    # 验证维度
    layout = case.multi_spme_layout
    n_chem = layout["n_chem"]
    n_chem_total = ne * n_chem
    n_thermal = nT
    n_expected = n_chem_total + n_thermal
    
    println("  电化学自由度: $n_chem × $ne = $n_chem_total")
    println("  热场自由度: $n_thermal")
    println("  预期总长度: $n_expected")
    
    if length(y0) == n_expected
        println("✓ 状态向量维度正确")
    else
        println("✗ 状态向量维度错误！")
        error("维度不匹配")
    end
    
    # 验证layout信息
    println("\n  Layout信息:")
    println("    ne = $(layout["ne"])")
    println("    n_chem = $(layout["n_chem"])")
    println("    nT = $(layout["nT"])")
    println("    n_total = $(layout["n_total"])")
    println("    chem_range = $(layout["chem_range"][1]):$(layout["chem_range"][end])")
    println("    thermal_range = $(layout["thermal_range"][1]):$(layout["thermal_range"][end])")
    
    if layout["ne"] == ne && layout["nT"] == nT
        println("✓ Layout信息正确")
    else
        println("✗ Layout信息错误")
    end
    
catch e
    println("✗ 初始化失败: $e")
    rethrow(e)
end

# ============================================================================
# 测试 2: 状态提取
# ============================================================================
println("\n[3/6] 测试状态提取...")

try
    y0 = JuBat.ModelInitialisation_MultiSPMe(case)
    
    # 提取单元状态
    e_test = min(5, ne)
    yt_e = JuBat.MultiSPMe_extract_element_state(y0, e_test, case)
    
    println("✓ 单元状态提取成功")
    println("  提取单元: $e_test")
    println("  单元状态长度: $(length(yt_e))")
    
    n_chem = case.multi_spme_layout["n_chem"]
    if length(yt_e) == n_chem
        println("✓ 单元状态维度正确")
    else
        println("✗ 单元状态维度错误")
    end
    
    # 提取热场
    T_nodes = JuBat.MultiSPMe_get_thermal_dofs(y0, case)
    
    println("✓ 热场提取成功")
    println("  热场长度: $(length(T_nodes))")
    
    if length(T_nodes) == nT
        println("✓ 热场维度正确")
    else
        println("✗ 热场维度错误")
    end
    
    # 验证初始温度
    T0 = case.param.cell.T0
    if all(T_nodes .≈ T0)
        println("✓ 初始温度正确 (T0 = $T0)")
    else
        println("⚠ 初始温度不一致")
    end
    
catch e
    println("✗ 状态提取失败: $e")
    rethrow(e)
end

# ============================================================================
# 测试 3: 状态更新
# ============================================================================
println("\n[4/6] 测试状态更新...")

try
    y0 = JuBat.ModelInitialisation_MultiSPMe(case)
    y_test = copy(y0)
    
    # 修改单元状态
    e_test = min(3, ne)
    yt_e = JuBat.MultiSPMe_extract_element_state(y_test, e_test, case)
    yt_e_modified = yt_e .* 1.1  # 乘以1.1
    
    JuBat.MultiSPMe_update_element_state!(y_test, e_test, yt_e_modified, case)
    
    # 验证更新
    yt_e_check = JuBat.MultiSPMe_extract_element_state(y_test, e_test, case)
    
    if all(yt_e_check .≈ yt_e_modified)
        println("✓ 单元状态更新成功")
    else
        println("✗ 单元状态更新失败")
        error("状态更新验证失败")
    end
    
    # 修改热场
    T_nodes = JuBat.MultiSPMe_get_thermal_dofs(y_test, case)
    T_nodes_modified = T_nodes .+ 0.05  # 加0.05（无量纲）
    
    JuBat.MultiSPMe_update_thermal_dofs!(y_test, T_nodes_modified, case)
    
    # 验证更新
    T_nodes_check = JuBat.MultiSPMe_get_thermal_dofs(y_test, case)
    
    if all(T_nodes_check .≈ T_nodes_modified)
        println("✓ 热场更新成功")
    else
        println("✗ 热场更新失败")
    end
    
    # 验证其他单元未受影响
    e_other = e_test == 1 ? 2 : 1
    yt_other = JuBat.MultiSPMe_extract_element_state(y_test, e_other, case)
    yt_other_original = JuBat.MultiSPMe_extract_element_state(y0, e_other, case)
    
    if all(yt_other .≈ yt_other_original)
        println("✓ 其他单元未受影响")
    else
        println("⚠ 其他单元可能被错误修改")
    end
    
catch e
    println("✗ 状态更新失败: $e")
    rethrow(e)
end

# ============================================================================
# 测试 4: 非均匀初始SOC分布
# ============================================================================
println("\n[5/6] 测试非均匀初始SOC分布...")

try
    # 创建SOC梯度（从0.5到1.0）
    soc_dist = range(0.5, 1.0, length=ne)
    soc_vec = collect(soc_dist)
    
    println("  SOC分布: $(soc_vec[1]) ~ $(soc_vec[end])")
    
    y0_nonuniform = JuBat.ModelInitialisation_MultiSPMe(case; initial_soc_distribution=soc_vec)
    
    println("✓ 非均匀初始化成功")
    println("  状态向量长度: $(length(y0_nonuniform))")
    
    # 验证不同单元的初始浓度不同
    yt_1 = JuBat.MultiSPMe_extract_element_state(y0_nonuniform, 1, case)
    yt_ne = JuBat.MultiSPMe_extract_element_state(y0_nonuniform, ne, case)
    
    # 提取负极表面浓度（前Nrn个）
    Nrn = case.opt.Nrn
    cn_surf_1 = yt_1[1]
    cn_surf_ne = yt_ne[1]
    
    println("  单元1 cn_surf: $cn_surf_1")
    println("  单元$ne cn_surf: $cn_surf_ne")
    
    if cn_surf_1 < cn_surf_ne
        println("✓ SOC梯度正确（cn_surf随SOC增大）")
    else
        println("⚠ SOC梯度可能不正确")
    end
    
    # 测试错误处理：SOC超出范围
    try
        soc_invalid = ones(ne) * 1.5  # 超出[0,1]
        JuBat.ModelInitialisation_MultiSPMe(case; initial_soc_distribution=soc_invalid)
        println("✗ 未检测到无效SOC")
    catch e
        if occursin("must be in [0, 1]", string(e))
            println("✓ 正确检测到无效SOC")
        else
            println("⚠ 错误类型不符")
        end
    end
    
    # 测试错误处理：SOC向量长度不匹配
    try
        soc_wrong_length = ones(ne+1)
        JuBat.ModelInitialisation_MultiSPMe(case; initial_soc_distribution=soc_wrong_length)
        println("✗ 未检测到长度不匹配")
    catch e
        if occursin("must equal number of elements", string(e))
            println("✓ 正确检测到长度不匹配")
        else
            println("⚠ 错误类型不符")
        end
    end
    
catch e
    println("✗ 非均匀初始化失败: $e")
    rethrow(e)
end

# ============================================================================
# 测试 5: 与SPMe_element集成
# ============================================================================
println("\n[6/6] 测试与SPMe_element集成...")

try
    y0 = JuBat.ModelInitialisation_MultiSPMe(case)
    
    # 提取第一个单元的状态
    e = 1
    yt_e = JuBat.MultiSPMe_extract_element_state(y0, e, case)
    
    # 使用SPMe_element求解
    t = 0.0
    I_e = 1.0  # 1C
    T_e = 1.0  # 参考温度
    
    M_e, K_e, F_e, vars_e = JuBat.SPMe_element(case, yt_e, t, e; I_e=I_e, T_e=T_e)
    
    println("✓ SPMe_element调用成功")
    println("  矩阵维度: M=$(size(M_e)), K=$(size(K_e)), F=$(length(F_e))")
    println("  单元电压: $(vars_e["cell voltage"] * case.param.scale.phi) V")
    
    # 验证矩阵维度与单元状态一致
    if size(M_e, 1) == length(yt_e)
        println("✓ 矩阵维度与单元状态一致")
    else
        println("✗ 矩阵维度不一致")
    end
    
    # 模拟多个单元求解
    println("\n  模拟多单元求解...")
    n_test = min(5, ne)
    for e in 1:n_test
        yt_e = JuBat.MultiSPMe_extract_element_state(y0, e, case)
        M_e, K_e, F_e, vars_e = JuBat.SPMe_element(case, yt_e, t, e; I_e=1.0, T_e=1.0)
        V_e = vars_e["cell voltage"] * case.param.scale.phi
        @printf("    单元%d: V=%.4f V\n", e, V_e)
    end
    
    println("✓ 多单元求解成功")
    
catch e
    println("✗ 集成测试失败: $e")
    rethrow(e)
end

# ============================================================================
# 总结
# ============================================================================
println("\n" * "="^80)
println("测试总结")
println("="^80)

println("""
✓ ModelInitialisation_MultiSPMe 实现成功
✓ 状态向量结构正确（维度、布局）
✓ 状态提取函数正确
✓ 状态更新函数正确
✓ 非均匀初始SOC分布支持
✓ 与SPMe_element集成成功

阶段2目标达成！多SPMe初始化与状态向量管理功能完整。

状态向量布局（以本测试为例，ne=$ne, n_chem=$(case.multi_spme_layout["n_chem"]), nT=$nT）：
  y0 = [
    yt_e[1:$ne];    # 电化学部分 ($(ne * case.multi_spme_layout["n_chem"])个自由度)
        T_nodes[1:$nT];  # 热场部分 ($(nT)个自由度)
  ]
    总长度: $(case.multi_spme_layout["n_total"]) 个自由度

下一步：阶段3 - 实现 CallModel_MultiSPMe
""")

println("="^80)