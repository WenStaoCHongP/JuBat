# 测试极片模量修复（颗粒 vs 极片 vs 全叠合有效模量）
#
# 运行：julia --project=. example/力学模块验证/test_electrode_coat_modulus.jl

include(joinpath(@__DIR__, "../../src/JuBat.jl"))
using .JuBat
using Printf

println("="^60)
println("TEST 1: Separator/CurrentCollector 结构体字段")
println("="^60)

param_dim = JuBat.ChooseCell("Jellyroll")

@assert param_dim.SP.E   == 1e9   "SP.E 应为 1e9，实际 $(param_dim.SP.E)"
@assert param_dim.SP.nu  == 0.3   "SP.nu 应为 0.3，实际 $(param_dim.SP.nu)"
@assert param_dim.PCC.E  == 70e9  "PCC.E 应为 70e9，实际 $(param_dim.PCC.E)"
@assert param_dim.PCC.nu == 0.3   "PCC.nu 应为 0.3，实际 $(param_dim.PCC.nu)"
@assert param_dim.NCC.E  == 69e9  "NCC.E 应为 69e9，实际 $(param_dim.NCC.E)"
@assert param_dim.NCC.nu == 0.3   "NCC.nu 应为 0.3，实际 $(param_dim.NCC.nu)"

@assert hasproperty(param_dim.SP, :alphaT)  "SP 缺少 alphaT 字段"
@assert hasproperty(param_dim.PCC, :alphaT) "PCC 缺少 alphaT 字段"
@assert hasproperty(param_dim.NCC, :alphaT) "NCC 缺少 alphaT 字段"

println("  PASS: SP/PCC/NCC 结构体字段全部就位")

println("\n" * "="^60)
println("TEST 2: Electrode.E_coat/nu_coat 字段 + Scale.E_coat")
println("="^60)

@assert !hasproperty(param_dim.PE, :E_p)  "PE.E_p 应已重命名为 PE.E_coat"
@assert !hasproperty(param_dim.PE, :nu_p) "PE.nu_p 应已重命名为 PE.nu_coat"
@assert !hasproperty(param_dim.NE, :E_n)  "NE.E_n 应已重命名为 NE.E_coat"
@assert !hasproperty(param_dim.NE, :nu_n) "NE.nu_n 应已重命名为 NE.nu_coat"

param_lgm = JuBat.ChooseCell("LG M50")
@assert param_lgm.PE.E_coat == 0 "LGM50 PE.E_coat 应默认 0"
@assert param_lgm.NE.E_coat == 0 "LGM50 NE.E_coat 应默认 0"

@assert hasproperty(param_dim.scale, :E_coat) "Scale 缺少 E_coat 字段"
@assert hasproperty(param_lgm.scale, :E_coat) "Scale 缺少 E_coat 字段（LGM50）"

println("  PASS: 字段定义符合规格")

println("\n" * "="^60)
println("TEST 3: scale.E_coat 计算正确性 + NaN 防御")
println("="^60)

p = param_dim
expected_E_coat = (
    p.PE.E_coat * p.PE.thickness +
    p.NE.E_coat * p.NE.thickness +
    p.SP.E      * p.SP.thickness  +
    p.PCC.E     * p.PCC.thickness +
    p.NCC.E     * p.NCC.thickness
) / (p.PE.thickness + p.NE.thickness + p.SP.thickness + p.PCC.thickness + p.NCC.thickness)

@assert abs(param_dim.scale.E_coat - expected_E_coat) / expected_E_coat < 1e-12 "scale.E_coat 与手算不一致：$(param_dim.scale.E_coat) vs $expected_E_coat"

@printf("  scale.E_coat = %.3e Pa (预期 %.3e)\n", param_dim.scale.E_coat, expected_E_coat)
@assert 1e9 < param_dim.scale.E_coat < 1e11 "scale.E_coat 量级应在 1e9-1e11 Pa（集流体主导），实际 $(param_dim.scale.E_coat)"

@assert param_lgm.scale.E_coat == 0 "LGM50 scale.E_coat 应为 0（未定义），实际 $(param_lgm.scale.E_coat)"
@assert !isnan(param_lgm.scale.E_coat) "LGM50 scale.E_coat 不能是 NaN"

println("  PASS: scale.E_coat 计算正确，NaN 防御有效")

println("\n" * "="^60)
println("TEST 4: NormaliseParam 归一化一致性")
println("="^60)

opt = JuBat.Option()
opt.model = "SPMe"
opt.per_element_spme = true
case = JuBat.SetCase(param_dim, opt)

expected_PE_E_coat_norm = param_dim.PE.E_coat / param_dim.scale.E_coat
expected_NE_E_coat_norm = param_dim.NE.E_coat / param_dim.scale.E_coat

@assert abs(case.param.PE.E_coat - expected_PE_E_coat_norm) / expected_PE_E_coat_norm < 1e-8 "PE.E_coat 归一化错位：$(case.param.PE.E_coat) vs $expected_PE_E_coat_norm"
@assert abs(case.param.NE.E_coat - expected_NE_E_coat_norm) / expected_NE_E_coat_norm < 1e-8 "NE.E_coat 归一化错位：$(case.param.NE.E_coat) vs $expected_NE_E_coat_norm"
@assert case.param.PE.nu_coat == param_dim.PE.nu_coat "nu_coat 不应被归一化"
@assert case.param.NE.nu_coat == param_dim.NE.nu_coat "nu_coat 不应被归一化"

# SP/PCC/NCC.E 也须按 scale.E_coat 归一化（spec §4.4 增补）
for (comp, label) in ((:SP, "SP"), (:PCC, "PCC"), (:NCC, "NCC"))
    dim_E = getfield(param_dim, comp).E
    norm_E = getfield(case.param, comp).E
    @assert abs(norm_E - dim_E / param_dim.scale.E_coat) / (dim_E / param_dim.scale.E_coat) < 1e-8 "$label.E 归一化错位：$norm_E vs $(dim_E / param_dim.scale.E_coat)"
end

@printf("  PE.E_coat: 物理=%.3e Pa, 归一化=%.3f, 还原=%.3e Pa\n",
        param_dim.PE.E_coat, case.param.PE.E_coat,
        case.param.PE.E_coat * param_dim.scale.E_coat)
@printf("  NE.E_coat: 物理=%.3e Pa, 归一化=%.3f, 还原=%.3e Pa\n",
        param_dim.NE.E_coat, case.param.NE.E_coat,
        case.param.NE.E_coat * param_dim.scale.E_coat)
@printf("  SP/PCC/NCC.E 已按 scale.E_coat 归一化（spec §4.4）\n")

println("  PASS: 归一化一致性（容差 1e-8）")

println("\n" * "="^60)
println("TEST 7: thermal_diffusion_stress_2D 使用极片模量（@assert 防御）")
println("="^60)

mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=20, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)

variables = Dict{String, Union{Array{Float64},Float64}}()
thermal_mesh = case.mesh["thermal2D"]
nnodes = size(thermal_mesh.node, 1)
nelems = size(thermal_mesh.element, 1)
# 注意：函数内部读取归一化 T_nodes 与归一化 param.cell.T0 比较，
# 因此零扰动情景需用归一化温度（case.param.cell.T0），而非物理温度。
variables["T_nodes"] = fill(case.param.cell.T0, nnodes)
variables["thermal2D element soc_n"] = fill(case.param.NE.cs0, nelems)
variables["thermal2D element soc_p"] = fill(case.param.PE.cs0, nelems)

new_vars = JuBat.thermal_diffusion_stress_2D(case, variables)

@assert haskey(new_vars, "diffusion stress vonMises") "thermal_diffusion_stress_2D 未输出 vonMises"
@assert haskey(new_vars, "displacement x") "thermal_diffusion_stress_2D 未输出 displacement"
@printf("  max vonMises = %.3e Pa\n", maximum(new_vars["diffusion stress vonMises"]))
@printf("  max disp_x   = %.3e m\n",  maximum(abs.(new_vars["displacement x"])))

@assert maximum(abs.(new_vars["diffusion stress vonMises"])) < 1e-3 "零扰动情景应力应≈0"

println("  PASS: thermal_diffusion_stress_2D 正常运行，零扰动情景应力≈0")

println("\n" * "="^60)
println("TEST 8: LGM50（无 E_coat）启用力学应被拦截")
println("="^60)

param_lgm = JuBat.ChooseCell("LG M50")
@assert param_lgm.scale.E_coat == 0 "LGM50 scale.E_coat 应保持 0"
@assert !isnan(param_lgm.scale.E_coat) "scale.E_coat 不能是 NaN"

opt_lgm = JuBat.Option()
opt_lgm.model = "SPMe"
case_lgm = JuBat.SetCase(param_lgm, opt_lgm)

variables_lgm = Dict{String, Union{Array{Float64},Float64}}()

try
    JuBat.thermal_diffusion_stress_2D(case_lgm, variables_lgm)
    error("FAIL: 期望抛 AssertionError 但未抛")
catch e
    @assert e isa AssertionError "期望 AssertionError，实际抛 $(typeof(e)): $(e.msg)"
    @assert occursin("E_coat", e.msg) "AssertionError 消息应包含 E_coat，实际：$(e.msg)"
    println("  拦截成功（thermal_diffusion_stress_2D）：$(e.msg)")
end

println("  PASS: LGM50 启用宏观力学被入口 @assert 拦截，未产生 NaN 污染")
println("\n" * "="^60)
println("ALL TESTS PASSED")
println("="^60)
