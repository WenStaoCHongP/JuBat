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

@printf("  PE.E_coat: 物理=%.3e Pa, 归一化=%.3f, 还原=%.3e Pa\n",
        param_dim.PE.E_coat, case.param.PE.E_coat,
        case.param.PE.E_coat * param_dim.scale.E_coat)
@printf("  NE.E_coat: 物理=%.3e Pa, 归一化=%.3f, 还原=%.3e Pa\n",
        param_dim.NE.E_coat, case.param.NE.E_coat,
        case.param.NE.E_coat * param_dim.scale.E_coat)

println("  PASS: 归一化一致性（容差 1e-8）")

println("\n" * "="^60)
println("TEST 5: compute_effective_coating_modulus 全叠合加权")
println("="^60)

E_eff, nu_eff, alpha_eff = JuBat.compute_effective_coating_modulus(case)

@assert abs(E_eff - 1.0) < 1e-8 "归一化 E_eff 应≈1.0（同尺度），实际 $E_eff"

expected_nu = (
    case.param.PE.nu_coat * case.param.PE.thickness +
    case.param.NE.nu_coat * case.param.NE.thickness +
    case.param.SP.nu      * case.param.SP.thickness  +
    case.param.PCC.nu     * case.param.PCC.thickness +
    case.param.NCC.nu     * case.param.NCC.thickness
) / (case.param.PE.thickness + case.param.NE.thickness +
    case.param.SP.thickness  + case.param.PCC.thickness + case.param.NCC.thickness)
@assert abs(nu_eff - expected_nu) < 1e-12 "nu_eff 与手算不符：$nu_eff vs $expected_nu"

@assert alpha_eff >= 0 "alpha_eff 应非负"

@printf("  E_eff (归一化) = %.6f\n", E_eff)
@printf("  nu_eff         = %.6f (预期 %.6f)\n", nu_eff, expected_nu)
@printf("  alpha_eff      = %.6e\n", alpha_eff)
println("  PASS: compute_effective_coating_modulus 输出正确")

println("\n" * "="^60)
println("TEST 6: compute_czm_effective_params 与共享函数尺度换算一致")
println("="^60)

E_eff_czm, nu_eff_czm, alpha_eff_czm, beta_n, beta_p = JuBat.compute_czm_effective_params(case)

E_eff_coat_norm, _, _ = JuBat.compute_effective_coating_modulus(case)
expected_E_eff_czm = E_eff_coat_norm * param_dim.scale.E_coat / param_dim.scale.σ_czm

@assert abs(E_eff_czm - expected_E_eff_czm) / expected_E_eff_czm < 1e-12 "CZM E_eff 与共享函数尺度换算不一致：$E_eff_czm vs $expected_E_eff_czm"

_, expected_nu_czm, _ = JuBat.compute_effective_coating_modulus(case)
@assert nu_eff_czm == expected_nu_czm "nu_eff 应与共享函数完全一致"

@assert beta_n == case.param.NE.Omega / 3.0
@assert beta_p == case.param.PE.Omega / 3.0

@printf("  E_eff (σ_czm 归一化) = %.6e\n", E_eff_czm)
@printf("  E_eff (E_coat 归一化→σ_czm) = %.6e\n", expected_E_eff_czm)
println("  PASS: compute_czm_effective_params 与共享函数一致")
