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
