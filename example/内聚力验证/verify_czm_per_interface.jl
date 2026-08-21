include(joinpath(@__DIR__, "../../src/JuBat.jl"))
using .JuBat
using Printf

"""
单元级验证：单 PE-PCC cohesive 单元单轴拉开
按 spec v2 §8.1，断言与 param_cache 自洽（不硬编码占位数值）

注：CohesiveElement 实际有 8 字段（id, nodes, nodes_bottom, nodes_top,
    length, interface_type, host_outer_elem, host_inner_elem）—— 单元级
    验证只需前 6 个语义字段，后两个 host_elem 在全网格装配时才用到，
    此处置为 0。
"""
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
case = JuBat.SetCase(param_dim, opt)
param_cache = JuBat.compute_czm_params_per_interface(case)
pe = param_cache.by_interface[:PE_PCC]

# 构造单 cohesive 单元（手工 4 节点，无 mesh）
coh = JuBat.CohesiveElement(
    1,                       # id
    [1, 2, 3, 4],            # nodes [n1, n2, n3, n4]
    [1, 2],                  # nodes_bottom
    [4, 3],                  # nodes_top
    1.0,                     # length
    :PE_PCC,                 # interface_type
    0,                       # host_outer_elem（单元级测试未使用）
    0                        # host_inner_elem（单元级测试未使用）
)
damage = JuBat.DamageState()

# 位移控制加载：δ_n 从 0 到 1.5*δ_c_n
# δ_0_n/δ_c_n 极小，必须显式注入 δ_0_n 附近样本
n_steps = 100
δ_n_history = sort(unique(vcat(
    collect(range(0, 1.5 * pe.δ_c_n; length=n_steps)),
    [0.5 * pe.δ_0_n, pe.δ_0_n, 2.0 * pe.δ_0_n]
)))
n_samples = length(δ_n_history)
σ_n_history = zeros(n_samples)
D_history = zeros(n_samples)

for (i, δ_n) in enumerate(δ_n_history)
    T_n, _, D_i = JuBat.bilinear_traction(δ_n, 0.0, damage, pe; update=true)
    σ_n_history[i] = T_n
    D_history[i] = D_i
end

# 打印概要
println("=" ^ 60)
println("PE-PCC 单元级验证")
println("=" ^ 60)
println("σ_max (param_cache) : $(pe.σ_max)")
println("δ_0_n (param_cache) : $(pe.δ_0_n)")
println("δ_c_n (param_cache) : $(pe.δ_c_n)")
println("G_c   (param_cache) : $(pe.G_c)")
println()
println("峰值 σ_n (实测)     : $(maximum(σ_n_history))")
println("δ @ 峰值            : $(δ_n_history[argmax(σ_n_history)])")
println("D 终值              : $(D_history[end])")
println("=" ^ 60)

# 与 param_cache 自洽的断言（spec §8.1，不硬编码）
@assert maximum(σ_n_history) ≈ pe.σ_max rtol=1e-6 "峰值牵引应等于 σ_max"
@assert argmax(σ_n_history) <= cld(n_samples, 2) + 1 "峰值应在 δ_0_n 附近（前半段）"
@assert D_history[end] ≈ 1.0 atol=1e-6 "完全断裂时 D=1"
@assert D_history[1] ≈ 0.0 "未加载时 D=0"

# 健全性检查
@assert all(isfinite, σ_n_history) "σ_n 应有限"
@assert all(isfinite, D_history) "D 应有限"
@assert all(0 .<= D_history .<= 1.0) "D 应在 [0, 1]"

# 损伤单调
nonzero_D = findall(D_history .> 0)
if length(nonzero_D) > 1
    @assert all(diff(D_history[nonzero_D]) .>= -1e-10) "D 应单调非减"
end

println("✓ 单元级验证通过")

# 可选：保存 σ-δ 曲线 CSV
if haskey(ENV, "SAVE_CSV")
    using DelimitedFiles
    out_dir = joinpath(@__DIR__, "..", "..", "output", "verify_czm_per_interface")
    mkpath(out_dir)
    writedlm(joinpath(out_dir, "czm_pe_pcc_sigma_delta.csv"),
             [δ_n_history σ_n_history D_history], ',')
    println("σ-δ 曲线已保存到 ", joinpath(out_dir, "czm_pe_pcc_sigma_delta.csv"))
end
