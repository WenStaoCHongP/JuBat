# 多 SPMe 输出按热单元顺序保存；场量在历史矩阵中按单元连续排列行。
# PostProcessing 将场量恢复为（内部自由度，单元，时间），标量为（单元，时间）。
const SPME_ELEMENT_CHEMICAL_OUTPUT_KEYS = (
    "negative particle lithium concentration",
    "positive particle lithium concentration",
    "negative particle surface lithium concentration",
    "positive particle surface lithium concentration",
    "electrolyte lithium concentration",
    "electrolyte lithium concentration in negative electrode",
    "electrolyte lithium concentration in positive electrode",
    "electrolyte lithium concentration in separator",
    "electrolyte lithium concentration at negative electrode Gauss point",
    "electrolyte lithium concentration at positive electrode Gauss point",
    "electrolyte lithium concentration at separator Gauss point",
    "negative electrode exchange current density",
    "positive electrode exchange current density",
    "negative electrode interfacial current density",
    "positive electrode interfacial current density",
    "negative electrode overpotential",
    "positive electrode overpotential",
    "negative electrode open circuit potential",
    "positive electrode open circuit potential",
)

const SPME_ELEMENT_MECHANICAL_OUTPUT_KEYS = (
    "negative particle center radial stress",
    "positive particle center radial stress",
    "negative particle surface tangential stress",
    "positive particle surface tangential stress",
    "negative particle surface displacement",
    "positive particle surface displacement",
    "negative particle stress coupling diffusion coefficient",
    "positive particle stress coupling diffusion coefficient",
    "negative particle concentration at gauss point",
    "positive particle concentration at gauss point",
)

function spme_element_output_keys(case::Case)
    if case.opt.mechanicalmodel == "full"
        return (SPME_ELEMENT_CHEMICAL_OUTPUT_KEYS..., SPME_ELEMENT_MECHANICAL_OUTPUT_KEYS...)
    end
    return SPME_ELEMENT_CHEMICAL_OUTPUT_KEYS
end

function StandardVariables(case::Case, num::Int64)
    n_spme = 1
    if case.opt.model == "SPMe" && case.opt.thermalmodel == "distributed2D"
        n_spme = size(case.mesh["thermal2D"].element, 1)
    end
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        Nn = 1
        Np = 1
        Ne_ngs = case.opt.Nn * case.opt.gsorder
        Ne_pgs = case.opt.Np * case.opt.gsorder
    elseif case.opt.model == "P2D"
        Nn = case.mesh["negative electrode"].nlen
        Np = case.mesh["positive electrode"].nlen
        Ne_ngs = case.opt.Nn * case.opt.gsorder
        Ne_pgs = case.opt.Np * case.opt.gsorder
    end

    variables = Dict{String, Union{Array{Float64}, Float64}}(
        "negative particle lithium concentration" => zeros(Float64, Nrn * n_spme, num),
        "positive particle lithium concentration" => zeros(Float64, Nrp * n_spme, num),
        "negative particle averaged lithium concentration" => zeros(Float64, Nn, num),
        "positive particle averaged lithium concentration" => zeros(Float64, Np, num),
        "negative particle surface lithium concentration" => zeros(Float64, Nn * n_spme, num),
        "positive particle surface lithium concentration" => zeros(Float64, Np * n_spme, num),
        "negative electrode porosity" => zeros(Float64, Nn, num),
        "positive electrode porosity" => zeros(Float64, Np, num),
        "negative electrode temperature" => zeros(Float64, Nn, num),
        "positive electrode temperature" => zeros(Float64, Np, num),
        "negative electrode exchange current density" => zeros(Float64, Nn * n_spme, num),
        "positive electrode exchange current density" => zeros(Float64, Np * n_spme, num),
        "negative electrode interfacial current density" => zeros(Float64, Nn * n_spme, num),
        "positive electrode interfacial current density" => zeros(Float64, Np * n_spme, num),
        "negative electrode overpotential" => zeros(Float64, Nn * n_spme, num),
        "positive electrode overpotential" => zeros(Float64, Np * n_spme, num),
        "negative electrode open circuit potential" => zeros(Float64, Nn * n_spme, num),
        "positive electrode open circuit potential" => zeros(Float64, Np * n_spme, num),
        "negative particle center radial stress" => zeros(Float64, Nn * n_spme, num),
        "positive particle center radial stress" => zeros(Float64, Np * n_spme, num),
        "negative particle surface tangential stress" => zeros(Float64, Nn * n_spme, num),
        "positive particle surface tangential stress" => zeros(Float64, Np * n_spme, num),
        "negative particle surface displacement" => zeros(Float64, Nn * n_spme, num),
        "positive particle surface displacement" => zeros(Float64, Np * n_spme, num),
        "negative particle concentration at gauss point" => zeros(Float64, Nn * case.opt.Nrn * case.opt.gsorder * n_spme, num),
        "positive particle concentration at gauss point" => zeros(Float64, Np * case.opt.Nrp * case.opt.gsorder * n_spme, num),
        "negative particle surface tangential stress at gauss point" => zeros(Float64, Ne_ngs, num),
        "positive particle surface tangential stress at gauss point" => zeros(Float64, Ne_pgs, num),
        "negative particle stress coupling diffusion coefficient" => zeros(Float64, Nn * n_spme, num),
        "positive particle stress coupling diffusion coefficient" => zeros(Float64, Np * n_spme, num),
        "cell voltage" => zeros(Float64, 1, num),
        "time" => zeros(Float64, 1, num),
        "cell current" => zeros(Float64, 1, num),
    )
    # additional variables for SPMe and P2D
    if case.opt.model == "SPMe" || case.opt.model == "P2D"
        Ne = case.mesh["electrolyte"].nlen
        Ne_n = case.mesh["negative electrode"].nlen
        Ne_p = case.mesh["positive electrode"].nlen
        Ne_sp = case.mesh["separator"].nlen
        Ne_ngs = case.opt.Nn * case.opt.gsorder
        Ne_pgs = case.opt.Np * case.opt.gsorder
        Ne_spgs = case.opt.Ns * case.opt.gsorder
        variables["electrolyte lithium concentration"] = zeros(Float64, Ne * n_spme, num)
        variables["electrolyte lithium concentration in negative electrode"] = zeros(Float64, Ne_n * n_spme, num)
        variables["electrolyte lithium concentration in positive electrode"] = zeros(Float64, Ne_p * n_spme, num)
        variables["electrolyte lithium concentration in separator"] = zeros(Float64, Ne_sp * n_spme, num)
        variables["electrolyte lithium concentration at negative electrode Gauss point"] = zeros(Float64, Ne_ngs * n_spme, num)
        variables["electrolyte lithium concentration at positive electrode Gauss point"] = zeros(Float64, Ne_pgs * n_spme, num)
        variables["electrolyte lithium concentration at separator Gauss point"] = zeros(Float64, Ne_spgs * n_spme, num)
    end
    # extra variables for P2D
    if case.opt.model == "P2D"
        variables["negative electrode potential"] = zeros(Float64, Nn, num)
        variables["positive electrode potential"] = zeros(Float64, Np, num)
        variables["electrolyte potential"] = zeros(Float64, Ne, num)
        variables["electrolyte potential in negative electrode"] = zeros(Float64, Ne_n, num)
        variables["electrolyte potential in positive electrode"] = zeros(Float64, Ne_p, num)
        variables["electrolyte potential in separator"] = zeros(Float64, Ne_sp, num)
        variables["negative electrode interfacial current at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode interfacial current at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode open circuit potential at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode open circuit potential at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode overpotential at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode overpotential at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative electrode exchange current density at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive electrode exchange current density at Gauss point"] = zeros(Float64, Ne_pgs, num)
        variables["negative particle surface lithium concentration at Gauss point"] = zeros(Float64, Ne_ngs, num)
        variables["positive particle surface lithium concentration at Gauss point"] = zeros(Float64, Ne_pgs, num)
    end

    variables["temperature"] = zeros(Float64, 1, num)

    if case.opt.thermalmodel == "lumped"
        variables["thermal lumped internal heat"] = zeros(Float64, 1, num)
    end

    if case.opt.thermalmodel == "distributed2D"
        ne = size(case.mesh["thermal2D"].element, 1)
        nT = case.mesh["thermal2D"].nlen
        variables["thermal2D temperature at nodes"] = zeros(Float64, nT, num)
        variables["heat_source_fields"] = zeros(Float64, ne, num)
        variables["thermal2D q_rxn_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_rev_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_s_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_e_ne"] = zeros(Float64, ne, num)
        variables["thermal2D q_sp"] = zeros(Float64, ne, num)
        variables["thermal2D q_rxn_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_rev_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_s_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_ohm_e_pe"] = zeros(Float64, ne, num)
        variables["thermal2D q_pcc"] = zeros(Float64, ne, num)
        variables["thermal2D q_ncc"] = zeros(Float64, ne, num)
        variables["thermal2D q_el"] = zeros(Float64, ne, num)
        variables["thermal2D element current"] = zeros(Float64, ne, num)
        if case.opt.per_element_spme
            variables["thermal2D branch area factor"] = zeros(Float64, ne, num)
            variables["thermal2D total effective area factor"] = zeros(Float64, 1, num)
            variables["thermal2D element branch current [A]"] = zeros(Float64, ne, num)
        end
        variables["thermal2D eta_n_e"] = zeros(Float64, ne, num)
        variables["thermal2D eta_p_e"] = zeros(Float64, ne, num)
        variables["thermal2D dUdT_n_e"] = zeros(Float64, ne, num)
        variables["thermal2D dUdT_p_e"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_n"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_n_vol"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_p_vol"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_n_nodeavg"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_p_nodeavg"] = zeros(Float64, ne, num)
        variables["thermal2D element soc_p"] = zeros(Float64, ne, num)
        variables["thermal2D element voltages"] = zeros(Float64, ne, num)
        variables["thermal2D element OCV"] = zeros(Float64, ne, num)
        variables["thermal2D active_mask"] = zeros(Float64, ne, num)
        variables["thermal2D deactivated_mask"] = zeros(Float64, ne, num)
        variables["thermal2D inactive_reason"] = zeros(Float64, ne, num)
        variables["thermal2D n_active_elements"] = zeros(Float64, 1, num)
        variables["thermal2D n_deactivated_elements"] = zeros(Float64, 1, num)
        variables["thermal2D Vsolve status"] = zeros(Float64, 1, num)
        variables["thermal2D Vsolve iters"] = zeros(Float64, 1, num)
        variables["thermal2D Vsolve converged"] = zeros(Float64, 1, num)
        variables["thermal2D n_cutoff_elements"] = zeros(Float64, 1, num)
        variables["thermal2D nearest_cutoff_element"] = zeros(Float64, 1, num)
        variables["thermal2D nearest_cutoff_ocv"] = zeros(Float64, 1, num)
        variables["thermal2D margin_to_cutoff"] = zeros(Float64, 1, num)
        variables["total heat source"] = zeros(Float64, 1, num)
    end
    if case.opt.czm.enabled == true
        variables["negative electrode cohesive zone damage"] = zeros(Float64, Nn, num)
        variables["positive electrode cohesive zone damage"] = zeros(Float64, Np, num)
        variables["czm D_max"] = zeros(Float64, 1, num)
        variables["czm D_mean"] = zeros(Float64, 1, num)
        variables["czm δ_max_n"] = zeros(Float64, 1, num)
        variables["czm δ_max_eff"] = zeros(Float64, 1, num)
        variables["czm δ_mean_n"] = zeros(Float64, 1, num)
        variables["czm n_fractured"] = zeros(Float64, 1, num)
        if case.opt.czm.winding_prestress
            variables["winding prestress"] = zeros(Float64, 1, num)
        end
        if haskey(case.mesh, "czm")
            n_coh = length(case.mesh["czm"].cohesive_elements)
            n_czm_nodes = size(case.mesh["czm"].node, 1)
            variables["czm damage"] = zeros(Float64, n_coh, num)
            variables["czm displacement x"] = zeros(Float64, n_czm_nodes, num)
            variables["czm displacement y"] = zeros(Float64, n_czm_nodes, num)
            variables["czm traction normal"] = zeros(Float64, n_coh, num)
            variables["czm traction tangent"] = zeros(Float64, n_coh, num)
            variables["czm separation normal"] = zeros(Float64, n_coh, num)
            variables["czm separation tangent"] = zeros(Float64, n_coh, num)
            variables["czm separation effective"] = zeros(Float64, n_coh, num)
            ne_czm = size(case.mesh["czm"].bulk_element, 1)
            if case.opt.czm.j2_plasticity
                variables["diffusion stress xx"] = fill(NaN, ne_czm, num)
                variables["diffusion stress yy"] = fill(NaN, ne_czm, num)
                variables["diffusion stress xy"] = fill(NaN, ne_czm, num)
                variables["diffusion stress vonMises"] = fill(NaN, ne_czm, num)
                variables["diffusion stress max vonMises"] = fill(NaN, 1, num)
                variables["equivalent plastic strain max"] = fill(NaN, 1, num)
            elseif !(case.opt.czm.geo_nonlinear || case.opt.czm.winding_prestress)
                variables["diffusion stress xx"] = fill(NaN, ne_czm, num)
                variables["diffusion stress yy"] = fill(NaN, ne_czm, num)
                variables["diffusion stress xy"] = fill(NaN, ne_czm, num)
                variables["diffusion stress vonMises"] = fill(NaN, ne_czm, num)
            else
                @warn "层分辨应力恢复与无 J2 状态的几何非线性/预应力路径不相容，本运行不导出 diffusion stress 历史" maxlog=1
            end
        end
    end
    return variables
end

"""
    create_element_workspace(case)

创建精简型单元工作区 Dict，仅包含 SPMe 分布式单元调用链实际需要的键。
排除 distributed2D 专用的 ~30 个 thermal2D 键（由 CallModel_MultiSPMe
在单元循环外独立管理）。比 StandardVariables(case, 1) 减少约 60% 数组分配。
"""
function create_element_workspace(case::Case)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    Nn = 1; Np = 1  # SPMe 模式
    Ne_ngs = case.opt.Nn * case.opt.gsorder
    Ne_pgs = case.opt.Np * case.opt.gsorder

    ws = Dict{String, Union{Array{Float64}, Float64}}()

    # ── 状态提取键（case.index 对应，SPMe_variables 的 workspace 模式从 yt 原位覆写）──
    ws["negative particle lithium concentration"] = zeros(Float64, Nrn, 1)
    ws["positive particle lithium concentration"] = zeros(Float64, Nrp, 1)
    ws["negative particle surface lithium concentration"] = zeros(Float64, Nn, 1)
    ws["positive particle surface lithium concentration"] = zeros(Float64, Np, 1)

    # ── SPMe 专用键 ──
    if case.opt.model == "SPMe"
        Ne_n = case.mesh["negative electrode"].nlen
        Ne_p = case.mesh["positive electrode"].nlen
        Ne_sp = case.mesh["separator"].nlen
        Ne_spgs = case.opt.Ns * case.opt.gsorder
        ws["electrolyte lithium concentration in negative electrode"] = zeros(Float64, Ne_n, 1)
        ws["electrolyte lithium concentration in positive electrode"] = zeros(Float64, Ne_p, 1)
        ws["electrolyte lithium concentration in separator"] = zeros(Float64, Ne_sp, 1)
        # Gauss 点计算结果（ElectrolyteDiffusion 读取）
        ws["electrolyte lithium concentration at negative electrode Gauss point"] = zeros(Float64, Ne_ngs, 1)
        ws["electrolyte lithium concentration at positive electrode Gauss point"] = zeros(Float64, Ne_pgs, 1)
        ws["electrolyte lithium concentration at separator Gauss point"] = zeros(Float64, Ne_spgs, 1)
    end

    ws["temperature"] = 0.0

    # ── SPMe_variables 的 workspace 模式计算结果键 ──
    ws["cell voltage"] = 0.0
    ws["time"] = 0.0
    ws["cell current"] = 0.0
    ws["negative electrode exchange current density"] = zeros(Float64, Nn, 1)
    ws["positive electrode exchange current density"] = zeros(Float64, Np, 1)
    ws["negative electrode interfacial current density"] = zeros(Float64, Nn, 1)
    ws["positive electrode interfacial current density"] = zeros(Float64, Np, 1)
    ws["negative electrode overpotential"] = zeros(Float64, Nn, 1)
    ws["positive electrode overpotential"] = zeros(Float64, Np, 1)
    ws["negative electrode open circuit potential"] = zeros(Float64, Nn, 1)
    ws["positive electrode open circuit potential"] = zeros(Float64, Np, 1)

    # ── 高斯点浓度键（SPMe 装配无条件读取）──
    # 注意：用 opt.Nrn/Nrp（单元数）而非 mesh.nlen（节点数）来计算高斯点数
    ws["negative particle concentration at gauss point"] = zeros(Float64, Nn * case.opt.Nrn * case.opt.gsorder, 1)
    ws["positive particle concentration at gauss point"] = zeros(Float64, Np * case.opt.Nrp * case.opt.gsorder, 1)

    # ── Mechanicaloutput 结果键（条件）──
    if case.opt.mechanicalmodel == "full"
        ws["negative particle center radial stress"] = zeros(Float64, Nn, 1)
        ws["positive particle center radial stress"] = zeros(Float64, Np, 1)
        ws["negative particle surface tangential stress"] = zeros(Float64, Nn, 1)
        ws["positive particle surface tangential stress"] = zeros(Float64, Np, 1)
        ws["negative particle surface displacement"] = zeros(Float64, Nn, 1)
        ws["positive particle surface displacement"] = zeros(Float64, Np, 1)
        ws["negative particle surface tangential stress at gauss point"] = zeros(Float64, Ne_ngs, 1)
        ws["positive particle surface tangential stress at gauss point"] = zeros(Float64, Ne_pgs, 1)
        ws["negative particle stress coupling diffusion coefficient"] = zeros(Float64, Nn, 1)
        ws["positive particle stress coupling diffusion coefficient"] = zeros(Float64, Np, 1)
    end

    # ── CZM 键（条件）──
    if case.opt.czm.enabled
        ws["negative electrode cohesive zone damage"] = zeros(Float64, Nn, 1)
        ws["positive electrode cohesive zone damage"] = zeros(Float64, Np, 1)
        ws["czm D_max"] = 0.0
        ws["czm D_mean"] = 0.0
        ws["czm n_fractured"] = 0.0
        if haskey(case.mesh, "czm")
            ws["czm damage"] = zeros(Float64, length(case.mesh["czm"].cohesive_elements), 1)
        end
    end

    return ws
end

function Variable_update!(variables_hist::Dict{String, Union{Array{Float64},Float64}}, variables::Dict{String, Union{Array{Float64},Float64}}, v::Int64; debug_coupling::Bool=false)
    v >= 1 || throw(BoundsError(variables_hist,v))
    # Resolve a complete column and all growth before touching any existing history.
    updates = Tuple{String,Union{Matrix{Float64},Float64},Union{Vector{Float64},Float64}}[]
    for (k,hist_val) in variables_hist
        val = variables[k]
        if hist_val isa Array{Float64}
            ndims(hist_val) == 2 || throw(DimensionMismatch("history variable \"$k\" must be a matrix"))
            nrows,ncols = size(hist_val)
            column = if val isa Array{Float64}
                length(val) == nrows || throw(DimensionMismatch(
                    "variable \"$k\" length $(length(val)) does not match history row count $nrows"))
                copy(vec(val))
            else
                nrows == 1 || throw(DimensionMismatch(
                    "scalar variable \"$k\" cannot populate history with $nrows rows"))
                [val]
            end
            target = hist_val
            if v > ncols
                expansion = max(1000,ncols÷2)
                v <= ncols+expansion || throw(BoundsError(hist_val,(:,v)))
                @warn "时间步 $(v) 超过预分配 $(ncols)，扩展 $(expansion) 步（变量: $(k)）"
                target = hcat(hist_val,zeros(nrows,expansion))
            end
            push!(updates,(k,target,column))
        else
            val isa Float64 || throw(ArgumentError(
                "scalar history variable \"$k\" requires a Float64 value, got $(typeof(val))"))
            push!(updates,(k,hist_val,val))
        end
    end
    # Only checked scalar assignments remain; no callbacks or allocations in publication.
    for (k,target,column) in updates
        if target isa Matrix{Float64}
            for i in eachindex(column)
                target[i,v] = column[i]
            end
            variables_hist[k] = target
        else
            variables_hist[k] = column
        end
    end
    return variables_hist
end
