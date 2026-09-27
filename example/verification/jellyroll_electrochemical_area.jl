#=
Independent geometry audit of the double-sided electrochemical area.

Run with Julia 1.11.2 --startup-file=no --project=. -t1.
Uses the production Jellyroll mesh at nθ=80, 160, 320, without solving a cycle.
Outputs are relative to this file: output/jellyroll_electrochemical_area/.
No geometric area is rescaled to cell.area before comparison.
=#
using Printf, TOML
include(joinpath(@__DIR__, "..", "..", "src", "JuBat.jl"))
using .JuBat

const AREA_OUTDIR = normpath(joinpath(@__DIR__, "..", "..", "output",
                                    "jellyroll_electrochemical_area"))

# Integral of sqrt(r(theta)^2 + b^2) dtheta for r = a + offset + b*theta.
function spiral_length(a, b, offset, theta_start, theta_end)
    primitive(r) = (r * hypot(r, b) + b^2 * asinh(r / b)) / (2b)
    return primitive(a + offset + b * theta_end) -
           primitive(a + offset + b * theta_start)
end

function audit_mesh(pd, param, ntheta)
    data = JuBat.jellyroll_collector_seed_mesh(param; nθ=ntheta,
                                              gsorder=2, czm_enabled=true)
    thermal = data.thermal2D
    czm = JuBat.create_czm_mesh(data.czm_submesh, data.thermal2D_merged, param)
    ne = size(thermal.element, 1)
    L, height, pitch = pd.scale.L, pd.cell.width, pd.cell.layer
    a, b = pd.cell.Rin, pitch / (2pi)
    cross_areas, _ = JuBat.jellyroll_element_properties(thermal, param)
    cross_areas .*= L^2
    @assert all(>(0.0), cross_areas)
    @assert czm.n_cohesive == 4ne

    first_inner = thermal.element[1, 1]
    last_inner = thermal.element[end, 4]
    r_start = hypot(thermal.node[first_inner, 1], thermal.node[first_inner, 2]) * L
    r_end = hypot(thermal.node[last_inner, 1], thermal.node[last_inner, 2]) * L
    theta_start, theta_end = (r_start - a) / b, (r_end - a) / b
    theta_cap = (pd.cell.Rout - a - pitch) / b
    @assert isapprox(theta_end - theta_start, ne * 2pi / ntheta; rtol=1e-12)

    face_names = ("PE_inner", "PE_outer", "NE_inner", "NE_outer")
    offsets = (pd.PE.thickness,
               pd.PE.thickness + pd.PCC.thickness,
               2pd.PE.thickness + pd.PCC.thickness + pd.SP.thickness + pd.NE.thickness,
               2pd.PE.thickness + pd.PCC.thickness + pd.SP.thickness + pd.NE.thickness + pd.NCC.thickness)
    face_areas = zeros(ne, 4)
    face_counts = zeros(Int, ne, 4)
    for (j, elem) in enumerate(czm.cohesive_elements)
        e = czm.cohesive_to_thermal[j]
        inner_material = data.czm_submesh.material_type[elem.host_inner_elem]
        column = if elem.interface_type == :PE_PCC
            inner_material == :PE ? 1 : inner_material == :PCC ? 2 : error("Invalid PE face")
        elseif elem.interface_type == :NE_NCC
            inner_material == :NE ? 3 : inner_material == :NCC ? 4 : error("Invalid NE face")
        else
            error("Unexpected interface $(elem.interface_type)")
        end
        n1, n2 = elem.nodes_bottom
        length_m = hypot(czm.node[n2, 1] - czm.node[n1, 1],
                         czm.node[n2, 2] - czm.node[n1, 2]) * L
        @assert isapprox(length_m, elem.length * L; rtol=1e-12)
        # Both ends have the same radial offset from their parent inner spiral.
        base_r = minimum(hypot(thermal.node[n, 1], thermal.node[n, 2])
                         for n in thermal.element[e, [1, 4]]) * L
        face_r = min(hypot(czm.node[n1, 1], czm.node[n1, 2]),
                     hypot(czm.node[n2, 1], czm.node[n2, 2])) * L
        @assert isapprox(face_r - base_r, offsets[column]; atol=1e-12, rtol=1e-10)
        face_areas[e, column] += length_m * height
        face_counts[e, column] += 1
    end
    @assert all(==(1), face_counts) "Each thermal segment must own four distinct faces"
    totals = vec(sum(face_areas; dims=1))
    exact = [height * spiral_length(a, b, offset, theta_start, theta_end) for offset in offsets]
    exact_unclipped = [height * spiral_length(a, b, offset, 0.0, theta_cap) for offset in offsets]
    @assert all(totals .< exact) "Straight chords must be shorter than the spiral arcs"

    midline_length = 0.0
    for e in 1:ne
        n1, n2, n3, n4 = thermal.element[e, :]
        dx = (thermal.node[n3, 1] + thermal.node[n4, 1] -
              thermal.node[n1, 1] - thermal.node[n2, 1]) / 2
        dy = (thermal.node[n3, 2] + thermal.node[n4, 2] -
              thermal.node[n1, 2] - thermal.node[n2, 2]) / 2
        midline_length += hypot(dx, dy) * L
    end

    open(joinpath(AREA_OUTDIR, "element_areas_n$(ntheta).csv"), "w") do io
        println(io, "thermal_element,cross_section_m2,PE_inner_m2,PE_outer_m2,NE_inner_m2,NE_outer_m2")
        for e in 1:ne
            println(io, join((e, cross_areas[e], face_areas[e, :]...), ","))
        end
    end
    return (; ntheta, ne, theta_start, theta_end, theta_cap, face_names, totals, exact,
            exact_unclipped, cross_area=sum(cross_areas), midline_length,
            midline_exact=spiral_length(a, b, pitch / 2, theta_start, theta_end),
            max_chord_error=maximum(abs.(totals ./ exact .- 1)))
end

function main()
    mkpath(AREA_OUTDIR)
    pd = JuBat.ChooseCell("Jellyroll")
    param = JuBat.NormaliseParam(pd)
    reference = pd.cell.width * pd.cell.length * pd.cell.no_layers
    @assert reference == pd.cell.area
    results = [audit_mesh(pd, param, ntheta) for ntheta in (80, 160, 320)]
    @assert results[3].max_chord_error < results[2].max_chord_error < results[1].max_chord_error
    base = first(results)
    rows = NamedTuple[]
    addrow(label, area) = push!(rows, (metric=label, area_m2=area,
                                      relative_difference_percent=100 * (area / reference - 1)))
    addrow("cell.area reference", reference)
    addrow("2 x thermal cross-section (NOT unfolded area)", 2base.cross_area)
    addrow("2 x width x midline chord length", 2pd.cell.width * base.midline_length)
    addrow("2 x width x midline exact arc length", 2pd.cell.width * base.midline_exact)
    addrow("PE two faces: mesh", sum(base.totals[1:2]))
    addrow("NE two faces: mesh", sum(base.totals[3:4]))
    addrow("PE two faces: exact same angular extent", sum(base.exact[1:2]))
    addrow("NE two faces: exact same angular extent", sum(base.exact[3:4]))
    addrow("PE two faces: before angular grid clipping", sum(base.exact_unclipped[1:2]))
    addrow("NE two faces: before angular grid clipping", sum(base.exact_unclipped[3:4]))
    addrow("2 x width x mesh cross-section / repeat thickness", 2pd.cell.width * base.cross_area / pd.cell.layer)
    addrow("full-annulus volume proxy: 2 x volume / repeat thickness", 2pd.cell.volume / pd.cell.layer)

    open(joinpath(AREA_OUTDIR, "area_comparison.csv"), "w") do io
        println(io, "metric,area_m2,relative_difference_percent")
        for row in rows
            println(io, join((row.metric, row.area_m2, row.relative_difference_percent), ","))
        end
    end
    open(joinpath(AREA_OUTDIR, "metrics.toml"), "w") do io
        TOML.print(io, Dict(
            "parameters" => Dict("cell_area_m2"=>reference, "width_m"=>pd.cell.width,
                "length_m"=>pd.cell.length, "no_layers"=>pd.cell.no_layers,
                "Rin_m"=>pd.cell.Rin, "Rout_m"=>pd.cell.Rout,
                "repeat_thickness_m"=>pd.cell.layer, "normalization_length_m"=>pd.scale.L),
            "comparison" => [Dict(string(k)=>v for (k,v) in pairs(row)) for row in rows],
            "mesh_checks" => [Dict("ntheta"=>r.ntheta, "thermal_elements"=>r.ne,
                "theta_start_rad"=>r.theta_start, "theta_end_rad"=>r.theta_end,
                "theta_cap_rad"=>r.theta_cap, "max_chord_error_percent"=>100r.max_chord_error,
                "PE_area_m2"=>sum(r.totals[1:2]), "NE_area_m2"=>sum(r.totals[3:4])) for r in results]))
    end

    report = IOBuffer()
    println(report, "# Jellyroll 电化学面积独立核查\n")
    println(report, "使用当前参数和真实生产网格；不运行电化学/力学求解，不按 cell.area 校准几何面积。\n")
    @printf(report, "cell.area = %.6g × %.6g × %.6g = **%.9f m²**。\n\n",
            pd.cell.width, pd.cell.length, pd.cell.no_layers, reference)
    @printf(report, "Rin=%.6f mm，Rout=%.6f mm，完整重复层厚度=%.3f μm，width=%.6f m。\n\n",
            pd.cell.Rin*1e3, pd.cell.Rout*1e3, pd.cell.layer*1e6, pd.cell.width)
    @printf(report, "主核查 nθ=80：%d 个热单元，%d 个真实 cohesive 面；θ 范围 %.9f–%.9f rad（%.9f 匝）。\n\n",
            base.ne, 4base.ne, base.theta_start, base.theta_end, (base.theta_end-base.theta_start)/(2pi))
    println(report, "| 面积定义 | 面积 [m²] | 相对 cell.area 偏差 |\n|---|---:|---:|")
    for row in rows
        @printf(report, "| %s | %.9f | %+.6f%% |\n", row.metric, row.area_m2, row.relative_difference_percent)
    end
    println(report, "\n## 四个真实面的贡献\n\n| 面 | 单面面积 [m²] |\n|---|---:|")
    for i in 1:4
        @printf(report, "| %s | %.9f |\n", base.face_names[i], base.totals[i])
    end
    println(report, "\n## 离散误差复核\n\n| nθ | 热单元数 | PE 双面 [m²] | NE 双面 [m²] | 面弦长对同域解析弧长最大误差 |\n|---:|---:|---:|---:|---:|")
    for r in results
        @printf(report, "| %d | %d | %.9f | %.9f | %.6f%% |\n",
                r.ntheta, r.ne, sum(r.totals[1:2]), sum(r.totals[3:4]), 100r.max_chord_error)
    end
    println(report, "\n## 解释与验证边界\n")
    println(report, "- 热横截面积 Q4 积分乘 L²；真实面展开面积为 cohesive 边长乘 L，再乘电池高度 width。每极分别加两面，不能正负极相加成四倍面积。")
    println(report, "- 2×热横截面积并非双面展开面积。中线近似为 2×width×卷绕中线长度；横截面积代理须额外乘 width/完整重复层厚度。")
    println(report, "- PE 与 NE 位于不同径向偏置，所以实际长度和双面面积略有差别；同极两面也非严格等面积。面积损失宜保留逐面 A_j 权重，f=Σ[A_j(1-D_j)]/ΣA_j。")
    println(report, "- 解析弧长与实际网格使用相同 θ 范围；另列未作角节点裁剪的结果。完整圆环体积代理不等于当前螺旋网格实际覆盖域。")
    println(report, "- 检查包括四面拓扑与归属、径向偏置、坐标与长度单位、正面积、弦长小于弧长及加密后误差降低。相对 cell.area 的偏差如实报告，不以预设容差强行判定标定通过。")
    text = String(take!(report))
    write(joinpath(AREA_OUTDIR, "report.md"), text)
    print(text)
    println("\nOUTPUT: ", AREA_OUTDIR)
end

main()
