using LinearAlgebra, Statistics, Plots, SparseArrays
using Printf
import PyPlot

if !isdefined(Main, :JuBat)
    include(joinpath(@__DIR__, "..", "..", "src", "JuBat.jl"))
end

function compute_q_elem(mesh, q_func, t)
    ne = size(mesh.element, 1)
    q_elem = zeros(Float64, ne)
    w_elem = zeros(Float64, ne)

    wJ = mesh.gs.weight .* mesh.gs.detJ
    xg = mesh.gs.x[:, 1]
    yg = mesh.gs.x[:, 2]
    rg = hypot.(xg, yg)
    thetag = atan.(yg, xg)
    ele = mesh.gs.ele

    @inbounds for g in 1:length(wJ)
        e = ele[g]
        qg = q_func(rg[g], thetag[g], t)
        q_elem[e] += qg * wJ[g]
        w_elem[e] += wJ[g]
    end

    @inbounds for e in 1:ne
        w_elem[e] > 0.0 && (q_elem[e] /= w_elem[e])
    end

    return q_elem
end

function radial_profile(mesh, T)
    r = hypot.(mesh.node[:, 1], mesh.node[:, 2])
    r_unique = unique(sort(r))
    T_avg = zeros(Float64, length(r_unique))
    for (i, rv) in enumerate(r_unique)
        idx = findall(abs.(r .- rv) .< 1e-10)
        T_avg[i] = mean(T[idx])
    end
    return r_unique, T_avg
end

function angular_profile(mesh, T)
    theta = atan.(mesh.node[:, 2], mesh.node[:, 1])
    theta = mod.(theta, 2.0 * pi)
    theta_unique = unique(sort(theta))
    T_avg = zeros(Float64, length(theta_unique))
    for (i, tv) in enumerate(theta_unique)
        idx = findall(abs.(theta .- tv) .< 1e-10)
        T_avg[i] = mean(T[idx])
    end
    return theta_unique, T_avg
end

function analytical_solution_ring(r, q, k, h, r_i, r_o, T_f)
    term1 = (q / (4.0 * k)) * (r_o^2 - r^2)
    term2 = (q / (2.0 * k)) * r_i^2 * log(r / r_o)
    term3 = (q / (2.0 * h)) * (r_o - (r_i^2 / r_o))
    return T_f + term1 + term2 + term3
end

function rt_q4_shape(xi, eta)
    N = 0.25 * [(1 - xi) * (1 - eta);
                (1 + xi) * (1 - eta);
                (1 + xi) * (1 + eta);
                (1 - xi) * (1 + eta)]
    dN_dxi = 0.25 * [-(1 - eta)    -(1 - xi);
                      (1 - eta)    -(1 + xi);
                      (1 + eta)     (1 + xi);
                     -(1 + eta)     (1 - xi)]
    return N, dN_dxi
end

function rt_mesh(Rin, Rout, ntheta, dr)
    nr = max(1, round(Int, (Rout - Rin) / dr))
    r_nodes = collect(range(Rin, Rout; length=nr + 1))
    theta_nodes = collect(range(0.0, 2.0 * pi; length=ntheta + 1))[1:ntheta]
    return r_nodes, theta_nodes, nr
end

function rt_connectivity(nr, ntheta)
    ne = nr * ntheta
    element = zeros(Int64, ne, 4)
    e = 0
    for ir in 1:nr
        for it in 1:ntheta
            it_next = it == ntheta ? 1 : it + 1
            e += 1
            n1 = (ir - 1) * ntheta + it
            n2 = ir * ntheta + it
            n3 = ir * ntheta + it_next
            n4 = (ir - 1) * ntheta + it_next
            element[e, 1] = n1
            element[e, 2] = n2
            element[e, 3] = n3
            element[e, 4] = n4
        end
    end
    return element
end

function rt_fem_assemble(r_nodes, theta_nodes, element, k_r, k_t, q0, Bi, T_amb)
    ntheta = length(theta_nodes)
    nr = length(r_nodes) - 1
    nnode = length(r_nodes) * ntheta
    I = Int[]
    J = Int[]
    V = Float64[]
    F = zeros(Float64, nnode)
    M = zeros(Float64, nnode)

    gp = [-0.577350269189626, 0.577350269189626]
    gw = [1.0, 1.0]

    # Volume integrals
    for e in 1:size(element, 1)
        nodes = element[e, :]
        r_e = [r_nodes[ceil(Int, n / ntheta)] for n in nodes]
        t_e = [theta_nodes[mod(n - 1, ntheta) + 1] for n in nodes]

        for (xi, wx) in zip(gp, gw)
            for (eta, wy) in zip(gp, gw)
                N, dN_dxi = rt_q4_shape(xi, eta)
                Jmat = transpose(dN_dxi) * hcat(r_e, t_e)
                detJ = det(Jmat)
                invJ = inv(Jmat)
                dN_drt = dN_dxi * invJ
                dN_dr = dN_drt[:, 1]
                dN_dt = dN_drt[:, 2]

                r_g = dot(N, r_e)
                w = wx * wy * detJ

                for a in 1:4
                    ia = nodes[a]
                    Na = N[a]
                    M[ia] += Na * Na * r_g * w
                    F[ia] += q0 * Na * r_g * w
                    for b in 1:4
                        ib = nodes[b]
                        kb = k_r * dN_dr[a] * dN_dr[b] * r_g
                        kt = k_t * dN_dt[a] * dN_dt[b] / r_g
                        push!(I, ia); push!(J, ib); push!(V, -(kb + kt) * w)
                    end
                end
            end
        end
    end

    # Outer convection (r = Rout)
    for it in 1:ntheta
        it_next = it == ntheta ? 1 : it + 1
        n1 = nr * ntheta + it
        n2 = nr * ntheta + it_next
        r_edge = r_nodes[end]
        t1 = theta_nodes[it]
        t2 = theta_nodes[it_next]
        for (s, w) in zip(gp, gw)
            N1 = 0.5 * (1 - s)
            N2 = 0.5 * (1 + s)
            theta_g = 0.5 * (1 - s) * t1 + 0.5 * (1 + s) * t2
            dtheta_ds = 0.5 * (t2 - t1)
            Jedge = r_edge * dtheta_ds
            wt = w * Jedge
            for (a, Na, ia) in ((1, N1, n1), (2, N2, n2))
                F[ia] += Bi * T_amb * Na * wt
                for (b, Nb, ib) in ((1, N1, n1), (2, N2, n2))
                    push!(I, ia); push!(J, ib); push!(V, -Bi * Na * Nb * wt)
                end
            end
        end
    end

    MT = spdiagm(0 => M)
    KT = sparse(I, J, V, nnode, nnode)
    return MT, KT, F
end

function run_rtheta_fem_case(param_dim, ntheta, dr, label)
    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout
    scale = param_dim.scale
    param = JuBat.NormaliseParam(param_dim)

    Rin_nd = param.cell.Rin
    Rout_nd = param.cell.Rout
    dr_nd = dr / scale.L
    r_nodes, theta_nodes, nr = rt_mesh(Rin_nd, Rout_nd, ntheta, dr_nd)
    element = rt_connectivity(nr, ntheta)

    k_r = param.cell.lambda_r
    k_t = param.cell.lambda_t
    Bi = scale.h
    T_f_nd = param.cell.T_amb
    q0_nd = 2.0e5 / scale.q

    MT, KT, F = rt_fem_assemble(r_nodes, theta_nodes, element, k_r, k_t, q0_nd, Bi, T_f_nd)
    T_nd = -(KT \ F)
    T = T_nd .* scale.T_ref

    # Build cartesian coords for plotting
    nnode = length(r_nodes) * ntheta
    x = zeros(Float64, nnode)
    y = zeros(Float64, nnode)
    for ir in 1:length(r_nodes)
        for it in 1:ntheta
            idx = (ir - 1) * ntheta + it
            x[idx] = r_nodes[ir] * scale.L * cos(theta_nodes[it])
            y[idx] = r_nodes[ir] * scale.L * sin(theta_nodes[it])
        end
    end

    ang_range, ang_std = angular_variation((node = hcat(x, y),), T)
    return (ang_range=ang_range, ang_std=ang_std)
end


function plot_mesh_outline(mesh, out_dir)
    x = mesh.node[:, 1]
    y = mesh.node[:, 2]
    plt = Plots.plot(size=(700, 700), aspect_ratio=1, legend=false,
        xlabel="x [m]", ylabel="y [m]", title="Ring Mesh Outline")
    ne = size(mesh.element, 1)
    for e in 1:ne
        n = mesh.element[e, :]
        xs = [x[n[1]], x[n[2]], x[n[3]], x[n[4]], x[n[1]]]
        ys = [y[n[1]], y[n[2]], y[n[3]], y[n[4]], y[n[1]]]
        Plots.plot!(plt, xs, ys, lw=0.4, color=:black, alpha=0.7)
    end
    savefig(plt, joinpath(out_dir, "ring_mesh_outline.png"))
end

function quad_to_triangles(element)
    ne = size(element, 1)
    tris = Array{Int64}(undef, ne * 2, 3)
    for e in 1:ne
        n1, n2, n3, n4 = element[e, 1], element[e, 2], element[e, 3], element[e, 4]
        tris[2 * e - 1, :] = [n1, n2, n3]
        tris[2 * e, :] = [n1, n3, n4]
    end
    return tris
end

function tripcolor_field(x, y, tris, values, out_path; title="", cmap="hot", Rin=0.0, Rout=Inf)
    fig = PyPlot.figure(figsize=(6.5, 5.0))
    ax = fig.add_subplot(1, 1, 1)
    tri_mod = PyPlot.pyimport("matplotlib.tri")
    tri = tri_mod.Triangulation(x, y, tris .- 1)
    cx = (x[tris[:, 1]] + x[tris[:, 2]] + x[tris[:, 3]]) / 3.0
    cy = (y[tris[:, 1]] + y[tris[:, 2]] + y[tris[:, 3]]) / 3.0
    rc = hypot.(cx, cy)
    mask = (rc .< Rin) .| (rc .> Rout)
    tri.set_mask(mask)
    tpc = ax.tripcolor(tri, values, shading="gouraud", cmap=cmap)
    fig.colorbar(tpc, ax=ax)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    fig.savefig(out_path, dpi=200)
    PyPlot.close(fig)
end

function mesh_diagnostics(mesh, out_dir)
    return nothing
end

function angular_variation(mesh, T)
    r = hypot.(mesh.node[:, 1], mesh.node[:, 2])
    r_unique = unique(sort(r))
    max_range = 0.0
    max_std = 0.0
    for rv in r_unique
        idx = findall(abs.(r .- rv) .< 1e-10)
        Tv = T[idx]
        max_range = max(max_range, maximum(Tv) - minimum(Tv))
        max_std = max(max_std, std(Tv))
    end
    return max_range, max_std
end

function radius_variation_report(mesh, T; topn=8, tol=1e-10)
    r = hypot.(mesh.node[:, 1], mesh.node[:, 2])
    r_unique = unique(sort(r))
    rows = Vector{Tuple{Float64, Float64, Float64, Int}}()
    for rv in r_unique
        idx = findall(abs.(r .- rv) .< tol)
        Tv = T[idx]
        push!(rows, (rv, maximum(Tv) - minimum(Tv), std(Tv), length(idx)))
    end
    sort!(rows, by = x -> x[2], rev = true)
    return rows[1:min(topn, length(rows))]
end

function steady_residual(case, vars, T_vec, model, mesh_data)
    if model == "ring2D_polar"
        MT, KT, FT = JuBat.ThermalPolar2D_Ring(case, vars, mesh_data)
    else
        MT, KT, FT = JuBat.ThermalDistributed2D_Ring(case, vars)
        outer_nodes = get(vars, "thermal2D outer_nodes", Int[])
        isempty(outer_nodes) || (KT, FT = JuBat.ThermalRing2D_BC(KT, FT, case, outer_nodes, 0.0))
    end
    r = KT * T_vec + FT
    r_l2 = sqrt(mean(r .^ 2))
    r_linf = maximum(abs.(r))
    return r_l2, r_linf
end

function steady_solve(case, vars, model, mesh_data)
    if model == "ring2D_polar"
        MT, KT, FT = JuBat.ThermalPolar2D_Ring(case, vars, mesh_data)
    else
        MT, KT, FT = JuBat.ThermalDistributed2D_Ring(case, vars)
        outer_nodes = get(vars, "thermal2D outer_nodes", Int[])
        isempty(outer_nodes) || (KT, FT = JuBat.ThermalRing2D_BC(KT, FT, case, outer_nodes, 0.0))
    end
    return -(KT \ FT)
end

function solve_for_ntheta(param_dim, ntheta, dr, model)
    opt = JuBat.Option()
    opt.model = "thermal"
    opt.thermal_enabled = true
    opt.thermalmodel = model
    opt.time = [0.0, 3600]
    opt.dt = [1.0, 10]
    case = JuBat.SetCase(param_dim, opt)

    scale = param_dim.scale
    T_ref = scale.T_ref
    q_ref = scale.q  # 统一能量尺度热源参考 (P_ref / L^3)

    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout
    nr = round(Int, (Rout - Rin) / dr)
    mesh_data = JuBat.ring_mesh(case.param, ntheta=ntheta, nr=nr, gsorder=2)
    mesh = mesh_data.mesh
    case.mesh["thermal2D"] = mesh

    q0 = 2.0e5
    q0_nd = q0 / q_ref
    q_func = (r, theta, t) -> q0

    variables = Dict{String,Any}()
    variables["T_nodes"] = fill(param_dim.cell.T0 / T_ref, mesh.nlen)
    variables["thermal2D outer_nodes"] = mesh_data.outer_nodes
    if model == "ring2D_polar"
        ne = size(mesh.element, 1)
        variables["heat_source_fields"] = fill(q0_nd, ne)
    end

    update_fn = (t, vars) -> begin
        if model == "ring2D_polar"
            vars["heat_source_fields"] = fill(q0_nd, ne)
        else
            vars["heat_source_fields"] = compute_q_elem(mesh, q_func, t) ./ q_ref
        end
    end

    case.multi_spme_layout["thermal_variables"] = variables
    case.multi_spme_layout["thermal_update_fn"] = update_fn
    case.multi_spme_layout["thermal_record"] = false
    case.multi_spme_layout["polar_mesh_data"] = mesh_data

    result = JuBat.Solve(case)
    T = result.T_nodes .* T_ref

    # Convert normalized coordinates back to dimensional for exact solution comparison
    r_nodes_dim = hypot.(mesh.node[:, 1], mesh.node[:, 2]) .* scale.L
    k_r = param_dim.cell.lambda_r
    h = param_dim.cell.h
    T_f = param_dim.cell.T_amb
    T_exact = analytical_solution_ring.(r_nodes_dim, q0, k_r, h, Rin, Rout, T_f)
    err = T .- T_exact
    err_l2 = sqrt(mean(err .^ 2))
    err_linf = maximum(abs.(err))
    err_rel = err_l2 / max(1e-12, maximum(abs.(T_exact)))
    ang_range, ang_std = angular_variation(mesh, T)

    return (ang_range=ang_range, ang_std=ang_std, err_l2=err_l2, err_linf=err_linf, err_rel=err_rel)
end

function run_case(param_dim, ntheta, dr, model, label)
    opt = JuBat.Option()
    opt.model = "thermal"
    opt.thermal_enabled = true
    opt.thermalmodel = model
    opt.time = [0.0, 3600]
    opt.dt = [1.0, 10]

    case = JuBat.SetCase(param_dim, opt)

    scale = param_dim.scale
    T_ref = scale.T_ref
    q_ref = scale.q  # 统一能量尺度热源参考 (P_ref / L^3)

    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout
    nr = round(Int, (Rout - Rin) / dr)
    mesh_data = JuBat.ring_mesh(case.param, ntheta=ntheta, nr=nr, gsorder=2)
    mesh = mesh_data.mesh
    case.mesh["thermal2D"] = mesh

    q0 = 2.0e5
    q0_nd = q0 / q_ref
    q_func = (r, theta, t) -> q0

    variables = Dict{String,Any}()
    variables["T_nodes"] = fill(param_dim.cell.T0 / T_ref, mesh.nlen)
    variables["thermal2D outer_nodes"] = mesh_data.outer_nodes
    if model == "ring2D_polar"
        ne = size(mesh.element, 1)
        variables["heat_source_fields"] = fill(q0_nd, ne)
    end

    update_fn = (t, vars) -> begin
        if model == "ring2D_polar"
            vars["heat_source_fields"] = fill(q0_nd, ne)
        else
            vars["heat_source_fields"] = compute_q_elem(mesh, q_func, t) ./ q_ref
        end
    end

    case.multi_spme_layout["thermal_variables"] = variables
    case.multi_spme_layout["thermal_update_fn"] = update_fn
    case.multi_spme_layout["thermal_record"] = false
    case.multi_spme_layout["polar_mesh_data"] = mesh_data

    result = JuBat.Solve(case)
    T = result.T_nodes .* T_ref

    out_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "thermal_verify", label))
    isdir(out_dir) || mkpath(out_dir)
    plot_mesh_outline(mesh, out_dir)

    r_prof, T_r = radial_profile(mesh, T)
    r_prof_dim = r_prof .* scale.L  # Convert to dimensional for exact solution
    k_r = param_dim.cell.lambda_r
    h = param_dim.cell.h
    T_f = param_dim.cell.T_amb
    T_r_exact = analytical_solution_ring.(r_prof_dim, q0, k_r, h, Rin, Rout, T_f)
    p_r = Plots.plot(r_prof_dim, T_r, xlabel="r [m]", ylabel="T [K]", lw=2, label="FEM",
        title="Radial Temperature Profile")
    Plots.plot!(p_r, r_prof_dim, T_r_exact, lw=2, linestyle=:dash, label="Exact")
    savefig(p_r, joinpath(out_dir, "ring_temperature_radial.png"))

    x = mesh.node[:, 1]
    y = mesh.node[:, 2]
    tris = quad_to_triangles(mesh.element)
    Rin_nd = Rin / scale.L  # Normalized for mask calculation
    Rout_nd = Rout / scale.L
    tripcolor_field(x, y, tris, T, joinpath(out_dir, "ring_temperature_field.png");
        title="Ring Temperature Field", cmap="hot", Rin=Rin_nd, Rout=Rout_nd)

    theta_prof, T_theta = angular_profile(mesh, T)
    p_theta = Plots.plot(theta_prof, T_theta, xlabel="theta [rad]", ylabel="T [K]", lw=2,
        title="Angular Temperature Profile")
    savefig(p_theta, joinpath(out_dir, "ring_temperature_angular.png"))

    r_nodes_dim = hypot.(mesh.node[:, 1], mesh.node[:, 2]) .* scale.L
    T_exact = analytical_solution_ring.(r_nodes_dim, q0, k_r, h, Rin, Rout, T_f)
    err = T .- T_exact
    err_l2 = sqrt(mean(err .^ 2))
    err_linf = maximum(abs.(err))
    err_rel = err_l2 / max(1e-12, maximum(abs.(T_exact)))

    update_fn !== nothing && update_fn(opt.time[end], variables)
    T_exact_nd = T_exact ./ T_ref
    res_l2, res_linf = steady_residual(case, variables, T_exact_nd, model, mesh_data)
    T_steady = steady_solve(case, variables, model, mesh_data) .* T_ref
    err_s = T_steady .- T_exact
    err_s_l2 = sqrt(mean(err_s .^ 2))
    err_s_linf = maximum(abs.(err_s))
    err_s_rel = err_s_l2 / max(1e-12, maximum(abs.(T_exact)))
    diff_ts = T .- T_steady
    diff_ts_l2 = sqrt(mean(diff_ts .^ 2))
    ang_range, ang_std = angular_variation(mesh, T)
    worst_r = radius_variation_report(mesh, T)

    tripcolor_field(x, y, tris, err, joinpath(out_dir, "ring_temperature_error.png");
        title="Temperature Error (T - T_exact)", cmap="coolwarm", Rin=Rin_nd, Rout=Rout_nd)


    return (out_dir=out_dir, ang_range=ang_range, ang_std=ang_std)
end

function run_case_data(param_dim, ntheta, dr, model)
    opt = JuBat.Option()
    opt.model = "thermal"
    opt.thermal_enabled = true
    opt.thermalmodel = model
    opt.time = [0.0, 3600]
    opt.dt = [1.0, 10]

    case = JuBat.SetCase(param_dim, opt)

    scale = param_dim.scale
    T_ref = scale.T_ref
    q_ref = scale.q  # 统一能量尺度热源参考 (P_ref / L^3)

    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout
    nr = round(Int, (Rout - Rin) / dr)
    mesh_data = JuBat.ring_mesh(case.param, ntheta=ntheta, nr=nr, gsorder=2)
    mesh = mesh_data.mesh
    case.mesh["thermal2D"] = mesh

    q0 = 2.0e5
    q0_nd = q0 / q_ref
    q_func = (r, theta, t) -> q0

    variables = Dict{String,Any}()
    variables["T_nodes"] = fill(param_dim.cell.T0 / T_ref, mesh.nlen)
    variables["thermal2D outer_nodes"] = mesh_data.outer_nodes
    if model == "ring2D_polar"
        ne = size(mesh.element, 1)
        variables["heat_source_fields"] = fill(q0_nd, ne)
    end

    update_fn = (t, vars) -> begin
        if model == "ring2D_polar"
            vars["heat_source_fields"] = fill(q0_nd, ne)
        else
            vars["heat_source_fields"] = compute_q_elem(mesh, q_func, t) ./ q_ref
        end
    end

    case.multi_spme_layout["thermal_variables"] = variables
    case.multi_spme_layout["thermal_update_fn"] = update_fn
    case.multi_spme_layout["thermal_record"] = false
    case.multi_spme_layout["polar_mesh_data"] = mesh_data

    result = JuBat.Solve(case)
    T = result.T_nodes .* T_ref

    r_prof, T_r = radial_profile(mesh, T)
    r_prof_dim = r_prof .* scale.L  # Convert to dimensional for exact solution
    k_r = param_dim.cell.lambda_r
    h = param_dim.cell.h
    T_f = param_dim.cell.T_amb
    T_r_exact = analytical_solution_ring.(r_prof_dim, q0, k_r, h, Rin, Rout, T_f)
    ang_range, ang_std = angular_variation(mesh, T)

    r_nodes_dim = hypot.(mesh.node[:, 1], mesh.node[:, 2]) .* scale.L
    T_exact = analytical_solution_ring.(r_nodes_dim, q0, k_r, h, Rin, Rout, T_f)
    err = T .- T_exact
    err_l2 = sqrt(mean(err .^ 2))
    err_linf = maximum(abs.(err))
    err_rel = err_l2 / max(1e-12, maximum(abs.(T_exact)))

    update_fn !== nothing && update_fn(opt.time[end], variables)
    T_exact_nd = T_exact ./ T_ref
    res_l2, res_linf = steady_residual(case, variables, T_exact_nd, model, mesh_data)

    return (mesh=mesh, mesh_data=mesh_data, T=T, r_prof=r_prof_dim, T_r=T_r,
        T_r_exact=T_r_exact, ang_range=ang_range, ang_std=ang_std,
        err_l2=err_l2, err_linf=err_linf, err_rel=err_rel,
        res_l2=res_l2, res_linf=res_linf,
        Tmin=minimum(T), Tmax=maximum(T))
end

function run_rtheta_fem_data(param_dim, ntheta, dr)
    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout
    scale = param_dim.scale
    param = JuBat.NormaliseParam(param_dim)

    Rin_nd = param.cell.Rin
    Rout_nd = param.cell.Rout
    dr_nd = dr / scale.L
    r_nodes, theta_nodes, nr = rt_mesh(Rin_nd, Rout_nd, ntheta, dr_nd)
    element = rt_connectivity(nr, ntheta)

    k_r_nd = param.cell.lambda_r
    k_t_nd = param.cell.lambda_t
    Bi = scale.h
    T_f_nd = param.cell.T_amb
    q0 = 2.0e5
    q0_nd = q0 / scale.q

    MT, KT, F = rt_fem_assemble(r_nodes, theta_nodes, element, k_r_nd, k_t_nd, q0_nd, Bi, T_f_nd)
    T_nd = -(KT \ F)
    T = T_nd .* scale.T_ref

    nnode = length(r_nodes) * ntheta
    x = zeros(Float64, nnode)
    y = zeros(Float64, nnode)
    for ir in 1:length(r_nodes)
        for it in 1:ntheta
            idx = (ir - 1) * ntheta + it
            x[idx] = r_nodes[ir] * scale.L * cos(theta_nodes[it])
            y[idx] = r_nodes[ir] * scale.L * sin(theta_nodes[it])
        end
    end

    r_prof = r_nodes .* scale.L
    T_r = [mean(T[(ir - 1) * ntheta + 1:ir * ntheta]) for ir in 1:length(r_nodes)]
    k_r = param_dim.cell.lambda_r
    h = param_dim.cell.h
    T_f = param_dim.cell.T_amb
    T_r_exact = analytical_solution_ring.(r_prof, q0, k_r, h, Rin, Rout, T_f)
    ang_range, ang_std = angular_variation((node = hcat(x, y),), T)

    r_nodes_dim = hypot.(x, y)
    T_exact = analytical_solution_ring.(r_nodes_dim, q0, k_r, h, Rin, Rout, T_f)
    err = T .- T_exact
    err_l2 = sqrt(mean(err .^ 2))
    err_linf = maximum(abs.(err))
    err_rel = err_l2 / max(1e-12, maximum(abs.(T_exact)))

    res = KT * T_nd + F
    res_l2 = sqrt(mean(res .^ 2))
    res_linf = maximum(abs.(res))

    return (x=x, y=y, element=element, T=T, r_prof=r_prof, T_r=T_r,
        T_r_exact=T_r_exact, ang_range=ang_range, ang_std=ang_std,
        err_l2=err_l2, err_linf=err_linf, err_rel=err_rel,
        res_l2=res_l2, res_linf=res_linf,
        Tmin=minimum(T), Tmax=maximum(T))
end

function print_check_summary(tag, data)
    println("\n[Check] " * tag)
    @printf("  T range [K]         : [%.6f, %.6f]\n", data.Tmin, data.Tmax)
    @printf("  Error L2 [K]        : %.6e\n", data.err_l2)
    @printf("  Error Linf [K]      : %.6e\n", data.err_linf)
    @printf("  Error Rel [-]       : %.6e\n", data.err_rel)
    @printf("  Angular range [K]   : %.6e\n", data.ang_range)
    @printf("  Angular std [K]     : %.6e\n", data.ang_std)
    @printf("  Residual L2 [-]     : %.6e\n", data.res_l2)
    @printf("  Residual Linf [-]   : %.6e\n", data.res_linf)
end

function print_consistency_checks(fem, polar, fem_rt)
    println("\n[Check] Cross-method consistency")
    l2_polar_rt = sqrt(mean((polar.T .- fem_rt.T).^2))
    linf_polar_rt = maximum(abs.(polar.T .- fem_rt.T))
    l2_fem_polar = sqrt(mean((fem.T .- polar.T).^2))
    linf_fem_polar = maximum(abs.(fem.T .- polar.T))

    @printf("  Polar vs r-theta FEM L2 [K]   : %.6e\n", l2_polar_rt)
    @printf("  Polar vs r-theta FEM Linf [K] : %.6e\n", linf_polar_rt)
    @printf("  Q4 FEM vs Polar L2 [K]        : %.6e\n", l2_fem_polar)
    @printf("  Q4 FEM vs Polar Linf [K]      : %.6e\n", linf_fem_polar)

    if l2_polar_rt < 1e-2 && polar.err_rel < 0.2 && fem_rt.err_rel < 0.2
        println("  Status: PASS (polar and r-theta FEM are mutually consistent)")
    else
        println("  Status: WARN (check heat source / boundary scaling)")
    end

    if fem.err_rel > 0.2 && l2_fem_polar > 5.0
        println("  Q4 FEM Note: large deviation from polar/r-theta baseline, likely assembly/model mismatch instead of unit scaling.")
    end
end


function main()
    param_dim = JuBat.ChooseCell("Ring")
    scale = param_dim.scale
    Rin = param_dim.cell.Rin
    Rout = param_dim.cell.Rout
    Rin_nd = Rin / scale.L  # Normalized for mask calculation
    Rout_nd = Rout / scale.L
    ntheta = 40
    dr = (Rout - Rin) / 20
    out_root = normpath(joinpath(@__DIR__, "..", "..", "output", "thermal_verify"))
    isdir(out_root) || mkpath(out_root)

    fem = run_case_data(param_dim, ntheta, dr, "ring2D")
    polar = run_case_data(param_dim, ntheta, dr, "ring2D_polar")
    fem_rt = run_rtheta_fem_data(param_dim, ntheta, dr)

    print_check_summary("FEM (Q4)", fem)
    print_check_summary("Polar FVM", polar)
    print_check_summary("r-theta FEM", fem_rt)
    print_consistency_checks(fem, polar, fem_rt)

    # Mesh outline (FEM mesh)
    plot_mesh_outline(fem.mesh, out_root)

    # Temperature fields (three methods) - use normalized Rin/Rout for mask
    tris_fem = quad_to_triangles(fem.mesh.element)
    tripcolor_field(fem.mesh.node[:, 1], fem.mesh.node[:, 2], tris_fem, fem.T,
        joinpath(out_root, "ring_temperature_field_fem.png"); title="FEM (Q4) Temperature", cmap="hot", Rin=Rin_nd, Rout=Rout_nd)
    tripcolor_field(fem.mesh.node[:, 1], fem.mesh.node[:, 2], tris_fem, polar.T,
        joinpath(out_root, "ring_temperature_field_polar.png"); title="Polar FVM Temperature", cmap="hot", Rin=Rin_nd, Rout=Rout_nd)
    tris_rt = quad_to_triangles(fem_rt.element)
    tripcolor_field(fem_rt.x, fem_rt.y, tris_rt, fem_rt.T,
        joinpath(out_root, "ring_temperature_field_fem_rt.png"); title="r-θ FEM Temperature", cmap="hot", Rin=Rin, Rout=Rout)

    # Error comparison curve (radial)
    p_err = Plots.plot(fem.r_prof, abs.(fem.T_r .- fem.T_r_exact), lw=2, label="FEM")
    Plots.plot!(p_err, polar.r_prof, abs.(polar.T_r .- polar.T_r_exact), lw=2, label="Polar FVM")
    Plots.plot!(p_err, fem_rt.r_prof, abs.(fem_rt.T_r .- fem_rt.T_r_exact), lw=2, label="r-θ FEM")
    Plots.plot!(p_err, xlabel="r [m]", ylabel="|T - T_exact| [K]", title="Radial Error Comparison")
    savefig(p_err, joinpath(out_root, "ring_error_radial_compare.png"))

    # Convergence study for angular variation vs ntheta
    ntheta_list = [20, 40, 80, 120]
    conv_dir = out_root

    # Convergence plot
    fem_ang = Float64[]
    polar_ang = Float64[]
    femrt_ang = Float64[]
    for nθ in ntheta_list
        push!(fem_ang, solve_for_ntheta(param_dim, nθ, dr, "ring2D").ang_range)
        push!(polar_ang, solve_for_ntheta(param_dim, nθ, dr, "ring2D_polar").ang_range)
        push!(femrt_ang, run_rtheta_fem_data(param_dim, nθ, dr).ang_range)
    end
    p_conv = Plots.plot(ntheta_list, fem_ang, marker=:o, label="FEM")
    Plots.plot!(p_conv, ntheta_list, femrt_ang, marker=:o, label="r-θ FEM")
    Plots.plot!(p_conv, ntheta_list, polar_ang, marker=:o, label="Polar FVM")
    Plots.plot!(p_conv, xlabel="ntheta", ylabel="Angular range [K]", title="Angular Symmetry Convergence")
    savefig(p_conv, joinpath(out_root, "ring_convergence_angular.png"))
end

main()
