"""
    ThermalPolar2D_Ring(case, variables, mesh_data)

Polar finite-volume assembly for ring geometry (r-theta grid).
Uses nondimensional units with KT containing the negative sign convention:
    M dT/dt = KT T + F

Expected:
- mesh_data from ring_mesh (r, ntheta, nr, etc.)
- variables["heat_source_fields"] element source (W/m^3)
"""
function ThermalPolar2D_Ring(case::Case, variables::Dict{String,Any}, mesh_data)
    mesh = case.mesh["thermal2D"]
    r_nodes = mesh_data.r
    ntheta = mesh_data.ntheta
    scale = case.param_dim.scale
    L_th = scale.L_th
    k_th = scale.k_th
    rho_c_th = scale.rho_c_th
    Bi = scale.h_th
    T_ref = scale.T_ref

    nnode = length(r_nodes) * ntheta
    dtheta = 2.0 * pi / ntheta

    rho = case.param_dim.cell.rho
    cp = case.param_dim.cell.heat_Q
    k_r = case.param_dim.cell.lambda_r
    k_t = case.param_dim.cell.lambda_t
    T_amb = case.param_dim.cell.T_amb / T_ref
    rho_c_nd = (rho * cp) / rho_c_th
    k_r_nd = k_r / k_th
    k_t_nd = k_t / k_th
    r_nodes_nd = r_nodes ./ L_th

    I = Int[]
    J = Int[]
    V = Float64[]
    F = zeros(Float64, nnode)
    M = zeros(Float64, nnode)

    q_node = zeros(Float64, nnode)
    q_elem = get(variables, "heat_source_fields", nothing)
    q_nodes = get(variables, "heat_source_nodes", nothing)
    if q_elem !== nothing
        counts = zeros(Int, nnode)
        for e in 1:size(mesh.element, 1)
            for n in mesh.element[e, :]
                q_node[n] += q_elem[e]
                counts[n] += 1
            end
        end
        for i in 1:nnode
            counts[i] > 0 && (q_node[i] /= counts[i])
        end
    elseif q_nodes !== nothing
        q_node = q_nodes
        length(q_node) == nnode || error("heat_source_nodes length mismatch")
    else
        error("Missing heat source: set heat_source_fields or heat_source_nodes")
    end

    for ir in 1:length(r_nodes)
        r_i = r_nodes_nd[ir]
        r_imh = ir == 1 ? r_nodes_nd[1] : 0.5 * (r_nodes_nd[ir - 1] + r_i)
        r_iph = ir == length(r_nodes) ? r_nodes_nd[end] : 0.5 * (r_i + r_nodes_nd[ir + 1])
        dr_im = ir == 1 ? (r_i - r_nodes_nd[1]) : (r_i - r_nodes_nd[ir - 1])
        dr_ip = ir == length(r_nodes) ? (r_nodes_nd[end] - r_i) : (r_nodes_nd[ir + 1] - r_i)
        vol = 0.5 * (r_iph^2 - r_imh^2) * dtheta
        area_theta = r_iph - r_imh

        for it in 1:ntheta
            idx = (ir - 1) * ntheta + it
            M[idx] = rho_c_nd * vol
            F[idx] += q_node[idx] * vol

            if ir > 1
                a_rm = k_r_nd * r_imh * dtheta / dr_im
                push!(I, idx); push!(J, idx); push!(V, -a_rm)
                idxm = (ir - 2) * ntheta + it
                push!(I, idx); push!(J, idxm); push!(V, a_rm)
            end
            if ir < length(r_nodes)
                a_rp = k_r_nd * r_iph * dtheta / dr_ip
                push!(I, idx); push!(J, idx); push!(V, -a_rp)
                idxp = ir * ntheta + it
                push!(I, idx); push!(J, idxp); push!(V, a_rp)
            else
                area = r_iph * dtheta
                push!(I, idx); push!(J, idx); push!(V, -Bi * area)
                F[idx] += Bi * area * T_amb
            end

            itp = it == ntheta ? 1 : it + 1
            itm = it == 1 ? ntheta : it - 1
            a_t = k_t_nd * area_theta / (r_i * dtheta)
            push!(I, idx); push!(J, idx); push!(V, -2.0 * a_t)
            push!(I, idx); push!(J, (ir - 1) * ntheta + itp); push!(V, a_t)
            push!(I, idx); push!(J, (ir - 1) * ntheta + itm); push!(V, a_t)
        end
    end

    MT = spdiagm(0 => M)
    KT = sparse(I, J, V, nnode, nnode)
    return MT, KT, F
end
