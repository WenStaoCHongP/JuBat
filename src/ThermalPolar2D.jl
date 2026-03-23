function ThermalPolar2D_Ring(case::Case, variables::Dict{String,Any}, mesh_data)
    mesh = case.mesh["thermal2D"]
    param = case.param  # Use normalized parameters
    scale = case.param.scale

    r_nodes = mesh_data.r  # Already normalized (ring_mesh called with L=scale.L)
    ntheta = mesh_data.ntheta

    # Grid spacing (in radians)
    dtheta = 2.0 * pi / ntheta

    # --- Normalized thermal properties (from case.param.cell) ---
    # Volumetric heat capacity: (ρc)* = C* / V* = heat_Q / volume
    # where C* is total heat capacity and V* is normalized volume
    V_nd = param.cell.volume  # V* = V / L³
    C_nd = param.cell.heat_Q  # C* = m·c·T_ref/(t0·P_ref)
    V_nd > 0 || error("ThermalPolar2D_Ring: param.cell.volume must be positive, got $V_nd")
    rho_c_nd = C_nd / V_nd  # (ρc)* = C*/V*

    # Anisotropic conductivities (already normalized)
    k_r_nd = param.cell.lambda_r  # k_r* = k_r · L · T_ref / P_ref
    k_t_nd = param.cell.lambda_t  # k_t* = k_t · L · T_ref / P_ref

    # Boundary condition parameters
    Bi = scale.h * param.cell.lambda_r
    T_amb_nd = param.cell.T_amb  # Normalized ambient temperature

    nnode = length(r_nodes) * ntheta

    # Sparse matrix assembly
    I = Int[]
    J = Int[]
    V = Float64[]
    F = zeros(Float64, nnode)
    M = zeros(Float64, nnode)

    # --- Heat source mapping: element -> node ---
    q_elem = variables["heat_source_fields"]
    q_node = zeros(Float64, nnode)
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

    # --- FVM assembly over radial rings ---
    for ir in 1:length(r_nodes)
        r_i = r_nodes[ir]

        # Cell face positions (halfway between nodes)
        r_imh = ir == 1 ? r_nodes[1] : 0.5 * (r_nodes[ir - 1] + r_i)
        r_iph = ir == length(r_nodes) ? r_nodes[end] : 0.5 * (r_i + r_nodes[ir + 1])

        # Radial distances to neighbors
        dr_im = ir == 1 ? (r_i - r_nodes[1]) : (r_i - r_nodes[ir - 1])
        dr_ip = ir == length(r_nodes) ? (r_nodes[end] - r_i) : (r_nodes[ir + 1] - r_i)

        # Cell volume (in 2D polar: area per unit height, normalized)
        # V* = 0.5 · (r_iph² - r_imh²) · Δθ
        vol_nd = 0.5 * (r_iph^2 - r_imh^2) * dtheta

        # Radial face length for tangential conduction
        area_theta_nd = r_iph - r_imh

        for it in 1:ntheta
            idx = (ir - 1) * ntheta + it

            # --- Mass matrix: M[idx] = (ρc)* · V* ---
            M[idx] = rho_c_nd * vol_nd

            # --- Load vector: F[idx] += q* · V* ---
            F[idx] += q_node[idx] * vol_nd

            # --- Radial conduction (inner neighbor) ---
            if ir > 1
                # Conductance: a = k* · r_face · Δθ / Δr
                dr_im > 0 || error("ThermalPolar2D_Ring: radial spacing dr_im must be positive at ir=$ir, got $dr_im")
                a_rm = k_r_nd * r_imh * dtheta / dr_im
                push!(I, idx); push!(J, idx); push!(V, -a_rm)
                idxm = (ir - 2) * ntheta + it
                push!(I, idx); push!(J, idxm); push!(V, a_rm)
            end

            # --- Radial conduction (outer neighbor) ---
            if ir < length(r_nodes)
                dr_ip > 0 || error("ThermalPolar2D_Ring: radial spacing dr_ip must be positive at ir=$ir, got $dr_ip")
                a_rp = k_r_nd * r_iph * dtheta / dr_ip
                push!(I, idx); push!(J, idx); push!(V, -a_rp)
                idxp = ir * ntheta + it
                push!(I, idx); push!(J, idxp); push!(V, a_rp)
            else
                # --- Convection BC at outer boundary (r = R_out) ---
                # Boundary area: A* = r_iph · Δθ (normalized arc length)
                area_bc_nd = r_iph * dtheta
                push!(I, idx); push!(J, idx); push!(V, -Bi * area_bc_nd)
                F[idx] += Bi * area_bc_nd * T_amb_nd
            end

            # --- Tangential conduction (periodic in θ) ---
            itp = it == ntheta ? 1 : it + 1
            itm = it == 1 ? ntheta : it - 1
            # Conductance: a = k_t* · Δr_face / (r · Δθ)
            r_i > 0 || error("ThermalPolar2D_Ring: radial position r_i must be positive at ir=$ir, got $r_i")
            a_t = k_t_nd * area_theta_nd / (r_i * dtheta)
            push!(I, idx); push!(J, idx); push!(V, -2.0 * a_t)
            push!(I, idx); push!(J, (ir - 1) * ntheta + itp); push!(V, a_t)
            push!(I, idx); push!(J, (ir - 1) * ntheta + itm); push!(V, a_t)
        end
    end

    MT = spdiagm(0 => M)
    KT = sparse(I, J, V, nnode, nnode)

    return MT, KT, F
end
