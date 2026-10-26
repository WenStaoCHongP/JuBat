"""Validate the native radial-major, uniformly spaced periodic ring layout."""
function _polar_ring_grid(mesh::Mesh, mesh_data, param::Params)
    bad(message) = throw(ArgumentError("ThermalPolar2D_Ring: " * message))
    r, nt = mesh_data.r, mesh_data.ntheta
    r isa AbstractVector{<:Real} && length(r) >= 2 && all(isfinite, r) &&
        all(>(0), r) && all(>(0), diff(r)) || bad("r must be finite, positive and strictly increasing")
    nr = length(r) - 1
    tol = 128 * eps(Float64)
    close(a, b) = isapprox(a, b; rtol=tol, atol=tol * maximum(r))
    close(first(r), param.cell.Rin) && close(last(r), param.cell.Rout) ||
        bad("normalized mesh radii must match cell.Rin/Rout")
    phase = atan(mesh.node[1, 2], mesh.node[1, 1])
    dtheta = 2pi / nt
    for ir in eachindex(r), it in 1:nt
        n = (ir - 1) * nt + it
        angle = phase + (it - 1) * dtheta
        close(mesh.node[n, 1], r[ir] * cos(angle)) &&
            close(mesh.node[n, 2], r[ir] * sin(angle)) ||
            bad("nodes must follow the same uniform angular order on every radius")
    end
    for ir in 1:nr, it in 1:nt
        itp = it == nt ? 1 : it + 1
        e = (ir - 1) * nt + it
        expected = ((ir - 1) * nt + it, ir * nt + it, ir * nt + itp, (ir - 1) * nt + itp)
        all(mesh.element[e, j] == expected[j] for j in 1:4) ||
            bad("element numbering must match the native periodic ring connectivity")
    end
    return r, nt
end

"""
    ThermalPolar2D_Ring(case, variables, mesh_data)

Assemble the native annular-sector FVM. `heat_source_fields[e]` is the signed
average dimensionless volumetric source on curved sector e, with area
`Asector[e] = (r_outer^2-r_inner^2)*pi/ntheta`. Mesh coordinates and radii are
normalized, e.g. from `ring_mesh(case.param)`.

A source from the same-numbered straight Q4 ring has a different integration
domain. Its caller must first convert EACH cell's power:
`qsector = qQ4 .* AQ4 ./ Asector`, using the actual Q4 areas. A global power
rescaling, or treating a Jellyroll partition as this ring, is not a valid bridge.
"""
function ThermalPolar2D_Ring(case::Case, variables::AbstractDict{String,<:Any}, mesh_data)
    mesh = case.mesh["thermal2D"]
    param = case.param  # Use normalized parameters
    scale = case.param.scale

    r_nodes, ntheta = _polar_ring_grid(mesh, mesh_data, param)

    # Grid spacing (in radians)
    dtheta = 2.0 * pi / ntheta

    # --- Normalized thermal properties (from case.param.cell) ---
    # Volumetric heat capacity: (ρc)* = C* / V* = heat_Q / volume
    # where C* is total heat capacity and V* is normalized volume
    V_nd = param.cell.volume  # V* = V / L³
    C_nd = param.cell.heat_Q  # C* = m·c·T_ref/(t0·P_ref)
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

    # Each sector contributes its four true overlaps with nodal control volumes.
    # Column sums equal sector areas; row sums equal the control areas below.
    q_elem = variables["heat_source_fields"]
    for ir in 1:(length(r_nodes) - 1)
        rmid = 0.5 * (r_nodes[ir] + r_nodes[ir + 1])
        win = 0.25 * (rmid^2 - r_nodes[ir]^2) * dtheta
        wout = 0.25 * (r_nodes[ir + 1]^2 - rmid^2) * dtheta
        for it in 1:ntheta
            e = (ir - 1) * ntheta + it
            n1, n2, n3, n4 = mesh.element[e, :]
            F[n1] += win * q_elem[e]
            F[n4] += win * q_elem[e]
            F[n2] += wout * q_elem[e]
            F[n3] += wout * q_elem[e]
        end
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

            # --- Radial conduction (inner neighbor) ---
            if ir > 1
                # Conductance: a = k* · r_face · Δθ / Δr
                a_rm = k_r_nd * r_imh * dtheta / dr_im
                push!(I, idx); push!(J, idx); push!(V, -a_rm)
                idxm = (ir - 2) * ntheta + it
                push!(I, idx); push!(J, idxm); push!(V, a_rm)
            end

            # --- Radial conduction (outer neighbor) ---
            if ir < length(r_nodes)
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
