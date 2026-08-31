function Arrhenius(Eac::Float64, T::Union{Float64, Array{Float64}})
"""
    Arrhenius function 
        Input: Eac::Float64 -- activation energy (normalised with Eac/(R T0))
                T::Float64 -- temperature (normalised with T/T0)
        Output: the results of Arrhenius function
"""
    return exp.(Eac * (1 .- 1 ./ T))
end

function IntV(x::Array{Float64}, mesh::Mesh)
"""
    Integration function over a domain (mesh) 
        Input: x::Array{Float64} -- the coefficients to integrate
                mesh::Mesh -- mesh structure
        Output: the integration of (x) over "mesh"
"""
    return sum(x .* mesh.gs.weight .* mesh.gs.detJ)
end

"""
    IntQ4(x_e, y_e; order=2) do ξ, η, w, dNdx, dNdy, detJ
        # contribution
    end

Generic Q4 element quadrature driver. The callback receives local coordinates,
integration weight `w`, spatial gradients, and `detJ`.
"""
function IntQ4(f::Function,x_e::AbstractVector{<:Real},y_e::AbstractVector{<:Real};order::Int=2)
    w, q = GSweight(order, 2)
    node_e = hcat(collect(Float64, x_e), collect(Float64, y_e))
    for i in eachindex(w)
        ξ = q[i, 1]
        η = q[i, 2]
        dNdxi = last(LagrangeBasis("Q4", 2, [ξ, η]))
        J = dNdxi * node_e
        detJ = det(J)
        detJ > 1e-15 || error("invalid Q4 Jacobian determinant $detJ at (ξ=$ξ, η=$η)")
        dNdxy = transpose(inv(J) * dNdxi)
        dNdx = @view dNdxy[:, 1]
        dNdy = @view dNdxy[:, 2]
        f(ξ, η, w[i], dNdx, dNdy, detJ)
    end
    return nothing
end

function identify_boundary_nodes(mesh, param, opt=nothing)
    if mesh isa Mesh
        nnode = mesh.nlen
    elseif mesh isa CohesiveMesh
        nnode = mesh.nnode
    else
        error("identify_boundary_nodes: unsupported mesh type $(typeof(mesh))")
    end
    Rin, Rout = param.cell.Rin, param.cell.Rout
    t_repeat = param.PCC.thickness + 2 * (param.PE.thickness + param.SP.thickness + param.NE.thickness) + param.NCC.thickness
    a = Rin
    b = t_repeat / (2 * pi)

    theta0_mesh = max(0.0, (Rin - a) / b)
    theta1_mesh = min((Rout - a - t_repeat) / b, (Rout - a) / b)
    theta_in_range = (theta0_mesh, min(theta0_mesh + 2.0 * pi, theta1_mesh))
    theta_out_range = (max(theta1_mesh - 2.0 * pi, theta0_mesh), theta1_mesh)
    tol = 1e-4

    is_inner = [edge_boundary(mesh, i, param; which=:inner, theta_range=theta_in_range, tol=tol) for i in 1:nnode]
    is_outer = [edge_boundary(mesh, i, param; which=:outer, theta_range=theta_out_range, tol=tol) for i in 1:nnode]
    return is_inner, is_outer
end

"""
    compute_separation(czm_mesh, elem, u)

计算单个 cohesive 单元的平均法/切向分离 `(δ_n, δ_t)`。法向由 `cohesive_local_frame`
按 host-inner → host-outer 拓扑定向，避免逆时针边序把物理张开记为负值。
"""
function compute_separation(
    czm_mesh::CohesiveMesh,
    elem::AbstractCohesiveElement,
    u::Vector{Float64},
)
    _, _, _, R = cohesive_local_frame(czm_mesh, elem)
    n1, n2 = elem.nodes_bottom
    n4, n3 = elem.nodes_top
    wts, pts = NCweight(2)
    delta_n_avg = 0.0
    delta_t_avg = 0.0
    total_w = 0.0

    for (xi, w) in zip(pts, wts)
        N1 = 0.5 * (1.0 - xi)
        N2 = 0.5 * (1.0 + xi)
        B_global = zeros(Float64, 2, 8)
        B_global[1, 1] = -N1; B_global[2, 2] = -N1
        B_global[1, 3] = -N2; B_global[2, 4] = -N2
        B_global[1, 5] = N2;  B_global[2, 6] = N2
        B_global[1, 7] = N1;  B_global[2, 8] = N1

        u_e = [
            u[2n1 - 1], u[2n1], u[2n2 - 1], u[2n2],
            u[2n3 - 1], u[2n3], u[2n4 - 1], u[2n4],
        ]
        delta_local = R * B_global * u_e
        delta_n_avg += w * delta_local[1]
        delta_t_avg += w * delta_local[2]
        total_w += w
    end

    return delta_n_avg / total_w, delta_t_avg / total_w
end

"""
    element_nodal_mean(mesh, nodal_values)

Compute element-wise arithmetic mean from nodal field values.
"""
function element_nodal_mean(mesh::Mesh, nodal_values::AbstractVector{<:Real})
    ne = size(mesh.element, 1)
    out = zeros(Float64, ne)
    @inbounds for e in 1:ne
        nodes = mesh.element[e, :]
        out[e] = sum(@view nodal_values[nodes]) / length(nodes)
    end
    return out
end

function thermal2D_volume_average_temperature(mesh::Mesh, T_nodes::AbstractVector{<:Real})
    ngs = length(mesh.gs.detJ)
    nloc = size(mesh.element, 2)
    num = 0.0
    den = 0.0
    @inbounds for g in 1:ngs
        e = mesh.gs.ele[g]
        wJ = mesh.gs.weight[g] * mesh.gs.detJ[g]
        nodes = mesh.element[e, :]
        for i in 1:nloc
            wi = wJ * mesh.gs.Ni[g, i]
            num += wi * T_nodes[nodes[i]]
            den += wi
        end
    end
    return num / den
end

"""
    q4_center_gradients(node, elem_nodes)

Return Q4 shape-function spatial gradients at element center (xi=0, eta=0).

# Returns
- `dNdx`: length-4 vector
- `dNdy`: length-4 vector
- `detJ`: Jacobian determinant at center
"""
function q4_center_gradients(node::AbstractMatrix{<:Real}, elem_nodes::AbstractVector{<:Integer})
    dNdxi = last(LagrangeBasis("Q4", 2, [0.0, 0.0]))
    x_nodes = collect(Float64, node[elem_nodes, 1])
    y_nodes = collect(Float64, node[elem_nodes, 2])
    node_e = hcat(x_nodes, y_nodes)
    J = dNdxi * node_e
    detJ = det(J)
    detJ > 1e-12 || error("invalid Q4 center Jacobian determinant $detJ")
    dNdxy = transpose(inv(J) * dNdxi)
    dNdx = @view dNdxy[:, 1]
    dNdy = @view dNdxy[:, 2]
    return dNdx, dNdy, detJ
end
