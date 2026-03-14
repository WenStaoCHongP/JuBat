using LinearAlgebra

if !isdefined(@__MODULE__, :ShapeFunction2D)
    include(joinpath(@__DIR__, "SetMesh.jl"))
end

"""
    ring_mesh(param; ntheta, dr, nr, phase, gsorder)

Build a Q4 ring mesh using equal-angle segmentation and radial seeding.
Returns a NamedTuple with mesh and boundary node sets.

# Arguments
- `param`: Normalized parameter set (`case.param`) that provides
           `param.cell.Rin` and `param.cell.Rout` in nondimensional form
- `ntheta`: Number of angular divisions
- `dr`: Nondimensional radial step size (used if nr is nothing)
- `nr`: Number of radial divisions (overrides dr)
- `phase`: Phase offset for theta [radians]
- `gsorder`: Gauss integration order

# Returns
- `mesh`: Mesh struct with nondimensional node coordinates
- `inner_nodes`, `outer_nodes`: Boundary node indices
- `r`: Nondimensional radial coordinates
- `theta`, `nr`, `ntheta`: Mesh metadata
"""
function ring_mesh(param;
    ntheta::Int=40,
    dr::Union{Nothing,Float64}=nothing,
    nr::Union{Nothing,Int}=nothing,
    phase::Float64=0.0,
    gsorder::Int=2)

    Rin = param.cell.Rin
    Rout = param.cell.Rout

    (Rout > Rin) || error("Rout must be greater than Rin")
    (ntheta >= 3) || error("ntheta must be >= 3")

    dr_eff = dr === nothing ? (Rout - Rin) / 20 : dr
    (dr_eff > 0.0) || error("dr must be positive")

    nr_eff = nr === nothing ? max(1, round(Int, (Rout - Rin) / dr_eff)) : max(1, nr)
    r = collect(range(Rin, Rout; length=nr_eff + 1))
    theta = phase .+ collect(range(0.0, 2.0 * pi; length=ntheta + 1))[1:ntheta]

    nnode = (nr_eff + 1) * ntheta
    node = zeros(Float64, nnode, 2)

    idx(ir, it) = (ir - 1) * ntheta + it

    for ir in 1:(nr_eff + 1)
        for it in 1:ntheta
            n = idx(ir, it)
            node[n, 1] = r[ir] * cos(theta[it])
            node[n, 2] = r[ir] * sin(theta[it])
        end
    end

    ne = nr_eff * ntheta
    element = zeros(Int64, ne, 4)
    e = 0
    for ir in 1:nr_eff
        for it in 1:ntheta
            it_next = it == ntheta ? 1 : it + 1
            e += 1
            n1 = idx(ir, it)
            n2 = idx(ir + 1, it)
            n3 = idx(ir + 1, it_next)
            n4 = idx(ir, it_next)
            element[e, 1] = n1
            element[e, 2] = n2
            element[e, 3] = n3
            element[e, 4] = n4
        end
    end

    # Ensure consistent counterclockwise orientation to keep detJ positive
    for e in 1:ne
        n1, n2, n3, n4 = element[e, 1], element[e, 2], element[e, 3], element[e, 4]
        x1, y1 = node[n1, 1], node[n1, 2]
        x2, y2 = node[n2, 1], node[n2, 2]
        x3, y3 = node[n3, 1], node[n3, 2]
        x4, y4 = node[n4, 1], node[n4, 2]
        area = 0.5 * ((x1*y2 - x2*y1) + (x2*y3 - x3*y2) + (x3*y4 - x4*y3) + (x4*y1 - x1*y4))
        if area < 0.0
            element[e, 2], element[e, 4] = element[e, 4], element[e, 2]
        end
    end

    gs = GetGS(element, node, gsorder, "Q4")
    mesh = Mesh("Q4", 2, node, nnode, element, gs)

    inner_nodes = [idx(1, it) for it in 1:ntheta]
    outer_nodes = [idx(nr_eff + 1, it) for it in 1:ntheta]

    return (
        mesh = mesh,
        inner_nodes = inner_nodes,
        outer_nodes = outer_nodes,
        r = r,
        theta = vcat(theta, theta[1] + 2.0 * pi),
        nr = nr_eff,
        ntheta = ntheta
    )
end
