using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
import PyPlot
root_dir = abspath(joinpath(@__DIR__, ".."))
example_dir = @__DIR__
cd(example_dir)
src_dir = joinpath(root_dir, "src")
include(joinpath(src_dir, "JuBat.jl"))
using .JuBat

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

function tripcolor_field(x, y, tris, values, out_path; title="r-θ FEM Temperature")
	fig = PyPlot.figure(figsize=(6.5, 5.0))
	ax = fig.add_subplot(1, 1, 1)
	tri_mod = PyPlot.pyimport("matplotlib.tri")
	tri = tri_mod.Triangulation(x, y, tris .- 1)
	tpc = ax.tripcolor(tri, values, shading="gouraud", cmap="hot")
	fig.colorbar(tpc, ax=ax)
	ax.set_aspect("equal")
	ax.set_title(title)
	ax.set_xlabel("x [m]")
	ax.set_ylabel("y [m]")
	fig.savefig(out_path, dpi=200)
	PyPlot.close(fig)
end

function run_case(param_dim, crate, t_end; ntheta=80)
	opt = JuBat.Option()
	opt.model = "SPMe"
	opt.dimension = 1
	opt.Nn = 10
	opt.Ns = 5
	opt.Np = 10
	opt.Nrn = 10
	opt.Nrp = 10
	opt.gsorder = 2
	opt.solveType = "Crank-Nicolson"
	opt.dtType = "auto"
	opt.dt = [0.5, 10.0]
	opt.time = [0.0, t_end]

	opt.thermal_enabled = true
	opt.thermalmodel = "distributed2D"
	opt.thermal_dim = "2D"
	opt.cool_method = "surface"

	opt.per_element_spme = true
	opt.mechanicalmodel = "none"
	opt.czm_enabled = false

	I1C = param_dim.cell.I1C
	I = I1C * crate
	opt.Current = x -> I

	case = JuBat.SetCase(param_dim, opt)

	mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=ntheta, gsorder=2)
	case = JuBat.setup_thermal2D_mesh(case, mesh_data)
	mesh_th = case.mesh["thermal2D"]

	result = JuBat.Solve(case)
	return (; result, mesh=mesh_th)
end

function main()
	param_dim = JuBat.ChooseCell("Jellyroll")
	cases = [
		(crate=1.0, t_end=3600.0)
	]

	out_dir = joinpath(root_dir, "output")
	isdir(out_dir) || mkpath(out_dir)

	for c in cases
		run = run_case(param_dim, c.crate, c.t_end; ntheta=80)
		result = run.result
		mesh = run.mesh

		T_nodes_hist = result["thermal2D T_nodes [K]"]
		t = result["time [s]"]
		_, idx = findmin(abs.(t .- 60.0))
		T_nodes_60 = ndims(T_nodes_hist) == 2 ? T_nodes_hist[:, idx] : T_nodes_hist
		out_path = joinpath(out_dir, @sprintf("spme_thermal_T60_%.1fC.png", c.crate))
		x = mesh.node[:, 1]
		y = mesh.node[:, 2]
		tris = quad_to_triangles(mesh.element)
		tripcolor_field(x, y, tris, T_nodes_60, out_path)
	end
end

main()
