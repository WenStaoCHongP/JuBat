using LinearAlgebra, SparseArrays, Statistics, Plots, Printf
root_dir = abspath(joinpath(@__DIR__, ".."))
example_dir = @__DIR__
cd(example_dir)
src_dir = joinpath(root_dir, "src")
include(joinpath(src_dir, "JuBat.jl"))
using .JuBat

function plot_mesh_outline(mesh, out_path)
	x = mesh.node[:, 1]
	y = mesh.node[:, 2]
	plt = plot(size=(700, 700), aspect_ratio=1, legend=false,
		xlabel="x [m]", ylabel="y [m]", title="Jellyroll Mesh Outline")
	ne = size(mesh.element, 1)
	for e in 1:ne
		n = mesh.element[e, :]
		xs = [x[n[1]], x[n[2]], x[n[3]], x[n[4]], x[n[1]]]
		ys = [y[n[1]], y[n[2]], y[n[3]], y[n[4]], y[n[1]]]
		plot!(plt, xs, ys, lw=0.4, color=:black, alpha=0.7)
	end
	savefig(plt, out_path)
end

function mean_temperature(result)
	return result["temperature [K]"]
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

	mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=ntheta, gsorder=2)
	case = JuBat.setup_thermal2D_mesh(case, mesh_data)
	mesh_th = case.mesh["thermal2D"]

	result = JuBat.Solve(case)
	return (; result, mesh=mesh_th)
end

function mesh_convergence(param_dim, t_end)
	ntheta_list = [40, 60, 80, 120]
	@printf("\n[Mesh convergence] 1C, t_end=%.1f s\n", t_end)
	@printf("  ntheta | T_mean_end [K] | V_end [V]\n")
	for nθ in ntheta_list
		run = run_case(param_dim, 1.0, t_end; ntheta=nθ)
		T_mean = mean_temperature(run.result)
		V = run.result["cell voltage [V]"]
		@printf("  %6d | %13.4f | %9.4f\n", nθ, T_mean[end], V[end])
		if nθ == ntheta_list[1]
			out_dir = joinpath(root_dir, "output")
			isdir(out_dir) || mkpath(out_dir)
			plot_mesh_outline(run.mesh, joinpath(out_dir, "spme_thermal_mesh.png"))
		end
	end
end

function main()
	param_dim = JuBat.ChooseCell("Jellyroll")
	mesh_convergence(param_dim, 3600.0)
end

main()
