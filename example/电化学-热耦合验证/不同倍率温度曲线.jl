using LinearAlgebra, SparseArrays, Statistics, Plots, Printf, CSV
root_dir = abspath(joinpath(@__DIR__, "..", ".."))
example_dir = @__DIR__
cd(example_dir)
src_dir = joinpath(root_dir, "src")
include(joinpath(src_dir, "JuBat.jl"))
using .JuBat

function read_pybamm_curve(path)
	rows = CSV.File(path)
	cap = Float64[]
	volt = Float64[]
	temp = Float64[]
	for row in rows
		push!(cap, row[Symbol("capacity [A.h]")])
		push!(volt, row[Symbol("voltage [V]")])
		push!(temp, row[Symbol("temperature [K]")])
	end
	return cap, volt, temp
end

function read_sim_curve(path)
	rows = CSV.File(path)
	cap = Float64[]
	volt = Float64[]
	temp = Float64[]
	temp_max = Float64[]
	has_temp_max = Symbol("temperature max [K]") in propertynames(rows)
	for row in rows
		push!(cap, row[Symbol("capacity [A.h]")])
		push!(volt, row[Symbol("voltage [V]")])
		push!(temp, row[Symbol("temperature [K]")])
		if has_temp_max
			push!(temp_max, row[Symbol("temperature max [K]")])
		else
			push!(temp_max, row[Symbol("temperature [K]")])
		end
	end
	return cap, volt, temp, temp_max
end

function load_or_run_sim(param_dim, crate, t_end, sim_path)
	if isfile(sim_path)
		cap, volt, temp, temp_max = read_sim_curve(sim_path)
		if length(cap) > 1 && length(volt) > 1 && length(temp) > 1
			return cap, volt, temp, temp_max
		end
		rm(sim_path; force=true)
	end

	run = run_case(param_dim, crate, t_end; ntheta=80)
	I = abs(run.I)
	capacity_Ah = I .* run.time ./ 3600.0
	write_sim_curve(sim_path, capacity_Ah, run.voltage, run.temp, run.temp_max)
	return capacity_Ah, run.voltage, run.temp, run.temp_max
end

function write_sim_curve(path, cap, volt, temp, temp_max)
	mkpath(dirname(path))
	n = minimum(length.([cap, volt, temp, temp_max]))
	open(path, "w") do io
		println(io, "capacity [A.h],voltage [V],temperature [K],temperature max [K]")
		for i in 1:n
			println(io, cap[i], ",", volt[i], ",", temp[i], ",", temp_max[i])
		end
	end
end

function interp1(x, y, xi)
	j = searchsortedfirst(x, xi)
	if j <= 1
		return y[1]
	elseif j > length(x)
		return y[end]
	end
	x0 = x[j - 1]
	x1 = x[j]
	y0 = y[j - 1]
	y1 = y[j]
	if abs(x1 - x0) < 1e-12
		return y0
	end
	t = (xi - x0) / (x1 - x0)
	return y0 + t * (y1 - y0)
end

function sample_by_capacity(cap, y; n=100)
	order = sortperm(cap)
	cap_sorted = cap[order]
	y_sorted = y[order]
	cmin, cmax = extrema(cap_sorted)
	cap_samples = collect(range(cmin, cmax; length=n))
	y_samples = [interp1(cap_sorted, y_sorted, c) for c in cap_samples]
	return cap_samples, y_samples
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

	initial_state = Dict{String, Any}()
	_, timestep_data = JuBat.solve_phase_with_export(
		case, JuBat.PHASE_DISCHARGE, t_end, I, param_dim.cell.v_l, initial_state;
		dt_range=opt.dt,
		export_interval=1
	)

	times = [ts.time for ts in timestep_data]
	voltage = [ts.V for ts in timestep_data]
	temp = [ts.T_mean for ts in timestep_data]
	temp_max = [ts.T_max for ts in timestep_data]

	# Add explicit initial state to avoid misreading first computed step as 0Ah.
	T0_raw = param_dim.cell.T0
	T0_K = T0_raw < 10.0 ? T0_raw * param_dim.scale.T_ref : T0_raw
	if !isempty(times)
		pushfirst!(times, 0.0)
		pushfirst!(voltage, voltage[1])
		pushfirst!(temp, T0_K)
		pushfirst!(temp_max, T0_K)
	end

	return (; time=times, voltage=voltage, temp=temp, temp_max=temp_max, I=I, mesh=mesh_th)
end

function main()
	param_dim = JuBat.ChooseCell("Jellyroll")
	cases = [
		(crate=0.3, t_end=12000.0),
		(crate=0.7, t_end=4800.0),
		(crate=1.0, t_end=3600.0)
	]

	out_dir = joinpath(root_dir, "output")
	csv_dir = joinpath(out_dir, "csv")
	isdir(out_dir) || mkpath(out_dir)

	sim_curves = Dict{Float64, NamedTuple{(:cap, :volt, :temp, :temp_max), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()
	for c in cases
		sim_path = joinpath(csv_dir, @sprintf("spme_thermal_sim_%.1fC.csv", c.crate))
		cap, volt, temp, temp_max = load_or_run_sim(param_dim, c.crate, c.t_end, sim_path)
		sim_curves[c.crate] = (cap=cap, volt=volt, temp=temp, temp_max=temp_max)
	end

	data_curves = Dict{Float64, NamedTuple{(:cap, :volt, :temp), Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()
	data_dir = joinpath(root_dir, "src", "data")
	for c in cases
		path = joinpath(data_dir, @sprintf("pybamm_SPMe_LGM50_%.1fC.csv", c.crate))
		cap, volt, temp = read_pybamm_curve(path)
		data_curves[c.crate] = (cap=cap, volt=volt, temp=temp)
	end

	crates_sorted = sort(collect(keys(sim_curves)))
	colors = [:blue, :red, :green, :orange, :purple]
	markers = [:circle, :square, :diamond, :utriangle, :pentagon]

	p_v = plot(xlabel="Capacity [Ah]", ylabel="Voltage [V]",
		title="Voltage vs Capacity", legend=:topright)
	p_t = plot(xlabel="Capacity [Ah]", ylabel="Temperature [K]",
		title="Temperature vs Capacity", legend=:topleft)

	for (i, c) in enumerate(crates_sorted)
		sim = sim_curves[c]
		data = data_curves[c]

		plot!(p_v, data.cap, data.volt, lw=2, color=colors[i], label=@sprintf("Data %.1fC", c))
		plot!(p_t, data.cap, data.temp, lw=2, color=colors[i], label=@sprintf("Data %.1fC", c))

		cap_s, volt_s = sample_by_capacity(sim.cap, sim.volt; n=30)
		cap_t, temp_s = sample_by_capacity(sim.cap, sim.temp; n=30)
		scatter!(p_v, cap_s, volt_s, color=colors[i], marker=markers[i],
			ms=4, msw=0.6, label=@sprintf("Sim %.1fC", c))
		scatter!(p_t, cap_t, temp_s, color=colors[i], marker=markers[i],
			ms=4, msw=0.6, label=@sprintf("Sim %.1fC", c))
	end

	savefig(p_v, joinpath(out_dir, "spme_thermal_voltage_capacity_compare.png"))
	savefig(p_t, joinpath(out_dir, "spme_thermal_temperature_capacity_compare.png"))
end

main()
