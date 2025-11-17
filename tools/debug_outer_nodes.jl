using Statistics
include("../src/JuBat.jl")

param_dim = JuBat.ChooseCell("Jellyroll")
mesh_th = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=160, gsorder=2)
p = JuBat.jellyroll_spiral_params(param_dim)
b = p.b
t_repeat = p.t_repeat
N = p.n_wind

# classify with default tol (1e-4)
is_outer_default = [JuBat.edge_boundary(:node_on, mesh_th, i, param_dim; which=:outer) for i in 1:mesh_th.nlen]
count_default = count(identity, is_outer_default)

# classify with tight tol (1e-8)
is_outer_tight = [JuBat.edge_boundary(:node_on, mesh_th, i, param_dim; which=:outer, tol=1e-8) for i in 1:mesh_th.nlen]
count_tight = count(identity, is_outer_tight)

inner_curve_range = 1:(160+1)                # first nθ+1 nodes are inner curve
outer_curve_range = (160+2):(2*(160+1))      # second half
misclassified_inner_default = count(i -> i in inner_curve_range, findall(is_outer_default))
misclassified_inner_tight   = count(i -> i in inner_curve_range, findall(is_outer_tight))

eps_theta_default = 1e-4 / max(b, 1e-12)
eps_theta_tight   = 1e-8 / max(b, 1e-12)

println("b = ", b)
println("t_repeat = ", t_repeat)
println("n_wind (N) = ", N)
println("Default tol classification count = ", count_default, ", inner misclassified = ", misclassified_inner_default)
println("Tight tol classification count   = ", count_tight, ", inner misclassified = ", misclassified_inner_tight)
println("eps_theta_default = ", eps_theta_default, ", eps_theta_tight = ", eps_theta_tight)

# Show first 10 node cumulative theta values for inner and outer curves for inspection
function cumulative_theta_values(indices)
    vals = Float64[]
    for i in indices
        x = mesh_th.node[i,1]; y = mesh_th.node[i,2]
        r = hypot(x,y)
        push!(vals, (r - p.a - p.t_repeat)/b)  # theta_cum_out formula
    end
    return vals
end
println("Sample theta_cum_out inner[1:5] = ", cumulative_theta_values(1:5))
println("Sample theta_cum_out outer[1:5] (shifted index) = ", cumulative_theta_values( (160+2):(160+6) ))
