using Statistics
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
param=JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.czm_enabled = false
case = JuBat.SetCase(param, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param; nθ=10, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
mesh = case.mesh["thermal2D"]
ne=size(mesh.element,1)
println("ne=$ne, nnode=$(mesh.nlen)")
ngs = length(mesh.gs.detJ)
A=zeros(Float64,ne)
for g=1:ngs
  A[mesh.gs.ele[g]] += mesh.gs.weight[g]*mesh.gs.detJ[g]
end
println("area: min=$(minimum(A)), max=$(maximum(A)), mean=$(mean(A))")
println("detJ stats: min=$(minimum(mesh.gs.detJ)), max=$(maximum(mesh.gs.detJ))")
lw = JuBat.jellyroll_get_layer_weights(mesh)
println("layer_weights bound=$(lw !== nothing)")
if lw !== nothing
  println("layer_weights shape: ", size(lw))
  println("sample row 1: ", lw[1,:])
end
println("first 6 nodes:")
for i in 1:min(6,size(mesh.node,1))
  println(i, " ", mesh.node[i,:])
end
println("first 6 elements:")
for i in 1:min(6,size(mesh.element,1))
  println(i, " ", mesh.element[i,:])
end

# print outer spiral r_out(θ) = a + b θ + s_out endpoint angle
begin
  cell = param.cell
  Rin, Rout = cell.Rin, cell.Rout
  t_PE = param.PE.thickness
  t_NE = param.NE.thickness
  t_SP = param.SP.thickness
  t_PCC = param.PCC.thickness
  t_NCC = param.NCC.thickness
  t_repeat = t_PCC + 2 * (t_PE + t_SP + t_NE) + t_NCC
  a = Rin
  b = t_repeat / (2 * pi)
  s_out = t_repeat
  # solve for theta_end from r_out(theta_end) = Rout => a + b*theta_end + s_out = Rout
  theta_end = (Rout - a - s_out) / b
  println("outer spiral endpoint theta (rad) = ", theta_end)
  println("outer spiral endpoint theta (deg) = ", theta_end * 180.0 / pi)
end
