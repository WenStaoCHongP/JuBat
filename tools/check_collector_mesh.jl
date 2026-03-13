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
areas, lw = JuBat.jellyroll_element_properties(mesh, param)
println("layer_weights shape: ", size(lw))
println("sample row 1: ", lw[1,:])
println("areas: min=$(minimum(areas)), max=$(maximum(areas)), sum=$(sum(areas))")
println("first 6 nodes:")
for i in 1:min(6,size(mesh.node,1))
  println(i, " ", mesh.node[i,:])
end
println("first 6 elements:")
for i in 1:min(6,size(mesh.element,1))
  println(i, " ", mesh.element[i,:])
end

# print outer spiral r_out(θ) = a + b θ + s_out endpoint angle and arc lengths
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
  s_in = 0.0
  s_out = t_repeat

  # Arc length integral for Archimedean spiral: r(θ) = a + b*θ + s
  # L = ∫₀^θ sqrt(r² + (dr/dθ)²) dθ = ∫₀^θ sqrt((a + b*θ + s)² + b²) dθ
  # Using substitution u = a + b*θ + s:
  # F(u) = (u * sqrt(u² + b²) + b² * asinh(u/b)) / (2b)
  function spiral_arc_length(u::Float64, b::Float64)
    if b <= 0
      return 0.0
    end
    return (u * sqrt(u^2 + b^2) + b^2 * asinh(u / b)) / (2 * b)
  end

  # Inner spiral: s_in = 0, r_in(θ) = a + b*θ
  # Outer spiral: s_out = t_repeat, r_out(θ) = a + b*θ + t_repeat

  # Truncation logic (consistent with Jellyrollmodel.jl):
  # theta1 = min((Rout - a - s_out) / b, (Rout - a) / b)
  # Since s_out = t_repeat > 0, theta1 = (Rout - a - s_out) / b (outer spiral endpoint)
  # Both spirals use the SAME angle range [0, theta1]
  theta1 = min((Rout - a - s_out) / b, (Rout - a) / b)

  # Calculate arc lengths using the actual theta range [0, theta1]
  # At theta1:
  #   - r_in(theta1) = a + b*theta1 = Rout - t_repeat
  #   - r_out(theta1) = a + b*theta1 + t_repeat = Rout

  # Inner spiral: from u=a to u=a + b*theta1 = Rout - t_repeat
  u_start_in = a + s_in
  u_end_in = a + b * theta1 + s_in
  L_in = spiral_arc_length(u_end_in, b) - spiral_arc_length(u_start_in, b)

  # Outer spiral: from u=a+t_repeat to u=Rout
  u_start_out = a + s_out
  u_end_out = a + b * theta1 + s_out
  L_out = spiral_arc_length(u_end_out, b) - spiral_arc_length(u_start_out, b)

  println("\n=== Spiral Arc Lengths ===")
  println("Grid theta range: [0, $theta1] rad = [0, $(theta1 * 180.0 / pi)] deg")
  println("At theta1: r_in = $(a + b * theta1) m, r_out = $(a + b * theta1 + s_out) m")
  println("Inner spiral arc length: $L_in m")
  println("Outer spiral arc length: $L_out m")
  println("Total spiral length (inner + outer): $(L_in + L_out) m")
end
