using Test, LinearAlgebra, Serialization
include(joinpath(@__DIR__, "..", "src", "JuBat.jl"))
using .JuBat
const OUT_B = normpath(joinpath(@__DIR__, "..", "output", "test_case_batch_b"))
mkpath(OUT_B)

p = JuBat.SetCase(JuBat.ChooseCell("Jellyroll"),JuBat.Option()).param
mesh, meta = JuBat.create_unit_czm_strip(p)
workspace = JuBat.CZMWorkspace(mesh,p)
kw = (;czm_workspace=workspace)
ndof = 2size(mesh.node,1)
u = collect(range(0.0,1e-6;length=ndof))
eig = (dT=fill(0.001,size(mesh.bulk_element,1)),Δsn=zeros(size(mesh.bulk_element,1)),Δsp=zeros(size(mesh.bulk_element,1)))
K = JuBat.bulk_stiffness(mesh,p;kw...)
f, tangent = JuBat.assemble_bulk_residual_tangent(mesh,u,p;geo_nl=true,eigenstrain=eig,kw...)
state = JuBat.MechState(mesh)
coupled = JuBat.assemble_coupled_system(mesh,u,p;damage_states=state.damage_states,geo_nl=true,eigenstrain=eig,kw...)
mapping = [mesh.czm_submesh.thermal_elem_map[c.host_outer_elem] for c in mesh.cohesive_elements]
numeric = (;K=Matrix(K),f,tangent=Matrix(tangent),coupled_K=Matrix(coupled[1]),coupled_f=coupled[2],
    separations=coupled[3],tractions=coupled[4],mapping,node=mesh.node,bulk_element=mesh.bulk_element)
@testset "B numeric identity with A" begin
    ref = deserialize(joinpath(OUT_B,"reference.bin"))
    for key in keys(numeric)
        @test getproperty(numeric,key) == getproperty(ref,key)
    end
end

@testset "Boundary diagnostic helper receives its session workspace" begin
    # Load only this pure helper's actual definition; the full diagnostic is a long-cycle launcher.
    tree = Meta.parseall(read(joinpath(@__DIR__,"..","example","check_mechanical_boundary_branches.jl"),String))
    definitions = filter(x -> x isa Expr && x.head == :function &&
        x.args[1] isa Expr && x.args[1].head == :call && x.args[1].args[1] == :compact_GP_metrics, tree.args)
    @test length(definitions) == 1
    Core.eval(Main,only(definitions))
    state.gp_damage_states = JuBat.gp_history_for_step(state,mesh)
    metrics = compact_GP_metrics(mesh,state,[0.0,1.0];czm_workspace=workspace)
    @test metrics.d_max == 0.0
    @test metrics.mean_length_weighted == 0.0
end

@testset "B session workspaces, topology and physical-state isolation" begin
    @test fieldnames(JuBat.CohesiveMesh) == (:node,:bulk_element,:cohesive_elements,:czm_submesh)
    @test !hasfield(JuBat.Case,:czm_mesh)
    other = JuBat.CZMWorkspace(mesh,p)
    @test JuBat.bulk_stiffness(mesh,p;czm_workspace=other) == K
    @test other.K_bulk !== workspace.K_bulk
    @test JuBat.cohesive_geometry(mesh;czm_workspace=other) !== JuBat.cohesive_geometry(mesh;kw...)
    @test JuBat.assembly_workspace(mesh;czm_workspace=other) !== JuBat.assembly_workspace(mesh;kw...)
    @test JuBat.bulk_gauss_geometry(mesh;czm_workspace=other) !== JuBat.bulk_gauss_geometry(mesh;kw...)
    JuBat.assemble_bulk_residual_tangent(mesh,u,p;geo_nl=true,eigenstrain=eig,czm_workspace=other)
    @test other.bulk_asm !== workspace.bulk_asm
    saved_a = deepcopy((workspace.bulk_asm.ke,workspace.bulk_asm.fe,workspace.bulk_asm.nzval))
    saved_mech = deepcopy(state)
    JuBat.assemble_coupled_system(mesh,2u,p;damage_states=state.damage_states,geo_nl=true,
        eigenstrain=eig,czm_workspace=other)
    @test (workspace.bulk_asm.ke,workspace.bulk_asm.fe,workspace.bulk_asm.nzval) == saved_a
    @test state.u_prev == saved_mech.u_prev
    @test [s.D for s in state.damage_states] == [s.D for s in saved_mech.damage_states]
    @test other.bulk_asm.ke !== workspace.bulk_asm.ke
    @test other.bulk_asm.fe !== workspace.bulk_asm.fe
    @test other.bulk_asm.nzval !== workspace.bulk_asm.nzval
    # Reuse A after B overwrote its own buffers; all A results remain exact.
    again = JuBat.assemble_coupled_system(mesh,u,p;damage_states=state.damage_states,geo_nl=true,eigenstrain=eig,kw...)
    @test Matrix(again[1]) == numeric.coupled_K
    @test again[2] == numeric.coupled_f
    @test_throws ArgumentError JuBat.bulk_stiffness(deepcopy(mesh),p;kw...)
    opt = JuBat.CzmOptions()
    @test_throws ArgumentError JuBat.solve_czm_step(deepcopy(mesh),state,p,zeros(ndof),opt;kw...)
    @test_throws ArgumentError JuBat.solve_czm_step(mesh,state,deepcopy(p),zeros(ndof),opt;kw...)

    strain = zeros(size(mesh.bulk_element,1))
    stress_a = JuBat.recover_bulk_stress_geo(mesh,u,p,strain;kw...)
    stress_b = JuBat.recover_bulk_stress_geo(mesh,u,p,strain;czm_workspace=other)
    @test stress_a == stress_b
    @test all(isfinite,vcat(stress_a...))

    hist = Dict{String,Union{Array{Float64},Float64}}("force"=>zeros(ndof,1))
    vars = Dict{String,Union{Array{Float64},Float64}}("force"=>copy(numeric.f))
    JuBat.Variable_update!(hist,vars,1)
    saved_hist = copy(hist["force"])
    fill!(other.bulk_asm.fe,123.0)
    vars["force"] .= 456.0
    @test hist["force"] == saved_hist

    path = joinpath(OUT_B,"snapshot.bin")
    opt.iter_method = "arc_length"
    JuBat.save_czm_step_snapshot(path,mesh,deepcopy(state),p,zeros(ndof),opt; dT_elem=zeros(8),Δsoc_n_elem=zeros(8),Δsoc_p_elem=zeros(8))
    payload = deserialize(path)
    @test fieldnames(typeof(payload.czm_mesh)) == (:node,:bulk_element,:cohesive_elements,:czm_submesh)
    @test !hasproperty(payload,:czm_workspace)
    @test payload.mech.u_prev !== state.u_prev
    @test payload.mech.damage_states !== state.damage_states
    @test payload.czm_mesh.node !== mesh.node
    replayed, replay_ms = JuBat.replay_czm_step_snapshot(path)
    @test replayed.converged
    @test replay_ms.u_prev !== payload.mech.u_prev
    @test payload.mech.u_prev == saved_mech.u_prev
    @test state.u_prev == saved_mech.u_prev

    damage = [JuBat.DamageState() for _ in mesh.cohesive_elements]
    damage[2].D = 0.4
    @test JuBat.map_czm_damage_to_thermal(mesh,damage,1) == [0.4]
    @test_throws ErrorException JuBat.map_czm_damage_to_thermal(mesh,damage[1:end-1],1)
    damage[2].D = NaN
    @test_throws ErrorException JuBat.map_czm_damage_to_thermal(mesh,damage,1)
    damage[2].D = 0.4
    bad = deepcopy(mesh); bad.czm_submesh.thermal_elem_map .= 2
    @test_throws ErrorException JuBat.map_czm_damage_to_thermal(bad,damage,1)
    bad.czm_submesh = nothing
    @test_throws ErrorException JuBat.map_czm_damage_to_thermal(bad,damage,1)
end
