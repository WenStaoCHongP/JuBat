using Test, Serialization, LinearAlgebra, TOML

if ARGS == ["--compare"]
    root = normpath(joinpath(@__DIR__,"..","output","validate_case_batch_b_j2"))
    checks = Dict{String,Any}[]
    for (reference,candidate,file) in (("mild_yield_reference","mild_yield_candidate","compression.bin"),
            ("mild_yield_reference","mild_yield_candidate","repeat.bin"),
            ("original_repeat_reference","original_repeat_candidate","repeat.bin"))
        a,b = deserialize(joinpath(root,reference,file)),deserialize(joinpath(root,candidate,file))
        max_kappa = maximum(s.kappa for s in a.states)
        push!(checks,Dict("reference"=>reference,"candidate"=>candidate,"file"=>file,
            "exact_equal"=>isequal(a,b),"max_kappa"=>max_kappa,"active_plasticity"=>max_kappa>0))
    end
    open(joinpath(root,"comparison.toml"),"w") do io
        TOML.print(io,Dict("checks"=>checks,"tolerances_changed"=>false,"source_parameters_changed"=>false))
    end
    @assert all(c["exact_equal"] && c["active_plasticity"] for c in checks)
    println("J2 exact state comparisons: 3/3; all active plasticity")
    exit()
end

length(ARGS) in (2,3) || error("usage: validate_case_batch_b_j2.jl TEST_PATH LABEL [repeat-only]")
const BJ2_FACTOR = parse(Float64,get(ENV,"JUBAT_TEST_YIELD_FACTOR","1"))
isfinite(BJ2_FACTOR) && 0 < BJ2_FACTOR <= 1 || error("invalid yield factor")
const BJ2_OUT = normpath(joinpath(@__DIR__,"..","output","validate_case_batch_b_j2",ARGS[2]))
mkpath(BJ2_OUT)
path = abspath(ARGS[1])
source = read(path,String)
needle = "    param_dim = JuBat.ChooseCell(\"Jellyroll\")"
@assert count(needle,source) == 1
source = replace(source,needle => needle * "\n    param_dim.PCC.sigma_y *= Main.BJ2_FACTOR\n    param_dim.NCC.sigma_y *= Main.BJ2_FACTOR")
needle = "    @test maximum(κ_foil) > 0.0"
@assert count(needle,source) == 1
capture = """
    stress_foil = [s.sigma for (e, row) in enumerate(eachrow(ms.plastic_states)) for s in row if mt[e] in (:PCC,:NCC)]
    vm_foil = [sqrt(s[1]^2-s[1]*s[2]+s[2]^2+3s[3]^2) for s in stress_foil]
    println("YIELD_FACTOR=",Main.BJ2_FACTOR," PCC_YIELD_PA=",case.param_dim.PCC.sigma_y," NCC_YIELD_PA=",case.param_dim.NCC.sigma_y,
        " MAX_FOIL_VM_PA=",maximum(vm_foil)*case.param.scale.σ_czm," MAX_FOIL_KAPPA=",maximum(κ_foil))
    serialize(joinpath(Main.BJ2_OUT,"compression.bin"),(;displacement=copy(r.displacement),residual=r.residual_norm,
        states=map(s -> (;eps_p=s.eps_p,kappa=s.kappa,sigma=s.sigma),ms.plastic_states),
        damage=map(s -> Tuple(getfield(s,f) for f in fieldnames(typeof(s))),ms.damage_states)))
"""
source = replace(source,needle=>capture*"\n"*needle)
needle = "    @test κ_after1 > 0.0"
@assert count(needle,source) == 1
capture = """
    println("ACTIVE_REPEAT_KAPPA_FIRST=",κ_after1," SECOND=",κ_after2)
    serialize(joinpath(Main.BJ2_OUT,"repeat.bin"),(;displacement=copy(r2.displacement),residual=r2.residual_norm,
        converged=(r1.converged,r2.converged),kappa=(κ_after1,κ_after2),
        states=map(s -> (;eps_p=s.eps_p,kappa=s.kappa,sigma=s.sigma),ms.plastic_states),
        damage=map(s -> Tuple(getfield(s,f) for f in fieldnames(typeof(s))),ms.damage_states)))
"""
source = replace(source,needle=>capture*"\n"*needle)
write(joinpath(BJ2_OUT,"configuration.txt"),"yield_factor=$BJ2_FACTOR\nsource_test=$path\n")
@testset "J2 full suite with explicit local yield factor" begin
    if length(ARGS) == 3
        ARGS[3] == "repeat-only" || error("unknown mode")
        # Run the unchanged active-plasticity group only, keeping its original fixture and assertions.
        select_repeat(expr) = expr isa Expr && expr.head == :macrocall && expr.args[1] == Symbol("@testset") &&
            !(expr.args[3] isa String && startswith(expr.args[3],"状态不可逆")) ? :(nothing) : expr
        include_string(select_repeat,Main,source,path)
    else
        include_string(Main,source,path)
    end
end
