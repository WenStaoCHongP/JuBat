using Test, LinearAlgebra, SparseArrays, DelimitedFiles, Serialization, SHA

# Production entry points; the frozen source override is for consumer A/B only.
const ROOT66 = normpath(joinpath(@__DIR__, ".."))
const SOURCE66 = get(ENV, "TASK66_SOURCE", joinpath(ROOT66, "src"))
include(joinpath(SOURCE66, "JuBat.jl"))
using .JuBat
const RUN66 = get(ENV, "TASK66_RUN", "green")
const OUT66 = joinpath(ROOT66, "output", "test_thermal_fem_rotation", RUN66)
mkpath(OUT66)
if "--capture-baseline" in ARGS
    # Execute the exact production example configuration, adding a read-only result sink.
    example_path=joinpath(ROOT66,"example","testexample.jl")
    original=read(example_path,String)
    loading="include(joinpath(@__DIR__, \"../src/JuBat.jl\"))"
    needle="result = JuBat.Solve(case)"
    @assert count(loading,original)==1 && count(needle,original)==1
    capture="chemistry66 = NamedTuple[]\n    result = JuBat.Solve(case; return_final_state=true, chemical_step_observer=ev -> begin\n        ev.kind == :accepted_commit || return\n        push!(chemistry66, (;t_eval=ev.t_eval, y_before=copy(ev.y_before), y_after=copy(ev.y_after)))\n    end)\n    serialize(joinpath(OUT66, \"baseline-result.bin\"), (;result, chemistry=chemistry66, final=deepcopy(case.mech)))\n    serialize(joinpath(OUT66, \"baseline-geometry.bin\"), (thermal_nodes=case.mesh[\"thermal2D\"].node, thermal_elements=case.mesh[\"thermal2D\"].element, mechanical_nodes=case.mesh[\"czm\"].node, mechanical_elements=case.mesh[\"czm\"].bulk_element))"
    instrumented=replace(replace(original,loading=>""),needle=>capture)
    write(joinpath(OUT66,"baseline-original-sha256.txt"),bytes2hex(sha256(original)))
    write(joinpath(OUT66,"baseline-instrumented.jl"),instrumented)
    include_string(Main,instrumented,example_path)
    exit()
end
if "--analyze-baseline" in ARGS
    prior=deserialize(joinpath(ROOT66,"output","task65_scientific_review","after-production.bin"))
    current=deserialize(joinpath(ROOT66,"output","test_thermal_fem_rotation","baseline-after-capture","baseline-result.bin"))
    b,a=prior.result,current.result
    rows=Any[]; finite_rows=Any[]; nonfinite_locations=Any[]
    for key in sort(collect(intersect(keys(a),keys(b))))
        startswith(key,"timing ") && continue
        aa,bb=a[key],b[key]
        if aa isa AbstractArray{<:Number} && bb isa AbstractArray{<:Number} && size(aa)==size(bb) && !isempty(aa)
            finite=isfinite.(aa)
            @assert finite==isfinite.(bb) "nonfinite history mask changed at $key"
            @assert isequal(aa[.!finite],bb[.!finite]) "nonfinite history values changed at $key"
            allowed_initial = startswith(key,"diffusion stress ") || key=="equivalent plastic strain max [-]"
            if any(.!finite)
                @assert allowed_initial && ndims(aa)==2 && size(aa,2)==19 "unexpected nonfinite required data at $key"
                @assert all(isnan,aa[:,1]) && all(isnan,bb[:,1]) "initial placeholder convention differs at $key"
                @assert all(isfinite,aa[:,2:end]) && all(isfinite,bb[:,2:end]) "nonfinite committed mechanical history at $key"
                for index in findall(.!finite)
                    push!(nonfinite_locations,(key,index[1],index[2],"NaN initial uncommitted column; both sources"))
                end
            end
            push!(finite_rows,(key,join(size(aa),'x'),length(aa),count(finite),count(!,finite),any(.!finite) ? "only initial column 1 NaN; columns 2:19 finite" : "all finite"))
            if any(finite)
                av,bv=aa[finite],bb[finite]
                push!(rows,(key,join(size(aa),'x'),maximum(abs,av-bv),norm(av-bv),maximum(abs,bv),maximum(abs,av),count(!,finite)))
            end
        elseif aa isa Number && bb isa Number
            @assert isfinite(aa) && isfinite(bb) "nonfinite scalar result at $key"
            push!(finite_rows,(key,"scalar",1,1,0,"all finite"))
            push!(rows,(key,"scalar",abs(aa-bb),abs(aa-bb),abs(bb),abs(aa),0))
        end
    end
    for (filename,header,values) in (
        ("full-result-finiteness.csv","key,shape,total_count,finite_count,nonfinite_count,location_contract",finite_rows),
        ("nonfinite-locations.csv","key,row,column,source",nonfinite_locations))
        open(joinpath(OUT66,filename),"w") do io
            println(io,header)
            for row in values; println(io,join(('"'*replace(string(v),'"'=>"\"\"")*'"' for v in row),',')); end
        end
    end
    open(joinpath(OUT66,"full-result-differences.csv"),"w") do io
        println(io,"key,shape,max_abs_difference,L2_difference,before_max_abs,after_max_abs,matching_nonfinite_count")
        for row in rows
            println(io,join(('"'*replace(string(v),'"'=>"\"\"")*'"' for v in row),','))
        end
    end
    @assert a["time [s]"] == b["time [s]"]
    @assert length(prior.chemistry)==length(current.chemistry)
    ch=Any[]
    for (i,(bc,ac)) in enumerate(zip(prior.chemistry,current.chemistry))
        @assert bc.t_eval==ac.t_eval
        @assert all(isfinite,bc.y_before) && all(isfinite,bc.y_after) &&
            all(isfinite,ac.y_before) && all(isfinite,ac.y_after) "nonfinite actual accepted chemical state at event $i"
        push!(ch,(i,bc.t_eval,maximum(abs,ac.y_before-bc.y_before),maximum(abs,ac.y_after-bc.y_after)))
    end
    open(joinpath(OUT66,"actual-state-differences.csv"),"w") do io
        println(io,"event,t_eval,max_y_before_difference,max_y_after_difference")
        for row in ch; println(io,join(row,',')); end
    end
    # Independently confirm the incidental frozen-before capture matches the reused oracle.
    extra_packet=deserialize(joinpath(ROOT66,"output","test_thermal_fem_rotation","baseline-before-capture","baseline-result.bin"))
    extra=extra_packet isa NamedTuple ? extra_packet.result : extra_packet
    reused=Any[]
    for key in sort(collect(intersect(keys(extra),keys(b))))
        startswith(key,"timing ") && continue
        if extra[key] isa AbstractArray{<:Number} && b[key] isa AbstractArray{<:Number} && size(extra[key])==size(b[key])
            @assert isequal(extra[key],b[key]) "incidental before capture differs at $key"
            push!(reused,(key,true))
        end
    end
    @assert !isempty(reused)
    @assert length(reused)==count(key -> !startswith(key,"timing ") && extra[key] isa AbstractArray{<:Number} &&
        b[key] isa AbstractArray{<:Number} && size(extra[key])==size(b[key]),intersect(keys(extra),keys(b)))
    open(joinpath(OUT66,"reuse-oracle-numeric-equality.csv"),"w") do io
        println(io,"key,exact_equal")
        for row in reused; println(io,'"'*row[1]*"\",true"); end
    end
    println("Compared $(length(rows)) numeric result fields and $(length(ch)) accepted state events; reused Task65 before oracle agrees exactly with incidental before capture on $(length(reused)) fields.")
    exit()
end
const ROT66 = (0.37, 0.71, pi/4)
rot66(a) = [cos(a) -sin(a); sin(a) cos(a)]
rel66(a, b) = norm(a-b)/norm(b)
csv66(name, header, rows) = open(joinpath(OUT66, name), "w") do io
    println(io, join(header, ','))
    for row in rows
        println(io, join(row, ','))
    end
end

function params66(kind=:ratio)
    p = JuBat.ChooseCell("Jellyroll")
    if kind === :ratio
        # Equal fractions: arithmetic mean=30 and harmonic mean=1 W/(m K).
        p.NE.lambda = 30-sqrt(870.0)
        p.SP.lambda = 30+sqrt(870.0)
        p.cell.lambda_r, p.cell.lambda_t = 1.0, 30.0
    end
    # This dimensional fixture is configured BEFORE SetCase; no frozen params edits.
    p.cell.h = 20.0
    p.scale.h = p.cell.h*p.scale.L/p.cell.lambda_r
    return p
end

function weights66(p, ne, kind)
    f = kind === :ratio ? [0.5,0.5,0.,0.,0.] :
        [2p.NE.thickness,2p.SP.thickness,2p.PE.thickness,p.PCC.thickness,p.NCC.thickness]
    f ./= sum(f)
    return repeat(permutedims(f), ne, 1)
end

function case66(X, E; kind=:ratio)
    p = params66(kind)
    c = JuBat.SetCase(p, JuBat.Option(model="SPMe", thermal_enabled=true,
        thermalmodel="distributed2D", cool_method="none"))
    gs = JuBat.GetGS(E, X, 2, "Q4")
    m = JuBat.Mesh("Q4",2,X,size(X,1),E,gs)
    c.mesh["thermal2D"] = m
    ne = size(E,1)
    c.mesh["thermal2D layer weights"] = weights66(p,ne,kind)
    c.mesh["thermal2D boundary edges"] = JuBat.BoundaryEdgeCache(Tuple{Int,Int}[],Float64[])
    return c
end

function assemble66(c, q; branch=:cartesian)
    v = Dict{String,Any}("heat_source_fields"=>q)
    return branch === :cartesian ? JuBat.ThermalDistributed2D(c,v) :
        JuBat.ThermalDistributed2D_Ring(c,v)
end

# Independent test-only algebra permits four mechanism combinations; production has no switch.
function ref66(c; gradient=true, tensor=true)
    m = c.mesh["thermal2D"]; f = c.mesh["thermal2D layer weights"][1,:]
    ls = [c.param.NE.lambda,c.param.SP.lambda,c.param.PE.lambda,
          c.param.PCC.lambda,c.param.NCC.lambda]
    lr, lt = inv(sum(f./ls)), sum(f.*ls)
    K = zeros(m.nlen,m.nlen)
    for g in eachindex(m.gs.ele)
        ids = m.element[m.gs.ele[g],:]
        A = transpose(last(JuBat.LagrangeBasis("Q4",2,vec(m.gs.xi[g,:]))))
        J = transpose(A)*m.node[ids,:]
        G = gradient ? A/transpose(J) : A/J
        x,y = m.gs.x[g,:]; er = [x,y]/hypot(x,y); et = [-er[2],er[1]]
        k = lr*(er*er')+lt*(et*et')
        if !tensor; k[1,2] *= -1; k[2,1] *= -1; end
        K[ids,ids] -= G*k*G'*(m.gs.weight[g]*m.gs.detJ[g])
    end
    return K
end

function robin66(K,F,m,edges,hnd)
    # Independent edge weak form with GIVEN SI h converted by L/lambda_scale.
    # Only the operator is under test; Task68 owns production boundary scales.
    H = -copy(K); b = copy(F)
    for (i,j) in edges
        len = norm(m.node[i,:]-m.node[j,:])
        H[[i,j],[i,j]] += hnd*len/6 .* [2. 1.;1. 2.]
    end
    return H,b
end

function annulus66(nr,nt,a=0.)
    # Independent polygon annulus; radius groups are integer generation indices.
    p = params66(); r = collect(range(p.cell.Rin/p.scale.L,p.cell.Rout/p.scale.L;length=nr+1))
    X = zeros((nr+1)*nt,2)
    for ir in 1:nr+1, it in 1:nt
        theta = 2pi*(it-1)/nt+a
        X[(ir-1)*nt+it,:] = r[ir].*[cos(theta),sin(theta)]
    end
    idx(ir,it) = (ir-1)*nt+mod1(it,nt)
    E = reduce(vcat,[permutedims([idx(ir,it),idx(ir+1,it),idx(ir+1,it+1),idx(ir,it+1)]) for ir in 1:nr for it in 1:nt])
    edges = [(idx(nr+1,it),idx(nr+1,it+1)) for it in 1:nt]
    return case66(X,E), r, edges
end

function consumer66()
    c = JuBat.SetCase(JuBat.ChooseCell("Jellyroll"), JuBat.Option(model="SPMe",
        thermal_enabled=true,thermalmodel="distributed2D",per_element_spme=true))
    md = JuBat.jellyroll_collector_seed_mesh(c.param;nθ=12,czm_enabled=true,gsorder=2)
    c = JuBat.setup_thermal2D_mesh(c,md)
    c.mesh["czm"] = JuBat.create_czm_mesh(md.czm_submesh,c.mesh["thermal2D"],c.param)
    c.mech = JuBat.MechState(c.mesh["czm"])
    m = c.mesh["thermal2D"]; ne = size(m.element,1)
    v = Dict{String,Union{Array{Float64},Float64}}("T_nodes"=>fill(c.param.cell.T0+0.001,m.nlen),
        "thermal2D element soc_n"=>fill(c.param.NE.cs0-0.02,ne),
        "thermal2D element soc_p"=>fill(c.param.PE.cs0+0.02,ne))
    K,F = JuBat.assemble_bonded_elastic(c,v)
    bulk = JuBat.assemble_bulk_stiffness(c.mesh["czm"],c.param)
    oneD = Dict(k=>(copy(mm.gs.Ni),copy(mm.gs.dNidx),copy(mm.gs.detJ)) for (k,mm) in c.mesh if mm isa JuBat.Mesh && mm.dimension==1)
    for typ in ("L2","L3")
        mm=JuBat.SetMesh([0.,1.],4,typ,2)
        oneD[typ]=(mm.gs.Ni,mm.gs.dNidx,mm.gs.detJ)
        for g in eachindex(mm.gs.ele)
            @test abs(sum(mm.gs.dNidx[g,:]))<=1e-12
            @test abs(dot(mm.node[mm.element[mm.gs.ele[g],:]],mm.gs.dNidx[g,:])-1)<=1e-12
        end
    end
    Xpatch=[2. .3;3.1 .55;3.3 .92;2.2 .64]
    JuBat.IntQ4(Xpatch[:,1],Xpatch[:,2]) do xi,eta,w,gx,gy,detj
        G=hcat(gx,gy)
        @test norm(Xpatch'*G-I)<=1e-12
        # Small strain admits translation and infinitesimal rotation.
        @test norm(ones(4)'*G)<=1e-12
        gradu=([0. -1.;1. 0.]*Xpatch')*G
        @test norm(gradu+gradu')<=1e-12
    end
    # COH2D4 derivatives are a separate arc-length contract and not changed by Q4.
    X = [0. 0.; 2. 1.; 2.1 1.2; 0.1 0.2]; E = reshape(Int64[1,2,3,4],1,4)
    co = JuBat.GetGS(E,X,2,"COH2D4")
    MT,_,FT = assemble66(c,fill(1.0,ne))
    data = (K=K,F=F,bulk=bulk,oneD=oneD,co=(co.Ni,co.dNidx,co.detJ),
        thermal_unchanged=(MT,FT,m.gs.Ni,m.gs.x,m.gs.detJ,m.gs.weight))
    serialize(joinpath(OUT66,"consumer-state.bin"),data)
    before = joinpath(ROOT66,"output","test_thermal_fem_rotation","red","consumer-state.bin")
    if RUN66 != "red" && isfile(before)
        old = deserialize(before)
        @test data == old
    end
    return nothing
end

const PATCH66 = Any[]; const COMB66 = Any[]; const ANN66 = Any[]
const TENSOR66=Any[]; const WHOLE66=Any[]
try
@testset "Task66 real production thermal FEM" begin
    @testset "fixed-state unaffected consumers" begin
        consumer66()
    end
    @testset "Q4 patch, covariance and material axes" begin
        X0 = [2. 0.3;3.1 0.55;3.3 0.92;2.2 0.64]
        E0 = reshape(Int64[1,2,3,4],1,4)
        for (shape,X) in ((:moderate,X0),(:thin,[2. 0.3;3.1 0.55;3.10002 0.55007;2.00002 0.30007]))
            c = case66(X,E0); m = c.mesh["thermal2D"]
            for a in ROT66
                R = rot66(a); cr = case66(X*R',E0); mr = cr.mesh["thermal2D"]
                for g in eachindex(m.gs.ele)
                    G = hcat(m.gs.dNidx[g,1:4],m.gs.dNidx[g,5:8])
                    Gr = hcat(mr.gs.dNidx[g,1:4],mr.gs.dNidx[g,5:8])
                    A = last(JuBat.LagrangeBasis("Q4",2,vec(m.gs.xi[g,:])))
                    condJ = cond(A*X)
                    patch = norm(X'*G-I); cov = rel66(Gr,G*R')
                    push!(PATCH66,(shape,a,g,condJ,patch,cov))
                    # Thin condition-dependent budget is separate from fixed 1e-12 moderate gate.
                    budget = shape===:moderate ? 1e-12 : 20eps()*condJ
                    @test patch <= budget
                    @test cov <= budget
                    @test norm(vec(sum(G;dims=1))) <= budget
                    @test all(>(0.), m.gs.detJ)
                end
            end
        end
        for kind in (:ratio,:default)
            c = case66(X0,E0;kind=kind); p = c.param
            f = c.mesh["thermal2D layer weights"][1,:]
            ls = [p.NE.lambda,p.SP.lambda,p.PE.lambda,p.PCC.lambda,p.NCC.lambda]
            lr,lt = inv(sum(f./ls)),sum(f.*ls)
            for theta in (0.,pi/8,pi/4,.63)
                kx,kxy,ky = JuBat.thermal_anisotropic_conductivity_2d(p,c.mesh["thermal2D layer weights"],[1],[cos(theta)],[sin(theta)])
                k = [kx[1] kxy[1];kxy[1] ky[1]]
                er,et = [cos(theta),sin(theta)],[-sin(theta),cos(theta)]
                push!(TENSOR66,(kind,theta,lr,lt,dot(er,k*er),dot(et,k*et),dot(er,k*et)))
                @test abs(dot(er,k*er)-lr)/lt <= 1e-12
                @test abs(dot(et,k*et)-lt)/lt <= 1e-12
                @test abs(dot(er,k*et))/lt <= 1e-12
                for a in ROT66
                    x,y,z = JuBat.thermal_anisotropic_conductivity_2d(p,c.mesh["thermal2D layer weights"],[1],[cos(theta+a)],[sin(theta+a)])
                    @test rel66([x[1] y[1];y[1] z[1]],rot66(a)*k*rot66(a)') <= 1e-12
                end
            end
            M,K,F = assemble66(c,[1.])
            for a in ROT66
                cr = case66(X0*rot66(a)',E0;kind=kind)
                Mr,Kr,Fr = assemble66(cr,[1.])
                @test rel66(Kr,K) <= 1e-12
                @test rel66(Mr,M) <= 1e-12
                @test rel66(Fr,F) <= 1e-12
                for grad in (false,true), tensor in (false,true)
                    kb,kbr = ref66(c;gradient=grad,tensor=tensor),ref66(cr;gradient=grad,tensor=tensor)
                    push!(COMB66,(kind,a,grad,tensor,rel66(kbr,kb),rel66(kb,K)))
                end
            end
            H = -K
            @test rel66(H,H') <= 1e-12
            @test norm(H*ones(4))/norm(H) <= 1e-12
            @test dot(X0[:,1],H*X0[:,1])>0
        end
    end
    @testset "two distorted elements, rotated boundary and solved temperature rise" begin
        X = [2. .3;3.1 .55;4.2 .4;2.2 .64;3.3 .92;4.4 .95]
        E = Int64[1 2 5 4;2 3 6 5]
        for kind in (:ratio,:default),branch in (:cartesian,:ring)
            c = case66(X,E;kind=kind); M,K,F = assemble66(c,[1.,2.];branch=branch)
            H,b = robin66(K,F,c.mesh["thermal2D"],[(3,6)],0.1)
            u = H\b
            @test norm(H*u-b)/norm(b) <= 1e-10
            writedlm(joinpath(OUT66,"whole-$(kind)-$(branch)-M.csv"),Matrix(M),',')
            writedlm(joinpath(OUT66,"whole-$(kind)-$(branch)-K.csv"),Matrix(K),',')
            writedlm(joinpath(OUT66,"whole-$(kind)-$(branch)-F-u.csv"),hcat(F,u),',')
            for a in ROT66
                cr = case66(X*rot66(a)',E;kind=kind)
                Mr,Kr,Fr = assemble66(cr,[1.,2.];branch=branch)
                Hr,br = robin66(Kr,Fr,cr.mesh["thermal2D"],[(3,6)],0.1)
                ur = Hr\br
                push!(WHOLE66,(kind,branch,a,rel66(Mr,M),rel66(Kr,K),rel66(Fr,F),rel66(Hr,H),rel66(br,b),rel66(ur,u),norm(Hr*ur-br)/norm(br)))
                @test rel66(Mr,M)<=1e-12
                @test rel66(Kr,K)<=1e-12
                @test rel66(Fr,F)<=1e-12
                @test rel66(Hr,H)<=1e-12
                @test rel66(br,b)<=1e-12
                @test rel66(ur,u)<=1e-10
                @test norm(Hr*ur-br)/norm(br)<=1e-10
            end
        end
    end
    @testset "independent annulus SI radial analytic convergence" begin
        qSI=10000.; hSI=20.
        for branch in (:cartesian,:ring), bc in (:robin,:dirichlet)
            previous = Inf
            previousmax = Inf
            for (nr,nt) in ((4,24),(8,48),(16,96))
                c,r,edges = annulus66(nr,nt); p = c.param; scale = p.scale
                qnd = qSI/scale.q; m = c.mesh["thermal2D"]
                M,K,F = assemble66(c,fill(qnd,size(m.element,1));branch=branch)
                H,b = robin66(K,F,m,edges,hSI*scale.L/scale.lambda)
                outer = nr*nt .+ (1:nt)
                if bc===:dirichlet
                    H,b = -K, F
                    free = setdiff(1:m.nlen,outer)
                    u = zeros(m.nlen); u[free] = H[free,free]\b[free]
                    residual = norm((H*u-b)[free])/norm(b[free])
                else
                    u = H\b; residual = norm(H*u-b)/norm(b)
                end
                Ri,Ro = c.param_dim.cell.Rin,c.param_dim.cell.Rout
                rr = r*scale.L; lrSI = c.param_dim.cell.lambda_r
                exact = qSI/(4lrSI).*(Ro^2 .- rr.^2) .+ qSI*Ri^2/(2lrSI).*log.(rr/Ro)
                if bc===:robin; exact .+= qSI/(2hSI)*(Ro-Ri^2/Ro); end
                exactnodes = repeat(exact;inner=nt)
                usi = u*scale.T_ref
                # Area lumping only defines the error norm, not the reference solution.
                areaweights = vec(sum(M;dims=2)); areaweights ./=sum(areaweights)
                l2 = sqrt(sum(areaweights.*(usi-exactnodes).^2))/sqrt(sum(areaweights.*exactnodes.^2))
                maxerr = maximum(abs,usi-exactnodes)/maximum(exactnodes)
                dispersion = maximum((maximum(u[(ir-1)*nt .+ (1:nt)])-minimum(u[(ir-1)*nt .+ (1:nt)]))/abs(sum(u[(ir-1)*nt .+ (1:nt)])/nt) for ir in 1:nr)
                # SI source power per metre and external loss from the same discrete weak form.
                source = sum(F)*scale.q*scale.L^2
                loss = bc===:robin ? sum((H+K)*u)*scale.q*scale.L^2 : sum((K*u+F)[outer])*scale.q*scale.L^2
                push!(ANN66,(branch,bc,nr,nt,l2,maxerr,dispersion,residual,source,loss))
                @test l2 < previous
                @test maxerr < previousmax
                @test dispersion <= 1e-11
                @test residual <= 1e-10
                @test abs(loss-source)/source <= 1e-10
                previous,previousmax = l2,maxerr
                writedlm(joinpath(OUT66,"annulus-$(branch)-$(bc)-$(nr)-field.csv"),hcat(m.node,usi,exactnodes),',')
                if nr==4
                    for a in ROT66
                        cr,_,er = annulus66(nr,nt,a)
                        Mr,Kr,Fr = assemble66(cr,fill(qnd,size(m.element,1));branch=branch)
                        @test rel66(Kr,K)<=1e-12
                        @test rel66(Mr,M)<=1e-12
                        @test rel66(Fr,F)<=1e-12
                        Hr,br = robin66(Kr,Fr,cr.mesh["thermal2D"],er,hSI*scale.L/scale.lambda)
                        if bc===:dirichlet
                            ur=zeros(m.nlen); free=setdiff(1:m.nlen,outer)
                            ur[free]=(-Kr)[free,free]\Fr[free]
                        else
                            ur=Hr\br
                        end
                        @test rel66(ur,u)<=1e-10
                    end
                end
            end
        end
        # One physical FVM reference uses the actual entry, constant q only.
        c,r,_ = annulus66(8,48)
        qnd=qSI/c.param.scale.q; md=(r=r,ntheta=48)
        Mf,Kf,Ff = JuBat.ThermalPolar2D_Ring(c,Dict{String,Any}("heat_source_fields"=>fill(qnd,8*48)),md)
        uf=(-Kf)\Ff .- c.param.cell.T_amb
        Ri,Ro=c.param_dim.cell.Rin,c.param_dim.cell.Rout
        rr=r*c.param.scale.L
        exact=qSI/4 .* (Ro^2 .- rr.^2) .+ qSI*Ri^2/2 .* log.(rr/Ro) .+ qSI/(2hSI)*(Ro-Ri^2/Ro)
        fvmerr=norm(uf*c.param.scale.T_ref-repeat(exact;inner=48))/norm(repeat(exact;inner=48))
        @test fvmerr < .01
        @test norm(Kf*(uf .+ c.param.cell.T_amb)+Ff)/norm(Ff)<=1e-10
        csv66("fvm-reference.csv",["nr","nt","relative_radial_error"],[(8,48,fvmerr)])
    end
end
finally
    csv66("rotation-and-patch.csv",["shape","alpha","gauss","condJ","patch_absolute","gradient_covariance_relative"],PATCH66)
    csv66("four-combinations.csv",["material","alpha","gradient_correct","tensor_correct","rotation_error","production_difference"],COMB66)
    csv66("annulus-convergence.csv",["branch","bc","nr","nt","relative_area_L2","relative_max","same_radius_relative_dispersion","solve_residual","source_W_per_m","loss_W_per_m"],ANN66)
    csv66("material-projection.csv",["material","theta","lambda_r","lambda_t","radial_projection","tangent_projection","mixed_projection"],TENSOR66)
    csv66("wholemesh-rotation.csv",["material","branch","alpha","M_relative","K_relative","F_relative","H_BC_relative","b_BC_relative","rise_solution_relative","solve_residual"],WHOLE66)
    open(joinpath(OUT66,"fixture.json"),"w") do io
        print(io,"{\"nodes\":[[2,0.3],[3.1,0.55],[3.3,0.92],[2.2,0.64]],\"elements\":[[1,2,3,4]],\"angles\":[0.37,0.71,0.7853981633974483],\"SI_lambda_r_t\":[1,30],\"SI_q\":10000,\"SI_h\":20,\"annulus_grids\":[[4,24],[8,48],[16,96]],\"parameter_contract\":\"configured before SetCase, frozen afterwards\",\"Robin_boundary\":\"independent given SI h weak form, Task68 boundary not validated\"}")
    end
end
