using Test
using LinearAlgebra
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

# Fixed-old-history differentiation: every perturbation starts from the same state.
statekey(s) = (s.D, s.D_visc, s.δ_max_n, s.δ_max_t, s.δ_max_eff,
               s.fractured, s.accumulated_damage)
function fixture(; equal=false, swapped=false)
    ip = deepcopy(JuBat.SetCase(JuBat.ChooseCell("Jellyroll"), JuBat.Option()).param.PCC)
    ip.K_n = ip.K_t = 1.0
    ip.δ_0, ip.δ_c = swapped ? (0.1, 0.5) : (0.4, 1.0)
    ip.δ_0_t, ip.δ_c_t = equal ? (ip.δ_0, ip.δ_c) : swapped ? (0.4, 1.0) : (0.1, 0.5)
    ip.σ_max = ip.K_n * ip.δ_0
    ip.τ_max = ip.K_t * ip.δ_0_t
    ip.G_c = ip.σ_max * ip.δ_c / 2
    ip.G_c_t = ip.τ_max * ip.δ_c_t / 2
    ip.eta = 1.45
    ip
end
function history(D, Dv=D; fractured=(D == Dv == 1.0))
    s = JuBat.DamageState()
    s.D, s.D_visc, s.fractured = D, Dv, fractured
    s.accumulated_damage = D
    s
end
trial(n,t,s,ip; model="mix", b=1.0) = JuBat.bilinear_traction_state(n,t,s,ip,model; visc_beta=b)
tangent(n,t,s,ip; model="mix", b=1.0) = JuBat.bilinear_tangent(n,t,s,ip,model; visc_beta=b)
function fd(n,t,s,ip; model="mix", b=1.0, h=1e-6)
    J = zeros(2,2)
    for j in 1:2
        p = j==1 ? (n+h,t) : (n,t+h)
        m = j==1 ? (n-h,t) : (n,t-h)
        J[:,j] = (collect(trial(p...,s,ip; model,b)[1:2]) -
                  collect(trial(m...,s,ip; model,b)[1:2])) / (2h)
    end
    J
end

@testset "72 constitutive physics" begin
    ip = fixture()
    @testset "direction envelope replaces the geometric peak gate" begin
        fresh = JuBat.DamageState()
        first = last(trial(0.3,0.0,fresh,ip))
        second = last(trial(0.0,0.2,first,ip))
        @info "nonproportional counterexample" old_peak=first.δ_max_eff new_norm=0.2 actual=second.D target=0.625
        @test second.D ≈ 0.625
        @test second.δ_max_eff == 0.3
        for sign_t in (-1.0,1.0)
            @test last(trial(0.0,sign_t*0.2,first,ip)).D ≈ 0.625
            swapped = fixture(swapped=true)
            s = last(trial(0.0,sign_t*0.3,fresh,swapped))
            @test last(trial(0.2,0.0,s,swapped)).D ≈ 0.625
        end
        eq = fixture(equal=true)
        for (n,t) in ((0.6,0.0),(0.0,0.6),(0.36,0.48),(0.36,-0.48))
            @test last(trial(n,t,fresh,eq)).D ≈ 5/9
        end
        s = last(trial(0.7,0.0,fresh,ip))
        for (n,t) in ((0.4,0.0),(-0.2,0.0),(0.6,0.0),(0.0,0.2))
            @test last(trial(n,t,s,ip)).D >= s.D
        end
        @test last(trial(0.8,0.0,s,ip)).D > s.D
        closed=last(trial(-0.2,0.0,s,ip))
        reopened=last(trial(0.6,0.0,closed,ip))
        @test reopened.D == closed.D == s.D
        @test reopened.δ_max_n == s.δ_max_n
    end
    @testset "pure modes, model1 and closure" begin
        fresh = JuBat.DamageState()
        for (a,c,mode) in ((ip.δ_0,ip.δ_c,1),(ip.δ_0_t,ip.δ_c_t,2))
            for e in (0.0,a/2,a,(a+c)/2,c,1.2c)
                n,t = mode==1 ? (e,0.0) : (0.0,e)
                d = e<=a ? 0.0 : e>=c ? 1.0 : c*(e-a)/(e*(c-a))
                r = trial(n,t,fresh,ip)
                @test r[3] ≈ d
                @test collect(r[1:2]) ≈ (1-d)*[n,t]
            end
        end
        @test last(trial(0.0,2.0,fresh,ip; model="model1")).D == 0.0
        full = last(trial(1.1,2.0,fresh,ip; model="model1"))
        @test full.fractured
        @test trial(-0.2,0.3,full,ip; model="model1")[1:2] == (-0.2,0.3)
        @test tangent(-0.2,0.3,full,ip; model="model1") == Matrix{Float64}(I,2,2)
        @test trial(0.2,0.3,full,ip; model="model1")[1:2] == (0.0,0.3)
        @test tangent(0.2,0.3,full,ip; model="model1") == [0.0 0.0;0.0 1.0]
        @test trial(-0.2,0.3,full,ip)[1:2] == (-0.2,0.0)
        @test tangent(-0.2,0.3,full,ip) == [1.0 0.0;0.0 0.0]
    end
    @testset "same accepted history, viscous time and canonical saturation" begin
        s = history(0.7,0.2)
        before = statekey(s)
        for b in (0.0,0.25,1.0)
            r = trial(0.65,0.1,s,ip;b)
            @test r[4].D_visc ≈ s.D_visc+b*(r[3]-s.D_visc)
            @test statekey(last(trial(0.65,0.1,s,ip;b))) == statekey(r[4])
            @test statekey(s) == before
        end
        lag = history(1.0,0.5; fractured=false)
        for _ in 1:4
            next = last(trial(0.1,0.0,lag,ip;b=0.2))
            @test next.D == 1.0
            @test next.D_visc > lag.D_visc
            lag = next
        end
        # Use a canonical old state below the saturation band; control raw prediction via b.
        s = history(1.0,0.0;fractured=false)
        threshold = 1.0-1e-10
        for raw in (prevfloat(threshold),threshold,nextfloat(threshold))
            r = trial(0.2,0.1,s,ip;b=raw)
            saturated = raw>=threshold
            @test r[4].fractured == saturated
            @test r[4].D_visc == (saturated ? 1.0 : raw)
            @test r[1] == (1-r[4].D_visc)*0.2
            @test tangent(0.2,0.1,s,ip;b=raw) == (1-r[4].D_visc)*Matrix{Float64}(I,2,2)
            @test statekey(last(trial(0.2,0.1,r[4],ip;b=0.0))) == statekey(r[4])
        end
        # D<1 remains on the true softening law even inside the numeric saturation band.
        e=1.0-1e-11
        r=trial(e,0.0,JuBat.DamageState(),ip)
        @test r[3] < 1.0
        @test !r[4].fractured
        @test tangent(e,0.0,JuBat.DamageState(),ip)[1,1] ≈ -ip.δ_0/(ip.δ_c-ip.δ_0)
        @test tangent(e,0.0,JuBat.DamageState(),ip)[2,2] == 1-r[4].D_visc
        frozen=last(trial(e,0.0,r[4],ip;b=0.0))
        @test statekey(frozen)==statekey(r[4])
        full=history(1.0)
        @test tangent(0.4,0.2,full,ip) == zeros(2,2)
        # No silent migration or reconstruction of illegal preexisting histories.
        for bad in (history(0.8,0.8;fractured=true),history(1.0,0.5;fractured=true),
                    history(1.0,threshold;fractured=false),history(1.0,1.0;fractured=false),
                    history(0.5,0.6;fractured=false))
            @test_throws ArgumentError trial(0.2,0.1,bad,ip)
            @test_throws ArgumentError tangent(0.2,0.1,bad,ip)
        end
    end
    @testset "fixed history consistent derivative and one sided branches" begin
        for mat in (ip,fixture(equal=true)), b in (0.0,0.3,1.0),
            (n,t,s,model) in ((0.6,0.08,JuBat.DamageState(),"mix"),
                (0.1,0.25,JuBat.DamageState(),"mix"),
                (0.1,-0.25,JuBat.DamageState(),"mix"),
                (-0.2,0.25,JuBat.DamageState(),"mix"),
                (0.2,0.05,history(0.8,0.2),"mix"),
                (0.6,0.2,JuBat.DamageState(),"model1"))
            J=tangent(n,t,s,mat;model,b)
            @test J ≈ fd(n,t,s,mat;model,b) rtol=2e-5 atol=2e-7
        end
        first=last(trial(0.3,0.0,JuBat.DamageState(),ip))
        @test tangent(0.02,0.2,first,ip) ≈ fd(0.02,0.2,first,ip) rtol=2e-5
        @test !issymmetric(tangent(0.2,0.2,JuBat.DamageState(),ip))
        # Side samples approach each fold without crossing it.
        h=1e-6
        for (e,expected) in ((ip.δ_0-10h,1.0),(ip.δ_0+10h,-2/3),
                             (ip.δ_c-10h,-2/3),(ip.δ_c+10h,0.0))
            @test fd(e,0.0,JuBat.DamageState(),ip;h)[1,1] ≈ expected atol=1e-7
        end
        # Exact folds use their defined inactive derivative; evaluate the matching side only.
        @test tangent(ip.δ_0,0.0,JuBat.DamageState(),ip)[1,1] == 1.0
        @test tangent(ip.δ_c,0.0,JuBat.DamageState(),ip)[1,1] == 0.0
    end
    @testset "work and dissipation of the stated opening interpolation" begin
        for mat in (ip,fixture(equal=true)), angle in (0.0,0.3,pi/2,-pi/2)
            nn,nt=cos(angle),sin(angle)
            beta=abs(nt)
            a=sqrt(mat.δ_0^2+(mat.δ_0_t^2-mat.δ_0^2)*beta^mat.eta)
            c=sqrt(mat.δ_c^2+(mat.δ_c_t^2-mat.δ_c^2)*beta^mat.eta)
            k=mat.K_n*nn^2+mat.K_t*nt^2
            # Exact trapezoidal integration of the two linear traction segments.
            s=JuBat.DamageState(); W=0.0; old=[0.0,0.0]; Told=[0.0,0.0]
            for e in vcat(collect(range(0,a,length=251))[2:end], collect(range(a,c,length=751))[2:end])
                x=e*[nn,nt]; r=trial(x...,s,mat)
                T=collect(r[1:2]); W+=dot((T+Told)/2,x-old)
                s=r[4];old=x;Told=T
            end
            @test W ≈ 0.5*k*a*c rtol=1e-10
            if angle==0.0
                @test W ≈ mat.G_c
            elseif abs(angle)==pi/2
                @test W ≈ mat.G_c_t
            end
            if mat.δ_0==mat.δ_0_t
                @test W ≈ mat.G_c
            end
        end
        # For arbitrary nonproportional loading, degradation dissipation is
        # (psi0_n+psi0_t)*dD >= 0; no fixed nonproportional Gc is asserted.
        s=JuBat.DamageState(); dissipated=0.0
        for x in ((0.3,0.0),(0.0,0.2),(-0.1,-0.25),(0.2,-0.3),
                  (0.1,0.0),(0.8,0.1),(0.0,0.6))
            next=last(trial(x...,s,ip))
            increment=0.5*(ip.K_n*max(x[1],0)^2+ip.K_t*x[2]^2)*(next.D-s.D)
            @test increment>=0.0
            dissipated+=increment;s=next
        end
        @test dissipated>0.0
    end
    @testset "actual COH GP assembly from unchanged accepted history" begin
        opt=JuBat.Option();opt.thermal_enabled=true;opt.thermalmodel="distributed2D";opt.per_element_spme=true
        case=JuBat.SetCase(JuBat.ChooseCell("Jellyroll"),opt)
        @test case.param.PCC.δ_0==case.param.PCC.δ_0_t
        @test case.param.PCC.δ_c==case.param.PCC.δ_c_t
        @test case.param.NCC.δ_0==case.param.NCC.δ_0_t
        @test case.param.NCC.δ_c==case.param.NCC.δ_c_t
        md=JuBat.jellyroll_collector_seed_mesh(case.param;nθ=8,czm_enabled=true,gsorder=2)
        case=JuBat.setup_thermal2D_mesh(case,md)
        case.mesh["czm"]=JuBat.create_czm_mesh(md.czm_submesh,case.mesh["thermal2D"],case.param)
        cm=case.mesh["czm"];ms=JuBat.MechState(cm);ws=JuBat.CZMWorkspace(cm,case.param)
        bc_dofs,bc_vals=JuBat.extract_bc_dofs(cm,case.param)
        selected=findfirst(elem->all(node->!(2node-1 in bc_dofs)&&!(2node in bc_dofs),elem.nodes_top),cm.cohesive_elements)
        @test selected !== nothing
        elem=cm.cohesive_elements[selected]
        _,normal,_,_=JuBat.cohesive_local_frame(cm,elem)
        mip=JuBat.collector_params(case.param,elem.interface_type)
        opening=mip.δ_0*case.param.scale.δ_czm/case.param.scale.L
        u=zeros(2size(cm.node,1))
        for (node,factor) in zip(elem.nodes_top,(10.0,50.0))
            u[2node-1]=factor*opening*normal[1];u[2node]=factor*opening*normal[2]
        end
        old=JuBat.gp_history_for_step(ms,cm);keys=statekey.(old);force=zeros(length(u))
        run(v,s; b=1.0)=JuBat.arc_equilibrium_trial(cm,v,case.param,s,force,bc_dofs,bc_vals;czm_workspace=ws,czm_model="mix",visc_beta=b)
        peak=run(u,old);again=run(u,old)
        @test statekey.(old)==keys
        @test statekey.(peak.gp_damage_states)==statekey.(again.gp_damage_states)
        @test maximum(s.D for s in peak.gp_damage_states[selected,:])>0.9
        @test maximum(s.D for s in peak.gp_damage_states[selected,:])-minimum(s.D for s in peak.gp_damage_states[selected,:])>1e-3
        # Independent pure-normal interpolation at the actual quadrature points:
        # the prescribed top-node gaps are 10a and 50a, bottom-node gaps zero.
        weights,points=JuBat.NCweight(cm.czm_submesh.mesh_bonded.gs.order)
        expectedT=0.0
        for (g,xi) in enumerate(points)
            gap=((1-xi)/2*10+(1+xi)/2*50)*mip.δ_0
            d=gap>=mip.δ_c ? 1.0 : mip.δ_c*(gap-mip.δ_0)/(gap*(mip.δ_c-mip.δ_0))
            @test peak.gp_damage_states[selected,g].D ≈ d
            expectedT+=weights[g]*(1-d)*mip.K_n*gap/sum(weights)
        end
        @test peak.tractions[selected][1] ≈ expectedT rtol=1e-10
        @test abs(peak.tractions[selected][2])<1e-10*max(abs(expectedT),1.0)
        unload=run(0.9u,peak.gp_damage_states;b=0.0)
        @test [s.D for s in unload.gp_damage_states]==[s.D for s in peak.gp_damage_states]
        @test [s.D_visc for s in unload.gp_damage_states]==[s.D_visc for s in peak.gp_damage_states]
        h=1e-6;plus=run((1+h)*u,old);minus=run((1-h)*u,old)
        free=setdiff(1:length(u),bc_dofs)
        numerical=-((plus.R-minus.R)/(2h))[free]
        @test norm((peak.K_total*u)[free]-numerical)/norm(numerical)<1e-4
        @test statekey.(old)==keys
    end
end
