using Serialization, LinearAlgebra
include(joinpath(@__DIR__,"..","src","JuBat.jl"))
using .JuBat
const C_OUT = normpath(joinpath(@__DIR__,"..","output","validate_case_batch_c_trace"))
mkpath(C_OUT)

function mech_numbers(ms)
    ms === nothing && return nothing
    return (;u=copy(ms.u_prev), damage=map(s->Tuple(getfield(s,f) for f in fieldnames(typeof(s))),ms.damage_states),
        gp=ms.gp_damage_states === nothing ? nothing : map(s->Tuple(getfield(s,f) for f in fieldnames(typeof(s))),ms.gp_damage_states),
        plastic=ms.plastic_states === nothing ? nothing : map(s->(;eps_p=s.eps_p,kappa=s.kappa,sigma=s.sigma),ms.plastic_states),
        prestress=deepcopy(ms.winding_prestress),node_ref=deepcopy(ms.node_ref),initialized=ms.mechanics_initialized)
end
function numerical_result(result)
    Dict(k=>deepcopy(v) for (k,v) in result if !startswith(k,"timing ") && k!="final_state")
end
function fixture_c(kind)
    p=JuBat.ChooseCell("Jellyroll")
    opt=JuBat.Option(model="SPMe",Nn=3,Ns=2,Np=3,Nrn=4,Nrp=4,
        thermalmodel="distributed2D",thermal_enabled=true,per_element_spme=true)
    opt.mechanicalmodel="none";opt.Current=t->0.5;opt.time=[0.0,0.03];opt.dt=[0.01,0.01];opt.dtType="constant"
    opt.czm.enabled=kind in ("gp","geo")
    opt.czm.area_loss_enabled=true;opt.czm.update_interval=2
    opt.czm.iter_method=kind=="gp" ? "gp_basic" : "basic"
    opt.czm.geo_nonlinear=kind=="geo";opt.czm.j2_plasticity=kind=="geo";opt.czm.winding_prestress=kind=="geo"
    c=JuBat.SetCase(p,opt)
    md=JuBat.jellyroll_collector_seed_mesh(c.param;nθ=8,gsorder=2,czm_enabled=opt.czm.enabled)
    c=JuBat.setup_thermal2D_mesh(c,md)
    if opt.czm.enabled
        c.mesh["czm"]=JuBat.create_czm_mesh(md.czm_submesh,c.mesh["thermal2D"],c.param)
    end
    return c
end

function phase_trace(kind)
    c=fixture_c(kind);state=(y=copy((c.opt.per_element_spme ? JuBat.ModelInitialisation_MultiSPMe(c) : JuBat.ModelInitialisation(c))),mech=deepcopy((haskey(c.mesh,"czm") ? JuBat.MechState(c.mesh["czm"]) : nothing)),t_global=0.0)
    events=Any[];mechanics=Any[];branches=Any[];phases=Any[]
    JuBat.mechanical_step_observer[]=ev->push!(mechanics,(;phase=ev.phase,dt=ev.dt_seconds,
        mech=mech_numbers(ev.ms),result=ev.result===nothing ? nothing :
        (;u=copy(ev.result.displacement),damage=copy(ev.result.damage),converged=ev.result.converged,residual=ev.result.residual_norm)))
    JuBat.branch_consumption_observer[]=ev->push!(branches,(;t=ev.t,D=deepcopy(ev.D_elem),w=copy(ev.w),I=copy(ev.I_e)))
    try
        for (phase,I) in ((JuBat.PHASE_DISCHARGE,0.5),(JuBat.PHASE_REST,0.0),(JuBat.PHASE_CHARGE,-0.5),(JuBat.PHASE_REST,0.0))
            result=JuBat.solve_phase(c,phase,0.03,state;dt_range=[0.01,0.01],
                chemical_step_observer=ev->push!(events,deepcopy(ev)))
            state=result.final_state
            push!(phases,(;result=numerical_result(result.solve_result),y=copy(state.y),
                t_global=state.t_global,mech=mech_numbers(state.mech)))
        end
    finally
        JuBat.mechanical_step_observer[]=nothing;JuBat.branch_consumption_observer[]=nothing
    end
    return (;events,mechanics,branches,phases)
end

function export_trace()
    c=fixture_c("export");state=nothing;phases=Any[]
    for (phase,I) in ((JuBat.PHASE_DISCHARGE,0.5),(JuBat.PHASE_REST,0.0))
        result,data=JuBat.solve_phase_with_export(c,phase,0.03,state;dt_range=[0.01,0.01])
        state=result.final_state
        push!(phases,(;y=copy(state.y),T=JuBat.get_thermal_dofs(state.y,c),V=result.V_end,t=state.t_global,
            samples=[Tuple(getfield(s,f) for f in fieldnames(typeof(s))) for s in data]))
    end
    return phases
end

function reject_trace()
    p=JuBat.ChooseCell("Jellyroll");o=JuBat.Option(model="SPMe",Nn=3,Ns=2,Np=3,Nrn=4,Nrp=4,
        thermalmodel="none",per_element_spme=false)
    o.Current=t->t<0.006 ? 0.0 : 0.5;o.time=[0.0,0.025];o.dt=[0.001,0.008];o.dtType="auto";o.dtThreshold=1e-10
    c=JuBat.SetCase(p,o);events=Any[]
    r=JuBat.Solve(c;return_final_state=true,chemical_step_observer=ev->push!(events,deepcopy(ev)))
    return (;result=numerical_result(r),y=r["final_state"].y,events)
end

if abspath(PROGRAM_FILE)==@__FILE__
    println("Historical B-to-C comparison retired: current temperature/current contract intentionally changes its exported trace. Reusable fixtures remain available via include; historical output is preserved.")
end
