#= Case 2 of czm_cycling_n20: one discharge/rest half-cycle, no charge phase.

Usage: julia -t1 --startup-file=no --project=. example/check_czm_n20_halfcycle.jl RUN_ID
Original setup/default parameters are shared; this entry only chooses two phases
and saves their real states/results. It neither calibrates nor repairs the solver.
=#
using Dates, SHA, Serialization, TOML, Statistics
include(joinpath(@__DIR__, "czm_cycling_n20.jl"))

savebin(path, x) = open(io -> serialize(io,x), path, "w")
function source_identity(outdir)
    root=normpath(joinpath(@__DIR__,".."))
    paths=String[]
    for (dir,_,files) in walkdir(joinpath(root,"src")), file in files
        push!(paths,relpath(joinpath(dir,file),root))
    end
    append!(paths,["Project.toml","Manifest.toml","example/czm_cycling_n20.jl",
        "example/cycling_csv_time.jl","example/check_czm_n20_halfcycle.jl"])
    open(joinpath(outdir,"source_manifest.csv"),"w") do io
        println(io,"relative_path,sha256")
        for p in sort!(paths)
            println(io,replace(p,'\\'=>'/'),",",bytes2hex(sha256(read(joinpath(root,p)))))
        end
    end
    write(joinpath(outdir,"source.diff"),read(`git -C $root diff -- src example`,String))
    write(joinpath(outdir,"git_status.txt"),read(`git -C $root status --short`,String))
    write(joinpath(outdir,"parent_cycling_script.jl"),read(`git -C $root show HEAD:example/czm_cycling_n20.jl`))
    return readchomp(`git -C $root rev-parse HEAD`)
end

function export_state(outdir,label,state,case)
    savebin(joinpath(outdir,label*"_accepted.bin"),state)
    ms=state.mech
    open(joinpath(outdir,label*"_state.txt"),"w") do io
        println(io,"t_global_s=",state.t_global)
        println(io,"mechanics_initialized=",ms.mechanics_initialized)
        println(io,"gp_history_present=",ms.gp_damage_states !== nothing)
        println(io,"temperature_K_extrema=",extrema(JuBat.get_thermal_dofs(state.y,case).*case.param.scale.T_ref))
        println(io,"element_D_extrema=",extrema(s.D for s in ms.damage_states))
        println(io,"element_Dv_extrema=",extrema(s.D_visc for s in ms.damage_states))
        println(io,"element_fractured=",count(s->s.fractured,ms.damage_states))
    end
    # A fresh MechState has no GP history; record that fact, never invent GP zeros.
    if ms.gp_damage_states !== nothing
        open(joinpath(outdir,label*"_gp.csv"),"w") do io
            println(io,"cohesive,gp,interface,D,D_visc,fractured,delta_max_n,delta_max_t,delta_max_eff,accumulated_damage")
            for k in axes(ms.gp_damage_states,1), g in axes(ms.gp_damage_states,2)
                s=ms.gp_damage_states[k,g]
                println(io,join((k,g,case.mesh["czm"].cohesive_elements[k].interface_type,
                    s.D,s.D_visc,s.fractured,s.δ_max_n,s.δ_max_t,s.δ_max_eff,s.accumulated_damage),","))
            end
        end
    end
end

function phase_metrics(case,label,r)
    sr=r.solve_result
    xx=sr["diffusion stress xx [Pa]"][:,end]
    yy=sr["diffusion stress yy [Pa]"][:,end]
    xy=sr["diffusion stress xy [Pa]"][:,end]
    cm=case.mesh["czm"];hoop=Float64[];shear=Float64[]
    for e in axes(cm.bulk_element,1)
        nodes=cm.bulk_element[e,:]
        x=mean(cm.node[nodes,1]);y=mean(cm.node[nodes,2]);r2=x*x+y*y
        push!(hoop,(y*y*xx[e]-2x*y*xy[e]+x*x*yy[e])/r2)
        push!(shear,(x*y*(yy[e]-xx[e])+(x*x-y*y)*xy[e])/r2)
    end
    gp=r.final_state.mech.gp_damage_states
    gp === nothing && error("Completed gp_basic phase has no real GP history")
    # geo/J2 are off: this is recovered element stress history, not a J2 GP maximum.
    vm=vec(sr["diffusion stress vonMises [Pa]"])
    first_valid=findfirst(isfinite,vm)
    first_valid===nothing && error("No committed stress history")
    all(isnan,@view vm[1:first_valid-1]) || error("Invalid leading stress history")
    all(isfinite,@view vm[first_valid:end]) || error("Nonfinite committed stress history")
    return (label,r.t_start,r.t_end,r.duration,r.V_start,r.V_end,r.capacity,
        sr["termination_reason"],r.T_max,r.T_mean_end,r.D_max,r.D_mean,
        maximum(s.D for s in gp),maximum(s.D_visc for s in gp),count(s->s.fractured,gp),
        minimum(hoop),maximum(hoop),minimum(shear),maximum(shear),maximum(@view vm[first_valid:end]),length(sr["time [s]"]))
end

function main()
    length(ARGS)==1 || error("Supply one unique RUN_ID")
    outdir=normpath(joinpath(@__DIR__,"..","output","check_czm_n20_halfcycle",ARGS[1]))
    mkpath(outdir)
    # Enumerate and remove external JUBAT overrides. The shared setup then uses
    # its original defaults, including max_iter=100; all choices are saved below.
    external=Dict(k=>v for (k,v) in ENV if startswith(k,"JUBAT"))
    for key in keys(external)
        delete!(ENV,key)
    end
    head=source_identity(outdir)
    metadata=Dict("head"=>head,"julia_version"=>string(VERSION),
        "julia_threads"=>Threads.nthreads(),"external_JUBAT_removed"=>external,
        "OPENBLAS_NUM_THREADS"=>ENV["OPENBLAS_NUM_THREADS"],"GKSwstype"=>ENV["GKSwstype"],
        "phase_order"=>["discharge","rest1"],"requested_seconds"=>[1800.0,600.0],
        "case"=>2,"strength_scale"=>0.085,"stiffness_scale"=>0.1,"gc_scale"=>2.0,
        "dt_cycle_seconds"=>[1.0,10.0],"SOC_init"=>0.65,"current_discharge_A"=>5.0)
    open(io->TOML.print(io,metadata),joinpath(outdir,"run_metadata.toml"),"w")
    probe=JuBat.ChooseCell("Jellyroll");tc=TAB_CASES[2]
    h_tab=probe.cell.h*probe.cell.cooling_surface/(2*4e-3*L_BASE)
    setup=build_tab_case(tc,theta_pos_for(tc.theta_kind,spiral_theta_range(probe)),h_tab,outdir,nothing;
        n_theta=20,n_cycles=1,strength_scale=0.085,stiffness_scale=0.1,gc_scale=2.0,
        iter_method="gp_basic",tau_seconds="5")
    case=setup.case
    JuBat.apply_initial_soc!(case,case.param_dim,0.65)
    savebin(joinpath(outdir,"configuration.bin"),(param_dim=case.param_dim,param=case.param,opt=case.opt,
        cycle_opt=setup.cycle_opt,mesh=case.mesh,tab_case=tc))
    open(joinpath(outdir,"configuration.txt"),"w") do io
        println(io,"Baseline Case 2; tab.h=",h_tab)
        println(io,"Parameters/option values are stored without truncation in configuration.bin.")
        println(io,"Options:");dump(io,case.opt;maxdepth=8)
        println(io,"Dimensional parameters:");dump(io,case.param_dim;maxdepth=8)
    end
    # Exactly the fresh solve_cycling initial-state construction, without running
    # either of its later charge/rest2 phases.
    initial=(y=copy(JuBat.ModelInitialisation_MultiSPMe(case)),
        mech=deepcopy(setup.mechanical_seed),t_global=0.0)
    runtime=Dict{String,Any}("accepted"=>deepcopy(initial))
    export_state(outdir,"initial",initial,case)
    snaps=JuBat.CZMSnapshot[]
    current_label=Ref("discharge")
    pending=Ref{Any}(nothing)
    previous_observer=JuBat.mechanical_step_observer[]
    JuBat.mechanical_step_observer[]=function(event)
        if event.phase==:before
            pending[]=(phase=current_label[],ms=deepcopy(event.ms),F_ext=copy(event.F_ext),
                dT=deepcopy(event.dT_elem),delta_soc_n=deepcopy(event.Δsoc_n_elem),
                delta_soc_p=deepcopy(event.Δsoc_p_elem),dt_seconds=event.dt_seconds)
        elseif event.phase==:error
            savebin(joinpath(outdir,"failed_mechanical_input.bin"),pending[])
        end
    end
    summary=open(joinpath(outdir,"phase_summary.csv"),"w")
    println(summary,"phase,t_start_s,t_end_s,duration_s,V_start_V,V_end_V,capacity_Ah,termination_reason,T_max_K,T_mean_end_K,D_max_element,D_mean_element,D_max_gp,Dv_max_gp,n_fractured_gp,hoop_min_Pa,hoop_max_Pa,shear_min_Pa,shear_max_Pa,recovered_element_vm_max_history_Pa,n_recorded_times")
    try
        state=initial
        for (label,kind,duration) in (("discharge",JuBat.PHASE_DISCHARGE,1800.0),("rest1",JuBat.PHASE_REST,600.0))
            current_label[]=label
            println("HALFCYCLE phase=",label," start=",state.t_global," requested_s=",duration)
            flush(stdout)
            r=JuBat.solve_phase(case,kind,duration,state;dt_range=[1.0,10.0],
                czm_snapshots=snaps,czm_cycle=1,czm_phase=label,variables=runtime)
            savebin(joinpath(outdir,label*"_phase_result.bin"),r)
            state=r.final_state
            export_state(outdir,label,state,case)
            metrics=phase_metrics(case,label,r)
            println(summary,join(metrics,","));flush(summary)
            println("HALFCYCLE completed ",metrics);flush(stdout)
        end
        savebin(joinpath(outdir,"czm_snapshots.bin"),snaps)
        write(joinpath(outdir,"completion.txt"),"SUCCESS: discharge and rest1 only; no charge/rest2 run.\n")
    catch err
        savebin(joinpath(outdir,"last_accepted.bin"),runtime["accepted"])
        savebin(joinpath(outdir,"czm_snapshots.bin"),snaps)
        pending[]===nothing || savebin(joinpath(outdir,"last_mechanical_input.bin"),pending[])
        open(joinpath(outdir,"failure.txt"),"w") do io
            println(io,"phase=",current_label[])
            showerror(io,err,catch_backtrace());println(io)
        end
        rethrow()
    finally
        close(summary)
        JuBat.mechanical_step_observer[]=previous_observer
    end
end

main()
