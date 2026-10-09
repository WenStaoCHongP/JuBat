#=
任务62：独立进程并行重放已冻结的真实机械输入包。

正式入口（只重放，不重新采集，不修改原串行输出）：
  julia -t1 --startup-file=no --project=. example/run_mechanical_boundary_replays_parallel.jl
烟测入口（从完整包派生真实第1步，明确0圈）：
  julia -t1 --startup-file=no --project=. example/run_mechanical_boundary_replays_parallel.jl --smoke
可选：--packet=<已有完整输入包> --workers=4 --worker-threads=2

上限：4个独立worker，每个最多2 Julia线程、2 BLAS线程；无共享演化状态或缓存。
正式输出 output/run_mechanical_boundary_replays_parallel/；烟测输出其 parallel_smoke/。
冻结的 check_mechanical_boundary_branches.jl 和 src/ 均保持原样。

每个进程在反序列化输入包前，先用与采集会话相同的 build_collection_case 序列
重建一次案例（约6s）。输入包的 param 字段携带 JuBat 模块内匿名编号闭包
（#407#408 等），这些编号只在创建它们的方法体首次 lower 时分配；新进程若
不重放同序调用，Serialization 无法解析对应符号（UndefVarError）。
=#
include(joinpath(@__DIR__, "check_mechanical_boundary_branches.jl"))

const PARALLEL_REPLAY_OUTPUT = normpath(joinpath(@__DIR__, "..", "output", "run_mechanical_boundary_replays_parallel"))
const PARALLEL_ORIGINAL_PACKET = normpath(joinpath(@__DIR__, "..", "output", "check_mechanical_boundary_branches", "input_packet.bin"))
const PARALLEL_LAUNCHER = abspath(@__FILE__)

file_sha256(path) = bytes2hex(open(sha256, path))

function canonical_text!(io,text)
    bytes = codeunits(String(text))
    write(io,htol(UInt64(length(bytes))))
    write(io,bytes)
end

function canonical_value!(io,value)
    # Explicit scalar bytes avoid the padding of mixed isbits tuple arrays.
    if value===nothing
        write(io,UInt8(0))
    elseif value isa Bool
        write(io,UInt8(1));write(io,UInt8(value))
    elseif value isa Float64
        write(io,UInt8(2))
        bits = isnan(value) ? reinterpret(UInt64,NaN) : reinterpret(UInt64,value)
        write(io,htol(bits))
    elseif value isa Float32
        write(io,UInt8(3))
        bits = isnan(value) ? reinterpret(UInt32,Float32(NaN)) : reinterpret(UInt32,value)
        write(io,htol(bits))
    elseif value isa Integer && isbitstype(typeof(value))
        write(io,UInt8(4));canonical_text!(io,string(typeof(value)));write(io,htol(value))
    elseif value isa Symbol
        write(io,UInt8(5));canonical_text!(io,string(value))
    elseif value isa AbstractString
        write(io,UInt8(6));canonical_text!(io,value)
    elseif value isa NamedTuple
        write(io,UInt8(7));write(io,htol(UInt64(length(value))))
        for (name,item) in pairs(value)
            canonical_text!(io,string(name));canonical_value!(io,item)
        end
    elseif value isa Tuple
        write(io,UInt8(8));write(io,htol(UInt64(length(value))))
        for item in value
            canonical_value!(io,item)
        end
    elseif value isa AbstractArray
        write(io,UInt8(9));canonical_text!(io,string(typeof(value)))
        write(io,htol(UInt64(ndims(value))))
        for dimension in size(value)
            write(io,htol(UInt64(dimension)))
        end
        for item in value
            canonical_value!(io,item)
        end
    elseif isstructtype(typeof(value))
        write(io,UInt8(10));canonical_text!(io,string(typeof(value)))
        write(io,htol(UInt64(fieldcount(typeof(value)))))
        for name in fieldnames(typeof(value))
            canonical_text!(io,string(name));canonical_value!(io,getfield(value,name))
        end
    else
        error("unsupported canonical mechanical-state type $(typeof(value))")
    end
    return nothing
end

function canonical_state_sha256(ms)
    io = IOBuffer()
    canonical_text!(io,"JuBat mechanical state scalar encoding v1")
    canonical_value!(io,mechanical_state_identity(ms))
    return bytes2hex(sha256(take!(io)))
end

function parallel_canonical_audit(packet,outdir)
    rows = Any[]
    previous_hash = nothing
    equal_hash_pairs,equal_value_pairs = 0,0
    for (k,input) in enumerate(packet.inputs)
        before_hash = canonical_state_sha256(input.ms_before)
        after_hash = canonical_state_sha256(packet.baseline_after[k])
        values_equal = k==1 ? nothing : isequal(mechanical_state_identity(input.ms_before),
            mechanical_state_identity(packet.baseline_after[k-1]))
        hashes_equal = k==1 ? nothing : before_hash==previous_hash
        if k>1
            equal_hash_pairs += hashes_equal
            equal_value_pairs += values_equal
            hashes_equal && values_equal || error("canonical committed-state chain differs at input $k")
        end
        push!(rows,(k,input.phase,input.time_s,before_hash,after_hash,hashes_equal,values_equal))
        previous_hash = after_hash
    end
    write_csv(joinpath(outdir,"canonical_input_state_hashes.csv"),("step","phase","time_s","before_canonical_sha256",
        "committed_after_canonical_sha256","before_matches_previous_committed_canonical_hash","before_matches_previous_committed_fieldwise"),rows)
    return (;pairs=length(packet.inputs)-1,equal_hash_pairs,equal_value_pairs)
end

function parallel_arguments(args)
    packet, smoke, workers, worker_threads, worker_manifest, worker_id = PARALLEL_ORIGINAL_PACKET, false, 4, 2, nothing, nothing
    for arg in args
        if arg=="--smoke"
            smoke = true
        elseif startswith(arg,"--packet=")
            packet = abspath(split(arg,"=";limit=2)[2])
        elseif startswith(arg,"--workers=")
            workers = parse(Int,split(arg,"=";limit=2)[2])
        elseif startswith(arg,"--worker-threads=")
            worker_threads = parse(Int,split(arg,"=";limit=2)[2])
        elseif startswith(arg,"--worker-manifest=")
            worker_manifest = abspath(split(arg,"=";limit=2)[2])
        elseif startswith(arg,"--worker-id=")
            worker_id = parse(Int,split(arg,"=";limit=2)[2])
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    1<=workers<=4 || throw(ArgumentError("worker process limit is 1:4"))
    1<=worker_threads<=2 || throw(ArgumentError("each worker is limited to 1:2 Julia/BLAS threads"))
    (worker_manifest===nothing)==(worker_id===nothing) || error("worker manifest/id must be provided together")
    return (;packet,smoke,workers,worker_threads,worker_manifest,worker_id)
end

function parallel_plan(packet)
    candidates = candidate_matrix(packet)
    boundaries = [JuBat.resolve_mechanical_bc(packet.mesh,packet.param,c.options;endpoint_variant=Symbol(c.variant)) for c in candidates]
    keys = [serialized_hash((;operator=boundary_operator_identity(b),enforcement=enforcement_tag(c))) for (c,b) in zip(candidates,boundaries)]
    seen,owners = Dict{String,String}(),String[]
    for (c,key) in zip(candidates,keys)
        push!(owners,get!(seen,key,c.id))
    end
    unique_ids = [c.id for (c,owner) in zip(candidates,owners) if c.id==owner]
    length(candidates)==33 && length(unique_ids)==25 || error("expected 33 candidates and 25 distinct enforcement/constraint cases")
    return (;candidates,boundaries,keys,owners,unique_ids)
end

function parallel_derive_smoke(packet,parent_sha,outdir)
    packet.completed_cycles==1 && length(packet.completed_phases)==4 || error("smoke derivation requires a complete parent cycle")
    input = first(packet.inputs)
    config = merge(packet.config,(;max_accepted_steps=1,requested_window_s=nothing))
    derived = merge(packet,(;created_at=string(now()),config,config_hash=serialized_hash(config),
        inputs=Any[input],baseline_after=Any[first(packet.baseline_after)],baseline_results=Any[first(packet.baseline_results)],
        failed_input=nothing,last_committed=first(packet.baseline_after),actual_end_time_s=input.time_s,
        completed_cycles=0,completed_phases=Any[],final_chemical_thermal_state=nothing,
        stop_reason="derived_parallel_smoke_parent_steps_1_to_1",collection_error="",
        derivation=(;parent_packet_sha256=parent_sha,range_start=1,range_end=1,original_configuration_hash=packet.config_hash)))
    path = joinpath(outdir,"input_packet.bin")
    save_binary(path,derived)
    validate_packet_for_replay(derived,source_identity())
    derived.derivation.parent_packet_sha256 == parent_sha || error("derived packet parent SHA anchor mismatch")
    return derived,path
end

function parallel_export_input_metadata(packet,outdir)
    write_csv(joinpath(outdir,"source_identity.csv"),("source_relative_path","sha256"),((s.path,s.sha256) for s in packet.sources))
    write_csv(joinpath(outdir,"capture_phases.csv"),("phase","start_s","end_s","duration_s","current_A","terminated_by",
        "initial_chemical_thermal_hash","final_chemical_thermal_hash","case_initial_SOC_setting","reset_T_before_phase"),
        (Tuple(p) for p in packet.completed_phases))
    write_csv(joinpath(outdir,"input_steps.csv"),("step","cycle","phase","time_s","phase_time_s","dt_seconds",
        "F_ext_hash","dT_hash","delta_soc_n_hash","delta_soc_p_hash","before_state_hash","committed_state_hash",
        "GP_count","D_max_GP","D_max_element"),
        ((x.index,x.cycle,x.phase,x.time_s,x.phase_time_s,x.dt_seconds,serialized_hash(x.F_ext),serialized_hash(x.dT_elem),
            serialized_hash(x.Δsoc_n_elem),serialized_hash(x.Δsoc_p_elem),serialized_hash(mechanical_state_identity(x.ms_before)),
            serialized_hash(mechanical_state_identity(packet.baseline_after[k])),length(packet.baseline_after[k].gp_damage_states),
            maximum(s.D for s in packet.baseline_after[k].gp_damage_states),maximum(s.D for s in packet.baseline_after[k].damage_states))
            for (k,x) in enumerate(packet.inputs)))
end

function warmup_jubat_closures(outdir)
    scratch = joinpath(outdir,"warmup_case_scratch")
    mkpath(scratch)
    build_collection_case(scratch)
    return nothing
end

function validate_packet_for_replay(packet, sources)
    # 与冻结版 validate_packet 逐项相同，仅豁免第303行 config 字节哈希一项。
    # 该哈希覆盖 czm_options::Dict{String,String}，序列化字节取决于字典布局；
    # 反序列化进程按流中迭代序重建，布局与采集进程的 fieldnames 插入序不同，
    # 跨进程不可复现（探针：param/mesh 字节稳定、dict 内容相等但字节不等）。
    # 内容级替代：czm_options 与 options_identity(packet.options) 字典相等；
    # 包级身份由文件 SHA256 锚定（controller 对照 identity.toml，worker 对照 manifest）。
    packet.schema == "JuBat mechanical boundary packet v1" || error("unsupported packet schema")
    packet.source_hash == serialized_hash(sources) || error("packet source identity differs from current source")
    packet.mesh_hash == serialized_hash(static_mesh_identity(packet.mesh)) || error("packet mesh hash mismatch")
    packet.param_hash == serialized_hash(packet.param) || error("packet param hash mismatch")
    packet.config.czm_options == options_identity(packet.options) ||
        error("packet config czm_options is inconsistent with packet options")
    length(packet.inputs) == length(packet.baseline_after) == length(packet.baseline_results) ||
        error("packet input/state/result lengths differ")
    isempty(packet.inputs) && error("no accepted mechanical inputs were captured")
    for (k,x) in enumerate(packet.inputs)
        x.index == k || error("noncontiguous packet indices")
        k == 1 || isequal(mechanical_state_identity(x.ms_before),
            mechanical_state_identity(packet.baseline_after[k-1])) || error("baseline committed-state chain differs at step $k")
        isfinite(x.dt_seconds) && x.dt_seconds > 0 || error("invalid packet dt at step $k")
        all(isfinite, x.F_ext) && all(isfinite, x.dT_elem) &&
            all(isfinite, x.Δsoc_n_elem) && all(isfinite, x.Δsoc_p_elem) || error("nonfinite packet input at step $k")
        packet.baseline_after[k].gp_damage_states === nothing && error("accepted gp_basic state lacks actual GP history")
    end
    return nothing
end

function verify_packet_file_anchor(packet_path)
    identity_path = joinpath(dirname(packet_path), "identity.toml")
    isfile(identity_path) || error("identity.toml not found next to packet; cannot anchor frozen file identity")
    recorded = TOML.parsefile(identity_path)
    recorded["packet_sha256"] == file_sha256(packet_path) ||
        error("packet file SHA256 differs from identity.toml record")
    return recorded
end

function parallel_worker(manifest_path,worker_id)
    manifest = TOML.parsefile(manifest_path)
    file_sha256(PARALLEL_LAUNCHER)==manifest["launcher_sha256"] || error("parallel launcher changed after manifest creation")
    Threads.nthreads()==manifest["worker_threads"] || error("worker Julia thread count differs from manifest")
    BLAS.set_num_threads(manifest["worker_threads"])
    BLAS.get_num_threads()==manifest["worker_threads"] || error("worker BLAS thread count differs from manifest")
    packet_path,outdir = manifest["packet_path"],manifest["output_directory"]
    file_sha256(packet_path)==manifest["packet_sha256"] || error("worker input packet SHA mismatch")
    warmup_jubat_closures(outdir)
    packet = open(deserialize,packet_path)
    validate_packet_for_replay(packet,source_identity())
    plan = parallel_plan(packet)
    ids = manifest["assignments"][string(worker_id)]
    all(id->id in plan.unique_ids,ids) && length(unique(ids))==length(ids) || error("invalid worker case assignment")
    theta = unfolded_angles(packet.mesh,packet.param)
    by_id = Dict(c.id=>i for (i,c) in enumerate(plan.candidates))
    results,unexpected = Any[],Any[]
    for id in ids
        i = by_id[id]
        try
            push!(results,replay_branch(plan.candidates[i],plan.boundaries[i],packet,theta,outdir))
        catch err
            message = sprint(showerror,err,catch_backtrace())
            push!(unexpected,(id,message))
            branch_dir = joinpath(outdir,"branches",id)
            mkpath(branch_dir)
            write(joinpath(branch_dir,"worker_exception.txt"),message*"\n")
            println("WORKER=$worker_id CASE=$id UNEXPECTED_EXCEPTION=$message")
        end
    end
    result_path = joinpath(outdir,"workers","worker_$(worker_id)_result.bin")
    save_binary(result_path,(;worker_id,pid=getpid(),execution_id=manifest["execution_id"],launcher_sha256=manifest["launcher_sha256"],
        julia_threads=Threads.nthreads(),blas_threads=BLAS.get_num_threads(),
        packet_sha256=manifest["packet_sha256"],source_hash=packet.source_hash,assigned_ids=ids,results,unexpected))
    success = isempty(unexpected) && length(results)==length(ids) && all(r->r.status=="completed",results)
    println("WORKER=$worker_id FINISHED=$(length(results))/$(length(ids)) SUCCESS=$success")
    return success
end

function parallel_summary(packet,plan,results,outdir)
    by_id = Dict(r.id=>r for r in results)
    rows = Any[]
    for (candidate,owner,key) in zip(plan.candidates,plan.owners,plan.keys)
        if !haskey(by_id,owner)
            push!(rows,(candidate.id,owner,candidate.id!=owner,key,candidate.outer,candidate.start,candidate.finish,
                candidate.variant,"missing_worker_result",nothing,length(packet.inputs),nothing,nothing,nothing,nothing,nothing,
                nothing,nothing,nothing,nothing,nothing,nothing,nothing,nothing,"worker did not return a valid committed result"))
            continue
        end
        r = by_id[owner]
        ds = r.diagnostics
        fs = r.fields
        gpmax = isempty(fs) ? nothing : maximum(f.D for f in fs)
        dvmax = isempty(fs) ? nothing : maximum(f.D_visc for f in fs)
        weighted = isempty(fs) ? nothing : sum(f.D*f.length_weight_m for f in fs)/sum(f.length_weight_m for f in fs)
        push!(rows,(candidate.id,owner,candidate.id!=owner,key,candidate.outer,candidate.start,candidate.finish,candidate.variant,
            r.status,r.accepted,length(packet.inputs),r.end_time,gpmax,dvmax,maximum(s.D for s in r.ms.damage_states),weighted,
            isempty(fs) ? nothing : count(f->f.D>1e-8,fs),isempty(fs) ? nothing : count(f->f.D>=0.99,fs),
            isempty(fs) ? nothing : count(f->f.fractured,fs),r.last_residual,ds===nothing ? nothing : ds.free_residual,
            ds===nothing ? nothing : ds.constraint_violation,ds===nothing ? nothing : ds.gauge_moment,r.wall_seconds,r.error))
    end
    write_csv(joinpath(outdir,"summary.csv"),("candidate","unique_owner","duplicate_solve_case","solve_identity_hash",
        "outer_bc","fix_start_raw","fix_end_raw","endpoint_variant","status","accepted_steps","requested_input_steps",
        "last_committed_time_s","D_max_GP","D_visc_max_GP","D_max_element_mean","D_length_weighted_mean",
        "damaged_GP_count_D_gt_1e_8","GP_count_D_ge_0_99","fractured_GP_count","last_solver_residual",
        "last_committed_free_residual","constraint_violation","gauge_moment","unique_replay_wall_seconds","error"),rows)
    return rows
end

function parallel_controller(settings)
    outdir = settings.smoke ? joinpath(PARALLEL_REPLAY_OUTPUT,"parallel_smoke") : PARALLEL_REPLAY_OUTPUT
    mkpath(joinpath(outdir,"workers"))
    warmup_jubat_closures(outdir)
    parent_sha = file_sha256(settings.packet)
    packet = open(deserialize,settings.packet)
    validate_packet_for_replay(packet,source_identity())
    identity_record = verify_packet_file_anchor(settings.packet)
    original_configuration_hash = packet.config_hash
    packet_path = settings.packet
    if settings.smoke
        packet,packet_path = parallel_derive_smoke(packet,parent_sha,outdir)
        GC.gc()
    else
        packet.completed_cycles==1 && length(packet.completed_phases)==4 &&
            isempty(packet.collection_error) && !isempty(packet.inputs) ||
            error("formal parallel validation requires one complete collected cycle with all four phases and no collection error")
    end
    packet_sha = file_sha256(packet_path)
    plan = parallel_plan(packet)
    audit = parallel_canonical_audit(packet,outdir)
    assignments = Dict(string(i)=>plan.unique_ids[i:settings.workers:end] for i in 1:settings.workers)
    # Conservative launch gate, not a claim of an exact peak-memory prediction.
    reserve_bytes = settings.workers*(3*filesize(packet_path)+512*1024^2)
    Sys.free_memory()>=reserve_bytes || error("insufficient free RAM for bounded worker launch: need conservative reserve $reserve_bytes bytes")
    execution_id = string(now())*"-pid"*string(getpid())
    manifest = Dict("schema"=>"JuBat independent mechanical replay execution v1","created_at"=>string(now()),"execution_id"=>execution_id,
        "launcher_path"=>PARALLEL_LAUNCHER,"launcher_sha256"=>file_sha256(PARALLEL_LAUNCHER),
        "packet_path"=>abspath(packet_path),"packet_sha256"=>packet_sha,"parent_packet_path"=>abspath(settings.packet),
        "parent_packet_sha256"=>parent_sha,"identity_toml_packet_sha256"=>identity_record["packet_sha256"],
        "source_sha256"=>packet.source_hash,
        "configuration_sha256"=>packet.config_hash,"parent_configuration_sha256"=>original_configuration_hash,
        "mesh_sha256"=>packet.mesh_hash,"param_sha256"=>packet.param_hash,
        "smoke"=>settings.smoke,"derived_range_start"=>1,"derived_range_end"=>settings.smoke ? 1 : length(packet.inputs),
        "accepted_input_steps"=>length(packet.inputs),"actual_end_time_s"=>packet.actual_end_time_s,"completed_cycles"=>packet.completed_cycles,
        "canonical_encoding_schema"=>"JuBat mechanical state scalar encoding v1: type tags, container shape, individual little-endian scalars, Bool byte, length-prefixed UTF8",
        "canonical_committed_chain_pairs"=>audit.pairs,"canonical_committed_chain_hash_equal_pairs"=>audit.equal_hash_pairs,
        "direct_fieldwise_committed_chain_equal_pairs"=>audit.equal_value_pairs,
        "workers"=>settings.workers,"worker_threads"=>settings.worker_threads,"worker_blas_threads"=>settings.worker_threads,
        "max_simultaneous_worker_julia_threads"=>settings.workers*settings.worker_threads,
        "max_simultaneous_worker_blas_threads"=>settings.workers*settings.worker_threads,
        "conservative_memory_reserve_bytes"=>reserve_bytes,"free_memory_before_launch_bytes"=>Sys.free_memory(),
        "output_directory"=>outdir,"assignments"=>assignments)
    manifest_path = joinpath(outdir,"execution_manifest.toml")
    open(io->TOML.print(io,manifest),manifest_path,"w")
    parallel_export_input_metadata(packet,outdir)
    theta = unfolded_angles(packet.mesh,packet.param)
    export_geometry(packet.mesh,packet.param,theta,outdir)
    export_candidate_boundaries(plan.candidates,plan.boundaries,plan.owners,plan.keys,packet.mesh,packet.param,outdir)
    processes,logs,pids = Any[],IO[],Int[]
    executable = joinpath(Sys.BINDIR,Base.julia_exename())
    for i in 1:settings.workers
        log_path = joinpath(outdir,"workers","worker_$i.log")
        io = open(log_path,"w")
        cmd = Cmd([executable,"-t$(settings.worker_threads)","--startup-file=no","--project=$BC_DIAGNOSTIC_ROOT",
            PARALLEL_LAUNCHER,"--worker-manifest=$manifest_path","--worker-id=$i"])
        cmd = addenv(cmd,"OPENBLAS_NUM_THREADS"=>string(settings.worker_threads),"GKSwstype"=>"100")
        process = run(pipeline(cmd;stdout=io,stderr=io);wait=false)
        pid = getpid(process)
        push!(processes,process);push!(logs,io);push!(pids,pid)
        println("WORKER_STARTED=$i PID=$pid CASES=$(length(assignments[string(i)])) JULIA_THREADS=$(settings.worker_threads) BLAS_THREADS=$(settings.worker_threads)")
    end
    done = falses(settings.workers)
    while !all(done)
        for i in 1:settings.workers
            if !done[i] && process_exited(processes[i])
                done[i] = true
                close(logs[i])
                println("WORKER_EXITED=$i EXIT_CODE=$(processes[i].exitcode)")
            end
        end
        all(done) || sleep(1)
    end
    results,status_rows = Any[],Any[]
    for i in 1:settings.workers
        result_path = joinpath(outdir,"workers","worker_$(i)_result.bin")
        returned,invalid_reason = nothing,""
        try
            candidate = isfile(result_path) ? open(deserialize,result_path) : nothing
            if candidate!==nothing
                candidate.worker_id==i && candidate.pid==pids[i] && candidate.execution_id==execution_id &&
                    candidate.launcher_sha256==manifest["launcher_sha256"] && candidate.packet_sha256==packet_sha &&
                    candidate.source_hash==packet.source_hash && candidate.assigned_ids==assignments[string(i)] &&
                    candidate.julia_threads==settings.worker_threads && candidate.blas_threads==settings.worker_threads ||
                    error("worker $i result identity/pid/thread/assignment mismatch")
                returned = candidate
            end
        catch err
            invalid_reason = sprint(showerror,err)
        end
        if returned!==nothing
            append!(results,returned.results)
        end
        push!(status_rows,(i,pids[i],processes[i].exitcode,length(assignments[string(i)]),
            returned===nothing ? 0 : length(returned.results),returned===nothing ? "missing_or_invalid_result: "*invalid_reason :
                join((string(x[1])*": "*string(x[2]) for x in returned.unexpected)," | "),file_sha256(joinpath(outdir,"workers","worker_$i.log"))))
    end
    length(unique(r.id for r in results))==length(results) || error("duplicate independent worker result")
    order = Dict(id=>i for (i,id) in enumerate(plan.unique_ids))
    sort!(results;by=r->get(order,r.id,typemax(Int)))
    all(r->r.id in plan.unique_ids,results) || error("worker returned an unassigned case")
    write_csv(joinpath(outdir,"worker_status.csv"),("worker","pid","exit_code","assigned_cases","returned_cases","unexpected_errors","log_sha256"),status_rows)
    rows = parallel_summary(packet,plan,results,outdir)
    shared_field_plots(results,outdir)
    complete = length(results)==length(plan.unique_ids)
    report_ok = if complete
        diagnostic_report(packet,packet_sha,plan.candidates,rows,results,outdir)
    else
        write(joinpath(outdir,"report.md"),"# Incomplete parallel replay\n\nReturned $(length(results))/$(length(plan.unique_ids)) independent cases. Missing worker results are marked in summary.csv; no committed state was invented.\n")
        false
    end
    open(joinpath(outdir,"report.md"),"a") do io
        println(io,"\nEvery process (controller and each worker) replayed the collection case setup once before deserialization, so the packet's anonymously numbered JuBat closures (e.g. #407#408, first observed as a fresh-process deserialization UndefVarError) exist with identical numbering. The warmup writes only to warmup_case_scratch/ and does not touch frozen collection outputs or source identity.")
        println(io,"\nPacket validation in this launcher keeps the frozen source/mesh/param byte checks (probe-verified cross-process stable) and replaces the config byte-hash with content-level checks: czm_options must equal options_identity(packet.options), and packet identity is anchored by the packet file SHA256 against identity.toml (controller) or the execution manifest (workers). Reason: the packet's stored config_hash covers czm_options::Dict{String,String}, whose serialized bytes depend on the dict layout; a fresh process rebuilds the dict from stream iteration order and cannot reproduce the collecting process's layout even though the content is equal.")
        println(io,"\nIndependent process execution: $(settings.workers) workers, each $(settings.worker_threads) Julia/BLAS threads. Core and original diagnostic source were not edited. The launcher has its own execution_manifest.toml SHA identity; the input packet's original source identity remains unchanged.")
        println(io,"\nInput packet: [saved input](<$(abspath(packet_path))>). Parent baseline packet: [original complete packet](<$(abspath(settings.packet))>). Parent SHA-256: `$parent_sha`.")
        println(io,"\nDerivation range: 1:$(settings.smoke ? 1 : length(packet.inputs)); smoke=$(settings.smoke). Original configuration hash: `$original_configuration_hash`. Smoke contains only the actual first input and its actual committed GP state, no future chemical/thermal state; it proves launcher/aggregation contracts and reports 0 complete cycles.")
        println(io,"\nThe input packet is referenced at the path above; a full packet is not copied into this output directory. All case outputs belong to this independent run. Worker exits/log hashes/assignments are in worker_status.csv and execution_manifest.toml. Shared plots are generated after all workers return using the original common-range plotting function.")
        println(io,"\nState verification uses validate_packet_for_replay's direct fieldwise committed-state chain comparison (identical to the frozen collector's validate_packet except the exempted config byte-hash) and the legacy replay's per-step exact comparison. The original before/after state CSV SHA values hash Julia Serialization byte representations, which can contain padding; unequal state byte hashes alone do not establish unequal engineering values. The saved input packet's file SHA-256 is the authoritative byte identity for the common input.")
        println(io,"\nSupplemental canonical_input_state_hashes.csv encodes every mechanical state field using explicit type tags, container shape and individual scalar bytes. Canonical hash chain: $(audit.equal_hash_pairs)/$(audit.pairs) equal pairs; direct fieldwise chain: $(audit.equal_value_pairs)/$(audit.pairs) equal pairs. The launcher never serializes a whole tuple, struct or array to obtain these canonical hashes.")
        println(io,"\nThe collector uses Newton-Cotes endpoint GP at xi=-1,+1. Adjacent elements can place distinct GP histories at the same physical location. GP counts and physical hotspot-location counts must remain separate; quadrature-weighted damaged support length is not a measured or resolved crack-front extent.")
    end
    success = report_ok && complete && all(p->p.exitcode==0,processes) && all(r->r.status=="completed",results)
    completed = count(r->r.status=="completed",results)
    println("PARALLEL_REPLAY_SUCCESS=$success CANDIDATES=$(length(plan.candidates)) UNIQUE_CASES=$(length(plan.unique_ids)) RETURNED=$(length(results)) COMPLETED=$completed")
    return success
end

function parallel_main(args=ARGS)
    settings = parallel_arguments(args)
    return settings.worker_manifest===nothing ? parallel_controller(settings) : parallel_worker(settings.worker_manifest,settings.worker_id)
end

if abspath(PROGRAM_FILE)==PARALLEL_LAUNCHER
    parallel_main() || exit(1)
end
