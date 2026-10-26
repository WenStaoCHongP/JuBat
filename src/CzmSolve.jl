include("CzmSolveCommon.jl")
include("CzmSolveBasic.jl")
include("CzmSolveLoadSubstep.jl")
include("CzmSolveArcLength.jl")

"""
    solve_czm_step(czm_mesh, ms, param, load::MechanicalLoad, czm_opt; czm_workspace, ...)

load.F_ext 为二维单位高度节点力的无量纲向量：物理 N/m 除 σ_ref*L；若给定高度 H
上的物理总节点力 N，则除 σ_ref*L*H。坐标与位移除 L。load.dt_seconds 及 viscous_tau
均为实际物理秒；load.sample_time_seconds 为全局物理采样秒数，不是候选发布时刻。
此接口不按数值大小猜测或转换单位，纯外力调用须显式提供三组零热化场。
"""
function solve_czm_step(czm_mesh::CohesiveMesh, ms::MechState, param,
        load::MechanicalLoad, czm_opt::CzmOptions; czm_workspace::CZMWorkspace,
        prestress=nothing, boundary=nothing)
    length(ms.u_prev)==2size(czm_mesh.node,1) || throw(DimensionMismatch("CZM displacement length"))
    length(ms.damage_states)==length(czm_mesh.cohesive_elements) || throw(DimensionMismatch("element damage history length"))
    czm_opt.winding_prestress && prestress===nothing && throw(ArgumentError(
        "winding_prestress=true requires explicit prestress input"))
    check_czm_workspace(czm_mesh, czm_workspace)
    czm_workspace.param === param || throw(ArgumentError("CZM workspace belongs to different frozen parameters"))
    czm_opt.outer_bc === :radial_slide && czm_opt.geo_nonlinear && throw(ArgumentError(
        "radial_slide requires geo_nonlinear=false"))
    if boundary === nothing && (czm_opt.outer_bc !== :fixed_xy ||
            czm_opt.fix_start !== nothing || czm_opt.fix_end !== nothing ||
            czm_opt.endpoint_variant !== :keep)
        boundary = resolve_mechanical_bc(czm_mesh,param,czm_opt)
    end
    if boundary !== nothing
        boundary.coordinates == czm_mesh.node || throw(ArgumentError(
            "mechanical boundary coordinates differ from solver mesh"))
        czm_opt.geo_nonlinear && any(==(:radial),values(boundary.nodes)) &&
            throw(ArgumentError("radial_slide requires geo_nonlinear=false"))
    end
    observer = mechanical_step_observer[]
    if observer === nothing
        result = _solve_czm_step_impl(czm_mesh,ms,param,load,czm_opt; czm_workspace=czm_workspace,
            prestress=prestress,boundary=boundary)
        result.converged && (ms.mechanics_initialized=true)
        return result
    end
    event = (;czm_mesh,ms,param,load,czm_opt,prestress)
    observer(merge(event,(phase=:before,result=nothing,error=nothing)))
    result = try
        _solve_czm_step_impl(czm_mesh,ms,param,load,czm_opt; czm_workspace=czm_workspace,
            prestress=prestress,boundary=boundary)
    catch err
        observer(merge(event,(phase=:error,result=nothing,error=err)))
        rethrow()
    end
    result.converged && (ms.mechanics_initialized=true)
    observer(merge(event,(phase=:after,result=result,error=nothing)))
    return result
end

"""
    czm_viscous_beta(czm_opt, dt_seconds)

按两次 CZM 更新间实际经过的物理秒数计算后向 Euler 粘性松弛系数。
"""
function czm_viscous_beta(czm_opt::CzmOptions, dt_seconds::Union{Nothing, Real})
    czm_opt.viscous_enabled || return 1.0
    τ = czm_opt.viscous_tau
    isfinite(τ) && τ >= 0.0 || throw(ArgumentError(
        "CZM viscous_tau must be finite and nonnegative seconds, got $τ"))
    τ == 0.0 && return 1.0
    Δt = Float64(dt_seconds)
    isfinite(Δt) && Δt > 0.0 || throw(ArgumentError(
        "CZM viscosity requires finite positive elapsed seconds, got $Δt"))
    return Δt / (τ + Δt)
end

function czm_snapshot_source_state()
    repo_root = normpath(joinpath(@__DIR__, ".."))
    return (
        git_head=readchomp(`git -c core.excludesFile= -C $repo_root rev-parse HEAD`),
        source_status=read(`git -c core.excludesFile= -c core.autocrlf=false -C $repo_root status --short -- src`, String),
        source_diff=read(`git -c core.excludesFile= -c core.autocrlf=false -C $repo_root diff -- src`, String),
        project_toml=read(joinpath(repo_root, "Project.toml"), String),
        julia_version=string(VERSION),
    )
end

function save_czm_step_snapshot(path::String, czm_mesh::CohesiveMesh,
        ms_before::MechState, param::Params, load::MechanicalLoad,
        czm_opt::CzmOptions)
    payload = (
        format_version=2,
        source=czm_snapshot_source_state(),
        czm_mesh=czm_mesh,
        mech=ms_before,
        param=param,
        load=deepcopy(load),
        czm_opt=deepcopy(czm_opt),
    )
    mkpath(dirname(path))
    temp_path = path * ".tmp"
    open(temp_path, "w") do io
        serialize(io, payload)
    end
    mv(temp_path, path; force=true)
    return path
end

function replay_czm_step_snapshot(path::String)
    snapshot = open(deserialize, path)
    snapshot.format_version == 2 || throw(ArgumentError(
        "unsupported CZM snapshot format $(snapshot.format_version)"))
    snapshot.source == czm_snapshot_source_state() || throw(ArgumentError(
        "CZM snapshot source differs from the current source tree"))
    ms = deepcopy(snapshot.mech)
    czm_opt = deepcopy(snapshot.czm_opt)
    czm_opt.snapshot_path = ""
    czm_workspace = CZMWorkspace(snapshot.czm_mesh, snapshot.param)
    result = solve_czm_step(snapshot.czm_mesh, ms, snapshot.param,
        deepcopy(snapshot.load), czm_opt; czm_workspace=czm_workspace)
    return result, ms
end

"""
    solve_czm_step(czm_mesh, ms, param, load, czm_opt; czm_workspace, ...) -> CZMResult

CZM 单步统一入口。演化状态在 `ms::MechState` 上收敛提交；失败不触碰。
"""
function _solve_czm_step_impl(czm_mesh::CohesiveMesh, ms::MechState, param, load::MechanicalLoad,
        czm_opt::CzmOptions; czm_workspace::CZMWorkspace,
        prestress=nothing, boundary=nothing)
    method = lowercase(czm_opt.iter_method)
    gp_arc = method in ("arc_length", "arclength", "arc-length") &&
        !czm_opt.geo_nonlinear
    if czm_opt.viscous_enabled && czm_opt.viscous_tau > 0.0 && method != "gp_basic"
        throw(ArgumentError("Physical-time CZM viscosity requires gp_basic so damage relaxes once per accepted time step"))
    end
    max_iter = czm_opt.max_iter
    tol = czm_opt.tol
    n_load_steps = czm_opt.load_steps
    arc_length_alpha = czm_opt.arc_length_alpha
    visc_beta = czm_viscous_beta(czm_opt, load.dt_seconds)
    czm_model = czm_opt.model
    fix_inner = czm_opt.fix_inner
    geo_nl = czm_opt.geo_nonlinear

    snapshot_enabled = !isempty(czm_opt.snapshot_path)
    snapshot_enabled && !gp_arc && throw(ArgumentError(
        "CZM step snapshots currently require nongeometric arc_length"))
    ms_before = if snapshot_enabled
        deepcopy(ms)
    else
        nothing
    end
    load_before = if snapshot_enabled
        deepcopy(load)
    else
        nothing
    end

    result = try
        if method == "load_substep"
            newton_raphson_czm(czm_mesh, load, param, ms; czm_workspace=czm_workspace,
                max_iter=max_iter, tol=tol, n_load_steps=n_load_steps,
                visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, boundary=boundary, geo_nl=geo_nl,
                plasticity=czm_opt.j2_plasticity, prestress=prestress)
        elseif method == "basic" || method == "gp_basic"
            solve_czm_basic_step(czm_mesh, load, param, ms; czm_workspace=czm_workspace,
                max_iter=max_iter, tol=tol,
                visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, boundary=boundary, geo_nl=geo_nl,
                plasticity=czm_opt.j2_plasticity, prestress=prestress)
        elseif (method == "arc_length" || method == "arclength" || method == "arc-length") && geo_nl
            solve_czm_arc_geo_step(czm_mesh, load, param, ms; czm_workspace=czm_workspace,
                max_iter=max_iter, tol=tol, n_load_steps=n_load_steps,
                arc_length_alpha=arc_length_alpha,
                visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, boundary=boundary,
                plasticity=czm_opt.j2_plasticity, prestress=prestress)
        elseif method == "arc_length" || method == "arclength" || method == "arc-length"
            solve_czm_arc_length_step(czm_mesh, load, param, ms; czm_workspace=czm_workspace,
                max_iter=max_iter, tol=tol, n_load_steps=n_load_steps, arc_length_alpha=arc_length_alpha,
                visc_beta=visc_beta, czm_model=czm_model, fix_inner=fix_inner, boundary=boundary,
                wall_time_limit_seconds=czm_opt.wall_time_limit_seconds)
        else
            error("Unknown CZM iteration method: $(czm_opt.iter_method). Use 'basic', 'gp_basic', 'load_substep', or 'arc_length'.")
        end
    catch err
        if snapshot_enabled
            try
                save_czm_step_snapshot(czm_opt.snapshot_path, czm_mesh,
                    ms_before, param, load_before, czm_opt)
            catch snapshot_err
                error("CZM solve failed: $(sprint(showerror, err)); snapshot failed: $(sprint(showerror, snapshot_err))")
            end
        end
        rethrow()
    end
    if snapshot_enabled && !result.converged
        save_czm_step_snapshot(czm_opt.snapshot_path, czm_mesh,
            ms_before, param, load_before, czm_opt)
    end
    return result
end
