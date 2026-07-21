# tools/czm_convergence_diag.jl
# CZM 不收敛根因诊断
#
# 运行方式: julia tools/czm_convergence_diag.jl

using LinearAlgebra
using SparseArrays
using Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function run_diagnostics()
    println("="^60)
    println("CZM Convergence Root-Cause Diagnostics")
    println("="^60)

    # 1. 参数和网格
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.gsorder = 2
    opt.czm_model = "model1"

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=40, gsorder=2)
    czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)

    ne = size(czm_mesh.bulk_element, 1)
    ndof = 2 * czm_mesh.nnode
    n_coh = czm_mesh.n_cohesive

    println("\n--- Mesh Info ---")
    println("  Nodes: $(czm_mesh.nnode), DOFs: $ndof")
    println("  Bulk Elements: $ne, Cohesive Elements: $n_coh")

    E_eff, ν_eff, α_eff, β_n, β_p = JuBat.compute_czm_effective_params(case)
    @printf("  E_eff = %.4e, ν_eff = %.4f\n", E_eff, ν_eff)
    @printf("  α_eff = %.4e, β_n = %.4e, β_p = %.4e\n", α_eff, β_n, β_p)

    # 2. 检查归一化尺度
    println("\n--- Normalization Scales ---")
    sc = param_dim.scale
    @printf("  σ_czm (stress scale)  = %.4e Pa\n", sc.σ_czm)
    @printf("  δ_czm (length scale)  = %.4e m\n", sc.δ_czm)
    @printf("  K_czm (stiffness scale)= %.4e Pa/m\n", sc.K_czm)
    @printf("  E_n (NE energy scale) = %.4e Pa\n", sc.E_n)
    @printf("  E_p (PE energy scale) = %.4e Pa\n", sc.E_p)
    @printf("  L   (length scale)    = %.4e m\n", sc.L)

    # 3. 检查 cohesive 参数
    coh = param_dim.cohesive
    # δ_0_n = coh.δ_0_n                                                    # TODO Chunk 2 Task 2.1
    # δ_c_n = coh.δ_c_n                                                    # TODO Chunk 2 Task 2.1
    # K_n = coh.K_n                                                        # TODO Chunk 2 Task 2.1
    @printf("\n--- Cohesive Parameters (physical) ---\n")
    # @printf("  δ_0_n = %.4e m  (damage initiation)\n", δ_0_n)             # TODO Chunk 2 Task 2.1
    # @printf("  δ_c_n = %.4e m  (complete failure)\n", δ_c_n)              # TODO Chunk 2 Task 2.1
    # @printf("  K_n   = %.4e Pa/m (penalty stiffness)\n", K_n)             # TODO Chunk 2 Task 2.1
    # @printf("  σ_max = %.4e Pa\n", coh.σ_max_n)                           # TODO Chunk 2 Task 2.1

    # 4. 检查归一化后参数
    coh_norm = case.param.cohesive
    @printf("\n--- Cohesive Parameters (normalized) ---\n")
    # @printf("  δ_0_n* = %.4e\n", coh_norm.δ_0_n)                          # TODO Chunk 2 Task 2.1
    # @printf("  δ_c_n* = %.4e\n", coh_norm.δ_c_n)                          # TODO Chunk 2 Task 2.1
    # @printf("  K_n*   = %.4e\n", coh_norm.K_n)                            # TODO Chunk 2 Task 2.1
    # @printf("  σ_max* = %.4e\n", coh_norm.σ_max_n)                        # TODO Chunk 2 Task 2.1

    # 5. 检查 NE/PE 归一化弹性模量
    p = case.param
    @printf("\n--- Electrode Elastic Moduli ---\n")
    @printf("  NE.E* = %.4e (normalized by E_n=%.4e)\n", p.NE.E, sc.E_n)
    @printf("  PE.E* = %.4e (normalized by E_p=%.4e)\n", p.PE.E, sc.E_p)
    @printf("  E_eff* = %.4e\n", E_eff)
    @printf("  E_eff (physical estimate) = %.4e Pa\n", E_eff * sc.σ_czm)

    # 6. 小扰动载荷
    dT_elem = fill(0.05 / param_dim.cell.T0, ne)
    Δsoc_n_elem = fill(0.0002, ne)
    Δsoc_p_elem = fill(0.0002, ne)

    # 计算初始应变
    ε_0_per_elem = α_eff .* dT_elem .+ β_n .* Δsoc_n_elem .+ β_p .* Δsoc_p_elem
    @printf("\n--- Initial Strain (per element) ---\n")
    @printf("  ε_0_thermal  = %.4e\n", α_eff * dT_elem[1])
    @printf("  ε_0_diff_n   = %.4e\n", β_n * Δsoc_n_elem[1])
    @printf("  ε_0_diff_p   = %.4e\n", β_p * Δsoc_p_elem[1])
    @printf("  ε_0_total    = %.4e\n", ε_0_per_elem[1])

    # 7. 组装系统并检查（单次线性求解，D=0）
    cache = JuBat.build_czm_cache(czm_mesh, E_eff, ν_eff, case.param)
    bc_dofs = cache.bc_dofs
    bc_vals = cache.bc_vals

    F_thermo_chem = JuBat.assemble_thermal_chemical_load(
        czm_mesh, E_eff, ν_eff, α_eff, β_n, β_p,
        dT_elem, Δsoc_n_elem, Δsoc_p_elem)
    F_ext = zeros(Float64, ndof)

    u0 = zeros(Float64, ndof)
    damage_states = czm_mesh.damage_states

    K_total, f_int_total, separations, tractions = JuBat.assemble_coupled_system(
        czm_mesh, u0, E_eff, ν_eff, coh_norm;
        damage_states=damage_states, K_bulk_cached=cache.K_bulk,
        geom_cache=cache.cohesive_geom, ws=cache.ws)

    R = F_ext + F_thermo_chem - f_int_total
    for (dof, val) in zip(bc_dofs, bc_vals)
        R[dof] = val - u0[dof]
    end

    @printf("\n--- Initial System (u=0, D=0) ---\n")
    @printf("  ||F_thermo_chem|| = %.6e\n", norm(F_thermo_chem))
    @printf("  ||f_int(u=0)||    = %.6e\n", norm(f_int_total))
    @printf("  ||R|| (initial)   = %.6e\n", norm(R))

    # 8. 刚度矩阵条件数估计
    K_bc, R_bc = JuBat.apply_bc_czm(K_total, R; bc_dofs=bc_dofs, bc_vals=bc_vals)
    @printf("\n--- Stiffness Matrix ---\n")
    @printf("  K_total size: %d×%d, nnz: %d\n", size(K_total, 1), size(K_total, 2), nnz(K_total))
    @printf("  K_bulk  nnz: %d\n", nnz(cache.K_bulk))
    @printf("  K cohesive contribution norm: %.6e\n", norm(K_total - cache.K_bulk, 1))

    # 检查刚度矩阵的对称性和对角占优
    K_dense = Matrix(K_bc)
    K_asym = norm(K_dense - K_dense', 1) / norm(K_dense, 1)
    @printf("  K asymmetry (relative): %.6e\n", K_asym)

    # 9. 求解线性系统并检查位移
    u_linear = try
        K_bc \ R_bc
    catch e
        @warn "Linear solve failed" e
        zeros(ndof)
    end
    @printf("\n--- Linear Solution (elastic, D=0) ---\n")
    @printf("  ||u_linear||   = %.6e\n", norm(u_linear))
    @printf("  max(|u_linear|) = %.6e\n", maximum(abs, u_linear))

    # 10. 检查 cohesive 分离位移
    max_sep_n = 0.0
    max_sep_t = 0.0
    for i in 1:n_coh
        geom = cache.cohesive_geom[i]
        dofs = geom.dofs
        L = geom.length
        R_mat = geom.R

        u_e = u_linear[dofs]
        B_global = zeros(2, 8)
        B_global[1, 1] = -0.5; B_global[2, 2] = -0.5
        B_global[1, 3] = -0.5; B_global[2, 4] = -0.5
        B_global[1, 5] = 0.5;  B_global[2, 6] = 0.5
        B_global[1, 7] = 0.5;  B_global[2, 8] = 0.5

        B_local = R_mat * B_global
        δ_local = B_local * u_e
        δ_n = δ_local[1]
        δ_t = δ_local[2]

        if abs(δ_n) > abs(max_sep_n)
            max_sep_n = δ_n
        end
        if abs(δ_t) > abs(max_sep_t)
            max_sep_t = δ_t
        end
    end
    @printf("  max(|δ_n|) at cohesive = %.6e\n", max_sep_n)
    @printf("  max(|δ_t|) at cohesive = %.6e\n", max_sep_t)
    # @printf("  δ_0_n (damage initiation)  = %.6e\n", coh_norm.δ_0_n)      # TODO Chunk 2 Task 2.1
    # if max_sep_n > coh_norm.δ_0_n                                        # TODO Chunk 2 Task 2.1
    #     @printf("  ⚠ WARNING: max(|δ_n|) > δ_0_n! Damage will be triggered!\n")  # TODO Chunk 2 Task 2.1
    # else                                                                 # TODO Chunk 2 Task 2.1
    #     @printf("  ✓ max(|δ_n|) < δ_0_n, elastic regime\n")              # TODO Chunk 2 Task 2.1
    # end                                                                  # TODO Chunk 2 Task 2.1

    # 11. 手动执行一次 Newton 迭代跟踪
    println("\n--- Manual Newton Iteration Tracking (basic solver) ---")
    u = zeros(Float64, ndof)
    states = [begin
        ns = JuBat.DamageState()
        ns.D = s.D; ns.D_visc = s.D_visc
        ns.δ_max_n = s.δ_max_n; ns.δ_max_t = s.δ_max_t
        ns.δ_max_eff = s.δ_max_eff; ns.fractured = s.fractured
        ns.accumulated_damage = s.accumulated_damage
        ns
    end for s in czm_mesh.damage_states]

    for iter in 1:5
        K_tot, f_int, seps, tracts = JuBat.assemble_coupled_system(
            czm_mesh, u, E_eff, ν_eff, coh_norm;
            damage_states=states, K_bulk_cached=cache.K_bulk,
            geom_cache=cache.cohesive_geom, ws=cache.ws)

        R_cur = F_ext + F_thermo_chem - f_int
        for (dof, val) in zip(bc_dofs, bc_vals)
            R_cur[dof] = val - u[dof]
        end
        R_norm = norm(R_cur)

        # 检查损伤状态
        D_vals = [s.D for s in states]
        max_D = maximum(D_vals)
        loading = count(s -> !s.fractured && s.δ_max_eff > 1e-15, states)

        @printf("  iter=%d: ||R||=%.6e, max(D)=%.6f, n_loading=%d\n",
            iter, R_norm, max_D, loading)

        if iter == 1
            R0 = max(R_norm, 1e-10)
        end
        if R_norm < 1e-4 || R_norm / R0 < 1e-4
            println("  → Converged!")
            break
        end

        K_bc_i, R_bc_i = JuBat.apply_bc_czm(K_tot, R_cur; bc_dofs=bc_dofs, bc_vals=bc_vals)
        Δu = try
            K_bc_i \ R_bc_i
        catch e
            println("  → Linear solve failed at iter $iter: $e")
            break
        end

        if any(isnan, Δu) || any(isinf, Δu)
            println("  → NaN/Inf in Δu at iter $iter")
            break
        end

        @printf("    ||Δu||=%.6e, max(|Δu|)=%.6e\n", norm(Δu), maximum(abs, Δu))

        # No line search, just full Newton step
        u .= u .+ Δu
        for (dof, val) in zip(bc_dofs, bc_vals)
            u[dof] = val
        end
    end

    # 12. 检查刚度矩阵对角项
    println("\n--- K_bulk diagnostic ---")
    Kb_dense = Matrix(cache.K_bulk)
    diag_vals = diag(Kb_dense)
    sorted_diag = sort(diag_vals)
    mid_val = sorted_diag[length(sorted_diag) ÷ 2 + 1]
    @printf("  K_bulk diag: min=%.4e, max=%.4e, median=%.4e\n",
        minimum(diag_vals), maximum(diag_vals), mid_val)
    n_zero_diag = count(x -> abs(x) < 1e-15, diag_vals)
    @printf("  Zero-diagonal entries: %d / %d\n", n_zero_diag, length(diag_vals))

    # 13. 深入追踪 basic solver 的线搜索行为
    println("\n--- Line Search Trace (simulating solve_czm_basic_step logic) ---")
    u2 = zeros(Float64, ndof)
    states2 = [JuBat.DamageState() for _ in 1:n_coh]
    ws2 = JuBat.CZMAssemblyWorkspace(ndof, n_coh)
    R_norm_02 = 1.0  # initialize before loop

    for iter in 1:5
        K_tot2, f_int2, seps2, tracts2 = JuBat.assemble_coupled_system(
            czm_mesh, u2, E_eff, ν_eff, coh_norm;
            damage_states=states2, K_bulk_cached=cache.K_bulk,
            geom_cache=cache.cohesive_geom, ws=ws2)

        R2 = F_ext + F_thermo_chem - f_int2
        for (dof, val) in zip(bc_dofs, bc_vals)
            R2[dof] = val - u2[dof]
        end
        R_norm2 = norm(R2)
        @printf("  iter=%d: ||R||=%.6e", iter, R_norm2)

        if iter == 1
            R_norm_02 = max(R_norm2, 1e-10)
        end
        rel_norm2 = R_norm2 / R_norm_02

        if R_norm2 < 1e-4 || rel_norm2 < 1e-4
            @printf("  → converged (||R||=%.2e, rel=%.2e)\n", R_norm2, rel_norm2)
            break
        end
        println()

        K_bc2, R_bc2 = JuBat.apply_bc_czm(K_tot2, R2; bc_dofs=bc_dofs, bc_vals=bc_vals)
        Δu2 = K_bc2 \ R_bc2
        @printf("    ||Δu||=%.6e, max(|Δu|)=%.6e\n", norm(Δu2), maximum(abs, Δu2))

        # Simulate line search with detailed trace
        α = 1.0
        ls_accepted = false
        for halving in 1:8
            u_trial = u2 + α * Δu2
            JuBat.apply_czm_dirichlet!(u_trial, bc_dofs, bc_vals)

            _, f_int_trial, _, _ = JuBat.assemble_coupled_system(
                czm_mesh, u_trial, E_eff, ν_eff, coh_norm;
                damage_states=states2, K_bulk_cached=cache.K_bulk,
                geom_cache=cache.cohesive_geom, ws=ws2)

            R_trial = F_ext + F_thermo_chem - f_int_trial
            for (dof, val) in zip(bc_dofs, bc_vals)
                R_trial[dof] = val - u_trial[dof]
            end
            R_trial_norm = norm(R_trial)
            @printf("    α=%.4f: ||R_trial||=%.6e (target < %.6e)\n", α, R_trial_norm, R_norm2)

            if !isnan(R_trial_norm) && R_trial_norm < R_norm2
                @printf("    → line search accepted at α=%.4f\n", α)
                ls_accepted = true
                break
            end
            α *= 0.5
        end

        if !ls_accepted
            println("    → line search FAILED, would exit solver")
            break
        end

        u2 .= u2 .+ α .* Δu2
        for (dof, val) in zip(bc_dofs, bc_vals)
            u2[dof] = val
        end
    end

    # 14. 直接调用 solve_czm_step 检查实际行为
    println("\n--- Actual solve_czm_step call (method=basic) ---")
    println("  Calling JuBat.solve_czm_step with method='basic'...")
    result_basic, mesh_basic = JuBat.solve_czm_step(
        czm_mesh, F_ext, E_eff, ν_eff, coh_norm, case.param, u0;
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
        max_iter=200, tol=1e-4, n_load_steps=50,
        iter_method="basic", cache=cache)
    @printf("  Result: converged=%s, iterations=%d, ||R||=%.6e\n",
        result_basic.converged, result_basic.iterations, result_basic.residual_norm)

    # 15. 不使用缓存的情况下调用
    println("\n--- Actual solve_czm_step call (method=basic, NO cache) ---")
    result_basic2, mesh_basic2 = JuBat.solve_czm_step(
        czm_mesh, F_ext, E_eff, ν_eff, coh_norm, case.param, u0;
        α_eff=α_eff, β_n=β_n, β_p=β_p,
        dT_elem=dT_elem, Δsoc_n_elem=Δsoc_n_elem, Δsoc_p_elem=Δsoc_p_elem,
        max_iter=200, tol=1e-4, n_load_steps=50,
        iter_method="basic", cache=nothing)
    @printf("  Result: converged=%s, iterations=%d, ||R||=%.6e\n",
        result_basic2.converged, result_basic2.iterations, result_basic2.residual_norm)

    println("\n" * "="^60)
    println("Diagnostics complete.")
    println("="^60)
end

run_diagnostics()
