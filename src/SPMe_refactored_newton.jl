# SPMe.jl - solve_branch_currents_newton 函数精简版
# 
# 精简策略：
# 1. 提取调试代码到独立函数
# 2. 提取电化学系数计算
# 3. 提取牛顿迭代主循环
# 4. 统一边界检查逻辑
# 5. 简化条件分支

# ========================================================================
# 辅助函数：调试输出
# ========================================================================

"""检查并报告 NaN/Inf 值（调试用）"""
function _debug_check_prefactors(prefactor_n, prefactor_p, csn_av, csp_av, 
                                  u_n_ref_val, u_p_ref_val, du_n_dT_val, du_p_dT_val, 
                                  c_sigma, cn_surf, cp_surf, ce_n_gs, ce_p_gs)
    has_nan = !isfinite(prefactor_n) || !isfinite(prefactor_p) || 
              !isfinite(csn_av) || !isfinite(csp_av) || 
              !isfinite(u_n_ref_val) || !isfinite(u_p_ref_val) || 
              !isfinite(c_sigma)
    
    if !has_nan
        return false
    end
    
    println("\n" * "="^80)
    println("❌ [DEBUG] 预计算值包含 NaN/Inf")
    println("="^80)
    println("📊 预因子: prefactor_n=$prefactor_n, prefactor_p=$prefactor_p")
    println("📊 平均浓度: csn_av=$csn_av, csp_av=$csp_av")
    println("📊 开路电位: u_n=$u_n_ref_val, u_p=$u_p_ref_val")
    println("📊 温度系数: du_n_dT=$du_n_dT_val, du_p_dT=$du_p_dT_val")
    println("📊 电导: c_sigma=$c_sigma")
    println("\n💡 输入浓度前5个值:")
    println("  cn_surf: $(cn_surf[1:min(5,end)])")
    println("  cp_surf: $(cp_surf[1:min(5,end)])")
    println("  ce_n_gs: $(ce_n_gs[1:min(5,end)])")
    println("  ce_p_gs: $(ce_p_gs[1:min(5,end)])")
    println("="^80 * "\n")
    
    return true
end

"""检查单元系数是否有效"""
function _debug_check_coefficients(e, has_nan_prefactor, C1, C2, alpha_p, alpha_n, C5, 
                                    T_e, j0_n, j0_p, u_n_val_T, u_p_val_T)
    if has_nan_prefactor || e > 1
        return  # 只打印第一个异常单元，且预计算值正常时
    end
    
    has_nan_coeff = !isfinite(C1) || !isfinite(C2) || 
                    !isfinite(alpha_p) || !isfinite(alpha_n) || !isfinite(C5)
    
    if !has_nan_coeff
        return
    end
    
    println("\n" * "="^80)
    println("❌ [DEBUG] 单元 $e 系数异常")
    println("="^80)
    println("📊 T_e=$T_e, j0_n=$j0_n, j0_p=$j0_p")
    println("📊 u_n(T)=$u_n_val_T, u_p(T)=$u_p_val_T")
    println("📊 系数: C1=$C1, C2=$C2")
    println("   alpha_p=$alpha_p, alpha_n=$alpha_n, C5=$C5")
    println("="^80 * "\n")
end

"""检查初始电压是否有效"""
function _debug_check_initial_voltage(has_nan_prefactor, V, V_branches, I_e, 
                                      coeffs, I_total, ne)
    if has_nan_prefactor || isfinite(V)
        return
    end
    
    println("\n" * "="^80)
    println("❌ [DEBUG] 初始电压异常: V=$V")
    println("="^80)
    println("  I_total=$I_total, ne=$ne")
    println("  前3个异常单元:")
    
    count = 0
    for e in 1:ne
        if !isfinite(V_branches[e]) && count < 3
            count += 1
            println("  单元 $e: V=$(V_branches[e]), I=$(I_e[e])")
            println("    C1=$(coeffs[e].C1), C2=$(coeffs[e].C2)")
            println("    α_p=$(coeffs[e].alpha_p), α_n=$(coeffs[e].alpha_n)")
        end
    end
    println("="^80 * "\n")
end

# ========================================================================
# 辅助函数：电化学计算
# ========================================================================

"""标量化：将数组转为标量（取第一个元素）"""
_scalarize(x) = isa(x, Number) ? Float64(x) : Float64(x[1])

"""计算电化学预因子"""
function _compute_electrochemical_prefactors(variables, param, mesh_ne, mesh_pe)
    cn_surf = variables["negative particle surface lithium concentration"]
    cp_surf = variables["positive particle surface lithium concentration"]
    ce_n_gs = variables["electrolyte lithium concentration at negative electrode Gauss point"]
    ce_p_gs = variables["electrolyte lithium concentration at positive electrode Gauss point"]
    
    # 预因子计算
    prefactor_n = IntV(abs.(cn_surf .* (1.0 .- cn_surf) .* ce_n_gs) .^ 0.5, mesh_ne) / param.NE.thickness
    prefactor_p = IntV(abs.(cp_surf .* (1.0 .- cp_surf) .* ce_p_gs) .^ 0.5, mesh_pe) / param.PE.thickness
    csn_av = IntV(ce_n_gs, mesh_ne) / param.NE.thickness
    csp_av = IntV(ce_p_gs, mesh_pe) / param.PE.thickness
    
    # 开路电位
    u_n_ref = param.NE.U(cn_surf)
    u_p_ref = param.PE.U(cp_surf)
    du_n_dT = param.NE.dUdT(cn_surf)
    du_p_dT = param.PE.dUdT(cp_surf)
    
    # 标量化
    u_n_ref_val = _scalarize(u_n_ref)
    u_p_ref_val = _scalarize(u_p_ref)
    du_n_dT_val = _scalarize(du_n_dT)
    du_p_dT_val = _scalarize(du_p_dT)
    
    # 固相电导
    c_sigma = (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig) / 3.0
    
    return (prefactor_n=prefactor_n, prefactor_p=prefactor_p,
            csn_av=csn_av, csp_av=csp_av,
            u_n_ref_val=u_n_ref_val, u_p_ref_val=u_p_ref_val,
            du_n_dT_val=du_n_dT_val, du_p_dT_val=du_p_dT_val,
            c_sigma=c_sigma,
            cn_surf=cn_surf, cp_surf=cp_surf, ce_n_gs=ce_n_gs, ce_p_gs=ce_p_gs)
end

"""计算单个单元的电化学系数"""
function _compute_element_coefficients(e, T_e, param, prefactors, T_ref, debug_mode=false)
    # 交换电流密度
    arr_n = Arrhenius(param.NE.Eac_k, T_e)
    arr_p = Arrhenius(param.PE.Eac_k, T_e)
    j0_n = max(param.NE.k * arr_n * prefactors.prefactor_n, 1e-16)
    j0_p = max(param.PE.k * arr_p * prefactors.prefactor_p, 1e-16)
    
    # 电解液电导率
    kappa_ne = max(param.EL.kappa(param.EL.ce0, T_e) * param.NE.eps^param.NE.brugg, 1e-16)
    kappa_pe = max(param.EL.kappa(param.EL.ce0, T_e) * param.PE.eps^param.PE.brugg, 1e-16)
    kappa_sp = max(param.EL.kappa(param.EL.ce0, T_e) * param.SP.eps^param.SP.brugg, 1e-16)
    R_EL = param.NE.thickness / (3.0 * kappa_ne) + 
           param.SP.thickness / kappa_sp + 
           param.PE.thickness / (3.0 * kappa_pe)
    
    # 开路电位（温度相关）
    u_n = prefactors.u_n_ref_val + (T_e - T_ref) * prefactors.du_n_dT_val
    u_p = prefactors.u_p_ref_val + (T_e - T_ref) * prefactors.du_p_dT_val
    
    # 计算系数
    C1 = (u_p - u_n) + 2.0 * T_e * (1.0 - param.EL.tplus) * 
         (prefactors.csp_av - prefactors.csn_av) / param.EL.ce0
    C2 = 2.0 * T_e
    alpha_p = -1.0 / (2.0 * j0_p * param.PE.as * param.PE.thickness)
    alpha_n = 1.0 / (2.0 * j0_n * param.NE.as * param.NE.thickness)
    C5 = R_EL + prefactors.c_sigma
    
    # 调试检查（仅第一个异常单元）
    if debug_mode
        _debug_check_coefficients(e, false, C1, C2, alpha_p, alpha_n, C5, 
                                  T_e, j0_n, j0_p, u_n, u_p)
    end
    
    return (C1=C1, C2=C2, alpha_p=alpha_p, alpha_n=alpha_n, C5=C5)
end

"""批量计算所有单元的系数"""
function _compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref, debug_mode=false)
    coeffs = Vector{NamedTuple{(:C1,:C2,:alpha_p,:alpha_n,:C5)}}(undef, ne)
    for e in 1:ne
        coeffs[e] = _compute_element_coefficients(e, Te_prev[e], param, prefactors, T_ref, debug_mode)
    end
    return coeffs
end

# ========================================================================
# 辅助函数：分支电压模型
# ========================================================================

"""计算分支电压 V = C1 + C2*(asinh(α_p*I) - asinh(α_n*I)) - C5*I"""
function _branch_voltage(coeff, I::Float64)
    apI = coeff.alpha_p * I
    anI = coeff.alpha_n * I
    return coeff.C1 + coeff.C2 * (asinh(apI) - asinh(anI)) - coeff.C5 * I
end

"""计算分支电压对电流的导数 dV/dI"""
function _branch_dVdI(coeff, I::Float64)
    apI = coeff.alpha_p * I
    anI = coeff.alpha_n * I
    denom_p = sqrt(1.0 + apI * apI)
    denom_n = sqrt(1.0 + anI * anI)
    return coeff.C2 * (coeff.alpha_p / denom_p - coeff.alpha_n / denom_n) - coeff.C5
end

# ========================================================================
# 辅助函数：初始化和边界检查
# ========================================================================

"""初始化单元电流猜测"""
function _initialize_currents(ne, w, I_total, x_prev)
    I_e = x_prev !== nothing && length(x_prev) == ne ? copy(x_prev) : (w .* I_total)
    
    # 归一化到满足总电流约束
    sx = sum(w .* I_e)
    if sx != 0.0
        I_e .*= (I_total / sx)
    else
        I_e .= w .* I_total
    end
    
    return I_e
end

"""检查电压边界"""
function _check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e, context="")
    if V_MIN <= V <= V_MAX
        return true
    end
    
    V_phys = V * phi_scale
    V_MIN_phys = V_MIN * phi_scale
    V_MAX_phys = V_MAX * phi_scale
    
    error_msg = "thermal2D common voltage out of bounds$context: " *
                "V(nd)=$V, V(V)=$V_phys, " *
                "allowed [$V_MIN, $V_MAX] nd -> [$V_MIN_phys, $V_MAX_phys] V; " *
                "I_total_nd=$I_total, sum(w.*I_e)=$(sum(w .* I_e))"
    
    throw(ErrorException(error_msg))
end

# ========================================================================
# 辅助函数：牛顿迭代
# ========================================================================

"""牛顿迭代主循环"""
function _newton_iteration!(I_e, V, ne, w, I_total, coeffs; 
                           tol_V=1e-8, tol_I=1e-10, max_iters=25)
    converged = false
    last_iter = 0
    F = zeros(Float64, ne)
    dFdI = similar(F)
    I_trial = similar(I_e)
    
    for iter in 1:max_iters
        last_iter = iter
        
        # 计算残差和雅可比
        for e in 1:ne
            V_e = _branch_voltage(coeffs[e], I_e[e])
            F[e] = V_e - V
            dFdI[e] = _branch_dVdI(coeffs[e], I_e[e])
            
            # 防止奇异雅可比
            if abs(dFdI[e]) < 1e-12
                dFdI[e] = sign(dFdI[e]) != 0.0 ? sign(dFdI[e]) * 1e-12 : -coeffs[e].C5
            end
        end
        
        # 检查收敛
        res_V = maximum(abs.(F))
        res_I = sum(w .* I_e) - I_total
        
        if res_V <= tol_V && abs(res_I) <= tol_I
            converged = true
            break
        end
        
        # 牛顿步
        denom = sum(w ./ dFdI)
        abs(denom) < 1e-12 && break
        
        num = -res_I + sum(w .* F ./ dFdI)
        ΔV = num / denom
        ΔI = ((-F) .+ ΔV) ./ dFdI
        
        # 线搜索
        λ, V_trial = _line_search(I_e, V, ΔI, ΔV, I_trial, ne)
        λ == 0.0 && break
        
        # 更新
        I_e .= I_trial
        V = V_trial
    end
    
    return V, converged, last_iter
end

"""线搜索：确保更新后的值有效"""
function _line_search(I_e, V, ΔI, ΔV, I_trial, ne; max_attempts=12)
    λ = 1.0
    
    for attempt in 1:max_attempts
        # 试探电压
        V_trial = V + λ * ΔV
        !isfinite(V_trial) && (λ *= 0.5; continue)
        
        # 试探电流
        ok = true
        @inbounds for e in 1:ne
            val = I_e[e] + λ * ΔI[e]
            if !isfinite(val) || abs(val) > 1e12
                ok = false
                break
            end
            I_trial[e] = val
        end
        
        ok && return λ, V_trial
        λ *= 0.5
    end
    
    return 0.0, V  # 失败
end

# ========================================================================
# 主函数：精简版 solve_branch_currents_newton
# ========================================================================

"""
    solve_branch_currents_newton(case, variables, yt, t, I_total, areas, Te_prev, x_prev)

非线性分流求解器（精简版）

使用牛顿法求解电流分配问题：
- 未知量：各单元电流 I_e 和公共端电压 V
- 约束：V_e(I_e, T_e) = V (每个单元) 和 Σ(I_e*A_e) = I_total

# 参数
- `case`: 案例对象
- `variables`: 变量字典
- `yt`: 状态向量
- `t`: 时间
- `I_total`: 总电流（无量纲）
- `areas`: 各单元面积
- `Te_prev`: 各单元温度
- `x_prev`: 上一步的电流分配（可选）

# 返回
- `variables`: 更新后的变量字典
- `I_e`: 各单元电流
- `V`: 公共端电压
"""
function solve_branch_currents_newton(case::Case, 
                                     variables::Dict{String,Union{Array{Float64},Float64}}, 
                                     yt::Array{Float64}, t::Float64, I_total::Float64, 
                                     areas::Vector{Float64}, Te_prev::Vector{Float64}, 
                                     x_prev::Union{Nothing,Vector{Float64}}=nothing)
    
    # 1. 快速退出：无热网格时按面积分配
    if !haskey(case.mesh, "thermal2D")
        return _fallback_solution(variables, areas, I_total, case.param.scale.I_typ)
    end
    
    # 2. 初始化
    ne = length(areas)
    A_global = sum(areas)
    w = areas ./ A_global  # 面积权重
    phi_scale = case.param.scale.phi
    V_MIN = hasproperty(case.param_dim.cell, :v_l) ? case.param_dim.cell.v_l / phi_scale : -Inf
    V_MAX = hasproperty(case.param_dim.cell, :v_h) ? case.param_dim.cell.v_h / phi_scale : Inf
    
    # 3. 计算电化学预因子
    param = case.param
    mesh_ne = case.mesh["negative electrode"]
    mesh_pe = case.mesh["positive electrode"]
    prefactors = _compute_electrochemical_prefactors(variables, param, mesh_ne, mesh_pe)
    
    # 调试：检查预因子
    debug_mode = hasproperty(case.opt, :debug_coupling) && case.opt.debug_coupling
    has_nan_prefactor = _debug_check_prefactors(
        prefactors.prefactor_n, prefactors.prefactor_p, 
        prefactors.csn_av, prefactors.csp_av,
        prefactors.u_n_ref_val, prefactors.u_p_ref_val,
        prefactors.du_n_dT_val, prefactors.du_p_dT_val,
        prefactors.c_sigma,
        prefactors.cn_surf, prefactors.cp_surf, 
        prefactors.ce_n_gs, prefactors.ce_p_gs
    )
    
    # 4. 计算各单元系数
    T_ref = case.param.cell.T0
    coeffs = _compute_all_coefficients(ne, Te_prev, param, prefactors, T_ref, debug_mode)
    
    # 5. 零电流特殊处理
    if abs(I_total) <= 1e-14
        return _zero_current_solution(variables, coeffs, ne, V_MIN, V_MAX, phi_scale, w)
    end
    
    # 6. 初始化电流猜测
    I_e = _initialize_currents(ne, w, I_total, x_prev)
    
    # 计算初始电压
    V_branches = [_branch_voltage(coeffs[e], I_e[e]) for e in 1:ne]
    V = sum(V_branches) / ne
    
    # 调试：检查初始电压
    _debug_check_initial_voltage(has_nan_prefactor, V, V_branches, I_e, coeffs, I_total, ne)
    
    # 7. 牛顿迭代求解
    V, converged, last_iter = _newton_iteration!(I_e, V, ne, w, I_total, coeffs)
    
    # 8. 未收敛时回退到面积分配
    if !converged
        I_e .= w .* I_total
        V = sum(coeffs[e].C1 for e in 1:ne) / ne
    end
    
    # 9. 归一化确保总电流约束
    sx = sum(w .* I_e)
    if sx != 0.0
        I_e .*= (I_total / sx)
    end
    
    # 10. 边界检查
    _check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, I_total, w, I_e)
    
    # 11. 写入结果
    variables["thermal2D element current"] = I_e
    variables["thermal2D element current A"] = case.param.scale.I_typ .* I_e
    variables["thermal2D common voltage"] = V
    variables["thermal2D Vsolve status"] = converged ? 3.0 : 3.5
    variables["thermal2D Vsolve iters"] = Float64(last_iter)
    variables["thermal2D Vsolve converged"] = converged ? 1.0 : 0.0
    
    return variables, I_e, V
end

# ========================================================================
# 特殊情况处理
# ========================================================================

"""无热网格时的回退解法：按面积均匀分配"""
function _fallback_solution(variables, areas, I_total, I_typ)
    A_tot = sum(areas)
    w = areas ./ A_tot
    I_e = w .* I_total
    
    variables["thermal2D element current"] = I_e
    variables["thermal2D element current A"] = I_typ .* I_e
    variables["thermal2D common voltage"] = 0.0
    variables["thermal2D Vsolve status"] = 1.0
    variables["thermal2D Vsolve iters"] = 0.0
    variables["thermal2D Vsolve converged"] = 0.0
    
    return variables, I_e, 0.0
end

"""零电流特殊解法"""
function _zero_current_solution(variables, coeffs, ne, V_MIN, V_MAX, phi_scale, w)
    I_e = zeros(Float64, ne)
    V = sum(coeffs[e].C1 for e in 1:ne) / ne
    
    # 边界检查
    _check_voltage_bounds(V, V_MIN, V_MAX, phi_scale, 0.0, w, I_e, " at zero-current")
    
    variables["thermal2D element current"] = I_e
    variables["thermal2D common voltage"] = V
    variables["thermal2D Vsolve status"] = 0.5
    variables["thermal2D Vsolve iters"] = 0.0
    variables["thermal2D Vsolve converged"] = 1.0
    
    return variables, I_e, V
end
