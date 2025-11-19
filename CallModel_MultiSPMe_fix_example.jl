# CallModel_MultiSPMe 修复示例
# 展示错误用法和正确用法的对比

using SparseArrays

# ============================================================
# 错误示例：试图给 Case 添加新字段（会报错）
# ============================================================

function CallModel_MultiSPMe_WRONG(case::Case, yt::Vector{Float64}, t::Float64; jacobi::String="update")
    """
    ❌ 这个版本会报错！
    """
    variables = SPMe_variables(case, yt, t)
    
    if haskey(case.mesh, "thermal2D")
        mesh_th = case.mesh["thermal2D"]
        
        # ❌ 错误 1：试图给 case 添加不存在的字段
        if !hasproperty(case, :thermal2D_element_area_cache)
            ne = size(mesh_th.element, 1)
            areas = zeros(Float64, ne)
            ngs = length(mesh_th.gs.detJ)
            for g in 1:ngs
                e = mesh_th.gs.ele[g]
                areas[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
            end
            # ❌ 这里会报错：type Case has no field thermal2D_element_area_cache
            case.thermal2D_element_area_cache = areas  # ❌ 错误！
        end
        
        # ❌ 错误 2：从 case 读取不存在的字段
        areas = case.thermal2D_element_area_cache  # ❌ 错误！
    end
    
    # ... 其他代码 ...
    return M, K, F, variables, y_phi
end

# ============================================================
# ✅ 正确示例 1：基本修复
# ============================================================

function CallModel_MultiSPMe_FIXED_BASIC(case::Case, yt::Vector{Float64}, t::Float64; jacobi::String="update")
    """
    ✅ 将缓存存储在 variables 字典中
    """
    variables = SPMe_variables(case, yt, t)
    
    if haskey(case.mesh, "thermal2D")
        mesh_th = case.mesh["thermal2D"]
        
        # ✅ 正确：检查并缓存在 variables 中
        if !haskey(variables, "thermal2D element area")
            ne = size(mesh_th.element, 1)
            areas = zeros(Float64, ne)
            ngs = length(mesh_th.gs.detJ)
            
            @inbounds for g in 1:ngs
                e = mesh_th.gs.ele[g]
                areas[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
            end
            
            # ✅ 存储在 variables 字典中
            variables["thermal2D element area"] = areas
        end
        
        # ✅ 从 variables 读取
        areas = variables["thermal2D element area"]
    end
    
    # ... 其他代码 ...
    return M, K, F, variables, y_phi
end

# ============================================================
# ✅ 正确示例 2：完整的 CallModel_MultiSPMe 实现
# ============================================================

function CallModel_MultiSPMe_FIXED_COMPLETE(case::Case, yt::Vector{Float64}, t::Float64; jacobi::String="update")
    """
    ✅ 参考 Solve.jl 中 CallModel 函数的正确实现
    """
    
    # 1. 获取电化学变量
    variables = SPMe_variables(case, yt, t)
    
    # 2. 计算电化学部分的 M, K, F
    param = case.param
    
    # 负极粒子扩散
    mesh_np = case.mesh["negative particle"]
    csn_gs = variables["negative particle concentration at gauss point"]
    M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, 0.0)
    M_np = M_np .* param.scale.ts_n / case.param_dim.scale.t0
    
    # 正极粒子扩散
    mesh_pp = case.mesh["positive particle"]
    csp_gs = variables["positive particle concentration at gauss point"]
    M_pp, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, 0.0)
    M_pp = M_pp .* param.scale.ts_p / case.param_dim.scale.t0
    
    # 电解质扩散
    mesh_el = case.mesh["electrolyte"]
    M_el, K_el = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, variables)
    M_el = M_el .* param.scale.te / case.param_dim.scale.t0
    
    # 边界条件
    F = SPMe_BC(case, variables)
    
    # 组装矩阵
    M = blockdiag(M_np, M_pp, M_el)
    K = blockdiag(K_np, K_pp, K_el)
    y_phi = Float64[]
    
    # 3. 热模块（如果启用）
    if case.opt.thermal_enabled && case.opt.thermalmodel == "distributed2D"
        if haskey(case.mesh, "thermal2D")
            mesh_th = case.mesh["thermal2D"]
            
            # ✅ 正确：确保 T_nodes 存在（在 variables 中）
            if !haskey(variables, "T_nodes") || isempty(variables["T_nodes"])
                nT = mesh_th.nlen
                variables["T_nodes"] = fill(case.param.cell.T0, nT)
            end
            
            # ✅ 正确：计算并缓存元素面积（在 variables 中）
            if !haskey(variables, "thermal2D element area")
                ne_loc = size(mesh_th.element, 1)
                A_loc = zeros(Float64, ne_loc)
                ngs_loc = length(mesh_th.gs.detJ)
                
                @inbounds for g in 1:ngs_loc
                    e = mesh_th.gs.ele[g]
                    A_loc[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
                end
                
                # ✅ 存储在 variables 中
                variables["thermal2D element area"] = A_loc
            end
            
            # ✅ 从 variables 读取缓存的面积
            areas = variables["thermal2D element area"]
            ne_loc = length(areas)
            
            # ✅ 计算元素平均温度（使用 variables 中的 T_nodes）
            T_nodes_loc = variables["T_nodes"]
            Te_prev = zeros(Float64, ne_loc)
            
            @inbounds for e in 1:ne_loc
                nds = mesh_th.element[e, :]
                Te_prev[e] = sum(T_nodes_loc[nds]) / length(nds)
            end
            
            # ✅ 获取总电流（从 variables 中）
            I_total = 0.0
            if haskey(variables, "cell current")
                I_total = Float64(variables["cell current"])
            end
            
            # ✅ 求解分支电流（结果存储在 variables 中）
            try
                variables, _Ie, _Vc = solve_branch_currents_newton(
                    case, variables, yt, t, I_total, areas, Te_prev, nothing
                )
            catch err
                error("solve_branch_currents_newton failed: $(err)")
            end
            
            # ✅ 计算热源（结果存储在 variables 中）
            try
                variables = heatQ_Source(case, variables, t, yt)
            catch err
                @warn "heatQ_Source failed, using zero heat" err
                variables["heat_source_fields"] = zeros(Float64, ne_loc)
            end
            
            # ✅ 装配热学矩阵
            MT, KT, FT = ThermalDistributed2D(case, variables)
            
            # 时间尺度匹配
            t_ratio = case.param_dim.scale.t0 / case.param_dim.scale.t_th
            MT = MT .* t_ratio
            
            # 边界条件
            ThermalDistributed2D_BC(KT, FT, case, t)
            
            # 拼接热学部分
            M = blockdiag(M, sparse(MT))
            K = blockdiag(K, sparse(KT))
            F = [F; FT]
        end
    end
    
    return M, K, F, variables, y_phi
end

# ============================================================
# 对比总结
# ============================================================

"""
关键区别：

❌ 错误做法：
   case.xxx_cache = value              # 试图给 Case 添加字段
   if hasproperty(case, :xxx_cache)    # 检查 case 的字段
   value = case.xxx_cache              # 从 case 读取

✅ 正确做法：
   variables["xxx"] = value            # 存储在 variables 字典中
   if haskey(variables, "xxx")         # 检查 variables 的键
   value = variables["xxx"]            # 从 variables 读取

记住：
- Case 结构体 → 只有5个固定字段，不能添加新字段
- variables 字典 → 可以动态添加任意键值对
- 所有计算结果和缓存 → 都应该存储在 variables 中
"""

# ============================================================
# 常见的需要缓存的数据（都应该存储在 variables 中）
# ============================================================

"""
所有这些都应该存储在 variables 字典中：

✅ variables["thermal2D element area"]        # 元素面积
✅ variables["T_nodes"]                        # 节点温度
✅ variables["T_prev"]                         # 前一步温度
✅ variables["heat_source_fields"]             # 热源场
✅ variables["thermal2D element current"]      # 元素电流
✅ variables["thermal2D element current A"]    # 元素电流（安培）
✅ variables["thermal2D common voltage"]       # 公共电压
✅ variables["thermal2D layer_weights"]        # 层权重
✅ variables["thermal2D Vsolve status"]        # 求解状态
✅ variables["thermal2D Vsolve iters"]         # 迭代次数
✅ variables["thermal2D Vsolve converged"]     # 是否收敛

❌ case.thermal2D_element_area_cache           # 错误！
❌ case.T_nodes_cache                          # 错误！
❌ case.heat_source_cache                      # 错误！
"""

println("✅ 示例代码加载完成")
println("请参考 CallModel_MultiSPMe_FIXED_BASIC 或 CallModel_MultiSPMe_FIXED_COMPLETE")
println("将所有 case.xxx_cache 改为 variables[\"xxx\"]")
