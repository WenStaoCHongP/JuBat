# test3.jl 错误修复示例
# 展示常见的错误模式和正确的修复方法

# ============================================================
# 示例 1: 检查 case 是否有某个字段
# ============================================================

# ❌ 错误写法
# if haskey(case, :mesh)
#     println("case has mesh")
# end

# ✅ 正确写法
if hasproperty(case, :mesh)
    println("case has mesh")
end

# ============================================================
# 示例 2: 检查 case.opt 的字段
# ============================================================

# ❌ 错误写法
# if haskey(case.opt, :thermal_enabled)
#     if case.opt.thermal_enabled
#         println("Thermal enabled")
#     end
# end

# ✅ 正确写法
if hasproperty(case.opt, :thermal_enabled)
    if case.opt.thermal_enabled
        println("Thermal enabled")
    end
end

# ============================================================
# 示例 3: 检查多个字段的组合条件
# ============================================================

# ❌ 错误写法
# thermal_distributed = haskey(case.opt, :thermal_enabled) && case.opt.thermal_enabled &&
#                      haskey(case.opt, :thermalmodel) && case.opt.thermalmodel == "distributed2D"

# ✅ 正确写法（参考 SPMe.jl 第143-144行）
thermal_distributed = hasproperty(case.opt, :thermal_enabled) && case.opt.thermal_enabled &&
                     hasproperty(case.opt, :thermalmodel) && case.opt.thermalmodel == "distributed2D"

# ============================================================
# 示例 4: 检查字典的键（这些是正确的，保持不变）
# ============================================================

# ✅ 正确 - case.mesh 是 Dict，应该用 haskey
if haskey(case.mesh, "thermal2D")
    mesh_th = case.mesh["thermal2D"]
    println("Thermal mesh exists")
end

# ✅ 正确 - case.index 是 Dict，应该用 haskey
if haskey(case.index, "temperature")
    temp_index = case.index["temperature"]
end

# ✅ 正确 - variables 是 Dict，应该用 haskey
if haskey(variables, "T_nodes")
    T_nodes = variables["T_nodes"]
end

# ============================================================
# 示例 5: CallModel_MultiSPMe 中可能的错误模式
# ============================================================

function CallModel_MultiSPMe_FIXED(case::Case, yt::Vector{Float64}, t::Float64; jacobi::String="update")
    """
    这是一个假设的修复示例
    根据错误信息，原函数可能在第277行附近使用了 haskey(case, ...)
    """
    
    # ❌ 可能的错误代码（需要修复）：
    # if haskey(case, :mesh)
    #     ...
    # end
    
    # ✅ 修复后的代码：
    if hasproperty(case, :mesh)
        # 处理逻辑
    end
    
    # ❌ 可能的错误代码：
    # if haskey(case, :opt) && haskey(case.opt, :thermal_enabled)
    #     ...
    # end
    
    # ✅ 修复后的代码：
    if hasproperty(case, :opt) && hasproperty(case.opt, :thermal_enabled)
        if case.opt.thermal_enabled
            # 热模块处理
        end
    end
    
    # ✅ 正确的字典检查（不需要改）：
    if hasproperty(case, :mesh) && haskey(case.mesh, "thermal2D")
        mesh_th = case.mesh["thermal2D"]
    end
    
    # 返回示例值（实际函数返回值会不同）
    return nothing
end

# ============================================================
# 示例 6: 安全的字段检查辅助函数
# ============================================================

"""
安全地检查并获取 case 的字段值
"""
function safe_get_field(obj, field::Symbol, default=nothing)
    if hasproperty(obj, field)
        return getproperty(obj, field)
    else
        return default
    end
end

# 使用示例：
# thermal_enabled = safe_get_field(case.opt, :thermal_enabled, false)
# thermalmodel = safe_get_field(case.opt, :thermalmodel, "none")

# ============================================================
# 示例 7: 类型检查与字段检查结合
# ============================================================

function check_thermal_setup(case::Case)
    """检查热模块配置的完整示例"""
    
    # 检查 opt 字段是否存在
    if !hasproperty(case, :opt)
        @warn "case 没有 opt 字段"
        return false
    end
    
    # 检查 thermal_enabled 字段
    if !hasproperty(case.opt, :thermal_enabled)
        @warn "case.opt 没有 thermal_enabled 字段"
        return false
    end
    
    if !case.opt.thermal_enabled
        println("热模块未启用")
        return false
    end
    
    # 检查 thermalmodel 字段
    if !hasproperty(case.opt, :thermalmodel)
        @warn "case.opt 没有 thermalmodel 字段"
        return false
    end
    
    # 检查 mesh 字典
    if !hasproperty(case, :mesh)
        @warn "case 没有 mesh 字段"
        return false
    end
    
    # 检查 thermal2D mesh（字典操作）
    if case.opt.thermalmodel == "distributed2D" && !haskey(case.mesh, "thermal2D")
        @warn "thermalmodel 设置为 distributed2D 但缺少 thermal2D mesh"
        return false
    end
    
    println("✅ 热模块配置检查通过")
    return true
end

# ============================================================
# 总结
# ============================================================

"""
记忆规则：
1. 结构体（struct）→ hasproperty(obj, :field) 或 hasfield(typeof(obj), :field)
2. 字典（Dict）→ haskey(dict, key)
3. Case 是 struct → hasproperty
4. case.mesh 是 Dict → haskey  
5. case.index 是 Dict → haskey
6. case.opt 是 Option struct → hasproperty
7. case.param 是 Params struct → hasproperty
8. variables 通常是 Dict → haskey
"""

println("✅ 示例代码加载完成")
println("请参考上述示例修复您的 test3.jl 和 CallModel_MultiSPMe 函数")
