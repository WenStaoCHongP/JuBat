"""
诊断工具：检查 Scale 结构中 T_ref 的初始化

该工具用于追踪 T_ref 在整个初始化流程中的变化：
1. 参数文件中的初始值
2. ChooseCell 后的值
3. SetCase 后的值（归一化参数）
"""

using LinearAlgebra, SparseArrays, Statistics, Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

function check_T_ref_initialization()
    println("="^80)
    println("T_ref 初始化诊断工具")
    println("="^80)
    
    # ========================================================================
    # 第 1 步：直接检查参数文件中的 Scale 默认值
    # ========================================================================
    println("\n[步骤 1] 检查 Scale 结构的默认值:")
    default_scale = JuBat.Scale()
    println("  T_ref (默认值) = $(default_scale.T_ref) K")
    println("  R (默认值) = $(default_scale.R) J/(mol·K)")
    println("  F (默认值) = $(default_scale.F) C/mol")
    
    if default_scale.T_ref <= 0 || isnan(default_scale.T_ref) || isinf(default_scale.T_ref)
        @warn "⚠️  默认 T_ref 值异常！" value=default_scale.T_ref
    else
        println("  ✓ 默认 T_ref 值正常")
    end
    
    # ========================================================================
    # 第 2 步：检查 ChooseCell 后的 param_dim.scale.T_ref
    # ========================================================================
    println("\n[步骤 2] 检查 ChooseCell('Jellyroll') 后的 param_dim.scale:")
    param_dim = JuBat.ChooseCell("Jellyroll")
    
    println("  T_ref = $(param_dim.scale.T_ref) K")
    println("  T_amb = $(param_dim.cell.T_amb) K")
    println("  T0 = $(param_dim.cell.T0) K")
    println("  L_th = $(param_dim.scale.L_th) m")
    println("  k_th = $(param_dim.scale.k_th) W/(m·K)")
    println("  q_th = $(param_dim.scale.q_th) W/m³")
    println("  phi = $(param_dim.scale.phi) V")
    
    if param_dim.scale.T_ref <= 0 || isnan(param_dim.scale.T_ref) || isinf(param_dim.scale.T_ref)
        @error "❌ ChooseCell 后 T_ref 值异常！" value=param_dim.scale.T_ref
        println("\n可能的原因分析:")
        println("  1. Jellyroll.jl 中 Scale() 初始化被覆盖")
        println("  2. ChooseCell 中的衍生参数计算修改了 T_ref")
        return
    else
        println("  ✓ ChooseCell 后 T_ref 值正常")
    end
    
    # 检查相关衍生参数
    println("\n  检查相关衍生参数的合理性:")
    if param_dim.cell.Rout <= 0
        @warn "  ⚠️  外半径 Rout 异常" value=param_dim.cell.Rout
    else
        println("    Rout = $(param_dim.cell.Rout) m ✓")
    end
    
    if param_dim.scale.L_th <= 0
        @warn "  ⚠️  热特征长度 L_th 异常" value=param_dim.scale.L_th
    else
        println("    L_th = $(param_dim.scale.L_th) m ✓")
    end
    
    if param_dim.scale.k_th <= 0
        @warn "  ⚠️  参考热导率 k_th 异常" value=param_dim.scale.k_th
    else
        println("    k_th = $(param_dim.scale.k_th) W/(m·K) ✓")
    end
    
    # ========================================================================
    # 第 3 步：检查 SetCase 后的归一化参数
    # ========================================================================
    println("\n[步骤 3] 检查 SetCase 后的 param.scale (归一化参数):")
    
    opt = JuBat.Option()
    opt.model = "SPMe"
    opt.Nn = 10
    opt.Ns = 5
    opt.Np = 10
    opt.Nrn = 10
    opt.Nrp = 10
    opt.Current = x -> 5.0
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.time = [0.0, 100.0]
    
    case = JuBat.SetCase(param_dim, opt)
    
    println("  case.param_dim.scale.T_ref = $(case.param_dim.scale.T_ref) K (原始)")
    println("  case.param.scale.T_ref = $(case.param.scale.T_ref) K (归一化后)")
    
    if case.param_dim.scale.T_ref != param_dim.scale.T_ref
        @warn "⚠️  SetCase 后 param_dim.scale.T_ref 值发生了变化！"
    end
    
    if case.param.scale.T_ref <= 0 || isnan(case.param.scale.T_ref) || isinf(case.param.scale.T_ref)
        @error "❌ 归一化参数中 T_ref 值异常！" value=case.param.scale.T_ref
    else
        println("  ✓ 归一化参数中 T_ref 值正常")
    end
    
    # ========================================================================
    # 第 4 步：模拟极耳边界条件计算
    # ========================================================================
    println("\n[步骤 4] 模拟极耳边界条件中的 T_ref 使用:")
    println("  这是 _apply_tab_bc! 函数中使用 T_ref 的方式：")
    
    scale = case.param_dim.scale
    T_amb = case.param_dim.cell.T_amb
    t = 0.0
    rate_Ks = 0.1
    
    println("  scale.T_ref = $(scale.T_ref) K")
    println("  T_amb = $(T_amb) K")
    println("  rate_Ks = $(rate_Ks) K/s")
    println("  t = $(t) s")
    
    println("\n  计算过程:")
    println("    T_amb_nd = T_amb / scale.T_ref")
    T_amb_nd = T_amb / scale.T_ref
    println("             = $(T_amb) / $(scale.T_ref)")
    println("             = $(T_amb_nd)")
    
    if isnan(T_amb_nd) || isinf(T_amb_nd)
        @error "❌ T_amb_nd 计算结果异常！"
        println("\n根因分析:")
        if scale.T_ref == 0
            println("  scale.T_ref 为 0，导致除零错误")
        elseif isnan(scale.T_ref)
            println("  scale.T_ref 为 NaN")
        elseif isnan(T_amb)
            println("  T_amb 为 NaN")
        end
    else
        println("    ✓ T_amb_nd 计算正常")
    end
    
    println("\n    T_tab_nd = T_amb_nd + (rate_Ks * t) / scale.T_ref")
    T_tab_nd = T_amb_nd + (rate_Ks * t) / scale.T_ref
    println("             = $(T_amb_nd) + ($(rate_Ks) * $(t)) / $(scale.T_ref)")
    println("             = $(T_tab_nd)")
    
    if isnan(T_tab_nd) || isinf(T_tab_nd)
        @error "❌ T_tab_nd 计算结果异常！"
    else
        println("    ✓ T_tab_nd 计算正常")
    end
    
    # ========================================================================
    # 第 5 步：检查极耳数组
    # ========================================================================
    println("\n[步骤 5] 检查极耳角度数组:")
    println("  默认值（Jellyroll.jl）：")
    println("    theta_pos = $(param_dim.tab.theta_pos)")
    println("    theta_neg = $(param_dim.tab.theta_neg)")
    
    if isempty(param_dim.tab.theta_pos) && isempty(param_dim.tab.theta_neg)
        println("  ✓ 极耳数组为空（默认状态）")
        println("  ⚠️  注意：当添加极耳角度值后，才会触发极耳边界条件")
        println("\n  如果您添加了极耳角度（例如：param_dim.tab.theta_pos = [0.0]）")
        println("  请重新运行此诊断工具，并传入修改后的 param_dim")
    else
        println("  ℹ️  极耳数组已设置，将应用极耳边界条件")
    end
    
    # ========================================================================
    # 总结
    # ========================================================================
    println("\n" * "="^80)
    println("诊断总结")
    println("="^80)
    println("T_ref 初始化流程检查：")
    println("  1. Scale() 默认值: $(default_scale.T_ref) K")
    println("  2. ChooseCell 后: $(param_dim.scale.T_ref) K")
    println("  3. SetCase 后 (param_dim): $(case.param_dim.scale.T_ref) K")
    println("  4. SetCase 后 (param 归一化): $(case.param.scale.T_ref) K")
    println("  5. 极耳边界条件计算: T_amb_nd=$(T_amb_nd), T_tab_nd=$(T_tab_nd)")
    
    if all([
        default_scale.T_ref > 0,
        param_dim.scale.T_ref > 0,
        case.param_dim.scale.T_ref > 0,
        case.param.scale.T_ref > 0,
        !isnan(T_amb_nd),
        !isnan(T_tab_nd)
    ])
        println("\n✅ 所有检查通过！T_ref 初始化正常。")
        println("\n如果您仍然遇到 NaN 错误，可能的原因：")
        println("  1. 在添加极耳角度后，某个后续计算步骤修改了 param_dim")
        println("  2. 热网格创建过程中出现了问题")
        println("  3. 初始化函数中存在其他数值问题")
    else
        println("\n⚠️  发现异常值，请检查上述步骤中标记为 ❌ 或 ⚠️  的项")
    end
    
    println("="^80)
end

# 运行诊断
check_T_ref_initialization()
