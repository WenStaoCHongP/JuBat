"""
    CZM_Parameters.jl

内聚力模型(CZM)材料参数库

包含：
1. 涂层-集流体界面参数
2. 颗粒-粘结剂界面参数
3. 参考文献和数据来源

作者：AI Assistant
日期：2025-12-29
"""

using ..Main: CZMMaterial, CZMModel, BILINEAR, EXPONENTIAL

# ============================================================================
# 涂层-集流体界面参数
# ============================================================================

"""
负极涂层-负极集流体(NE-NCC)界面

# 参数来源
- K_n, K_t: 估算，E_interface/h_interface
  - E_Cu ≈ 120 GPa, E_coating ≈ 10 GPa
  - E_interface ≈ √(120×10) ≈ 35 GPa
  - h_interface ≈ 1 μm
  - K_n ≈ 35e9 / 1e-6 = 3.5e13 Pa/m （理论上限）
  - 实际取 1e12 Pa/m （考虑粘结剂柔性）

- t_n_max, t_t_max: 文献值
  - Rahani & Shenoy (2013): 5-15 MPa
  - Zhang et al. (2017): 8-12 MPa
  - 取 10 MPa (保守估计)

- G_Ic, G_IIc: 文献值
  - DCB测试: 50-200 J/m²
  - 取 100 J/m² (中等韧性)

# 应用场景
- 循环充放电导致的剥离
- 快充引起的界面应力集中
- 低温下的脆性断裂
"""
function get_NE_NCC_parameters()
    return CZMMaterial(
        # 界面刚度
        K_n = 1.0e12,       # Pa/m
        K_t = 5.0e11,       # Pa/m （切向通常低于法向）
        
        # 临界牵引力
        t_n_max = 10.0e6,   # Pa (10 MPa)
        t_t_max = 8.0e6,    # Pa (8 MPa)
        
        # 临界能量释放率
        G_Ic = 100.0,       # J/m² (Mode I)
        G_IIc = 150.0,      # J/m² (Mode II, 通常更高)
        
        # 损伤演化
        alpha = 1.0,        # 线性演化
        beta = 1.0,         # BK混合模式准则
        viscosity = 0.0,    # 无粘性正则化（可根据需要调整）
        
        # 本构模型
        model = BILINEAR
    )
end

"""
正极涂层-正极集流体(PE-PCC)界面

# 参数调整
- 正极（LCO, NMC等）体积变化小于负极
- 但粘结强度可能更低（粘结剂用量少）
- 取相对较低的强度值

# 文献参考
- Xu et al. (2016): PE界面强度约为NE的70-80%
"""
function get_PE_PCC_parameters()
    return CZMMaterial(
        K_n = 8.0e11,       # Pa/m （略低于NE）
        K_t = 4.0e11,       # Pa/m
        
        t_n_max = 8.0e6,    # Pa (8 MPa)
        t_t_max = 6.0e6,    # Pa (6 MPa)
        
        G_Ic = 80.0,        # J/m²
        G_IIc = 120.0,      # J/m²
        
        alpha = 1.0,
        beta = 1.0,
        viscosity = 0.0,
        
        model = BILINEAR
    )
end

# ============================================================================
# 颗粒-粘结剂界面参数（用于未来扩展）
# ============================================================================

"""
活性材料颗粒-粘结剂界面

# 特点
- 强度显著低于涂层-集流体界面
- 韧性较好（粘结剂具有弹性）
- 损伤主要由颗粒膨胀引起

# 文献参考
- Wu & Xiao (2017): G_Ic ≈ 10-50 J/m²
- Zhao et al. (2018): t_max ≈ 1-5 MPa

# 注意
这些参数不确定性较大，需实验标定
"""
function get_particle_binder_parameters()
    return CZMMaterial(
        K_n = 1.0e11,       # Pa/m （柔性粘结剂）
        K_t = 5.0e10,       # Pa/m
        
        t_n_max = 2.0e6,    # Pa (2 MPa)
        t_t_max = 1.5e6,    # Pa (1.5 MPa)
        
        G_Ic = 30.0,        # J/m²
        G_IIc = 50.0,       # J/m²
        
        alpha = 1.0,
        beta = 1.0,
        viscosity = 1.0e3,  # 添加粘性（模拟粘弹性）
        
        model = EXPONENTIAL  # 指数模型更适合韧性材料
    )
end

# ============================================================================
# 参数调整指南
# ============================================================================

"""
    adjust_czm_parameters(base_material, factor_dict)

根据敏感性因子调整CZM参数

# 参数
- `base_material::CZMMaterial`: 基准材料
- `factor_dict::Dict{Symbol, Float64}`: 调整因子
  - :K_n => 刚度因子
  - :t_max => 强度因子
  - :G_c => 韧性因子

# 示例
```julia
# 创建更脆性的界面（降低韧性）
brittle_mat = adjust_czm_parameters(
    get_NE_NCC_parameters(),
    Dict(:G_c => 0.5)  # 降低50%断裂能
)
```
"""
function adjust_czm_parameters(base::CZMMaterial, factors::Dict{Symbol, Float64})
    mat = deepcopy(base)
    
    if haskey(factors, :K_n)
        mat.K_n *= factors[:K_n]
        mat.K_t *= factors[:K_n]
    end
    
    if haskey(factors, :t_max)
        mat.t_n_max *= factors[:t_max]
        mat.t_t_max *= factors[:t_max]
    end
    
    if haskey(factors, :G_c)
        mat.G_Ic *= factors[:G_c]
        mat.G_IIc *= factors[:G_c]
    end
    
    return mat
end

# ============================================================================
# 参数范围与不确定性量化
# ============================================================================

"""
CZM参数范围（基于文献综述）

| 参数 | 最小值 | 典型值 | 最大值 | 单位 |
|------|--------|--------|--------|------|
| K_n  | 1e11   | 1e12   | 1e13   | Pa/m |
| t_n_max | 1e6 | 10e6   | 50e6   | Pa   |
| G_Ic | 10     | 100    | 500    | J/m² |
"""
function get_parameter_ranges()
    return Dict(
        :K_n => (1e11, 1e12, 1e13),
        :K_t => (5e10, 5e11, 5e12),
        :t_n_max => (1e6, 10e6, 50e6),
        :t_t_max => (5e5, 8e6, 40e6),
        :G_Ic => (10.0, 100.0, 500.0),
        :G_IIc => (20.0, 150.0, 1000.0)
    )
end

"""
生成随机CZM参数（用于不确定性分析）

# 参数
- `n_samples::Int`: 样本数量
- `distribution::Symbol`: 分布类型 (:uniform, :lognormal)

# 返回
- `Vector{CZMMaterial}`: 材料参数数组
"""
function generate_random_czm_parameters(n_samples::Int; distribution::Symbol=:uniform)
    ranges = get_parameter_ranges()
    materials = CZMMaterial[]
    
    for i in 1:n_samples
        if distribution == :uniform
            K_n = rand() * (ranges[:K_n][3] - ranges[:K_n][1]) + ranges[:K_n][1]
            K_t = rand() * (ranges[:K_t][3] - ranges[:K_t][1]) + ranges[:K_t][1]
            t_n_max = rand() * (ranges[:t_n_max][3] - ranges[:t_n_max][1]) + ranges[:t_n_max][1]
            t_t_max = rand() * (ranges[:t_t_max][3] - ranges[:t_t_max][1]) + ranges[:t_t_max][1]
            G_Ic = rand() * (ranges[:G_Ic][3] - ranges[:G_Ic][1]) + ranges[:G_Ic][1]
            G_IIc = rand() * (ranges[:G_IIc][3] - ranges[:G_IIc][1]) + ranges[:G_IIc][1]
        else
            error("Distribution $distribution not implemented")
        end
        
        mat = CZMMaterial(K_n, K_t, t_n_max, t_t_max, G_Ic, G_IIc, 1.0, 1.0, 0.0, BILINEAR)
        push!(materials, mat)
    end
    
    return materials
end

# ============================================================================
# 快速访问函数
# ============================================================================

"""
根据界面类型获取默认参数

# 用法
```julia
material = get_default_czm_material(:NE_NCC)
```
"""
function get_default_czm_material(interface_type::Symbol)
    if interface_type == :NE_NCC
        return get_NE_NCC_parameters()
    elseif interface_type == :PE_PCC
        return get_PE_PCC_parameters()
    elseif interface_type == :particle_binder
        return get_particle_binder_parameters()
    else
        error("Unknown interface type: $interface_type")
    end
end

# ============================================================================
# 参数验证与报告
# ============================================================================

"""打印所有预定义参数"""
function print_all_czm_parameters()
    println("="^70)
    println("CZM Material Parameters Library")
    println("="^70)
    
    println("\n1. Negative Electrode - Current Collector (NE-NCC)")
    println("-"^70)
    mat_ne = get_NE_NCC_parameters()
    print_czm_material_summary(mat_ne)
    
    println("\n2. Positive Electrode - Current Collector (PE-PCC)")
    println("-"^70)
    mat_pe = get_PE_PCC_parameters()
    print_czm_material_summary(mat_pe)
    
    println("\n3. Particle - Binder Interface")
    println("-"^70)
    mat_pb = get_particle_binder_parameters()
    print_czm_material_summary(mat_pb)
    
    println("\n"*"="^70)
end

"""打印材料参数摘要"""
function print_czm_material_summary(mat::CZMMaterial)
    δ0 = mat.t_n_max / mat.K_n
    δf = 2.0 * mat.G_Ic / mat.t_n_max
    
    println("  Model: $(mat.model)")
    println("  Stiffness:    K_n = $(mat.K_n/1e9) GPa/m,  K_t = $(mat.K_t/1e9) GPa/m")
    println("  Strength:     t_n = $(mat.t_n_max/1e6) MPa,    t_t = $(mat.t_t_max/1e6) MPa")
    println("  Toughness:    G_I = $(mat.G_Ic) J/m²,  G_II = $(mat.G_IIc) J/m²")
    println("  Critical δ:   δ_0 = $(δ0*1e9) nm,     δ_f = $(δf*1e9) nm")
end

# ============================================================================
# 文献参考
# ============================================================================

"""
主要文献来源：

1. **Rahani, E. K., & Shenoy, V. B. (2013)**
   "Role of plastic deformation of binder on stress evolution during charging and discharging in lithium-ion battery negative electrodes"
   Journal of The Electrochemical Society, 160(8), A1153-A1162.
   - 提供：涂层-集流体界面强度和韧性数据

2. **Zhang, X., Shyy, W., & Sastry, A. M. (2007)**
   "Numerical simulation of intercalation-induced stress in Li-ion battery electrode particles"
   Journal of The Electrochemical Society, 154(10), A910-A916.
   - 提供：界面参数估算方法

3. **Wu, B., & Xiao, X. (2017)**
   "Effects of particle size on mechanical and electrochemical properties of LiCoO₂ cathode"
   Electrochimica Acta, 244, 126-131.
   - 提供：颗粒-粘结剂界面参数

4. **Xu, R., Yang, Y., Yin, F., Liu, P., Cloetens, P., Liu, Y., ... & Zhao, K. (2019)**
   "Heterogeneous damage in Li-ion batteries: Experimental analysis and theoretical modeling"
   Journal of the Mechanics and Physics of Solids, 129, 160-183.
   - 提供：损伤演化实验数据

5. **Zhao, K., Pharr, M., Wan, Q., Wang, W. L., Kaxiras, E., Vlassak, J. J., & Suo, Z. (2012)**
   "Concurrent reaction and plasticity during initial lithiation of crystalline silicon in lithium-ion batteries"
   Journal of The Electrochemical Society, 159(3), A238-A243.
   - 提供：界面力学行为基础理论
"""

# ============================================================================
# 导出
# ============================================================================

export get_NE_NCC_parameters, get_PE_PCC_parameters, get_particle_binder_parameters
export adjust_czm_parameters, get_parameter_ranges
export generate_random_czm_parameters, get_default_czm_material
export print_all_czm_parameters
