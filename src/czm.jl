"""
    czm.jl - 内聚力模型(Cohesive Zone Model)

用于模拟锂离子电池极片脱粘问题

理论基础：
- Dugdale-Barenblatt内聚力模型
- 牵引-分离律(Traction-Separation Law)
- 渐进损伤力学

应用场景：
1. 涂层-集流体界面脱粘
2. 颗粒-粘结剂界面失效
3. 电极内部分层

作者：AI Assistant
日期：2025-12-29
"""

# ============================================================================
# 数据结构定义
# ============================================================================

"""界面类型枚举"""
@enum InterfaceType begin
    COATING_COLLECTOR    # 涂层-集流体界面
    PARTICLE_BINDER      # 颗粒-粘结剂界面
    PARTICLE_PARTICLE    # 颗粒-颗粒界面
    LAYER_LAYER          # 层间界面
end

"""CZM本构模型类型"""
@enum CZMModel begin
    BILINEAR      # 双线性模型
    EXPONENTIAL   # 指数软化模型
    PPR           # Park-Paulino-Roesler模型
    TRAPEZOIDAL   # 梯形模型
end

"""
    CZMMaterial

CZM材料参数

# 字段
- `K_n::Float64`: 法向界面刚度 [Pa/m]
- `K_t::Float64`: 切向界面刚度 [Pa/m]
- `t_n_max::Float64`: 法向临界牵引力（强度） [Pa]
- `t_t_max::Float64`: 切向临界牵引力（强度） [Pa]
- `G_Ic::Float64`: 模式I临界能量释放率 [J/m²]
- `G_IIc::Float64`: 模式II临界能量释放率 [J/m²]
- `alpha::Float64`: 损伤演化指数 [-]
- `beta::Float64`: 混合模式指数 [-]
- `viscosity::Float64`: 粘性正则化参数 [Pa·s/m]
- `model::CZMModel`: 本构模型类型

# 参数获取来源
- 实验：剥离测试(t_max)、DCB测试(G_Ic)、ENF测试(G_IIc)
- 估算：K_n ≈ E/h_interface
- 文献：典型范围见文档

# 示例
```julia
material = CZMMaterial(
    K_n = 1e12,
    K_t = 5e11,
    t_n_max = 10e6,
    t_t_max = 8e6,
    G_Ic = 100.0,
    G_IIc = 150.0,
    alpha = 1.0,
    beta = 1.0,
    viscosity = 0.0,
    model = BILINEAR
)
```
"""
mutable struct CZMMaterial
    # 界面刚度
    K_n::Float64          # 法向刚度 [Pa/m]
    K_t::Float64          # 切向刚度 [Pa/m]
    
    # 临界牵引力（强度）
    t_n_max::Float64      # 法向强度 [Pa]
    t_t_max::Float64      # 切向强度 [Pa]
    
    # 临界能量释放率（断裂韧性）
    G_Ic::Float64         # 模式I [J/m²]
    G_IIc::Float64        # 模式II [J/m²]
    
    # 损伤演化参数
    alpha::Float64        # 损伤指数 [-]
    beta::Float64         # 混合模式指数 [-]
    viscosity::Float64    # 粘性正则化 [Pa·s/m]
    
    # 本构模型选择
    model::CZMModel
end

"""
    CZMInterface

内聚界面对象

# 字段
## 几何信息
- `id::Int64`: 界面编号
- `type::InterfaceType`: 界面类型
- `node_plus::Int64`: "+"侧节点（上侧/外侧）
- `node_minus::Int64`: "-"侧节点（下侧/内侧）
- `normal::Vector{Float64}`: 法向向量 [nx, ny]（指向+侧）
- `tangent::Vector{Float64}`: 切向向量 [tx, ty]
- `area::Float64`: 界面面积（2D中为长度×厚度） [m²]

## 材料参数
- `material::CZMMaterial`: 材料参数

## 状态变量（当前时刻）
- `D::Float64`: 损伤变量 [-], 0 ≤ D ≤ 1
- `delta_n::Float64`: 法向分离量 [m]
- `delta_t::Float64`: 切向分离量 [m]
- `t_n::Float64`: 法向牵引力 [Pa]
- `t_t::Float64`: 切向牵引力 [Pa]

## 历史变量
- `delta_max::Float64`: 最大历史等效分离量 [m]
- `G_dissipated::Float64`: 累积耗散能 [J/m²]
- `failed::Bool`: 失效标志
"""
mutable struct CZMInterface
    # 几何信息
    id::Int64
    type::InterfaceType
    node_plus::Int64
    node_minus::Int64
    normal::Vector{Float64}
    tangent::Vector{Float64}
    area::Float64
    
    # 材料参数
    material::CZMMaterial
    
    # 状态变量
    D::Float64
    delta_n::Float64
    delta_t::Float64
    t_n::Float64
    t_t::Float64
    
    # 历史变量
    delta_max::Float64
    G_dissipated::Float64
    failed::Bool
end

"""构造函数：创建新界面"""
function CZMInterface(id::Int64, type::InterfaceType, 
                     node_plus::Int64, node_minus::Int64,
                     normal::Vector{Float64}, tangent::Vector{Float64},
                     area::Float64, material::CZMMaterial)
    return CZMInterface(
        id, type, node_plus, node_minus, normal, tangent, area, material,
        0.0, 0.0, 0.0, 0.0, 0.0,  # 初始状态变量
        0.0, 0.0, false           # 初始历史变量
    )
end

# ============================================================================
# CZM本构模型
# ============================================================================

"""
    bilinear_cohesive_law(delta_n, delta_t, material, D_old)

双线性内聚力模型

# 理论
牵引-分离曲线分为三个阶段：

1. **弹性阶段** (δ < δ_0)：
   ```
   t = K·δ
   D = 0
   ```

2. **软化阶段** (δ_0 ≤ δ < δ_f)：
   ```
   t = t_max·(δ_f - δ)/(δ_f - δ_0)
   D = (δ - δ_0)/(δ_f - δ_0)
   ```

3. **完全失效** (δ ≥ δ_f)：
   ```
   t = 0
   D = 1
   ```

其中临界分离量：
- δ_0 = t_max / K （损伤起始）
- δ_f = 2·G_c / t_max （完全失效）

# 参数
- `delta_n::Float64`: 当前法向分离量 [m]
- `delta_t::Float64`: 当前切向分离量 [m]
- `material::CZMMaterial`: 材料参数
- `D_old::Float64`: 上一步损伤变量 [-]

# 返回
- `t_n::Float64`: 法向牵引力 [Pa]
- `t_t::Float64`: 切向牵引力 [Pa]
- `D::Float64`: 更新后损伤变量 [-]
- `dG::Float64`: 当前步耗散能增量 [J/m²]

# 混合模式准则
使用等效分离量：
```
δ_eff = √(⟨δ_n⟩²+ δ_t²)
```
其中 ⟨·⟩ 为Macaulay括号（仅考虑拉伸）
"""
function bilinear_cohesive_law(delta_n::Float64, delta_t::Float64, 
                               material::CZMMaterial, D_old::Float64)
    # Macaulay括号：仅拉伸产生损伤
    delta_n_pos = max(delta_n, 0.0)
    
    # 等效分离量（混合模式）
    delta_eff = sqrt(delta_n_pos^2 + delta_t^2)
    
    # 临界分离量
    delta_0 = material.t_n_max / material.K_n  # 损伤起始
    delta_f = 2.0 * material.G_Ic / material.t_n_max  # 完全失效
    
    # 防止退化
    if delta_f <= delta_0
        @warn "Invalid CZM parameters: δ_f ≤ δ_0" delta_0 delta_f
        delta_f = 2.0 * delta_0
    end
    
    # 计算损伤变量
    D = 0.0
    if delta_eff < delta_0
        # 阶段1：弹性
        D = 0.0
    elseif delta_eff < delta_f
        # 阶段2：软化
        D = (delta_eff - delta_0) / (delta_f - delta_0)
        D = max(D, D_old)  # 单调增长（不可恢复损伤）
    else
        # 阶段3：完全失效
        D = 1.0
    end
    
    # 限制范围
    D = clamp(D, 0.0, 1.0)
    
    # 计算牵引力
    if D < 1.0
        # 有效刚度（考虑损伤）
        K_n_eff = (1.0 - D) * material.K_n
        K_t_eff = (1.0 - D) * material.K_t
        
        # 牵引力
        t_n = K_n_eff * delta_n  # 允许压缩（负值）
        t_t = K_t_eff * delta_t
    else
        # 完全失效，无法承载
        t_n = 0.0
        t_t = 0.0
    end
    
    # 能量耗散增量
    if D > D_old && delta_eff > 0
        # 近似：当前步能量耗散
        dG = 0.5 * (abs(t_n * delta_n) + abs(t_t * delta_t)) * (D - D_old) / max(D, 1e-10)
    else
        dG = 0.0
    end
    
    return t_n, t_t, D, dG
end

"""
    exponential_cohesive_law(delta_n, delta_t, material, D_old)

指数软化内聚力模型

# 理论
牵引力随分离量指数衰减：
```
t = t_max·(δ/δ_0)·exp(1 - δ/δ_0)
```

损伤变量：
```
D = 1 - (t/t_max) / (δ/δ_0)
```

# 特点
- 更平滑的软化曲线
- 适合韧性材料
- 无明显的"完全失效"点
"""
function exponential_cohesive_law(delta_n::Float64, delta_t::Float64,
                                 material::CZMMaterial, D_old::Float64)
    # Macaulay括号
    delta_n_pos = max(delta_n, 0.0)
    
    # 等效分离量
    delta_eff = sqrt(delta_n_pos^2 + delta_t^2)
    
    # 特征长度
    delta_0 = material.G_Ic / material.t_n_max
    
    if delta_eff < 1e-12
        # 初始状态
        return 0.0, 0.0, 0.0, 0.0
    end
    
    # 无量纲分离
    xi = delta_eff / delta_0
    
    # 等效牵引力
    t_eff = material.t_n_max * xi * exp(1.0 - xi)
    
    # 损伤变量
    if xi > 0
        D = 1.0 - t_eff / (material.t_n_max * xi)
    else
        D = 0.0
    end
    D = max(D, D_old)  # 单调增长
    D = clamp(D, 0.0, 1.0)
    
    # 分解到法向和切向
    if delta_eff > 0
        t_n = t_eff * delta_n_pos / delta_eff
        t_t = t_eff * delta_t / delta_eff
    else
        t_n = 0.0
        t_t = 0.0
    end
    
    # 能量耗散
    if D > D_old
        dG = t_eff * delta_eff * (D - D_old) / max(D, 1e-10)
    else
        dG = 0.0
    end
    
    return t_n, t_t, D, dG
end

# ============================================================================
# 界面状态更新
# ============================================================================

"""
    compute_interface_separation(interface, U_global)

从全局位移向量计算界面分离量

# 参数
- `interface::CZMInterface`: 界面对象
- `U_global::Vector{Float64}`: 全局位移向量 [u1_x, u1_y, u2_x, u2_y, ...]

# 返回
- `delta_n::Float64`: 法向分离量 [m]
- `delta_t::Float64`: 切向分离量 [m]

# 原理
1. 提取两侧节点的位移
2. 计算位移跳跃 Δu = u⁺ - u⁻
3. 投影到局部坐标系（法向n，切向t）
"""
function compute_interface_separation(interface::CZMInterface, U_global::Vector{Float64})
    # 提取两侧节点位移（全局坐标系）
    node_plus = interface.node_plus
    node_minus = interface.node_minus
    
    u_plus = [U_global[2*node_plus-1], U_global[2*node_plus]]
    u_minus = [U_global[2*node_minus-1], U_global[2*node_minus]]
    
    # 位移跳跃
    Delta_u = u_plus - u_minus
    
    # 投影到局部坐标系
    n = interface.normal
    t = interface.tangent
    
    delta_n = dot(Delta_u, n)  # 法向分离（可正可负）
    delta_t = dot(Delta_u, t)  # 切向分离
    
    return delta_n, delta_t
end

"""
    update_czm_interface!(interface, U_global, dt)

更新单个界面的CZM状态

# 参数
- `interface::CZMInterface`: 待更新的界面（原地修改）
- `U_global::Vector{Float64}`: 全局位移向量
- `dt::Float64`: 时间步长 [s]（用于粘性正则化）

# 更新内容
- 分离量 (delta_n, delta_t)
- 牵引力 (t_n, t_t)
- 损伤变量 D
- 累积耗散能 G_dissipated
- 失效标志 failed
"""
function update_czm_interface!(interface::CZMInterface, U_global::Vector{Float64}, dt::Float64)
    # 计算当前分离量
    delta_n, delta_t = compute_interface_separation(interface, U_global)
    
    # 保存旧损伤状态
    D_old = interface.D
    
    # 调用本构模型
    if interface.material.model == BILINEAR
        t_n, t_t, D, dG = bilinear_cohesive_law(
            delta_n, delta_t, interface.material, D_old
        )
    elseif interface.material.model == EXPONENTIAL
        t_n, t_t, D, dG = exponential_cohesive_law(
            delta_n, delta_t, interface.material, D_old
        )
    else
        error("Unknown CZM model: $(interface.material.model)")
    end
    
    # 粘性正则化（可选，防止数值振荡）
    if interface.material.viscosity > 0
        visc = interface.material.viscosity
        t_n += visc * (delta_n - interface.delta_n) / dt
        t_t += visc * (delta_t - interface.delta_t) / dt
    end
    
    # 更新状态变量
    interface.delta_n = delta_n
    interface.delta_t = delta_t
    interface.t_n = t_n
    interface.t_t = t_t
    interface.D = D
    interface.G_dissipated += dG
    
    # 更新最大历史分离
    delta_n_pos = max(delta_n, 0.0)
    delta_eff = sqrt(delta_n_pos^2 + delta_t^2)
    interface.delta_max = max(interface.delta_max, delta_eff)
    
    # 失效判定
    if D ≥ 0.99
        interface.failed = true
    end
end

"""
    update_all_czm_interfaces!(interfaces, U_global, dt)

批量更新所有界面的CZM状态

# 参数
- `interfaces::Vector{CZMInterface}`: 界面列表
- `U_global::Vector{Float64}`: 全局位移向量
- `dt::Float64`: 时间步长 [s]

# 并行化
可以并行处理各界面（无相互依赖）
"""
function update_all_czm_interfaces!(interfaces::Vector{CZMInterface}, 
                                   U_global::Vector{Float64}, 
                                   dt::Float64)
    for interface in interfaces
        update_czm_interface!(interface, U_global, dt)
    end
end

# ============================================================================
# 与力学模块的接口
# ============================================================================

"""
    extract_displacement_from_variables(variables)

从variables字典提取全局位移向量

# 参数
- `variables::Dict`: variables字典（来自mechanical.jl）

# 返回
- `U_global::Vector{Float64}`: 交错排列的位移向量
  格式：[u1_x, u1_y, u2_x, u2_y, ..., un_x, un_y]

# 输入格式
variables应包含：
- "displacement x": x方向位移数组
- "displacement y": y方向位移数组
"""
function extract_displacement_from_variables(variables::Dict{String, <:Union{Array{Float64},Float64}})
    if !haskey(variables, "displacement x") || !haskey(variables, "displacement y")
        error("Missing displacement fields in variables")
    end
    
    u_x = variables["displacement x"]
    u_y = variables["displacement y"]
    
    nnode = length(u_x)
    U_global = zeros(Float64, 2*nnode)
    
    for i in 1:nnode
        U_global[2*i-1] = u_x[i]
        U_global[2*i] = u_y[i]
    end
    
    return U_global
end

"""
    extract_stress_from_variables(variables)

从variables字典提取应力场

# 返回
Dict包含：
- "xx": σ_xx
- "yy": σ_yy
- "xy": σ_xy
- "vm": Von Mises应力
"""
function extract_stress_from_variables(variables::Dict{String, <:Union{Array{Float64},Float64}})
    stress = Dict{String, Vector{Float64}}()
    
    if haskey(variables, "diffusion stress xx")
        stress["xx"] = variables["diffusion stress xx"]
        stress["yy"] = variables["diffusion stress yy"]
        stress["xy"] = variables["diffusion stress xy"]
        stress["vm"] = variables["diffusion stress vonMises"]
    else
        @warn "Stress fields not found in variables"
    end
    
    return stress
end

"""
    compute_interface_traction_from_stress(interface, stress, mesh)

从单元应力场计算界面牵引力（用于验证）

# 原理
在界面处：
- 法向牵引力：t_n = n·σ·n
- 切向牵引力：t_t = t·σ·n

# 注意
这是从连续体应力外推到界面，用于损伤起始判据。
与CZM本构计算的牵引力不同（CZM是界面内禀性质）。
"""
function compute_interface_traction_from_stress(interface::CZMInterface, 
                                               stress::Dict{String, Vector{Float64}},
                                               element_id::Int64)
    # 提取单元应力分量
    σ_xx = stress["xx"][element_id]
    σ_yy = stress["yy"][element_id]
    σ_xy = stress["xy"][element_id]
    
    # 界面法向和切向
    n = interface.normal
    t = interface.tangent
    nx, ny = n[1], n[2]
    tx, ty = t[1], t[2]
    
    # 应力张量
    # σ = [σ_xx  σ_xy]
    #     [σ_xy  σ_yy]
    
    # 法向牵引力：t_n = n·σ·n
    t_n = σ_xx * nx^2 + σ_yy * ny^2 + 2*σ_xy * nx * ny
    
    # 切向牵引力：t_t = t·σ·n
    t_t = σ_xx * tx * nx + σ_yy * ty * ny + σ_xy * (tx*ny + ty*nx)
    
    return t_n, t_t
end

# ============================================================================
# 后处理与统计
# ============================================================================

"""
    compute_czm_statistics(interfaces)

计算CZM统计信息

# 返回
Dict包含：
- n_total: 总界面数
- n_failed: 失效界面数
- D_max: 最大损伤
- D_mean: 平均损伤
- G_total: 总耗散能
"""
function compute_czm_statistics(interfaces::Vector{CZMInterface})
    n_total = length(interfaces)
    n_failed = count(iface -> iface.failed, interfaces)
    
    D_values = [iface.D for iface in interfaces]
    D_max = maximum(D_values)
    D_mean = mean(D_values)
    
    G_total = sum(iface.G_dissipated * iface.area for iface in interfaces)
    
    return Dict(
        "n_total" => n_total,
        "n_failed" => n_failed,
        "D_max" => D_max,
        "D_mean" => D_mean,
        "G_total" => G_total,
        "failure_percentage" => 100.0 * n_failed / n_total
    )
end

"""
    write_czm_results!(variables, interfaces)

将CZM结果写入variables字典

# 添加字段
- "czm damage": 损伤变量数组
- "czm traction normal": 法向牵引力数组
- "czm traction tangent": 切向牵引力数组
- "czm failed": 失效标志数组
- "czm dissipated energy": 累积耗散能数组
"""
function write_czm_results!(variables::Dict{String, <:Union{Array{Float64},Float64}}, 
                           interfaces::Vector{CZMInterface})
    n = length(interfaces)
    
    variables["czm damage"] = [iface.D for iface in interfaces]
    variables["czm traction normal"] = [iface.t_n for iface in interfaces]
    variables["czm traction tangent"] = [iface.t_t for iface in interfaces]
    variables["czm failed"] = Float64.([iface.failed for iface in interfaces])
    variables["czm dissipated energy"] = [iface.G_dissipated for iface in interfaces]
    variables["czm separation normal"] = [iface.delta_n for iface in interfaces]
    variables["czm separation tangent"] = [iface.delta_t for iface in interfaces]
    
    # 统计信息
    stats = compute_czm_statistics(interfaces)
    for (key, value) in stats
        variables["czm_stat_$key"] = value
    end
end

# ============================================================================
# 辅助函数
# ============================================================================

"""检查CZM参数的物理合理性"""
function validate_czm_material(material::CZMMaterial)
    valid = true
    
    if material.K_n ≤ 0 || material.K_t ≤ 0
        @warn "Interface stiffness must be positive" material.K_n material.K_t
        valid = false
    end
    
    if material.t_n_max ≤ 0 || material.t_t_max ≤ 0
        @warn "Strength must be positive" material.t_n_max material.t_t_max
        valid = false
    end
    
    if material.G_Ic ≤ 0 || material.G_IIc ≤ 0
        @warn "Fracture energy must be positive" material.G_Ic material.G_IIc
        valid = false
    end
    
    # 检查一致性
    delta_0 = material.t_n_max / material.K_n
    delta_f = 2.0 * material.G_Ic / material.t_n_max
    
    if delta_f ≤ delta_0
        @warn "Inconsistent parameters: δ_f ≤ δ_0" delta_0 delta_f
        @warn "This implies G_Ic < 0.5*t_max²/K_n"
        valid = false
    end
    
    return valid
end

"""打印CZM材料参数"""
function print_czm_material(material::CZMMaterial, name::String="")
    if name != ""
        println("CZM Material: $name")
    else
        println("CZM Material Parameters:")
    end
    println("  Model: $(material.model)")
    println("  Interface Stiffness:")
    println("    K_n = $(material.K_n) Pa/m")
    println("    K_t = $(material.K_t) Pa/m")
    println("  Strength:")
    println("    t_n_max = $(material.t_n_max*1e-6) MPa")
    println("    t_t_max = $(material.t_t_max*1e-6) MPa")
    println("  Fracture Energy:")
    println("    G_Ic = $(material.G_Ic) J/m²")
    println("    G_IIc = $(material.G_IIc) J/m²")
    
    # 计算派生量
    delta_0 = material.t_n_max / material.K_n
    delta_f = 2.0 * material.G_Ic / material.t_n_max
    println("  Critical Separations:")
    println("    δ_0 = $(delta_0*1e9) nm (damage initiation)")
    println("    δ_f = $(delta_f*1e9) nm (complete failure)")
end

"""打印CZM界面状态"""
function print_czm_interface(interface::CZMInterface)
    println("CZM Interface #$(interface.id):")
    println("  Type: $(interface.type)")
    println("  Nodes: $(interface.node_plus) ← → $(interface.node_minus)")
    println("  Normal: [$(interface.normal[1]), $(interface.normal[2])]")
    println("  State:")
    println("    D = $(interface.D)")
    println("    δ_n = $(interface.delta_n*1e9) nm")
    println("    δ_t = $(interface.delta_t*1e9) nm")
    println("    t_n = $(interface.t_n*1e-6) MPa")
    println("    t_t = $(interface.t_t*1e-6) MPa")
    println("    G_diss = $(interface.G_dissipated) J/m²")
    println("    Failed: $(interface.failed)")
end

# ============================================================================
# 导出符号
# ============================================================================

export InterfaceType, CZMModel
export COATING_COLLECTOR, PARTICLE_BINDER, PARTICLE_PARTICLE, LAYER_LAYER
export BILINEAR, EXPONENTIAL, PPR, TRAPEZOIDAL
export CZMMaterial, CZMInterface
export bilinear_cohesive_law, exponential_cohesive_law
export compute_interface_separation
export update_czm_interface!, update_all_czm_interfaces!
export extract_displacement_from_variables, extract_stress_from_variables
export compute_interface_traction_from_stress
export compute_czm_statistics, write_czm_results!
export validate_czm_material, print_czm_material, print_czm_interface
