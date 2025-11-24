# ================================================================================
# 内聚力模型与界面接触理论模块
# Cohesive Zone Model (CZM) and Interface Contact Mechanics
#
# 功能：
# 1. 内聚力模型：双线性牵引-分离定律，混合模式断裂准则
# 2. 界面接触：罚函数法，Hertz接触，Coulomb摩擦
# 3. 多物理场耦合：力学-电化学-热耦合
#
# 作者：AI Assistant
# 日期：2025-11-24
# ================================================================================

using LinearAlgebra, Statistics

# ================================================================================
# 数据结构定义
# ================================================================================

"""
    CohesiveInterface

内聚力界面数据结构，描述材料界面的渐进损伤和脱粘过程

# 字段
## 几何拓扑
- `element_pairs::Vector{Tuple{Int,Int}}` - 界面单元对 (负极单元ID, 正极单元ID)
- `interface_area::Vector{Float64}` - 界面面积 [m²]
- `normals::Matrix{Float64}` - 法向量 (n_pairs × 3)

## 材料参数
- `K_n::Float64` - 法向界面刚度 [Pa/m]
- `K_t::Float64` - 切向界面刚度 [Pa/m]
- `T_n_max::Float64` - 法向最大牵引应力 [Pa]
- `T_t_max::Float64` - 切向最大牵引应力 [Pa]
- `δ_n_0::Float64` - 法向损伤起始位移 [m]
- `δ_t_0::Float64` - 切向损伤起始位移 [m]
- `δ_n_f::Float64` - 法向完全失效位移 [m]
- `δ_t_f::Float64` - 切向完全失效位移 [m]
- `Γ_n::Float64` - 法向断裂能（I型） [J/m²]
- `Γ_t::Float64` - 切向断裂能（II型） [J/m²]
- `η::Float64` - Benzeggagh-Kenane准则指数 (通常1.5-2.5)

## 状态变量（历史相关）
- `damage::Vector{Float64}` - 损伤变量 D ∈ [0, 1]
- `δ_eq_max::Vector{Float64}` - 历史最大等效分离位移 [m]
- `is_failed::Vector{Bool}` - 完全失效标志
- `energy_dissipated::Vector{Float64}` - 累积耗散能量 [J]
"""
mutable struct CohesiveInterface
    # 几何拓扑
    element_pairs::Vector{Tuple{Int,Int}}
    interface_area::Vector{Float64}
    normals::Matrix{Float64}
    
    # 材料参数
    K_n::Float64
    K_t::Float64
    T_n_max::Float64
    T_t_max::Float64
    δ_n_0::Float64
    δ_t_0::Float64
    δ_n_f::Float64
    δ_t_f::Float64
    Γ_n::Float64
    Γ_t::Float64
    η::Float64
    
    # 状态变量
    damage::Vector{Float64}
    δ_eq_max::Vector{Float64}
    is_failed::Vector{Bool}
    energy_dissipated::Vector{Float64}
end

"""
    ContactInterface

界面接触数据结构，描述材料界面的法向接触和切向摩擦行为

# 字段
## 几何拓扑
- `master_surface::Vector{Int}` - 主表面节点ID列表
- `slave_nodes::Vector{Int}` - 从节点ID列表
- `contact_pairs::Vector{Tuple{Int,Int}}` - 潜在接触对 (主表面ID, 从节点ID)
- `normals::Matrix{Float64}` - 接触法向量 (n_pairs × 3)

## 接触参数
- `K_n_contact::Float64` - 法向罚刚度 [Pa/m]
- `K_t_contact::Float64` - 切向罚刚度 [Pa/m]
- `p_max::Float64` - 最大接触压力 [Pa]
- `g_penetration::Float64` - 允许穿透深度（特征间隙） [m]

## 摩擦参数
- `μ_s::Float64` - 静摩擦系数
- `μ_k::Float64` - 动摩擦系数 (μ_k ≤ μ_s)
- `v_critical::Float64` - 临界滑移速度 [m/s]

## 状态变量
- `gap::Vector{Float64}` - 间隙函数 g [m] (< 0 表示接触)
- `p_contact::Vector{Float64}` - 接触压力 [Pa]
- `τ_friction::Vector{Float64}` - 摩擦应力 [Pa]
- `is_contact::Vector{Bool}` - 接触激活标志
- `is_sliding::Vector{Bool}` - 滑移标志
- `contact_area::Vector{Float64}` - 接触面积 [m²]
- `cumulative_slip::Vector{Float64}` - 累积滑移距离 [m]
"""
mutable struct ContactInterface
    # 几何拓扑
    master_surface::Vector{Int}
    slave_nodes::Vector{Int}
    contact_pairs::Vector{Tuple{Int,Int}}
    normals::Matrix{Float64}
    
    # 接触参数
    K_n_contact::Float64
    K_t_contact::Float64
    p_max::Float64
    g_penetration::Float64
    
    # 摩擦参数
    μ_s::Float64
    μ_k::Float64
    v_critical::Float64
    
    # 状态变量
    gap::Vector{Float64}
    p_contact::Vector{Float64}
    τ_friction::Vector{Float64}
    is_contact::Vector{Bool}
    is_sliding::Vector{Bool}
    contact_area::Vector{Float64}
    cumulative_slip::Vector{Float64}
end

# ================================================================================
# 内聚力模型 (CZM) 核心函数
# ================================================================================

"""
    init_cohesive_interface(case::Case, interface_type::Symbol; kwargs...)

初始化内聚力界面

# 参数
- `case::Case` - JuBat案例对象
- `interface_type::Symbol` - 界面类型
  - `:particle_binder` - 颗粒-粘结剂界面
  - `:electrode_separator` - 电极-隔膜界面
  - `:internal_crack` - 内部裂纹

# 关键字参数
可选覆盖材料参数，例如：
- `T_n_max::Float64` - 法向强度 [Pa]
- `Γ_n::Float64` - 法向断裂能 [J/m²]
- `η::Float64` - BK准则指数

# 返回
- `CohesiveInterface` - 初始化的内聚力界面对象
"""
function init_cohesive_interface(case, interface_type::Symbol; kwargs...)
    param = case.param_dim
    
    # 默认参数（根据界面类型）
    if interface_type == :particle_binder
        # 石墨负极-粘结剂典型参数
        K_n = get(kwargs, :K_n, 1e12)          # Pa/m
        K_t = get(kwargs, :K_t, 5e11)          # Pa/m
        T_n_max = get(kwargs, :T_n_max, 10e6)  # 10 MPa
        T_t_max = get(kwargs, :T_t_max, 5e6)   # 5 MPa
        Γ_n = get(kwargs, :Γ_n, 100.0)         # J/m²
        Γ_t = get(kwargs, :Γ_t, 150.0)         # J/m²
        η = get(kwargs, :η, 2.0)
        
    elseif interface_type == :electrode_separator
        # 电极-隔膜界面参数
        K_n = get(kwargs, :K_n, 5e11)
        K_t = get(kwargs, :K_t, 2e11)
        T_n_max = get(kwargs, :T_n_max, 15e6)  # 15 MPa
        T_t_max = get(kwargs, :T_t_max, 8e6)   # 8 MPa
        Γ_n = get(kwargs, :Γ_n, 150.0)
        Γ_t = get(kwargs, :Γ_t, 200.0)
        η = get(kwargs, :η, 2.0)
        
    else  # :internal_crack 或其他
        K_n = get(kwargs, :K_n, 1e12)
        K_t = get(kwargs, :K_t, 5e11)
        T_n_max = get(kwargs, :T_n_max, 20e6)  # 20 MPa
        T_t_max = get(kwargs, :T_t_max, 10e6)  # 10 MPa
        Γ_n = get(kwargs, :Γ_n, 200.0)
        Γ_t = get(kwargs, :Γ_t, 300.0)
        η = get(kwargs, :η, 2.0)
    end
    
    # 从断裂能推导失效位移
    δ_n_f = 2.0 * Γ_n / T_n_max
    δ_t_f = 2.0 * Γ_t / T_t_max
    
    # 从刚度推导损伤起始位移
    δ_n_0 = T_n_max / K_n
    δ_t_0 = T_t_max / K_t
    
    # 构建界面单元对（简化：假设负极和正极颗粒一一对应）
    # 实际应用中需要根据网格拓扑构建
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        n_particles = 1  # SPM/SPMe只有一个等效颗粒
    else  # P2D
        mesh_ne = case.mesh["negative electrode"]
        n_particles = mesh_ne.nlen
    end
    
    element_pairs = [(i, i) for i in 1:n_particles]
    
    # 界面面积（简化：基于颗粒表面积）
    rs_n = param.NE.rs
    rs_p = param.PE.rs
    interface_area = fill(4π * (rs_n^2 + rs_p^2) / 2, n_particles)
    
    # 法向量（简化：径向方向）
    normals = zeros(Float64, n_particles, 3)
    normals[:, 1] .= 1.0  # x方向
    
    # 初始化状态变量
    damage = zeros(Float64, n_particles)
    δ_eq_max = zeros(Float64, n_particles)
    is_failed = falses(n_particles)
    energy_dissipated = zeros(Float64, n_particles)
    
    return CohesiveInterface(
        element_pairs, interface_area, normals,
        K_n, K_t, T_n_max, T_t_max,
        δ_n_0, δ_t_0, δ_n_f, δ_t_f,
        Γ_n, Γ_t, η,
        damage, δ_eq_max, is_failed, energy_dissipated
    )
end

"""
    compute_cohesive_traction(coh::CohesiveInterface, δ_n, δ_t)

计算内聚力牵引力（双线性牵引-分离定律）

# 理论
双线性定律分三个阶段：
1. 弹性加载 (0 ≤ δ ≤ δ₀): T = K·δ, D = 0
2. 损伤软化 (δ₀ < δ < δ_f): T = T_max·(δ_f - δ)/(δ_f - δ₀), D ∈ (0,1)
3. 完全失效 (δ ≥ δ_f): T = 0, D = 1

# 参数
- `coh::CohesiveInterface` - 内聚力界面对象
- `δ_n::Vector{Float64}` - 法向分离位移 [m]
- `δ_t::Vector{Float64}` - 切向分离位移 [m]

# 返回
- `T_n::Vector{Float64}` - 法向牵引力 [Pa]
- `T_t::Vector{Float64}` - 切向牵引力 [Pa]
- `D_new::Vector{Float64}` - 更新后的损伤变量
"""
function compute_cohesive_traction(
    coh::CohesiveInterface,
    δ_n::Vector{Float64},
    δ_t::Vector{Float64}
)
    ne = length(δ_n)
    T_n = zeros(Float64, ne)
    T_t = zeros(Float64, ne)
    D_new = copy(coh.damage)
    
    @inbounds for i in 1:ne
        # 等效分离位移（仅受拉有效，Macaulay括号）
        δ_n_eff = max(0.0, δ_n[i])
        δ_eq = sqrt(δ_n_eff^2 + δ_t[i]^2)
        
        # 更新历史最大值（不可逆）
        δ_eq_max_i = max(coh.δ_eq_max[i], δ_eq)
        
        # 判断损伤阶段
        if δ_eq_max_i <= coh.δ_n_0
            # ============= 阶段1：弹性加载 =============
            T_n[i] = coh.K_n * δ_n_eff
            T_t[i] = coh.K_t * δ_t[i]
            D_new[i] = 0.0
            
        elseif δ_eq_max_i < coh.δ_n_f
            # ============= 阶段2：损伤软化 =============
            # 损伤变量 (Camacho-Ortiz 定义)
            D_new[i] = (coh.δ_n_f / δ_eq_max_i) * 
                       (δ_eq_max_i - coh.δ_n_0) / (coh.δ_n_f - coh.δ_n_0)
            
            # 有效牵引力（降刚度）
            if δ_eq > 1e-16
                # 混合模式分配
                T_n[i] = (1.0 - D_new[i]) * coh.T_n_max * (δ_n_eff / δ_eq)
                T_t[i] = (1.0 - D_new[i]) * coh.T_t_max * (δ_t[i] / δ_eq)
            else
                T_n[i] = 0.0
                T_t[i] = 0.0
            end
            
        else
            # ============= 阶段3：完全失效 =============
            T_n[i] = 0.0
            T_t[i] = 0.0
            D_new[i] = 1.0
            coh.is_failed[i] = true
        end
    end
    
    return T_n, T_t, D_new
end

"""
    update_cohesive_damage!(coh::CohesiveInterface, δ_n, δ_t)

更新内聚力损伤状态（用于时间推进）
"""
function update_cohesive_damage!(coh::CohesiveInterface, δ_n, δ_t)
    T_n, T_t, D_new = compute_cohesive_traction(coh, δ_n, δ_t)
    
    # 更新状态变量
    @inbounds for i in eachindex(D_new)
        δ_n_eff = max(0.0, δ_n[i])
        δ_eq = sqrt(δ_n_eff^2 + δ_t[i]^2)
        δ_eq_old = coh.δ_eq_max[i]
        
        # 更新历史最大值
        coh.δ_eq_max[i] = max(coh.δ_eq_max[i], δ_eq)
        
        # 计算增量耗散能（梯形积分）
        if D_new[i] > coh.damage[i] && δ_eq > δ_eq_old
            dΔ = δ_eq - δ_eq_old
            T_avg = 0.5 * (sqrt(T_n[i]^2 + T_t[i]^2) + 
                          coh.K_n * δ_eq_old * (1.0 - coh.damage[i]))
            coh.energy_dissipated[i] += T_avg * dΔ * coh.interface_area[i]
        end
        
        coh.damage[i] = D_new[i]
    end
    
    return coh
end

"""
    cohesive_stiffness_matrix(coh::CohesiveInterface, δ_n, δ_t)

计算内聚力切线刚度矩阵
"""
function cohesive_stiffness_matrix(coh::CohesiveInterface, δ_n, δ_t)
    ne = length(δ_n)
    K_nn = zeros(Float64, ne)
    K_tt = zeros(Float64, ne)
    
    @inbounds for i in 1:ne
        D = coh.damage[i]
        δ_n_eff = max(0.0, δ_n[i])
        δ_eq = sqrt(δ_n_eff^2 + δ_t[i]^2)
        
        if δ_eq <= coh.δ_n_0
            # 弹性刚度
            K_nn[i] = coh.K_n
            K_tt[i] = coh.K_t
            
        elseif δ_eq < coh.δ_n_f && δ_eq > 1e-16
            # 降刚度（切线模量）
            factor = (1.0 - D) * coh.T_n_max / δ_eq
            K_nn[i] = factor * (1.0 - δ_n_eff^2 / δ_eq^2)
            K_tt[i] = factor * (1.0 - δ_t[i]^2 / δ_eq^2)
            
        else
            # 完全失效（零刚度）
            K_nn[i] = 0.0
            K_tt[i] = 0.0
        end
    end
    
    return K_nn, K_tt
end

"""
    compute_BK_fracture_energy(coh::CohesiveInterface, δ_n, δ_t)

计算混合模式断裂能（Benzeggagh-Kenane准则）

# 理论
Γ_c(β) = Γ_n + (Γ_t - Γ_n)·β^η
其中 β = G_t/(G_n + G_t) 为模式比
"""
function compute_BK_fracture_energy(
    coh::CohesiveInterface,
    δ_n::Float64,
    δ_t::Float64
)
    δ_n_eff = max(0.0, δ_n)
    
    # 能量释放率
    G_n = 0.5 * coh.K_n * δ_n_eff^2
    G_t = 0.5 * coh.K_t * δ_t^2
    G_total = G_n + G_t
    
    if G_total < 1e-16
        return coh.Γ_n
    end
    
    # 模式比
    β = G_t / G_total
    
    # BK准则
    Γ_c = coh.Γ_n + (coh.Γ_t - coh.Γ_n) * β^coh.η
    
    return Γ_c
end

# ================================================================================
# 界面接触理论 核心函数
# ================================================================================

"""
    init_contact_interface(case::Case, contact_type::Symbol; kwargs...)

初始化界面接触

# 参数
- `case::Case` - JuBat案例对象
- `contact_type::Symbol` - 接触类型
  - `:particle_particle` - 颗粒-颗粒接触
  - `:electrode_collector` - 电极-集流体接触

# 关键字参数
- `μ_s::Float64` - 静摩擦系数
- `μ_k::Float64` - 动摩擦系数
- `p_max::Float64` - 最大接触压力 [Pa]

# 返回
- `ContactInterface` - 初始化的接触界面对象
"""
function init_contact_interface(case, contact_type::Symbol; kwargs...)
    param = case.param_dim
    
    # 默认参数
    if contact_type == :particle_particle
        # 颗粒-颗粒接触参数
        E_n = hasproperty(param.NE, :E) ? param.NE.E : 1e9
        ν_n = hasproperty(param.NE, :nu) ? param.NE.nu : 0.3
        E_star = E_n / (1.0 - ν_n^2)
        
        K_n_contact = get(kwargs, :K_n_contact, E_star)
        K_t_contact = get(kwargs, :K_t_contact, 0.5 * E_star)
        p_max = get(kwargs, :p_max, 50e6)       # 50 MPa
        g_penetration = get(kwargs, :g_penetration, 1e-8)  # 10 nm
        μ_s = get(kwargs, :μ_s, 0.3)
        μ_k = get(kwargs, :μ_k, 0.2)
        v_critical = get(kwargs, :v_critical, 0.01)  # 1 cm/s
        
    else  # :electrode_collector
        # 电极-集流体接触参数
        K_n_contact = get(kwargs, :K_n_contact, 1e11)
        K_t_contact = get(kwargs, :K_t_contact, 5e10)
        p_max = get(kwargs, :p_max, 100e6)      # 100 MPa
        g_penetration = get(kwargs, :g_penetration, 1e-8)
        μ_s = get(kwargs, :μ_s, 0.2)
        μ_k = get(kwargs, :μ_k, 0.15)
        v_critical = get(kwargs, :v_critical, 0.01)
    end
    
    # 构建接触对（简化）
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        n_contacts = 1
    else
        mesh_ne = case.mesh["negative electrode"]
        n_contacts = mesh_ne.nlen
    end
    
    master_surface = collect(1:n_contacts)
    slave_nodes = collect(1:n_contacts)
    contact_pairs = [(i, i) for i in 1:n_contacts]
    
    # 法向量
    normals = zeros(Float64, n_contacts, 3)
    normals[:, 1] .= 1.0
    
    # 初始化状态变量
    gap = fill(1e-9, n_contacts)  # 初始小间隙
    p_contact = zeros(Float64, n_contacts)
    τ_friction = zeros(Float64, n_contacts)
    is_contact = falses(n_contacts)
    is_sliding = falses(n_contacts)
    
    # 接触面积（简化）
    rs_n = param.NE.rs
    contact_area = fill(π * rs_n^2, n_contacts)
    
    cumulative_slip = zeros(Float64, n_contacts)
    
    return ContactInterface(
        master_surface, slave_nodes, contact_pairs, normals,
        K_n_contact, K_t_contact, p_max, g_penetration,
        μ_s, μ_k, v_critical,
        gap, p_contact, τ_friction,
        is_contact, is_sliding, contact_area, cumulative_slip
    )
end

"""
    detect_contact(contact::ContactInterface, u_master, u_slave, x_master, x_slave)

检测接触状态

# 返回
- `gap::Vector{Float64}` - 间隙函数 g [m]
- `is_contact::Vector{Bool}` - 接触激活标志
"""
function detect_contact(
    contact::ContactInterface,
    u_master::Vector{Float64},
    u_slave::Vector{Float64},
    x_master::Matrix{Float64},
    x_slave::Matrix{Float64}
)
    n_pairs = length(contact.contact_pairs)
    gap = zeros(Float64, n_pairs)
    is_contact = falses(n_pairs)
    
    @inbounds for i in 1:n_pairs
        (master_id, slave_id) = contact.contact_pairs[i]
        
        # 当前位置
        x_m = x_master[master_id, :] + u_master[master_id, :]
        x_s = x_slave[slave_id, :] + u_slave[slave_id, :]
        
        # 间隙函数（投影到法向）
        n = contact.normals[i, :]
        gap[i] = dot(x_s - x_m, n)
        
        # 接触判断
        is_contact[i] = gap[i] < 0.0
    end
    
    return gap, is_contact
end

"""
    compute_contact_pressure(contact::ContactInterface)

计算接触压力（指数软化罚函数法）

# 理论
p = p_max·[1 - exp(g/g₀)]  , g < 0 (接触)
p = 0                       , g ≥ 0 (分离)
"""
function compute_contact_pressure(contact::ContactInterface)
    n_pairs = length(contact.gap)
    p = zeros(Float64, n_pairs)
    
    @inbounds for i in 1:n_pairs
        if contact.is_contact[i]
            g = contact.gap[i]
            g0 = contact.g_penetration
            
            # 指数软化（避免刚度突变）
            p[i] = contact.p_max * (1.0 - exp(g / g0))
            
            # 限制最大值
            p[i] = min(p[i], contact.p_max)
        else
            p[i] = 0.0
        end
    end
    
    return p
end

"""
    compute_friction_stress(contact::ContactInterface, v_tangential)

计算摩擦应力（Coulomb定律 + 正则化）

# 理论
τ = μ(v_t)·p·v_t/v_reg
μ(v_t) = μ_k + (μ_s - μ_k)·exp(-|v_t|/v_c)
v_reg = √(v_t² + v_c²)
"""
function compute_friction_stress(
    contact::ContactInterface,
    v_tangential::Vector{Float64}
)
    n_pairs = length(contact.p_contact)
    τ = zeros(Float64, n_pairs)
    is_sliding = falses(n_pairs)
    
    @inbounds for i in 1:n_pairs
        if !contact.is_contact[i]
            τ[i] = 0.0
            continue
        end
        
        p = contact.p_contact[i]
        v_t = v_tangential[i]
        v_c = contact.v_critical
        
        # 速度依赖摩擦系数
        μ = contact.μ_k + (contact.μ_s - contact.μ_k) * exp(-abs(v_t) / v_c)
        
        # 正则化速度（避免除零）
        v_reg = sqrt(v_t^2 + v_c^2)
        
        # 摩擦应力
        τ[i] = μ * p * v_t / v_reg
        
        # 判断滑移状态
        is_sliding[i] = abs(v_t) > v_c
    end
    
    return τ, is_sliding
end

"""
    compute_friction_dissipation(contact::ContactInterface, Δu_tangential, dt)

计算摩擦耗散能

# 理论
W_friction = ∫ τ·v_t dt = Σ τ_avg·|Δu_t|·A_contact
"""
function compute_friction_dissipation(
    contact::ContactInterface,
    Δu_tangential::Vector{Float64},
    dt::Float64
)
    W_friction = 0.0
    
    @inbounds for i in eachindex(contact.p_contact)
        if contact.is_contact[i] && contact.is_sliding[i]
            # 平均摩擦力
            τ_avg = contact.τ_friction[i]
            
            # 功增量
            dW = τ_avg * abs(Δu_tangential[i]) * contact.contact_area[i]
            W_friction += dW
            
            # 累积滑移
            contact.cumulative_slip[i] += abs(Δu_tangential[i])
        end
    end
    
    return W_friction
end

"""
    hertz_contact_pressure(R, F_normal, E_star)

计算Hertz接触压力分布（球-平面）

# 理论
接触半径：a = (3FR/4E*)^(1/3)
最大压力：p_max = (6FE*²/π³R²)^(1/3)
压力分布：p(r) = p_max·√(1 - (r/a)²)  , r ≤ a

# 返回
- `a::Float64` - 接触半径 [m]
- `p_max::Float64` - 最大压力 [Pa]
- `p_hertz::Function` - 压力分布函数 p(r)
"""
function hertz_contact_pressure(
    R::Float64,
    F_normal::Float64,
    E_star::Float64
)
    # 接触半径
    a = (3.0 * F_normal * R / (4.0 * E_star))^(1.0/3.0)
    
    # 最大压力
    p_max = (6.0 * F_normal * E_star^2 / (π^3 * R^2))^(1.0/3.0)
    
    # 压力分布函数
    function p_hertz(r)
        if r <= a
            return p_max * sqrt(1.0 - (r/a)^2)
        else
            return 0.0
        end
    end
    
    return a, p_max, p_hertz
end

"""
    effective_modulus(E1, ν1, E2, ν2)

计算有效模量（Hertz接触理论）

# 理论
1/E* = (1-ν₁²)/E₁ + (1-ν₂²)/E₂
"""
function effective_modulus(E1, ν1, E2, ν2)
    return 1.0 / ((1.0 - ν1^2)/E1 + (1.0 - ν2^2)/E2)
end

# ================================================================================
# 多物理场耦合接口函数
# ================================================================================

"""
    cohesive_output(case::Case, variables::Dict)

内聚力模型输出和更新（集成到主求解器）
"""
function cohesive_output(case, variables::Dict)
    if !case.opt.cohesive_enabled
        return variables
    end
    
    # 获取界面位移
    if haskey(variables, "negative particle surface displacement") &&
       haskey(variables, "positive particle surface displacement")
        
        disp_n = variables["negative particle surface displacement"]
        disp_p = variables["positive particle surface displacement"]
        
        # 计算相对分离
        δ_n, δ_t = compute_interface_separation(case, disp_n, disp_p)
        
        # 计算内聚力
        coh = case.cohesive_interface
        T_n, T_t, D = compute_cohesive_traction(coh, δ_n, δ_t)
        
        # 更新损伤状态
        update_cohesive_damage!(coh, δ_n, δ_t)
        
        # 存储结果
        variables["cohesive traction normal"] = T_n
        variables["cohesive traction tangential"] = T_t
        variables["cohesive damage"] = D
        variables["cohesive separation normal"] = δ_n
        variables["cohesive separation tangential"] = δ_t
        variables["cohesive energy dissipation"] = sum(coh.energy_dissipated)
        variables["cohesive failed fraction"] = 100.0 * count(coh.is_failed) / length(coh.is_failed)
    end
    
    return variables
end

"""
    contact_output(case::Case, variables::Dict)

接触模型输出和更新（集成到主求解器）
"""
function contact_output(case, variables::Dict)
    if !case.opt.contact_enabled
        return variables
    end
    
    # 获取节点位移
    u_total = assemble_total_displacement(case, variables)
    
    # 简化：假设主从节点位移相同（需根据实际网格修改）
    contact = case.contact_interface
    n_pairs = length(contact.contact_pairs)
    
    # 简化的接触检测
    gap = zeros(Float64, n_pairs)
    for i in 1:n_pairs
        # 基于颗粒表面位移估算间隙
        if haskey(variables, "negative particle surface displacement")
            disp_avg = mean(variables["negative particle surface displacement"])
            gap[i] = contact.gap[i] - disp_avg
        end
    end
    
    contact.gap = gap
    contact.is_contact = gap .< 0.0
    
    # 计算接触压力
    if any(contact.is_contact)
        p_n = compute_contact_pressure(contact)
        contact.p_contact = p_n
        
        # 计算摩擦力（需要速度信息）
        v_t = zeros(Float64, n_pairs)  # 简化：零速度
        τ_t, is_sliding = compute_friction_stress(contact, v_t)
        contact.τ_friction = τ_t
        contact.is_sliding = is_sliding
        
        # 存储结果
        variables["contact gap"] = gap
        variables["contact pressure"] = p_n
        variables["friction stress"] = τ_t
        variables["contact area"] = sum(contact.contact_area[contact.is_contact])
        variables["friction dissipation"] = 0.0  # 需要时间步长信息
    end
    
    return variables
end

"""
    compute_interface_separation(case::Case, disp_negative, disp_positive)

计算界面分离（法向和切向）
"""
function compute_interface_separation(
    case,
    disp_negative::Vector{Float64},
    disp_positive::Vector{Float64}
)
    coh = case.cohesive_interface
    n_pairs = length(coh.element_pairs)
    
    δ_n = zeros(Float64, n_pairs)
    δ_t = zeros(Float64, n_pairs)
    
    @inbounds for i in 1:n_pairs
        (neg_id, pos_id) = coh.element_pairs[i]
        
        # 相对位移
        if neg_id <= length(disp_negative) && pos_id <= length(disp_positive)
            Δu = disp_positive[pos_id] - disp_negative[neg_id]
            
            # 投影到法向（简化：假设法向为正）
            δ_n[i] = Δu
            δ_t[i] = 0.0  # 简化：忽略切向
        end
    end
    
    return δ_n, δ_t
end

"""
    assemble_total_displacement(case::Case, variables::Dict)

组装总位移（扩散位移 + 热位移）
"""
function assemble_total_displacement(case, variables::Dict)
    # 简化实现：返回颗粒表面位移
    if haskey(variables, "negative particle surface displacement")
        return variables["negative particle surface displacement"]
    else
        return Float64[]
    end
end

# ================================================================================
# 导出函数
# ================================================================================

# 导出到父模块（在 JuBat.jl 中添加）
# export CohesiveInterface, ContactInterface
# export init_cohesive_interface, compute_cohesive_traction, update_cohesive_damage!
# export init_contact_interface, detect_contact, compute_contact_pressure
# export cohesive_output, contact_output, hertz_contact_pressure
