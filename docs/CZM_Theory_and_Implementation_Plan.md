# 内聚力模型(CZM)在极片脱粘中的应用

## 一、背景与动机

### 1.1 极片脱粘问题

在锂离子电池循环过程中，由于化学应变、热应变和机械应力的累积作用，极片可能发生多种界面失效：

```
循环充放电 → 体积变化 → 界面应力 → 裂纹萌生 → 裂纹扩展 → 脱粘失效
     ↓            ↓           ↓           ↓           ↓           ↓
  电化学      化学应变    应力集中    损伤累积    能量释放    容量衰减
```

#### 关键界面
1. **活性材料-粘结剂界面**
   - 颗粒膨胀/收缩
   - 粘结剂断裂
   - 电接触损失

2. **涂层-集流体界面**
   - 剪切应力
   - 法向拉应力
   - 分层脱粘

3. **颗粒-颗粒界面**
   - 接触损失
   - 裂纹扩展
   - 网络断裂

### 1.2 为什么需要CZM？

传统线弹性断裂力学(LEFM)的局限：
- ❌ 假设裂纹尖端应力奇异
- ❌ 需要预先存在裂纹
- ❌ 不能描述渐进损伤过程

内聚力模型(CZM)的优势：
- ✅ 描述渐进失效过程
- ✅ 无需预置裂纹
- ✅ 自然考虑损伤累积
- ✅ 能量平衡明确

---

## 二、CZM理论基础

### 2.1 核心思想

**界面用"内聚力-位移关系"描述，而非完美粘接**

```
完美界面（传统）:
  u⁺ = u⁻  (位移连续)
  t = ∞    (可承受任意牵引力)

内聚界面（CZM）:
  δ = u⁺ - u⁻  (允许位移跳跃)
  t = t(δ, D)   (牵引力-分离关系)
  D ∈ [0,1]     (损伤变量)
```

### 2.2 牵引-分离律(Traction-Separation Law)

#### 基本形式
```
牵引力：t = t_n·n + t_t·t  (法向+切向)
分离量：δ = δ_n·n + δ_t·t

关系：t = f(δ, D)
```

#### 典型本构模型

**1. 双线性模型(Bilinear)**
```
阶段1：弹性加载 (D=0)
  t = K·δ,  δ < δ_0

阶段2：损伤软化 (0<D<1)
  t = (1-D)·K·δ,  δ_0 ≤ δ < δ_f

阶段3：完全失效 (D=1)
  t = 0,  δ ≥ δ_f
```

**2. 指数软化模型(Exponential)**
```
t = t_max·(δ/δ_0)·exp[1 - δ/δ_0]
```

**3. Park-Paulino-Roesler (PPR)模型**
```
t = t_max·[(δ/δ_f)·exp(-δ/δ_f)]^α
```

### 2.3 关键参数

#### 材料参数
1. **界面刚度** `K_n`, `K_t` [Pa/m]
   - 法向刚度（抗拉）
   - 切向刚度（抗剪）

2. **临界牵引力** `t_n^max`, `t_t^max` [Pa]
   - 法向强度
   - 切向强度

3. **临界能量释放率** `G_Ic`, `G_IIc` [J/m²]
   - 模式I（张开）
   - 模式II（剪切）

4. **临界分离量** `δ_n^c`, `δ_t^c` [m]
   ```
   δ_c = 2·G_c / t_max
   ```

#### 损伤演化参数
1. **损伤起始准则**
   ```
   最大应力：max(t_n/t_n^max, t_t/t_t^max) = 1
   二次应力：(t_n/t_n^max)² + (t_t/t_t^max)² = 1
   ```

2. **损伤演化律**
   ```
   线性：D = (δ - δ_0) / (δ_f - δ_0)
   指数：D = 1 - exp[-α·(δ - δ_0)]
   幂律：D = [(δ - δ_0) / (δ_f - δ_0)]^β
   ```

### 2.4 能量平衡

**断裂能（内聚能）**：
```
G = ∫₀^δ_f t(δ) dδ
```

对于双线性模型：
```
G = (1/2)·t_max·δ_f
```

**能量释放条件**：
```
G ≥ G_c  → 裂纹扩展
```

---

## 三、多物理场耦合

### 3.1 电化学-力学-CZM耦合

```
电化学模型 (SPMe/P2D)
    ↓
锂浓度 cs(r,x,t), SOC(x,t)
    ↓
化学应变 ε_chem = eps_s × Ω × cs_av
    ↓
热应变 ε_thermal = α × ΔT
    ↓
总应变 ε_total = ε_elastic + ε_chem + ε_thermal
    ↓
应力场 σ(x,y,t)  [来自mechanical.jl]
    ↓
界面牵引力 t_n, t_t  [投影到界面]
    ↓
CZM损伤演化 D(t)
    ↓
界面刚度退化 K_eff = (1-D)·K
    ↓
反馈到力学求解
```

### 3.2 关键变量传递

#### 从力学模块获取
```julia
# 宏观应力场（2D）
σ_xx, σ_yy, σ_xy  [Pa]

# 位移场
u_x, u_y  [m]

# 应变能密度
U_strain  [J/m³]

# Von Mises应力
σ_vm  [Pa]

# 界面法向应力（需投影）
σ_n = n·σ·n

# 界面剪切应力
τ = t·σ·n
```

#### 输出到系统
```julia
# 损伤变量
D  [-], 0 ≤ D ≤ 1

# 有效界面刚度
K_eff = (1-D)·K  [Pa/m]

# 界面完整性
integrity = 1 - D

# 累积损伤能
G_dissipated  [J/m²]

# 失效位置
failed_interfaces  [boolean array]
```

---

## 四、技术路线

### 4.1 实施层次

#### Level 1：界面定义与参数（基础）
- [x] 定义界面类型（涂层-集流体，颗粒-粘结剂）
- [x] 输入材料参数（K, t_max, G_c）
- [x] 选择本构模型（双线性/指数/PPR）

#### Level 2：界面单元与几何（核心）
- [x] 标识界面单元/节点
- [x] 计算界面法向和切向
- [x] 提取界面位移跳跃

#### Level 3：CZM本构计算（算法）
- [x] 计算当前牵引力 t(δ, D)
- [x] 检查损伤起始准则
- [x] 更新损伤变量 D
- [x] 计算有效刚度 K_eff

#### Level 4：与力学模块耦合（集成）
- [x] 从力学模块读取应力/位移
- [x] 投影到界面坐标系
- [x] 反馈界面刚度退化
- [x] 更新边界条件

#### Level 5：时间演化与后处理（分析）
- [x] 损伤历史追踪
- [x] 界面失效判定
- [x] 能量平衡验证
- [x] 可视化

### 4.2 实施步骤（10步）

#### 步骤1：定义数据结构
```julia
# 文件：src/czm.jl

"""内聚界面类型"""
@enum InterfaceType begin
    COATING_COLLECTOR    # 涂层-集流体
    PARTICLE_BINDER      # 颗粒-粘结剂
    PARTICLE_PARTICLE    # 颗粒-颗粒
end

"""CZM本构模型类型"""
@enum CZMModel begin
    BILINEAR      # 双线性
    EXPONENTIAL   # 指数软化
    PPR           # Park-Paulino-Roesler
end

"""CZM材料参数"""
mutable struct CZMMaterial
    # 界面刚度
    K_n::Float64          # 法向刚度 [Pa/m]
    K_t::Float64          # 切向刚度 [Pa/m]
    
    # 临界牵引力
    t_n_max::Float64      # 法向强度 [Pa]
    t_t_max::Float64      # 切向强度 [Pa]
    
    # 临界能量释放率
    G_Ic::Float64         # 模式I [J/m²]
    G_IIc::Float64        # 模式II [J/m²]
    
    # 损伤演化参数
    alpha::Float64        # 损伤指数 [-]
    beta::Float64         # 幂律指数 [-]
    
    # 本构模型选择
    model::CZMModel
end

"""CZM界面元素"""
mutable struct CZMInterface
    # 几何信息
    id::Int64
    type::InterfaceType
    node_plus::Int64      # "+"侧节点
    node_minus::Int64     # "-"侧节点
    normal::Vector{Float64}    # 法向向量 [nx, ny]
    tangent::Vector{Float64}   # 切向向量 [tx, ty]
    area::Float64         # 界面面积 [m²]
    
    # 材料参数
    material::CZMMaterial
    
    # 状态变量
    D::Float64            # 损伤变量 [-]
    delta_n::Float64      # 法向分离 [m]
    delta_t::Float64      # 切向分离 [m]
    t_n::Float64          # 法向牵引力 [Pa]
    t_t::Float64          # 切向牵引力 [Pa]
    
    # 历史变量
    delta_max::Float64    # 最大历史分离
    G_dissipated::Float64 # 累积耗散能 [J/m²]
    failed::Bool          # 失效标志
end
```

#### 步骤2：界面识别算法
```julia
"""
识别涂层-集流体界面

基于几何：
- 在thermal2D网格中
- 单元属于涂层（NE或PE）
- 且贴近集流体边界
"""
function identify_coating_collector_interfaces(mesh, param_dim)
    interfaces = CZMInterface[]
    
    ne = size(mesh.element, 1)
    pgeo = jellyroll_spiral_params(param_dim)
    
    # 计算单元中心
    centers = jellyroll_element_centers(mesh)
    
    for e in 1:ne
        x, y = centers[e, 1], centers[e, 2]
        r = sqrt(x^2 + y^2)
        θ = atan(y, x)
        
        # 判断材料类型
        layer = material_at(r, θ, param_dim; logic=:spiral)
        
        # 检查是否接近集流体边界
        if layer == :NE
            # 负极涂层，检查是否接近NCC内侧
            boundary = check_near_NCC_boundary(r, θ, pgeo)
            if boundary
                # 创建界面单元
                interface = create_interface(e, :COATING_COLLECTOR, :NE, mesh)
                push!(interfaces, interface)
            end
        elseif layer == :PE
            # 正极涂层，检查是否接近PCC外侧
            boundary = check_near_PCC_boundary(r, θ, pgeo)
            if boundary
                interface = create_interface(e, :COATING_COLLECTOR, :PE, mesh)
                push!(interfaces, interface)
            end
        end
    end
    
    return interfaces
end
```

#### 步骤3：提取界面位移跳跃
```julia
"""
从位移场计算界面分离量

输入：
- interface: CZM界面对象
- u_global: 全局位移向量 [u_x, u_y] for all nodes

输出：
- delta_n: 法向分离
- delta_t: 切向分离
"""
function compute_interface_separation(interface::CZMInterface, u_global::Vector{Float64})
    # 提取两侧节点位移
    node_plus = interface.node_plus
    node_minus = interface.node_minus
    
    u_plus = [u_global[2*node_plus-1], u_global[2*node_plus]]
    u_minus = [u_global[2*node_minus-1], u_global[2*node_minus]]
    
    # 位移跳跃（全局坐标系）
    Δu = u_plus - u_minus
    
    # 投影到局部坐标系
    n = interface.normal
    t = interface.tangent
    
    delta_n = dot(Δu, n)  # 法向分离
    delta_t = dot(Δu, t)  # 切向分离
    
    return delta_n, delta_t
end
```

#### 步骤4：CZM本构计算
```julia
"""
双线性内聚力模型

阶段1：弹性 (δ < δ_0)
  t = K·δ
  D = 0

阶段2：软化 (δ_0 ≤ δ < δ_f)
  t = t_max·(δ_f - δ)/(δ_f - δ_0)
  D = (δ - δ_0)/(δ_f - δ_0)

阶段3：失效 (δ ≥ δ_f)
  t = 0
  D = 1
"""
function bilinear_cohesive_law(delta_n, delta_t, material::CZMMaterial, D_old::Float64)
    # 等效分离（混合模式）
    delta_eff = sqrt(delta_n^2 + delta_t^2)
    
    # 临界分离量
    delta_0 = material.t_n_max / material.K_n  # 损伤起始
    delta_f = 2 * material.G_Ic / material.t_n_max  # 完全失效
    
    # 计算损伤变量
    if delta_eff < delta_0
        # 阶段1：弹性
        D = 0.0
    elseif delta_eff < delta_f
        # 阶段2：软化
        D = (delta_eff - delta_0) / (delta_f - delta_0)
        D = max(D, D_old)  # 单调增长（不可恢复）
    else
        # 阶段3：完全失效
        D = 1.0
    end
    
    # 计算牵引力
    if D < 1.0
        K_eff = (1.0 - D) * material.K_n
        t_n = K_eff * delta_n
        t_t = (1.0 - D) * material.K_t * delta_t
    else
        t_n = 0.0
        t_t = 0.0
    end
    
    # 能量耗散
    if D > D_old
        dG = 0.5 * (t_n * delta_n + t_t * delta_t) * (D - D_old) / D
    else
        dG = 0.0
    end
    
    return t_n, t_t, D, dG
end
```

#### 步骤5：损伤演化算法
```julia
"""
更新所有界面的CZM状态

时间推进：t^n → t^{n+1}
"""
function update_czm_state!(interfaces::Vector{CZMInterface}, 
                          u_global::Vector{Float64}, 
                          dt::Float64)
    for interface in interfaces
        # 计算当前分离量
        delta_n, delta_t = compute_interface_separation(interface, u_global)
        
        # 保存旧状态
        D_old = interface.D
        
        # 调用本构模型
        if interface.material.model == BILINEAR
            t_n, t_t, D, dG = bilinear_cohesive_law(
                delta_n, delta_t, interface.material, D_old
            )
        elseif interface.material.model == EXPONENTIAL
            t_n, t_t, D, dG = exponential_cohesive_law(...)
        else
            error("Unknown CZM model")
        end
        
        # 更新状态
        interface.delta_n = delta_n
        interface.delta_t = delta_t
        interface.t_n = t_n
        interface.t_t = t_t
        interface.D = D
        interface.G_dissipated += dG
        
        # 更新最大历史分离
        delta_eff = sqrt(delta_n^2 + delta_t^2)
        interface.delta_max = max(interface.delta_max, delta_eff)
        
        # 失效判定
        if D ≥ 0.99
            interface.failed = true
        end
    end
end
```

#### 步骤6：与力学模块耦合
```julia
"""
在力学求解中考虑CZM界面

修改方式：
1. 界面单元的刚度矩阵添加内聚力贡献
2. 或在边界条件中施加内聚牵引力
"""
function apply_czm_to_mechanics!(K_mech, F_mech, interfaces, mesh)
    for interface in interfaces
        if interface.failed
            continue  # 已失效界面不贡献
        end
        
        # 获取界面刚度
        K_n = (1.0 - interface.D) * interface.material.K_n
        K_t = (1.0 - interface.D) * interface.material.K_t
        
        # 构造局部刚度矩阵（2节点，4自由度）
        # K_local = [K_n*n⊗n + K_t*t⊗t]
        
        # 组装到全局刚度矩阵
        # K_mech[dofs, dofs] += K_local
        
        # 或者：施加等效节点力
        # F_mech[dofs] += t_n*n + t_t*t
    end
end
```

#### 步骤7：输入参数接口
```julia
"""
从参数文件读取CZM材料参数

添加到 parameters/Jellyroll.jl 或单独文件
"""
function load_czm_parameters(param_dim)
    # 涂层-集流体界面（负极）
    czm_NE_NCC = CZMMaterial(
        K_n = 1e12,        # Pa/m (估算：E/h_interface)
        K_t = 5e11,
        t_n_max = 10e6,    # Pa (粘结强度)
        t_t_max = 8e6,
        G_Ic = 100.0,      # J/m² (断裂能)
        G_IIc = 150.0,
        alpha = 1.0,
        beta = 1.0,
        model = BILINEAR
    )
    
    # 涂层-集流体界面（正极）
    czm_PE_PCC = CZMMaterial(
        K_n = 8e11,
        K_t = 4e11,
        t_n_max = 8e6,
        t_t_max = 6e6,
        G_Ic = 80.0,
        G_IIc = 120.0,
        alpha = 1.0,
        beta = 1.0,
        model = BILINEAR
    )
    
    return Dict(
        :NE_NCC => czm_NE_NCC,
        :PE_PCC => czm_PE_PCC
    )
end
```

#### 步骤8：时间积分策略
```julia
"""
CZM时间推进算法

隐式（稳定但需迭代）：
  t^{n+1} = f(δ^{n+1}, D^{n+1})
  D^{n+1} = g(δ^{n+1}, D^n)
  → 需要在力学迭代中同步更新

显式（简单但可能不稳定）：
  D^{n+1} = g(δ^n, D^n)
  t^{n+1} = f(δ^{n+1}, D^{n+1})
  → 先更新损伤，再计算牵引力
"""
function time_integration_czm!(case, variables, dt)
    # 获取界面列表
    interfaces = case.czm_interfaces
    
    # 获取当前位移场
    u_global = extract_displacement(variables)
    
    # 方法1：显式更新（简单，推荐初期使用）
    update_czm_state!(interfaces, u_global, dt)
    
    # 方法2：隐式更新（需要与力学迭代耦合）
    # 在mechanical.jl的thermal_diffusion_stress_2D中嵌入
end
```

#### 步骤9：后处理与可视化
```julia
"""
输出CZM结果

添加到 result 字典
"""
function postprocess_czm(interfaces, result)
    n_interfaces = length(interfaces)
    
    # 初始化数组
    damage_array = zeros(Float64, n_interfaces)
    t_n_array = zeros(Float64, n_interfaces)
    t_t_array = zeros(Float64, n_interfaces)
    failed_array = zeros(Bool, n_interfaces)
    
    for (i, interface) in enumerate(interfaces)
        damage_array[i] = interface.D
        t_n_array[i] = interface.t_n
        t_t_array[i] = interface.t_t
        failed_array[i] = interface.failed
    end
    
    # 存储到result
    result["czm damage"] = damage_array
    result["czm traction normal"] = t_n_array
    result["czm traction tangent"] = t_t_array
    result["czm failed interfaces"] = failed_array
    result["czm number failed"] = sum(failed_array)
    
    # 统计信息
    println("CZM Summary:")
    println("  Total interfaces: $(n_interfaces)")
    println("  Failed interfaces: $(sum(failed_array))")
    println("  Max damage: $(maximum(damage_array))")
    println("  Avg damage: $(mean(damage_array))")
end
```

#### 步骤10：验证案例
```julia
"""
CZM验证案例

双悬臂梁(DCB)测试
"""
function test_czm_dcb()
    # 创建简单几何
    # 施加上下拉伸载荷
    # 观察裂纹扩展
    # 对比理论解
    
    # 理论断裂能：G = ∫ t(δ) dδ
    # 数值断裂能：从CZM计算的G_dissipated
    
    @test abs(G_theory - G_numerical) / G_theory < 0.05
end
```

---

## 五、需要的初始参数

### 5.1 CZM材料参数（必需）

#### 涂层-集流体界面
```julia
# 负极：NE-NCC界面
CZM_NE_NCC = Dict(
    # 界面刚度（估算：E/h_interface）
    :K_n => 1e12,         # Pa/m
    :K_t => 5e11,         # Pa/m
    
    # 临界牵引力（来自粘结强度测试）
    :t_n_max => 10e6,     # Pa (10 MPa)
    :t_t_max => 8e6,      # Pa (8 MPa)
    
    # 临界能量释放率（断裂韧性）
    :G_Ic => 100.0,       # J/m²
    :G_IIc => 150.0,      # J/m²
    
    # 本构模型
    :model => :bilinear
)

# 正极：PE-PCC界面
CZM_PE_PCC = Dict(
    :K_n => 8e11,
    :K_t => 4e11,
    :t_n_max => 8e6,
    :t_t_max => 6e6,
    :G_Ic => 80.0,
    :G_IIc => 120.0,
    :model => :bilinear
)
```

#### 典型参数范围（文献值）
| 参数 | 范围 | 来源 |
|------|------|------|
| K_n | 1e11 ~ 1e13 Pa/m | 估算：E/h |
| t_n_max | 1 ~ 50 MPa | 剥离测试 |
| G_Ic | 10 ~ 500 J/m² | DCB测试 |
| G_IIc | 50 ~ 1000 J/m² | ENF测试 |

### 5.2 界面几何参数

```julia
# 界面识别参数
Interface_Geometry = Dict(
    # 界面类型
    :type => :coating_collector,  # or :particle_binder
    
    # 空间定位
    :location => :inner_spiral,    # NCC内侧
    # or :outer_spiral  # PCC外侧
    
    # 几何容差
    :tolerance => 1e-6,  # m
    
    # 是否考虑分层
    :use_layers => true
)
```

### 5.3 数值参数

```julia
# CZM求解参数
CZM_Numerics = Dict(
    # 时间积分
    :integration_scheme => :explicit,  # or :implicit
    
    # 损伤容差
    :damage_tolerance => 1e-6,
    
    # 失效判据
    :failure_criterion => 0.99,  # D ≥ 0.99 → 失效
    
    # 正则化参数（防止数值振荡）
    :viscous_regularization => 1e-5,
    
    # 最大迭代次数（隐式）
    :max_iterations => 50,
    
    # 收敛容差
    :convergence_tol => 1e-8
)
```

### 5.4 控制开关

```julia
# CZM模块开关
opt.czm_enabled = true
opt.czm_interface_type = "coating_collector"  # or "particle_binder"
opt.czm_model = "bilinear"                    # or "exponential", "PPR"
opt.czm_output_frequency = 10                 # 每10步输出一次
```

---

## 六、需要从力学模块获取的变量

### 6.1 必需变量（Level 1）

#### 位移场
```julia
# 全局位移向量
u_x::Vector{Float64}  # x方向位移 [m]
u_y::Vector{Float64}  # y方向位移 [m]

# 或合并为
U_global::Vector{Float64}  # [u1_x, u1_y, u2_x, u2_y, ...]

# 来源：thermal_diffusion_stress_2D → U_M
```

**用途**：
- 计算界面位移跳跃
- 投影到法向/切向分量

#### 应力场
```julia
# 单元应力
σ_xx::Vector{Float64}  # [Pa]
σ_yy::Vector{Float64}  # [Pa]
σ_xy::Vector{Float64}  # [Pa]

# 来源：_recover_stress_2D
```

**用途**：
- 投影到界面法向应力
- 验证损伤准则

### 6.2 推荐变量（Level 2）

#### 应变能
```julia
# 单元应变能密度
U_strain::Vector{Float64}  # [J/m³]

# 计算方式
U_strain[e] = 0.5 * (σ_xx[e]*ε_xx + σ_yy[e]*ε_yy + σ_xy[e]*γ_xy)
```

**用途**：
- 能量平衡验证
- 损伤驱动力评估

#### Von Mises应力
```julia
σ_vm::Vector{Float64}  # [Pa]

# 来源：已有输出
```

**用途**：
- 损伤萌生判据
- 应力集中识别

### 6.3 可选变量（Level 3）

#### 应力分量历史
```julia
# 时间历程
σ_xx_history::Matrix{Float64}  # [ne × num_steps]
σ_yy_history::Matrix{Float64}
σ_xy_history::Matrix{Float64}
```

**用途**：
- 疲劳损伤累积
- 循环加载分析

#### 塑性应变（如有）
```julia
ε_plastic::Vector{Float64}  # [-]
```

**用途**：
- 考虑塑性耗散
- 改进损伤模型

### 6.4 接口设计

#### 建议在 `mechanical.jl` 中添加
```julia
function extract_variables_for_czm(variables, case)
    """
    从variables字典提取CZM所需变量
    
    返回：
    - u_global: 全局位移向量
    - sigma: 应力场 (σ_xx, σ_yy, σ_xy)
    - strain_energy: 应变能密度
    """
    
    # 提取位移
    u_x = variables["displacement x"]
    u_y = variables["displacement y"]
    u_global = interleave_xy(u_x, u_y)
    
    # 提取应力
    sigma = Dict(
        "xx" => variables["diffusion stress xx"],
        "yy" => variables["diffusion stress yy"],
        "xy" => variables["diffusion stress xy"],
        "vm" => variables["diffusion stress vonMises"]
    )
    
    # 计算应变能（如果未存储）
    if !haskey(variables, "strain energy density")
        strain_energy = compute_strain_energy(sigma, case)
        variables["strain energy density"] = strain_energy
    else
        strain_energy = variables["strain energy density"]
    end
    
    return u_global, sigma, strain_energy
end
```

---

## 七、实施优先级

### Phase 1：基础框架（1周）
**目标**：建立CZM基本数据结构和接口

- [ ] 创建 `src/czm.jl`
- [ ] 定义数据结构（CZMMaterial, CZMInterface）
- [ ] 实现双线性本构模型
- [ ] 编写参数加载函数
- [ ] 单元测试（简单拉伸）

**交付物**：
- `src/czm.jl` (500行)
- `test/test_czm_basic.jl`
- 文档：接口说明

### Phase 2：界面识别（3天）
**目标**：在现有网格中识别CZM界面

- [ ] 实现涂层-集流体界面识别
- [ ] 计算界面法向和切向
- [ ] 可视化界面位置
- [ ] 验证几何正确性

**交付物**：
- 界面识别函数
- 可视化脚本
- 验证报告

### Phase 3：力学耦合（1周）
**目标**：将CZM集成到力学求解

- [ ] 从mechanical.jl提取变量
- [ ] 计算界面分离量
- [ ] 更新损伤状态
- [ ] 反馈到力学边界条件

**交付物**：
- 耦合函数
- 示例案例
- 对比验证

### Phase 4：时间演化（3天）
**目标**：实现多步时间推进

- [ ] 实现显式时间积分
- [ ] 损伤历史追踪
- [ ] 失效判定
- [ ] 能量平衡检查

**交付物**：
- 时间积分算法
- 后处理工具
- 能量验证

### Phase 5：验证与应用（1周）
**目标**：验证模型并应用到实际问题

- [ ] DCB测试（对比理论）
- [ ] 循环加载测试
- [ ] Jellyroll脱粘仿真
- [ ] 参数敏感性分析

**交付物**：
- 验证案例集
- 应用示例
- 技术报告

---

## 八、参数获取来源

### 8.1 实验测试

#### 粘结强度测试
- **90度剥离测试**：测量 t_n_max
- **剪切强度测试**：测量 t_t_max
- 设备：万能试验机 + 剥离夹具

#### 断裂韧性测试
- **双悬臂梁(DCB)**：测量 G_Ic
- **端部切口弯曲(ENF)**：测量 G_IIc
- 数据处理：载荷-位移曲线积分

### 8.2 文献数据

#### 涂层-集流体界面
- Rahani & Shenoy (2013): G_Ic ≈ 50-200 J/m²
- Zhang et al. (2017): t_max ≈ 5-15 MPa
- Xu et al. (2016): K_n ≈ 1e12 Pa/m

#### 颗粒-粘结剂界面
- Wu & Xiao (2017): G_Ic ≈ 10-50 J/m²
- Zhao et al. (2018): t_max ≈ 1-5 MPa

### 8.3 估算方法

#### 界面刚度
```
K_n ≈ E_interface / h_interface
```
其中：
- E_interface：界面等效模量（几何平均：√(E1·E2)）
- h_interface：界面厚度（通常取 0.1-1 μm）

#### 临界牵引力
```
t_max ≈ σ_yield_interface
```
使用粘结剂的屈服强度作为上限

#### 断裂能
```
G_c ≈ K_Ic² / E
```
从应力强度因子转换

---

## 九、预期挑战与解决方案

### 挑战1：界面单元定义
**问题**：2D网格中没有显式的界面单元

**解决方案**：
- 基于节点对定义"虚拟界面"
- 或插入零厚度内聚单元
- 使用penalty方法施加界面约束

### 挑战2：多时间尺度
**问题**：损伤演化可能比电化学快/慢

**解决方案**：
- 子循环：在一个电化学步内多次更新CZM
- 或仅在关键时刻激活CZM
- 自适应时间步长

### 挑战3：参数不确定性
**问题**：CZM参数难以准确测量

**解决方案**：
- 参数识别：反向拟合实验曲线
- 敏感性分析：识别关键参数
- 不确定性量化：区间/概率分析

### 挑战4：计算成本
**问题**：界面单元多，计算慢

**解决方案**：
- 仅在关键区域激活CZM
- 显式积分（避免迭代）
- 并行计算界面更新

---

## 十、总结

### 核心要素
1. **理论**：内聚力模型描述渐进损伤
2. **几何**：识别界面位置和方向
3. **本构**：牵引-分离律（双线性/指数）
4. **耦合**：从力学获取应力/位移
5. **演化**：时间推进损伤变量

### 关键参数
| 类别 | 参数 | 来源 |
|------|------|------|
| 刚度 | K_n, K_t | 估算或DIC |
| 强度 | t_max | 剥离测试 |
| 韧性 | G_c | DCB/ENF |
| 模型 | 双线性/指数 | 根据材料选择 |

### 需要的变量
| 变量 | 来源 | 用途 |
|------|------|------|
| u_x, u_y | thermal_diffusion_stress_2D | 计算分离 |
| σ_xx, σ_yy, σ_xy | _recover_stress_2D | 投影牵引力 |
| U_strain | 计算 | 能量平衡 |

### 实施路径
```
Phase 1: 数据结构 → Phase 2: 界面识别 → 
Phase 3: 力学耦合 → Phase 4: 时间演化 → 
Phase 5: 验证应用
```

---

**文档版本**: v1.0  
**创建日期**: 2025-12-29  
**状态**: 技术路线已完成，待实施  
**预计工时**: 3-4周
