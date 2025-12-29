# 基于颗粒体积膨胀的化学应变计算新路线

## 一、新路线核心思想

### 1.1 物理本质
**化学应变的根源是颗粒的真实体积膨胀，而非SOC的变化**

```
SOC变化 → 颗粒内浓度分布变化 → 颗粒体积膨胀 → 宏观应变
   ↑                                      ↓
   └──────── 不应跳过中间过程 ──────────┘
```

### 1.2 当前方法的问题
**旧路线**（SOC代理法）：
```julia
# 直接用SOC变化计算宏观应变
ε_macro = β * ΔSOC = (Ω/3) * eps_s * ΔSOC
```

**问题**：
1. ❌ 假设颗粒内浓度均匀（实际有梯度）
2. ❌ 忽略颗粒尺度的力学响应
3. ❌ SOC仅是平均值，丢失空间分布信息
4. ❌ 未利用已有的颗粒力学计算（`Calstressdisp`）

### 1.3 新路线优势
**新路线**（颗粒位移法）：
```julia
# 步骤1：颗粒力学计算 → 表面位移
disp_surf = Calstressdisp(cs(r))  # 从浓度分布计算位移

# 步骤2：颗粒体积膨胀
ε_particle = 3 * disp_surf / rs  # 球形颗粒体积应变

# 步骤3：宏观应变
ε_macro = eps_s * ε_particle  # 均质化
```

**优势**：
1. ✅ 考虑颗粒内真实浓度梯度
2. ✅ 利用已有的力学计算
3. ✅ 物理路径完整（无跳跃）
4. ✅ 自然考虑应力-扩散耦合

---

## 二、理论基础

### 2.1 颗粒体积膨胀计算

#### 球形颗粒的体积应变
对于半径为 `rs` 的球形颗粒，表面径向位移为 `disp_surf` 时：

```
新半径: rs_new = rs + disp_surf
新体积: V_new = (4/3)π(rs + disp_surf)³
原体积: V_old = (4/3)πrs³

体积应变:
ΔV/V = (V_new - V_old) / V_old
     = [(rs + disp_surf)³ - rs³] / rs³
     = [1 + disp_surf/rs]³ - 1
     ≈ 3 * disp_surf/rs  (小变形假设)
```

**关键公式**：
```julia
ε_particle_volumetric = 3 * disp_surf / rs
```

#### 与代码中的disp_surf关系
在 `mechanical.jl::Calstressdisp` 中（第136行）：
```julia
disp_surf = (Omega * rs * cs_av) / 3
```

其中 `cs_av` 是颗粒平均浓度（第133行）：
```julia
cs_av = (3/(4*π*rs³)) * ∫ cs(r) * 4π*r² dr
```

代入得：
```julia
ε_particle = 3 * disp_surf / rs
          = 3 * (Omega * rs * cs_av / 3) / rs
          = Omega * cs_av
```

**这与经典理论一致！**（Christensen & Newman 2006）

### 2.2 宏观均质化

#### 单元总体积膨胀
假设单元体积为 `V_elem`，其中：
- 固相体积：`V_solid = eps_s * V_elem`（包含所有颗粒）
- 液相体积：`V_liquid = eps_e * V_elem`（电解液，不膨胀）
- 粘结剂体积：`V_binder = eps_fi * V_elem`（假设不膨胀）

**关键假设**：
> 单元所有的体积膨胀都来自于颗粒的体积膨胀

则：
```
ΔV_elem / V_elem = ΔV_solid / V_elem
                 = (V_solid * ε_particle) / V_elem
                 = eps_s * ε_particle
```

**这正是均质化理论！**

#### 宏观应变张量
对于各向同性体积膨胀：
```
ε_macro = (ΔV/V) / 3 * I
        = (eps_s * ε_particle) / 3 * I
```

但在2D平面应力问题中，厚度方向可以自由膨胀，因此：
```
ε_macro_2D = (eps_s * ε_particle) * [1, 1, 0]^T
```

### 2.3 与SOC方法的对比

#### SOC方法（旧）
```
ε_particle = Ω * Δc = Ω * c_max * ΔSOC
ε_macro = eps_s * Ω * c_max * ΔSOC
```
**假设**：颗粒内浓度均匀变化

#### 位移方法（新）
```
ε_particle = 3 * disp_surf / rs
           = 3 * (Ω * rs * cs_av) / (3 * rs)
           = Ω * cs_av
ε_macro = eps_s * Ω * cs_av
```
**优势**：`cs_av` 由颗粒内真实浓度分布积分得到，考虑梯度

#### 差异分析
当颗粒内浓度均匀时（`cs(r) = const`）：
```
cs_av = cs = c_max * SOC
→ 两种方法等价
```

当颗粒内浓度有梯度时（扩散受限）：
```
cs_av ≠ cs_surf
→ 位移方法更准确（基于真实分布）
```

---

## 三、技术路线：从颗粒到宏观

### 3.1 数据流

```
电化学求解 → 颗粒浓度分布 cs(r,x,t)
                    ↓
颗粒力学计算 Calstressdisp(cs) → disp_surf
                    ↓
颗粒体积应变 ε_particle = 3 * disp_surf / rs
                    ↓
宏观应变 ε_macro = eps_s * ε_particle
                    ↓
2D有限元应力求解 → σ(x,y)
```

### 3.2 实施位置

#### 当前架构（SPMe + Thermal2D）

**颗粒尺度计算**（已有）：
- 文件：`src/mechanical.jl::Calstressdisp`
- 输入：`cs(r)` - 颗粒内浓度分布
- 输出：`disp_surf` - 颗粒表面位移

**宏观应力计算**（需修改）：
- 文件：`src/mechanical.jl::thermal_diffusion_stress_2D`
- 当前输入：`SOC_n`, `SOC_p` - 单元平均SOC
- **新输入**：`disp_surf_n`, `disp_surf_p` - 颗粒表面位移

### 3.3 代码修改方案

#### 方案A：在宏观应力计算前计算disp_surf（推荐）

**位置**：`thermal_diffusion_stress_2D` 函数内部

**步骤**：
1. 获取单元的颗粒浓度数据
2. 调用 `Calstressdisp` 计算 `disp_surf`
3. 计算颗粒体积应变：`ε_particle = 3 * disp_surf / rs`
4. 计算宏观应变：`ε_macro = eps_s * ε_particle`

**伪代码**：
```julia
function thermal_diffusion_stress_2D(case, variables)
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    
    # 初始化体积应变数组
    ε_particle_n = zeros(Float64, ne)
    ε_particle_p = zeros(Float64, ne)
    
    # 对每个单元
    for e in 1:ne
        # 获取该单元的颗粒浓度分布
        cs_n = get_particle_concentration_in_element(e, variables, :negative)
        cs_p = get_particle_concentration_in_element(e, variables, :positive)
        
        # 调用颗粒力学计算
        _, _, disp_surf_n, _, _ = Calstressdisp(param.NE, mesh_particle_n, cs_n, T[e])
        _, _, disp_surf_p, _, _ = Calstressdisp(param.PE, mesh_particle_p, cs_p, T[e])
        
        # 计算颗粒体积应变
        ε_particle_n[e] = 3.0 * disp_surf_n / param.NE.rs
        ε_particle_p[e] = 3.0 * disp_surf_p / param.PE.rs
    end
    
    # 计算宏观初始应变（含体积分数修正）
    for e in 1:ne
        epsilon_0_elem[e] = α_eff * dT_elem[e] + 
                           param.NE.eps_s * ε_particle_n[e] +
                           param.PE.eps_s * ε_particle_p[e]
    end
    
    # 后续的有限元求解保持不变
    ...
end
```

#### 方案B：在电化学求解时同步计算并存储（更高效）

**位置**：`Mechanicaloutput` 或 `ThermalDistributed` 中已有颗粒力学计算的地方

**优势**：
- 避免重复计算 `Calstressdisp`
- 利用已有的数据流
- 与颗粒应力计算共享结果

**实施**：
1. 在 `Mechanicaloutput` 中计算 `disp_surf` 时，同时存储：
   ```julia
   variables["negative particle surface displacement"] = disp_surf_n
   variables["positive particle surface displacement"] = disp_surf_p
   ```

2. 在 `thermal_diffusion_stress_2D` 中读取：
   ```julia
   if haskey(variables, "negative particle surface displacement")
       # 使用位移方法（新路线）
       disp_surf_n = variables["negative particle surface displacement"]
       ε_particle_n = 3.0 .* disp_surf_n ./ param.NE.rs
       ε_macro_n = param.NE.eps_s .* ε_particle_n
   else
       # 回退到SOC方法（旧路线，向后兼容）
       ε_macro_n = (param.NE.Omega / 3.0) * param.NE.eps_s * Δsoc_n
   end
   ```

---

## 四、详细实施步骤

### 步骤1：理解当前数据流

#### 1.1 SPMe模型中的颗粒数据
在 `SPMe.jl` 和 `P2D.jl` 中：
- 每个节点有独立的颗粒浓度分布 `cs(r)`
- 颗粒网格：`case.mesh["negative particle"]` 和 `case.mesh["positive particle"]`

#### 1.2 Thermal2D模型中的单元
在 `thermal_diffusion_stress_2D` 中：
- 每个热单元有平均SOC：`variables["thermal2D element soc_n"]`
- **问题**：颗粒数据与热单元如何对应？

#### 1.3 映射关系（关键）
- **单SPMe模式**：整体一个颗粒 → 直接用
- **多SPMe模式**（`per_element_spme=true`）：每个热单元对应一个子SPMe → 逐单元提取颗粒数据

### 步骤2：扩展变量传递

#### 2.1 修改 `Mechanicaloutput` 函数

**文件**：`src/mechanical.jl`

**在第12-26行（SPM/SPMe分支）添加颗粒体积应变计算**：

```julia
function Mechanicaloutput(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        mesh_n = case.mesh["negative particle"]
        mesh_p = case.mesh["positive particle"]
        c_n = variables["negative particle lithium concentration"]
        c_p = variables["positive particle lithium concentration"]
        ...
        T = variables["temperature"]
        
        # 计算颗粒应力和位移
        stress_rn_center,stress_theta_n_surf,disp_surf_n,theta_Mn,csn_gs = Calstressdisp(param.NE, mesh_n, c_n, T)
        stress_rp_center,stress_theta_p_surf,disp_surf_p,theta_Mp,csp_gs = Calstressdisp(param.PE, mesh_p, c_p, T)
        
        # ✨ 新增：计算颗粒体积应变
        volumetric_strain_n = 3.0 * disp_surf_n / param.NE.rs
        volumetric_strain_p = 3.0 * disp_surf_p / param.PE.rs
        
        # 存储所有变量
        variables["negative particle center radial stress"] = stress_rn_center
        variables["positive particle center radial stress"] = stress_rp_center
        variables["negative particle surface tangential stress"] = stress_theta_n_surf
        variables["positive particle surface tangential stress"] = stress_theta_p_surf
        variables["negative particle surface displacement"] = disp_surf_n
        variables["positive particle surface displacement"] = disp_surf_p
        
        # ✨ 新增：存储体积应变
        variables["negative particle volumetric strain"] = volumetric_strain_n
        variables["positive particle volumetric strain"] = volumetric_strain_p
        
        ...
    end
end
```

#### 2.2 对P2D模型（逐节点颗粒）
在第30-107行（P2D分支）中，需要对每个节点的颗粒计算：

```julia
elseif case.opt.model == "P2D" || case.opt.model == "sP2D"
    ...
    volumetric_strain_n = zeros(Float64, mesh_ne.nlen)
    volumetric_strain_p = zeros(Float64, mesh_pe.nlen)
    
    for i = 1:mesh_ne.nlen
        mesh = PickElement(mesh_n, ...)
        cs = cs_n[...]
        stress_rn_center[i],stress_theta_n_surf[i],disp_surf_n[i],theta_Mn[i],csn_gs[...] = Calstressdisp(...)
        
        # ✨ 新增：计算体积应变
        volumetric_strain_n[i] = 3.0 * disp_surf_n[i] / param.NE.rs
    end
    
    for i = 1:mesh_pe.nlen
        ...
        volumetric_strain_p[i] = 3.0 * disp_surf_p[i] / param.PE.rs
    end
    
    # 存储
    variables["negative particle volumetric strain"] = volumetric_strain_n
    variables["positive particle volumetric strain"] = volumetric_strain_p
    ...
end
```

### 步骤3：修改宏观应力计算

#### 3.1 在 `thermal_diffusion_stress_2D` 中使用新变量

**文件**：`src/mechanical.jl`  
**位置**：第165-223行

**当前代码**（第174-200行）：
```julia
# 提取温度场和SOC分布
T_nodes = ...
soc_n_elem = variables["thermal2D element soc_n"]
soc_p_elem = variables["thermal2D element soc_p"]
soc_ref_n = param.NE.cs0 
soc_ref_p = param.PE.cs0

# 计算单元级别的温度和SOC
...
Δsoc_n_elem[e] = soc_n_elem[e] - soc_ref_n
Δsoc_p_elem[e] = soc_p_elem[e] - soc_ref_p
```

**修改为**（使用颗粒体积应变）：

```julia
function thermal_diffusion_stress_2D(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    ...
    
    # 提取温度场
    T_nodes = haskey(variables, "T_nodes") ? variables["T_nodes"] : fill(T0, mesh.nlen)
    T_nodes = isa(T_nodes, AbstractVector) ? T_nodes : T_nodes[:, end]
    
    # ✨ 新路线：优先使用颗粒体积应变
    if haskey(variables, "negative particle volumetric strain") && 
       haskey(variables, "positive particle volumetric strain")
        
        # 路径1：从颗粒位移计算（新方法）
        ε_particle_n = variables["negative particle volumetric strain"]
        ε_particle_p = variables["positive particle volumetric strain"]
        
        # 映射到热单元
        # (需要根据具体模型处理，见步骤4)
        ε_particle_n_elem = map_particle_to_thermal_elements(ε_particle_n, case)
        ε_particle_p_elem = map_particle_to_thermal_elements(ε_particle_p, case)
        
        println("  [化学应变] 使用颗粒体积膨胀法（新路线）")
        
    elseif haskey(variables, "thermal2D element soc_n")
        
        # 路径2：从SOC计算（旧方法，向后兼容）
        soc_n_elem = variables["thermal2D element soc_n"]
        soc_p_elem = variables["thermal2D element soc_p"]
        soc_ref_n = param.NE.cs0 
        soc_ref_p = param.PE.cs0
        
        Δsoc_n_elem = soc_n_elem .- soc_ref_n
        Δsoc_p_elem = soc_p_elem .- soc_ref_p
        
        # 转换为等效体积应变
        ε_particle_n_elem = param.NE.Omega .* param.NE.cs_max .* Δsoc_n_elem
        ε_particle_p_elem = param.PE.Omega .* param.PE.cs_max .* Δsoc_p_elem
        
        println("  [化学应变] 使用SOC代理法（旧路线）")
        
    else
        error("无法获取化学应变数据：缺少颗粒体积应变或SOC信息")
    end
    
    # 获取材料参数
    E_eff = ...
    ν_eff = ...
    α_eff = ...
    
    # ✨ 不再需要β（颗粒体积应变已直接计算）
    # 但保留eps_s用于均质化
    eps_s_n = param.NE.eps_s
    eps_s_p = param.PE.eps_s
    
    # 计算单元级别的初始应变
    ne = size(mesh.element, 1)
    epsilon_0_elem = zeros(Float64, ne)
    
    @inbounds for e in 1:ne
        # 热应变
        ε_thermal = α_eff * dT_elem[e]
        
        # 化学应变（均质化）
        ε_chem_n = eps_s_n * ε_particle_n_elem[e]
        ε_chem_p = eps_s_p * ε_particle_p_elem[e]
        
        # 总初始应变
        epsilon_0_elem[e] = ε_thermal + ε_chem_n + ε_chem_p
    end
    
    # 后续有限元求解保持不变
    K_mech = _assemble_mechanical_stiffness_2D(...)
    F_mech = _assemble_thermal_diffusion_load_2D(mesh, E_eff, ν_eff, α_eff, 1.0, 1.0, dT_elem, 
                                                  ε_chem_n_elem, ε_chem_p_elem)
    ...
end
```

### 步骤4：颗粒数据到热单元的映射

这是**关键步骤**，取决于模型架构：

#### 情况A：单SPMe模型
```julia
function map_particle_to_thermal_elements(ε_particle_scalar, case)
    ne = size(case.mesh["thermal2D"].element, 1)
    # 所有单元使用相同的颗粒体积应变
    return fill(ε_particle_scalar, ne)
end
```

#### 情况B：多SPMe模型（每热单元一个SPMe）
```julia
function map_particle_to_thermal_elements(ε_particle_array, case)
    # 如果已经是按热单元存储的，直接返回
    if length(ε_particle_array) == size(case.mesh["thermal2D"].element, 1)
        return ε_particle_array
    else
        # 需要从电化学节点映射到热单元
        # (具体实现取决于数据存储方式)
        error("需要实现颗粒→热单元映射")
    end
end
```

#### 情况C：P2D模型（空间分布的颗粒）
```julia
function map_particle_to_thermal_elements(ε_particle_nodes, case)
    # 方法1：单元平均
    ne = size(case.mesh["thermal2D"].element, 1)
    ε_elem = zeros(Float64, ne)
    
    # 从电化学1D网格节点映射到2D热单元
    # 需要利用layer weights和几何映射
    for e in 1:ne
        # 获取单元中心坐标
        center = jellyroll_element_centers(mesh)[e,:]
        
        # 判断属于哪一层
        layer = material_at(center[1], center[2], param_dim)
        
        if layer == :NE
            # 从负极节点插值
            ε_elem[e] = interpolate_from_1D_to_element(ε_particle_nodes, e, :negative)
        elseif layer == :PE
            ε_elem[e] = interpolate_from_1D_to_element(ε_particle_nodes, e, :positive)
        else
            ε_elem[e] = 0.0  # SP/PCC/NCC无颗粒
        end
    end
    
    return ε_elem
end
```

---

## 五、唯一性与可解性分析（新路线）

### 5.1 问题是否唯一可解？

**答案：是的** ✅

#### 输入量（已知）
1. **颗粒浓度分布** `cs_n(r)`, `cs_p(r)` - 来自电化学求解
2. **温度场** `T(x,y)` - 来自热求解
3. **材料参数** `E, ν, α, Ω, eps_s, rs` - 常数

#### 中间量（可计算）
通过 `Calstressdisp` 函数：
```
cs(r) → disp_surf
```
这是**唯一确定的**（线弹性球形颗粒的解析解）

#### 颗粒体积应变（唯一确定）
```
ε_particle = 3 * disp_surf / rs
```

#### 宏观应变（唯一确定）
```
ε_macro = eps_s * ε_particle
```

#### 宏观应力（唯一确定）
给定 `ε_thermal` 和 `ε_chemical` 后，2D平面应力问题有唯一解（见前文分析）

### 5.2 与SOC方法的比较

|  | SOC方法（旧） | 位移方法（新） |
|---|--------------|--------------|
| **输入** | SOC（标量） | cs(r)（分布） |
| **假设** | 颗粒内均匀 | 考虑梯度 |
| **计算路径** | 直接公式 | 颗粒力学 |
| **精度** | 近似 | 更准确 |
| **计算量** | 小 | 略大 |
| **唯一性** | ✅ | ✅ |

### 5.3 物理一致性检验

#### 检验1：极限情况
当颗粒内浓度均匀（无梯度）时：
```
cs(r) = cs_av = const
disp_surf = Ω * rs * cs_av / 3
ε_particle = 3 * disp_surf / rs = Ω * cs_av
ε_macro = eps_s * Ω * cs_av

→ 与SOC方法一致 ✅
```

#### 检验2：能量守恒
颗粒弹性势能 + 宏观弹性势能 = 化学势能变化
（可通过变分原理验证）

#### 检验3：尺度一致性
- 颗粒尺度：Ω的单位 [m³/mol]
- 宏观尺度：ε的单位 [-]（无量纲）
- 转换关系：通过cs（[mol/m³]）联系

---

## 六、实施计划

### 阶段一：核心功能实现（1-2天）

#### 任务1：扩展变量存储
- [x] 理论分析完成
- [ ] 修改 `Mechanicaloutput` - 存储 `volumetric_strain`
- [ ] 测试：打印验证 `ε_particle` 的计算

#### 任务2：修改宏观应力计算
- [ ] 实现双路径逻辑（优先新方法，回退旧方法）
- [ ] 处理不同模型（SPM/SPMe/P2D）的映射
- [ ] 测试：单元测试

### 阶段二：映射与验证（2-3天）

#### 任务3：实现颗粒→热单元映射
- [ ] 单SPMe：全局均匀映射
- [ ] 多SPMe：逐单元映射
- [ ] P2D：1D→2D空间映射

#### 任务4：创建对比验证
- [ ] 运行新旧方法对比脚本
- [ ] 分析差异（均匀vs梯度场景）
- [ ] 可视化对比

### 阶段三：扩展与优化（可选，1周）

#### 任务5：非均匀场景测试
- [ ] 大电流（强梯度）
- [ ] 温度非均匀
- [ ] 循环累积

#### 任务6：性能优化
- [ ] 缓存 `Calstressdisp` 结果
- [ ] 并行计算
- [ ] 减少重复调用

---

## 七、代码实现示例

### 示例1：修改 `Mechanicaloutput`（SPMe分支）

```julia
# 文件: src/mechanical.jl
# 位置: 第12-26行

function Mechanicaloutput(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    param = case.param
    if case.opt.model == "SPM" || case.opt.model == "SPMe"
        mesh_n = case.mesh["negative particle"]
        mesh_p = case.mesh["positive particle"]
        c_n = variables["negative particle lithium concentration"]
        c_p = variables["positive particle lithium concentration"]
        eta_n = variables["negative electrode overpotential"]
        eta_p = variables["positive electrode overpotential"]
        V_cell = variables["cell voltage"] 
        T = variables["temperature"]
        
        # 调用颗粒力学计算（已有）
        stress_rn_center,stress_theta_n_surf,disp_surf_n,theta_Mn,csn_gs = Calstressdisp(param.NE, mesh_n, c_n, T)
        stress_rp_center,stress_theta_p_surf,disp_surf_p,theta_Mp,csp_gs = Calstressdisp(param.PE, mesh_p, c_p, T)
        
        # ✨✨✨ 新增：计算颗粒体积应变 ✨✨✨
        volumetric_strain_n = 3.0 * disp_surf_n / param.NE.rs
        volumetric_strain_p = 3.0 * disp_surf_p / param.PE.rs
        
        # 修正过电位（已有）
        eta_p_new = eta_p - (2/3) * stress_theta_p_surf * param.PE.Omega 
        eta_n_new = eta_n - (2/3) * stress_theta_n_surf * param.NE.Omega
        V_cell_new = V_cell  - (2/3) * stress_theta_p_surf * param.PE.Omega + (2/3) * stress_theta_n_surf * param.NE.Omega
        
        # 存储所有变量
        variables["negative particle center radial stress"] = stress_rn_center
        variables["positive particle center radial stress"] = stress_rp_center
        variables["negative particle surface tangential stress"] = stress_theta_n_surf
        variables["positive particle surface tangential stress"] = stress_theta_p_surf
        variables["negative particle surface displacement"] = disp_surf_n
        variables["positive particle surface displacement"] = disp_surf_p
        variables["negative particle concentration at gauss point"] = csn_gs
        variables["positive particle concentration at gauss point"] = csp_gs
        variables["negative particle stress coupling diffusion coefficient"] = theta_Mn
        variables["positive particle stress coupling diffusion coefficient"] = theta_Mp
        
        # ✨✨✨ 新增：存储体积应变 ✨✨✨
        variables["negative particle volumetric strain"] = volumetric_strain_n
        variables["positive particle volumetric strain"] = volumetric_strain_p
        
        variables["negative electrode overpotential"] = eta_n_new
        variables["positive electrode overpotential"] = eta_p_new
        variables["cell voltage"] = V_cell_new[1]
    elseif case.opt.model == "P2D" || case.opt.model == "sP2D"
        # P2D分支（类似处理，略）
        ...
    end
    return variables
end
```

### 示例2：修改 `thermal_diffusion_stress_2D`（核心修改）

```julia
# 文件: src/mechanical.jl
# 位置: 第165-223行

function thermal_diffusion_stress_2D(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})
    @assert haskey(case.mesh, "thermal2D") "thermal2D mesh is required for 2D diffusion stress"
    mesh = case.mesh["thermal2D"]
    @assert mesh.type == "Q4" "diffusion_stress_2D requires Q4 mesh"
    
    param = case.param
    Tref = param.scale.T_ref
    T0 = hasproperty(param.cell, :T0) ? param.cell.T0 : 298.0 / Tref
    
    # 提取温度场
    T_nodes = haskey(variables, "T_nodes") ? variables["T_nodes"] : fill(T0, mesh.nlen)
    T_nodes = isa(T_nodes, AbstractVector) ? T_nodes : T_nodes[:, end]
    
    # 获取材料参数
    E_eff = (param.NE.E * param.NE.thickness + param.PE.E * param.PE.thickness) / (param.NE.thickness + param.PE.thickness)
    ν_eff = (param.NE.nu * param.NE.thickness + param.PE.nu * param.PE.thickness) / (param.NE.thickness + param.PE.thickness)
    α_eff = (param.NE.alphaT * param.NE.thickness + param.PE.alphaT * param.PE.thickness) / (param.NE.thickness + param.PE.thickness)
    eps_s_n = param.NE.eps_s
    eps_s_p = param.PE.eps_s
    
    # 计算单元级别的温度
    ne = size(mesh.element, 1)
    T_elem = zeros(Float64, ne)
    dT_elem = zeros(Float64, ne)
    
    @inbounds for e in 1:ne
        nodes = mesh.element[e, :]
        T_elem[e] = sum(T_nodes[nodes]) / length(nodes)
        dT_elem[e] = T_elem[e] - T0
    end
    
    # ✨✨✨ 关键修改：双路径化学应变计算 ✨✨✨
    ε_chem_n_elem = zeros(Float64, ne)
    ε_chem_p_elem = zeros(Float64, ne)
    
    if haskey(variables, "negative particle volumetric strain") && 
       haskey(variables, "positive particle volumetric strain")
        
        # ═══════════════════════════════════════
        # 路径1：颗粒体积膨胀法（新路线）
        # ═══════════════════════════════════════
        println("  [化学应变] 使用颗粒体积膨胀法（基于disp_surf）")
        
        ε_particle_n = variables["negative particle volumetric strain"]
        ε_particle_p = variables["positive particle volumetric strain"]
        
        # 映射到热单元
        if isa(ε_particle_n, Number)
            # 标量：所有单元相同（单SPMe模式）
            ε_particle_n_elem = fill(ε_particle_n, ne)
            ε_particle_p_elem = fill(ε_particle_p, ne)
        elseif length(ε_particle_n) == ne
            # 向量且长度匹配：直接使用（多SPMe模式）
            ε_particle_n_elem = ε_particle_n
            ε_particle_p_elem = ε_particle_p
        else
            # 需要映射（P2D模式）
            ε_particle_n_elem = map_particle_to_thermal_elements(ε_particle_n, case, :negative)
            ε_particle_p_elem = map_particle_to_thermal_elements(ε_particle_p, case, :positive)
        end
        
        # 均质化：宏观应变 = 固相体积分数 × 颗粒体积应变
        @inbounds for e in 1:ne
            ε_chem_n_elem[e] = eps_s_n * ε_particle_n_elem[e]
            ε_chem_p_elem[e] = eps_s_p * ε_particle_p_elem[e]
        end
        
    elseif haskey(variables, "thermal2D element soc_n") && 
           haskey(variables, "thermal2D element soc_p")
        
        # ═══════════════════════════════════════
        # 路径2：SOC代理法（旧路线，向后兼容）
        # ═══════════════════════════════════════
        println("  [化学应变] 使用SOC代理法（向后兼容）")
        
        soc_n_elem = variables["thermal2D element soc_n"]
        soc_p_elem = variables["thermal2D element soc_p"]
        soc_ref_n = param.NE.cs0 
        soc_ref_p = param.PE.cs0
        
        # 计算SOC变化
        Δsoc_n_elem = zeros(Float64, ne)
        Δsoc_p_elem = zeros(Float64, ne)
        
        @inbounds for e in 1:ne
            Δsoc_n_elem[e] = soc_n_elem[e] - soc_ref_n
            Δsoc_p_elem[e] = soc_p_elem[e] - soc_ref_p
        end
        
        # 转换为体积应变（假设均匀浓度）
        # ε_particle = Ω * Δc = Ω * c_max * ΔSOC
        ε_particle_n_elem = param.NE.Omega .* Δsoc_n_elem .* param.NE.cs_max
        ε_particle_p_elem = param.PE.Omega .* Δsoc_p_elem .* param.PE.cs_max
        
        # 均质化
        @inbounds for e in 1:ne
            ε_chem_n_elem[e] = eps_s_n * ε_particle_n_elem[e]
            ε_chem_p_elem[e] = eps_s_p * ε_particle_p_elem[e]
        end
        
    else
        error("无法计算化学应变：缺少颗粒体积应变或SOC数据")
    end
    
    # ═══════════════════════════════════════
    # 计算总初始应变（热 + 化学）
    # ═══════════════════════════════════════
    epsilon_0_elem = zeros(Float64, ne)
    
    @inbounds for e in 1:ne
        epsilon_0_elem[e] = α_eff * dT_elem[e] + ε_chem_n_elem[e] + ε_chem_p_elem[e]
    end
    
    # ═══════════════════════════════════════
    # 后续有限元求解（保持不变）
    # ═══════════════════════════════════════
    
    # 装配力学刚度矩阵
    K_mech = _assemble_mechanical_stiffness_2D(mesh, E_eff, ν_eff)
    
    # 装配热-扩散载荷向量（现在直接使用 ε_chem_elem）
    F_mech = _assemble_thermal_diffusion_load_2D_new(mesh, E_eff, ν_eff, α_eff, 
                                                     dT_elem, ε_chem_n_elem, ε_chem_p_elem)
    
    # 施加边界条件
    K_mech, F_mech = _apply_mechanical_BC_2D(K_mech, F_mech, mesh, case)
    
    # 求解位移场
    U_M = _solve_mechanical_displacement_2D(K_mech, F_mech, mesh.nlen)
    
    # 恢复应力场
    σ_xx, σ_yy, σ_xy, σ_vm, σ_thermal, σ_diffusion = _recover_stress_2D_new(U_M, mesh, E_eff, ν_eff, α_eff, 
                                                                              Tref, dT_elem, ε_chem_n_elem, ε_chem_p_elem)
    
    # 写入结果
    L_ref = hasproperty(param.scale, :L_th) ? param.scale.L_th : 1.0
    _write_mechanical_results!(variables, U_M, σ_xx, σ_yy, σ_xy, σ_vm, σ_thermal, σ_diffusion, L_ref)
    
    return variables
end
```

### 示例3：辅助函数 - 映射颗粒数据到热单元

```julia
# 文件: src/mechanical.jl（在文件末尾添加）

"""
    map_particle_to_thermal_elements(ε_particle, case, electrode)

将颗粒尺度的体积应变映射到热单元。

# 参数
- `ε_particle`: 颗粒体积应变（标量或向量）
- `case`: 案例对象
- `electrode`: :negative 或 :positive

# 返回
- 长度为热单元数的向量
"""
function map_particle_to_thermal_elements(ε_particle::Union{Float64,Vector{Float64}}, 
                                         case::Case, 
                                         electrode::Symbol)
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)
    
    # 情况1：标量输入（单SPMe）
    if isa(ε_particle, Float64)
        return fill(ε_particle, ne)
    end
    
    # 情况2：向量且长度匹配（多SPMe）
    if length(ε_particle) == ne
        return ε_particle
    end
    
    # 情况3：P2D模式 - 需要从1D电化学网格映射到2D热网格
    if electrode == :negative
        mesh_1d = case.mesh["negative electrode"]
    elseif electrode == :positive
        mesh_1d = case.mesh["positive electrode"]
    else
        error("Unknown electrode: $electrode")
    end
    
    # 简化方法：基于层判定
    ε_elem = zeros(Float64, ne)
    param_dim = case.param_dim
    centers = jellyroll_element_centers(mesh_th)
    
    @inbounds for e in 1:ne
        x, y = centers[e, 1], centers[e, 2]
        r = sqrt(x^2 + y^2)
        θ = atan(y, x)
        
        # 判断单元属于哪一层
        layer = material_at(r, θ, param_dim; logic=:spiral)
        
        if (electrode == :negative && layer == :NE) || 
           (electrode == :positive && layer == :PE)
            # 该单元属于目标电极
            # 简化：使用平均值（更精确的方法需要空间插值）
            ε_elem[e] = mean(ε_particle)
        else
            # 其他层无贡献
            ε_elem[e] = 0.0
        end
    end
    
    return ε_elem
end
```

---

## 八、验证与对比

### 8.1 理论验证

#### 验证1：均匀浓度极限
**场景**：颗粒内浓度均匀（`cs(r) = const`）

**预期**：
```
新方法 ≈ 旧方法
```

**验证脚本**：
```julia
# 创建均匀浓度分布
cs_uniform = fill(0.5, Nrn)

# 新方法
_, _, disp_surf, _, _ = Calstressdisp(param.NE, mesh, cs_uniform, T)
ε_new = eps_s * 3.0 * disp_surf / rs

# 旧方法
SOC = 0.5
ε_old = eps_s * Omega * cs_max * SOC

# 对比
@test isapprox(ε_new, ε_old, rtol=1e-6)
```

#### 验证2：梯度场差异
**场景**：大电流放电，颗粒表面浓度低、中心浓度高

**预期**：
```
cs_av < cs_center
→ disp_surf_new < disp_surf_uniform
→ ε_new < ε_old
```

**物理意义**：
- 浓度梯度导致等效平均浓度偏低
- 旧方法用表面SOC（或平均SOC）会高估体积膨胀

### 8.2 数值对比

#### 对比脚本框架
```julia
# example/chemical_strain_particle_vs_soc_comparison.jl

using JuBat, Plots

function main()
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    
    # 测试不同C-rate
    C_rates = [0.5, 1.0, 2.0, 5.0]
    
    results_new = []
    results_old = []
    
    for C in C_rates
        println("\n测试 $(C)C ...")
        
        # 运行新方法
        opt_new = deepcopy(opt)
        opt_new.use_particle_displacement = true  # 新增标志
        result_new = run_simulation(param_dim, opt_new, C)
        push!(results_new, result_new)
        
        # 运行旧方法
        opt_old = deepcopy(opt)
        opt_old.use_particle_displacement = false
        result_old = run_simulation(param_dim, opt_old, C)
        push!(results_old, result_old)
    end
    
    # 对比绘图
    plot_comparison(C_rates, results_new, results_old)
end
```

### 8.3 预期差异

| C-rate | 浓度梯度 | 新vs旧（应力峰值） | 物理解释 |
|--------|---------|-------------------|---------|
| 0.5C | 小 | ~99% | 接近均匀，差异小 |
| 1C | 中 | ~95% | 开始出现梯度 |
| 2C | 大 | ~85% | 梯度显著，旧方法高估 |
| 5C | 很大 | ~70% | 强梯度，差异明显 |

---

## 九、总结

### 9.1 新路线核心

**从颗粒位移到宏观应变的完整物理路径**：
```
cs(r) → disp_surf → ε_particle → ε_macro → σ(x,y)
```

### 9.2 关键优势

1. ✅ **物理路径完整**：无跳跃，每步有明确物理意义
2. ✅ **考虑浓度梯度**：利用颗粒内真实分布
3. ✅ **与现有架构融合**：复用 `Calstressdisp`
4. ✅ **双路径兼容**：保留旧方法作为回退
5. ✅ **唯一性保证**：理论严格，数值稳定

### 9.3 实施优先级

**立即执行**（1-2天）：
1. ✅ 理论文档完成
2. 修改 `Mechanicaloutput` - 存储 `volumetric_strain`
3. 修改 `thermal_diffusion_stress_2D` - 双路径逻辑
4. 单元测试

**短期完善**（1周）：
5. 实现P2D模式的映射函数
6. 创建对比验证脚本
7. 可视化新旧方法差异

**长期优化**（可选）：
8. 性能优化（缓存、并行）
9. 多循环累积效应分析
10. 与实验数据对比

---

## 十、参考文献

### 核心理论
1. **Christensen & Newman (2006)**: 颗粒应力的经典理论
2. **Bower et al. (2011)**: 有限应变与电化学耦合
3. **Zhang et al. (2007)**: 颗粒内应力的数值模拟

### 均质化方法
4. **Salvadori et al. (2014)**: 计算均质化框架
5. **Zuo & Zhao (2015)**: 多尺度模型中的体积分数效应

---

**文档版本**: v2.0 (新路线)  
**创建日期**: 2025-12-29  
**技术路线**: 颗粒位移法（基于disp_surf）  
**状态**: 理论完备，待代码实现 ✅
