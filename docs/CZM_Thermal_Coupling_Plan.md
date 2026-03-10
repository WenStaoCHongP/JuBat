# 电化学-热-内聚力全耦合模型修改计划

## 一、目标概述

实现损伤驱动的电化学-热-力学全耦合模型：
- **损伤 → 热阻增加 → 局部热集中**
- **断裂(D=1) → 内圈单元退出 → 电流/温度重分配**
- **容量衰减 → SOH监控 → 循环终止(SOH=80%)**

---

## 二、设计决策汇总

| 项目 | 决策 |
|------|------|
| 断裂判据 | D = 1（完全断裂） |
| 退出方式 | 突然退出（非渐进） |
| 退出范围 | 仅内圈单元退出电化学反应 |
| 电流处理 | 失效单元电流=0，由其他单元承担 |
| 失效性质 | 永久失效（不可恢复） |
| CZM模式 | Mode I only（只考虑法向） |
| 间隙导热 | 串联热阻模型（接触+空气间隙） |
| 节点温度 | 内聚力单元与相邻体积单元共享 |
| 网格策略 | 热网格划分时统一创建CZM网格 |
| 耦合策略 | 带固定点迭代的强交错求解 |
| 损伤更新 | 每个时间步更新一次 |
| 终止条件 | SOH ≤ 80%（容量保持率） |

---

## 三、间隙导热模型

### 3.1 物理模型

考虑脱粘产生的物理间隙（空气）与残余接触面的串联热阻：

$$h_{eff} = \left[ \frac{1}{h_0 \cdot (1 - D)} + \frac{\delta_n \cdot H(\delta_n)}{k_{air}} \right]^{-1}$$

其中：
- $h_0$：完好界面的接触换热系数 [W/(m²·K)]
- $D$：损伤变量 [0, 1]
- $\delta_n$：法向分离位移 [m]
- $H(\delta_n)$：Heaviside函数，$H(x) = 1$ 当 $x > 0$，否则为 $0$
- $k_{air}$：空气热导率 ≈ 0.026 W/(m·K)

### 3.2 物理含义

```
完好状态 (D=0, δ_n=0):
    h_eff = h_0  （纯接触传热）

部分损伤 (0<D<1, δ_n>0):
    接触面积减少 → 接触热阻增加
    出现间隙 → 空气热阻增加
    两者串联

完全断裂 (D=1, δ_n>0):
    h_eff → 0  （近似绝热）
```

### 3.3 参数建议值

| 参数 | 符号 | 建议值 | 单位 |
|------|------|--------|------|
| 完好接触换热系数 | h_0 | 1000~5000 | W/(m²·K) |
| 空气热导率 | k_air | 0.026 | W/(m·K) |
| 临界间隙（近似绝热） | δ_crit | ~100 μm | m |

---

## 四、代码修改计划

### 阶段1：数据结构扩展

#### 1.1 新增间隙导热参数 (`SetParams.jl`)

```julia
@with_kw mutable struct GapConductance
    h_0::Float64 = 2000.0       # 完好接触换热系数 [W/(m²·K)]
    k_air::Float64 = 0.026      # 空气热导率 [W/(m·K)]
    enabled::Bool = true        # 是否启用间隙导热模型
end
```

**修改文件**: `src/SetParams.jl`
- 在 `Cohesive` 结构后添加 `GapConductance` 结构
- 在 `Params` 结构中添加 `gap_conductance` 字段

#### 1.2 扩展单元状态 (`czm.jl`)

```julia
# 在 CohesiveMesh 中添加
mutable struct CohesiveMesh
    # ... 现有字段 ...
    
    # 新增：相邻体积单元映射
    inner_elements::Vector{Int64}   # 内圈体积单元索引
    outer_elements::Vector{Int64}   # 外圈体积单元索引
    
    # 新增：单元活跃状态
    element_active::Vector{Bool}    # 体积单元是否活跃
    fractured_interface::Vector{Bool}  # 内聚力单元是否断裂
end
```

**修改文件**: `src/czm.jl`

#### 1.3 扩展Option (`Option.jl`)

```julia
# 在 Option 中添加
czm_enabled::Bool = false              # 是否启用CZM
czm_mode_I_only::Bool = false          # 只考虑法向（已有）
czm_update_interval::Int = 1           # 损伤更新间隔（改为1）
czm_soh_threshold::Float64 = 0.8       # SOH终止阈值
czm_gap_conductance_enabled::Bool = true  # 启用间隙导热
```

**修改文件**: `src/Option.jl`

---

### 阶段2：网格统一

#### 2.1 修改热网格生成 (`Jellyrollmodel.jl`)

将CZM网格创建逻辑集成到热网格生成中：

```julia
function jellyroll_collector_seed_mesh(param_dim; nθ=60, gsorder=2, 
                                        create_czm::Bool=false)
    # ... 现有热网格生成代码 ...
    
    if create_czm
        # 识别界面节点对
        interface_pairs = _find_coincident_node_pairs(mesh, inner_nodes, outer_nodes, tol)
        
        # 创建内聚力单元
        cohesive_elements = _create_cohesive_elements_direct(mesh, interface_pairs)
        
        # 建立体积单元-内聚力单元映射
        inner_elements, outer_elements = _map_elements_to_interface(mesh, interface_pairs)
        
        # 将CZM信息附加到mesh
        mesh.czm_data = (
            cohesive_elements = cohesive_elements,
            interface_pairs = interface_pairs,
            inner_elements = inner_elements,
            outer_elements = outer_elements
        )
    end
    
    return mesh
end
```

**修改文件**: `src/Jellyrollmodel.jl`

#### 2.2 简化CZM网格创建 (`czm.jl`)

```julia
function create_czm_mesh_from_thermal(thermal_mesh::Mesh, param_dim; tol=1e-8)
    # 直接使用热网格中预先创建的CZM数据
    if haskey(thermal_mesh, :czm_data) && thermal_mesh.czm_data !== nothing
        czm_data = thermal_mesh.czm_data
        # 基于czm_data构建CohesiveMesh
        ...
    else
        # 回退到原有逻辑
        return create_czm_mesh(thermal_mesh, param_dim; tol=tol)
    end
end
```

**修改文件**: `src/czm.jl`

---

### 阶段3：间隙导热实现

#### 3.1 计算有效换热系数 (`czm.jl`)

```julia
"""
    compute_gap_conductance(D, δ_n, gap_params)

计算间隙有效换热系数。

# 公式
h_eff = [1/(h_0·(1-D)) + δ_n·H(δ_n)/k_air]^(-1)
"""
function compute_gap_conductance(D::Float64, δ_n::Float64, gap_params::GapConductance)
    if !gap_params.enabled
        return gap_params.h_0  # 不启用时返回完好值
    end
    
    # 避免除零
    contact_term = 1.0 / (gap_params.h_0 * max(1.0 - D, 1e-10))
    
    # Heaviside函数：只有张开时才有间隙
    gap_term = δ_n > 0 ? δ_n / gap_params.k_air : 0.0
    
    # 串联热阻
    h_eff = 1.0 / (contact_term + gap_term)
    
    return h_eff
end
```

**新增文件或添加到**: `src/czm.jl`

#### 3.2 修改热传导装配 (`ThermalDistributed.jl`)

```julia
function ThermalDistributed2D_with_CZM(case, variables, czm_mesh)
    # 1. 基础热矩阵（体积单元）
    MT, KT, FT = ThermalDistributed2D(case, variables)
    
    # 2. 添加界面热阻贡献
    if czm_mesh !== nothing && czm_mesh.n_cohesive > 0
        gap_params = case.param_dim.gap_conductance
        
        for (i, coh_elem) in enumerate(czm_mesh.cohesive_elements)
            D = czm_mesh.damage_states[i].D
            δ_n = czm_mesh.damage_states[i].δ_max_n
            
            # 计算有效换热系数
            h_eff = compute_gap_conductance(D, δ_n, gap_params)
            
            # 在界面节点间添加热耦合
            _add_interface_thermal_coupling!(KT, FT, coh_elem, h_eff, case)
        end
    end
    
    return MT, KT, FT
end

"""
在界面节点间添加热耦合（共享温度或热阻）
"""
function _add_interface_thermal_coupling!(KT, FT, coh_elem, h_eff, case)
    # 获取界面节点
    n_bottom = coh_elem.nodes_bottom  # 外圈节点
    n_top = coh_elem.nodes_top        # 内圈节点
    
    L = coh_elem.length  # 单元长度
    H = case.param_dim.cell.width  # 厚度方向
    
    # 界面面积
    A_interface = L * H
    
    # 热导贡献 = h_eff * A
    k_interface = h_eff * A_interface
    
    # 组装到全局矩阵
    # 对于每对界面节点：
    for (nb, nt) in zip(n_bottom, n_top)
        # KT[nb, nb] += k_interface
        # KT[nt, nt] += k_interface
        # KT[nb, nt] -= k_interface
        # KT[nt, nb] -= k_interface
        # 这实现了 q = h_eff * A * (T_bottom - T_top)
    end
end
```

**修改文件**: `src/ThermalDistributed.jl`

---

### 阶段4：电流重分配

#### 4.1 修改分流求解器 (`Solve.jl`)

```julia
function solve_branch_currents_with_damage(case, variables, czm_mesh, ...)
    ne = size(case.mesh["thermal2D"].element, 1)
    
    # 获取单元活跃状态
    element_active = get_element_active_state(czm_mesh, ne)
    
    # 计算活跃单元数量
    n_active = sum(element_active)
    
    if n_active == 0
        @warn "所有单元都已失效！"
        return zeros(ne), 0.0
    end
    
    # 只在活跃单元间分配电流
    I_total = case.opt.Current(t * case.param.scale.t0) / case.param.scale.I_typ
    
    # 重新计算面积权重（仅活跃单元）
    areas_active = [element_active[e] ? areas[e] : 0.0 for e in 1:ne]
    total_active_area = sum(areas_active)
    
    # 分配电流
    I_e = zeros(ne)
    for e in 1:ne
        if element_active[e]
            I_e[e] = I_total * areas_active[e] / total_active_area
        end
    end
    
    return I_e, Vc
end

"""
根据CZM断裂状态获取单元活跃状态。
断裂的内聚力单元会导致其内圈相邻单元失效。
"""
function get_element_active_state(czm_mesh, ne)
    element_active = fill(true, ne)
    
    if czm_mesh === nothing
        return element_active
    end
    
    for (i, coh_elem) in enumerate(czm_mesh.cohesive_elements)
        if czm_mesh.damage_states[i].fractured
            # 内圈单元失效
            inner_elem_idx = czm_mesh.inner_elements[i]
            if inner_elem_idx > 0 && inner_elem_idx <= ne
                element_active[inner_elem_idx] = false
            end
        end
    end
    
    return element_active
end
```

**修改文件**: `src/Solve.jl`

#### 4.2 修改热源计算 (`ThermalDistributed.jl`)

```julia
function compute_element_heat_sources_with_damage(case, variables, I_e, T_e, ne, 
                                                   element_active)
    q_elem = zeros(Float64, ne)
    
    for e in 1:ne
        if element_active[e]
            # 正常计算热源
            q_elem[e] = _compute_single_element_heat(case, variables, I_e[e], T_e[e], e)
        else
            # 失效单元无热源
            q_elem[e] = 0.0
        end
    end
    
    return q_elem
end
```

**修改文件**: `src/ThermalDistributed.jl`

---

### 阶段5：SOH监控与终止条件

#### 5.1 容量计算与SOH估算 (`CycleSolver.jl`)

```julia
"""
计算当前SOH（基于放电容量）
"""
function compute_soh(result::CyclingResult)
    if length(result.capacity_discharge) < 2
        return 1.0
    end
    
    initial_capacity = result.capacity_discharge[1]
    current_capacity = result.capacity_discharge[end]
    
    soh = current_capacity / initial_capacity
    return soh
end

"""
检查是否应该终止循环（SOH阈值）
"""
function should_terminate_cycling(result::CyclingResult, soh_threshold::Float64)
    soh = compute_soh(result)
    return soh <= soh_threshold
end
```

#### 5.2 修改循环求解器 (`CycleSolver.jl`)

```julia
function solve_cycling(case, cycle_opt, czm_mesh; ...)
    # ... 现有初始化代码 ...
    
    # 获取SOH阈值
    soh_threshold = hasproperty(case.opt, :czm_soh_threshold) ? 
                    case.opt.czm_soh_threshold : 0.8
    
    for cycle in 1:n_cycles
        # ... 现有循环代码 ...
        
        # 检查SOH终止条件
        if should_terminate_cycling(result, soh_threshold)
            if verbose
                soh = compute_soh(result)
                @printf("  ⚠ SOH = %.1f%% ≤ %.1f%%，达到终止条件\n", 
                        soh * 100, soh_threshold * 100)
            end
            break
        end
        
        # 检查断裂单元比例
        if czm_mesh !== nothing
            n_fractured = count(s -> s.fractured, czm_mesh.damage_states)
            if n_fractured >= czm_mesh.n_cohesive
                if verbose
                    println("  ⚠ 所有内聚力单元断裂，终止循环")
                end
                break
            end
        end
    end
    
    return result
end
```

**修改文件**: `src/CycleSolver.jl`

---

### 阶段6：损伤更新频率调整

#### 6.1 修改CZM更新间隔 (`CycleSolver.jl`)

```julia
function _solve_phase_internal(...)
    # ... 现有代码 ...
    
    # 获取CZM更新间隔（从opt读取，默认改为1）
    czm_update_interval = hasproperty(case.opt, :czm_update_interval) ? 
                          case.opt.czm_update_interval : 1
    
    while t < t_end_nd
        # ... 电化学求解 ...
        
        # CZM损伤计算（每个时间步）
        if czm_mesh !== nothing && czm_params !== nothing && 
           (step_count % czm_update_interval == 0)
            # ... CZM更新代码 ...
        end
        
        # ... 其余代码 ...
    end
end
```

**修改文件**: `src/CycleSolver.jl`

---

## 五、修改文件清单

| 序号 | 文件 | 修改内容 | 优先级 |
|------|------|----------|--------|
| 1 | `src/SetParams.jl` | 添加GapConductance结构 | P0 |
| 2 | `src/Option.jl` | 添加CZM相关选项 | P0 |
| 3 | `src/parameters/Jellyroll.jl` | 添加间隙导热默认参数 | P0 |
| 4 | `src/Jellyrollmodel.jl` | 统一热网格和CZM网格创建 | P1 |
| 5 | `src/czm.jl` | 添加间隙导热计算、简化网格创建 | P1 |
| 6 | `src/ThermalDistributed.jl` | 添加界面热阻、修改热源计算 | P2 |
| 7 | `src/Solve.jl` | 修改分流求解器考虑失效单元 | P2 |
| 8 | `src/CycleSolver.jl` | SOH监控、终止条件、更新频率 | P3 |
| 9 | `src/JuBat.jl` | 导出新函数和结构 | P3 |
| 10 | `example/czm_cycle_example.jl` | 更新示例 | P4 |

---

## 六、实现顺序

```
阶段1: 数据结构 (1-2天)
├── SetParams.jl: GapConductance
├── Option.jl: CZM选项
└── Jellyroll.jl: 默认参数

阶段2: 网格统一 (2-3天)
├── Jellyrollmodel.jl: 集成CZM网格创建
└── czm.jl: 简化接口

阶段3: 间隙导热 (2-3天)
├── czm.jl: compute_gap_conductance
└── ThermalDistributed.jl: 界面热阻

阶段4: 电流重分配 (2天)
├── Solve.jl: 分流求解器
└── ThermalDistributed.jl: 热源计算

阶段5: 循环控制 (1-2天)
├── CycleSolver.jl: SOH监控
└── CycleSolver.jl: 终止条件

阶段6: 测试验证 (2-3天)
├── 单元测试
└── 集成测试
```

---

## 七、关键接口变更

### 7.1 新增函数

```julia
# czm.jl
compute_gap_conductance(D, δ_n, gap_params) → h_eff

# ThermalDistributed.jl
ThermalDistributed2D_with_CZM(case, variables, czm_mesh) → MT, KT, FT
compute_element_heat_sources_with_damage(...) → q_elem

# Solve.jl
get_element_active_state(czm_mesh, ne) → element_active
solve_branch_currents_with_damage(...) → I_e, Vc

# CycleSolver.jl
compute_soh(result) → soh
should_terminate_cycling(result, threshold) → Bool
```

### 7.2 修改的函数签名

```julia
# 原：jellyroll_collector_seed_mesh(param_dim; nθ, gsorder)
# 新：jellyroll_collector_seed_mesh(param_dim; nθ, gsorder, create_czm=false)

# 原：create_czm_mesh(thermal_mesh, param_dim; tol)
# 新：create_czm_mesh_from_thermal(thermal_mesh, param_dim; tol)
```

---

## 八、预期结果

### 8.1 功能验证

1. **间隙导热**：损伤增加 → 界面热阻增加 → 局部温升
2. **电流重分配**：断裂 → 内圈失效 → 其他单元电流增加
3. **热集中**：断裂区域温度梯度增大
4. **容量衰减**：失效单元增多 → 有效容量下降 → SOH降低
5. **循环终止**：SOH ≤ 80% 时自动终止

### 8.2 输出数据

```julia
result = solve_cycling(case, cycle_opt, czm_mesh)

# 可获取的数据
result.capacity_discharge  # 每循环放电容量
result.D_max              # 最大损伤
result.n_fractured        # 断裂单元数
result.T_max              # 最高温度

# 新增
result.soh                # SOH历史
result.n_active_elements  # 活跃单元数历史
result.termination_reason # 终止原因
```

---

## 九、风险与注意事项

1. **数值稳定性**：间隙导热系数在D→1时趋近于零，需要设置最小值避免除零
2. **收敛性**：电流重分配可能导致局部电流密度过高，需要监控
3. **网格依赖性**：CZM结果可能依赖网格密度，需要进行网格收敛性研究
4. **计算效率**：每步更新CZM会增加计算量，需要评估性能影响
5. **边界情况**：如果所有单元都失效，需要优雅地处理

---

## 十、后续扩展（暂不实现）

1. **可视化**：损伤场、温度场、电流分布的动态可视化
2. **渐进退出**：基于损伤程度的渐进电流衰减
3. **双向脱粘**：外圈也可能退出（当前只有内圈）
4. **热失控预警**：基于温升速率的安全预警
5. **参数敏感性分析**：间隙导热参数对结果的影响