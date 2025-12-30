# 极耳边界条件使用指南

## 快速开始

### 方案1：极耳强化散热（推荐）✅

适用于极耳区域有额外冷却措施的场景。

```julia
# 在 testexample.jl 中添加
opt.tab_bc_type = "tab_convection"
opt.h_tab = 100.0  # W/(m²·K)，根据实际情况调整
```

**数值稳定，无矩阵病态问题！**

### 方案2：一般表面散热

适用于整体冷却，所有节点均匀散热。

```julia
opt.tab_bc_type = "surface_convection"
opt.h_surface = 10.0  # W/(m²·K)
```

### 方案3：惩罚法（不推荐）

仅用于必须强制固定温度的场景。

```julia
opt.tab_bc_type = "penalty"
opt.tab_penalty = 1e6  # ⚠️ 必须降低到 1e6-1e8
```

## 详细配置

### 1. 修改 testexample.jl

在创建 `Option` 后添加：

```julia
opt = JuBat.Option()
opt.model = "SPMe"
# ... 其他设置 ...

# ✨ 新增：选择极耳边界类型
opt.tab_bc_type = "tab_convection"  # 或 "surface_convection" 或 "penalty"

# 如果选择 tab_convection
opt.h_tab = 100.0  # 极耳换热系数 [W/(m²·K)]

# 如果选择 surface_convection
opt.h_surface = 10.0  # 表面换热系数 [W/(m²·K)]

# 如果选择 penalty（不推荐）
opt.tab_penalty = 1e6  # 降低惩罚值
```

### 2. 换热系数取值建议

| 冷却条件 | h [W/(m²·K)] | 适用场景 |
|---------|--------------|---------|
| 自然对流 | 5 - 15 | 空气自然冷却 |
| 强制风冷 | 20 - 100 | 风扇冷却 |
| 液冷接触 | 100 - 1000 | 极耳与冷板接触 |
| 相变冷却 | 1000 - 5000 | PCM或蒸发冷却 |

### 3. 参数敏感性分析

运行参数扫描：

```julia
# 对比不同换热系数的影响
h_values = [10.0, 50.0, 100.0, 200.0, 500.0]

for h in h_values
    opt.h_tab = h
    result = JuBat.Solve(case)
    # 记录温度分布...
end
```

## 理论对照

### 对流边界条件

**物理方程**：
$$
-k \frac{\partial T}{\partial n} = h(T - T_{\text{amb}})
$$

**弱形式贡献**：
$$
\int_{\Gamma} h(T - T_{\text{amb}}) N_i \, d\Gamma
$$

**离散后**：
- 刚度矩阵：$K_{ij} += \int h N_i N_j \, dA$
- 载荷向量：$F_i += \int h T_{\text{amb}} N_i \, dA$

### Z方向冷却的2D实现

由于2D模型只有 x-y 平面，z方向的表面积需要"投影"到节点：

**节点投影面积**：
$$
A_z(i) = A_{\text{voronoi}}(i) \times H
$$

其中：
- $A_{\text{voronoi}}(i)$: 节点在 x-y 平面的影响面积
- $H$: 电池高度（z方向）

**等效对流项**：
$$
K_{ii}^{\text{conv,z}} = h \cdot A_z(i) / (k_{\text{th}} \cdot L_{\text{th}})
$$

$$
F_i^{\text{conv,z}} = h \cdot T_{\text{amb}} \cdot A_z(i) / (k_{\text{th}} \cdot L_{\text{th}})
$$

### 无量纲化

引入 Biot 数：
$$
Bi_z = \frac{h \cdot H}{k_{\text{th}}}
$$

无量纲刚度矩阵：
$$
K_{ii}^* = Bi_z \cdot \frac{A_{\text{voronoi}}(i)}{L_{\text{th}}^2}
$$

## 代码修改位置

### 原始代码（ThermalDistributed.jl）

**修改前**：
```julia
function _apply_tab_bc!(KT, FT, mesh, case, t)
    # ... 惩罚法实现 ...
    penalty = 1e12  # ❌ 数值不稳定
    for n in tab_nodes
        KT[n, n] += penalty
        FT[n] += penalty * T_tab_nd
    end
end
```

**修改后**：
```julia
function _apply_tab_bc!(KT, FT, mesh, case, t)
    # 根据 opt.tab_bc_type 选择
    if bc_type == "tab_convection"
        _apply_tab_convection_bc!(KT, FT, mesh, case, t)
    elseif bc_type == "surface_convection"
        _apply_surface_convection_bc!(KT, FT, mesh, case, t)
    else
        _apply_tab_bc_penalty!(KT, FT, mesh, case, t)
    end
end
```

### 新增函数

**1. 极耳强化散热**（仅极耳节点）：
```julia
function _apply_tab_convection_bc!(KT, FT, mesh, case, t)
    # 1. 识别极耳节点
    tab_nodes = jellyroll_tab_node_indices(...)
    
    # 2. 计算节点面积
    node_areas = _compute_node_areas(mesh)
    
    # 3. 施加对流边界
    for n in tab_nodes
        A_z_nd = node_areas[n] * H / L_th^2
        KT[n, n] += Bi_z_tab * A_z_nd
        FT[n] += Bi_z_tab * T_amb_nd * A_z_nd
    end
end
```

**2. 一般表面散热**（所有节点）：
```julia
function _apply_surface_convection_bc!(KT, FT, mesh, case, t)
    # 对所有节点，通过高斯积分施加对流项
    for g in 1:ngs
        e = ele[g]
        nodes = mesh.element[e, :]
        wt = conv_factor * wJ[g]
        
        for i in 1:nn_per_elem
            ni = nodes[i]
            for j in 1:nn_per_elem
                nj = nodes[j]
                KT[ni, nj] += wt * Ni[g,i] * Ni[g,j]
            end
            FT[ni] += wt * T_amb_nd * Ni[g,i]
        end
    end
end
```

## 验证方法

### 1. 能量守恒检查

```julia
# 计算总产热
Q_gen = sum(q_elements .* volumes)

# 计算对流散热
Q_conv = sum(h * (T_nodes - T_amb) .* areas_z)

# 检查守恒
energy_balance = Q_gen - Q_conv - Q_storage
@assert abs(energy_balance) < tolerance
```

### 2. 温度分布合理性

- 极耳区域温度应略低（如果 h_tab 较大）
- 温度场应连续光滑
- 最高温度通常在内部远离极耳的位置

### 3. 对比实验数据

如果有实验测温数据：
```julia
T_exp = [...]  # 实验温度
T_sim = [...]  # 仿真温度

RMSE = sqrt(mean((T_exp - T_sim).^2))
println("温度RMSE: $RMSE K")
```

## 故障排除

### 问题1：仍然出现 NaN

**原因**：可能还在使用旧的惩罚法，或参数设置错误。

**解决**：
```julia
# 检查设置
@show opt.tab_bc_type
@show opt.h_tab  # 或 opt.h_surface

# 如果未定义，手动设置
opt.tab_bc_type = "tab_convection"
opt.h_tab = 100.0
```

### 问题2：温度过高或过低

**原因**：换热系数 h 不合理。

**解决**：
- 检查 h 的数量级（应在 5-500 范围内）
- 参考上表中的典型值
- 进行敏感性分析

### 问题3：求解很慢

**原因**：对流项可能导致时间步长减小。

**解决**：
```julia
# 适当放宽时间步长限制
opt.dt = [1.0, 20.0]  # [dt_min, dt_max]
```

## 示例脚本

运行对比测试：
```bash
julia example/testexample_tab_convection.jl
```

该脚本会自动对比三种边界条件，生成对比图表。

## 参考文献

详见 `docs/Tab_Convection_BC_Theory.md` 完整理论推导。
