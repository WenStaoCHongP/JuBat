# 极耳冷却实现的再次修正（V2）

## 核心纠正

### 修正前（V1，错误）

❌ **错误理解**：
- 极耳冷却的散热面积 = 极耳节点的xy平面投影面积之和
- 每个节点贡献自己的投影面积

```julia
// 错误实现
node_areas = compute_node_areas(mesh)
for n in tab_nodes
    A_node = node_areas[n]  // ❌ 节点投影面积
    K[n,n] += (2h_tab/H) * A_node / L_th^2
end
```

### 修正后（V2，正确）

✅ **正确理解**：
- 极耳冷却的散热面积 = `tab.area`（在 Jellyroll.jl 中定义）
- `jellyroll_tab_node_indices` 识别受影响的节点（不是计算散热面积）
- 将 `tab.area` 对应的散热功率按权重分配到节点

```julia
// 正确实现
tab_area = case.param_dim.tab.area  // ✅ 实际极耳面积
tab_nodes = jellyroll_tab_node_indices(...)

node_areas = compute_node_areas(mesh)
total_node_area = sum(node_areas[tab_nodes])

for n in tab_nodes
    weight = node_areas[n] / total_node_area  // 权重
    A_eff = tab_area * weight  // 分配的散热面积
    
    coeff = h_tab * A_eff / (H * k_th * L_th)
    K[n,n] += coeff
    F[n] += coeff * T_amb
end
```

## 物理意义对比

### 表面冷却（surface）

**散热机制**：电池上下表面（z方向）与环境对流

**散热面积**：整个网格的xy平面投影面积
$$
A_{\text{surface}} = \sum_{\text{all elements}} A_{\text{elem}} = \text{网格总面积}
$$

**实现方式**：对所有单元高斯积分
```julia
for g in 1:ngs  // 所有高斯点
    wt = (2h_surface/H) * wJ[g] / L_th^2
    K[i,j] += wt * N_i * N_j  // 面积分
end
```

**数值特征**：
- 分布式散热，作用于整个域
- 体积散热率：$q = 2h_{\text{surface}}/H$ [W/(m³·K)]

### 极耳冷却（tab）

**散热机制**：极耳与外界（如冷板）接触散热

**散热面积**：`tab.area`（实际极耳接触面积）
$$
A_{\text{tab}} = \text{tab.width} \times \text{tab.length} \times 2
$$

在 Jellyroll.jl 中：
```julia
tab.width = 40e-3          // 40 mm
tab.length = 0.75 * 99.06e-3  // 约 74 mm
tab.area = 40e-3 * 74e-3 * 2  // 约 5.9e-3 m²
```

**实现方式**：识别受影响的节点，按权重分配
```julia
tab_area = case.param_dim.tab.area  // 固定值
tab_nodes = jellyroll_tab_node_indices(...)  // 受影响的节点

// 按节点面积权重分配
for n in tab_nodes
    weight = node_areas[n] / total_node_area
    A_eff = tab_area * weight
    coeff = h_tab * A_eff / (H * k_th * L_th)
    K[n,n] += coeff
end
```

**数值特征**：
- 集中散热，仅作用于极耳节点
- 总散热功率：$\dot{Q} = h_{\text{tab}} \cdot A_{\text{tab}} \cdot (T - T_{\text{amb}})$

## 关键区别

| 特性 | 表面冷却 | 极耳冷却 |
|------|---------|---------|
| **散热面积来源** | 网格面积（计算得到）| `tab.area`（参数定义）|
| **作用范围** | 整个域 | 极耳节点邻域 |
| **节点识别作用** | 无（所有节点）| 识别受影响的节点 |
| **实现方式** | 高斯积分 | 权重分配 |
| **散热强度** | 均匀分布 | 集中在极耳 |
| **物理模型** | 体积热汇（2h/H）| 实际接触散热 |

## 为什么不能用节点投影面积？

### 问题1：物理意义不对应

- 节点投影面积 = xy平面的几何面积
- 极耳散热面积 = z方向的接触面积
- 两者的物理意义完全不同！

### 问题2：面积量级错误

**典型数值**：
- 网格总投影面积：$\sim 3.5 \times 10^{-4}$ m² (π R_out²)
- 极耳实际散热面积：$\sim 5.9 \times 10^{-3}$ m² (tab.area)
- 极耳面积 > 网格总面积的 **17 倍**！

如果用节点投影面积，极耳散热会被严重低估。

### 问题3：网格依赖性

如果用节点投影面积：
- 网格加密 → 极耳节点数增加 → 每个节点面积减小
- 但总散热功率应该是固定的（由 `tab.area` 决定）
- 会导致网格依赖性问题

使用 `tab.area` 并按权重分配：
- 网格加密 → 极耳节点数增加 → 每个节点分配的面积自动减小
- 总散热功率保持 $h_{\text{tab}} \cdot A_{\text{tab}}$
- 无网格依赖性 ✓

## 节点识别函数的真正作用

### `jellyroll_tab_node_indices(mesh, param_dim)`

**功能**：
识别受极耳冷却**影响**的节点

**输入**：
- `mesh`：网格数据
- `param_dim.tab.theta_pos`：正极耳角度位置
- `param_dim.tab.theta_neg`：负极耳角度位置
- `param_dim.tab.width`：极耳宽度（用于确定影响范围）

**输出**：
- `pos_idx`：正极耳影响的节点索引
- `neg_idx`：负极耳影响的节点索引

**物理解释**：
- 极耳位于螺旋线的某个位置（由 `theta_pos/neg` 确定）
- 极耳的散热会影响附近的节点
- `width` 决定影响范围的大小
- 返回的节点是螺旋线上的离散点（以直代曲）

**不是**：
- ❌ 计算散热面积
- ❌ 确定极耳的几何形状
- ❌ 网格剖分的一部分

## 分配策略

### 为什么按节点面积权重分配？

总散热功率固定：
$$
\dot{Q}_{\text{total}} = h_{\text{tab}} \cdot A_{\text{tab}} \cdot (T_{\text{avg}} - T_{\text{amb}})
$$

需要分配到 n 个节点：
$$
\dot{Q}_{\text{total}} = \sum_{i=1}^{n} \dot{Q}_i
$$

**方案1：均匀分配**
$$
\dot{Q}_i = \frac{\dot{Q}_{\text{total}}}{n}
$$

问题：与节点的"代表性"无关，不合理。

**方案2：按节点面积权重分配**（推荐）✅
$$
\dot{Q}_i = \frac{A_{\text{node},i}}{\sum_j A_{\text{node},j}} \cdot \dot{Q}_{\text{total}}
$$

优点：
- 面积大的节点"代表"更多区域，应该分配更多散热
- 无网格依赖性
- 物理更合理

### 权重计算

```julia
node_areas = compute_node_areas(mesh)  // 每个节点的xy平面面积
total_node_area = sum(node_areas[tab_nodes])  // 极耳节点总面积

for n in tab_nodes
    weight[n] = node_areas[n] / total_node_area  // 归一化权重
end

// 验证：sum(weight) = 1.0
```

然后：
```julia
A_eff[n] = tab_area * weight[n]  // 节点 n 分配到的散热面积
```

## 无量纲化

### 物理系数

对于节点 n：
$$
\text{coeff}_n = \frac{h_{\text{tab}} \cdot A_{\text{eff},n}}{H \cdot k_{\text{th}} \cdot L_{\text{th}}}
$$

其中：
- $A_{\text{eff},n} = w_n \cdot A_{\text{tab}}$
- $w_n = A_{\text{node},n} / \sum_j A_{\text{node},j}$

### 量纲检查

```
[h_tab] = W/(m²·K)
[A_eff] = m²
[H] = m
[k_th] = W/(m·K)
[L_th] = m

coeff = [W/(m²·K)] · [m²] / ([m] · [W/(m·K)] · [m])
      = [W/K] / [W/K]
      = 无量纲 ✓
```

### 数值量级

假设：
- $h_{\text{tab}} = 100$ W/(m²·K)
- $A_{\text{tab}} = 5.9 \times 10^{-3}$ m²
- $H = 0.07$ m
- $k_{\text{th}} = 1.0$ W/(m·K)
- $L_{\text{th}} = 0.01$ m
- $n = 27$ 个节点

平均每个节点：
$$
\text{coeff}_{\text{avg}} = \frac{100 \times (5.9 \times 10^{-3}/27)}{0.07 \times 1.0 \times 0.01} \approx 0.3
$$

数值稳定 ✓

## 代码实现

### 完整函数（最终版本）

```julia
function _apply_cool_tab!(KT, FT, mesh, case, t)
    # 1. 识别受影响的节点
    pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
    tab_nodes = unique(vcat(pos_idx, neg_idx))
    
    isempty(tab_nodes) && return
    
    # 2. 获取参数
    h_tab = case.opt.h_tab
    tab_area = case.param_dim.tab.area  # ✅ 实际极耳散热面积
    H = case.param_dim.cell.width
    
    scale = case.param_dim.scale
    k_th, L_th = scale.k_th, scale.L_th
    T_amb_nd = case.param_dim.cell.T_amb / scale.T_ref
    
    # 3. 计算节点面积权重
    node_areas = _compute_node_areas(mesh)
    total_node_area = sum(node_areas[tab_nodes])
    
    # 4. 按权重分配散热功率
    for n in tab_nodes
        # 节点权重
        weight = node_areas[n] / total_node_area
        
        # 分配的散热面积
        A_eff = tab_area * weight
        
        # 无量纲系数
        coeff_nd = h_tab * A_eff / (H * k_th * L_th)
        
        # 装配到矩阵
        KT[n, n] += coeff_nd
        FT[n] += coeff_nd * T_amb_nd
    end
end
```

## 验证方法

### 1. 权重和检查

```julia
weights = [node_areas[n] / total_node_area for n in tab_nodes]
@assert abs(sum(weights) - 1.0) < 1e-10 "权重和应为1"
```

### 2. 总散热功率

如果所有极耳节点温度相同为 $T$：
$$
\dot{Q}_{\text{total}} = h_{\text{tab}} \cdot A_{\text{tab}} \cdot (T - T_{\text{amb}})
$$

应该等于各节点散热功率之和。

### 3. 网格无关性

加密网格，极耳节点数增加，但总散热功率应保持不变。

## 总结

### 核心要点

1. ✅ **散热面积 = `tab.area`**：在 Jellyroll.jl 中定义的实际极耳面积

2. ✅ **节点识别 ≠ 面积计算**：`jellyroll_tab_node_indices` 只识别受影响的节点

3. ✅ **分配策略**：按节点面积权重分配 `tab.area` 对应的散热功率

4. ✅ **无网格依赖**：网格加密时，总散热功率保持不变

5. ✅ **数值稳定**：系数量级 $O(10^{-1})$ to $O(1)$

### 与表面冷却的区别

| 特性 | 表面冷却 | 极耳冷却 |
|------|---------|---------|
| 散热面积 | 网格总面积（计算）| `tab.area`（参数）|
| 作用范围 | 整个域 | 极耳节点 |
| 实现方式 | 高斯积分 | 权重分配 |
| 散热模型 | 体积热汇（2h/H）| 接触散热 |

### 修正历史

- **V0**：惩罚法（删除，数值不稳定）
- **V1**：使用节点投影面积（错误）
- **V2**：使用 `tab.area` + 权重分配（正确）✅
