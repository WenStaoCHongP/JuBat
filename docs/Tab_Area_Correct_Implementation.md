# 极耳冷却的正确实现：使用 tab.area

## 核心纠正

### 错误理解（之前）

❌ **错误**：极耳冷却的散热面积 = 极耳节点的xy平面投影面积之和

```julia
# 错误的实现
node_areas = compute_node_areas(mesh)
for n in tab_nodes
    A_node = node_areas[n]  # ❌ 用节点投影面积
    K[n,n] += (2h/H) * A_node / L_th^2
end
```

### 正确理解（现在）

✅ **正确**：极耳冷却的散热面积 = `tab.area`（在 Jellyroll.jl 中定义）

```julia
# 正确的实现
tab_area = case.param_dim.tab.area  # ✅ 实际极耳面积
tab_nodes = jellyroll_tab_node_indices(...)  # 识别受影响的节点

# 将 tab_area 对应的散热功率分配到节点
```

## 物理机制

### 1. 极耳的实际散热面积

在 `Jellyroll.jl` 中定义：

```88:90:src/parameters/Jellyroll.jl
tab.width = 40e-3
tab.length = 0.75 * 99.06e-3
tab.area = tab.width * tab.length * 2
```

- `tab.width`：极耳宽度 [m]
- `tab.length`：极耳长度 [m]
- `tab.area`：极耳总散热面积（正负极耳之和）[m²]

**物理意义**：极耳与外界（如冷板）接触的实际面积。

### 2. 极耳节点识别函数的作用

`jellyroll_tab_node_indices(mesh, param_dim)` 的作用：

- ✅ 识别受极耳冷却**影响**的节点
- ✅ 这些节点位于螺旋线上（以直代曲）
- ❌ **不是**计算散热面积

**类比**：
- 极耳就像一个"热交换器"，有固定的散热面积 `tab.area`
- `tab_nodes` 是这个热交换器"影响"到的FEM节点
- 需要将热交换器的散热功率合理分配到这些节点

### 3. 总散热功率

极耳通过 `tab.area` 与环境对流换热：

$$
\dot{Q}_{\text{total}} = h_{\text{tab}} \cdot A_{\text{tab}} \cdot (T_{\text{tab}} - T_{\text{amb}})
$$

其中 $T_{\text{tab}}$ 是极耳的某种平均温度（或认为各节点温度不同）。

## 分配策略

### 策略1：均匀分配（简单）

将总散热功率均匀分配到所有极耳节点：

$$
\dot{Q}_i = \frac{h_{\text{tab}} \cdot A_{\text{tab}}}{n} (T_i - T_{\text{amb}})
$$

其中 $n = $ `length(tab_nodes)`。

**实现**：
```julia
n_nodes = length(tab_nodes)
coeff = h_tab * tab_area / (n_nodes * H)

for n in tab_nodes
    K[n,n] += coeff / L_th
    F[n] += coeff * T_amb / L_th
end
```

### 策略2：按节点面积权重分配（推荐）✅

按节点的xy平面面积权重分配，使得面积较大的节点分配更多散热：

$$
w_i = \frac{A_{\text{node},i}}{\sum_j A_{\text{node},j}}
$$

$$
\dot{Q}_i = h_{\text{tab}} \cdot (w_i \cdot A_{\text{tab}}) \cdot (T_i - T_{\text{amb}})
$$

**实现**：
```julia
node_areas = compute_node_areas(mesh)
total_node_area = sum(node_areas[tab_nodes])

for n in tab_nodes
    weight = node_areas[n] / total_node_area
    A_eff = tab_area * weight
    coeff = h_tab * A_eff / H
    
    K[n,n] += coeff / L_th
    F[n] += coeff * T_amb / L_th
end
```

**优点**：
- 物理更合理（面积大的节点应该承担更多散热）
- 避免节点数依赖性（如果网格加密，每个节点分配的散热会自动减少）

## 无量纲化

### 物理量

- 散热功率：$\dot{Q} = h_{\text{tab}} \cdot A_{\text{tab}} \cdot (T - T_{\text{amb}})$ [W]
- 体积散热率：$q_{\text{vol}} = \frac{\dot{Q}}{V} = \frac{h_{\text{tab}} \cdot A_{\text{tab}}}{H \cdot A_{\text{domain}}} (T - T_{\text{amb}})$

### 无量纲系数

对于单个节点 $i$，分配到的散热面积为 $A_{\text{eff},i} = w_i \cdot A_{\text{tab}}$：

$$
\text{coeff}_i = \frac{h_{\text{tab}} \cdot A_{\text{eff},i}}{H \cdot k_{\text{th}} \cdot L_{\text{th}}}
$$

**量纲检查**：
```
[h] = W/(m²·K)
[A] = m²
[H] = m
[k] = W/(m·K)
[L] = m

coeff = [W/(m²·K)] * [m²] / ([m] * [W/(m·K)] * [m])
      = [W/K] / [W/K]
      = 无量纲 ✓
```

## 与表面冷却的对比

### 表面冷却（surface）

**散热面积**：整个网格的xy平面投影面积

$$
A_{\text{surface}} = \sum_{\text{all elements}} A_{\text{elem}}
$$

**实现**：对所有单元高斯积分

```julia
for g in 1:ngs
    wt = (2h_surface/H) * wJ[g] / L_th^2
    for i, j in element_nodes
        K[i,j] += wt * N_i * N_j
    end
end
```

### 极耳冷却（tab）

**散热面积**：`tab.area`（实际极耳面积，通常远小于网格总面积）

$$
A_{\text{tab}} \ll A_{\text{surface}}
$$

**实现**：识别极耳节点，按权重分配

```julia
tab_area = case.param_dim.tab.area
tab_nodes = jellyroll_tab_node_indices(...)

for n in tab_nodes
    weight = node_areas[n] / total_node_area
    A_eff = tab_area * weight
    coeff = h_tab * A_eff / (H * k_th * L_th)
    
    K[n,n] += coeff
    F[n] += coeff * T_amb
end
```

## 典型数值

### Jellyroll.jl 中的参数

```julia
tab.width = 40e-3          # 40 mm
tab.length = 0.75 * 99.06e-3  # 约 74 mm
tab.area = 40e-3 * 74e-3 * 2  # 约 5.9e-3 m² (正负极耳之和)
```

### 网格总面积（估算）

```julia
cell.Rout = 0.0105 m
A_total ≈ π * R_out² ≈ 3.5e-4 m²
```

### 面积比

```julia
tab.area / A_total ≈ 5.9e-3 / 3.5e-4 ≈ 17
```

**注意**：极耳面积可能大于网格总面积！这是合理的，因为：
- 网格面积是xy平面投影
- 极耳面积是z方向的实际接触面积
- 它们的物理意义不同

## 验证方法

### 1. 能量守恒

总散热功率应等于：

$$
\dot{Q}_{\text{tab}} = h_{\text{tab}} \cdot A_{\text{tab}} \cdot (T_{\text{avg}} - T_{\text{amb}})
$$

其中 $T_{\text{avg}}$ 是极耳节点的平均温度。

### 2. 节点权重检查

所有节点权重之和应为1：

$$
\sum_{i \in \text{tab\_nodes}} w_i = 1
$$

### 3. 与表面冷却对比

如果 $h_{\text{tab}} = h_{\text{surface}}$ 且 $A_{\text{tab}}$ 足够大，极耳冷却应比表面冷却更强。

## 代码实现

### 完整函数

```julia
function _apply_cool_tab!(KT, FT, mesh, case, t)
    # 识别受影响的节点
    pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
    tab_nodes = unique(vcat(pos_idx, neg_idx))
    
    isempty(tab_nodes) && return
    
    # 参数
    h_tab = case.opt.h_tab
    tab_area = case.param_dim.tab.area  # ✅ 关键：使用实际极耳面积
    H = case.param_dim.cell.width
    
    scale = case.param_dim.scale
    k_th, L_th, T_ref = scale.k_th, scale.L_th, scale.T_ref
    T_amb_nd = case.param_dim.cell.T_amb / T_ref
    
    # 计算节点面积权重
    node_areas = _compute_node_areas(mesh)
    total_node_area = sum(node_areas[tab_nodes])
    
    # 按权重分配
    for n in tab_nodes
        weight = node_areas[n] / total_node_area
        A_eff = tab_area * weight
        
        # 无量纲系数
        coeff_nd = h_tab * A_eff / (H * k_th * L_th)
        
        KT[n, n] += coeff_nd
        FT[n] += coeff_nd * T_amb_nd
    end
end
```

## 关键要点总结

1. ✅ **散热面积 = `tab.area`**：在 Jellyroll.jl 中定义的实际极耳面积

2. ✅ **节点识别 ≠ 面积计算**：`jellyroll_tab_node_indices` 只识别受影响的节点

3. ✅ **分配策略**：按节点面积权重分配 `tab.area` 对应的散热功率

4. ✅ **物理量纲**：
   - 表面冷却：$2h/H$ [W/(m³·K)]（体积散热率）
   - 极耳冷却：$h \cdot A_{\text{tab}} / H$ [W/(m·K)]（线散热率，分配后）

5. ✅ **数值稳定**：系数量级仍为 $O(10^{-3})$ to $O(1)$，数值稳定

## 修正前后对比

| 特性 | 修正前（错误）| 修正后（正确）|
|------|-------------|-------------|
| 散热面积 | 节点投影面积之和 | `tab.area` |
| 节点识别作用 | 计算散热面积 | 识别受影响的节点 |
| 分配策略 | 直接用节点面积 | 按权重分配 `tab.area` |
| 物理意义 | 不清晰 | 清晰（实际极耳接触面积）|
