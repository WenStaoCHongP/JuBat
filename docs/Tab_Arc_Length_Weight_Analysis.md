# 极耳冷却：弧长权重 vs 均匀分配

## 问题分析

### 物理机制

极耳节点是螺旋**线**上的离散点：
- 相邻节点之间有一段**弧长**
- 每个节点"代表"一段螺旋弧
- 不同节点代表的弧长**不相等**（外层>内层）

### 分配策略对比

#### 方案1：均匀分配（V3初版）

```julia
coeff_per_node = h_tab * tab_area / (n_nodes * H * k_th * L_th)
```

**假设**：每个节点代表相等的弧长

**问题**：
- 忽略了节点间距的差异
- 外层节点间距大，内层节点间距小
- 物理上不合理

#### 方案2：弧长权重分配（正确）✅

```julia
arc_lengths = compute_arc_lengths(tab_nodes, mesh)
total_arc_length = sum(arc_lengths)

for i in 1:length(tab_nodes)
    weight = arc_lengths[i] / total_arc_length
    coeff = h_tab * tab_area * weight / (H * k_th * L_th)
    K[n,n] += coeff
end
```

**假设**：每个节点按其代表的弧长分配散热

**优点**：
- 考虑节点间距差异
- 外层节点（间距大）分配更多散热
- 物理更合理

## 几何分析

### 节点分布特征

极耳节点识别逻辑：
```julia
if (Rin <= r <= Rout) &&           # 整个卷绕半径范围
   (θ_start <= θ_cum <= θ_end)    # 极耳角度范围
    push!(idx, i)
end
```

识别出的节点分布：
```
外层 (r ≈ Rout)  ○--------○--------○  间距大
中层 (r ≈ R_mid) ○------○------○------○  间距中等
内层 (r ≈ Rin)   ○----○----○----○----○  间距小
         │        │
         ↓        ↓
    角度范围 [θ_start, θ_end]
```

**关键**：节点间距随半径增大而增大！

### 弧长估算

对于阿基米德螺旋 $r = a + b\theta$：

弧长微元：
$$
ds = \sqrt{(dr)^2 + (r d\theta)^2} = \sqrt{b^2 + r^2} \, d\theta
$$

节点 i（半径 $r_i$）代表的弧长：
$$
s_i \approx \sqrt{b^2 + r_i^2} \cdot \Delta\theta_i
$$

其中 $\Delta\theta_i$ 是节点 i 的角度"影响范围"。

### 简化方法：直线距离

由于节点密集，螺旋弧 ≈ 直线段（以直代曲）：

```julia
# 节点i代表的弧长 ≈ 与相邻节点的直线距离平均
s[i] = (dist(node[i], node[i-1]) + dist(node[i], node[i+1])) / 2
```

边界节点：
```julia
s[1] = dist(node[1], node[2])
s[end] = dist(node[end], node[end-1])
```

## 数值影响分析

### 典型数值

假设：
- 内层半径：$R_{\text{in}} = 2$ mm
- 外层半径：$R_{\text{out}} = 10$ mm
- 极耳角度范围：$\Delta\theta = 0.1$ rad（约6°）

节点间距估算：
- 内层：$\Delta s_{\text{in}} \approx R_{\text{in}} \cdot \Delta\theta / n_{\text{layer}} \approx 2 \times 0.1 / 5 = 0.04$ mm
- 外层：$\Delta s_{\text{out}} \approx R_{\text{out}} \cdot \Delta\theta / n_{\text{layer}} \approx 10 \times 0.1 / 5 = 0.2$ mm

**间距比**：$\Delta s_{\text{out}} / \Delta s_{\text{in}} = 5$

### 系数差异

**均匀分配**：
```julia
coeff_in = coeff_out = h * A_tab / (n_nodes * H)
```

**弧长权重分配**：
```julia
coeff_in = h * A_tab * (s_in / s_total) / H
coeff_out = h * A_tab * (s_out / s_total) / H

coeff_out / coeff_in = s_out / s_in ≈ 5
```

**影响**：
- 外层节点散热强度是内层的 **5倍**（弧长权重）
- 均匀分配则完全相同

## 对组装的影响

### 刚度矩阵

**均匀分配**：
```julia
K[n_in, n_in] += c
K[n_out, n_out] += c  # 相同
```

**弧长权重**：
```julia
K[n_in, n_in] += c * w_in
K[n_out, n_out] += c * w_out  # w_out >> w_in
```

**影响**：
- 外层节点对环境耦合更强
- 矩阵对角占优性改变
- 可能影响收敛性

### 载荷向量

**均匀分配**：
```julia
F[n_in] += c * T_amb
F[n_out] += c * T_amb  # 相同
```

**弧长权重**：
```julia
F[n_in] += c * w_in * T_amb
F[n_out] += c * w_out * T_amb  # w_out >> w_in
```

**影响**：
- 外层节点受环境影响更大
- 温度分布改变

### 温度场

**物理预期**：
- 外层更接近环境温度（散热强）
- 内层温度更高（散热弱）

**均匀分配**：可能无法正确反映这一差异

**弧长权重**：正确反映物理机制 ✅

## 合理性判断

### 均匀分配是否合理？

**支持**：
1. 如果节点分布非常均匀（间距相近），影响小
2. 如果极耳很小（角度范围小），半径变化小
3. 简化计算，避免几何复杂性

**反对**：
1. 卷绕结构中，半径变化大（Rin=2mm, Rout=10mm）
2. 节点间距差异显著（可达5倍）
3. 物理机制不正确

### 建议

**推荐：使用弧长权重** ✅

理由：
1. **物理正确**：节点代表的弧长不同，应按弧长分配
2. **数值合理**：外层散热强，内层散热弱，符合实际
3. **计算简单**：用直线距离近似弧长（以直代曲）
4. **影响显著**：系数差异可达数倍

## 实现方案

### 计算节点弧长

```julia
"""
计算极耳节点代表的弧长（以直代曲）

方法：
对于节点i，其代表的弧长为与相邻节点的距离平均。

返回：
arc_lengths[i] = 节点 i 代表的弧长 [m]
"""
function _compute_tab_node_arc_lengths(mesh, tab_nodes)
    n_nodes = length(tab_nodes)
    arc_lengths = zeros(Float64, n_nodes)
    
    if n_nodes == 0
        return arc_lengths
    elseif n_nodes == 1
        # 单个节点：假设代表一个最小弧长
        arc_lengths[1] = 1.0  # 相对权重
        return arc_lengths
    end
    
    # 节点坐标
    coords = [mesh.node[n, :] for n in tab_nodes]
    
    # 计算相邻节点间距
    for i in 1:n_nodes
        if i == 1
            # 第一个节点：到下一个节点的距离
            arc_lengths[i] = norm(coords[2] - coords[1])
        elseif i == n_nodes
            # 最后一个节点：到前一个节点的距离
            arc_lengths[i] = norm(coords[i] - coords[i-1])
        else
            # 中间节点：前后距离的平均
            dist_prev = norm(coords[i] - coords[i-1])
            dist_next = norm(coords[i+1] - coords[i])
            arc_lengths[i] = (dist_prev + dist_next) / 2.0
        end
    end
    
    return arc_lengths
end
```

### 应用弧长权重

```julia
function _apply_cool_tab!(KT, FT, mesh, case, t)
    # 识别极耳节点
    pos_idx, neg_idx = jellyroll_tab_node_indices(mesh, case.param_dim)
    tab_nodes = unique(vcat(pos_idx, neg_idx))
    
    isempty(tab_nodes) && return
    
    # 参数
    h_tab = case.opt.h_tab
    tab_area = case.param_dim.tab.area
    H = case.param_dim.cell.width
    k_th, L_th = case.param_dim.scale.k_th, case.param_dim.scale.L_th
    T_amb_nd = case.param_dim.cell.T_amb / case.param_dim.scale.T_ref
    
    # 计算节点弧长（以直代曲）
    arc_lengths = _compute_tab_node_arc_lengths(mesh, tab_nodes)
    total_arc_length = sum(arc_lengths)
    
    if total_arc_length < 1e-12
        @warn "极耳总弧长过小" total_arc_length=total_arc_length
        return
    end
    
    # 按弧长权重分配
    for (i, n) in enumerate(tab_nodes)
        weight = arc_lengths[i] / total_arc_length
        coeff = h_tab * tab_area * weight / (H * k_th * L_th)
        
        KT[n, n] += coeff
        FT[n] += coeff * T_amb_nd
    end
end
```

## 验证

### 权重归一化

```julia
@assert abs(sum(arc_lengths) - total_arc_length) < 1e-10
@assert abs(sum(weights) - 1.0) < 1e-10
```

### 总弧长检查

```julia
# 总弧长应接近 tab.width
@assert abs(total_arc_length - tab.width) / tab.width < 0.5
```

### 能量守恒

```julia
Q_total = sum([coeff[i] * (T[n] - T_amb) for (i,n) in enumerate(tab_nodes)])
Q_expected = h_tab * tab_area * (T_avg - T_amb)
@assert abs(Q_total - Q_expected) / Q_expected < 0.01
```

## 总结

### 均匀分配的问题

1. ❌ **物理不正确**：忽略节点间距差异
2. ❌ **数值误差大**：系数差异可达5倍
3. ❌ **温度场失真**：无法正确反映外层散热强于内层

### 弧长权重的优势

1. ✅ **物理正确**：每个节点按其代表的弧长分配
2. ✅ **数值合理**：外层节点散热强，内层散热弱
3. ✅ **计算简单**：用直线距离近似（以直代曲）
4. ✅ **影响显著**：对温度场和收敛性有明显改善

### 建议

**必须使用弧长权重分配** ✅

理由：
- 物理机制要求
- 数值精度要求
- 对后续组装有显著影响

**不推荐均匀分配** ❌

除非：
- 极耳角度范围极小（<1°）
- 半径变化极小（Rout/Rin < 1.2）
- 仅用于定性分析

---

**结论**：弧长权重分配是**必要**的，不是过度优化！
