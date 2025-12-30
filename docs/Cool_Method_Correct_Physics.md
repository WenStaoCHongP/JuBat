# Z方向冷却的正确物理机制与实现

## 核心认识

### 关键问题
2D Jellyroll模型是xy平面的俯视图，只考虑xy方向的热传导。那么**z方向（上下表面）的对流散热如何影响xy平面的刚度矩阵和载荷向量？**

### 答案：体积热汇（Volume Heat Sink）

z方向的对流散热**不是**在xy平面边界上的线积分，而是以**体积热汇**的形式分布在整个xy域内。

## 物理推导

### 1. 三维真实情况

考虑一个微元体 $dV = dA \times H$：
- $dA$：xy平面微元面积
- $H$：电池高度（z方向）

上下表面的对流散热功率：
$$
\dot{Q}_{\text{conv}} = [h_{\text{top}} + h_{\text{bottom}}](T - T_{\text{amb}}) \cdot dA \approx 2h(T - T_{\text{amb}}) \cdot dA
$$

### 2. 单位体积散热率

$$
q_{\text{vol}} = \frac{\dot{Q}_{\text{conv}}}{dV} = \frac{2h(T - T_{\text{amb}}) \cdot dA}{dA \times H} = \frac{2h}{H}(T - T_{\text{amb}})
$$

**关键系数**：$\alpha = \frac{2h}{H}$ [W/(m³·K)]

### 3. 修正的2D热传导方程

$$
\rho c_p \frac{\partial T}{\partial t} = \nabla_{xy} \cdot (k \nabla_{xy} T) + q_{\text{gen}} - \frac{2h}{H}(T - T_{\text{amb}})
$$

体积热汇项：$-\frac{2h}{H}(T - T_{\text{amb}})$

## 弱形式与FEM离散

### 1. 乘以权函数并在xy域积分

$$
\int_{\Omega_{xy}} \rho c_p \frac{\partial T}{\partial t} N_i \, d\Omega + \int_{\Omega_{xy}} k \nabla T \cdot \nabla N_i \, d\Omega + \int_{\Omega_{xy}} \frac{2h}{H}(T - T_{\text{amb}}) N_i \, d\Omega = \int_{\Omega_{xy}} q_{\text{gen}} N_i \, d\Omega
$$

注意：**第三项是xy平面上的面积分，不是边界线积分！**

### 2. 离散化

温度场 $T = \sum_j N_j T_j$，代入对流项：

$$
\int_{\Omega_{xy}} \frac{2h}{H} \sum_j N_j T_j N_i \, d\Omega - \int_{\Omega_{xy}} \frac{2h}{H} T_{\text{amb}} N_i \, d\Omega
$$

**刚度矩阵贡献**：
$$
K_{ij}^{\text{conv}} = \int_{\Omega_{xy}} \frac{2h}{H} N_i N_j \, d\Omega
$$

**载荷向量贡献**：
$$
F_i^{\text{conv}} = \int_{\Omega_{xy}} \frac{2h}{H} T_{\text{amb}} N_i \, d\Omega
$$

## 无量纲化

### 特征尺度
- 温度：$T_{\text{ref}} = 298$ K
- 长度：$L_{\text{th}} = R_{\text{out}}$
- 热导率：$k_{\text{th}}$

### 修正的Biot数

$$
Bi_z = \frac{2h L_{\text{th}}^2}{H k_{\text{th}}}
$$

**物理意义**：
- 分子：z方向对流散热能力（按单位体积）
- 分母：xy平面传导能力（按单位体积）

### 无量纲刚度矩阵

$$
K_{ij}^{*, \text{conv}} = \int_{\Omega_{xy}^*} Bi_z \cdot N_i N_j \, d\Omega^*
$$

其中 $\Omega_{xy}^* = \Omega_{xy} / L_{\text{th}}^2$

## 实现方式对比

### 方式1：整体表面冷却（surface）

**适用范围**：整个xy域 $\Omega_{xy}$

**积分**：
$$
K_{ij}^{\text{surface}} = \int_{\Omega_{xy}} \frac{2h_{\text{surface}}}{H} N_i N_j \, d\Omega
$$

**实现**：对所有单元进行高斯积分求和

```julia
for g in 1:ngs  # 所有高斯点
    e = ele[g]
    wt = (2h/H) * wJ[g] / L_th^2
    
    for i, j in nodes_of_element(e)
        K[i,j] += wt * N_i(g) * N_j(g)
    end
    
    for i in nodes_of_element(e)
        F[i] += wt * T_amb * N_i(g)
    end
end
```

### 方式2：极耳强化冷却（tab）

**适用范围**：仅极耳节点邻域

**节点识别**：
通过 `jellyroll_tab_node_indices` 识别螺旋线上的离散节点（以直代曲）

**简化方法**：
假设极耳节点的"影响面积"为 $A_{\text{node}}$（通过单元面积平均分配）

$$
K_{ii}^{\text{tab}} = \frac{2h_{\text{tab}}}{H} \cdot A_{\text{node}}
$$

$$
F_i^{\text{tab}} = \frac{2h_{\text{tab}}}{H} \cdot T_{\text{amb}} \cdot A_{\text{node}}
$$

**实现**：
```julia
# 计算节点面积（以直代曲）
node_areas = compute_node_areas(mesh)

for n in tab_nodes
    A_node = node_areas[n]
    coeff = (2h_tab/H) * A_node / L_th^2
    
    K[n,n] += coeff
    F[n] += coeff * T_amb
end
```

## 关键洞察

### 1. 不是边界积分！

❌ **错误理解**：z方向对流在xy平面边界上的线积分
$$
\int_{\partial\Omega_{xy}} h(T - T_{\text{amb}}) N_i \, ds
$$

✅ **正确理解**：z方向对流在xy平面域内的面积分（体积热汇）
$$
\int_{\Omega_{xy}} \frac{2h}{H}(T - T_{\text{amb}}) N_i \, d\Omega
$$

### 2. 极耳节点是"螺旋线离散点"

通过 `jellyroll_tab_node_indices` 识别的极耳节点是：
- 螺旋线上的离散采样点
- 使用"以直代曲"思想
- 每个节点代表一段螺旋弧的影响范围

### 3. 数值尺度

$$
K^{\text{conv}} \sim \frac{2h L_{\text{th}}^2}{H k_{\text{th}}} \cdot A_{\text{node}} \sim Bi_z \cdot O(10^{-3}) \sim O(10^{-3}) \text{ to } O(1)
$$

远小于惩罚法的 $10^{12}$，因此数值稳定。

## 物理解释

### 热阻网络类比

将每个xy位置看作一个"热节点"，它有：

1. **xy方向传导**：与周围节点通过热阻连接
   $$R_{\text{cond}} \sim \frac{L}{kA}$$

2. **z方向对流**：与环境通过热阻连接
   $$R_{\text{conv}} \sim \frac{H}{2hA}$$

3. **并联工作**：热量可以通过xy传导或z散热

### 极耳强化的意义

- 普通位置：$h_{\text{surface}} \sim 10$ W/(m²·K)
- 极耳位置：$h_{\text{tab}} \sim 100$ W/(m²·K)

极耳位置的对流热阻降低10倍，形成"冷却通道"。

## 验证方法

### 1. 能量守恒

$$
\frac{dE}{dt} = \int_{\Omega_{xy}} q_{\text{gen}} \, dA - \int_{\Omega_{xy}} \frac{2h}{H}(T - T_{\text{amb}}) \, dA
$$

储能变化 = 产热 - z方向散热

### 2. Biot数检查

$$
Bi_z = \frac{2h L_{\text{th}}^2}{H k_{\text{th}}}
$$

典型值：
- 自然对流：$Bi_z \sim 10^{-3}$ to $10^{-2}$
- 强制对流：$Bi_z \sim 10^{-2}$ to $10^{-1}$
- 液冷接触：$Bi_z \sim 10^{-1}$ to $1$

### 3. 温度分布

**表面冷却**：整体温度下降，分布形状保持类似

**极耳冷却**：极耳邻域形成"冷斑"，温度梯度增大

## 数值实现细节

### 方案A：高斯积分（精确，用于surface）

```julia
function _apply_cool_surface!(KT, FT, mesh, case, t)
    vol_coeff = 2h_surface / H
    Bi_z = vol_coeff * L_th^2 / k_th
    conv_factor = Bi_z / L_th^2
    
    for g in 1:ngs
        e = ele[g]
        wt = conv_factor * wJ[g]
        
        for i, j in element_nodes(e)
            KT[i,j] += wt * Ni[g,i] * Ni[g,j]
        end
        
        for i in element_nodes(e)
            FT[i] += wt * T_amb * Ni[g,i]
        end
    end
end
```

### 方案B：节点面积法（简化，用于tab）

```julia
function _apply_cool_tab!(KT, FT, mesh, case, t)
    vol_coeff = 2h_tab / H
    Bi_z = vol_coeff * L_th^2 / k_th
    
    node_areas = compute_node_areas(mesh)
    
    for n in tab_nodes
        A_nd = node_areas[n] / L_th^2
        
        KT[n,n] += Bi_z * A_nd
        FT[n] += Bi_z * T_amb * A_nd
    end
end
```

## 总结

### 核心要点

1. ✅ **z方向冷却 = 体积热汇**：$q_{\text{vol}} = \frac{2h}{H}(T - T_{\text{amb}})$

2. ✅ **xy平面面积分**：$K_{ij} += \int_{\Omega_{xy}} \frac{2h}{H} N_i N_j \, dA$

3. ✅ **不是边界线积分**：这不是在xy平面边界上的积分

4. ✅ **极耳节点 = 螺旋线离散点**：以直代曲，每个节点代表一段弧

5. ✅ **数值稳定**：$Bi_z \sim O(10^{-3})$ to $O(1)$，远小于惩罚法

### 使用方法

**整体表面冷却**：
```julia
opt.cool_method = "surface"
opt.h_surface = 10.0  # W/(m²·K)
```

**极耳强化冷却**：
```julia
opt.cool_method = "tab"
opt.h_tab = 100.0  # W/(m²·K)
```

### 物理图景

想象电池俯视图（xy平面），每个位置都有一个"竖直的冷却管道"通往z方向的上下表面：

- 普通位置：细管道（$h_{\text{surface}}$ 小）
- 极耳位置：粗管道（$h_{\text{tab}}$ 大）

冷却强度 $\propto 2h/H$，以体积热汇形式分布在xy域内。
