# Z方向冷却对XY平面2D模型的贡献机制

## 核心问题

2D模型是xy平面的俯视图，热传导方程只涉及xy方向：
$$
\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x}\left(k \frac{\partial T}{\partial x}\right) + \frac{\partial}{\partial y}\left(k \frac{\partial T}{\partial y}\right) + q
$$

**问题**：z方向的对流散热如何影响这个2D方程的刚度矩阵和载荷向量？

## 物理机制：体积热汇

### 1. 三维真实情况

考虑一个微元体 $dV = dA \times H$，其中：
- $dA$：xy平面微元面积
- $H$：电池高度（z方向）

上下表面的对流散热功率：
$$
\dot{Q}_{\text{conv}} = h_{\text{top}}(T - T_{\text{amb}}) \cdot dA + h_{\text{bottom}}(T - T_{\text{amb}}) \cdot dA
$$

假设上下表面换热系数相同 $h_{\text{top}} = h_{\text{bottom}} = h$：
$$
\dot{Q}_{\text{conv}} = 2h(T - T_{\text{amb}}) \cdot dA
$$

### 2. 转换为体积热汇

单位体积的散热率（负热源）：
$$
q_{\text{vol}} = \frac{\dot{Q}_{\text{conv}}}{dV} = \frac{2h(T - T_{\text{amb}}) \cdot dA}{dA \times H} = \frac{2h}{H}(T - T_{\text{amb}})
$$

### 3. 修正的2D热传导方程

将z方向对流散热作为体积热汇项加入：
$$
\rho c_p \frac{\partial T}{\partial t} = \nabla_{xy} \cdot (k \nabla_{xy} T) + q_{\text{gen}} - \frac{2h}{H}(T - T_{\text{amb}})
$$

其中 $\nabla_{xy}$ 表示只在xy平面的梯度算子。

## 弱形式推导

### 1. 乘以权函数并积分

$$
\int_{\Omega_{xy}} \rho c_p \frac{\partial T}{\partial t} N_i \, d\Omega + \int_{\Omega_{xy}} k \nabla_{xy}T \cdot \nabla_{xy}N_i \, d\Omega + \int_{\Omega_{xy}} \frac{2h}{H}(T - T_{\text{amb}}) N_i \, d\Omega = \int_{\Omega_{xy}} q_{\text{gen}} N_i \, d\Omega
$$

### 2. 对流项展开

$$
\int_{\Omega_{xy}} \frac{2h}{H}T N_i \, d\Omega - \int_{\Omega_{xy}} \frac{2h}{H}T_{\text{amb}} N_i \, d\Omega
$$

### 3. 离散化

温度场离散为 $T = \sum_j N_j T_j$，代入对流项：

**刚度矩阵贡献**：
$$
K_{ij}^{\text{conv}} = \int_{\Omega_{xy}} \frac{2h}{H} N_i N_j \, d\Omega
$$

**载荷向量贡献**：
$$
F_i^{\text{conv}} = \int_{\Omega_{xy}} \frac{2h}{H} T_{\text{amb}} N_i \, d\Omega
$$

## 无量纲化

### 1. 特征量

- 温度尺度：$T_{\text{ref}}$
- 长度尺度：$L_{\text{th}}$（通常取 $R_{\text{out}}$）
- 热导率尺度：$k_{\text{th}}$

### 2. 无量纲参数

定义修正的Biot数：
$$
Bi_z = \frac{2h L_{\text{th}}^2}{H k_{\text{th}}}
$$

物理意义：
- 分子 $2hL_{\text{th}}^2$：对流散热能力
- 分母 $Hk_{\text{th}}$：体积热传导能力
- $Bi_z$ 衡量z方向散热与xy平面传导的相对强度

### 3. 无量纲刚度矩阵和载荷

$$
K_{ij}^{*, \text{conv}} = \int_{\Omega_{xy}^*} Bi_z \cdot N_i N_j \, d\Omega^*
$$

$$
F_i^{*, \text{conv}} = \int_{\Omega_{xy}^*} Bi_z \cdot T_{\text{amb}}^* N_i \, d\Omega^*
$$

其中 $\Omega_{xy}^* = \Omega_{xy} / L_{\text{th}}^2$，$T^* = T / T_{\text{ref}}$。

## 两种冷却方式的差异

### 表面冷却（surface）

**适用范围**：整个xy平面
$$
K_{ij}^{\text{surface}} = \int_{\Omega_{xy}} \frac{2h_{\text{surface}}}{H} N_i N_j \, d\Omega
$$

**实现**：对所有单元进行高斯积分求和

### 极耳冷却（tab）

**适用范围**：仅极耳节点邻域
$$
K_{ij}^{\text{tab}} = \int_{\Omega_{\text{tab}}} \frac{2h_{\text{tab}}}{H} N_i N_j \, d\Omega
$$

其中 $\Omega_{\text{tab}}$ 是极耳节点所在的单元区域。

**实现方法1**（精确）：
对包含极耳节点的单元进行高斯积分

**实现方法2**（简化）：
假设极耳节点的"影响面积"为 $A_{\text{node}}$（通过Voronoi分割或单元面积平均）：
$$
K_{ii}^{\text{tab}} \approx \frac{2h_{\text{tab}}}{H} \cdot A_{\text{node}}
$$

## 关键洞察

### 1. 不是边界积分！

z方向对流**不是**在xy平面边界上的线积分，而是在xy平面域内的**面积分**（体积热汇）。

### 2. 对角占优

对流项只增加刚度矩阵的对角元素（如果采用节点集中方式）或近对角元素（如果采用高斯积分）。

### 3. 数值稳定性

$$
K_{ii}^{\text{conv}} \sim \frac{2h}{H} \cdot A_i \sim O(10^{-3}) \text{ to } O(1)
$$

远小于惩罚法的 $10^{12}$，因此数值稳定。

## 物理解释

可以将z方向对流散热理解为：

1. **热阻网络**：
   - xy平面传导：热阻 $\sim L/(kA)$
   - z方向对流：热阻 $\sim 1/(hA)$
   - 并联工作

2. **等效体积热源**：
   - 每个xy位置都有一个"竖直的冷却管道"
   - 冷却强度 $\propto 2h/H$

3. **极耳强化**：
   - 极耳位置的 $h$ 更大（$h_{\text{tab}} \gg h_{\text{surface}}$）
   - 相当于"加粗的冷却管道"

## 数值实现细节

### 方案A：单元高斯积分（精确）

```julia
for each element e
    for each gauss point g in e
        wt = (2h/H) * weight[g] * detJ[g] / L_th^2
        for i in nodes_of_e
            for j in nodes_of_e
                K[i,j] += wt * N_i(g) * N_j(g)
            end
            F[i] += wt * T_amb * N_i(g)
        end
    end
end
```

### 方案B：节点面积法（简化，仅用于极耳节点）

```julia
for each tab node n
    A_node = compute_node_area(n)  # Voronoi或单元平均
    K[n,n] += (2h_tab/H) * A_node / L_th^2
    F[n] += (2h_tab/H) * T_amb * A_node / L_th^2
end
```

## 验证方法

### 1. 能量守恒

$$
\frac{dE}{dt} = \int_{\Omega_{xy}} q_{\text{gen}} \, dA - \int_{\Omega_{xy}} \frac{2h}{H}(T - T_{\text{amb}}) \, dA
$$

左边：储能变化率
右边：产热 - z方向散热

### 2. 温度分布合理性

- 表面冷却：整体温度降低，分布形状类似
- 极耳冷却：极耳邻域温度降低更明显，形成"冷斑"

### 3. Biot数检查

$$
Bi_z = \frac{2h L_{\text{th}}^2}{H k_{\text{th}}} \sim O(10^{-3}) \text{ to } O(1)
$$

如果 $Bi_z \ll 1$：z方向散热可忽略
如果 $Bi_z \sim O(1)$：z方向散热与xy传导同等重要
如果 $Bi_z \gg 1$：z方向散热主导（不太可能）

## 总结

**关键结论**：

1. z方向对流散热以**体积热汇**的形式进入2D方程
2. 贡献系数：$2h/H$（单位：W/(m³·K)）
3. 刚度矩阵：$K_{ij} += \int \frac{2h}{H} N_i N_j \, dA$（面积分，不是线积分！）
4. 物理尺度：$K^* \sim Bi_z = \frac{2hL^2}{Hk} \sim O(10^{-3})$ to $O(1)$
5. 数值稳定，无病态问题

**与原先理解的差异**：

- ❌ 不是在xy平面边界上的线积分
- ✅ 是在xy平面域内的面积分
- ❌ 不需要识别"边界边"
- ✅ 对所有节点（或极耳节点）施加体积散热
