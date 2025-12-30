# 极耳边界对流散热理论与实现

## 1. 物理模型

### 1.1 两种边界处理方式

极耳边界有两种冷却策略，都模拟电池在z方向（上下表面）的散热：

#### 方式1：一般表面散热
- **适用场景**：整个电池上下表面均匀冷却
- **物理意义**：电池通过上下表面与环境换热
- **特点**：所有节点使用相同的换热系数

#### 方式2：极耳强化散热
- **适用场景**：极耳区域有额外的冷却措施（如极耳与外壳接触）
- **物理意义**：极耳区域散热能力强于其他区域
- **特点**：极耳节点使用更大的换热系数

## 2. 理论推导

### 2.1 三维热传导方程

$$
\rho c_p \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + q
$$

### 2.2 对流边界条件

在边界 $\Gamma$ 上，对流换热：

$$
-k \frac{\partial T}{\partial n} = h(T - T_{\text{amb}})
$$

其中：
- $h$: 对流换热系数 [W/(m²·K)]
- $T_{\text{amb}}$: 环境温度 [K]
- $n$: 边界外法向

### 2.3 弱形式

将对流边界条件代入弱形式：

$$
\int_{\Omega} \rho c_p \frac{\partial T}{\partial t} N_i \, d\Omega + \int_{\Omega} k \nabla T \cdot \nabla N_i \, d\Omega + \int_{\Gamma} h(T - T_{\text{amb}}) N_i \, d\Gamma = \int_{\Omega} q N_i \, d\Omega
$$

边界项离散：

$$
\int_{\Gamma} h T N_i \, d\Gamma - \int_{\Gamma} h T_{\text{amb}} N_i \, d\Gamma
$$

对应矩阵形式：
- 刚度矩阵贡献：$K_{ij}^{\text{bc}} = \int_{\Gamma} h N_i N_j \, d\Gamma$
- 载荷向量贡献：$F_i^{\text{bc}} = \int_{\Gamma} h T_{\text{amb}} N_i \, d\Gamma$

## 3. 2D模型模拟z方向冷却

### 3.1 问题描述

Jellyroll 2D模型是圆柱电池的**俯视图**（x-y平面），忽略z方向变化。但物理上，电池通过**上下表面**（z方向）与环境换热。

### 3.2 等效处理

将z方向的表面积分"投影"到2D节点上：

**方法1：节点投影面积法**

每个节点代表的z方向投影面积：

$$
A_{\text{z}}(i) = A_{\text{voronoi}}(i) \times H
$$

其中：
- $A_{\text{voronoi}}(i)$: 节点在x-y平面的Voronoi面积
- $H$: 电池高度（z方向）

等效对流项：

$$
K_{ii}^{\text{conv,z}} = h \cdot A_{\text{z}}(i)
$$

$$
F_i^{\text{conv,z}} = h \cdot T_{\text{amb}} \cdot A_{\text{z}}(i)
$$

**方法2：单元积分法（推荐）**

对每个单元，将z方向表面积分分配到节点：

$$
K_{ij}^{\text{conv,z}} = \int_{\Omega_e} h \cdot \frac{H}{L_{\text{th}}^2} \cdot N_i N_j \, dx dy
$$

其中 $\frac{H}{L_{\text{th}}^2}$ 是将z方向面密度转换为2D积分的尺度因子。

### 3.3 无量纲化

使用 Scheme B 热尺度：

- 温度：$T^* = T / T_{\text{ref}}$
- Biot数：$Bi_z = h H / k_{\text{th}}$
- 特征长度：$L_{\text{th}}$ (通常取 $R_{\text{out}}$)

无量纲对流项：

$$
K_{ij}^{*, \text{conv,z}} = \int_{\Omega_e^*} Bi_z \cdot \frac{H}{L_{\text{th}}} \cdot N_i N_j \, dx^* dy^*
$$

## 4. 两种实现方式对比

| 特性 | 方式1：一般表面散热 | 方式2：极耳强化散热 |
|------|-------------------|-------------------|
| 作用范围 | 所有节点 | 仅极耳节点 |
| 换热系数 | $h_{\text{cell}}$ | $h_{\text{tab}}$ (通常更大) |
| 物理意义 | 整体冷却 | 局部强化冷却 |
| 实现位置 | 全局循环 | 极耳节点循环 |
| 典型值 | $h \sim 10$ W/(m²·K) | $h_{\text{tab}} \sim 100$ W/(m²·K) |

## 5. 数值稳定性

### 5.1 与惩罚法对比

| 方法 | 对角元素增量 | 条件数 | 稳定性 |
|------|------------|--------|--------|
| 惩罚法 | $10^{12}$ | $10^{12}$ | ❌ 不稳定 |
| 对流法 | $h A \sim 10^{-3}$ | $O(1)$ | ✅ 稳定 |

### 5.2 时间步长限制

对流边界条件不会引入额外的时间步长限制，因为：

$$
\Delta t < \frac{\rho c_p L^2}{k}
$$

对流项 $h$ 的量级远小于扩散项 $k/L$。

## 6. 参数设置指南

### 6.1 一般表面散热

```julia
opt.tab_bc_type = "surface_convection"  # 选择方式1
opt.h_surface = 10.0  # W/(m²·K)，典型值：5-50
```

### 6.2 极耳强化散热

```julia
opt.tab_bc_type = "tab_convection"  # 选择方式2
opt.h_tab = 100.0  # W/(m²·K)，典型值：50-500
```

### 6.3 惩罚法（不推荐，仅用于固定温度边界）

```julia
opt.tab_bc_type = "penalty"  # 传统方式
opt.tab_penalty = 1e6  # 降低到1e6-1e8
```

## 7. 验证方法

### 7.1 能量守恒检查

总热流应等于产热减去储热：

$$
Q_{\text{gen}} = Q_{\text{storage}} + Q_{\text{conv,outer}} + Q_{\text{conv,z}}
$$

### 7.2 温度场合理性

- 极耳区域温度应略低于内部（如果 $h_{\text{tab}}$ 较大）
- 温度分布应连续光滑
- 不应出现非物理的温度突变

## 8. 参考文献

1. Heat Transfer Textbook, Section on Convection Boundary Conditions
2. Finite Element Heat Transfer Analysis, Chapter on Natural Boundary Conditions
3. Battery Thermal Management Review Papers
