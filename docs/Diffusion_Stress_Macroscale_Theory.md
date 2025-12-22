# 宏观层面扩散应力计算理论

## 1. 引言

本文档描述如何在宏观电池单元层面计算锂浓度变化引起的扩散应力，参照二维热应力的计算方法。该方法适用于果冻卷电池的 2D 分布式热-力耦合模型。

## 2. 物理模型

### 2.1 锂浓度与SOC的关系

在宏观单元层面，将锂浓度变化等效为该单元的SOC变化：

```
SOC = (c_s - c_s_min) / (c_s_max - c_s_min)
```

其中：
- `c_s`: 固相锂浓度的平均值
- `c_s_max`: 最大锂浓度
- `c_s_min`: 最小锂浓度

### 2.2 扩散诱导的体积应变

锂嵌入/脱出会引起电极材料的体积变化，定义为：

```
ε_vol = β_c * Δc_s = β_c * c_s_max * ΔSOC
```

其中：
- `β_c`: 锂浓度变化引起的体积膨胀系数 (类似热膨胀系数)
- `Δc_s`: 锂浓度变化量
- `ΔSOC`: SOC 变化量

对于各向同性材料，线性应变分量：

```
ε_xx = ε_yy = ε_zz = ε_vol / 3 = (β_c * c_s_max * ΔSOC) / 3
```

### 2.3 总应变的分解

总应变可以分解为：
1. 弹性应变 `ε_e`
2. 热应变 `ε_th = α * ΔT`
3. 扩散应变 `ε_diff = β_c * c_s_max * ΔSOC`

```
ε_total = ε_e + ε_th + ε_diff
```

## 3. 控制方程

### 3.1 应力平衡方程

在准静态条件下（忽略惯性项），应力平衡方程为：

```
∂σ_xx/∂x + ∂σ_xy/∂y = 0
∂σ_xy/∂x + ∂σ_yy/∂y = 0
```

边界条件：
- 内边界：绝热边界，法向应力为零 `σ_n = 0`
- 外边界：自由边界 `σ_n = 0`，或约束边界 `u = 0`

### 3.2 变形协调方程

应变-位移关系（小变形假设）：

```
ε_xx = ∂u/∂x
ε_yy = ∂v/∂y
ε_xy = (1/2) * (∂u/∂y + ∂v/∂x)
```

其中 `u`, `v` 是 x, y 方向的位移。

### 3.3 本构关系

对于平面应力问题（2D电池截面，垂直方向应力 σ_zz ≈ 0）：

弹性应变与应力的关系：

```
ε_e_xx = (1/E) * (σ_xx - ν*σ_yy) + α*ΔT + β_c*c_s_max*ΔSOC
ε_e_yy = (1/E) * (σ_yy - ν*σ_xx) + α*ΔT + β_c*c_s_max*ΔSOC
ε_e_xy = (1+ν)/E * σ_xy
```

反过来，应力-应变关系：

```
σ_xx = E/(1-ν²) * [(ε_xx - α*ΔT - β_c*c_s_max*ΔSOC) + ν*(ε_yy - α*ΔT - β_c*c_s_max*ΔSOC)]
σ_yy = E/(1-ν²) * [(ε_yy - α*ΔT - β_c*c_s_max*ΔSOC) + ν*(ε_xx - α*ΔT - β_c*c_s_max*ΔSOC)]
σ_xy = E/(2*(1+ν)) * ε_xy
```

简化后：

```
σ_xx = E/(1-ν²) * [(ε_xx + ν*ε_yy) - (1+ν)*(α*ΔT + β_c*c_s_max*ΔSOC)]
σ_yy = E/(1-ν²) * [(ε_yy + ν*ε_xx) - (1+ν)*(α*ΔT + β_c*c_s_max*ΔSOC)]
σ_xy = E/(2*(1+ν)) * ε_xy
```

### 3.4 弱形式（虚功原理）

将应力平衡方程乘以虚位移 `δu`, `δv` 并在域 Ω 上积分：

```
∫_Ω [σ_xx * ∂(δu)/∂x + σ_yy * ∂(δv)/∂y + σ_xy * (∂(δu)/∂y + ∂(δv)/∂x)] dΩ = 0
```

引入应力-应变关系后，得到以位移 `u`, `v` 为未知量的弱形式方程。

## 4. 有限元离散

### 4.1 位移插值

对于 Q4 单元（四节点四边形单元），位移场采用双线性插值：

```
u(x,y) = Σ N_i(x,y) * u_i
v(x,y) = Σ N_i(x,y) * v_i
```

其中 `N_i` 是形函数，`u_i`, `v_i` 是节点位移。

### 4.2 应变-位移矩阵

应变可以表示为：

```
{ε} = [B] {d}
```

其中：
- `{ε} = [ε_xx, ε_yy, γ_xy]^T`
- `{d} = [u_1, v_1, u_2, v_2, ..., u_4, v_4]^T` 是单元节点位移向量
- `[B]` 是应变-位移矩阵

对于节点 i：

```
[B_i] = [
  ∂N_i/∂x    0
  0          ∂N_i/∂y
  ∂N_i/∂y    ∂N_i/∂x
]
```

### 4.3 本构矩阵

对于平面应力问题：

```
{σ} = [D] * ({ε} - {ε_0})
```

其中弹性矩阵：

```
[D] = E/(1-ν²) * [
  1    ν    0
  ν    1    0
  0    0    (1-ν)/2
]
```

初始应变（热应变 + 扩散应变）：

```
{ε_0} = [
  α*ΔT + β_c*c_s_max*ΔSOC
  α*ΔT + β_c*c_s_max*ΔSOC
  0
]
```

### 4.4 单元刚度矩阵

单元刚度矩阵：

```
[K_e] = ∫_Ω_e [B]^T [D] [B] dΩ
```

使用高斯积分：

```
[K_e] = Σ_g w_g * |J_g| * [B_g]^T [D] [B_g]
```

### 4.5 等效节点力（初始应变引起）

热应变和扩散应变引起的等效节点力：

```
{F_e} = ∫_Ω_e [B]^T [D] {ε_0} dΩ
```

使用高斯积分：

```
{F_e} = Σ_g w_g * |J_g| * [B_g]^T [D] {ε_0,g}
```

其中：

```
{ε_0,g} = [
  α*ΔT_g + β_c*c_s_max*ΔSOC_g
  α*ΔT_g + β_c*c_s_max*ΔSOC_g
  0
]
```

### 4.6 全局方程组

组装全局刚度矩阵 `[K]` 和载荷向量 `{F}`：

```
[K] {U} = {F}
```

其中：
- `{U}` 是全局节点位移向量
- `{F}` 是由热应变和扩散应变引起的等效节点力

### 4.7 边界条件

1. **固定约束**：某些节点位移为零（如对称轴）
   ```
   u_i = 0  或  v_i = 0
   ```

2. **自由边界**：不施加约束

3. **对称边界**：法向位移为零，切向自由
   ```
   u_n = 0,  τ_t = 0
   ```

## 5. 求解流程

### 5.1 输入数据

对于每个时间步：
1. 温度场 `T(x,y)` 从热传导方程求解得到
2. SOC分布 `SOC(x,y)` 从电化学模型计算
3. 材料参数：`E`, `ν`, `α`, `β_c`, `c_s_max`

### 5.2 计算步骤

1. **计算初始应变**：
   ```julia
   ΔT_elem = T_elem - T_ref
   ΔSOC_elem = SOC_elem - SOC_ref
   ε_0_elem = α * ΔT_elem + β_c * c_s_max * ΔSOC_elem
   ```

2. **装配刚度矩阵**：
   ```julia
   K = assemble_stiffness_matrix(mesh, material_props)
   ```

3. **装配载荷向量**：
   ```julia
   F = assemble_thermal_diffusion_load(mesh, ε_0_elem, material_props)
   ```

4. **施加边界条件**：
   ```julia
   K, F = apply_boundary_conditions(K, F, bc)
   ```

5. **求解位移**：
   ```julia
   U = K \ F
   ```

6. **恢复应力**：
   ```julia
   σ = recover_stress(U, ε_0_elem, mesh, material_props)
   ```

### 5.3 应力恢复

对于每个单元，在高斯点处：

```julia
ε_g = B_g * d_e              # 总应变
ε_elastic = ε_g - ε_0_g      # 弹性应变
σ_g = D * ε_elastic          # 应力
```

### 5.4 应力分量

最终得到的应力包括：
- 正应力：`σ_xx`, `σ_yy`
- 剪应力：`σ_xy`
- Von Mises等效应力：
  ```
  σ_vm = √(σ_xx² + σ_yy² - σ_xx*σ_yy + 3*σ_xy²)
  ```

## 6. 与热应力的类比

| 项目 | 热应力 | 扩散应力 |
|------|--------|----------|
| 驱动量 | 温度变化 ΔT | SOC变化 ΔSOC |
| 膨胀系数 | 热膨胀系数 α | 浓度膨胀系数 β_c |
| 初始应变 | ε_th = α*ΔT | ε_diff = β_c*c_s_max*ΔSOC |
| 本构关系 | σ = E/(1-ν)*(ε - α*ΔT) | σ = E/(1-ν)*(ε - β_c*c_s_max*ΔSOC) |
| 控制方程 | 应力平衡 ∇·σ = 0 | 应力平衡 ∇·σ = 0 |

## 7. 热-扩散应力耦合

当同时考虑热应力和扩散应力时，总的初始应变为：

```
{ε_0} = [
  α*ΔT + β_c*c_s_max*ΔSOC
  α*ΔT + β_c*c_s_max*ΔSOC
  0
]
```

这样就自然地将热应力和扩散应力统一在同一个有限元框架内。

## 8. 代码实现要点

### 8.1 主要函数

```julia
# 1. 装配力学刚度矩阵
function assemble_mechanical_stiffness_2D(mesh, E, ν)
    # 返回全局刚度矩阵 K
end

# 2. 计算热-扩散载荷
function assemble_thermal_diffusion_load_2D(mesh, T_elem, SOC_elem, α, β_c, c_s_max, E, ν)
    # 返回等效节点力 F
end

# 3. 求解位移
function solve_mechanical_displacement(K, F, bc)
    # 返回节点位移 U
end

# 4. 恢复应力
function recover_stress_2D(U, mesh, T_elem, SOC_elem, α, β_c, c_s_max, E, ν)
    # 返回单元应力 σ_elem
end

# 5. 主函数：热-扩散应力计算
function thermal_diffusion_stress_2D(case, variables)
    # 集成上述函数，返回更新的 variables
end
```

### 8.2 需要的材料参数

在参数文件中添加：
```julia
# 负极
NE.beta_c = 0.03        # 浓度膨胀系数 [-]
NE.c_s_max = 31507.0    # 最大锂浓度 [mol/m³]

# 正极
PE.beta_c = 0.02
PE.c_s_max = 51765.0
```

### 8.3 变量命名约定

输入：
- `T_nodes`: 节点温度 [K]
- `SOC_elem`: 单元SOC [-]

输出：
- `displacement_x`, `displacement_y`: 节点位移 [m]
- `stress_xx`, `stress_yy`, `stress_xy`: 单元应力 [Pa]
- `stress_vonMises`: Von Mises等效应力 [Pa]

## 9. 验证方法

### 9.1 纯热应力验证

设置 `β_c = 0`，仅考虑热应力，与现有 `thermal_stress()` 函数对比。

### 9.2 纯扩散应力验证

设置 `α = 0`，仅考虑扩散应力，与文献数据对比。

### 9.3 耦合验证

同时考虑热应力和扩散应力，检查应力叠加原理。

## 10. 物理参数估计

### 10.1 浓度膨胀系数

基于文献数据：
- 石墨负极：`β_c ≈ 0.01 - 0.04` (体积变化约 10-13%)
- LCO/NMC正极：`β_c ≈ 0.005 - 0.02` (体积变化约 2-5%)

### 10.2 弹性模量和泊松比

- 石墨：`E ≈ 10-15 GPa`, `ν ≈ 0.3`
- LCO/NMC：`E ≈ 100-200 GPa`, `ν ≈ 0.3`
- 复合电极（考虑孔隙率）：需要使用有效模量

## 11. 参考文献

1. Zhang, X., et al. (2007). "Numerical simulation of intercalation-induced stress in Li-ion battery electrode particles." Journal of the Electrochemical Society, 154(10), A910-A916.

2. Christensen, J., & Newman, J. (2006). "Stress generation and fracture in lithium insertion materials." Journal of Solid State Electrochemistry, 10(5), 293-319.

3. Zhao, K., et al. (2011). "Concurrent reaction and plasticity during initial lithiation of crystalline silicon in lithium-ion batteries." Journal of the Electrochemical Society, 159(3), A238-A243.

4. Dai, Y., et al. (2016). "On graded electrode porosity as a design tool for improving the energy density of batteries." Journal of The Electrochemical Society, 163(3), A406-A416.
