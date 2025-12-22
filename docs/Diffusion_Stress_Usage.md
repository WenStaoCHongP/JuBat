# 宏观扩散应力计算使用指南

## 概述

本文档介绍如何在 JuBat 中使用新增的宏观扩散应力计算功能。该功能基于有限元方法，计算锂浓度变化和温度变化共同引起的应力场和位移场。

## 快速开始

### 1. 基本用法

```julia
using JuBat

# 创建案例
param_dim = JuBat.ChooseCell("LG M50")
opt = JuBat.Option()
opt.model = "SPMe"
case = JuBat.SetCase(param_dim, opt)

# 创建2D网格（必须是Q4单元）
mesh_th = ... # 创建或加载Q4网格
case.mesh["thermal2D"] = mesh_th

# 准备变量字典
variables = Dict{String, Union{Array{Float64},Float64}}()
variables["T_nodes"] = T_nodes  # 节点温度
variables["SOC_elem"] = SOC_elem  # 单元SOC

# 计算应力
variables = JuBat.diffusion_stress_2D(case, variables)

# 提取结果
σ_xx = variables["diffusion stress xx"]  # [Pa]
σ_yy = variables["diffusion stress yy"]  # [Pa]
σ_xy = variables["diffusion stress xy"]  # [Pa]
σ_vm = variables["diffusion stress vonMises"]  # [Pa]
u_x = variables["displacement x"]  # [m]
u_y = variables["displacement y"]  # [m]
```

### 2. 与电化学模型耦合

```julia
# 先求解电化学模型
result = JuBat.Solve(case)

# 提取浓度信息
variables["negative electrode lithium concentration"] = result["negative particle lithium concentration"][:, end]
variables["positive electrode lithium concentration"] = result["positive particle lithium concentration"][:, end]

# 提取温度信息
T_final = result["temperature[K]"][end]
variables["T_nodes"] = fill(T_final / param.scale.T_ref, mesh_th.nlen)

# 计算应力
variables = JuBat.diffusion_stress_2D(case, variables)
```

## 理论基础

### 控制方程

1. **应力平衡方程**（准静态）：
   ```
   ∇·σ = 0
   ```

2. **应变-位移关系**：
   ```
   ε_xx = ∂u/∂x
   ε_yy = ∂v/∂y
   γ_xy = ∂u/∂y + ∂v/∂x
   ```

3. **本构关系**（平面应力）：
   ```
   σ = E/(1-ν²) * [(ε - ε_0) + ν(ε - ε_0)]
   ```
   
   其中初始应变：
   ```
   ε_0 = α*ΔT + β_c*c_s_max*ΔSOC
   ```

### 离散方程

使用有限元方法（Q4单元），得到：

```
[K] {U} = {F}
```

其中：
- `[K]`: 刚度矩阵，尺寸 (2N×2N)，N 为节点数
- `{U}`: 位移向量，包含所有节点的 u 和 v
- `{F}`: 载荷向量，由热应变和扩散应变引起

## 必需的材料参数

在参数文件中需要定义以下参数：

```julia
# 负极
NE.E = 15e9        # 弹性模量 [Pa]
NE.nu = 0.3        # 泊松比 [-]
NE.alphaT = 1.5e-5 # 热膨胀系数 [1/K]
NE.beta_c = 0.03   # 浓度膨胀系数 [-]
NE.cs_max = 33133  # 最大锂浓度 [mol/m³]

# 正极
PE.E = 375e9       # 弹性模量 [Pa]
PE.nu = 0.2        # 泊松比 [-]
PE.alphaT = 1.0e-5 # 热膨胀系数 [1/K]
PE.beta_c = 0.01   # 浓度膨胀系数 [-]
PE.cs_max = 63104  # 最大锂浓度 [mol/m³]

# 隔膜
SP.E = 1e9         # 弹性模量 [Pa]
SP.nu = 0.3        # 泊松比 [-]
SP.alphaT = 1.0e-5 # 热膨胀系数 [1/K]
```

### 参数说明

1. **弹性模量 E [Pa]**：
   - 石墨负极：10-15 GPa
   - NMC/LCO正极：100-200 GPa
   - 隔膜：0.5-2 GPa
   - 集流体（Al）：70 GPa
   - 集流体（Cu）：130 GPa

2. **泊松比 ν [-]**：
   - 大多数电极材料：0.2-0.3
   - 金属：0.3-0.35

3. **热膨胀系数 α [1/K]**：
   - 石墨：1-2 × 10⁻⁵ 1/K
   - NMC/LCO：1-1.5 × 10⁻⁵ 1/K
   - 铝：23 × 10⁻⁶ 1/K
   - 铜：17 × 10⁻⁶ 1/K

4. **浓度膨胀系数 β_c [-]**：
   - 定义：Δε_vol / Δc_s，体积应变相对于浓度变化的敏感度
   - 石墨：0.01-0.04（文献报道体积变化10-13%）
   - NMC/LCO：0.005-0.02（体积变化2-5%）

## 输入变量

### 必需输入

1. **温度场** `T_nodes`：
   - 类型：Vector{Float64}，长度 = 节点数
   - 单位：无量纲（除以 T_ref）或直接用 [K]
   - 说明：每个节点的温度

2. **SOC分布** `SOC_elem`：
   - 类型：Vector{Float64}，长度 = 单元数
   - 单位：无量纲 [0-1]
   - 说明：每个单元的锂化程度

### 可选输入

1. **层权重矩阵** `thermal2D layer_weights`：
   - 类型：Matrix{Float64}，尺寸 (ne × 5)
   - 说明：每个单元中各层（NE, SP, PE, PCC, NCC）的体积分数
   - 用途：计算等效材料参数

2. **锂浓度分布**：
   - `negative electrode lithium concentration`
   - `positive electrode lithium concentration`
   - 用途：自动估计 SOC 分布

## 输出变量

计算完成后，`variables` 字典中会添加以下字段：

1. **应力场**（单元中心值）：
   - `diffusion stress xx` [Pa]：x方向正应力
   - `diffusion stress yy` [Pa]：y方向正应力
   - `diffusion stress xy` [Pa]：剪应力
   - `diffusion stress vonMises` [Pa]：Von Mises等效应力

2. **位移场**（节点值）：
   - `displacement x` [m]：x方向位移
   - `displacement y` [m]：y方向位移

## 边界条件

当前实现的边界条件策略：

1. **最小约束**：固定一个节点（最接近原点的节点）以防止刚体位移
2. **自由边界**：所有其他节点不施加约束

### 自定义边界条件

如需修改边界条件，可以编辑 `_identify_mechanical_bc_nodes` 函数：

```julia
function _identify_mechanical_bc_nodes(mesh, case)
    bc_nodes = Dict{Int, Symbol}()
    
    # 例如：固定左边界
    x = mesh.node[:, 1]
    x_min = minimum(x)
    for i in 1:mesh.nlen
        if abs(x[i] - x_min) < 1e-6
            bc_nodes[i] = :fixed_x  # 固定x方向
        end
    end
    
    return bc_nodes
end
```

支持的边界条件类型：
- `:fixed_x`：固定x方向位移
- `:fixed_y`：固定y方向位移
- `:fixed_xy`：固定x和y方向位移

## 网格要求

1. **单元类型**：必须是 Q4（四节点四边形）单元
2. **维度**：2D平面
3. **网格质量**：
   - 避免过度扭曲的单元
   - 推荐长宽比 < 5
   - 单元尺寸应与物理尺度匹配

## 可视化示例

```julia
using Plots

# 计算单元中心坐标
ne = size(mesh.element, 1)
x_elem = [mean(mesh.node[mesh.element[e, :], 1]) for e in 1:ne]
y_elem = [mean(mesh.node[mesh.element[e, :], 2]) for e in 1:ne]

# Von Mises应力云图
σ_vm = variables["diffusion stress vonMises"]
scatter(x_elem, y_elem, marker_z=σ_vm./1e6,
        color=:plasma, markersize=3,
        xlabel="x [m]", ylabel="y [m]",
        title="Von Mises Stress [MPa]",
        colorbar=true, aspect_ratio=:equal)
savefig("stress_field.png")

# 位移矢量场
u_x = variables["displacement x"]
u_y = variables["displacement y"]
quiver(mesh.node[:, 1], mesh.node[:, 2],
       quiver=(u_x.*1e6, u_y.*1e6),
       xlabel="x [m]", ylabel="y [m]",
       title="Displacement Field [μm]",
       aspect_ratio=:equal)
savefig("displacement_field.png")
```

## 注意事项

1. **单位一致性**：
   - 温度可以是无量纲或有量纲，程序会自动处理
   - SOC 必须是无量纲 [0-1]
   - 输出应力单位：Pa
   - 输出位移单位：m

2. **计算效率**：
   - 刚度矩阵装配的复杂度：O(n_elem × n_gauss × n_node_per_elem²)
   - 求解线性方程组：O(n_dof^3)（直接法）或 O(n_dof)（迭代法）
   - 对于大规模问题（> 10,000 节点），考虑使用迭代求解器

3. **数值稳定性**：
   - 泊松比应 < 0.5（可压缩材料）
   - 弹性模量应 > 0
   - 避免极端的材料参数比值（如 E_ratio > 10⁶）

4. **物理合理性检查**：
   - Von Mises应力应 < 材料屈服强度
   - 位移应与单元尺寸可比
   - 检查应力平衡（∫σ·n dS ≈ 0）

## 与热应力的对比

| 项目 | 纯热应力 | 纯扩散应力 | 耦合 |
|------|---------|-----------|------|
| 初始应变 | α*ΔT | β_c*c_s_max*ΔSOC | α*ΔT + β_c*c_s_max*ΔSOC |
| 驱动因素 | 温度梯度 | SOC梯度 | 两者同时 |
| 典型量级 | 10-100 MPa | 10-500 MPa | 叠加 |

## 故障排除

### 问题1：求解失败
**症状**：`Mechanical solve failed` 警告

**可能原因**：
- 刚度矩阵奇异（缺少边界约束）
- 材料参数不合理

**解决方案**：
- 检查边界条件
- 验证材料参数
- 增加数值稳定性（如添加少量阻尼）

### 问题2：应力过大
**症状**：应力值 > GPa 级别

**可能原因**：
- 温度或SOC变化过大
- 膨胀系数设置不当

**解决方案**：
- 检查输入温度和SOC的物理合理性
- 重新校准 α 和 β_c 参数

### 问题3：位移为零
**症状**：所有位移 ≈ 0

**可能原因**：
- 温度和SOC都没有变化（ΔT=0, ΔSOC=0）
- 膨胀系数 α 或 β_c 设置为 0

**解决方案**：
- 确保存在温度或浓度梯度
- 检查材料参数是否正确加载

## 示例代码

完整的示例代码见：
- `example/thermal_diffusion_stress_example.jl`

## 参考文献

1. **理论推导**：
   - `docs/Diffusion_Stress_Macroscale_Theory.md`

2. **相关论文**：
   - Zhang, X., et al. (2007). "Numerical simulation of intercalation-induced stress in Li-ion battery electrode particles." J. Electrochem. Soc., 154(10), A910-A916.
   - Christensen, J., & Newman, J. (2006). "Stress generation and fracture in lithium insertion materials." J. Solid State Electrochem., 10(5), 293-319.
   - Zhao, K., et al. (2011). "Concurrent reaction and plasticity during initial lithiation of crystalline silicon in lithium-ion batteries." J. Electrochem. Soc., 159(3), A238-A243.

## 联系与支持

如有问题或建议，请：
1. 查阅理论文档
2. 检查示例代码
3. 提交 Issue 到代码仓库
