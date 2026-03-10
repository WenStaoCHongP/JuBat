# 层间热传导修复说明

## 问题描述

在使用果冻卷（Jellyroll）二维分布式热传导模型时，发现 `thermal2D`（未合并节点）和 `thermal2D_merged`（合并节点）网格的仿真结果存在显著差异。

### 根本原因

在 `jellyroll_collector_seed_mesh` 函数生成的网格中：

1. **thermal2D（未合并节点）**：相邻螺旋圈的界面节点在空间位置重合，但拥有不同的节点编号。这导致在刚度矩阵装配时，相邻圈的节点之间没有热传导耦合项，形成**层间热绝缘**。

2. **thermal2D_merged（合并节点）**：重合的节点被合并成单一节点，相邻圈自然通过共享节点实现热传导。

### 物理影响

| 网格类型 | 层间热传导 | 温度分布 |
|---------|-----------|---------|
| `thermal2D_merged` | ✅ 自动通过共享节点传导 | 连续 |
| `thermal2D`（修复前） | ❌ 界面节点独立，无耦合 | 层间热绝缘 |
| `thermal2D`（修复后） | ✅ 通过耦合项传导 | 连续 |

## 解决方案

### 1. 新增函数：`_apply_interlayer_thermal_conductance!`

在 `ThermalDistributed.jl` 中新增函数，用于在界面节点对之间添加热传导耦合项。

**物理模型**：
使用完美接触假设，在界面节点对 (n_out, n_in) 之间添加热传导耦合：

```
KT[n_out, n_out] -= h_interface * A_rep / (k_th * L_th)
KT[n_out, n_in]  += h_interface * A_rep / (k_th * L_th)
KT[n_in, n_out]  += h_interface * A_rep / (k_th * L_th)
KT[n_in, n_in]   -= h_interface * A_rep / (k_th * L_th)
```

其中：
- `h_interface = λ_r * 100 / t_repeat`：界面等效换热系数 [W/(m²·K)]
- `λ_r`：径向等效热导率 [W/(m·K)]
- `t_repeat`：完整周期厚度 [m]
- `A_rep`：每个节点对代表的界面弧长 [m]
- `k_th`, `L_th`：热传导参考值，用于无量纲化

### 2. 修改函数：`ThermalDistributed2D_BC`

在应用边界条件时自动调用层间热传导函数：

```julia
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64=0.0)
    # ... 原有代码 ...
    
    # 应用层间热传导（关键修改！）
    _apply_interlayer_thermal_conductance!(KT, mesh, case)
    
    return nothing
end
```

### 3. 新增辅助函数

- `_get_or_compute_interface_pairs`：获取或计算界面节点对
- `_compute_interface_pairs_from_mesh`：从网格自动识别界面节点对

### 4. 新增便捷函数：`setup_thermal2D_mesh!`

在 `Jellyrollmodel.jl` 中新增函数，用于设置热网格并自动保存界面信息：

```julia
# 推荐使用方式
param_dim = JuBat.ChooseCell("Jellyroll")
case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
JuBat.setup_thermal2D_mesh!(case, mesh_data)  # 自动保存interface_pairs
```

## 使用方法

### 方法1（推荐）：使用 `setup_thermal2D_mesh!`

```julia
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
JuBat.setup_thermal2D_mesh!(case, mesh_data)
```

### 方法2：手动设置（向后兼容）

```julia
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
case.mesh["thermal2D"] = mesh_data.thermal2D

# 可选：保存interface_pairs以提高效率
case.multi_spme_layout["interface_pairs"] = mesh_data.interface_pairs
```

即使不保存 `interface_pairs`，系统也会自动从网格计算。

### 方法3：使用合并节点网格（不推荐用于CZM）

```julia
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
JuBat.setup_thermal2D_mesh!(case, mesh_data; use_merged=true)
```

注意：使用合并节点网格时，CZM（内聚力区域模型）将无法正常工作。

## 验证

修复后，使用 `thermal2D` 和 `thermal2D_merged` 的温度场结果应该高度一致。

## 文件修改列表

1. `src/ThermalDistributed.jl`
   - 修改 `ThermalDistributed2D_BC`
   - 新增 `_apply_interlayer_thermal_conductance!`
   - 新增 `_get_or_compute_interface_pairs`
   - 新增 `_compute_interface_pairs_from_mesh`

2. `src/Jellyrollmodel.jl`
   - 新增 `setup_thermal2D_mesh!`
   - 更新导出列表

3. `src/JuBat.jl`
   - 新增导出 `setup_thermal2D_mesh!`
