# 层间热传导实现说明

本文档按当前代码更新，说明当前仓库中“层间热传导问题”是如何处理的，以及与早期修复方案文档的差异。

---

## 1. 问题背景

在果冻卷二维热模型中，相邻卷绕圈在几何上可能共享同一空间位置，但在未合并网格里仍然对应不同节点编号。

这会带来两类不同的热学建模方式：

1. 使用合并节点网格时，相邻圈通过共享节点自然导热。
2. 使用未合并节点网格时，必须显式添加界面耦合项，否则层间会表现为热绝缘。

---

## 2. 当前代码采用的总体策略

当前代码并不是单纯给 `thermal2D` 打一个“统一的层间导热补丁”，而是采用了两套路径：

### 路径 A：不启用 CZM 时

默认使用 `thermal2D_merged`。

此时层间传热通过共享节点自然存在，不需要额外界面耦合项。

### 路径 B：启用 CZM 时

默认使用未合并的 `thermal2D`，然后在 `ThermalDistributed2D_BC` 中根据 `czm_mesh` 的界面单元和损伤状态显式装配界面导热项。

因此，当前实现的真实口径应为：

- 非 CZM 场景优先用 merged 网格解决层间传热。
- CZM 场景用 unmerged 网格保留界面自由度，再显式加入界面热耦合。

---

## 3. 当前网格生成与选择方式

### 3.1 `jellyroll_collector_seed_mesh`

当前 `src/Jellyrollmodel.jl` 会一次性返回两套热网格：

```julia
mesh_data = jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)

mesh_data.thermal2D
mesh_data.thermal2D_merged
mesh_data.interface_pairs
mesh_data.czm_element_map
mesh_data.is_inner_layer
```

其中：

- `thermal2D` 是未合并节点网格。
- `thermal2D_merged` 是重合节点已合并的网格。
- `interface_pairs` 是在未合并网格中识别出的重合内外节点对。

### 3.2 `setup_thermal2D_mesh`

当前正式接口是：

```julia
case = setup_thermal2D_mesh(case, mesh_data; use_merged=nothing)
```

注意两点：

1. 该函数不是原地修改版本，因此没有 `setup_thermal2D_mesh!`。
2. 默认选择逻辑是：
   - `czm_enabled=false` 时使用 `thermal2D_merged`
   - `czm_enabled=true` 时使用 `thermal2D`

所以，当前最推荐的用法是：

```julia
param_dim = JuBat.ChooseCell("Jellyroll")
case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
```

---

## 4. 当前层间热耦合真正发生的位置

当前仓库中没有早期文档里提到的：

- `_apply_interlayer_thermal_conductance!`
- `_get_or_compute_interface_pairs`
- `_compute_interface_pairs_from_mesh`

当前真正生效的入口是：

```julia
ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)
```

其逻辑是：

1. 当 `case.opt.czm_enabled` 为真时，读取或即时创建 `czm_mesh`。
2. 遍历每个 `cohesive_element`。
3. 读取该界面的 `DamageState.D` 与 `DamageState.δ_max_n`。
4. 调用 `compute_gap_conductance(D, δ_n, case.param_dim.cohesive)`。
5. 将界面导热项直接装配进热刚度矩阵 `K`。

当前装配形式为：

```julia
coeff = h_eff * czm_elem.length / (k_th * L_th)

K[nb, nb] -= coeff
K[nb, nt] += coeff
K[nt, nb] += coeff
K[nt, nt] -= coeff
```

这就是当前代码里“未合并网格仍可层间导热”的来源。

---

## 5. 当前界面导热系数模型

早期文档中写的是固定 `h_interface` 的完美接触模型。当前代码已经改成基于 CZM 状态的退化导热模型。

当前使用：

```julia
compute_gap_conductance(D, δ_n, cohesive)
```

其参数来自 `Cohesive`：

```julia
h_c0
k_air
lambda_m
beta
threshold
δ_0_n
δ_c_n
```

因此，当前层间导热不是常数，而是随：

- 损伤变量 `D`
- 界面法向开裂位移 `δ_n`

共同变化。

---

## 6. 对“修复是否完成”的准确表述

如果把问题定义为“未合并网格在热学上完全绝缘”，那么当前代码已经提供了解决方案，但需要分场景理解：

### 6.1 非 CZM 场景

通过默认选择 `thermal2D_merged` 规避问题。

### 6.2 CZM 场景

通过 `ThermalDistributed2D_BC` 中的界面导热项解决问题。

### 6.3 不能再沿用的旧说法

以下表述已经不再准确：

- “系统会自动从任意未合并网格计算 interface_pairs 并注入层间导热补丁。”
- “修复后 thermal2D 与 thermal2D_merged 在所有场景下应高度一致。”

更准确的说法是：

- 在非 CZM 模式下，默认不直接使用未合并网格做层间导热。
- 在 CZM 模式下，未合并网格通过界面单元导热项恢复层间热耦合，但其结果会受到损伤状态影响，因此不必与 merged 网格完全一致。

---

## 7. 当前推荐使用方式

### 7.1 常规热仿真

```julia
opt = JuBat.Option(
    thermal_enabled = true,
    thermalmodel = "distributed2D",
    czm_enabled = false,
)

mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
```

此时默认会选用 `thermal2D_merged`。

### 7.2 启用 CZM 的热-力耦合场景

```julia
opt = JuBat.Option(
    thermal_enabled = true,
    thermalmodel = "distributed2D",
    czm_enabled = true,
)

mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
czm_mesh = JuBat.create_czm_mesh(case.mesh["thermal2D"], param_dim)
case.multi_spme_layout["czm_mesh"] = czm_mesh
```

此时默认会选用未合并网格，并在热边界条件装配时加入界面导热项。

---

## 8. 相关文件

当前与该问题直接相关的文件如下：

1. `src/Jellyrollmodel.jl`
   - 生成 `thermal2D` / `thermal2D_merged`
   - 计算 `interface_pairs`
   - 提供 `setup_thermal2D_mesh`

2. `src/ThermalDistributed.jl`
   - 在 `ThermalDistributed2D_BC` 中装配 CZM 界面导热项

3. `src/Materialmatrix.jl`
   - 定义 `compute_gap_conductance`

4. `src/czm.jl`
   - 定义 `create_czm_mesh`
   - 保存界面单元与损伤状态

5. `src/JuBat.jl`
   - 导出 `setup_thermal2D_mesh`
   - 导出 `compute_gap_conductance`

---

## 9. 当前结论

当前代码已经不再采用早期文档中描述的“固定界面导热补丁函数”方案，而是演进为：

- 用 merged 网格处理无损伤层间传热。
- 用 unmerged + CZM 界面导热耦合处理可退化界面。

这比早期方案更贴近当前代码结构，也更符合后续电化学-热-力耦合扩展方向。
