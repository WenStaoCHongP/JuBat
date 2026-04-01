# Initialisation.jl 优化方案

> 日期: 2026-04-01
> 文件: `src/Initialisation.jl`
> 状态: 修改 (43→239 行, +196%)
> main 分支行数: 43

---

## 1. main 分支现状

### 1.1 唯一函数: `ModelInitialisation(case::Case)`

43 行，流程清晰：
1. 检查 `case.opt.y0` 是否为空
2. 按 model 类型分支：
   - SPM → `y0 = [csn0; csp0]`
   - SPMe → `y0 = [csn0; csp0; ce0]`
   - P2D → `y0 = [csn0; csp0; ce0]` + 后追电位
3. 可选 lumped 热追加 `T0`
4. P2D 追加电位 `phis_n, phis_p, phie0`

**设计简洁，不做任何修改。**

---

## 2. 当前分支新增内容

### 2.1 函数清单

| 函数 | 行范围 | 行数 | 职责 |
|------|--------|------|------|
| `ModelInitialisation` | 1-48 | 48 | 标准 (仅追 distributed2D 分支) |
| `ModelInitialisation_MultiSPMe` | 54-119 | 66 | 多SPMe初始化 + layout 构建 |
| `MultiSPMe_extract_element_state` | 149-161 | 13 | 提取单元状态 |
| `MultiSPMe_get_thermal_dofs` | 183-193 | 11 | 提取热 DOF |
| `MultiSPMe_update_state` | 205-239 | 35 | 更新全局状态 |

### 2.2 `ModelInitialisation` 的变更

仅增加了一个 `elseif` 分支（行 35-40）：

```julia
elseif case.opt.thermalmodel == "distributed2D"
    nT = case.mesh["thermal2D"].nlen
    T0_nodes = fill(case.param.cell.T0, nT)
    y0 = [y0; T0_nodes]
```

**这部分保留不动**——它是在标准 SPMe + distributed2D 单元路径下的正确逻辑。

---

## 3. 优化方案

### 3.1 约束

- `ModelInitialisation()` 的 SPM/SPMe/P2D/lumped 分支不动
- `ModelInitialisation()` 的 distributed2D 分支保留
- 仅重构 4 个 `MultiSPMe_*` 函数

### 3.2 `ModelInitialisation_MultiSPMe` 简化

```julia
# ===== 旧 (66 行) =====
function ModelInitialisation_MultiSPMe(case::Case; initial_soc_distribution=nothing)
    ne = size(case.mesh["thermal2D"].element, 1)
    nT = case.mesh["thermal2D"].nlen
    # 临时修改 thermalmodel → "none" 来获取纯电化学部分
    original_thermalmodel = case.opt.thermalmodel
    case.opt.thermalmodel = "none"
    y0_single_chem = ModelInitialisation(case)
    case.opt.thermalmodel = original_thermalmodel
    n_chem = length(y0_single_chem)
    # ... 30+ 行手动初始化 ...
    # 最后 6 行构建 layout Dict
    empty!(case.multi_spme_layout)
    case.multi_spme_layout["ne"] = ne
    # ... 5 行 Dict 赋值 ...
end

# ===== 新 (~35 行) =====
function model_initialisation_multi_spme(case::Case;
        initial_soc_distribution::Union{Nothing, Vector{Float64}}=nothing)
    ne = size(case.mesh["thermal2D"].element, 1)
    nT = case.mesh["thermal2D"].nlen

    # 获取单单元电化学 DOF
    original_thermalmodel = case.opt.thermalmodel
    case.opt.thermalmodel = "none"
    y0_single = ModelInitialisation(case)
    case.opt.thermalmodel = original_thermalmodel
    n_chem = length(y0_single)

    # 设置 layout（替代 6 行 Dict 构建）
    case.layout = MultiSPMeLayout(ne, n_chem, nT)

    # 批量初始化电化学状态
    y0_chem = repeat(y0_single, ne)

    # 非均匀 SOC 处理（仅在提供了 distribution 时）
    if initial_soc_distribution !== nothing
        _apply_nonuniform_soc!(y0_chem, initial_soc_distribution, case)
    end

    # 热场 + 组装
    T0_nodes = fill(case.param.cell.T0, nT)
    return [y0_chem; T0_nodes]
end

# 辅助：非均匀 SOC 应用
function _apply_nonuniform_soc!(y0_chem::Vector{Float64},
                                 soc_dist::Vector{Float64}, case::Case)
    Nrn = case.mesh["negative particle"].nlen
    Nrp = case.mesh["positive particle"].nlen
    n_chem = case.layout.n_chem
    for e in 1:length(soc_dist)
        offset = (e - 1) * n_chem
        soc_e = soc_dist[e]
        # NE 浓度
        cn = case.param.NE.cs0 * soc_e
        y0_chem[(offset + 1):(offset + Nrn)] .= cn
        # PE 浓度
        cp = case.param.PE.cs0 * (1.0 - soc_e)
        y0_chem[(offset + Nrn + 1):(offset + Nrn + Nrp)] .= cp
    end
end
```

**改进点**：
- 用 `repeat(y0_single, ne)` 替代手动 for 循环填充
- layout 构建从 6 行 Dict → 1 行 struct
- 非均匀 SOC 逻辑提取为独立辅助函数

### 3.3 三个辅助函数精简

```julia
# ===== 旧 MultiSPMe_extract_element_state (13 行) =====
function MultiSPMe_extract_element_state(y::Array{Float64}, e::Int, case::Case)
    y_vec = vec(y)
    layout = case.multi_spme_layout
    ne = layout["ne"]
    n_chem = layout["n_chem"]
    offset = (e - 1) * n_chem
    yt_e = y_vec[(offset + 1):(offset + n_chem)]
    return yt_e
end

# ===== 新 (~5 行) =====
function extract_element_state(y::AbstractVector, e::Int, layout::MultiSPMeLayout)
    offset = (e - 1) * layout.n_chem
    return y[(offset + 1):(offset + layout.n_chem)]
end

# ===== 旧 MultiSPMe_get_thermal_dofs (11 行) =====
function MultiSPMe_get_thermal_dofs(y::Array{Float64}, case::Case)
    y_vec = vec(y)
    layout = case.multi_spme_layout
    thermal_range = layout["thermal_range"]
    T_nodes = y_vec[thermal_range]
    return T_nodes
end

# ===== 新 (~3 行) =====
function get_thermal_dofs(y::AbstractVector, layout::MultiSPMeLayout)
    return y[layout.thermal_range]
end

# ===== 旧 MultiSPMe_update_state (35 行) =====
# 含大量手动 Dict 检查、边界检查、错误消息

# ===== 新 (~12 行) =====
function update_state(y::AbstractVector, layout::MultiSPMeLayout;
                      element_index::Union{Nothing,Int}=nothing,
                      element_state::Union{Nothing,Vector{Float64}}=nothing,
                      thermal_nodes::Union{Nothing,Vector{Float64}}=nothing)
    y_new = copy(y)
    if element_index !== nothing
        offset = (element_index - 1) * layout.n_chem
        y_new[(offset + 1):(offset + layout.n_chem)] .= element_state
    end
    if thermal_nodes !== nothing
        y_new[layout.thermal_range] .= thermal_nodes
    end
    return y_new
end
```

### 3.4 调用点更新

所有调用 `MultiSPMe_extract_element_state(y, e, case)` 的地方改为：

```julia
# 旧:
yt_e = MultiSPMe_extract_element_state(yt, e, case)

# 新:
yt_e = extract_element_state(vec(yt), e, case.layout)
```

受影响位置：
- `Solve.jl` CallModel_MultiSPMe 内 ~2 处
- `CycleSolver.jl` ~1 处
- `CycleData.jl` ~1 处

---

## 4. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| 总行数 | 239 | ~120 |
| 函数数 | 5 | 6 (+1 辅助) |
| Dict 访问 | ~20 处 | 0 |
| 类型安全 | 无 | 全部有 |
| `ModelInitialisation` 行数 | 48 (不变) | 48 (不变) |
