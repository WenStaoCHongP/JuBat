# Initialisation.jl 优化方案

> 日期: 2026-04-01 (修订)
> 文件: `src/Initialisation.jl`
> 状态: **已实施** (2026-04-07)
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
- 重构 4 个 `MultiSPMe_*` 函数（改名 + 类型替换），不新增函数
- 所有逻辑保持内联，不提取子函数

### 3.2 函数改名

| 旧名 | 新名 | 说明 |
|------|------|------|
| `ModelInitialisation_MultiSPMe` | `model_initialisation_multi_spme` | snake_case |
| `MultiSPMe_extract_element_state` | `extract_element_state` | 去前缀 + snake_case |
| `MultiSPMe_get_thermal_dofs` | `get_thermal_dofs` | 去前缀 + snake_case |
| `MultiSPMe_update_state` | `update_state` | 去前缀 + snake_case |

### 3.3 类型替换

所有函数中 `case.multi_spme_layout["key"]` 替换为 `case.layout.key` / `case.geometry.key`：

```julia
# 旧:
layout = case.multi_spme_layout
ne = layout["ne"]
n_chem = layout["n_chem"]
thermal_range = layout["thermal_range"]

# 新:
layout = case.layout
ne = layout.ne
n_chem = layout.n_chem
thermal_range = layout.thermal_range
```

### 3.4 `model_initialisation_multi_spme` 内部简化

```julia
# 旧 (行 110-116):
empty!(case.multi_spme_layout)
case.multi_spme_layout["ne"] = ne
case.multi_spme_layout["n_chem"] = n_chem
case.multi_spme_layout["nT"] = nT
case.multi_spme_layout["n_total"] = length(y0)
case.multi_spme_layout["chem_range"] = 1:(ne * n_chem)
case.multi_spme_layout["thermal_range"] = (ne * n_chem + 1):(ne * n_chem + nT)

# 新 (1 行):
case.layout = MultiSPMeLayout(ne, n_chem, nT)
```

**非均匀 SOC 逻辑保持内联**，不提取为子函数。

### 3.5 三个辅助函数精简（仅改名 + 类型替换）

```julia
# 旧 MultiSPMe_extract_element_state → 新 extract_element_state:
function extract_element_state(y::AbstractVector, e::Int, layout::MultiSPMeLayout)
    offset = (e - 1) * layout.n_chem
    return y[(offset + 1):(offset + layout.n_chem)]
end

# 旧 MultiSPMe_get_thermal_dofs → 新 get_thermal_dofs:
function get_thermal_dofs(y::AbstractVector, layout::MultiSPMeLayout)
    return y[layout.thermal_range]
end

# 旧 MultiSPMe_update_state → 新 update_state:
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

### 3.6 调用点更新

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
| 总行数 | 239 | ~230 |
| 函数数 | 5 | 5（仅改名，不新增） |
| Dict 访问 | ~20 处 | 0 |
| 类型安全 | 无 | 全部有 |
| `ModelInitialisation` 行数 | 48 (不变) | 48 (不变) |

---

## 5. 实施记录 (2026-04-07)

### 实际偏差

| 原方案 | 实际实施 | 原因 |
|--------|---------|------|
| `ModelInitialisation_MultiSPMe` 改名为 `model_initialisation_multi_spme` | 保持原名 | 避免同时改名+改类型降低出错风险 |
| 3 个辅助函数改名 + 类型替换 | 已完成（`extract_element_state`/`get_thermal_dofs`/`update_state`） | 按 05 文档执行 |
| `MultiSPMe_update_state` 不含验证 | 保留 `@assert` 验证 | 审查建议采纳 |

### 实际效果

| 指标 | 旧 | 新 | 备注 |
|------|-----|-----|------|
| `multi_spme_layout` Dict 访问 | ~20 处 | **0 处** | 全局替换完成 |
| 布局构建重复 | 4 处 x 6 行 | **4 处 x 1 行** | `MultiSPMeLayout` 构造器 |
| `_ensure_multi_spme_layout!` | 28 行 | **已删除** | 不再需要 |
| `haskey` 检查 | 2 处 | **0 处** | fail-fast 原则 |
| `layer_weights` try/catch | 1 处 | **已删除** | 迁入 MeshGeometry |

### 新增文件
- `src/CouplingState.jl` — `MultiSPMeLayout` + `MeshGeometry` struct 定义

### 连锁修改文件
| 文件 | 改动点数 |
|------|---------|
| `SetCase.jl` | Case struct + 构造器 |
| `JuBat.jl` | include + export |
| `Solve.jl` | ~20 处 |
| `CycleSolver.jl` | 删除 28 行函数 |
| `CycleData.jl` | 4 处 |
| `PostProcessing.jl` | 3 处（规格遗漏，实施中发现） |
| `Jellyrollmodel.jl` | MeshGeometry 构造 |
| `ThermalDistributed.jl` | haskey 替换 |

### 验证
- 全局搜索 `multi_spme_layout`：仅在注释中出现（1 处）
- 旧函数名搜索：零匹配
- `haskey` 搜索：零匹配
- Julia 冒烟测试：Case 构造 + MultiSPMeLayout 构造 均通过
