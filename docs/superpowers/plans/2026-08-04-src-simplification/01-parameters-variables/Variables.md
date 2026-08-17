# Variables.jl Simplification Plan

**Status:** ⚠️ Blocked（评审后待重写）
**Layer:** 1 - 参数与变量
**桶:** Consolidate（D2 重复簇）
**依赖:** `VariableKeys.md` 必须先完成
**阻塞:** 下游所有用 `variables["..."]` 的文件

**Goal:** 删除 `StandardVariables` (line 1) 与 `create_element_workspace` (line 157) 两份独立硬编码键字面量表，统一引用 `VariableKeys` 模块。

**Architecture:** `StandardVariables` 内部改为 `Dict(name => zeros(...) for name in VariableKeys.GLOBAL_KEYS)` 形式（元组驱动）；`create_element_workspace` 同理。键名只在 `VariableKeys.jl` 定义一次。

---

## Files

- Modify: `src/Variables.jl:1-156`（`StandardVariables`）
- Modify: `src/Variables.jl:157-235`（`create_element_workspace`）
- Reference: `src/VariableKeys.jl`

---

## Tasks

### Task 1: Characterization Test（先写测试）

**Files:**
- Create: `test/unit_variables_keys.jl`

- [ ] **Step 1.1: 写键快照测试**

```julia
# test/unit_variables_keys.jl
include("../src/JuBat.jl")
using .JuBat

@testset "StandardVariables 键集合稳定" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    vars = JuBat.StandardVariables(case, 1)
    expected_keys = Set(VariableKeys.GLOBAL_KEYS)
    @test Set(keys(vars)) == expected_keys
end

@testset "create_element_workspace 键集合稳定" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    ws = JuBat.create_element_workspace(case)
    expected_keys = Set(VariableKeys.ELEMENT_WORKSPACE_KEYS)
    @test Set(keys(ws)) == expected_keys
end
```

- [ ] **Step 1.2: 跑测试验证当前实现通过**

Run: `julia test/unit_variables_keys.jl`
Expected: PASS（捕获当前行为作快照）

- [ ] **Step 1.3: Commit**

```bash
git add test/unit_variables_keys.jl
git commit -m "test(Variables): 添加键集合 characterization test"
```

### Task 2: 重构 `StandardVariables` 用常量元组

**Files:**
- Modify: `src/Variables.jl:1-156`

- [ ] **Step 2.1: 改写为元组驱动**

读 `Variables.jl:1-156`，按维度分类（标量/向量/矩阵），改成：

```julia
function StandardVariables(case::Case, num::Int64)
    variables = Dict{String, Union{Array{Float64},Float64}}()
    # 标量键（值为 0.0）
    for k in VariableKeys.SCALAR_GLOBAL_KEYS
        variables[k] = 0.0
    end
    # 向量键（值为 zeros(num)）
    for k in VariableKeys.VECTOR_GLOBAL_KEYS
        variables[k] = zeros(Float64, num)
    end
    # 矩阵键（按 case 维度）
    for (k, n) in VariableKeys.matrix_keys_with_dim(case)
        variables[k] = zeros(Float64, n, num)
    end
    return variables
end
```

**注**：Task 2 中如果 `VariableKeys` 不区分 SCALAR/VECTOR/MATRIX，需要先回到 `VariableKeys.md` 补充分类。

- [ ] **Step 2.2: 跑 characterization test**

Run: `julia test/unit_variables_keys.jl`
Expected: PASS

- [ ] **Step 2.3: 跑 example 验证**

Run: `julia example/minimal_example.jl`
Expected: 完整运行无错

- [ ] **Step 2.4: Commit**

```bash
git add src/Variables.jl
git commit -m "refactor(Variables): StandardVariables 改用 VariableKeys 常量元组"
```

### Task 3: 重构 `create_element_workspace`

**Files:**
- Modify: `src/Variables.jl:157-235`

- [ ] **Step 3.1: 改写为元组驱动**

```julia
function create_element_workspace(case::Case)
    ws = Dict{String, Union{Array{Float64},Float64}}()
    for k in VariableKeys.ELEMENT_WORKSPACE_KEYS
        # 按维度规则初始化（标量 0.0 / 向量 zeros(N) / 矩阵 zeros(N,Ngp)）
        ws[k] = VariableKeys.init_workspace_value(k, case)
    end
    return ws
end
```

- [ ] **Step 3.2: 跑测试 + example**

Run: `julia test/unit_variables_keys.jl && julia example/minimal_example.jl`
Expected: PASS

- [ ] **Step 3.3: Commit**

```bash
git add src/Variables.jl
git commit -m "refactor(Variables): create_element_workspace 改用 VariableKeys"
```

### Task 4: 更新 baseline

- [ ] **Step 4.1: 记录行数变化**

```bash
wc -l src/Variables.jl >> Simplify/baseline.md
```

Expected: 行数从 277 减少（预计 ~150-180）

---

## Validation

- [ ] `test/unit_variables_keys.jl` 通过
- [ ] `example/minimal_example.jl`、`example/SPMe_Thermal_example.jl` 跑通
- [ ] `Variables.jl` 不再含 `"negative particle lithium concentration"` 等键字面量（grep 验证）

## Risk

- **维度分类错误**：标量/向量/矩阵键若归类错，运行时形状不匹配。**缓解**：characterization test 先行；每个 example 跑通
- **键遗漏**：元组漏写一个键，下游访问返回 KeyError。**缓解**：Step 4.1 grep 验证
