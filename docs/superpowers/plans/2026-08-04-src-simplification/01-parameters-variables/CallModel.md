# CallModel.jl Simplification Plan

**Status:** ⬜ Pending
**Layer:** 1 - 参数与变量
**桶:** Consolidate（D2 成员）
**依赖:** `VariableKeys.md`、`Variables.md` 必须先完成
**阻塞:** 下游 MultiSPMe 调度

**Goal:** 删除 `copy_element_results` (line 247-268) 中 16 个键字面量，改用 `VariableKeys.ELEMENT_RESULT_KEYS` 元组循环构造。

**Architecture:** `copy_element_results(vars_e)` 改为 `Dict(k => vars_e[k] for k in VariableKeys.ELEMENT_RESULT_KEYS)`，一行代替 16 行字面量。

---

## Files

- Modify: `src/CallModel.jl:247-268`
- Reference: `src/VariableKeys.jl`

---

## Tasks

### Task 1: Characterization Test

**Files:**
- Create: `test/unit_callmodel_copy.jl`

- [ ] **Step 1.1: 写测试**

```julia
include("../src/JuBat.jl")
using .JuBat

@testset "copy_element_results 键集合稳定" begin
    # 构造一个有完整键的 vars_e
    keys = VariableKeys.ELEMENT_RESULT_KEYS
    vars_e = Dict{String, Union{Array{Float64},Float64}}(k => (length(k) % 2 == 0 ? zeros(2) : 0.0) for k in keys)
    copied = JuBat.CallModel.copy_element_results(vars_e)
    @test Set(keys(copied)) == Set(keys)
    # 引用语义验证：值应共享
    @test copied["cell voltage"] === vars_e["cell voltage"]
end
```

- [ ] **Step 1.2: 跑测试验证当前实现通过**

Run: `julia test/unit_callmodel_copy.jl`
Expected: PASS

- [ ] **Step 1.3: Commit**

```bash
git add test/unit_callmodel_copy.jl
git commit -m "test(CallModel): copy_element_results 键集合快照"
```

### Task 2: 重构 `copy_element_results`

**Files:**
- Modify: `src/CallModel.jl:247-268`

- [ ] **Step 2.1: 替换 16 行字面量为元组循环**

```julia
function copy_element_results(vars_e)
    Dict{String, Union{Array{Float64},Float64}}(
        k => vars_e[k] for k in VariableKeys.ELEMENT_RESULT_KEYS
    )
end
```

- [ ] **Step 2.2: 跑 characterization test + example**

Run: `julia test/unit_callmodel_copy.jl && julia example/minimal_example.jl`
Expected: PASS

- [ ] **Step 2.3: Commit**

```bash
git add src/CallModel.jl
git commit -m "refactor(CallModel): copy_element_results 改用 VariableKeys 元组"
```

---

## Validation

- [ ] `CallModel.jl` 不再含 `vars_e["..."]` 字面量（grep 验证）
- [ ] `example/czm_cycle_example.jl` 跑通（覆盖 MultiSPMe 路径）

## Risk

- 低；引用语义保持一致（值共享，不复制）
