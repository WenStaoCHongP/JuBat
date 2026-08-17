# CycleData.jl D1 Refactor Plan

**Status:** ⚠️ Blocked（评审后待重写） | **Layer:** 5 后处理/导出 | **桶:** Consolidate

**Goal:** D1 重复消除。`solve_phase_with_export` / `solve_cycle_with_export` 重复 CycleSolver ~400 行。两版都有外部 example 调用者 → **不能直接删 CycleData**。改为：CycleSolver 增加 `export_callback` 关键字，CycleData 改为 thin wrapper。

## 现状（623 行）

| 函数 | 行号 | 用途 | 外部调用者 |
|---|---|---|---|
| `solve_phase_with_export` | 32 | 单相位 + 数据采集 | `example/电化学-热耦合验证/不同倍率温度曲线.jl:132` |
| `solve_cycle_with_export` | 268 | 完整循环 + 数据采集 | `example/循环验证/export_cycle_data_example.jl:134` |

**重复来源**：内部复刻了 `solve_phase`（CycleSolver.jl:118）的整个时间步进循环，只在每步追加 `TimeStepData` push。

**D9 triage 结果**：两版都有外部调用者 → **保留两版**，但 wrapper 化。

---

## Files

- Modify: `src/CycleSolver.jl:118`（`solve_phase` 签名增加 `export_callback`）
- Modify: `src/CycleData.jl:32-260`（`solve_phase_with_export` 改 wrapper）
- Modify: `src/CycleData.jl:268-450`（`solve_cycle_with_export` 改 wrapper）
- Test: `test/smoke_cycledata_wrapper.jl`（新建）

---

## Tasks

### Task 1: CycleSolver 增加 `export_callback` 关键字

**Files:** Modify `src/CycleSolver.jl:118`

- [ ] **Step 1: 写失败测试**

Create `test/smoke_cycledata_wrapper.jl`:

```julia
include("../src/JuBat.jl")
using .JuBat

# 构造 case
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.model = "SPMe"
case = JuBat.SetCase(param_dim, opt)

cycle_opt = JuBat.CycleOption(n_cycles=1, I_charge=5.0, I_discharge=5.0,
                              t_charge=300, t_discharge=300,
                              V_upper=4.2, V_lower=2.5, SOC_init=0.5)

# 收集器
collected = JuBat.TimeStepData[]
cb = (step_data::JuBat.TimeStepData, cycle::Int, phase::JuBat.PhaseType) -> push!(collected, step_data)

# 调用带 callback 的 solve_cycling
result = JuBat.solve_cycling(case, cycle_opt; export_callback=cb)
@assert !isempty(collected) "export_callback 未被触发"
println("PASS: export_callback 触发 $(length(collected)) 次")
```

- [ ] **Step 2: 跑测试，确认 FAIL**

Run: `julia test/smoke_cycledata_wrapper.jl`
Expected: FAIL（`export_callback` 关键字不存在）

- [ ] **Step 3: 修改 `solve_phase` 签名**

`src/CycleSolver.jl:118` 当前签名：
```julia
function solve_phase(case::Case, phase_type::PhaseType, t_max::Float64, ...)
```

改为追加关键字：
```julia
function solve_phase(case::Case, phase_type::PhaseType, t_max::Float64, ...;
                     export_callback::Union{Function, Nothing}=nothing)
```

- [ ] **Step 4: 在时间步循环内调用 callback**

定位 `solve_phase` 内每个成功时间步后的位置（results push 附近），插入：

```julia
if export_callback !== nothing
    step_data = TimeStepData(
        t * case.param.scale.t0, phase_type,
        variables["cell voltage [V]"][1],
        JuBat.extract_current(variables),  # 辅助函数，已存在则用之
        copy(variables["thermal2D temperature [K]"]),
        maximum(variables["thermal2D temperature [K]"]),
        mean(variables["thermal2D temperature [K]"]),
        copy(get(variables, "thermal2D element soc_NE", Float64[])),
        copy(get(variables, "thermal2D element soc_PE", Float64[])),
        mean(get(variables, "thermal2D element soc_NE", [0.0])))
    )
    export_callback(step_data, 0, phase_type)  # cycle=0 由 solve_cycling 覆盖
end
```

- [ ] **Step 5: `solve_cycling` 也加 `export_callback` 并透传 cycle index**

`src/CycleSolver.jl:217`：签名加 `export_callback::Union{Function, Nothing}=nothing`，
内部调用 `solve_phase` 时包装：

```julia
phase_cb = export_callback === nothing ? nothing :
    (sd, _, ph) -> export_callback(sd, cycle_idx, ph)
solve_phase(..., export_callback=phase_cb)
```

- [ ] **Step 6: 跑 smoke test，确认 PASS**

Run: `julia test/smoke_cycledata_wrapper.jl`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/CycleSolver.jl test/smoke_cycledata_wrapper.jl
git commit -m "feat(cycle): solve_phase/solve_cycling 增加 export_callback 关键字"
```

---

### Task 2: CycleData.jl 改 wrapper

**Files:** Modify `src/CycleData.jl:32-260, 268-450`

- [ ] **Step 1: `solve_phase_with_export` 改写**

将整个 line 32-260（~228 行）替换为 ~30 行 wrapper：

```julia
function solve_phase_with_export(case::Case, phase_type::PhaseType, t_max::Float64,
                                  I_current::Float64, V_limit::Float64,
                                  initial_state::Dict;
                                  czm_mesh=nothing, czm_params=nothing,
                                  dt_range::Vector{Float64}=[1.0, 10.0],
                                  export_interval::Int=1)
    timestep_data = TimeStepData[]
    step_counter = Ref(0)

    cb = (sd, cycle, ph) -> begin
        step_counter[] += 1
        if step_counter[] % export_interval == 0
            push!(timestep_data, sd)
        end
    end

    result = JuBat.solve_phase(case, phase_type, t_max, I_current, V_limit, initial_state;
                                czm_mesh=czm_mesh, dt_range=dt_range,
                                export_callback=cb)

    node_coords, element_connectivity = JuBat.extract_mesh_geometry(case)
    ne = size(element_connectivity, 1)
    nT = size(node_coords, 1)
    export_data = CycleExportData(0, timestep_data, node_coords, element_connectivity, ne, nT)
    return result, export_data
end
```

- [ ] **Step 2: `solve_cycle_with_export` 同样改写**

将 line 268-450（~180 行）替换为 ~40 行 wrapper：

```julia
function solve_cycle_with_export(case::Case, cycle_opt::CycleOption;
                                  verbose::Bool=true, export_interval::Int=1)
    timestep_data = TimeStepData[]
    step_counter = Ref(0)

    cb = (sd, cycle, ph) -> begin
        step_counter[] += 1
        if step_counter[] % export_interval == 0
            push!(timestep_data, sd)
        end
    end

    result = JuBat.solve_cycling(case, cycle_opt; verbose=verbose, export_callback=cb)

    node_coords, element_connectivity = JuBat.extract_mesh_geometry(case)
    ne = size(element_connectivity, 1)
    nT = size(node_coords, 1)
    export_data = CycleExportData(cycle_opt.n_cycles, timestep_data,
                                   node_coords, element_connectivity, ne, nT)
    return result, export_data
end
```

- [ ] **Step 3: 跑两个 example 验证**

Run: `julia example/循环验证/export_cycle_data_example.jl`
Expected: PASS

Run: `julia example/电化学-热耦合验证/不同倍率温度曲线.jl`
Expected: PASS

- [ ] **Step 4: 跑 smoke test**

Run: `julia test/smoke_cycledata_wrapper.jl`
Expected: PASS

- [ ] **Step 5: 行数验证**

```bash
wc -l src/CycleData.jl
```

预期：~623 → ~150-200（减少 ~400-450 行）。

- [ ] **Step 6: Commit + baseline**

```bash
git add src/CycleData.jl
git commit -m "refactor(cycledata): D1 重复消除，wrapper 化（-400 行）"
echo "$(date +%F): CycleData.jl 决策=[Consolidate: wrapper 化，-400 行]" >> Simplify/baseline.md
```

---

## Risk

| 动作 | 风险 | 缓解 |
|---|---|---|
| `solve_phase` 签名变更 | 中（影响 CycleSolver 内部多处调用） | 仅加 optional kw，向后兼容 |
| `TimeStepData` 字段映射 | 中（字段必须与 example 期望一致） | smoke test 验证 |
| 几何提取 `extract_mesh_geometry` | 低（如不存在则新增 helper） | 检查现有 API |

## 不做的事

- 不删 `solve_phase_with_export` / `solve_cycle_with_export`（有外部调用者）
- 不改 `TimeStepData` 字段定义（外部依赖）
- 不动 `solve_cycling` 的循环逻辑
