# ThermalDistributed.jl Twin Elimination Plan

**Status:** ✅ Completed | **Layer:** 3 物理模型 | **桶:** Consolidate

**Goal:** D3/D4 双胞胎消除：保留 `apply_convection_bc` 与 `apply_cool_method` 名称，直接采用原位变体实现，删除两个 `!` 变体。

## 现状（557 行）

| 函数对 | non-! 行号 | ! 行号 | 重复行数 | 处理 |
|---|---|---|---|---|
| `apply_convection_bc` / `apply_convection_bc!` | 49-96 | 178-212 | ~35 | non-! → wrapper |
| `apply_cool_method` / `apply_cool_method!` | 99-174 | 214-260 | ~45 | non-! → wrapper |
| `compute_heat_sources` / `_with_czm` | 385 | 518 | 不是双胞胎（_with_czm 调用 _） | Leave alone |

**D9 triage 结果**：
- non-! 版被 export（JuBat.jl:57-58），是公开 API
- ! 版仅 `ThermalDistributed2D_BC`（line 320-321）内部使用
- → **保留两版**，但 non-! 版改为 wrapper（消除逻辑重复）

---

## Files

- Modify: `src/ThermalDistributed.jl:49-96`（`apply_convection_bc` 改 wrapper）
- Modify: `src/ThermalDistributed.jl:99-174`（`apply_cool_method` 改 wrapper）
- Test: `test/smoke_thermal_bc.jl`（新建，验证两版行为等价）

---

## Tasks

### Task 1: 写 smoke test（先验等价性）

**Files:** Create `test/smoke_thermal_bc.jl`

- [x] **Step 1: 写测试**

```julia
include("../src/JuBat.jl")
using .JuBat

# 构造小 case
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=20, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)

mesh = case.mesh["thermal2D"]
nnode = size(mesh.node, 1)
KT = sparse(zeros(nnode, nnode))
FT = zeros(nnode)

# non-! 版（拷贝）
K1, F1 = JuBat.apply_convection_bc(KT, FT, mesh, nothing, case)
K2, F2 = JuBat.apply_cool_method(KT, FT, mesh, case)

# ! 版（需先 copy，因为原 Non-! 也是 copy 进入）
K3 = copy(KT); F3 = copy(FT)
JuBat.apply_convection_bc!(K3, F3, case;)
K4 = copy(KT); F4 = copy(FT)
JuBat.apply_cool_method!(K4, F4, mesh, case)

@assert K1 ≈ K3 "convection K mismatch"
@assert F1 ≈ F3 "convection F mismatch"
@assert K2 ≈ K4 "cool K mismatch"
@assert F2 ≈ F4 "cool F mismatch"
println("PASS: twin behavior equivalent")
```

- [x] **Step 2: 跑测试，验证当前等价**

Run: `julia test/smoke_thermal_bc.jl`
Expected: PASS（若 FAIL，说明两版已有偏差，先记录差异，本 plan 暂停）

- [ ] **Step 3: Commit baseline test**（工作树已有用户修改，未擅自提交）

```bash
git add test/smoke_thermal_bc.jl
git commit -m "test(thermal): 添加 convection/cool 双胞胎行为 baseline"
```

---

### Task 2: `apply_convection_bc` 改 wrapper

**Files:** Modify `src/ThermalDistributed.jl:49-96`

- [x] **Step 1: 替换函数体**

将 line 49-96 整段替换为：

```julia
function apply_convection_bc(KT, FT, mesh, is_outer, case; edge_cache=nothing)
    K = copy(KT)
    F = copy(FT)
    # 当 edge_cache===nothing 时，! 版内部会回退调用本函数；
    # 为避免无限递归，这里直接构造 edge_cache 再委托
    if edge_cache === nothing
        if is_outer === nothing
            _, is_outer = identify_boundary_nodes(mesh, case.param)
        end
        edge_cache = compute_boundary_edge_cache(mesh, is_outer)
    end
    apply_convection_bc!(K, F, case; edge_cache=edge_cache)
end
```

**注意**：当前 `apply_convection_bc!` 在 `edge_cache===nothing` 时（line 189-191）会回退调用 `apply_convection_bc`；wrapper 必须避免无限递归，所以 wrapper 内部主动构造 `edge_cache`。

- [x] **Step 2: 修改 `apply_convection_bc!` line 189-191**

将：
```julia
    if edge_cache === nothing
        return apply_convection_bc(K, F, case.mesh["thermal2D"], nothing, case)
    end
```
改为（自身主动构造，不再回调 non-! 版）：
```julia
    if edge_cache === nothing
        _, is_outer = identify_boundary_nodes(case.mesh["thermal2D"], case.param)
        edge_cache = compute_boundary_edge_cache(case.mesh["thermal2D"], is_outer)
    end
```

- [x] **Step 3: 跑 smoke test**

Run: `julia test/smoke_thermal_bc.jl`
Expected: PASS

- [ ] **Step 4: Commit**（工作树已有用户修改，未擅自提交）

```bash
git add src/ThermalDistributed.jl
git commit -m "refactor(thermal): apply_convection_bc non-! 改为 wrapper（消除 ~35 行重复）"
```

---

### Task 3: `apply_cool_method` 改 wrapper

**Files:** Modify `src/ThermalDistributed.jl:99-174`

- [x] **Step 1: 替换函数体**

将 line 99-174 替换为：

```julia
function apply_cool_method(KT, FT, mesh, case)
    K = copy(KT)
    F = copy(FT)
    apply_cool_method!(K, F, mesh, case)
end
```

（`apply_cool_method!` 已处理 `cool_method="none"/"surface"/"tab"` 全分支，且 non-! 版原行为也是 copy 后逐分支相同，无递归风险。）

- [x] **Step 2: 跑 smoke test**

Run: `julia test/smoke_thermal_bc.jl`
Expected: PASS

- [x] **Step 3: 跑主线 example**（以更强的 `example/testexample.jl` 全耦合基线替代并通过）

Run: `julia example/SPMe_Thermal_example.jl`
Expected: PASS

- [x] **Step 4: baseline**（已更新；commit 因用户工作树已有修改而跳过）

```bash
git add src/ThermalDistributed.jl
git commit -m "refactor(thermal): apply_cool_method non-! 改为 wrapper（消除 ~45 行重复）"
echo "$(date +%F): ThermalDistributed.jl 决策=[Consolidate: twin wrapper 化]" >> Simplify/baseline.md
```

---

### Task 4: 审查 `compute_heat_sources_with_czm`

**Files:** 仅审查 `src/ThermalDistributed.jl:518-557`

- [x] 通读，确认 `_with_czm` 是「调用 `_` + 追加 CZM 热源」的扩展，不是字面重复。
- [x] 追加部分具有独立 CZM 热源职责；本简化轮不动。

---

## Risk

| 动作 | 风险 | 缓解 |
|---|---|---|
| `apply_convection_bc!` 回调消除 | 中（递归风险） | Step 2 必须同步改 line 189-191 |
| wrapper 行为偏差 | 低 | smoke test 先验等价 |
| `cool_method="tab"` 分支 | 低 | smoke test + 主线 example 双验 |

## 不做的事

- 不合并 `compute_heat_sources` / `_with_czm`（不是字面重复）
- 不动 `ThermalDistributed2D_BC`（line 320-321 的调用点，行为不变）
- 不删 `edge_cache` 缓存机制

## Final Execution Result (2026-08-05)

- 按用户修订，不保留 wrapper/variant 双层：两个非 `!` 函数直接原位修改传入的 K/F，两个 `!` 定义已删除。
- `ThermalDistributed2D_BC` 已改为调用保留的函数名；`edge_cache` 机制保持不变。
- `smoke_thermal_bc.jl` 18/18，通过强制 `testexample` 基线且 PNG SHA-256 完全一致。
- 最终源码 diff：`+6/-133`，净删 127 行。
