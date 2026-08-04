# czm.jl Simplification Plan

**Status:** ⬜ Pending | **Layer:** 3 物理模型 | **桶:** Mixed（Consolidate + Dead code）

**Goal:** D9 triage `assemble_coupled_system` vs `_full`；确认 `build_czm_cache`/`ensure_czm_cache` 单一入口；不动 `assemble_czm_system`/`assemble_bulk_stiffness` 数值核心。

## 现状（873 行）

| 函数 | 行号 | 用途 | 桶 |
|---|---|---|---|
| `create_czm_mesh` | — | CZM 网格生成 | Leave alone |
| `assemble_czm_system` | 261-454 | 内聚力组装核心 | High-risk-leave-alone |
| `assemble_bulk_stiffness` | — | 体刚度 | High-risk-leave-alone |
| `build_czm_cache` | 621 | 缓存构造 | Leave alone |
| `ensure_czm_cache` | 719 | 缓存入口（守卫失效） | Leave alone |
| `assemble_coupled_system` | 737 | 体+内聚力组装 | **Keep（active）** |
| `assemble_coupled_system_full` | 777 | 含热-化学载荷 | **Dead code 候选** |
| `assemble_thermal_chemical_load` | 395-449 | 热-化学载荷 | Leave alone |

---

## Tasks

### Task 1: D9 triage — `assemble_coupled_system_full` 死代码验证

**Files:**
- Modify: `src/czm.jl:777-820`（删除函数）
- Modify: `src/JuBat.jl:69`（移除 export）

- [ ] **Step 1: grep 全仓确认无调用**

```bash
grep -rn "assemble_coupled_system_full" --include="*.jl" src/ example/ test/
```

预期：仅 `src/czm.jl:777` 定义、`src/JuBat.jl:69` export、其余命中只在 `.md` 文档。
若 `example/` 或 `test/` 出现调用 → 立即停止，改为 Leave alone。

- [ ] **Step 2: 删除函数定义**

删除 `src/czm.jl:771-820`（含 docstring 与函数体）。

- [ ] **Step 3: 移除 export**

`src/JuBat.jl:69`：从 export 行删除 `assemble_coupled_system_full`，保留 `assemble_thermal_chemical_load`。

- [ ] **Step 4: 跑 unit test**

Run: `julia test/unit_czm_bilinear.jl`
Expected: PASS

- [ ] **Step 5: 跑主线 example**

Run: `julia example/czm_cycle_example.jl`
Expected: PASS（应不再依赖 `_full`）

- [ ] **Step 6: Commit**

```bash
git add src/czm.jl src/JuBat.jl
git commit -m "chore(czm): 删除无调用者的 assemble_coupled_system_full"
echo "$(date +%F): czm.jl 决策=[删除 assemble_coupled_system_full]" >> Simplify/baseline.md
```

---

### Task 2: 审查 `build_czm_cache` / `ensure_czm_cache` 单一入口

**Files:** 仅审查 `src/czm.jl:604-728`

- [ ] 通读 `ensure_czm_cache` 守卫逻辑（line 719-728），确认 `cache.valid` 失效判据完整：
  - `cache === nothing`
  - `!cache.valid`
  - `cache.czm_mesh_id != objectid(czm_mesh)`
  - `cache.param_cache_id != param_cache.id`
  - `cache.fix_inner != fix_inner`
- [ ] grep 确认外部仅通过 `ensure_czm_cache` 进入：

```bash
grep -rn "build_czm_cache\b" --include="*.jl" src/
```

预期：仅 `czm.jl:621` 定义、`czm.jl:604`（docstring 示例）、`czm.jl:725`（ensure 内部调用）、`JuBat.jl:72` export。
若 example/test 直接调用 `build_czm_cache` → 保留 export；否则保留 export 不变（公开 API）。

- [ ] baseline 记录：

```bash
echo "$(date +%F): czm.jl build_czm_cache/ensure_czm_cache 审查通过，桶=Leave alone" >> Simplify/baseline.md
```

---

### Task 3: 数值核心不动（High-risk-leave-alone）

**Files:** 仅审查 `src/czm.jl:261-454`（`assemble_czm_system`）

- [ ] 通读 194 行，记录是否还存在以下问题：
  - try/catch 兜底（应为 0 处，spec §7 扫描已确认）
  - 静默回退（fallback / silent degrade）
  - 重复几何计算（已通过 `geom_cache` 命中）
- [ ] 无修改；若发现问题，单独开新 plan 评估，不在本简化轮处理。

---

## Risk

| 动作 | 风险 | 缓解 |
|---|---|---|
| 删 `_full` | 低（grep 是可靠判据） | Step 1 必须先跑 |
| 审查 cache | 无（只读） | — |
| 审查 assemble_czm_system | 无（只读） | — |

## 不做的事

- 不动 `assemble_czm_system` 内部（High-risk-leave-alone）
- 不合并 `assemble_coupled_system` 与 `_full`（_full 直接删，不合并）
- 不重命名任何函数（数值核心稳定）
