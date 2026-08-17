# Czm.jl Simplification Plan

**Status:** ✅ Completed（审计保留） | **Layer:** 3 物理模型 | **桶:** Mixed（Consolidate + Dead code）

**Goal:** D9 triage `assemble_coupled_system` vs `_full`；确认 `build_czm_cache`/`ensure_czm_cache` 单一入口；不动 `assemble_czm_system`/`assemble_bulk_stiffness` 数值核心。

## 现状（已拆分：`Czm.jl` 629 行、`CzmMesh.jl` 182 行、`CzmBC.jl` 63 行）

| 函数 | 行号 | 用途 | 桶 |
|---|---|---|---|
| `create_czm_mesh` | `CzmMesh.jl:33` | CZM 网格生成 | Leave alone |
| `assemble_czm_system` | `Czm.jl:80-262` | 内聚力组装核心 | High-risk-leave-alone |
| `assemble_bulk_stiffness` | `Czm.jl:273-343` | 体刚度 | High-risk-leave-alone |
| `build_czm_cache` | `Czm.jl:440` | 缓存构造 | Leave alone |
| `ensure_czm_cache` | `Czm.jl:538` | 缓存入口（守卫失效） | Leave alone |
| `assemble_coupled_system` | `Czm.jl:556` | 体+内聚力组装 | **Keep（active）** |
| `assemble_coupled_system_full` | `Czm.jl:596` | 含热-化学载荷 | **Dead code 候选** |
| `assemble_thermal_chemical_load` | `Czm.jl:352-416` | 热-化学载荷 | Leave alone |

---

## Tasks

### Task 1: D9 triage — `assemble_coupled_system_full` 死代码验证

**Files:**
- Modify: `src/Czm.jl:590-629`（删除函数及 docstring）
- Modify: `src/JuBat.jl:72`（移除 export）

- [ ] **Step 1: grep 全仓确认无调用**

```bash
grep -rn "assemble_coupled_system_full" --include="*.jl" src/ example/ test/
```

预期：仅 `src/Czm.jl:596` 定义、`src/JuBat.jl:72` export、其余命中只在 `.md` 文档。
若 `example/` 或 `test/` 出现调用 → 立即停止，改为 Leave alone。

- [ ] **Step 2: 删除函数定义**

删除 `src/Czm.jl:590-629`（含 docstring 与函数体）。

- [ ] **Step 3: 移除 export**

`src/JuBat.jl:72`：从 export 行删除 `assemble_coupled_system_full`，保留 `assemble_thermal_chemical_load`。

- [ ] **Step 4: 跑 unit test**

Run: `julia test/unit_czm_bilinear.jl`
Expected: PASS

- [ ] **Step 5: 跑主线 example**

Run: `julia example/czm_cycle_example.jl`
Expected: PASS（应不再依赖 `_full`）

- [ ] **Step 6: Commit**

```bash
git add src/Czm.jl src/JuBat.jl
git commit -m "chore(czm): 删除无调用者的 assemble_coupled_system_full"
echo "$(date +%F): Czm.jl 决策=[删除 assemble_coupled_system_full]" >> Simplify/baseline.md
```

---

### Task 2: 审查 `build_czm_cache` / `ensure_czm_cache` 单一入口

**Files:** 仅审查 `src/Czm.jl:421-548`

- [ ] 通读 `ensure_czm_cache` 守卫逻辑（line 538-548），确认 `cache.valid` 失效判据完整：
  - `cache === nothing`
  - `!cache.valid`
  - `cache.czm_mesh_id != objectid(czm_mesh)`
  - `cache.param_cache_id != param_cache.id`
  - `cache.fix_inner != fix_inner`
- [ ] grep 确认外部仅通过 `ensure_czm_cache` 进入：

```bash
grep -rn "build_czm_cache\b" --include="*.jl" src/
```

预期：仅 `Czm.jl:440` 定义、`Czm.jl:421`（docstring 示例）、`Czm.jl:544`（ensure 内部调用）、`JuBat.jl:75` export。
若 example/test 直接调用 `build_czm_cache` → 保留 export；否则保留 export 不变（公开 API）。

- [ ] baseline 记录：

```bash
echo "$(date +%F): Czm.jl build_czm_cache/ensure_czm_cache 审查通过，桶=Leave alone" >> Simplify/baseline.md
```

---

### Task 3: 数值核心不动（High-risk-leave-alone）

**Files:** 仅审查 `src/Czm.jl:80-262`（`assemble_czm_system`）

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

## Execution Result (2026-08-05)

- `_full` 无仓库调用者但已由 `JuBat.jl` 明确导出，按公共 API 规则保留。
- `ensure_czm_cache` 是活跃内部入口；`build_czm_cache` 同样已导出，且缓存失效判据已有定向测试覆盖。
- 数值组装核心保持不动。
