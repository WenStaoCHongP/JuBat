# CouplingState.jl Consolidation Plan

**Status:** ⬜ Pending | **Layer:** 5 后处理/导出 | **桶:** Consolidate + Leave alone（D2 + category B）

**Goal:** D2 部分整合（`MultiSPMeLayout` double ctor、`update_czm_damage!` double method）+ category B `@warn` 标注。

## 现状（762 行）

| 函数/结构 | 行号 | 用途 | 桶 |
|---|---|---|---|
| `MultiSPMeLayout` ctor (5-arg) | 98 | 完整构造 | Keep（active） |
| `MultiSPMeLayout` ctor (轻量) | 109 | 内部辅助 | **D2 评估** |
| `compute_czm_params_per_interface` | 302 | CZM 参数 | Leave alone |
| `update_czm_damage!` (双方法) | 497, 613 | two signatures | **D9 triage** |
| `@warn` × 4 | 340, 532, 540, 593 | NaN/数值检查 | category B，保留 |

---

## Tasks

### Task 1: D9 triage `update_czm_damage!` 双方法

**Files:** Modify `src/CouplingState.jl:497, 613`

- [ ] **Step 1: 阅读两方法签名 + 调用者**

```bash
grep -rn "update_czm_damage!" --include="*.jl" src/
```

- [ ] **Step 2: 决策**

| 情况 | 决策 |
|---|---|
| 两签名都被外部调用（不同参数集） | 保留为 Julia 多方法（multi-method），加 docstring 区分 |
| 仅一个被调用，另一个是旧接口 | 删除未被调用者 |

- [ ] **Step 3: 若保留多方法 → 加 docstring**

```julia
"""
    update_czm_damage!(case, variables, T_nodes) -> (u_czm, converged)

主入口：从 thermal2D 节点温度 + 电化学 SOC 计算 CZM 损伤增量。
"""
function update_czm_damage!(case::Case, variables::Dict, T_nodes::Vector{Float64})
    ...
end

"""
    update_czm_damage!(case, variables, T_nodes, debug::Bool) -> ...

debug 入口（额外日志）。spec §7 类目 E（兼容入口）。
"""
function update_czm_damage!(case::Case, variables::Dict, T_nodes::Vector{Float64}, debug::Bool)
    ...
end
```

- [ ] **Step 4: 若删一个 → 删除未被调用者**

- [ ] **Step 5: 跑主线 example**

Run: `julia example/coupled_czm_thermal_example.jl`
Expected: PASS

- [ ] **Step 6: Commit + baseline**

```bash
git add src/CouplingState.jl
git commit -m "refactor(coupling): D9 update_czm_damage! 双方法决策=[...]"
echo "$(date +%F): CouplingState.jl update_czm_damage! 决策=[...]" >> Simplify/baseline.md
```

---

### Task 2: `MultiSPMeLayout` 双 ctor 评估

**Files:** 仅审查 `src/CouplingState.jl:98, 109`

- [ ] 通读两个 ctor，确认是否：
  - 一个是完整构造（外部调用）
  - 一个是内部辅助（被完整构造调用）
- [ ] 若是辅助 → 保留（不是重复，是分层）
- [ ] 若是相同功能不同名 → 合并

- [ ] baseline 记录

---

### Task 3: category B `@warn` 标注

**Files:** Modify `src/CouplingState.jl:340, 532, 540, 593`

- [ ] 每处 @warn 上方加注释：

```julia
# spec §7 类目 B：数值输入检查的运行时警告（非 silent swallow）
@warn "CZM inputs contain NaN" ...
```

- [ ] Commit：

```bash
git add src/CouplingState.jl
git commit -m "docs(coupling): 标注 4 处 @warn 为数值检查（spec §7 类目 B）"
```

---

## Risk

| 动作 | 风险 | 缓解 |
|---|---|---|
| `update_czm_damage!` 多方法合并 | 中（数值路径） | 必须先 grep 调用者 |
| `@warn` 注释 | 无 | — |

## 不做的事

- 不动 `compute_czm_params_per_interface`（CZM 参数核心）
- 不删任何 @warn（都是数值健康检查）
- 不动 `MultiSPMeLayout` struct 定义
