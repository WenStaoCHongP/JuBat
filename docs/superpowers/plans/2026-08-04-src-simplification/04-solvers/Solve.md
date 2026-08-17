# Solve.jl Defensive Code Audit Plan

**Status:** ✅ Completed（按评审修订） | **Layer:** 4 求解器 | **桶:** Mixed（Leave alone + category A/B 评估）

**Goal:** 主求解器入口。评估 2 处 try/catch（line 273, 418）+ 2 处 @warn（line 123, 148）；数值主循环不动。

## 现状（471 行）

| 区块 | 行号 | 用途 | 桶 |
|---|---|---|---|
| 主循环入口 | — | 时间步进、自适应 | High-risk-leave-alone |
| CZM 调用 try/catch | 273-312 | `update_czm_damage!` 容错 | 评估 |
| 热数据 try/catch | 418-426 | 「non-fatal」silent degrade | **category A 嫌疑** |
| @warn 状态长度不匹配 | 123 | 外部 y0 长度校验 | 评估（category B） |
| @warn 时间步超限 | 148-149 | 性能提示 | Leave alone（有用的提示） |

---

## Tasks

### Task 1: 评估 line 418-426 silent catch

**Files:** Modify `src/Solve.jl:418-426`

- [ ] **Step 1: 阅读上下文**

```julia
# src/Solve.jl:418-426
try
    if case.opt.per_element_spme && case.opt.thermalmodel == "distributed2D" && !isempty(T_nodes_carry)
        Tref = case.param_dim.scale.T_ref
        result["thermal2D final temperature at nodes [K]"] = T_nodes_carry .* Tref
        result["thermal2D nodes xy [m]"] = case.mesh["thermal2D"].node
    end
catch
    # non-fatal
end
```

**问题**：条件已在 if 内做防御，try/catch 多余。`T_nodes_carry` 在主循环已确保存在，无抛异常路径。

- [ ] **Step 2: 删除 try/catch**

替换为：

```julia
if case.opt.per_element_spme && case.opt.thermalmodel == "distributed2D" && !isempty(T_nodes_carry)
    Tref = case.param_dim.scale.T_ref
    result["thermal2D final temperature at nodes [K]"] = T_nodes_carry .* Tref
    result["thermal2D nodes xy [m]"] = case.mesh["thermal2D"].node
end
```

- [ ] **Step 3: 跑主线 example**

Run: `julia example/SPMe_Thermal_example.jl`
Expected: PASS（result 字典含 `thermal2D final temperature at nodes [K]`）

- [ ] **Step 4: Commit**

```bash
git add src/Solve.jl
git commit -m "refactor(solve): 删除 line 418 silent try/catch（条件已防御）"
```

---

### Task 2: 评估 line 273-312 CZM 调用 try/catch

**Files:** 仅审查 `src/Solve.jl:273-312`

- [ ] 通读，确认这是**设计内回退**（spec §7 类目 D「数值核心中的恢复路径」）：
  - CZM 失效不应导致电化学仿真崩溃
  - 已有 `czm_converged` 标志位记录失败
  - 失败时 `@debug` 记录（不是 silent）
- [ ] **保留 + 加注释**：

```julia
# CZM 损伤更新（设计内容错：CZM 不收敛时保留上一时刻损伤，仿真继续）
try
    u_czm_new, czm_converged = update_czm_damage!(...)
    ...
catch e
    @debug "CZM damage update failed at step $czm_step_count: $e"
    czm_converged = false
end
```

- [ ] Commit 注释：

```bash
git add src/Solve.jl
git commit -m "docs(solve): 标注 CZM try/catch 为设计内回退（spec §7 类目 D）"
```

---

### Task 3: 评估 line 123 @warn（外部 y0 长度不匹配）

**Files:** 仅审查 `src/Solve.jl:120-130`

- [ ] 通读，确认 `@warn` 后是否仍执行（fallback 到 `ModelInitialisation`）。
- [ ] 决策：
  - 若 fallback 是合理行为 → 保留 + 加注释「外部 y0 长度与多 SPMe 布局不匹配时回退到模型初始化」
  - 若应直接报错（用户传入错误 y0） → 改为 `@error` 或 `throw(DimensionMismatch(...))`
- [ ] baseline 记录决策

---

## Risk

| 动作 | 风险 | 缓解 |
|---|---|---|
| 删 line 418 try/catch | 低（条件已防御） | 主线 example 跑通 |
| 保留 line 273 try/catch | 无（加注释） | — |
| @warn line 123 决策 | 低 | 仅注释/级别变更 |

## 不做的事

- 不动主循环时间步进
- 不动 `CallModel_MultiSPMe` 数值路径
- 不重命名任何变量

## Execution Result (2026-08-05)

- 删除附加最终热节点结果时的 silent try/catch（−4 行）；现有条件守卫保持不变。
- 保留 CZM 损伤更新的可观察恢复路径，以及外部状态长度不匹配时的显式告警/重初始化。
- 未添加 spec taxonomy 注释，避免无行为价值的生产注释。
- 热边界 smoke 与强制 `testexample` 基线均通过，PNG SHA-256 完全一致。
