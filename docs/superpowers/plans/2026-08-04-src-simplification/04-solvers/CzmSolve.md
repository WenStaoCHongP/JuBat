# CzmSolve.jl Newton Recovery Path Annotation Plan

**Status:** ✅ Completed（审计保留） | **Layer:** 4 求解器 | **桶:** High-risk-leave-alone + category D 标注

**Goal:** 4 处 try/catch 全部保留（数值核心恢复路径），加注释说明为何不是 silent swallow。不动 Newton 迭代。

## 现状（675 行）

| 函数 | 行号 | 用途 |
|---|---|---|
| `solve_czm_step` | — | 单步 Newton（4 处 try/catch） |
| `solve_czm_arc_length` | — | 弧长法 |
| `solve_czm_adaptive_load` | — | 自适应载荷 |

**4 处 try/catch（spec §7 类目 D）**：

| 行号 | 上下文 | catch 行为 | 评估 |
|---|---|---|---|
| 220-224 | `K_bc \ R_bc` | `break`（退出 Newton，标 fail） | 保留 + 注释 |
| 319-323 | tangent 计算 | `break` | 保留 + 注释 |
| 372-376 | `delta_u_R` | `break` | 保留 + 注释 |
| 381-385 | `delta_u_F` | `break` | 保留 + 注释 |

**@warn（line 430, 603）**：arc-length stalled / adaptive stalled — 保留（设计内用户提示）。

---

## Tasks

### Task 1: 4 处 try/catch 加标准化注释

**Files:** Modify `src/CzmSolve.jl:220, 319, 372, 381`

- [ ] **Step 1: 模板注释**

每处 try/catch 上方加：

```julia
# spec §7 类目 D：Newton 迭代内的恢复路径
# 奇异矩阵时 break 退出迭代，外层标 czm_converged=false，不静默吞异常
Δu = try
    K_bc \ R_bc
catch
    break
end
```

- [ ] **Step 2: 跑 unit test**

Run: `julia test/unit_czm_newton.jl`
Expected: PASS

- [ ] **Step 3: 跑主线 example**

Run: `julia example/czm_cycle_example.jl`
Expected: PASS

- [ ] **Step 4: Commit + baseline**

```bash
git add src/CzmSolve.jl
git commit -m "docs(czmsolve): 标注 4 处 try/catch 为 Newton 恢复路径（spec §7 类目 D）"
echo "$(date +%F): CzmSolve.jl 4 处 try/catch 全部保留+注释，桶=High-risk-leave-alone" >> Simplify/baseline.md
```

---

## Risk

无；纯注释修改。

## 不做的事

- 不动 Newton 迭代算法
- 不动弧长法/自适应载荷算法
- 不合并 4 处 try/catch（数值稳定性需要）
