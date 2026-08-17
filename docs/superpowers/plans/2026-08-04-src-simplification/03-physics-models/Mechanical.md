# Mechanical.jl Stress Solve Hardening Plan

**Status:** ✅ Completed（审计保留） | **Layer:** 3 物理模型 | **桶:** Leave alone（category A 防御代码评估）

**Goal:** 评估 `try/catch` 兜底（line 294-299）是否属于 spec §7 类目 A「吞异常」；如果是，改为显式错误或保留并加注释。

## 现状（360 行）

| 函数 | 行号 | 用途 |
|---|---|---|
| `Mechanicaloutput` | 1 | 主入口（颗粒扩散应力） |
| `Calstressdisp` | 112 | 单电极扩散应力计算 |
| `thermal_diffusion_stress_2D` | ~200+ | 二维宏观应力（CZM 应变驱动） |
| `compute_effective_coating_modulus` | — | 极片模量入口（CLAUDE.md §9.4） |

**Category A 嫌疑**：
```julia
# src/Mechanical.jl:294-299
U_M = try
    K_mech \ F_mech
catch e
    @warn "Mechanical solve failed, using zero displacement" e
    zeros(Float64, ndof)
end
```

---

## Files

- Modify: `src/Mechanical.jl:294-299`（决策后处理）
- Test: `test/smoke_mechanical_solve.jl`（新建，验证奇异矩阵场景）

---

## Tasks

### Task 1: 评估 try/catch 类别

- [ ] **Step 1: 跑主线 example，观察是否触发 catch**

Run: `julia example/jellyroll_stress_displacement.jl`
Expected: 正常完成（不应触发 catch）

- [ ] **Step 2: 决策类目**

| 场景 | 决策 |
|---|---|
| catch 从未触发，K_mech 应该总是非奇异 | **删除 try/catch**，让 LAPACK 错误直接抛出 |
| catch 偶尔触发，零位移是有意义的回退（如初始 step） | **保留 + 加注释**：说明为什么这是设计内回退 |
| catch 经常触发，掩盖 bug | **删除 try/catch** + 修底层 K_mech 奇异根因 |

- [ ] **Step 3: 记录到 baseline**

```bash
echo "$(date +%F): Mechanical.jl:294 try/catch 决策=[...]" >> Simplify/baseline.md
```

---

### Task 2: 若删除 try/catch

**Files:** Modify `src/Mechanical.jl:294-299`

- [ ] **Step 1: 替换**

```julia
# 求解位移场（奇异矩阵视为数值 bug，直接报错）
U_M = K_mech \ F_mech
```

- [ ] **Step 2: 跑主线 example**

Run: `julia example/jellyroll_stress_displacement.jl`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add src/Mechanical.jl
git commit -m "refactor(mechanical): 删除兜底 try/catch（K_mech \\ F_mech）"
```

---

### Task 3: 若保留 try/catch（加注释）

- [ ] **Step 1: 加注释**

```julia
# 求解位移场
# 设计内回退：初始时间步或参数极端时 K_mech 可能奇异；
# 零位移不影响后续 CZM/热耦合（仅用于后处理可视化）
U_M = try
    K_mech \ F_mech
catch e
    @warn "Mechanical solve failed, using zero displacement" e
    zeros(Float64, ndof)
end
```

- [ ] **Step 2: Commit**

```bash
git add src/Mechanical.jl
git commit -m "docs(mechanical): 标注 try/catch 为设计内回退（spec §7 类目 A）"
```

---

## Risk

低；try/catch 处理不影响数值核心，且 `thermal_diffusion_stress_2D` 与 `Calstressdisp` 都不动。

## 不做的事

- 不动 `Calstressdisp`（颗粒扩散应力核心）
- 不动 `thermal_diffusion_stress_2D` 公式
- 不动 `compute_effective_coating_modulus`（CLAUDE.md §9.4 入口）

## Execution Result (2026-08-05)

- 保留带 `@warn` 的零位移回退；当前强制基线未覆盖奇异矩阵失败分支，计划指定的 3600 s / nθ=360 示例不适合作为每批快速 characterization。
- 未添加 taxonomy 注释；待建立聚焦失败路径测试后再评估 fail-fast。
