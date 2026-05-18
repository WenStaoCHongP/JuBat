# CZM 求解器统一设计规格

> 日期: 2026-05-18
> 状态: proposed
> 关联: `docs/planning-with-files/代码简化/`

## 1. 问题定义

`src/CzmSolve.jl` (664行) 包含三个独立求解器和若干辅助函数，彼此有约 60% 的重复代码：

| 函数 | 行数 | 核心逻辑 |
|------|------|----------|
| `solve_czm_basic_step` | 104 | 全量载荷 + 回溯线搜索 |
| `newton_raphson_czm` | 170 | 自适应载荷子步 + 内联线搜索 |
| `solve_czm_arc_length_step` | 197 | Crisfield 圆柱弧长法 |

### 重复模式

1. **参数/缓存/BC 提取** (~15行) — 三个函数完全一致
2. **热化学载荷计算** (~3行) — 完全一致
3. **事后组装+残差计算** (~12行) — 完全一致
4. **结果填充** (damage/separation/traction → result) — 三处重复
5. **线搜索** — `backtrack_line_search!` 和 `newton_raphson_czm` 内联实现逻辑相同
6. **BC 自由度提取** — `extract_bc_dofs` 与 `build_czm_cache` 内部的 BC 提取重复

### 根因

三个求解器共享相同的 Newton-Raphson 骨架，但各自的迭代策略被"复制粘贴"进了独立函数。差异仅在：
- 载荷推进方式（全量 / 自适应子步 / 弧长）
- 修正步计算（KΔu=R / 弧长约束二次方程）
- 收敛准则（纯残差 / 残差+弧长约数）

## 2. 设计目标

1. **消除求解器重复代码**：将三个求解器合并为一个统一框架，通过策略参数切换
2. **保持三种方法的独立行为**：basic / load_substep / arc_length 的计算路径不变
3. **`solve_czm_step` 公开 API 不变**：外部调用方无需修改
4. **消除辅助函数重复**：BC 提取、结果填充、线搜索各只留一份
5. **不改变数值结果**：重构前后计算结果 bit-exact 一致

## 3. 设计方案

### 3.1 统一求解器架构

```
solve_czm_step (公开 API，不变)
    │
    └── solve_czm_newton (统一入口)
            │
            ├── 公共前处理 (一次性)
            │   • 提取 ndof, n_coh, cache, bc_dofs/bc_vals
            │   • 计算 F_thermo_chem_total
            │   • 初始化 result, u, damage_states
            │
            ├── 迭代循环 (策略路由)
            │   ├── basic 策略: 单步 Newton + backtrack_line_search!
            │   ├── load_substep 策略: 自适应子步 + 线搜索
            │   └── arc_length 策略: Crisfield 弧长法
            │
            └── 公共后处理 (一次性)
                • 最终装配 + 残差计算
                • 填充 result (displacement, damage, separation, traction)
                • 构建 new_czm_mesh
```

### 3.2 策略分发

在统一函数内部，通过 `iter_method` 参数选择执行路径：

```julia
function solve_czm_newton(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param, u_prev; ...)
    # === 公共前处理 ===
    ndof, n_coh = 2 * czm_mesh.nnode, czm_mesh.n_cohesive
    result = CZMResult(ndof, n_coh)
    u = copy(u_prev)
    damage_states = czm_mesh.damage_states

    bc_dofs, bc_vals = extract_bc_dofs(czm_mesh, param; cache=cache)
    F_thermo_chem_total = assemble_thermal_chemical_load(...)

    # 提取缓存
    K_bulk_cached, geom_cache, ws = ...

    # === 策略路由 ===
    if method == "basic"
        converged, total_iter, R_final = newton_basic!(...)
    elseif method == "load_substep"
        converged, total_iter, R_final, load_progress = newton_load_substep!(...)
    elseif method == "arc_length"
        converged, total_iter, R_final, load_progress = newton_arc_length!(...)
    end

    # === 公共后处理 ===
    # 最终装配 + 残差
    K_total, f_int_total, separations, tractions = assemble_coupled_system(...)
    R = F_ext + F_thermo_chem - f_int_total
    apply_czm_dirichlet_bc!(R, bc_dofs, bc_vals)
    R_norm = norm(R)

    # 填充 result
    result.converged = converged
    result.iterations = total_iter
    result.residual_norm = ...
    result.displacement = u
    fill_czm_result!(result, u, damage_states, separations, tractions)

    new_czm_mesh = clone_czm_mesh_with_damage(czm_mesh, damage_states)
    return result, new_czm_mesh
end
```

### 3.3 三个策略函数

三个策略函数各自只包含迭代核心逻辑：

- **`newton_basic!`** (~50行)：单步 Newton 循环 + `backtrack_line_search!` 调用
- **`newton_load_substep!`** (~70行)：自适应子步循环 + 内联线搜索（当前 `newton_raphson_czm` 逻辑）
- **`newton_arc_length!`** (~100行)：Crisfield 弧长法子步循环（当前 `solve_czm_arc_length_step` 逻辑）

所有策略函数共享相同的签名约定：接收预提取的参数，通过修改传入的 `u`、`damage_states`、`result` 来返回状态。

### 3.4 辅助函数统一

| 重复项 | 当前状态 | 重构后 |
|--------|----------|--------|
| BC 提取 | `extract_bc_dofs` + `build_czm_cache` 内重复 | 只保留 `extract_bc_dofs`，`build_czm_cache` 调用它 |
| 结果填充 | `fill_czm_result!` + `newton_raphson_czm` 末尾手动填充 | 统一使用 `fill_czm_result!` |
| 线搜索 | `backtrack_line_search!` + `newton_raphson_czm` 内联 | 统一使用 `backtrack_line_search!` |
| 收敛后状态更新 | 三处重复 `update_damage(...)` | 统一在后处理中调用 |

## 4. 影响范围

### 修改文件

| 文件 | 改动 |
|------|------|
| `src/CzmSolve.jl` | 主要改动：合并为统一求解器 + 3 个策略函数 + 复用辅助函数 |
| `src/czm.jl` | 微调：`build_czm_cache` 调用 `extract_bc_dofs` 消除 BC 提取重复 |

### 不改文件

| 文件 | 原因 |
|------|------|
| `src/CouplingState.jl` | `update_czm_damage!` 调用 `solve_czm_step`，API 不变 |
| `src/Materialmatrix.jl` | 本构律不改 |
| `src/CzmPostProcess.jl` | 后处理不改 |
| `src/Mechanical.jl` | 不属于 CZM 求解器 |
| `src/SetCase.jl`, `src/Option.jl` | 配置不变 |

### API 兼容性

- `solve_czm_step(czm_mesh, F_ext, E_eff, ν_eff, cohesive_params, param, u_prev; ...)` — **签名和行为完全不变**
- `newton_raphson_czm(...)` — 内部函数，改为委托给统一入口
- `solve_czm_basic_step(...)` — 内部函数，改为委托给统一入口
- `solve_czm_arc_length_step(...)` — 内部函数，改为委托给统一入口

## 5. 预期收益

| 指标 | 当前 | 重构后 |
|------|------|--------|
| `CzmSolve.jl` 行数 | ~664 | ~400 |
| 求解器函数 | 4 个 (3 solver + 1 dispatch) | 2 个 (1 dispatch + 1 unified) |
| 策略函数 | 内联 | 3 个轻量级策略函数 |
| 结果填充点 | 4 处 | 1 处 |
| BC 提取实现 | 2 处 | 1 处 |
| 线搜索实现 | 2 处 | 1 处 |

## 6. 验证方案

### 6.1 单元验证

运行现有 CZM 探针脚本，对比重构前后的残差历史和收敛行为：

```
tools/czm_baseline_probe.jl          — basic 方法基线
tools/czm_convergence_diag.jl        — load_substep 收敛诊断
tools/verify_czm_standalone.jl       — 独立 CZM 验证
```

### 6.2 集成验证

```
example/testexample.jl               — 全耦合仿真
example/coupled_czm_thermal_example.jl — 耦合 CZM-热仿真
```

### 6.3 验收标准

- 三个验证探针输出与重构前一致（残差范数、最终 D_max、迭代次数）
- `testexample.jl` 跑通，无 NaN、无 `CZM solve issue` 警告
- 三种方法（basic / load_substep / arc_length）均可正常收敛

## 7. 风险与缓解

| 风险 | 缓解 |
|------|------|
| 合并后引入回归 | 先写策略函数，再逐步替换，每步验证探针 |
| 弧长法逻辑复杂易出错 | 策略函数直接复制原逻辑，不做逻辑改动，只删除重复外围代码 |
| 对耦合入口有隐含影响 | `solve_czm_step` API 不变，调用方无需修改 |

## 8. 不做的优化

以下改进不在本次范围内，避免范围蔓延：

- **不重组文件结构**：CZM 本构律仍留在 `Materialmatrix.jl`
- **不改文件命名**：不拆分/重命名 `czm.jl` / `CzmSolve.jl`
- **不向量化单元循环**：`assemble_czm_system` 的逐单元循环保持现状
- **不优化内存分配**：`CZMAssemblyWorkspace` 的预分配模式不变
