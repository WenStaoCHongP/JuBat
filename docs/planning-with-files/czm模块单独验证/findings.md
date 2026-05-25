# Findings & Decisions

## Requirements
- 用户要的是 CZM 模块的独立验证计划，而不是立刻把电化学-热耦合接回去。
- 独立验证要聚焦内聚力子模块本身。
- 缺少的扩散应力和热应力要通过函数输入提供。
- 计划需要覆盖验证层级、接口边界、判据和后续落地路径。

## Research Findings
- `example/内聚力验证/czm_example.jl` 已经是一个不依赖电化学-热耦合的 CZM 纯力学验证入口。
- 该示例当前使用 `cohesive.czm_model = "mix"`，更适合混合模式和 BK 准则检查，不适合直接代表纯 Mode I 计划。
- `tools/verify_czm_unit.jl` 和 `tools/verify_czm_analytical.jl` 已经覆盖了本构解析对比、单元级和简单系统级验证，可作为独立验证基线。
- `src/CzmSolve.jl` 的 `solve_czm_step` 和 `newton_raphson_czm` 已经是 CZM 的直接入口，且保留了 `dT_elem`、`Δsoc_n_elem`、`Δsoc_p_elem` 等可选输入形态。
- `src/czm.jl` 中的 `assemble_thermal_chemical_load` 和 `assemble_coupled_system_full` 反映的是耦合路径；如果做独立验证，不应把它们作为必须依赖。
- `src/Materialmatrix.jl` 的本构已经区分 `model1` 和 `mix`，所以独立验证可以拆成纯 Mode I 和混合模式两条线。
- 已新增 `tools/verify_czm_standalone.jl`：脚本层用 provider callback 显式构造 `dT_elem / Δsoc_n_elem / Δsoc_p_elem`，并直接调用 `solve_czm_step`，无需接入电化学-热主链。
- 该脚本已完成可执行验证：Mode I 本构 max_abs_error=0，uniform/gradient 输入驱动求解均收敛，mix/BK 曲线已生成并写入 `output/czm_standalone/`。
- 该脚本的时间历史路径也已完成：热/扩散输入随时间变化，standalone 里记录了 D(t)，且 uniform/gradient 两个场景都出现了清晰的损伤演化。

## 求解器收敛性对比发现 (2026-05-12)

### 脚本重构
`tools/verify_czm_standalone.jl` 已重构为：**纯 Mode I 法向验证 + 三种求解器（basic / load_substep / arc_length）收敛性对比**。只输出汇总表，不画图。

### 归一化参数关键发现
- **必须使用 `case.param.cohesive`（归一化后参数）**，而非 `param_dim.cohesive`（物理参数）。
- 归一化后 `σ_max_n=1.0`, `K_n≈176`, `δ_0_n≈5.69e-3`, `E_eff≈210`——量级一致，系统可解。
- 使用物理参数时 `K_n=1.2e17` vs `E_eff=210`，系统极端病态，所有方法都不收敛。
- `czm_baseline_probe.jl` 也存在同样问题（使用 `param_dim.cohesive`），所有方法均不收敛。

### 收敛性对比结果（nθ=40, tol=1e-4, n_load_steps=50）

#### 初始结果（球面弧长约束）

| Δsoc_n | basic | load_substep | arc_length |
|--------|-------|-------------|------------|
| 0.1 | OK 2it D=0.00 r=4.6e-08 | OK 100it D=0.00 r=9.1e-10 | OK 50it D=0.00 r=4.6e-08 |
| 1.0 | OK 2it D=0.00 r=4.6e-07 | OK 100it D=0.00 r=9.1e-09 | OK 50it D=0.00 r=4.6e-07 |
| 1.5 | OK 3it D=0.77 r=2.8e-10 | FAIL D=0.77 r=0.13 | FAIL D=0.00 |
| 3.0 | OK 5it D=0.97 r=1.7e-02 | FAIL D=0.97 r=0.30 | FAIL D=0.00 |
| 5.0 | OK 8it D=1.00 r=6.6e-01 | FAIL D=1.00 r=0.84 | FAIL D=0.00 |
| 10.0 | FAIL 2it D=0.00 | FAIL D=1.00 r=0.26 | FAIL D=0.00 |

#### Crisfield 圆柱弧长法修复后（2026-05-12）

| Δsoc_n | basic | load_substep | arc_length (圆柱) |
|--------|-------|-------------|------------|
| 0.1 | OK 2it D=0.00 r=4.6e-08 | OK 100it D=0.00 r=9.1e-10 | OK 50it D=0.00 r=4.6e-08 |
| 0.5 | OK 2it D=0.00 r=2.3e-07 | OK 100it D=0.00 r=4.6e-09 | OK 50it D=0.00 r=2.3e-07 |
| 1.0 | OK 2it D=0.00 r=4.6e-07 | OK 100it D=0.00 r=9.1e-09 | OK 50it D=0.00 r=4.6e-07 |
| 1.5 | OK 3it D=0.77 r=2.8e-10 | FAIL 103it D=0.77 r=0.13 | FAIL 70it D=0.77 r=0.13 |
| 2.0 | OK 3it D=0.91 r=9.8e-10 | FAIL 106it D=0.91 r=0.082 | **OK** 77it D=0.91 r=4.1e-04 |
| 3.0 | OK 5it D=0.97 r=1.7e-02 | FAIL 115it D=0.97 r=0.30 | **OK** 90it D=0.97 r=1.9e-04 |
| 5.0 | OK 8it D=1.00 r=6.6e-01 | FAIL 132it D=1.00 r=0.84 | **OK** 99it D=1.00 r=8.7e-04 |
| 10.0 | FAIL 2it D=0.00 r=2.2e+04 | FAIL 148it D=1.00 r=0.26 | FAIL 112it D=1.00 r=4.1e+02 |

**arc_length 改进汇总**：
- 收敛率：3/8 → 6/8（+3 个损伤区级别）
- 损伤区 D 值：全为 0 → D=0.77~1.0（与 basic 一致）
- elastic regime 不受影响
- basic / load_substep 结果完全不变

### 分析
1. **elastic regime（Δsoc≤1）**：三种方法都收敛，basic 最快（2 it）。
2. **损伤 regime（Δsoc 1.5-5）**：
   - basic 收敛但残差较大（rel_norm 下降足够但绝对值不小）。
   - load_substep 残差始终无法降到 tol 以下，D 值接近 basic 但报 FAIL。
   - arc_length（球面）完全不工作（D=0），原因是弹性预测子低估位移，球面约束 hypersphere 与平衡路径不相交。
   - arc_length（Crisfield 圆柱）在 Δsoc=2.0-5.0 收敛，D 值与 basic 一致。Δsoc=1.5 处仍失败（接近损伤起始，步长 halving 后 stall）。
3. **大载荷（Δsoc≥10）**：所有方法都失败。

### basic 收敛判断 Bug 修复
- **问题**：`solve_czm_basic_step` 收敛时调用 `update_damage` 改变 `damage_states`，循环后用新 damage_states 重新组装计算 `result.residual_norm`，导致报告残差与收敛判断时不一致。
- **修复**：新增 `converged_R_norm` 变量记录收敛时刻的残差，最终报告时使用该值。修改位置 `src/CzmSolve.jl:170-247`。

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| 独立验证不回接电化学-热主链 | 避免把 CZM 验证和耦合调度混在一起 |
| 扩散应力和热应力以函数输入注入 | 满足 standalone 要求，也方便构造边界工况 |
| 验证顺序按层级推进 | 先本构，再单元，再系统，最后再做混合模式 |
| Pure Mode I 与 mix/BK 分开验证 | 这样每个基准的理论预期更清晰 |
| standalone 入口采用 provider callback | 便于复用不同热/扩散输入而不修改求解器 |
| 时间历史损伤演化用 standalone 本构驱动 | 比全场非线性迭代更稳定，也更容易直接观察 D(t) |
| 使用 `case.param.cohesive` 而非 `param_dim.cohesive` | 归一化后参数量级一致（K_n≈176 vs E_eff≈210），系统可解 |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| 使用 `param_dim.cohesive`（物理参数）导致系统病态，所有方法不收敛 | 改用 `case.param.cohesive`（归一化后参数） |
| `solve_czm_basic_step` 报告的 `residual_norm` 与收敛判断时的残差不一致 | 新增 `converged_R_norm` 记录收敛时刻残差（`src/CzmSolve.jl`） |
| arc_length 方法在损伤 regime 完全不工作 | 已修复：将球面弧长约束替换为 Crisfield 圆柱弧长法（`src/CzmSolve.jl`），收敛率从 3/8 提升到 6/8 |

## Resources
- `example/内聚力验证/czm_example.jl`
- `tools/verify_czm_unit.jl`
- `tools/verify_czm_analytical.jl`
- `tools/verify_czm_standalone.jl`
- `tools/czm_baseline_probe.jl`
- `src/czm.jl`
- `src/CzmSolve.jl`
- `src/Materialmatrix.jl`

## Visual/Browser Findings
- None.
