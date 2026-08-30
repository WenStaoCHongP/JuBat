# Findings: CZM 求解器性能瓶颈分析

## 系统配置
| 参数 | 值 |
|------|------|
| n_theta | 16 |
| cohesive 单元数 | 320 |
| bulk 单元数 | 336 |
| 节点数 nnode | 674 |
| 自由度 ndof | 1348 |
| K_bulk nnz | 16144 |
| K_coh nnz | 15376 |

## 关键发现

### F1: CZM 占 99.7% 仿真时间

| 模块 | 累计 [s] | 占比 |
|------|---------|------|
| CZM | 252.3 | 99.7% |
| SPMe | 0.274 | 0.1% |
| Thermal | 0.106 | <0.1% |
| Branch | 0.007 | <0.1% |

- 18 个时间步，每步 CZM 平均 14s
- JIT 不是瓶颈：热身运行 269s，生产运行 253s（仅 6% JIT 开销）

### F2: Newton 每次跑满 2000 迭代，全部不收敛

所有 18 步 CZM 调用全部 `converged = false`，`total_iter = 2000`（20 substeps × 100 max_iter 跑满）。

典型步骤耗时分解：

| 操作 | 耗时 [s] | 占比 |
|------|---------|------|
| **sparse_solve (apply_bc + K\R)** | **14.1** | **90%** |
| assemble_coupled | 1.5 | 10% |
| cache / strain / thermal_load | ~0 | ~0% |

### F3: 缓存复用正常

- 第 1 步：`_t_cache = 0.030s`（首次构建）
- 第 2-18 步：`_t_cache ≈ 0.3-5e-7 s`（命中缓存）
- 缓存优化工作正常，不是问题所在

---

## 根因链（四层深入）

### 第一层：`max_Δu = 1e-6` 硬编码限幅（表面原因）

`newton_raphson_czm` 中原代码（CzmSolve.jl:112-114）：

```julia
max_Δu = 1e-6
if Δu_norm > max_Δu
    Δu = Δu * (max_Δu / Δu_norm)
end
```

**诊断实验结果**（`example/czm_bottleneck_diagnosis.jl`）：

| max_Δu 设置 | 总迭代 | 收敛? | 残差行为 |
|------------|--------|-------|---------|
| **无限制** | **10** | **✓ 全部收敛** | R: 3.8e1 → 1.9e-8（二次收敛） |
| 1e-6 | 150 | ✗ 全不收敛 | R 完全停滞（3.80→3.80→3.80） |
| 1e-4 | 150 | ✗ 全不收敛 | R 几乎不动 |
| 1e-2 | 150 | ✗ 全不收敛 | R 缓慢下降但太慢 |

无限制时：每个 substep 仅 2 次迭代即收敛，载荷 0~2x 范围内均如此。
有限幅时：Newton 步长 `du_norm=0.55` 被裁剪到 `1e-6`（**55 万倍**），位移几乎不动，残差永远不降。

### 第二层：去掉限幅后产生 NaN（直接原因）

去掉 `max_Δu` 限幅后，实际仿真中 Newton 第一步就产生 **NaN**：

```
[NR-DBG] ls=1 it=1 R=NaN du=NaN
[NR-DBG] ls=1 it=2 R=NaN du=NaN
[NR-DBG] ls=1 it=3 R=NaN
```

**诊断脚本中不出现 NaN，但实际仿真中出现**。两者的区别：
- 诊断脚本：用 `fill(0.01, n_bulk)` 构造载荷，Newton 直接调用 `assemble_coupled_system`
- 实际仿真：从热模型获取真实 `dT_elem`、`Δsoc_n/p_elem`，经过 `update_czm_damage!` → `solve_czm_step` → `newton_raphson_czm`

### 第三层：NaN 跨模块传播链（传播机制）

NaN 一旦在 CZM 中产生，会通过以下路径永久污染所有后续计算：

```
CZM Newton 发散
  → NaN damage_states (D, δ_max_n, δ_max_eff)
  → update_czm_damage! 第700行: czm_mesh.damage_states = updated_czm_mesh.damage_states
  → ThermalDistributed2D_BC 第295-309行:
      读取 damage_states[i].D 和 .δ_max_n
      → compute_gap_conductance(NaN, NaN, ...) → NaN h_eff
      → K[nb,nb] -= NaN → 热刚度矩阵含 NaN
  → 热求解产生 NaN 温度场
  → NaN T_nodes_carry → NaN dT_elem → NaN F_thermo_chem
  → 下次 CZM 调用从第1次迭代就 NaN（输入已坏）
```

**关键传播节点**（代码位置）：

| 文件 | 行号 | 操作 | 风险 |
|------|------|------|------|
| CzmSolve.jl | 700 | `czm_mesh.damage_states = updated_czm_mesh.damage_states` | 无条件写入，即使是 NaN |
| Solve.jl | 279 | `u_czm_prev = u_czm_new` | 无条件更新，NaN 持久化 |
| ThermalDistributed.jl | 296-299 | 读取 `state.D` / `state.δ_max_n` | NaN 传播到热刚度矩阵 |
| Materialmatrix.jl | 323 | 读取 `damage_states[elem_idx]` | NaN 传播到材料矩阵 |

### 第四层：根因已确认并闭环

**已知事实**：
- 诊断脚本（人造载荷 `fill(0.01, n_bulk)`, F_tc norm=35）→ Newton 2次/substep 收敛 ✓
- 真实耦合路径在固定子步 / 无回退 / 强限幅时会空转或进入 NaN ✗
- 加入线搜索、自适应 load stepping 和失败回滚后，主脚本稳定通过 ✓

**结论**：
1. 不是缓存问题。
2. 不是归一化问题。
3. 不是单纯的材料本构奇异性。
4. 根因是 CZM 求解策略缺少“自适应载荷推进 + 失败回滚 + 收敛后再提交”，导致非收敛状态被回写并污染下游。

**验证结果**：
- `example/testexample.jl` 最新复跑已无 NaN。
- `CZM solve issue` 警告消失。
- `D_max = 0`，说明当前验证用例未进入损伤演化，但数值链路已经稳定。

---

## 载荷归一化验证

已确认 CZM 参数归一化正确：
- `α_eff = α_phys × T_ref = 8.94e-4` ✓
- `E_eff = E_phys / E_n = 182.7` ✓
- `T_nodes_carry` 是归一化 T* = T/T_ref ✓
- 实际仿真 dT ~ 7e-5, Δsoc ~ 0.005 → F_tc norm ~ 35
- 归一化方案本身没有问题

---

## 已实施的修复（CzmSolve.jl）

### 1. 回溯线搜索（newton_raphson_czm，替代硬编码 max_Δu）

在 Newton 迭代中：
- 计算搜索方向 `Δu = K_bc \ R_bc`
- NaN/Inf 检测：`any(isnan, Δu)` → 立即 break
- 最多 4 次回溯减半：`α = 1.0 → 0.5 → 0.25 → 0.125`
- 每次回溯重新计算 `assemble_coupled_system` 验证残差下降
- 无法找到下降方向时 break（避免空转）

**成本分析**：每次回溯需要一次 assembly (~0.75ms)，最多 3 次额外 assembly。
正常情况下（诊断脚本验证）Newton 不需要回溯，成本为零。
只在数值不稳定时才触发，比 2000 次空转节省 99% 以上。

### 2. 诊断日志（update_czm_damage!）

- 输入 NaN 检测：`@warn "CZM inputs contain NaN"`
- u_czm_prev NaN 检测：`@warn "CZM u_czm_prev contains NaN, resetting to zeros"`
- 求解结果异常检测：`@warn "CZM solve issue"`

这些日志帮助定位 NaN 首次出现的具体位置。

### 3. 自适应载荷推进

- 在 `newton_raphson_czm` 中加入 load progress 递增、失败回滚和步长减半重试
- 非收敛 substep 不再推进 `load_progress`
- 仅在收敛时才更新 `damage_states`

### 设计决策：不做静默跳过

**刻意选择不添加 NaN 静默跳过**（不隐藏错误）：
- 输入含 NaN 时：仅 `@warn`，不跳过（让错误暴露）
- 求解结果含 NaN 时：仅 `@warn`，仍然写入 damage_states（让传播链可见）
- 目的：通过日志定位根因，而不是掩盖问题
- 非收敛 substep 直接回滚，不把部分状态写回到下一时间步

---

## 2026-04-21 最终复跑

- `example/testexample.jl` 成功完成 18 步。
- wall-clock 7.887 s。
- CallModel 中 CZM model 累计 0.978 s，平均 54.324 ms/call。
- 运行中未再出现 NaN 或 `CZM solve issue` 警告。
- `D_max` / `D_mean` 仍为 0，说明当前 case 验证的是稳定性而不是损伤演化。

## 最终结论

1. 早期的 99.7% 时间瓶颈，本质上是固定子步、硬限幅和无回退造成的大量空转。
2. NaN 传播链已被切断：非收敛 damage_states 和 `u_czm_prev` 不再回写到下游。
3. 当前 case 已验证数值稳定性；若要验证损伤演化，需要更高载荷或专门的 CZM 验证 case。

## 2026-04-21 内聚力参数软化试验

- 将 Jellyroll 的 cohesive 刚度临时软化 10 倍：`K_n = 2.4e16`，`K_t = 4.3e11`。
- 复跑 `example/testexample.jl` 后，CZM 仍未完全收敛，但性能明显改善：
  - CZM 平均耗时从约 `115.6 ms/call` 降到约 `33.2 ms/call`。
  - 总 wall-clock 从约 `61.0 s` 降到约 `31.2 s`。
  - 输出中已出现明显损伤演化，`D_max ≈ 96.84%`，不再是全零。
  - 但每步仍有 `CZM solve issue`，说明仅软化刚度还不足以让 basic 方法彻底收敛。

## 2026-04-21 再次降档试验

- 继续将 Jellyroll cohesive 刚度再降 10 倍：`K_n = 2.4e15`，`K_t = 4.3e10`。
- 复跑 `example/testexample.jl` 后，`CZM solve issue` 已消失，脚本稳定完成。
- 结果表明：当前不收敛的主要诱因确实是本构刚度过硬；降到这一档后，basic 方法已经可以稳定推进。
- 代价是损伤水平显著下降到 `D_max ≈ 72.10%`、`D_mean ≈ 0.80%`，说明这是一组更“软”的试验参数，而不一定是最终物理标定值。

## 后续建议

1. 在更高载荷下再跑一次，确认损伤演化路径也能稳定收敛。
2. 若后续需要对比性能，保留当前自适应 load stepping 作为默认策略。

## 2026-04-21 求解器拆分结论

- `basic` 已拆成独立单次 Newton 路径，失败时回滚到进入求解前的位移和损伤状态。
- `arc_length` 已拆成独立预测-校正路径，不再 fallback 到 `load_substep`。
- 三种方法现在都遵循同一语义：只有收敛后才提交损伤状态，失败步不会污染下游。
- 主例程 `example/testexample.jl` 仍可稳定跑通；弧长法短探针也可完成一段完整推进。

## 资源

- 诊断脚本：`example/czm_bottleneck_diagnosis.jl`
- CZM 求解器：`src/CzmSolve.jl`
- CZM 组装：`src/czm.jl`
- 热化学载荷：`czm.jl:assemble_thermal_chemical_load`
- 热模型 BC（读取 damage）：`ThermalDistributed.jl:286-310`
