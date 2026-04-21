# Progress: CZM 求解器性能优化

## Session Log

### 2026-04-21 优化实施与初步测试

**完成的工作：**
- 实现 CZMAssemblyCache / CZMAssemblyWorkspace 结构体 (CouplingState.jl)
- Case 结构体添加 czm_cache 字段 (SetCase.jl)
- build_czm_cache / ensure_czm_cache 函数 (czm.jl)
- assemble_czm_system 使用 mul! 和预分配稀疏矩阵 (czm.jl)
- newton_raphson_czm / solve_czm_step 接受 cache 参数 (CzmSolve.jl)
- update_czm_damage! 调用 ensure_czm_cache (CzmSolve.jl)
- JuBat.jl 导出新类型和函数

**Bug 修复：**
- CZMAssemblyWorkspace 必须在 CZMAssemblyCache 之前定义（类型依赖顺序）
- workspace 最初在 newton_raphson_czm 中每次创建，改为从 cache 复用
- A*B 矩阵乘法替换为 mul! 预分配写入

### 2026-04-21 微基准测试

**测试配置：** n_theta=16, 320 cohesive, 336 bulk, ndof=1348

```
Op                          x100 耗时    单次耗时
assemble_czm_system         0.046 s       0.46 ms
K_bulk*u                    0.001 s       0.01 ms
K_bulk+K_coh                0.008 s       0.08 ms
copy(K_total)               0.007 s       0.07 ms
apply_bc_czm                0.023 s       0.23 ms
sparse_solve                0.269 s       2.69 ms
assemble_coupled            0.059 s       0.59 ms
full newton_raphson_czm     0.024 s       24 ms (20 iters)
```

### 2026-04-21 双运行对比 (JIT 排除)

| 运行 | 总耗时 | CZM 耗时 | CZM 平均/步 |
|------|--------|---------|------------|
| Warmup | 269.0 s | 259.3 s | 14408 ms |
| Production | 252.9 s | 252.3 s | 14015 ms |
| **差异 (JIT)** | 16.1 s (6%) | 7.1 s (3%) | — |

### 2026-04-21 Newton 不收敛确认

所有 18 步 CZM 调用：`converged=false, total_iter=2000`（跑满 20×100）。

NR-DETAIL 典型数据（第5步）：
```
total_iter   = 2000
_t_thermo    = 0.00118 s
_t_assem     = 1.51 s  (2000次 × 0.75ms)
_t_solve     = 14.08 s (2000次 × 7.0ms)
```

缓存复用正常（第2步起 cache ≈ 0s）。

### 2026-04-21 诊断脚本实验

创建 `example/czm_bottleneck_diagnosis.jl`，4 组对照实验：

**实验 1: Δu 限幅**
| max_Δu | 迭代 | 收敛? |
|--------|------|-------|
| 无限制 | 10 | ✓ 全部 |
| 1e-6 | 150 | ✗ |
| 1e-4 | 150 | ✗ |
| 1e-2 | 150 | ✗ |

无限制时每个 substep 仅 2 次迭代即二次收敛。

**实验 2: 载荷子步数**（无限幅）
n_substeps 1~20 都能收敛，每个 substep 固定 2 次迭代。

**实验 3: 载荷大小**（无限幅）
scale 0.0~2.0 都能收敛。

### 2026-04-21 去掉限幅后 NaN 发现

移除 `max_Δu=1e-6` 后实际仿真 Newton 产生 NaN：
```
[NR-DBG] ls=1 it=1 R=NaN du=NaN
[NR-DBG] ls=1 it=2 R=NaN du=NaN
```

诊断脚本（人造载荷 `fill(0.01, n_bulk)`）无 NaN。
实际仿真（真实热模型 dT/SOC 载荷）产生 NaN。

**结论**：`max_Δu=1e-6` 是在掩盖 NaN 问题。真正根因是载荷量级（‖F_tc‖=195）与 Newton 求解策略不匹配。

### 2026-04-21 载荷归一化验证

确认 CZM 参数归一化正确：
- α_eff = α_phys × T_ref ✓
- E_eff = E_phys / E_n ✓
- T_nodes_carry 是归一化 T* ✓
- 归一化方案本身没有问题

### 2026-04-21 NaN 传播链确认

追踪代码发现 NaN 跨模块传播路径：
```
CZM NaN → damage_states → ThermalDistributed2D_BC → 热模型 NaN → T_nodes_carry NaN → 下次 CZM 输入 NaN
```

关键传播节点：
- CzmSolve.jl:700 — 无条件写入 NaN damage_states
- Solve.jl:279 — 无条件更新 NaN u_czm_prev
- ThermalDistributed.jl:296-299 — 读取 NaN D/δ_max_n
- Materialmatrix.jl:323 — 读取 NaN damage_states

### 2026-04-21 回溯线搜索实施

在 `newton_raphson_czm` 中实现回溯线搜索：
- 替代硬编码 `max_Δu=1e-6`
- NaN/Inf 检测 → 立即 break
- 最多 4 次减半回溯（α = 1.0 → 0.125）
- 每次回溯重新评估 `assemble_coupled_system` 验证残差下降
- 无法找到下降方向时 break

添加诊断日志：
- `@warn "CZM inputs contain NaN"` — 输入异常检测
- `@warn "CZM u_czm_prev contains NaN"` — 初始状态异常
- `@warn "CZM solve issue"` — 求解结果异常

**设计决策**：不做静默跳过，让错误暴露，通过日志定位根因。

### 2026-04-21 复跑验证与诊断增强

- 修复 `CzmSolve.jl` 中 `@warn` 多行语法错误，恢复模块加载。
- 运行 `example/czm_bottleneck_diagnosis.jl`：人工载荷下纯 CZM Newton 仍可收敛，未复现 NaN。
- 运行 `example/testexample.jl`：完整耦合流程可正常完成，只有 `CZM solve issue` 警告，没有 NaN 相关告警。
- 增强 `update_czm_damage!` 的输入检查，新增 `T_nodes_carry` / `thermal2D element soc_n` / `thermal2D element soc_p` 的 NaN 分类日志。
- 现阶段结论：当前主例程不再复现 NaN，剩余问题是 CZM 不收敛而非 NaN 传播。

### 2026-04-21 主验证脚本复跑

- 重新运行 `example/testexample.jl`，仿真完整完成，没有发散中断，也没有 NaN。
- 仍有 `CZM solve issue` 警告，但残差维持在有限值范围内，未出现 `Inf/NaN`。
- 目前状态：发散问题已压住，后续重点是继续降低 CZM 残差并恢复真正收敛。

### 2026-04-21 最终稳定复跑与文档同步

- 将 `newton_raphson_czm` 改为自适应载荷推进：失败回滚、步长减半重试、收敛后再提交。
- 恢复并整理 `solve_czm_step` 的分支结构，保留 `load_substep` 作为稳定默认路径。
- 重新运行 `example/testexample.jl`：18 步完整完成，未再出现 NaN 或 `CZM solve issue` 警告。
- 最终性能：CZM model 累计 0.978 s，平均 54.324 ms/call，wall-clock 7.887 s。
- 文档已同步到 `task_plan.md`、`findings.md` 和当前进度记录。

### 2026-04-21 弧长法独立实现与验证

- 抽出 `solve_czm_basic_step` 和 `solve_czm_arc_length_step`，把 `solve_czm_step` 收敛为纯调度入口。
- `basic` 与 `arc_length` 都改为失败回滚，且只在收敛后提交损伤状态。
- `arc_length` 不再 fallback 到 `load_substep`，而是独立执行增广预测-校正流程。
- 复跑 `example/testexample.jl` 仍然成功。
- 运行短探针脚本验证弧长法路径：可稳定完成 68 步推进并输出结果。

### Bug 修复

| Bug | 文件 | 状态 |
|-----|------|------|
| `Δu = K_bc \ R_bc` 重复计算（第107-109行） | CzmSolve.jl | 已修复 |
| `max_Δu=1e-6` 硬编码限幅导致不收敛 | CzmSolve.jl | 已移除，替换为线搜索 |
| profiling 计时代码残留 | CzmSolve.jl | 已清理 |
| debug println 残留 | CzmSolve.jl | 已移除 |

## 当前状态

**已确认**：
1. 缓存优化有效（K_bulk/geom/BC/workspace 复用正常）
2. 根因已确认：固定子步 + 硬限幅 + 无回退导致空转/NaN
3. 去掉限幅后暴露的 NaN 已被线搜索和自适应子步拦住
4. `Δu = K_bc \ R_bc` 重复计算已修复
5. NaN 通过 damage_states → 热模型 → CZM 的跨模块传播链已切断
6. 载荷归一化验证通过，不是归一化的问题
7. 回溯线搜索 + 自适应 load stepping 已通过主脚本验证

**当前结论**：
- 主验证脚本已稳定通过，没有 NaN 或收敛警告
- CZM 的时间占比已回落到合理数量级
- 仅剩可选验证：更高载荷或专门损伤 case 下的 CZM 演化

## 测试结果

| 测试 | 结果 |
|------|------|
| testexample.jl（原始） | 18步完成，但 CZM 全部不收敛 |
| 微基准 newton_raphson_czm | 0.024s (20 iters) |
| 诊断脚本（无限幅） | 10 iters / 5 substeps，全部收敛 |
| 诊断脚本（有限幅 1e-6） | 150 iters，全部不收敛 |
| 实际仿真（去掉限幅） | NaN 传播，仍不收敛 |
| 实际仿真（线搜索 + 自适应子步） | 18 步完成，无 NaN / 无警告 |

## 2026-04-21 内聚力参数试验

- 试验性降低 Jellyroll cohesive 刚度 10 倍：`K_n = 2.4e16`，`K_t = 4.3e11`。
- 重新运行 `example/testexample.jl` 后：
	- 仿真仍然完成，未崩溃。
	- `CZM solve issue` 仍存在，但残差增长明显放缓，CZM 平均耗时显著下降。
	- `D_max` 从全零变为约 `96.84%`，说明损伤演化已经真正被触发。
	- 结论：问题没有完全消失，但刚度软化确实缓解了收敛与传播问题，说明本构过硬是重要因素之一。
