# 仿真提速-SPMe热CZM findings

## §锚点-1（2026-08-19，HEAD=8ad8953）

- 命令：`GKSwstype=100 JULIA_NUM_THREADS=1 julia-1.11.2 --startup-file=no example/testexample.jl`（本机路径 `C:\Users\19303\AppData\Local\Programs\Julia-1.11.2`）
- 退出码：0
- 墙钟：**249.4 s**（real 4m9.443s）
- 网格/步数：ne=1682、nT=1763、CallModel calls=19（与档案 metrics.toml 一致）
- 科学结果（与档案一致，按打印精度）：初始电压 4.0367 V、最终电压 3.9438 V、容量 0.0833 Ah、温度 298.15~299.00 K、CZM D_max/D_mean=0.0000%、19 次 CZM 更新全收敛
- PNG SHA-256（output/testexample_results.png）：`e879c50bad62222c0f7b39da98f480db74e35c31320da9bdad24b3129251dfb9`
  （与跨机档案 SHA `4ba6207c...` 不同，符合"档案跨机不可复现"预期，本机锚点为准）

### timing 分解（CallModel 口径）

| 模块 | total [s] | ratio | avg [ms/call] |
|------|-----------|-------|---------------|
| SPMe 求解 | 5.586 | 2.66% | 293.975 |
| 分流求解器 | 1.626 | 0.77% | 85.604 |
| 热分布式模型 | 1.614 | 0.77% | 84.957 |
| **CZM 模型** | **201.382** | **95.80%** | **10599.072** |

**结论：CZM 单点瓶颈（95.8%），单次更新 ~10.6s × 19 次。** 其余三项合计 <4.3%。墙钟 249.4s 与 CallModel 计时合计 210.2s 的差 ~39s 为 Julia 启动编译、初始化、逐步求解与绘图。

## §锚点-2（2026-08-19）

- 退出码：0，墙钟 245.1s（real 4m5.109s）
- PNG SHA-256：`e879c50bad62222c0f7b39da98f480db74e35c31320da9bdad24b3129251dfb9`——**与锚点-1 完全一致**

## §判据裁决

PNG SHA 本机稳定：**是**（两次运行逐位一致）；后续批次采用**四判据**（退出码、网格/步数、metrics 打印精度、本机 PNG SHA，锚点 SHA=`e879c50b...`）。

## §实测占比

见 §锚点-1 表。占比采用基线 debug 固有口径（testexample.jl:56 已设 debug_coupling=true，timing 累计本就在 result 键中，开销为基线固有），Profile 未做（spec §2 允许的简化，理由：占比差距悬殊——95.8% vs <3%——交叉确认不改变任何门的结论）。

## §事实核查

### 3a/3b（probe_matrix_constancy.jl，2026-08-19）

- **M identical = true**（全局质量矩阵跨状态常量）
- **K identical = false，差异 127832 个非零全部在化学块，热块（含 BC）差异 = 0**：
  - 热矩阵 KT/MT/KT_bc 跨状态常量（`jacobi="update"` 只作用于化学块 K_chem，电解液扩散矩阵随 ce 状态更新）
- **CZM 单次更新：iterations=1、converged=true**，单次耗时 10.6s——瓶颈不在 Newton 迭代次数，而在每迭代的装配+LU 求解+线搜索组装（czm_load_steps=10 载荷子步 × 每子步装配/求解/线搜索）
- **等价性 1（lu(A)\b == A\b）：true**
- **等价性 2（手写同序加法 == SparseArrays +）：true**

### CZM 向量化计划（2026-04-20）残留核对

- 该计划 27 个步骤**全部未勾选**（planning-with-files/向量化CZM/progress.md 显示该轮止步于"产出 spec+plan"，未系统执行）
- 但代码对照确认以下项已通过后续工作落地：K_bulk 缓存（czm.jl:440/:561 `CZMAssemblyCache`）、K_coh pattern 复用+fill! 清零（czm.jl:101-131）、geom_cache（czm.jl:152-159）、mul! 无分配（czm.jl:195-225）
- **残留热点**（Task 7 范围，互补不冲突）：
  1. `K_total = K_bulk + K_coh` 每 Newton 迭代分配（czm.jl:572）
  2. `f_int_bulk = K_bulk*u`、`f_int_total = f_int_bulk + f_int_coh` 每次分配（czm.jl:569/:575）
  3. `K_coh[dofs[a],dofs[b]] +=` 标量稀疏索引累加 801×64 次/迭代（czm.jl:251）
  4. 线搜索路径全量组装 K_total（CzmSolve.jl:113，只用 f_int）
- 增补机制候选（不在本计划内）：`compute_czm_strain_inputs` 批量化、`clone_czm_mesh_with_damage` 每步克隆、`apply_bc_czm` 每迭代 copy

## §批次裁决（2026-08-19，依计划执行门速查）

| 任务 | 门条件 | 实测 | 裁决 |
|------|--------|------|------|
| Task 3 热矩阵缓存 | 热矩阵块跨步常量 | 常量成立，**但热模块占比仅 0.77%（1.6s/249s）** | **关**（投入产出不成立，记录不动） |
| Task 4 全局 M/K 缓存 | 3a 全局 true | K 化学块随状态变 | **关** |
| Task 5 单元循环分配消除 | SPMe ratio ≥10% | 2.66% | **关** |
| Task 6 按 dt 分解缓存 | 3a true 且等价性 1 true | 3a false（K 每步变） | **关** |
| **Task 7 CZM 装配优化** | CZM ratio ≥10% + 残留核对 + 等价性 2 | **95.80% + 完成 + true** | **开（唯一开门批次）** |
| Task 8 收尾 | 恒开 | — | 开 |

**执行顺序：Task 7 → Task 8。**

**预期目标**：CZM 201.4s（四项残留：装配分配、标量稀疏索引、线搜索冗余组装）。保守折减按计划取开门模块 timing × 0.5 ≈ **预期总墙钟 249s → ~150s**（实际以四判据 + timing 实测为准；若不足，Task 8 走增补机制）。

## §批次记录

### Task 7（CZM 装配优化，commit 2f7d129，2026-08-19；**src 改动已于同日 revert（commit 5812c79，用户决定）**）

- 改动：`CZMAssemblyWorkspace` +6 缓冲字段；`assemble_czm_system` nzval 直索引（`_nz_index` + `cohesive_nzidx` 首建）；`assemble_coupled_system` f_int/K_total 预装配 + `assemble_K` 开关；线搜索 `assemble_K=false`
- 四判据：**全部通过**（退出码 0、19 步、科学结果一致、PNG SHA `e879c50b...` 与锚点逐位一致）
- timing：CZM 201.382s → 197.982s（**-1.7%**），总墙钟 249.4s → 244.0s——**收益不及预期**（预期 ×0.5）

### 根因再定位（profile 证据链）

1. `probe_czm_profile.jl`：预热后单次 `update_czm_damage!` 仅 **0.24s**（缓存命中、t=0 零载荷 → 1 次 Newton 迭代收敛）；czm ndof=43758、n_coh=6728
2. `probe_gc_attribution.jl`（20s 短仿真）：CZM avg 10.08s/次 × 15 次真实存在；**GC 仅 1.25%**（排除 GC 归因）；全程序 123.65M 分配 / 20.97 GiB
3. 缓存失效排除：`update_czm_damage!` 不替换 `case.czm_mesh`（CouplingState.jl:637 只写回 damage_states），`param_cache.id` 为稳定内容哈希，无人置 `cache.valid=false`
4. `probe_solve_profile.jl`（@profile 全仿真采样）：**84% 采样（25723/30538）落在 `src/CzmSolve.jl:222` 的 `K_bc \ R_bc`**——Newton 迭代内（:195 for 循环）每次迭代对 43,758 DOF 稀疏矩阵做完整 UMFPACK 分解
5. 机制解释：t=0 零载荷 1 次迭代即收敛（探针快）；仿真中 dT/Δsoc 载荷非零 → Newton 多轮迭代 × 每轮 LU ≈ 10s/次。本场景 D=0、δ≈0，bilinear tangent 固定 → **K_total/K_bc 数值跨迭代跨步几乎不变，却在每次迭代被重新分解**

### 增补批次建议（回 spec 补充后另立，不在本计划实施）

**K_bc 分解因子按内容判据复用**：
- 位置：`solve_czm_basic_step`（CzmSolve.jl:219-225）及 arc_length/load_substep 同型调用点
- 改法：缓存上一次分解的 `(K_bc_nzval_copy, factorization)`；每次 Newton 迭代装配出 K_total 后，若 `nonzeros(K_bc_new) == K_bc_nzval_copy`（O(nnz) 比较，远便宜于 LU）则 `ldiv!(F, R_bc)` 复用因子，否则重新分解并更新缓存
- bit 一致论证：同矩阵同 rhs 下，`K_bc \ R_bc`（每次新分解，确定性）与复用因子的 `F \ R_bc` 解逐位一致（等价性 1 `lu(A)\b == A\b` 已验证；数值不等时走原路径不变）
- 预期收益：本场景每次更新 ~10 迭代仅 1 次真分解 + 9 次回代 → CZM 198s → 约 20-40s，总墙钟 244s → **约 90-130s**（达到 spec 预期 ×0.5 目标）
- 附带候选（次要）：`apply_bc_czm` 每迭代 copy K_total（43758² copy）可改为修改-恢复式或复用判据的一部分

### 回滚记录（2026-08-19，用户决定）

- commit 2f7d129（src/czm.jl、src/CzmSolve.jl、src/CouplingState.jl 三文件）已 revert（5812c79），src 完全还原到批次 0 锚点状态
- **保留**：批次 0 锚点与占比数据、矩阵常量性/等价性结论、profile 根因定位（84% 在 CzmSolve.jl:222 每迭代 LU 分解）、增补批次建议——这些是纯文档与探针脚本，不含 src 改动
- 增补批次（K_bc 因子内容判据复用）若获批，将在还原后的干净基线上实施

### Task 8（收尾，2026-08-19）

- 汇总：批次 0（锚点+裁决）→ Task 7（-1.7%，四判据通过）→ 根因再定位（84% 在每迭代 LU 分解）→ 增补建议如上
- 关门任务（Task 3/4/5/6）未实施，理由见 §批次裁决
- 目标达成情况：本计划内未达 ×0.5 预期（-2%）；瓶颈根因已定位并有 bit 一致可行的增补方案，待用户批准后另立批次
