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

（各优化批次追加）
