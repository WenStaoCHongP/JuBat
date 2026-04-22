# Task Plan: 内聚力模块拆分与简化

## Goal
把内聚力模块从当前的"大而全"实现拆成边界清晰、职责单一、便于后续维护的若干小模块，同时保持现有公开接口和数值行为尽量不变。

## Current Phase
Phase 3+4 已完成，Phase 5 (验证与文档) 待执行

## Phases

### Phase 1: 责任边界盘点
- [x] 确认当前核心入口主要集中在 `src/CzmSolve.jl` 和 `src/czm.jl`
- [x] 识别出网格、状态、组装、求解、后处理的混杂点
- [x] 记录已知归一化与后处理约束
- [x] 量化基线：CzmSolve.jl 1001 行/19 函数，czm.jl 659 行/11 函数
- [x] 重复代码统计：~314 行重复（BC 提取 60 行、clone 44 行、Newton 循环 180 行、末尾重组 30 行）
- **Status:** complete

### Phase 2: 模块拆分设计
- [x] 设计目标模块边界与文件划分（见下方 "Target File Layout"）
- [x] 明确哪些符号继续保留为兼容导出（见下方 "Protected Symbols"）
- [x] 确定 `JuBat.jl` 的 include 顺序和依赖顺序
- [x] 回答 Key Questions（见下方 "Key Questions - Resolved"）
- **Status:** complete

### Phase 3: 代码结构简化方案
- [x] 3.1 提取 `clone_czm_mesh` 公共 helper（消除 4 处重复，~44 行） — `bd9d8f2`
- [x] 3.2 提取 `extract_bc_dofs` 公共 helper（消除 3 处 BC 提取重复，~60 行） — `9499832`
- [x] 3.3 提取 Newton 迭代 + 线搜索公共 helper（`backtrack_line_search!`，仅服务 basic 求解器） — `56ce8e1`
- [x] 3.4 将后处理/统计函数从 CzmSolve.jl 中拆分到 `src/CzmPostProcess.jl` (117 行) — `5272682`
- [x] 3.5 将参数计算函数 (`compute_czm_effective_params`, `compute_czm_strain_inputs`, `update_czm_damage!`) 并入 `src/CouplingState.jl` — `7bd29ee`
- **Status:** complete
- **说明:** 3.3 的 `backtrack_line_search!` 限定只服务 `solve_czm_basic_step`，因为 `newton_raphson_czm` 使用惩罚式 BC 残差语义，不能复用零化式 helper。arc_length 和 newton_raphson 仍有部分 Newton 循环重复，但风险可控。

### Phase 4: 迁移与兼容
- [x] 先做无行为变化的内部重构（Phase 3 的 helper 提取）
- [x] 再做文件级拆分与入口收束
- [x] 保持旧 API 可用，必要时加薄兼容层
- **Status:** complete
- **结果:** CzmSolve.jl 从 1001 行瘦身至 647 行；CzmPostProcess.jl 117 行；CouplingState.jl 扩展 375 行。所有公开接口保持不变。

### Phase 5: 验证与文档
- [ ] 跑最小 CZM 探针和主耦合例程
- [ ] 检查是否引入 include 顺序或类型初始化错误
- [ ] 同步更新文档中的模块说明
- **Status:** pending

## Key Questions - Resolved

### Q1: 拆分后是否保留 `src/CzmSolve.jl` 作为调度入口，还是直接拆成多个并列文件？

**决定：保留 `CzmSolve.jl` 作为调度入口，同时拆出 `CzmPostProcess.jl`，并在现有 `CouplingState.jl` 中扩展耦合 helpers。**

理由：
- `CzmSolve.jl` 已在 `JuBat.jl` 中被 include，且有大量 export 绑定到它
- 外部调用点（`tools/`, `example/`）直接 import `solve_czm_step`, `newton_raphson_czm` 等符号
- 保留调度入口 + 拆出辅助文件是最小改动策略

目标文件布局：
```
src/CzmSolve.jl          → 调度层 + 求解器 (保留，瘦身到 ~600 行)
src/CzmPostProcess.jl    → 后处理/统计/损伤管理 (新文件, ~120 行)
src/CouplingState.jl     → 耦合状态 + 参数计算/应变输入/损伤更新 (现有文件, 新增 ~90 行)
```

### Q2: 哪些函数和类型必须保持原名，以避免影响现有例程和导出接口？

**决定：以下符号必须保持原名和签名。**

已导出符号（`JuBat.jl` 中显式 export）：
- `newton_raphson_czm`, `solve_czm_step`
- `get_damage_statistics`, `check_fracture_criterion`, `reset_damage_states`, `accumulate_cycle_damage`
- `czm_output_to_variables`
- `bilinear_traction`, `bilinear_tangent`, `update_damage` (Materialmatrix.jl)

被外部脚本直接调用的未导出符号：
- `compute_czm_effective_params` (5 个脚本)
- `compute_czm_strain_inputs` (4 个脚本)

被内部调用的入口：
- `update_czm_damage!` (Solve.jl:276)

**策略：** 后处理拆到 `CzmPostProcess.jl`，耦合函数并入现有 `CouplingState.jl`，因此 `CzmSolve.jl` 只需要保留 `include("CzmPostProcess.jl")`；`CouplingState.jl` 继续沿用现有 include 链，不额外新增文件。

### Q3: 先拆"纯逻辑"还是先拆"高频调用路径"，哪个对风险最小？

**决定：先拆"纯逻辑"（后处理/统计），再抽 helper，最后动求解器。**

风险排序（从低到高）：
1. **最低风险：后处理/统计函数** — `get_damage_statistics`, `check_fracture_criterion`, `reset_damage_states`, `accumulate_cycle_damage`, `czm_output_to_variables` — 无状态、纯计算、独立性强，可以直接搬到新文件
2. **低风险：参数计算函数** — `compute_czm_effective_params`, `compute_czm_strain_inputs` — 只被 `update_czm_damage!` 调用，搬迁不影响求解逻辑
3. **中等风险：公共 helper 提取** — `clone_czm_mesh`, `extract_bc_dofs`, Newton 循环 — 涉及求解器内部重构，但只改内部结构不改外部接口
4. **最高风险：求解器本体** — 三个求解器的结构变化 — Phase 3 中仅在 helper 提取后进行

## Target File Layout

### 实际完成后的文件结构

```
src/
├── CzmSolve.jl            (调度 + 求解器，647 行)
│   ├── CZMResult
│   ├── clone_damage_states
│   ├── clone_czm_mesh_with_damage
│   ├── apply_czm_dirichlet!
│   ├── zero_czm_bc_entries!
│   ├── fill_czm_result!
│   ├── build_arc_length_augmented_matrix
│   ├── backtrack_line_search!    ← 新增 helper (原 3.3)
│   ├── extract_bc_dofs           ← 新增 helper (原 3.2)
│   ├── solve_czm_basic_step
│   ├── solve_czm_arc_length_step
│   ├── newton_raphson_czm
│   └── solve_czm_step
│
├── CzmPostProcess.jl      (新文件，117 行)
│   ├── get_damage_statistics
│   ├── check_fracture_criterion
│   ├── reset_damage_states
│   ├── accumulate_cycle_damage
│   └── czm_output_to_variables
│
├── CouplingState.jl       (现有文件扩展，375 行)
│   ├── ... (原有类型和 helpers)
│   ├── compute_czm_effective_params  ← 从 CzmSolve.jl 迁入
│   ├── compute_czm_strain_inputs     ← 从 CzmSolve.jl 迁入
│   └── update_czm_damage!            ← 从 CzmSolve.jl 迁入
│
├── czm.jl                 (组装层，659 行不变)
│
└── JuBat.jl               (include 顺序已更新)
  ├── include("CouplingState.jl")
  ├── include("czm.jl")
  ├── include("CzmSolve.jl")
  ├── include("CzmPostProcess.jl")
  └── ...
```

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 先做职责切分，再做物理文件搬迁 | 先稳定边界，再减少后续反复改动 |
| 保留兼容入口层 | 避免一次性改动过大导致现有脚本失效 |
| 先抽公共 helper，再减少主函数体积 | 最低风险，能立刻降低重复和嵌套深度 |
| 保留 `CzmSolve.jl` 作为调度入口 | 外部调用点和 export 绑定最少的改动路径 |
| 先拆后处理/统计（最低风险） | 无状态、纯计算、独立性强 |
| 新文件通过 include 引入而非 module | Julia 的 module 级 include 保证符号在同一作用域，无需重命名 |

## Protected Constraints (不可破坏)
以下代码路径/语义在重构过程中必须保持不变：
1. **线搜索逻辑**：`for _ in 1:8; α *= 0.5; end` — 回溯线搜索，最多 8 次减半
2. **失败回滚**：`if !converged; u = u_start; damage_states = damage_start; end` — 未收敛时完整回滚
3. **损伤提交门控**：`if result.converged; czm_mesh.damage_states = updated_czm_mesh.damage_states; end` — 只在收敛时提交损伤
4. **载荷子步自适应**：`step_size *= 0.5` 和 `step_size *= 1.25` — 子步失败时缩小、成功时放大
5. **BC 处理**：惩罚法 (penalty=1e12) + Dirichlet 直接赋值两种方式并存

## Rollback Strategy
- 所有拆分工作在 `czm-refactor` feature branch 上进行
- 每完成一个 Phase 的子任务，在 feature branch 上提交，确保 main 始终可用
- 回归检测手段：
  - 运行 `example/testexample.jl` 或 `example/coupled_czm_thermal_example.jl` 检查端到端输出
  - 对比拆分前后的 `D_max`, `D_mean`, `cell voltage [V]`, `thermal2D temperature [K]` 数值
  - 允许的数值偏差：相对误差 < 1e-8（由求解器 tol 决定）
- 如果某个 Phase 引入回归，直接 `git revert` 该次提交

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| None yet | 1 | - |

## Notes
- 这次拆分的优先级不是改物理模型，而是把职责边界整理干净。
- 任何结构改动都要优先保住现有线搜索、失败回滚和损伤状态提交语义。
