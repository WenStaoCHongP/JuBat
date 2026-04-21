# CZM 求解器向量化与缓存优化实施计划

> **Spec Document:** `docs/superpowers/specs/2026-04-20-czm-vectorized-solver-design.md`
>
> **For agentic workers:** REQUIRED: Use superpowers:executing-plans or superpowers:subagent-driven-development (if available) to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **Context:** 基于 `docs/代码优化/15_CZM耦合修复与单步损伤演化.md` 的性能结论，CZM 现在是 CZM-active 场景中的主耗时项，`solve_czm_step` / `update_czm_damage!` 单次更新约 12s。
>
> **Updated:** 2026-04-20 — 基于代码审查修正：Task 7（冗余更新）已确认完成，K_bulk 缓存提升为独立高优先级任务。

**Goal:** 通过"缓存 + 批处理单元核"优化，将单次 CZM 更新时间显著压缩，同时保持损伤演化、位移、牵引力和终止行为与当前实现一致。

**Architecture:** 三个递进阶段：
1. 定义缓存契约 + 基线确认。
2. 缓存 K_bulk（最高 ROI）+ 几何/DOF 映射预计算 + 输入路径固化。
3. 重写 cohesive 单元核为无分配工作区模式。
4. 回归验证。

**Tech Stack:** Julia 1.x, sparse assembly, `@views`, `@inbounds`, `mul!`, `Threads.@threads`（可选后置）

**Baseline (representative CZM-active case):**
| Metric | Current | Target |
|--------|---------|--------|
| `update_czm_damage!` 单次耗时 | ~12s | < 8s (phase 1) / < 6s (phase 2) |
| CZM 在总耗时中的占比 | ~98.9% | 明显下降 |
| `assemble_czm_system` 临时分配 | ~80K 次/更新 | 显著降低 |
| 阶段末重复 CZM 更新 | **已修复** | N/A |

**性能瓶颈量化：**
- 801 cohesive 单元 × 2 Gauss 点 × ~50 Newton 总迭代 ≈ ~80K 次临时数组分配/更新
- `assemble_bulk_stiffness` 每次 Newton 迭代重建（20 load_steps × ~3 Newton = ~60 次冗余 bulk 重组）
- triplet 向量 `push!` 每轮重建：~801 × 64 entries × ~60 轮

---

## 设计约束

1. **不改物理模型**：不更改损伤准则、牵引-分离律或热-化学耦合关系。
2. **外层 Newton 保持串行**：只优化单元核和缓存，不把依赖迭代强行广播化。
3. **缓存必须可失效**：mesh 或材料参数变化时必须重建，不能把旧矩阵错误复用到新场景。
4. **先单线程稳定，再考虑并行**：线程优化只能作为后置步骤，不能掩盖单线程内存布局问题。
5. ~~去重更新路径~~：**已由 `15_CZM耦合修复` 完成**。CycleSolver.jl:157-160 已含修复注释。

---

## File Structure

| File | Change Type | Responsibility |
|------|------------|----------------|
| `src/CouplingState.jl` | Modify | 新增 `CZMAssemblyCache` / `CZMAssemblyWorkspace` 结构体 |
| `src/SetCase.jl` | Modify | `Case` 新增 `czm_cache::Union{Nothing, CZMAssemblyCache}` 字段 |
| `src/czm.jl` | Modify | K_bulk 缓存、element kernel 重写、triplet 预分配、装配逻辑批处理化 |
| `src/CzmSolve.jl` | Modify | 输入准备批量化、`update_czm_damage!` 接缓存/工作区 |
| `src/Solve.jl` | Modify | 传递缓存、保留计时、控制 CZM 更新节奏 |
| `src/Option.jl` | Read | `czm_update_interval`、`czm_load_steps` 参数定义（可能需调默认值） |
| `src/JuBat.jl` | Modify | 新增结构体/接口的 include/export 对齐 |

**不改动的文件（显式排除）：**
- `src/CycleSolver.jl` — 冗余更新已在上一轮修复，无需再改
- `src/czm.jl:assemble_thermal_chemical_load` — 仅在 `newton_raphson_czm` 入口调用一次且 load_factor 缩放已正确，无需优化

---

## Chunk 1: Cache Contract & Baseline

### Task 1: 定义 CZM 缓存与工作区结构体

**Files:**
- Modify: `src/CouplingState.jl`
- Modify: `src/SetCase.jl`
- Modify: `src/JuBat.jl`

- [ ] **Step 1: 在 `CouplingState.jl` 新增 `CZMAssemblyCache` 结构体**
  - `K_bulk::SparseMatrixCSC{Float64, Int64}` — bulk 刚度矩阵（E/ν 不变时复用）
  - `bulk_dofs::Vector{Vector{Int64}}` — 每个 bulk 元素的 DOF 映射
  - `cohesive_geom::Vector{CohesiveElementGeom}` — 预计算的 cohesive 单元几何
  - `bc_dofs::Vector{Int64}` — 边界条件 DOF
  - `bc_vals::Vector{Float64}` — 边界条件值
  - `E_eff::Float64`, `ν_eff::Float64` — 用于缓存失效判断
  - `valid::Bool` — 缓存有效性标志

- [ ] **Step 2: 新增 `CohesiveElementGeom` 结构体**（内嵌或独立）
  - `length::Float64` — 单元长度
  - `n_vec::SVector{2, Float64}` — 法向量
  - `t_vec::SVector{2, Float64}` — 切向量
  - `R::SMatrix{2,2, Float64}` — 旋转矩阵
  - `dofs::SVector{8, Int64}` — 全局 DOF 编号
  - `nodes::SVector{4, Int64}` — 节点编号 (n1,n2,n3,n4)
  - `gauss_wts::Vector{Float64}` — Gauss 权重
  - `gauss_pts::Vector{Float64}` — Gauss 点

- [ ] **Step 3: 在 `CouplingState.jl` 新增 `CZMAssemblyWorkspace` 结构体**
  - `u_e::MVector{8, Float64}` — 单元位移
  - `B_global::MMatrix{2, 8, Float64}` — 全局 B 矩阵
  - `K_e::MMatrix{8, 8, Float64}` — 单元刚度
  - `f_int_e::MVector{8, Float64}` — 单元内力
  - `I_idx::Vector{Int64}` — triplet 行索引（预分配容量，`empty!` 复用）
  - `J_idx::Vector{Int64}` — triplet 列索引
  - `K_vals::Vector{Float64}` — triplet 值
  - `f_int_coh::Vector{Float64}` — 全局内力向量

- [ ] **Step 4: 修改 `Case` 结构体（`SetCase.jl`）**
  - 新增字段 `czm_cache::Union{Nothing, CZMAssemblyCache}`
  - 更新 5 参数构造器：传入 `nothing`

- [ ] **Step 5: 更新 `JuBat.jl`**
  - 确保 `CZMAssemblyCache`、`CZMAssemblyWorkspace`、`CohesiveElementGeom` 的 export

---

### Task 2: 基线确认

验证当前代码行为，确认优化前的性能基准。

**Files:**
- Read/inspect: `src/CzmSolve.jl`, `src/czm.jl`, `src/Solve.jl`

- [ ] **Step 1: 验证无冗余更新**
  - 确认 CycleSolver.jl:157-160 只含注释，无实际调用
  - 确认 Solve.jl:270-285 的 `czm_update_interval` 逻辑正确

- [ ] **Step 2: 记录当前基线数值**
  - `assemble_bulk_stiffness` 每次调用耗时和分配数
  - `assemble_czm_system` 每次调用耗时和分配数
  - `newton_raphson_czm` 中总 Newton 迭代数和 load_steps 数
  - `assemble_thermal_chemical_load` 确认只在入口调用一次

- [ ] **Step 3: 定义缓存失效条件**
  - `E_eff` 或 `ν_eff` 变化 → `K_bulk` 失效
  - mesh 变更 → 所有缓存失效
  - cohesive 参数变更 → 仅影响 cohesive 几何缓存（不影响 K_bulk）

---

## Chunk 2: Static Precompute — K_bulk 缓存（最高 ROI）

### Task 3: 缓存 K_bulk 与 bulk DOF 映射

**这是整个优化中收益最高的单点。** 每次 CZM 更新中 `assemble_bulk_stiffness` 被调用 ~60 次（20 load_steps × ~3 Newton），而 E/ν 不变时 K_bulk 完全相同。

**Files:**
- Modify: `src/czm.jl`
- Modify: `src/CzmSolve.jl`

- [ ] **Step 1: 新增 `build_czm_cache` 函数**
  - 调用 `assemble_bulk_stiffness` 一次，存储结果
  - 预计算并存储每个 bulk 元素的 DOF 映射
  - 存储当前 `E_eff`、`ν_eff` 用于失效判断
  - 在 `update_czm_damage!` 入口处构建/复用缓存

- [ ] **Step 2: 修改 `newton_raphson_czm` 接受缓存**
  - 新增关键字参数 `cache::Union{Nothing, CZMAssemblyCache}=nothing`
  - 将 `assemble_coupled_system` 调用改为接受 `K_bulk` 缓存
  - 如果 cache 为 nothing，回退到每轮重算（兼容旧行为）

- [ ] **Step 3: 修改 `assemble_coupled_system` 接受 `K_bulk`**
  - 新增关键字参数 `K_bulk_cached::Union{Nothing, SparseMatrixCSC}=nothing`
  - 如果提供了缓存，跳过 `assemble_bulk_stiffness` 调用
  - 保持 `K_bulk * u` 计算 `f_int_bulk` 不变

- [ ] **Step 4: 缓存失效逻辑**
  - 在 `update_czm_damage!` 入口检查 `case.czm_cache`
  - 如果 `nothing` 或 `E_eff/ν_eff` 不匹配，重建缓存
  - 设置 `cache.valid = true`

---

### Task 4: 预计算 cohesive 单元几何 + 固化输入路径

将静态几何信息和输入准备合并处理。

**Files:**
- Modify: `src/czm.jl`
- Modify: `src/CzmSolve.jl`

- [ ] **Step 1: 预计算 `CohesiveElementGeom` 数组**
  - 遍历 `czm_mesh.cohesive_elements` 一次
  - 为每个元素计算并存储 `length`、`n_vec`、`t_vec`、`R`、`dofs`、Gauss 常量
  - 存入 `CZMAssemblyCache.cohesive_geom`

- [ ] **Step 2: 预计算边界条件缓存**
  - 将 `identify_bc_nodes_czm` 的结果缓存
  - 存入 `CZMAssemblyCache.bc_dofs` / `bc_vals`

- [ ] **Step 3: 批量化 `compute_czm_strain_inputs`**
  - 保持输出为单元级数组
  - `epsilon_0_elem` 在一次更新中只算一次（当前已是如此，确认无重复）

- [ ] **Step 4: 将几何缓存传入 `assemble_czm_system`**
  - 新增关键字参数 `geom_cache::Union{Nothing, Vector{CohesiveElementGeom}}=nothing`
  - 如果提供了缓存，跳过单元内几何重算

---

## Chunk 3: Batched Element Kernel

### Task 5: 重写 cohesive 单元核为无分配工作区模式

这是主要性能改造点。将 `assemble_czm_system` 从"每单元分配"改为"预分配工作区 + 无分配内核"。

**Files:**
- Modify: `src/czm.jl`

- [ ] **Step 1: 创建 `CZMAssemblyWorkspace` 实例化函数**
  - `create_czm_workspace(ndof, n_coh, avg_nnz_per_elem)` — 预分配所有工作数组
  - triplet 缓冲区容量：`n_coh * 64`（每个 cohesive 单元最多 8×8=64 个非零元）

- [ ] **Step 2: 消除单元内临时分配**
  - `u_e`：改为 `ws.u_e`，从全局 `u` 向量 `copyto!` 而非 `zeros`
  - `B_global`：改为 `ws.B_global`，每 Gauss 点 `fill!(ws.B_global, 0.0)` 后赋值
  - `K_e`：改为 `ws.K_e`，每单元开始时 `fill!(ws.K_e, 0.0)`
  - `f_int_e`：改为 `ws.f_int_e`，每单元开始时 `fill!(ws.f_int_e, 0.0)`

- [ ] **Step 3: 改造 triplet 装配路径**
  - 预分配 `ws.I_idx`、`ws.J_idx`、`ws.K_vals`（容量 = n_coh × 64）
  - 每轮 Newton 开始时 `empty!(ws.I_idx)` 等（保留内存但不重新分配）
  - 单元核内用 `push!` 填充（无需动态增长，容量已预留）
  - 最终用 `sparse()` 构造（`sparse` 内部会复制，这是不可避免的）

- [ ] **Step 4: 利用几何缓存驱动局部计算**
  - 从 `geom_cache[i]` 读取 `R`、`dofs`、`nodes`、`length`
  - 消除 `elem.nodes_bottom/top` 查找和 `dx/dy/L` 计算
  - Gauss 循环从缓存的 `gauss_wts/pts` 读取

- [ ] **Step 5: 保持数值结果不变**
  - 不改变损伤更新顺序
  - 不改变牵引-分离关系（`bilinear_traction_state` / `bilinear_tangent` 不改）
  - 不改变边界约束处理方式

---

### Task 6: 视情况引入线程并行（可选后置）

只在单线程缓存化完成且结果稳定后，再考虑多线程。

**Files:**
- Modify: `src/czm.jl`

- [ ] **Step 1: 评估线程收益**
  - 先运行优化后 benchmark，确认单线程瓶颈位置
  - 如果单线程已经 < 6s，可能不需要线程化

- [ ] **Step 2: 拆分线程局部工作区**
  - 每线程独立的 `K_e` / `f_e` / triplet 缓冲区
  - 线程间只在最终合并阶段交互

- [ ] **Step 3: 验证线程版不会破坏结果**
  - 与单线程输出逐项对比
  - 检查损伤状态单调性

---

## Chunk 4: Regression Verification

### Task 7: 回归与性能验收

**Files:**
- Run: `example/testexample.jl`
- Run: `example/coupled_czm_thermal_example.jl`

- [ ] **Step 1: 跑功能回归**
  - 关注电压、温度、损伤、牵引/分离输出
  - 对比优化前后的 `D_max`、`D_mean`、`δ_max_n` 值

- [ ] **Step 2: 跑性能对比**
  - 单次 CZM 更新时间（优化前后）
  - `assemble_czm_system` 分配数（优化前后）
  - 总体仿真时间

- [ ] **Step 3: 检查无回归**
  - 结果误差在可接受范围内（相对误差 < 1e-8）
  - 确认阶段末无重复更新（验证 CycleSolver 注释有效）
  - 代表性样例全部通过

---

## 验收标准

1. `update_czm_damage!` 单次耗时较当前实现显著下降，优先目标 < 8s。
2. `assemble_czm_system` 和 `assemble_coupled_system` 的临时分配数明显减少。
3. ~~`CycleSolver` 不再做重复的阶段末 CZM 更新。~~ **已由上一轮修复确认。**
4. `example/testexample.jl` 和 `example/coupled_czm_thermal_example.jl` 正常运行。
5. 损伤、位移、牵引力、终止行为与当前实现一致或在可接受误差内。
6. 单线程版本先稳定，再考虑线程并行。

---

## 备注

如果实施过程中发现"缓存 + 工作区"已经足以把时间压到目标范围，就不要为了形式上的向量化继续扩张改动范围。CZM 的第一原则是减少每次更新必须做的工作，而不是把每一行都改成广播。
