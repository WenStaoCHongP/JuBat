# CZM 求解器向量化与缓存优化设计规格

> 日期: 2026-04-20
> 状态: 实施中
> 基于文档: `docs/代码优化/15_CZM耦合修复与单步损伤演化.md`
> 关联计划: `docs/superpowers/plans/2026-04-20-czm-vectorized-solver-plan.md`
> **Updated:** 2026-04-20 — 基于代码审查修正，补充量化数据

---

## 1. 目标

针对 CZM（内聚力模型）在耦合求解中的性能瓶颈，设计一条可落地的优化路线，目标是在不改变物理模型和数值结果的前提下，显著降低 `update_czm_damage!` / `solve_czm_step` 的单次耗时。

**核心判断**：这里的"向量化"不应理解为把外层 Newton 循环简单改成广播，而应理解为把单元几何、局部算子、输入提取和临时数组分配改造成批处理、预分配、缓存化的数据流。外层 Newton / load substep 依赖关系必须保留串行。

**建议目标**：
- 第一阶段把单次 CZM 更新从约 12s 压到 6-8s 量级
- 第二阶段在缓存和工作区复用稳定后，进一步把耗时压到 6s 以下
- 若后续引入线程并行，再争取更低的单次更新时间

---

## 2. 热点拆解

### 2.1 调用链

当前最热路径是：

`Solve.jl` → `update_czm_damage!` → `solve_czm_step` → `newton_raphson_czm` → `assemble_coupled_system`

其中：
- `Solve.jl` (L270-285) 决定更新频率
- `update_czm_damage!` 负责把温度/SOC 输入转换成 CZM 求解所需参数
- `solve_czm_step` 负责一次损伤更新的外层求解
- `newton_raphson_czm` 负责 Newton + load substep 迭代
- `assemble_coupled_system` 负责每次迭代的刚度和内力组装

### 2.2 主要耗时来源

| 位置 | 类型 | 结论 |
|------|------|------|
| `assemble_bulk_stiffness` | **重复计算（最高 ROI 缓存项）** | 每次 Newton 迭代都重算，但 E/ν 不变时 K_bulk 相同。~20 load_steps × ~3 Newton = ~60 次冗余 bulk 重组 |
| `assemble_czm_system` | **单元级循环（核心热区）** | 801 cohesive 单元 × 2 Gauss 点 × ~50 总迭代 ≈ ~80K 次临时数组分配/更新 |
| `compute_czm_strain_inputs` | 输入准备 | 当前仍按单元逐个取平均，属于可批量化输入层 |
| ~~`CycleSolver.solve_phase`~~ | ~~冗余调用~~ | **已在上一轮修复（CycleSolver.jl:157-160）** |

### 2.3 性能特征

- CZM 计算的成本主要来自大量小矩阵、小向量和重复积分，而不是单个大矩阵运算。
- 纯 broadcast 并不能自动解决 sparse assembly 和 Newton 迭代依赖。
- 真正有效的优化应优先减少临时分配、缓存静态量、复用局部工作区。

**临时分配量化**（每次 `update_czm_damage!` 调用）：
- `assemble_bulk_stiffness`：~800 bulk 元素 × `zeros(8,8)` + `zeros(3,8)` × Gauss 点 × ~60 轮 ≈ ~96K 次
- `assemble_czm_system`：801 单元 × (`zeros(8,8)` + `zeros(8)` + `zeros(2,8)` × 2 Gauss) × ~60 轮 ≈ ~288K 次
- triplet `push!`：~801 × 64 entries × ~60 轮 ≈ ~3M 次
- 总计 ~400K 次临时分配/更新，是主要 GC 压力来源

---

## 3. 优化方向

### 3.1 静态缓存优先（最高 ROI：K_bulk 缓存）

**K_bulk 缓存是整个优化中收益最高的单点**，因为：
- `assemble_bulk_stiffness` 在每次 Newton 迭代都被调用
- 但 `K_bulk` 在同一次 CZM 更新中完全不变（E/ν 不变、bulk mesh 不变）
- 缓存后可消除 ~60 次完整的 Q4 FEM 组装

优先缓存以下静态或准静态内容：
- **bulk stiffness 矩阵**（最高优先）
- cohesive 单元到全局 DOF 的映射
- 单元几何信息（长度、切向/法向、节点对）
- Gauss 积分常量和形函数常量
- 边界条件节点/自由度列表

**判断**：这是最高收益、最低风险的第一步，因为这些量在一次 CZM 更新中通常不变。

### 3.2 批处理单元核

将 `assemble_czm_system` 从"单元内不断分配临时数组"的模式，改为"预分配工作区 + 批处理单元核"的模式：

- 预先准备 `u_e`、`B`、`K_e`、`f_e` 等工作区
- 用 `@views`、`@inbounds`、`mul!` 等手段减少临时对象
- 将每个单元的局部计算尽量写成无分配的紧凑内核
- sparse triplet 的构建也尽量改成预分配缓冲区（`empty!` 复用）

**判断**：这是最关键的性能改造点，也是"向量化"最应落地的地方。

### 3.3 输入提取批量化

`compute_czm_strain_inputs` 现在逐单元遍历节点并做平均。该步骤本身不是最大热区，但它是所有 CZM 更新的入口，适合一起整理：

- 温度输入统一为单元级数组
- SOC 输入统一为单元级数组
- `epsilon_0_elem` 在一次更新中只算一次
- 输出保持为数组，不做不必要的 flatten

**判断**：收益中等，但能显著简化后续核函数接口。

### 3.4 求解频率与调用链

- `czm_update_interval` 仍然保留，作为控制求解频率的低风险手段
- `czm_load_steps` 可作为精度/速度折中参数
- ~~必须清理 `CycleSolver.solve_phase` 的阶段末冗余更新~~ → **已完成**

**判断**：这是控制总成本的辅助杠杆，不应替代核心内核优化。

### 3.5 可选线程并行

如果单线程缓存化和批处理后仍然过慢，可以考虑对 cohesive 单元层引入 `Threads.@threads`。

**前提**：
- 必须先把共享的 `push!` / sparse triplet 逻辑改成线程安全模式
- 必须把工作区拆成线程局部
- 必须验证并行不会破坏损伤状态更新顺序

**判断**：这是后续选项，不应作为第一步。

---

## 4. 推荐架构

建议把优化拆成两个层次：

### 4.1 `CZMAssemblyCache`

挂载在 `Case.czm_cache::Union{Nothing, CZMAssemblyCache}`，用于保存静态或准静态数据：
- `K_bulk::SparseMatrixCSC{Float64, Int64}`
- `bulk_dofs::Vector{Vector{Int64}}` — 每个 bulk 元素的 DOF 映射
- `cohesive_geom::Vector{CohesiveElementGeom}` — 预计算几何
- `bc_dofs::Vector{Int64}`, `bc_vals::Vector{Float64}`
- `E_eff::Float64`, `ν_eff::Float64` — 缓存失效判断
- `valid::Bool`

### 4.2 `CZMAssemblyWorkspace`

每次 Newton 迭代复用，用于保存可复用的临时数组：
- `u_e::MVector{8, Float64}`
- `K_e::MMatrix{8, 8, Float64}`
- `f_int_e::MVector{8, Float64}`
- `B_global::MMatrix{2, 8, Float64}`
- `I_idx::Vector{Int64}`, `J_idx::Vector{Int64}`, `K_vals::Vector{Float64}` — triplet（预分配容量，`empty!` 复用）
- `f_int_coh::Vector{Float64}`

### 4.3 `CohesiveElementGeom`

每个 cohesive 单元的预计算几何：
- `length::Float64`
- `n_vec::SVector{2, Float64}`, `t_vec::SVector{2, Float64}`
- `R::SMatrix{2, 2, Float64}`
- `dofs::SVector{8, Int64}`, `nodes::SVector{4, Int64}`
- `gauss_wts::Vector{Float64}`, `gauss_pts::Vector{Float64}`

**设计原则**：
- 缓存负责"不变的东西"（`CZMAssemblyCache`）
- 工作区负责"会变但可以复用的东西"（`CZMAssemblyWorkspace`）
- 外层 Newton 迭代保留串行，内层单元核尽量无分配

---

## 5. 风险与约束

| 风险 | 说明 | 缓解 |
|------|------|------|
| 缓存失效 | mesh 或材料参数变动后缓存过期 | 比较 E/ν 值，不匹配则重建 |
| 过度向量化 | broadcast 生成大量临时数组，反而更慢 | 优先手写无分配循环和 workspace 复用 |
| 稀疏装配并行冲突 | 多线程写入 triplet 或矩阵会冲突 | 先单线程优化，再考虑线程局部缓冲区 |
| ~~结果重复更新~~ | ~~`Solve` 与 `CycleSolver` 双重更新损伤~~ | **已在上一轮修复** |
| 数值漂移 | 改写核函数后误差累积 | 以现有 example 和 damage 输出做回归比较 |

---

## 6. 验证标准

- `example/testexample.jl` 能正常运行
- `example/coupled_czm_thermal_example.jl` 能正常运行
- `update_czm_damage!` 单次耗时明显下降
- `czm D_max`、`czm mean damage`、位移和牵引力结果保持在可接受误差内（相对误差 < 1e-8）
- ~~阶段末不再出现重复 CZM 更新~~ → **已确认无重复**
- 热点函数分配数下降，尤其是 `assemble_czm_system` 相关分配

---

## 7. 结论

CZM 的性能瓶颈更像一个"小型 FEM 组装器"问题，而不是单纯的循环写法问题。最佳路线不是把外层求解器强行广播化，而是：

1. **先缓存 K_bulk**（最高 ROI 单点，消除 ~60 次冗余 bulk 重组）
2. 再缓存几何/DOF 映射
3. 再把单元核改成预分配批处理
4. 最后再考虑线程并行

这条路线风险最低，也最符合当前代码结构。

**不需要优化的部分**：
- `assemble_thermal_chemical_load`（仅在 Newton 入口调用一次，且 load_factor 缩放已正确）
- `bilinear_traction_state` / `bilinear_tangent`（本构模型，不改）
