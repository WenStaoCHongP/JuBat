# Findings & Decisions

## Requirements
- 保持实体（固体力学）网格节点与拓扑完全不变。
- 仅对界面内聚力网格做单独加密，不允许改动 bulk mesh 的节点布局。
- 加密后的边缘节点必须与原实体表面节点严格共享，不能错位。
- 仅允许在单个实体单元的表面内细分，不能跨实体单元边界做加密。
- 加密后单元长宽比应受控，建议目标小于 5~10。
- 一个 cohesive zone length 内至少覆盖 3~5 个内聚力单元。
- 需要能接入当前 Jellyroll / CZM 网格敏感性验证流程。

## Research Findings
- `src/czm.jl:create_czm_mesh(thermal_mesh, param_dim; tol)` 当前是基于热网格坐标重合来识别界面，并直接用原始节点生成内聚力单元。
- `src/SetMesh.jl:CohesiveMesh` 已经包含 `bulk_mesh`、`bulk_element`、`node_map`、`interface_nodes`、`damage_states` 等字段，适合作为局部加密后的承载结构。
- `src/CouplingState.jl` 目前是 CZM 耦合/状态的统一归属点，`CzmLayout` 和 `compute_czm_effective_params` 都在这里，后续不应额外引入薄包装层。
- `example/网格敏感性/4_czm_mesh_sensitivity.jl` 当前是按 cohesive zone length 反推候选 `nθ`，说明验证路径已经具备“按目标尺度做网格对比”的基础。
- 现有记忆显示 Jellyroll 预设是当前可用于 CZM 归一化检查的有效参数集，且 `chooseCell("Jellyroll")` 比 `LG M50` 更稳妥。
- 现有 CZM 收敛记忆提示 Newton/CZM 过程存在收敛下限与损伤切线不连续问题，因此新的加密方案需要把“几何正确性”和“求解稳健性”分开验证。

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| 采用”界面参数化 + 局部细分 + 节点共享”方案 | 最符合”实体不动、只细分 cohesive”的目标 |
| 把加密限制在单个界面段内部 | 避免跨实体单元边界导致拓扑被破坏 |
| 用 cohesive zone length 驱动目标单元尺寸 | 直接满足每个 CZL 至少 3~5 个单元的要求 |
| 继续复用现有 `CohesiveMesh` 作为输出容器 | 降低对求解器和后处理的扰动 |
| 新建 `CZMMesh.jl` 承载网格生成与加密 | 将网格拓扑与力学组装解耦 |
| 三种加密模式：default / multiple / auto | 覆盖从简单到复杂的全部使用场景 |

## CZM 加密模式设计

### Option 新增字段（`src/Option.jl`）

```julia
# CZM mesh refinement
czm_mesh_refine::String = “default”   # “default” | “multiple” | “auto”
czm_refine_factor::Int64 = 2          # “multiple” 模式的固定倍数（每个原始 CZM 段 → n 个子段）
```

### 三种模式详细算法

#### 模式 1：`default`（默认，当前实现）
- 行为：1 个热单元 ↔ 1 个 CZM 单元
- 算法：`create_czm_mesh` 当前逻辑不变
- 输入：`thermal_mesh` + `param_dim`
- 输出：CZM 单元数 = 界面节点对数 - 1

#### 模式 2：`multiple`（固定倍数加密）
- 行为：1 个原始 CZM 段 → `n = czm_refine_factor` 个子段
- 算法：
  1. 先调用 `create_czm_mesh` 生成原始 1:1 CZM 网格
  2. 对每个原始 CZM 单元：
     - 在 bottom edge 和 top edge 上各插入 `n-1` 个等距中间节点
     - 中间节点坐标 = 原始边缘节点的线性插值
     - 生成 `n` 个子 CZM 单元，每个子单元的长度 = 原始长度 / n
  3. 更新 `damage_states`：每个子单元一个独立的 `DamageState()`
  4. 更新 `node_map`：新增中间节点到热节点的映射关系
- 示例：`czm_refine_factor = 3`，原始 80 个 CZM → 加密后 240 个 CZM

#### 模式 3：`auto`（自适应加密）
- 行为：自动确定每个原始 CZM 段的等分数，使子单元长度 ≤ lc/3
- 算法：
  1. 先调用 `create_czm_mesh` 生成原始 1:1 CZM 网格
  2. 计算内聚力区特征长度：
     ```
     lc = Gc · E_eff / σ_max²
     ```
     参数来源：
     - `Gc = param.cohesive.G_c_n`（归一化断裂能）
     - `E_eff = (E_NE·t_NE + E_PE·t_PE) / (t_NE + t_PE)`（归一化有效模量）
     - `σ_max = param.cohesive.σ_max_n`（归一化最大牵引力）
     - **注意**：此处全部使用归一化参数，结果 `lc` 也是无量纲的
  3. 对每个原始 CZM 单元：
     ```
     target_length = lc / 3                         # 内聚力网格收敛准则：至少 3 个单元覆盖 lc
     n_i = max(1, ceil(L_i / target_length))       # 该单元的等分数
     ```
  4. 其余与 `multiple` 模式相同（等距插入节点、生成子单元）
- 特点：不同 CZM 段可能有不同的 `n_i`（因为各段长度可能不同）
- 最小 `n_i = 1`（当原始单元已经足够短时）

### 节点编号规则

- **边缘节点**：直接复用热网格原始节点编号（不新增）
- **中间节点**：追加到节点表末尾，编号从 `nnode + 1` 开始
- **节点坐标**：沿原始 CZM 段的 bottom/top edge 线性插值
- **CZM 网格节点总数**：`nnode_czm = nnode_thermal + Σ(n_i - 1) × 2`（每个等分段插入 n-1 个节点，上下各一）

### 文件结构规划

```
src/CZMMesh.jl   ← 新建
├── struct CohesiveElement          ← 从 czm.jl 迁移
├── struct DamageState              ← 从 czm.jl 迁移
├── function create_czm_mesh        ← 从 czm.jl 迁移（原始 1:1 网格生成）
├── function compute_cohesive_zone_length  ← 新增：计算 lc
├── function refine_czm_mesh_multiple!     ← 新增：固定倍数加密
├── function refine_czm_mesh_auto!         ← 新增：自适应加密
└── function create_czm_mesh_refined       ← 新增：统一入口（根据 opt 分发）

src/czm.jl   ← 保留力学组装
├── assemble_czm_system
├── assemble_bulk_stiffness
├── assemble_thermal_chemical_load
├── build_czm_cache
├── ensure_czm_cache
├── assemble_coupled_system
├── assemble_coupled_system_full
├── apply_bc_czm
└── identify_bc_nodes_czm
```

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| `memory view` 的 `view_range` 写成了超出实际行数的范围 | 改为全文读取，随后成功获取记忆内容 |

## Resources
- `src/SetMesh.jl`
- `src/czm.jl`
- `src/CouplingState.jl`
- `example/网格敏感性/4_czm_mesh_sensitivity.jl`
- `md/02_几何与网格.md`
- `md/06_内聚力模型_CZM.md`
- `md/08_逐单元算法.md`
- `md/13_耦合验证方案.md`

## CZM 加密对多物理场耦合的影响分析

### 当前架构：1:1 映射（一个 CZM 单元 ↔ 一个热单元）

#### 1. CZM 网格生成
- `create_czm_mesh` (`src/czm.jl:54`) 基于热网格坐标重合检测界面节点对
- 每对相邻界面节点之间生成一个 CZM 单元，直接复用热网格节点
- `nθ` 个热单元 → `nθ - 1` 个 CZM 单元

#### 2. 电化学反应停止机制（D=1 → 热源 = 0）

**路径 A：热源置零** (`compute_heat_sources_with_czm` → `get_active_elements`)
- `get_active_elements` (`src/Materialmatrix.jl:345`) 对每个内层热单元检查其关联的 CZM 单元
- 通过 `czm_element_map`（`src/CouplingState.jl:91`，`Dict{Int,Vector{Int}}`）查找关联
- 若任一 CZM 单元断裂（`D ≥ 0.99` 或 `fractured = true`）→ 该热单元标记为非活跃
- 非活跃热单元的 `heat_source_fields[e]` 被设为 0（`src/ThermalDistributed.jl:538`）

**路径 B：分流电流置零** (`solve_branch_currents` + `deactivated_elements`)
- `CallModel.jl:62-70` 从 `variables["deactivated_elements"]` 读取失效热单元列表
- `solve_branch_currents` (`src/Parallelsolution.jl:358`) 将失效热单元的电流强制为 0
- 总电流由剩余活跃热单元重新分配

**当前 1:1 行为**：一个 CZM 单元 D=1 → 其唯一对应的热单元完全停止电化学反应（q=0, I=0）

#### 3. 间隙热导模型
- `ThermalDistributed2D_BC` (`src/ThermalDistributed.jl:292-309`) 遍历每个 CZM 单元
- 对每个 CZM 单元：`coeff = h_eff × elem.length`
- 修改热刚度矩阵的节点对耦合项：`K[nb,nb] -= coeff`, `K[nb,nt] += coeff`, ...
- **关键**：CZM 单元的 `nodes_bottom` / `nodes_top` 就是热网格节点，贡献直接进入热系统矩阵

#### 4. CZM 应变输入
- `compute_czm_strain_inputs` (`src/CouplingState.jl:259`) 基于 `czm_mesh.bulk_element`（热单元）
- 温度和 SOC 按热单元平均 → 作为 CZM 载荷输入
- **不受 CZM 加密影响**，因为 bulk mesh 不变

---

### 加密后架构：N:1 映射（多个 CZM 单元 ↔ 一个热单元）的影响

#### 影响一：电化学反应停止策略（关键，需要修改）

| 维度 | 当前 1:1 | 加密后 N:1 |
|------|----------|------------|
| 粒度 | 1 个 CZM D=1 → 整个热单元停止 | 1/N 个 CZM D=1 → 整个热单元仍然停止 |
| 问题 | 物理一致 | **过于激进**：热单元只有部分界面受损 |
| 物理真实 | 全有/全无 | 应支持渐变：部分损伤 → 部分降低电化学活性 |

**具体分析**：
- `get_active_elements` 已支持 N:1 的 `czm_element_map` 遍历（`for czm_idx in get(mesh_data.czm_element_map, e, Int64[])`）
- 但判断逻辑是 `if czm_idx in fractured_czm → active[e] = false`，即任一 CZM 断裂则整个热单元失效
- 加密后，一个热单元可能对应 3~5 个 CZM 单元，其中只有 1 个断裂时就会导致整个热单元停止产热
- **这与物理不符**：只有热单元界面的一部分失去电接触，大部分界面仍应保持电化学活性

**需要引入的机制（建议三选一）**：
1. **面积加权热源衰减**：`q_eff[e] = q[e] × (1 - D_avg[e])`，其中 `D_avg` 是关联 CZM 单元损伤的长度加权平均
2. **面积加权电流衰减**：`I_max[e] = I_original[e] × active_fraction`，其中 `active_fraction = 1 - n_fractured / n_total`
3. **阈值策略**：只有当 `n_fractured / n_total ≥ threshold`（如 0.5）时才完全停止热单元

**推荐方案 1（面积加权热源衰减）**，因为：
- 保持分流求解器的 `active_mask` 逻辑简单
- 物理上最直观：界面损伤比例直接折减电化学产热
- 修改范围最小：只需修改 `compute_heat_sources_with_czm` 中的热源折减逻辑

#### 影响二：间隙热导模型（需要节点映射）

| 维度 | 当前 1:1 | 加密后 N:1 |
|------|----------|------------|
| 节点 | CZM 节点 = 热网格节点 | CZM 引入新的细分节点 |
| 贡献 | 直接写入热刚度矩阵 | 需要映射到热网格节点 |
| 总导热 | 1 个 `h_eff × L` | N 个 `h_eff_i × L_i`，理论上 `Σ(h_i × L_i) = h_total × L` |

**具体分析**：
- 加密后的 CZM 单元引入新的中间节点（不在热网格中）
- `ThermalDistributed2D_BC` 中 `czm_elem.nodes_bottom/top` 可能不再是热网格节点
- 需要将每个加密 CZM 单元的间隙导热贡献映射到最近的热网格节点对上
- 如果损伤均匀，加密前后的总导热守恒（`Σ L_i = L_original`）
- 加密后可以更精确地描述**损伤梯度**沿界面的变化

**需要引入的机制**：
- 为每个加密 CZM 单元建立到热网格节点对的映射关系
- 可以在 `CohesiveElement` 中增加 `thermal_node_pair` 字段，或使用父单元的信息
- 或者：在 `ThermalDistributed2D_BC` 中按热网格节点对聚合多个 CZM 单元的导热贡献

#### 影响三：`czm_element_map` 构建（需要更新）

| 维度 | 当前 | 加密后 |
|------|------|--------|
| 映射 | ~1:1 | 真正 N:1 |
| 构建位置 | `Jellyrollmodel.jl:117-144` | 需要更新以包含加密后的 CZM 索引 |

**具体分析**：
- 当前构建逻辑通过 CZM 单元的节点查找关联的热单元
- 加密后的 CZM 单元引入新节点，但边缘节点仍与热网格共享
- 可以通过边缘节点或通过父 CZM 单元索引建立映射
- `get_active_elements` 的遍历逻辑本身已支持 N:1，不需要修改

#### 影响四：CZM 应变输入（无需修改）

- `compute_czm_strain_inputs` 基于 `czm_mesh.bulk_element`（热单元）
- Bulk mesh 不变，因此温度和 SOC 输入不受影响
- 但加密后的 CZM 单元可能需要**亚单元级别**的应变梯度（当前模型未实现）
- 短期：无需修改；长期：可考虑将热单元温度场插值到 CZM 单元中心

#### 影响五：CZM 力学求解（自身兼容）

- `assemble_czm_system` 逐 CZM 单元组装，天然支持加密
- `CohesiveElementGeom` 缓存每个单元几何，无需修改
- CZM Newton-Raphson 求解器的 DOF 数量随加密而增加
- **需要注意**：加密后 DOF 数量增加 → 稀疏矩阵规模增大 → 计算成本上升

---

### 修改优先级与依赖关系

```
1. [必须] 更新 czm_element_map 构建逻辑
   ↓
2. [必须] CZM 单元 → 热网格节点对映射（影响间隙热导）
   ↓
3. [必须] 引入部分损伤折减机制（替代全有/全无策略）
   ↓
4. [验证] 验证加密前后物理一致性（均匀损伤时结果应一致）
```

### 代码影响清单

| 文件 | 函数/位置 | 修改类型 | 优先级 |
|------|-----------|----------|--------|
| `src/czm.jl` | `create_czm_mesh` 或新函数 | 加密后需新增节点和 CZM 单元 | P0 |
| `src/Jellyrollmodel.jl:117-144` | `czm_element_map` 构建 | 更新为 N:1 映射 | P0 |
| `src/ThermalDistributed.jl:292-309` | `ThermalDistributed2D_BC` | CZM 节点→热节点映射 | P1 |
| `src/ThermalDistributed.jl:512-550` | `compute_heat_sources_with_czm` | 引入渐变折减 | P1 |
| `src/Materialmatrix.jl:345-363` | `get_active_elements` | 支持渐变活跃度 | P1 |
| `src/Parallelsolution.jl:358-440` | `solve_branch_currents` | 可选：渐变电流上限 | P2 |
| `src/CouplingState.jl:259-307` | `compute_czm_strain_inputs` | 无需修改（长期可选） | P3 |

## Visual/Browser Findings
-