# Task Plan: 内聚力网格单独加密

## Goal
在保持实体（固体力学）网格节点与拓扑完全不变的前提下，实现仅对界面内聚力（Cohesive）网格的单独加密，同时保证节点一一共享、单元质量可控、且不跨实体单元边界。加密后需正确处理 N:1（多个 CZM → 一个热单元）的电化学反应停止和间隙热导映射。

## Current Phase
Phase 1

## Phases

### Phase 1: Requirements & Discovery
- [x] 理解用户目标与约束
- [x] 识别当前 CZM 网格生成入口与数据结构
- [x] 记录已知风险与验证约束
- [x] 分析 CZM 加密对多物理场耦合的影响（详见 findings.md）
- **Status:** completed

### Phase 2: Algorithm Design — CZM 加密模式
- [ ] 设计三种加密模式的具体算法：
  - `default`：当前 1:1 实现，直接从热网格界面节点生成 CZM 单元
  - `multiple`：固定倍数 n，将每个原始 CZM 段等分为 n 个子段
  - `auto`：自适应等分，计算 lc = Gc·E / σ_max²，然后 n = ceil(L_elem / (lc/3))
- [ ] 定义节点拆分规则：
  - 边缘节点必须与热网格节点严格共享
  - 中间节点在原 CZM 段的弧线上等距插入
  - 中间节点仅存在于 CZM 网格，不加入热网格
- [ ] 设计 CZM→热节点映射策略（间隙热导装配用）
- [ ] 设计部分损伤折减机制（替代全有/全无策略）
- **Status:** pending

### Phase 3: File Restructuring — 新建 CZMMesh.jl
- [ ] 新建 `src/CZMMesh.jl`
- [ ] 从 `src/czm.jl` 迁移以下内容至 `CZMMesh.jl`：
  - `CohesiveElement` 结构体定义
  - `DamageState` 结构体定义
  - `create_czm_mesh` 函数（原始 1:1 网格生成）
- [ ] 在 `CZMMesh.jl` 中新增：
  - `refine_czm_mesh_default!`：default 路径（不加密，直接返回原始 mesh）
  - `refine_czm_mesh_multiple!`：倍数路径（固定 n 等分每个 CZM 段）
  - `refine_czm_mesh_auto!`：自适应路径（基于 lc 自动确定 n）
  - `create_czm_mesh_refined`：统一入口，根据 opt 分发到上述三者
  - `compute_cohesive_zone_length`：计算 lc = Gc·E / σ_max²
- [ ] 更新 `src/JuBat.jl`：
  - 将 `include("czm.jl")` 改为先 `include("CZMMesh.jl")` 再 `include("czm.jl")`
  - 或将 `czm.jl` 中已迁移部分删除，仅保留力学组装函数
  - 更新 export 列表
- [ ] `src/czm.jl` 保留：`assemble_czm_system`、`assemble_bulk_stiffness`、`assemble_thermal_chemical_load`、`build_czm_cache`、`ensure_czm_cache`、边界条件函数等力学相关代码
- **Status:** pending

### Phase 4: Option 扩展
- [ ] 在 `src/Option.jl` 的 `Option` 结构体中新增字段：
  ```julia
  czm_mesh_refine::String = "default"  # "default" | "multiple" | "auto"
  czm_refine_factor::Int64 = 2         # "multiple" 模式的固定倍数
  ```
  - auto 模式的目标长度 `lc/3` 直接在程序中内联，不作为用户选项
- [ ] 在 `CLAUDE.md` 的配置选项表中补充说明
- **Status:** pending

### Phase 5: Data Model & API Plan
- [ ] 规划 `CohesiveMesh` 需要新增或调整的字段：
  - `parent_czm_idx::Vector{Int}` — 每个加密 CZM 单元的父单元索引（default 模式下 = 自身索引）
  - `thermal_node_pairs::Vector{Tuple{Int,Int}}` — 每个加密 CZM 单元对应的热网格节点对（用于间隙热导装配）
  - `refine_mode::String` — 记录使用的加密模式
- [ ] 规划 `czm_element_map` 更新策略：在 `Jellyrollmodel.jl` 中扩展以包含加密后的 CZM 索引
- [ ] 规划 `get_active_elements` → `get_element_activity` 渐变接口
- [ ] 约束与现有 `CouplingState.jl`、`CzmLayout` 的接口兼容
- **Status:** pending

### Phase 6: Verification Plan
- [ ] 定义节点重合、拓扑不变、局部细分不跨单元的检查
- [ ] 定义长宽比与单元数量的验收标准
- [ ] 加密前后一致性验证：均匀损伤下 default/multiple/auto 结果应一致
- [ ] auto 模式验证：确认 lc 计算正确，单元数满足 lc/3 约束
- [ ] 渐变损伤场景测试：部分 CZM 断裂时热源/电流应渐变折减
- [ ] 间隙热导聚合验证：多个加密 CZM 的导热贡献之和等于原始值
- [ ] 扩展 `example/网格敏感性/4_czm_mesh_sensitivity.jl`
- **Status:** pending

### Phase 7: Delivery Plan
- [ ] 输出推荐的实现顺序
- [ ] 输出风险点、回滚方案与后续落地步骤
- **Status:** pending

## Key Questions
1. ~~单独加密是否以每个界面独立的目标单元长度或 `nθ` 输入为主？~~ → 已确定三种模式：default / multiple / auto
2. 节点拆分后，`node_map` 是否需要支持原节点到多个派生共享节点的映射？
3. ~~现有 `create_czm_mesh` 是否直接扩展？~~ → 已确定拆分到 CZMMesh.jl
4. 部分损伤折减采用哪种策略？（面积加权热源衰减 / 面积加权电流衰减 / 阈值策略）
5. 加密后的 CZM 中间节点如何映射到热网格节点对？

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 保持 bulk mesh 完全不变 | 满足实体网格节点与拓扑不变的硬约束 |
| 只在界面内聚力段上做局部细分 | 避免重网格实体并保持边界共享 |
| 细分必须局限于单个实体单元表面内 | 防止跨实体单元加密破坏拓扑 |
| 以 cohesive zone length 和长宽比约束驱动目标尺寸 | 保证每个 CZL 内有足够单元且质量可控 |
| 引入渐变损伤折减替代全有/全无策略 | N:1 映射下部分 CZM 断裂不应杀死整个热单元 |
| 三种加密模式：default / multiple / auto | 覆盖从简单到复杂的全部使用场景 |
| 新建 CZMMesh.jl 承载网格生成与加密 | 将网格拓扑与力学组装解耦，职责清晰 |
| czm.jl 仅保留力学组装函数 | 遵循单一职责原则 |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| `memory view` 的 `view_range` 超出范围 | 1 | 改为读取记忆文件全文或使用更小范围 |

## Notes
- 现有实现锚点：`src/czm.jl`（网格生成 + 力学组装）、`src/SetMesh.jl`（CohesiveMesh 结构）、`src/Option.jl`（选项）
- **加密对多物理场的完整影响分析已写入 `findings.md`**
- `lc` 计算公式：`lc = Gc·E / σ_max²`，参数来源 `param.cohesive.G_c_n`（归一化）或 `param_dim.cohesive.G_c_n`（物理量）和 `E_eff`
