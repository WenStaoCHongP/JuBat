# Task Plan: 内聚力模块变量与状态重构

## Goal
把内聚力模块的计算结果和中间变量统一收敛到 Variables.jl 定义、PostProcessing.jl 还原的模式，对齐电化学变量状态的设计规范。不新增 CZMState 结构，扩展现有 CouplingState.jl 中的类型。

## Current Phase
Complete

## Phases

### Phase 1: Requirements & Discovery
- [x] 梳理当前 CZM 调用链、状态持有点和数据流
- [x] 确认问题范围：变量键定义分散、命名不一致、预分配缺失
- [x] 将发现记录到 findings.md
- **Status:** complete

### Phase 2: 变量键统一设计
- [x] 2.1 审计全部 CZM 变量键的使用点 — findings.md 对照表
- [x] 2.2 统一键名：czm max damage → czm D_max 等 — `7ec95dd`
- [x] 2.3 补齐 StandardVariables 中 7 个场变量预分配 — `4d2c4c8`
- [x] 2.4 补齐 create_element_workspace 中 CZM 工作区键 — `ab39179`
- [x] 2.5 更新 czm_output_to_variables 使用统一键名 — `7ec95dd`
- [x] 2.6 PostProcessing.jl 还原逻辑已使用正确键名（无需改动）
- **Status:** complete

### Phase 3: 状态传递收敛
- [x] 3.1 新增 CzmLayout (n_coh, ndof, u_prev) — `dd9febf`
- [x] 3.2 Case 新增 czm_layout 字段 — `dd9febf`
- [x] 3.3 update_czm_damage! 从 case.czm_layout 获取 u_prev — `9f4da8b`
- [x] 3.4 收敛后写回 case.czm_layout.u_prev — `9f4da8b`
- **Status:** complete

### Phase 4: 验证与文档
- [x] 模块加载验证通过
- [x] 导出符号验证通过（含 CzmLayout）
- [ ] CLAUDE.md CZM 变量说明待更新
- **Status:** mostly complete

## Key Questions — Resolved

### Q1: 是否需要新增 CZMState 结构？
**决定：不需要。** 扩展现有结构：
- `Case.czm_mesh` → 网格+损伤（已有）
- `Case.czm_cache` → 装配缓存（已有）
- `Case.czm_layout` → 布局+u_prev（新增，轻量 struct）
- 计算结果和中间变量命名在 Variables.jl 定义
- PostProcessing.jl 负责还原物理单位

### Q2: 状态分层？
| 层 | 归属 | 持有者 | 状态 |
|----|------|--------|------|
| 网格+损伤 | 长期 | `Case.czm_mesh` | 已有 |
| 装配缓存 | 长期 | `Case.czm_cache` | 已有 |
| 布局+u_prev | 长期 | `Case.czm_layout` | 新增 |
| 装配工作区 | 步内 | `CZMAssemblyCache.ws` | 已有 |
| 位移+收敛标记 | 步内 | `CZMResult` | 已有 |
| 统计+场输出 | 步内 | `variables` Dict | 已有 |

### Q3: 字符串键？
**决定：全部保留在 Dict 中，但统一命名并集中预分配。**
- Variables.jl 定义所有键名和维度
- czm_output_to_variables 写入 variables Dict
- PostProcessing.jl 从 Dict 读取并还原物理单位
- 不做 typed struct 迁移，保持与电化学变量一致的 Dict 模式

### Q4: 刚度矩阵重建触发条件？
**决定：E/ν 变化时重建。** 已有 `ensure_czm_cache` 处理。不需要额外设计。

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 不新增 CZMState，扩展 CouplingState.jl | 现有结构已覆盖大部分需求，只缺 CzmLayout |
| 计算结果在 Variables.jl 定义，PostProcessing.jl 还原 | 与电化学变量状态设计对齐 |
| 字符串键保持 Dict 模式，统一命名 | 与电化学 Dict 模式一致，避免 typed struct 迁移的破坏性 |
| u_prev 归入 CzmLayout 而非松散传递 | 收敛 update_czm_damage! 签名 |
| 保留现有求解入口，优先做适配 | 降低重构风险 |

## Protected Constraints (不可破坏)
1. PostProcessing.jl 中的物理单位还原公式（scale.L, scale.E_n, scale.E_p, scale.r0）
2. czm_output_to_variables 的输出语义（哪些量是 per-element、哪些是 scalar）
3. update_czm_damage! 的损伤提交门控（只在收敛时提交）
4. Variable_update! 的动态扩展机制

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| memory.view_range 超出实际行数 | 1 | 改为读取正确范围 |

## Notes
- 核心改动集中在 Variables.jl、CzmPostProcess.jl、PostProcessing.jl 和 CouplingState.jl
- CzmLayout 只需 ~5 行 struct 定义 + ~10 行初始化逻辑
- 变量键统一是一次性改动，不影响求解器内部逻辑
