# JuBat src/ 简化实施计划索引

> **For agentic workers:** REQUIRED: 使用 `superpowers:subagent-driven-development` 或 `superpowers:executing-plans` 执行本目录下的 plan。每个 plan 文件独立可执行，按下方"执行顺序"推进。

**Goal:** 依据 spec `docs/superpowers/specs/2026-08-04-src-simplification-design.md`，对 JuBat `src/*.jl` + `src/parameters/*.jl` 共 41 个文件执行 Reduce 模式简化。

**Architecture:** Strategy B 依赖层级自底向上（参数/变量 → 几何/网格 → 物理模型 → 求解器 → 后处理）。每个文件独立 plan，按层分目录。索引文件聚合状态与执行顺序。

**Tech Stack:** Julia 1.x；项目内部测试框架（`test/unit_czm_*.jl`）；grep 脚本（`Simplify/scripts/`）。

---

## 执行顺序（推荐）

依据 spec 附录 A，按"低风险预热 → 重复消除 → 兜底清理 → dead code 批量"分批。**每个文件独立 commit；每个 PR 不跨层**。

| 阶段 | 文件 | 优先级 |
|---|---|---|
| 0. 基线固化 | `Simplify/baseline.md`（脚本生成） | 必须先做 |
| 1. 低风险预热 | install / SetCase 兼容入口 / parameters dead 验证 | 风险★ |
| 2. 第 1 层完成 | VariableKeys 新建 → Variables / CallModel / CouplingState 键集中化 | 风险★★ |
| 3. 第 5 层主菜 | CycleData D1 删除 + CycleSolver callback 化 | 风险★★★ |
| 4. 第 3 层重复消除 | ThermalDistributed D3/D4 合并 | 风险★★★ |
| 5. 兜底清理 | 类别 A/B 跨文件处理（CsvExport / Mechanical / Solve / CouplingState） | 风险★★ |
| 6. Dead code 批量 | SPM / P2D / ThermalLumped / ring / ThermalPolar2D 验证后删除 | 风险★ |
| 7. D9 triage | CouplingState overload 评估 + czm.jl 双胞胎 | 风险★★ |
| 8. High-risk（仅在测试到位时） | czm.jl / CzmSolve.jl 内部 diff 确认项 | 风险★★★★ |

---

## 第 1 层：参数与变量（`01-parameters-variables/`）

| Plan 文件 | 源码 | 桶 | 主要动作 |
|---|---|---|---|
| `SetParams.md` | `SetParams.jl` (522) | Leave alone | 不动；归一化已统一 |
| `Variables.md` | `Variables.jl` (277) | Consolidate (D2) | 删两份硬编码键表；接入 `VariableKeys.jl` |
| `VariableKeys.md` | **新建** `src/VariableKeys.jl` | — | 集中所有变量键常量 |
| `Option.md` | `Option.jl` (88) | Leave alone | 验证 solveType forward/backward 是否被走 |
| `parameters-Jellyroll.md` | `parameters/Jellyroll.jl` (202) | Leave alone | 主线必保；TODO 占位标注 |
| `parameters-Enertech.md` | `parameters/Enertech.jl` (175) | Delete 候选 | grep `ChooseCell("Enertech")` 验证 |
| `parameters-LGM50.md` | `parameters/LGM50.jl` (136) | Delete 候选 | grep 验证 |
| `parameters-Northrop.md` | `parameters/Northrop.jl` (133) | Delete 候选 | grep 验证 |
| `parameters-Ring.md` | `parameters/Ring.jl` (65) | Delete 候选 | grep 验证 |
| `CallModel.md` | `CallModel.jl` (268) | Consolidate (D2) | 删 `copy_element_results` 键字面量 |
| `SetCase.md` | `SetCase.jl` (116) | Leave alone | 5 参数兼容构造器 grep 后处理 |

## 第 2 层：几何与网格（`02-geometry-mesh/`）

| Plan 文件 | 源码 | 桶 | 主要动作 |
|---|---|---|---|
| `SetMesh.md` | `SetMesh.jl` (767) | High-risk-leave-alone | 不动；记录已审查 |
| `Jellyrollmodel.md` | `Jellyrollmodel.jl` (690) | High-risk-leave-alone | 不动 |
| `CzmUnitMesh.md` | `CzmUnitMesh.jl` (108) | Leave alone | 保留独立（unit test 隔离合理） |
| `ring.md` | `ring.jl` (70) | Delete 候选 | grep `ring_mesh` 决定删/留 |

## 第 3 层：物理模型（`03-physics-models/`）

| Plan 文件 | 源码 | 桶 | 主要动作 |
|---|---|---|---|
| `SPMe.md` | `SPMe.jl` (290) | Leave alone | 主线核心电化学；仅审查 |
| `SPM.md` | `SPM.jl` (85) | Leave alone | Option 默认模型；活跃 API |
| `P2D.md` | `P2D.jl` (287) | Leave alone | 主线 API；minimal_example 用 |
| `ElectrodeDiffusion.md` | `ElectrodeDiffusion.jl` (18) | Leave alone | 保留独立 |
| `ElectrodePotential.md` | `ElectrodePotential.jl` (19) | Leave alone | 保留独立 |
| `ElectrolyteDiffusion.md` | `ElectrolyteDiffusion.jl` (30) | Leave alone | 保留独立 |
| `ElectrolytePotential.md` | `ElectrolytePotential.jl` (23) | Leave alone | 保留独立 |
| `czm.md` | `czm.jl` (873) | Mixed（Consolidate + Dead code） | 删无调用者 `assemble_coupled_system_full`；审查 cache；数值核心不动 |
| `Mechanical.md` | `Mechanical.jl` (360) | Leave alone（category A 评估） | line 294 try/catch 决策（删或加注释） |
| `Materialmatrix.md` | `Materialmatrix.jl` (427) | Leave alone | 仅审查 |
| `ThermalDistributed.md` | `ThermalDistributed.jl` (557) | Consolidate (D3 + D4) | 双胞胎 non-! 改 wrapper；! 内部回退消除 |
| `Thermal.md` | `Thermal.jl` (80) | Leave alone | Lumped 备用路径；活跃 API |
| `ThermalPolar2D.md` | `ThermalPolar2D.jl` (120) | High-risk-leave-alone | 数值核心；不动 |

## 第 4 层：求解器（`04-solvers/`）

| Plan 文件 | 源码 | 桶 | 主要动作 |
|---|---|---|---|
| `Solve.md` | `Solve.jl` (471) | Mixed（Leave alone + category A/B 评估） | 删 line 418 silent try/catch；保留 line 273 CZM 容错 + 注释 |
| `CzmSolve.md` | `CzmSolve.jl` (675) | High-risk + category D | 4 处 try/catch 全保留 + 标准化注释 |
| `Parallelsolution.md` | `Parallelsolution.jl` (453) | High-risk-leave-alone | 不动 |
| `CycleSolver.md` | `CycleSolver.jl` (546) | Leave alone（D1 keep 侧） | 与 CycleData.md 联动；本文件不动 |
| `Initialisation.md` | `Initialisation.jl` (152) | Leave alone | 仅审查 |
| `Assemble.md` | `Assemble.jl` (40) | High-risk-leave-alone | 数值基础设施；不动 |

## 第 5 层：后处理与导出（`05-postprocessing-export/`）

| Plan 文件 | 源码 | 桶 | 主要动作 |
|---|---|---|---|
| `CycleData.md` | `CycleData.jl` (623) | Consolidate (D1) | wrapper 化（依赖 CycleSolver `export_callback`）；-400 行 |
| `CsvExport.md` | `CsvExport.jl` (636) | Consolidate（category B） | 7 处 try/catch 抽 `safe_write_csv` helper（保留 @warn） |
| `PostProcessing.md` | `PostProcessing.jl` (350) | Leave alone | 仅审查 |
| `CouplingState.md` | `CouplingState.jl` (762) | Consolidate（D9 + category B） | `update_czm_damage!` 双方法 triage；4 处 @warn 加注释 |
| `CzmPostProcess.md` | `CzmPostProcess.jl` (117) | Leave alone | 仅审查 |
| `Tools.md` | `Tools.jl` (191) | Leave alone | 仅审查 |
| `install.md` | `install.jl` (5) | Leave alone | 仅审查 |
| `JuBat.md` | `JuBat.jl` (89) | Leave alone（联动 export 清理） | 死码删除后同步 export 表 |

---

## 跨文件约束

执行任何 plan 时必须遵守 spec 的不变量：

1. **活跃公共 API 不删减**（spec §0.3）
2. **High-risk-leave-alone 仅记录不动手**（spec §4.1）
3. **每步 commit + 测试**（spec §5）
4. **每文件完成后更新 `Simplify/baseline.md`**（spec §3.3）
5. **每个 PR 不跨层**（spec §12.2）

---

## 状态追踪

每个 plan 文件顶部应有状态标记。执行时更新：

- ⬜ Pending（未开始）
- 🔄 In Progress（执行中）
- ✅ Completed（已完成且测试通过）
- ⚠️ Blocked（被阻塞；说明原因）

---

**入口**：从上方"执行顺序"表阶段 0 开始（基线固化），按推荐顺序推进。
