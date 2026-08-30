# 层分辨应力求解：任务计划

## 目标

彻底删除粗网格抹平应力求解 `thermal_diffusion_stress_2D` 的旧实现（热网格单一材料 + 求和特征应变），替换为力学子网格上的层分辨应力，并以"耦合在线导出 + 固体按需工具函数"双流程落地：

1. **耦合流程（opt.czm_enabled=true）**：求解过程中于 Solve 主循环 CZM 更新块内收割本步收敛位移（`case.czm_layout.u_prev`）与当步载荷，经共享恢复核得 σ，写入 `diffusion stress xx/yy/xy/vonMises` 历史，PostProcessing 导出 `[Pa]` 结果键。
2. **固体流程（opt.czm_enabled=false）**：`thermal_diffusion_stress_2D(case, variables)` 保留原名同签名、原址实现替换，仅在显式调用时于 `czm_submesh.mesh_bonded` 上做层分辨线性静力求解。
3. 两流程共享恢复核 `recover_bulk_stress` 与输出键（Pa/m 有量纲），单一应力定义。

## 约束

- `example/testexample.jl` 冻结求解指标（Simplify/baseline）按记录精度不变；应力输出纯增量单向。
- 不改 `update_czm_damage!` / `solve_czm_step` / czm.jl 装配的数值行为。
- 不新增 Option 字段；新函数命名无 `_` 前缀、无 `!` 结尾。
- 新代码不新增防御性断言（唯一保留旧函数首行 E_coat 断言，原文不动）。
- 沿用跨界面 α_eff（:PE_PCC.α）；不新增 SP/PCC/NCC alphaT。
- 几何非线性/J2/预应力与线性恢复核不相容：分配与收割整体跳过 + 一条 @warn。

## 阶段

- Phase 1：规划文件与事实核实（尺度转换、执行顺序、历史写入机制）。
- Phase 2：src/Mechanical.jl 原址重写（recover_bulk_stress / export_macro_stress / thermal_diffusion_stress_2D）。
- Phase 3：管线接线（Variables.jl 分配、Solve.jl CZM 更新块内收割、PostProcessing 导出）。
- Phase 4：调用方迁移（testexample / jellyroll_stress_displacement / test_electrode_coat_modulus）。
- Phase 5：新测试 test/test_layer_resolved_stress.jl。
- Phase 6：验证门（基线指标、示例、全量测试、静态 grep）与文档同步。

## 验收门

1. testexample 完整运行 exit 0，冻结指标一致；应力图出现层间拉压交替、量级回落 ~90 倍至 O(1) MPa。
2. 两个示例脚本运行通过；新测试与全量 test 套件通过。
3. 旧键名（"thermal stress vonMises"、"diffusion stress vonMises only"）零残留。
4. AGENTS/CLAUDE/md 文档同步。
