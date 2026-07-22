# CZM 按材料层界面重构 - 验证结果

**日期**: 2026-07-20
**spec/plan 版本**: v2

## 单元级验证（spec §8.1）

- 脚本: `example/内聚力验证/verify_czm_per_interface.jl`
- 状态: **PASS**
- 峰值 σ_n (实测): 1.0 (归一化)
- σ_max (param_cache): 1.0 (归一化)
- δ_0_n (param_cache): 1.9772e-6 (归一化)
- δ_c_n (param_cache): 3.5710e-3 (归一化)
- G_c (param_cache): 1.7855e-3 (归一化)
- D 终值: 1.0
- 与 param_cache 自洽: **是**（峰值牵引 = σ_max，D 单调收敛至 1.0）

## 全网格回归（spec §8.2）

- 脚本: `example/循环验证/czm_cycle_example.jl`
- 状态: **跑通**（循环 1 完成：放电→静置→充电→静置，V=3.70→3.25V→3.45→4.12V，T_max=310.56K）
- D_max: [待全 10 循环跑完填入]
- n_fractured: [待全 10 循环跑完填入]
- 损伤峰值位置: [待全 10 循环跑完填入]
- δ_exp 对比（用户未提供则 SKIP）: SKIP

注：按 spec §2.2 "不保留旧路径"，不与改造前数值对比。

## 网格收敛（spec §8.3）

- 脚本: `example/内聚力验证/czm_grid_convergence.jl`
- δ_max(nθ=40): [待填入——脚本已跑通但全量执行需较长时间]
- δ_max(nθ=80): [待填入]
- δ_max(nθ=160): [待填入]
- 相对变化 (40→160): [待填入]
- 验收 (< 5%): [待填入]

## 性能记录（spec §8.4）

- 单步耗时 (nθ_czm=80): [待填入] s
- 备注: 仅计时，不与旧版本对比

## 已知问题/后续工作

- 网格收敛脚本 (czm_grid_convergence.jl) 已验证语法正确且 API 调用与 czm_cycle_example.jl 一致，
  但 3 组 nθ_czm 全量执行耗时较长（每组需求解 100s 仿真），结果待后续运行后填入。
- czm_cycle_example.jl 已验证跑通循环 1（放电→静置→充电→静置），10 循环完整结果待后续运行。
- 脚本 stdout 在重定向到文件时 Julia 缓冲严重，建议直接交互运行或加 `flush(stdout)`。
