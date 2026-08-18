# 仿真提速-SPMe热CZM task_plan

来源：docs/superpowers/plans/2026-08-19-perf-optimization-spme-thermal-czm.md

- [x] Task 1: 批次 0a — 本机基线锚点与 PNG SHA 稳定性验证
- [x] Task 2: 批次 0b — 矩阵常量性探针 + 等价性实验 + CZM 残留核对 + 批次裁决
- [ ] Task 3: [关] 热矩阵不变量缓存（热块常量成立但占比 0.77%，投入产出不成立）
- [ ] Task 4: [关] 全局 M/K 缓存（K 化学块随状态变）
- [ ] Task 5: [关] 单元循环分配消除（SPMe 2.66% < 10%）
- [ ] Task 6: [关] 按 dt 分解缓存（矩阵常量性 false）
- [ ] Task 7: [开] CZM 装配优化（CZM 95.80%、残留核对完成、等价性 2 通过）
- [ ] Task 8: 收尾汇总与文档同步
