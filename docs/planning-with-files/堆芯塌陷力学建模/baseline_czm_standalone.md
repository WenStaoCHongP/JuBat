# CZM 基线快照（verify_czm_standalone，Batch 1 门禁基准，方案 B）

- **入口**: `tools/verify_czm_standalone.jl`
- **命令**: `GKSwstype=100 JULIA_NUM_THREADS=1 julia --startup-file=no --project=. tools/verify_czm_standalone.jl`
- **环境**: Julia 1.11.2，1 thread，`--startup-file=no`，`GKSwstype=100`
- **Git HEAD**: `e117fd2f4d48ebb80eefd2616cc19370d16cb65b`
- **冻结日期**: 2026-08-21（计划 v1.1 修订时执行，用户决策方案 B）
- **背景**: 原 Batch 1 门禁工具 `tools/czm_baseline_probe.jl` 已被 `2bf2ac7` 删除，经评审（2026-08-21）由本工具替代；spec v1.3 已同步替换 §5/§7 引用。工具在 `2bf2ac7` 中修复 `build_czm_cache` 签名后实测通过，本快照冻结其当前行为。
- **求解器参数**: tol=1e-4，max_iter=200，n_load_steps=50，visc_beta=1.0；F_ext=0，dT_elem=0，Δsoc_p=0（仅 SOC 驱动）。

## 网格与有效参数（实测）

| 量 | 值 |
|---|---|
| 网格 | nθ=40，gsorder=2，czm_enabled=true |
| Nodes | 10946 |
| Bulk elements | 6728 |
| Cohesive elements | 3364 |
| E_eff | 6.097561e+00 |
| ν_eff | 0.3000 |
| α_eff | 4.470000e-03 |
| β_n | 3.423743e-02 |
| β_p | -1.531324e-02 |

## 收敛对比表（冻结）

```
Δsoc_n | basic                  | load_substep           | arc_length
 0.100 | OK   2it D=0.0000 r=8.0e-09 | OK 100it D=0.0000 r=5.8e-09 | OK  50it D=0.0000 r=4.6e-08
 0.500 | OK   2it D=0.0000 r=4.1e-08 | OK 100it D=0.0000 r=2.8e-08 | OK  50it D=0.0000 r=2.4e-07
 1.000 | OK   2it D=0.0000 r=8.1e-08 | OK 100it D=0.0000 r=5.7e-08 | OK  50it D=0.0000 r=4.7e-07
 1.500 | OK   2it D=0.0000 r=1.3e-07 | OK 100it D=0.0000 r=8.7e-08 | OK  50it D=0.0000 r=6.9e-07
 2.000 | OK   2it D=0.0000 r=1.6e-07 | OK 100it D=0.0000 r=1.1e-07 | OK  50it D=0.0000 r=9.4e-07
 3.000 | OK   2it D=0.0000 r=2.6e-07 | OK 100it D=0.0000 r=1.7e-07 | OK  50it D=0.0000 r=1.4e-06
 5.000 | OK   2it D=0.0000 r=4.1e-07 | OK 100it D=0.0000 r=2.9e-07 | OK  50it D=0.0000 r=2.3e-06
10.000 | FAIL 2it D=0.0000 r=1.6e+03 | FAIL 94it D=0.0000 r=7.7e-04 | FAIL 1837it D=0.0000 r=8.0e-04
```

（`10.000` 行伴随两条 stall 警告：`CzmSolve.jl:603` adaptive load stepping 与 `CzmSolve.jl:430` arc-length stepping，load_progress≈0.639、step_size 至下限——该警告属 10.0 水平的冻结行为。）

## Summary（冻结）

```
basic          : converged 7/8 levels, D_max=0.0000 (final) / 0.0000 (peak), total_iter=16
load_substep   : converged 7/8 levels, D_max=0.0000 (final) / 0.0000 (peak), total_iter=794
arc_length     : converged 7/8 levels, D_max=0.0000 (final) / 0.0000 (peak), total_iter=2187
```

## 冻结说明

1. **`FAIL` 条目是冻结行为**。三方法在 Δsoc_n=10.0 均不收敛（与 `2bf2ac7` 提交信息"三方法 7/8 载荷水平收敛"一致）；不得通过调参（tol/max_iter/n_load_steps）使其收敛后当作通过。
2. 全表 D=0.0000：本快照只行使线弹性 bulk + cohesive 装配与求解路径，损伤通路由 `unit_czm_bilinear.jl`/`unit_czm_eigenstrain.jl` 等单测覆盖。这恰好匹配 Batch 1 的改动面（bulk 装配接线）。
3. 工具既有瑕疵（已登记 findings，不修改、不影响门禁）：`:66` 模板网格传 `case.param` 而 `:134` 逐方法重建传 `param_dim`；两处均用未合并 `mesh_data.thermal2D`。CZM 求解只消费 bulk/cohesive 拓扑，不受影响；若未来 `create_czm_mesh` 开始实质消费 `param` 尺度，须先统一并**重冻结本快照**（同批声明）。

## 比较规则

- Batch 1 起每批完成后重跑本工具，上两节全部数值必须在打印精度下逐位一致（含 OK↔FAIL 状态与迭代数）。
- 不一致即视为行为漂移：停止该批，先定位或回退，不得以"数值接近"或"FAIL 仍 FAIL"以外任何宽松判据放行。
- 重跑环境必须与冻结时一致（Julia 1.11.2、单线程、`GKSwstype=100`、`--startup-file=no`）。
