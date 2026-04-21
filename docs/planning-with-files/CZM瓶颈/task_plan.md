# Task Plan: CZM 求解器性能瓶颈诊断与优化

## Goal
诊断并消除 CZM 求解器在实际仿真中的性能瓶颈与不收敛问题，确保主例程稳定完成。

## Current Phase
Phase 5 — 弧长法独立与回滚验证

## Phases

### Phase 1: 缓存优化实施
- [x] 实现 CZMAssemblyCache / CZMAssemblyWorkspace 结构体
- [x] K_bulk / 几何 / BC 缓存
- [x] mul! 替换 A*B、预分配稀疏矩阵
- [x] workspace 跨时间步复用
- **Status:** complete — 缓存复用验证通过

### Phase 2: 瓶颈定位
- [x] 微基准测试各操作耗时
- [x] 实际仿真路径计时（Newton 2000 iter/步，全部不收敛）
- [x] 诊断 `max_Δu` 限幅（表面原因）
- [x] 发现去掉限幅后产生 NaN（直接原因）
- [x] 确认 NaN 跨模块传播链
- **Status:** complete

### Phase 3: 根因修复
- [x] 修复 `Δu = K_bc \ R_bc` 重复计算
- [x] 实现回溯线搜索（替代 `max_Δu=1e-6` 硬编码）
- [x] 添加诊断日志（NaN 来源追踪）
- [x] 清理 debug println 和 profiling 代码
- [x] **运行仿真验证线搜索效果**
- [x] **根据日志定位 NaN 首次出现位置**
- [x] 验证修复后收敛性和仿真时间
- **Status:** complete
  - 备注：线搜索、自适应 load stepping、失败回滚和收敛后再提交已落地；主脚本验证通过。

### Phase 4: 清理与验证
- [x] 运行完整 testexample.jl 验证
- [x] 确认 CZM 收敛且仿真时间合理
- [x] 更新文档
- **Status:** complete
  - 备注：`example/testexample.jl` 最新复跑成功，未再出现 NaN 或 `CZM solve issue` 警告。

### Phase 5: 求解器接口拆分
- [x] 抽出 `basic` / `arc_length` 独立实现
- [x] 三种方法统一失败回滚与收敛后提交语义
- [x] 运行主例程和弧长法短探针验证
- **Status:** complete
  - 备注：`arc_length` 已不再 fallback 到 `load_substep`；主例程复跑成功，弧长法短探针可稳定完成。

## 根因链总结

```
表面原因: max_Δu = 1e-6 硬编码限幅
  → Newton 步长 0.55 被裁到 1e-6 (55万倍)
  → 位移几乎不动，残差永远不降，2000 iter 空转

直接原因: 固定子步 + 无回退的求解策略在真实耦合载荷下不稳定
  → Newton 第一步可能出现 NaN 或在有限残差下反复空转
  → 仅靠硬限幅会把问题掩盖成“几乎不动”

传播机制: 非收敛 damage_states / u_czm_prev 被回写到下游
  → 热模型读取异常损伤 → 热刚度矩阵异常 → 温度场与下次 CZM 输入被污染
  → 这就是 NaN 和不收敛被放大的原因

根本原因: 求解策略缺少“自适应载荷推进 + 失败回滚 + 收敛后再提交”
  → 通过线搜索、自适应子步和回写保护已经修复
```

## Key Questions
1. ~~为何微基准 0.024s 但实际仿真 14s？~~ → Newton 跑满 2000 iter
2. ~~缓存是否正确复用？~~ → 是，正常
3. ~~为何 Newton 不收敛？~~ → `max_Δu=1e-6` 裁剪步长 55 万倍
4. ~~去掉限幅后为何仍不收敛？~~ → NaN 传播，载荷量级与求解策略不匹配
5. ~~载荷归一化是否正确？~~ → 是，已验证
6. **NaN 首次出现在哪里？CZM 内部 vs 上游传入？** ← 已通过复跑确认，问题出在 CZM 非收敛状态回写后污染下游，不是缓存或归一化
7. **线搜索是否能阻止 NaN？** ← 已验证，需配合自适应子步和回滚一起使用

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| 先做微基准再做实际路径计时 | 隔离变量 |
| 诊断脚本 `example/czm_bottleneck_diagnosis.jl` | 可重复验证，排除其他模块干扰 |
| 回溯线搜索替代硬限幅 | 标准做法，自适应步长控制 |
| 不做 NaN 静默跳过 | 让错误暴露，通过日志定位根因 |
| 保留 `@warn` 诊断日志 | 帮助追踪 NaN 来源 |
| 非收敛状态不回写 | 防止 damage_states / u_czm_prev 污染后续步骤 |
| 自适应 load stepping + 失败回滚 | 让真实耦合载荷下也能稳定推进 |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| @printf 未定义 | 1 | 添加 using Printf |
| nnz 未定义 | 1 | 添加 using SparseArrays |
| include 路径错误 (/tmp) | 1 | 改用 @__DIR__ 相对路径 |
| Julia soft scope (`idx` undefined) | 1 | 用函数封装避免全局 scope 问题 |
| 顶层 docstring 导致 parse error | 1 | 改为注释 |
| `Δu = K_bc \ R_bc` 重复计算 | 1 | 删除重复行 |
| 移除限幅后 NaN | 1 | 实施回溯线搜索并加入自适应子步 |
| `solve_czm_step` ParseError | 1 | 恢复 `if/elseif/end` 结构，把自适应回退逻辑移回 `newton_raphson_czm` |

## Notes
- 诊断脚本确认：无 Δu 限幅时 Newton 2次/子步即收敛，任何限幅都阻止收敛
- 实际仿真已改为自适应 load stepping + 回滚，主脚本验证通过
- `apply_bc_czm` 使用罚函数法 penalty=1e12，当前未观察到新的不稳定来源
- 最新复跑：`example/testexample.jl` 18 步完成，wall-clock 7.887 s，CZM 累计 0.978 s，平均 54.324 ms/call，未再出现 NaN 或 `CZM solve issue` 警告
