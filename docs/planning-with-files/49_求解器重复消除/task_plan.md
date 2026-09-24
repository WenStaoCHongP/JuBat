# 49_求解器重复消除（R1）任务计划

> 立项 2026-09-11（用户批准 R1 计划，指示先做 R1-d）。vibe-coding 框架评估的 R1 建议：
> 消除 CzmSolve 四求解器的真重复与陈旧拷贝，不强行统一 λ 延续控制流。

## 子批
- **R1-d（已完成）**：load_substep geo 求解统一到 solve_equilibrated（消灭未均衡化 69× 填充慢路径；数值路径改变走容差门）
- R1-a：公共初始化段提取 czm_step_init（4 调用点 ×25 行）
- R1-b：残差+BC 三行收编 czm_residual_at（×7）
- R1-c：load_substep 内联线搜索改调 backtrack_line_search!（同算术，位级目标）
- 不做：arc 两变体 λ 延续统一（YAGNI；留待任务 30 Batch 6/7 触发）

## 行数测算：R1 全部完成 ≈ −110 行（1132→~1020）；R1-d 本身 +5 行（换来的性能/一致性远超行数代价）
