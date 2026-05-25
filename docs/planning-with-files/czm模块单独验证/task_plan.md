# Task Plan: CZM模块单独验证

## Goal
在不接入电化学-热耦合的前提下，完成 CZM 内聚力子模块的独立验证方案，并将缺失的扩散应力和热应力改为函数输入，形成可复用的 standalone 验证入口，同时支持时间历史下的 D(t) 追踪。

## Current Phase
Complete

## Phases

### Phase 1: Scope and interface freeze
- [x] 确认验证目标：只验证内聚力子模块，不联动 SPMe/热求解器
- [x] 确认缺失载荷由函数输入提供，而不是回接耦合链
- [x] 汇总现有可复用脚本和入口
- **Status:** complete

### Phase 2: Benchmark design
- [x] 定义纯本构、单单元、单界面、混合模式四类基准
- [x] 为每类基准定义输入、输出和判据
- [x] 明确 Mode I only 与 mix/BK 两条验证路径
- **Status:** complete

### Phase 3: Standalone driver
- [x] 设计独立验证脚本或 wrapper
- [x] 让扩散应力和热应力以函数参数注入
- [x] 保持与现有 CZM API 最小耦合
- **Status:** complete

### Phase 4: Numerical verification
- [x] 与解析解和基准脚本对比
- [x] 记录收敛、损伤、牵引-分离曲线
- [x] 检查不同 loading path 的稳定性
- **Status:** complete

### Phase 5: Documentation and handoff
- [x] 写清输入接口、使用方法、测试边界
- [x] 更新相关计划和验证说明
- [x] 输出最终结论和剩余风险
- **Status:** complete

## Key Questions
1. standalone validation 应该用应力向量输入，还是用 provider callback? -> 采用 provider callback，返回 `dT_elem / Δsoc_n_elem / Δsoc_p_elem` 数组，便于构造不同工况。
2. 哪些 cases 是必须通过的最小验收集? -> 纯 Mode I 本构、输入驱动的单步求解（uniform + gradient）、mix/BK 基线。
3. 新入口应放在 `tools/` 还是 `example/内聚力验证/`? -> 以 `tools/verify_czm_standalone.jl` 作为 canonical CLI 入口，后续如需演示再加 example wrapper。

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| Keep electrochemical-thermal coupling out of the standalone path | isolate CZM behavior |
| Inject diffusion and thermal stress via function input | satisfy standalone requirement and avoid hidden dependencies |
| Split Mode I and mix/BK validation | align with current constitutive model behavior and keep expected results clear |
| Reuse existing analytical and unit scripts | minimize implementation cost |
| Use a provider callback that returns input arrays | keeps the standalone driver flexible without touching core solvers |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| None yet | 1 | - |

## Notes
- 这份计划只定义 CZM 独立验证边界，不修改电化学-热耦合主流程。
- 后续实现优先复用 `tools/verify_czm_unit.jl`、`tools/verify_czm_analytical.jl` 和 `example/内聚力验证/czm_example.jl` 里的验证思路。
