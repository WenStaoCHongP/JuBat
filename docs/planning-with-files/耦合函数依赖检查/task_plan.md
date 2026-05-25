# 耦合函数依赖检查 - 任务计划

## 目标
梳理 Jellyroll 锂电池 SPMe + 二维分布式热 + CZM 耦合仿真涉及的 src 函数，
明确调用关系、分支条件、启用开关和数据类型，为后续代码简化提供指导。

## 范围
- 仅关注 SPMe + distributed2D + CZM 耦合路径
- 跳过 SPM/P2D/Ring/lumped 等无关分支

## 产出
1. `findings.md` - 主文档（Part A: 执行流程总览 + Part B: 模块函数清单）
2. `progress.md` - 进度跟踪

## 关键开关
| 开关 | 耦合路径值 | 含义 |
|------|-----------|------|
| `opt.model` | "SPMe" | 电化学模型 |
| `opt.thermal_enabled` | true | 启用热耦合 |
| `opt.thermalmodel` | "distributed2D" | 二维分布式热 |
| `opt.per_element_spme` | true | 逐单元SPMe |
| `opt.czm_enabled` | true | 启用CZM损伤 |
| `opt.mechanicalmodel` | "full" | 完整力学耦合 |
| `opt.cool_method` | "tab"/"surface" | 冷却方式 |
| `opt.czm_iter_method` | "basic"/"load_substep"/"arc_length" | CZM迭代方法 |
