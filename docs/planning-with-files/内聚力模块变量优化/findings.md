# Findings & Decisions

## Requirements
- 不新增 CZMState，扩展现有结构
- 计算结果和中间变量命名在 Variables.jl 定义
- PostProcessing.jl 负责还原物理单位
- 需要兼顾现有耦合求解、损伤演化和后处理输出

## 当前问题：变量键定义分散且不一致

### 键名不一致对照表

| 语义 | StandardVariables (预分配) | czm_output_to_variables (写入) | PostProcessing.jl (还原) |
|------|---------------------------|-------------------------------|-------------------------|
| 最大损伤 | `czm D_max` | `czm max damage` | `czm D_max` |
| 平均损伤 | `czm D_mean` | `czm mean damage` | `czm D_mean` |
| 断裂数 | `czm n_fractured` | `czm fractured elements` | `czm n_fractured` |
| 法向最大分离 | `czm δ_max_n` | — | `czm δ_max_n` |
| 法向平均分离 | `czm δ_mean_n` | — | `czm δ_mean_n` |
| 位移 x | — (未预分配) | `czm displacement x` | `czm displacement x [m]` |
| 位移 y | — (未预分配) | `czm displacement y` | `czm displacement y [m]` |
| 损伤场 | — (未预分配) | `czm damage` | `czm damage [0-1]` |
| 法向牵引 | — (未预分配) | `czm traction normal` | `czm traction normal [Pa]` |
| 切向牵引 | — (未预分配) | `czm traction tangent` | `czm traction tangent [Pa]` |
| 法向分离 | — (未预分配) | `czm separation normal` | `czm separation normal [m]` |
| 切向分离 | — (未预分配) | `czm separation tangent` | `czm separation tangent [m]` |
| 负极损伤 | `negative electrode cohesive zone damage` | — | — |
| 正极损伤 | `positive electrode cohesive zone damage` | — | — |

### 问题总结
1. **键名不一致**: D_max vs max damage, n_fractured vs fractured elements
2. **预分配缺失**: czm_output_to_variables 写入的 6 个场变量（displacement, damage, traction, separation）在 StandardVariables 中未预分配
3. **create_element_workspace 不完整**: 只有电极损伤键，缺少场变量键
4. **变量尺寸不匹配风险**: 未预分配的键通过 `copy(variables)` 动态添加，Variable_update! 可能遗漏

### PostProcessing.jl 还原公式（不可破坏）

| 变量 | 还原公式 | 尺度字段 |
|------|---------|---------|
| displacement x/y | × `case.param_dim.scale.L` | 长度尺度 |
| traction normal | × `case.param.scale.E_n` | 正极弹性模量 |
| traction tangent | × `case.param.scale.E_p` | 负极弹性模量 |
| separation normal/tangent | × `case.param.scale.r0` | 参考长度 |
| damage | 无量纲，直接传递 | — |
| D_max, D_mean | 无量纲，直接传递 | — |
| δ_max_n, δ_mean_n | × `case.param_dim.scale.L` | 长度尺度 |
| n_fractured | 整数→Float64，直接传递 | — |

## 现有状态结构（不新增，只扩展）

### Case 中 CZM 相关字段 (SetCase.jl:97-107)
```julia
mutable struct Case
    ...
    czm_mesh::Union{Nothing, CohesiveMesh}      # 网格+损伤状态
    czm_cache::Union{Nothing, CZMAssemblyCache}  # 装配缓存
    # czm_layout::Union{Nothing, CzmLayout}      ← 新增
end
```

### CouplingState.jl 中的 CZM 类型
| 类型 | 行号 | 职责 |
|------|------|------|
| `CohesiveElementGeom` | 106-116 | 预计算单元几何 |
| `CZMAssemblyWorkspace` | 124-161 | 步内装配工作区 |
| `CZMAssemblyCache` | 170-186 | 长期装配缓存 |

### update_czm_damage! 当前签名 (CouplingState.jl:291)
```julia
function update_czm_damage!(czm_mesh, czm_params, case, variables, T_nodes_carry, u_czm_prev)
```
→ 目标: `update_czm_damage!(case, variables)` — 从 case.czm_layout 取 u_prev

## 建议的统一键名规范

### StandardVariables 中 CZM 块（替换现有 lines 136-144）
```julia
if case.opt.czm_enabled == true
    n_coh = case.czm_mesh.n_cohesive
    n_nodes = size(case.czm_mesh.nodes, 1)

    # ── 标量统计 ──
    variables["czm D_max"]          = zeros(Float64, 1, num)
    variables["czm D_mean"]         = zeros(Float64, 1, num)
    variables["czm n_fractured"]    = zeros(Float64, 1, num)
    variables["czm δ_max_n"]        = zeros(Float64, 1, num)
    variables["czm δ_mean_n"]       = zeros(Float64, 1, num)

    # ── 场变量 (per cohesive element / per node) ──
    variables["czm damage"]         = zeros(Float64, n_coh, num)
    variables["czm displacement x"] = zeros(Float64, n_nodes, num)
    variables["czm displacement y"] = zeros(Float64, n_nodes, num)
    variables["czm traction normal"]  = zeros(Float64, n_coh, num)
    variables["czm traction tangent"] = zeros(Float64, n_coh, num)
    variables["czm separation normal"]  = zeros(Float64, n_coh, num)
    variables["czm separation tangent"] = zeros(Float64, n_coh, num)

    # ── 电极级损伤 (per thermal element) ──
    variables["negative electrode cohesive zone damage"] = zeros(Float64, Nn, num)
    variables["positive electrode cohesive zone damage"] = zeros(Float64, Np, num)
end
```

### czm_output_to_variables 统一键名（去掉 czm max damage 等旧名）
```julia
new_variables["czm D_max"]       = stats.max_D
new_variables["czm D_mean"]      = stats.mean_D
new_variables["czm n_fractured"] = Float64(stats.n_fractured)
new_variables["czm damage"]      = result.damage
new_variables["czm displacement x"] = u_x
# ... 同理
```

### CzmLayout (CouplingState.jl 新增)
```julia
struct CzmLayout
    n_coh::Int
    ndof::Int
    u_prev::Vector{Float64}
end
```

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| 不新增 CZMState | 现有 Case.czm_mesh + czm_cache 覆盖大部分需求 |
| CzmLayout 只放 n_coh, ndof, u_prev | 最小扩展，收敛签名 |
| 变量键统一为英文无空格下划线风格 | 与现有 variables Dict 中电化学键风格一致 |
| 在 StandardVariables 中预分配所有 CZM 键 | 避免动态添加导致的 Variable_update! 遗漏 |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| 初次读取 memory 文件时行范围写得过大 | 重新读取正确范围并继续整理 |

## Resources
- src/Variables.jl (变量定义)
- src/CzmPostProcess.jl (czm_output_to_variables)
- src/PostProcessing.jl (物理单位还原)
- src/CouplingState.jl (CZM 类型定义)
- src/SetCase.jl (Case 结构体)
- src/Solve.jl (调用入口)

## Visual/Browser Findings
- 无

---
*Update this file after every 2 view/browser/search operations*
