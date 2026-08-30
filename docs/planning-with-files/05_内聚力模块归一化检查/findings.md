# Findings & Decisions

## Requirements
- 检查内聚力模块的参数归一化是否一致
- 重点关注强度、位移、刚度及其对应尺度字段
- 需要能判断是否存在重复缩放、漏缩放或文档与代码不一致

## Research Findings
- 归一化入口位于 src/SetParams.jl 的 NormaliseParam(param_dim::Params)
- 原始 cohesive 参数定义位于 src/parameters/Jellyroll.jl，包括 σ_max_n、K_n、δ_c_n 等
- CZM 相关材料矩阵在 src/Materialmatrix.jl 中直接读取 cohesive_params.K_n、δ_c_n 等字段
- 文档 md/01_参数定义与归一化.md 已列出 cohesive.σ_max_n、cohesive.δ_c_n、cohesive.K_n 的物理量纲说明，需与代码保持一致
- 目前的检查重点是确认 NormaliseParam 是否对 cohesive 字段做了统一缩放，以及后续使用是否重复缩放
- 目前看到的实现是：NormaliseParam 将 cohesive 的法向/切向强度、起始位移、临界位移、断裂能和刚度都按 param_dim.scale.σ_czm、δ_czm、G_czm、K_czm 归一化
- Materialmatrix.jl 中的 bilinear_traction_state / bilinear_tangent 直接消费 cohesive_params.K_n、δ_0_n、δ_c_n，说明这里期望输入已经是归一化后的量
- 目前未发现 CZM 求解路径里再次对 cohesive 参数做二次缩放的迹象；下一步要核对 scale.σ_czm、scale.δ_czm、scale.K_czm 的定义来源
- scale.σ_czm = cohesive.σ_max_n，scale.δ_czm = cohesive.δ_c_n，scale.G_czm = scale.σ_czm × scale.δ_czm，scale.K_czm = scale.σ_czm / scale.δ_czm
- 这意味着归一化后的 cohesive.σ_max_n、cohesive.δ_c_n、cohesive.K_n 期望是 1 量级的无量纲数，而不是 SI 单位
- 需要确认其他代码路径没有把已归一化的 cohesive 值再次按物理单位处理
- 运行时检查 `ChooseCell("LG M50")` 返回的 scale.σ_czm、scale.δ_czm、cohesive.σ_max_n、cohesive.δ_c_n、cohesive.K_n 都是 0，导致 NormaliseParam 之后 cohesive.σ_max_n、cohesive.δ_c_n、cohesive.K_n 变成 NaN
- `ChooseCell("LG M50")` 只 include src/parameters/LGM50.jl，而该文件目前仅创建 cohesive = Cohesive()，没有设置 cohesive 数值；因此它不是用于 CZM 归一化验证的有效 preset
- `ChooseCell("Jellyroll")` 才会 include src/parameters/Jellyroll.jl，其中明确设置了 cohesive.σ_max_n、K_n、δ_c_n 等值，适合作为归一化验证样本
- Jellyroll 的切向参数满足 δ_0_t > δ_c_t（原始值和归一化值都如此），所以当前 bilinear traction law 在 shear 方向会退化为“损伤起点即全损伤”的跳变行为，而不是标准的软化区间
- 这属于参数一致性/模型设定问题，不是归一化公式本身的问题，但它会影响对 CZM 行为的物理解释
- CZM 后处理结果还原的实际实现是：位移 x/y 乘 param.scale.L，法向/切向牵引力分别乘 param.scale.E_n 和 param.scale.E_p，法向/切向分离乘 param.scale.r0
- 这与 md/01_参数定义与归一化.md 中“牵引力按 σ_czm 还原”的描述不完全一致，当前应以代码实现为准
- 如果要还原“节点合力”或残差向量，仓库里没有单独的力尺度字段，不能直接套用 σ_czm；必须回到积分/装配定义推回物理量纲

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| 以归一化入口和调用链为核心检查对象 | 能快速覆盖参数定义、传递和使用三层问题 |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| Julia 脚本中直接访问 Unicode 字段名时触发 type Scale has no field ? | 下一次改用 getfield(Symbol(...)) 读取 σ_czm / δ_czm / K_czm |

## Resources
- md/01_参数定义与归一化.md
- src/Option.jl
- src/czm.jl
- src/CzmSolve.jl
- src/mechanical.jl
- src/SetParams.jl
- src/parameters/Jellyroll.jl
- src/Materialmatrix.jl
- src/PostProcessing.jl

## Visual/Browser Findings
- None yet

---
*Update this file after every 2 view/browser/search operations*