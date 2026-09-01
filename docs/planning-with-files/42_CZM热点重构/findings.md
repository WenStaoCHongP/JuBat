# CZM 热点重构：发现与决策

## 需求

- 用户要求“使用 superpowers 做 CZM 热点重构任务”。
- Superpowers 的 14 个技能已从官方 `obra/superpowers` 安装，本轮已调用 `using-superpowers` 与 `brainstorming`。
- 仓库已有 `docs/superpowers/specs/`、`docs/superpowers/plans/` 和 `docs/superpowers/findings/`，可作为现有设计/执行证据。
- 复杂任务按 AGENTS.md 使用 `docs/planning-with-files/<中文任务名>/` 三件套，并同步总索引。
- 2026-09-01 用户明确批准实施方案 A，并确认当前内聚力参数是正确实验值；本任务不得再把参数标为占位值，也不得修改这些参数。

## 初步发现

- 内存记录指出既有性能计划明确要求先对照 `2026-04-20-czm-vectorized-solver-plan.md` 与 `docs/planning-with-files/09_向量化CZM/`，只做残留项。
- 记录中的历史残留热点包括：`K_bulk * u`、`f_int_bulk + f_int_coh`、`K_bulk + K_coh` 的分配，cohesive 稀疏矩阵标量累加，以及线搜索中的全量组装；这些必须以当前源码复核，不能直接视为现状。
- 按文件名搜索“CZM 热点/热点重构/hotspot”没有直接命中的任务文件；需要按正文、近期修改与计划勾选状态继续定位。
- 项目已有 CZM 拓扑 API、机械结构重构、性能优化等多轮历史变更，必须先排除已完成或被拒绝的工作。
- 内存索引只保留了旧快照中的 `K_bulk * u` 与 `K_bulk + K_coh` 证据；它不证明当前代码仍有该热点，后续必须逐行复核当前装配实现与近期提交。

## 2026-08-31 当前基线审计

- 当前分支 `codex/src-physics-modularization` 与远端同在 `39387f3`；最新相关提交 `98adb9a feat(tools): CZM 热点定位探针（行级采样聚合，复核四层重构后热点结论）` 已进入 HEAD。
- tracked 工作区没有改动；未跟踪项包括 `.mimosa/`、本任务三件套及若干 `output/` 结果。后续不得把这些无关输出纳入源码改动。
- `docs/planning-with-files/29_仿真提速-SPMe热CZM/` 记录：CZM 曾占总耗时 95.80%；Task 7 的四项装配优化满足四判据且位级一致，但只把 CZM 从 201.4 s 降到 198.0 s（-1.7%），源码随后按用户决定于 `5812c79` 回退。
- 同一历史任务把根因再定位为 Newton 迭代中的稀疏 LU，约占 CZM 时间 84%；候选方案是基于 `K_bc` 内容一致性的因子复用，而不是仅按对象或 pattern 复用。
- 历史机制假设：当前场景 `D=0`、分离近零、bilinear tangent 固定，使 `K_total/K_bc` 多轮数值可能不变，但每轮仍重新分解。该假设必须由最新 `98adb9a` 探针和当前源码验证。
- `98adb9a` 只新增 `tools/probe_czm_hotspot.jl`，不改变求解路径；提交说明明确以 `testexample` 同款 `nθ=80`、`mix`、逐单元 SPMe 完整求解进行行级采样，目标正是复核四层结构重构后 `K_bc \\ R_bc` 是否仍占主导。
- 当前 `src/CzmSolve.jl` 不止一个反斜杠求解点：basic、load-substep、arc-length、几何非线性/增广路径均有直接 `K_bc \\ rhs`。不能把单个缓存变量机械塞进一个行号；应先确认探针实际命中的 solver 分支，并决定复用边界是否覆盖一个或多个 RHS。
- `apply_bc_czm` 位于独立 `src/CzmBC.jl`；历史候选还提到它复制/重建大矩阵。因此应把“BC 施加后的矩阵内容稳定性”和“因子复用”视作同一调用链审计，而不是只替换 `\\` 语法。
- 最新探针固定 `opt.czm.iter_method="basic"`、`model="mix"`、`fix_inner=false`、`nθ=80`，所以它只能直接裁决 basic 路径；不能据此声称 arc-length 或几何非线性路径也获得相同收益。
- basic 路径每次 `solve_czm_basic_step` 创建/取得装配工作区，Newton 循环中先装配 `K_total`，再由 `apply_bc_czm` 复制 `K_total`/`R`、计算相对罚值并修改对角，随后执行 `K_bc \\ R_bc`。
- `apply_bc_czm` 的罚值依赖 `maximum(abs, diag(K))`；因此因子复用比较必须针对施加 BC 后的矩阵内容（至少结构、`nzval` 与尺寸），不能只比较未施加 BC 的 `K_total`。
- `geo_nl=true` 时 bulk tangent 明确依赖 `u` 且禁止 bulk 缓存；内容判据若严格逐项相等仍能安全拒绝变化矩阵，但本轮性能收益证据仅来自 `geo_nl=false` 的 basic 配置。
- 直接把因子缓存挂到演化状态会混淆状态所有权。后续需检查 `CZMAssemblyWorkspace`/mesh 缓存的现有生命周期，优先把纯性能缓存放在装配缓存层，并保证求解失败不会提交任何物理状态。
- 当前 `CZMAssemblyWorkspace` 只含单元级矩阵/向量、cohesive 全局内力、分离/牵引和 `K_coh` 稀疏缓冲；没有 BC 矩阵或 LU 因子字段。
- `CouplingState.jl` 的明确契约是：缓存随 `CohesiveMesh`，演化状态随 `MechState`。因此因子复用属于 mesh/workspace 缓存，不能并入 `MechState` 的 `u_prev`、damage/plastic/prestress 等提交语义。
- 四层重构后的 `assemble_coupled_system` 仍执行 `K_total = K_bulk + K_coh` 与 `f_int_total = f_int_bulk + f_int_coh`；这是旧 Task 7 已验证可位级优化但收益仅 -1.7% 的次要热点，除非新 profile 推翻旧结论，本轮不扩张到该回退范围。
- `solve_czm_step` 是统一分派入口：basic、load_substep、普通 arc-length、geo nonlinear arc-length 分支不同。最新 probe 只运行 basic，因此第一批实现应保持其他分支原样，避免无性能证据地扩大矩阵复用面。
- `assembly_workspace(czm_mesh)` 会把同一 `CZMAssemblyWorkspace` 惰性挂载到 `czm_mesh.ws` 并跨时间步复用；LU 内容快照与因子缓存放入该 workspace 可自然覆盖 Newton 内与相邻 CZM 更新，同时仍属于纯缓存。
- `bulk_stiffness`、`cohesive_geometry`、`assembly_workspace` 都以 `nothing` 为唯一惰性初始化门，依赖 `SetCase` 后参数冻结和 mesh 对象身份稳定；新的线性求解缓存应遵循同一失效模式，不另造参数回退或静默重建语义。
- 尚需确认 `CohesiveMesh` 构造字段、Julia 1.11 UMFPACK 因子类型和 `ldiv!`/`F \\ b` 位级等价，之后才能给出可实现设计。
- `CohesiveMesh.ws::Any` 已是 mesh 级惰性缓存槽，具体类型由 `assembly_workspace` 断言；扩充 `CZMAssemblyWorkspace` 不需要给 `MechState` 或公开 `CohesiveMesh` API 增加新缓存字段。
- Julia 1.11.2 小实验得到显式 `lu(A)` 类型 `SparseArrays.UMFPACK.UmfpackLU{Float64,Int64}`，但 `A \\ b == lu(A) \\ b` 和 `A \\ b == ldiv!(...,lu(A),b)` 均为 `false`（末位不同）。这直接否定“任意稀疏矩阵上显式 LU 复用必定位级等价”的历史泛化。
- 可能原因是 `A \\ b` 对具体矩阵选择的分解/排序路径与显式 `lu(A)` 不同；必须进一步确认生产 `K_bc` 的实际 dispatch，并把生产矩阵上的逐位等价作为开批硬门。若无法逐位一致，则不能以近似相等替代项目基线。
- 后续 Julia 1.11.2 实验确认：`factorize(A)` 对对称稀疏矩阵返回 `SparseArrays.CHOLMOD.Factor`、对非对称稀疏矩阵返回 `SparseArrays.UMFPACK.UmfpackLU`；两类均满足 `A \\ b == factorize(A) \\ b` 逐位相等。正确候选是缓存 `factorize(K_bc)` 的结果，而不是强制 `lu(K_bc)`。
- `CohesiveMesh.ws::Any` 已允许 workspace 内持有具体分解对象，无需扩大 `CohesiveMesh` 公开字段；但缓存命中仍须同时校验尺寸、CSC 结构（`colptr/rowval`）和 `nzval` 内容，不能只比 `nzval` 后假定 pattern 永远不变。
- 未跟踪的 `output/ab_probe.jl` 是四层重构前后 API A/B 草稿，仍调用旧 `solve_czm_basic_step(..., u0)[1]` 形态；当前 API 已改为 `MechState`，本轮不修改或执行该用户文件。
- 现有 `test/unit_czm_newton.jl` 是独立单元条带参考求解器，不直接覆盖生产 `solve_czm_basic_step` 的 workspace/缓存生命周期；新回归应针对生产 basic 入口和 mesh-cached workspace。
- 可复用现有测试契约：`test_czm_mech_core.jl` 已用 `===` 验证 `K_bulk` 对象缓存，并用 `==` 要求装配逐位一致。新测试应同样区分缓存命中（因子对象身份/计数）与科学输出逐位相等，并覆盖矩阵数值改变后的失效。

## 当前热点探针结果

- 2026-08-31 在 Julia 1.11.2、单线程、当前 HEAD `39387f3` 运行 `tools/probe_czm_hotspot.jl` 成功，完整 60 s 仿真退出码 0。
- 顶层 `probe_czm_hotspot.jl:66` 与 `Solve.jl:1` 各有 8,771 个采样，代表约 8,771 个 profile 样本栈；`CzmSolve.jl:234` 的 `K_bc \\ R_bc` 出现 6,831 次，按样本栈计约为 **77.88%**，确认四层重构后线性求解仍是绝对主热点。
- 探针当前把所有 JuBat 栈帧计数相加得到 69,464 并作为百分比分母，因此打印 `CzmSolve.jl:234 = 9.83%` 是“栈帧占比”，不是墙钟样本占比。若用于验收，应同时报告样本栈分母或改用 flat/exclusive 口径，避免低估热点。
- 当前 profile 与历史 84% 方向一致，但数值口径不同；可据此开启“basic 路径线性求解复用”设计，不能据此外推其他 solver 分支。

## 生产 `K_bc` 因子复用可行性

- 用当前 probe 同款 `nθ=80` mesh 构造真实 basic 线弹性 `K_bc`：尺寸 `40552×40552`、`nnz=620080`、`issymmetric=false`，`factorize` 选择 `SparseArrays.UMFPACK.UmfpackLU{Float64,Int64}`。
- 同一真实矩阵/RHS 上，`K_bc \\ R_bc` 与 `factorize(K_bc) \\ R_bc` 满足 `==`，`max_abs_diff=0.0`，通过生产矩阵逐位等价门。
- 单次非隔离 timing：直接 `\\` 约 3.36 s，显式 `factorize` 约 2.39 s，缓存因子回代约 0.024 s。绝对数受预热/机器影响，但约两个数量级的分解-回代差距足以证明热点值得开批。
- 因为当前真实 `K_bc` 非对称，缓存对象在该路径是 UMFPACK；实现仍应调用 `factorize` 而不写死该类型选择，以保留未来对称性/算法 dispatch。
- 尚缺“连续 Newton/时间步矩阵内容命中率”实测；严格内容判据保证未命中时退回原语义，因此命中率只影响收益、不影响正确性。可把命中/未命中计数仅放入测试辅助或探针，不写入科学结果。

## 逐 cohesive 元素数组路径

- 当前有两套分离位移入口：`src/Tools.jl` 的公开单元素查询 `compute_separation(czm_mesh, elem, u)`，以及生产装配 `src/czm.jl::assemble_czm_system` 内的批量逐元素实现。
- 生产路径遍历 `1:n_coh`，从预计算 `geom_cache[i].dofs/R/...` 取几何，手工把 8 个全局位移分量装入 `ws.u_e`，用 `mul!` 得到 `ws.δ_local`，再计算高斯点分离、牵引与切线，并写入预分配的 `ws.separations/ws.tractions/ws.f_int_coh/K_coh`。
- 因此“逐元素数组计算”确实存在，但旧向量化工作已经覆盖了 workspace、几何缓存和 `mul!`；当前 profile 中整个 `czm.jl` 仅占全部 JuBat 栈帧约 1.27%，远低于线性分解。是否继续优化应先做行级/分段微基准，不能仅凭循环外观判断热点。
- 公开 `compute_separation`（`src/Tools.jl:77-109`）每个高斯点都会分配 `B_global=zeros(2,8)`、`u_e=[...]`，并执行 `R * B_global * u_e`，存在明显小数组临时分配；仓库当前内部调用搜索只命中测试，生产 Newton 装配不走该函数。
- 生产 `assemble_czm_system`（`src/czm.jl:178-288`）把 `u_e/B_global/B_local/δ_local` 放入 workspace，使用 `fill!` 和 `mul!`；分离结果写入预分配 `ws.separations[i]`，牵引写入 `ws.tractions[i]`。
- 可选的进一步逐元素优化是直接按形函数计算上下表面位移跳量，再旋转得到 `δ_n/δ_t`，以及直接填 `B_local`，省去 `B_global` 清零和小矩阵乘法；但 `B_local` 后续仍参与切线与内力，收益必须通过专门微基准证明。
- 当前生产循环更明显的残留开销是 `K_coh[dofs[a],dofs[b]] += ...` 的稀疏标量索引，以及可能的本构小矩阵返回分配；前者正是旧 Task 7 已做过且总收益仅 -1.7% 的一部分。

## 热点根因与内聚力参数的关系

- 源码仍留有“TODO 用户提供实测值”的旧注释，但用户已确认当前 `PCC/NCC.σ_max/K_n/G_c` 是正确实验值；该注释属于后续文档一致性问题，本性能批次不改参数，也不以参数错误解释热点。
- 但现有多个完整运行日志均显示 basic 求解 `converged=true`、`D_max=0`、最大物理分离约 `1e-13–1e-12 m`，同时 CZM 仍占总耗时约 88%–96%。因此当前性能热点在没有损伤演化、没有后峰软化的近线性状态下仍稳定出现。
- 这反证“热点主要由内聚力参数触发强非线性/损伤迭代”不是当前根因。当前直接根因是：约 4 万 DOF 稀疏系统、`update_interval=1` 每个耦合步调用、basic Newton 中对内容常常不变的 `K_bc` 重复执行完整分解。
- 内聚力参数仍会影响二阶性能：`K_n/K_t` 影响条件数与数值主元，`σ_max/G_c/δ_0/δ_c` 影响进入损伤区的时机、切线变化频率、Newton 迭代数和缓存命中率。参数若触发损伤，通常会让矩阵更常变化、缓存收益下降并可能使总耗时更高；但它不会解释当前 D=0 时已经发生的重复分解。
- 需区分两类问题：性能上属于求解算法/调度缺少因子复用；当前实验参数视为正确输入，零损伤结果是否符合本载荷与边界条件不在本性能批次中重新解释或调参。
- Git 读取出现用户级 ignore 文件权限 warning；命令仍成功，且不影响仓库状态判定。后续使用 scoped Git 命令，不尝试修改用户全局配置。

## 技术决策

| 决策 | 理由 |
|---|---|
| 先读现有 CZM 向量化 spec/plan、相关 planning 三件套与当前源码 | 避免重复、冲突或依据过时行号重构 |
| 先做只读审计，再决定是否实施 | 用户要求的是重构，目标与验收需从仓库治理文档中精确恢复 |
| 不启用子代理 | 用户未明确要求委派/并行代理，当前协作规则禁止主动派生 |
| 优先审计 LU 分解热点，不复做已回退 Task 7 | Task 7 仅 -1.7%，最新历史结论与 HEAD 探针均指向 LU |
| 参数冻结，性能重构不得修改材料参数或损伤规律 | 用户确认参数为正确实验值；保持科学契约 |
| 在 `K_bc` 尺寸、CSC 结构及数值内容完全一致时才复用 `factorize(K_bc)` | 命中时消除重复分解；任何内容变化立即重新分解，保持原求解语义 |

## 2026-09-01 执行恢复

- 当前是普通 checkout，分支 `codex/src-physics-modularization`，不是 linked worktree。
- tracked 改动只有本任务总索引；未跟踪项包含本任务目录与用户输出。后续只对目标源码、测试和任务文档做 scoped diff，不纳入 `.mimosa/` 或 `output/`。
- Superpowers 执行顺序：writing-plans 固化 spec/plan；TDD 建立失败测试；实现 basic 非几何路径的精确内容门控因子缓存；verification-before-completion 运行聚焦门和 `testexample.jl`。
- 既有 2026-08-19 性能 spec 的最高约束仍适用：不改求解器类型、收敛容差、浮点算法选择或科学/API 契约；本批把直接 `K_bc \\ R_bc` 对应的 `factorize(K_bc)` 结果作为缓存对象。
- `solve_czm_basic_step` 在 Newton 循环内完成装配、残差检查、`apply_bc_czm`，随后直接求解；缓存接线点应位于 BC 后、回代前，且只在 `geo_nl=false` 时启用。
- `apply_bc_czm` 每次复制矩阵和右端，并根据当前矩阵对角最大值施加相对罚；缓存判据必须使用它返回的最终 `K_bc`，不能只验证未约束 `K_total`。
- 现有 `test_czm_mech_core.jl` 已建立对象身份缓存和逐位等价的测试风格；新测试应沿用真实 mesh/workspace，不用 mock，也不只检查源码文本。
- `CZMAssemblyWorkspace` 实际定义在 `src/CouplingState.jl`，由 `assembly_workspace(czm_mesh)` 惰性挂到 `czm_mesh.ws`；新增快照矩阵与分解对象字段可沿用当前生命周期，不修改 `CohesiveMesh` 或 `MechState` 接口。
- 当前 basic 接线精确位置为 `src/CzmSolve.jl:231-234`：`apply_bc_czm` 返回最终 `K_bc/R_bc` 后立即直接求解；非几何缓存 helper 只替代第 234 行的分解路径。
- 缓存命中必须比较 `size`、CSC `colptr`、`rowval` 和 `nzval`；数值数组采用 `isequal`，连 NaN 与正负零差异也不误判为同一内容。缓存未命中调用 `factorize(K_bc)`，不写死 UMFPACK/LU 类型。
- spec/plan 自检未发现占位步骤；测试代码覆盖首次建缓存、相同矩阵新 RHS 命中、数值变化失效、结构变化失效和 production basic 接线。修正了文档中的 Julia `\` 运算符转义，避免执行者照抄为双反斜杠。

## 方案 A 实施基线（源码修改前）

- 运行环境：Julia 1.11.2、`--startup-file=no --project=.`、`JULIA_NUM_THREADS=1`、`GKSwstype=100`。
- `example/testexample.jl` 退出码 0；`n_theta=80`、`ne=1682`、`nT=1763`、总时间步数/CallModel calls 均为 19。
- 墙钟 `55.976 s`；CZM `44.804 s`，占 CallModel 计时 `89.94%`，平均 `2358.126 ms/call`。
- 科学锚点：最终电压 `3.9438 V`，容量 `0.0833 Ah`，温度 `298.15–299.00 K`，`D_max=D_mean=0`，最大法向分离 `7.4202e-13 m`，断裂单元 0，环向应力 `-1.8433–4.1344 MPa`，切向剪应力 `-0.26184–0.21666 MPa`。
- 基线全过程 CZM 收敛，损伤保持 0；这为相同切线的跨调用复用提供当前场景证据，但实际命中仍由严格矩阵内容判据逐次裁决。

## TDD 实施结果

- 第一轮 RED：退出码 1；helper 与两个 workspace 字段不存在的 3 个预期失败，0 error。
- helper 实现后：首次创建/相同矩阵复用 8/8，通过；数值与 CSC 结构变化失效 5/5，通过。
- production 接线 RED：helper 共 13 个断言通过；真实 basic 非几何入口只在两个缓存非空断言失败，位移有限。
- 接线 GREEN：三个测试集分别 8/8、5/5、4/4，共 17/17，通过；只接入 basic 非几何路径。
- scoped 源码审查确认改动集中在 `CZMAssemblyWorkspace` 两个字段、两个内部 helper 和 basic 求解选择；没有参数、本构、BC、残差、线搜索或状态提交差异。
- 现有测试中直接调用 `solve_czm_basic_step` 的文件只有 `unit_czm_newton.jl`；另以 `test_czm_mech_core.jl` 覆盖 workspace/装配契约，并继续搜索统一 `solve_czm_step` 的间接调用范围。
- 聚焦回归一次运行 6 个文件并全部退出码 0：新缓存测试、`unit_czm_newton.jl`、机械核心、geo C1、J2 integration、arc-geo。由此同时覆盖 basic 命中、workspace 构造、几何非线性旁路以及未改动的 J2/弧长分支。

## 方案 A 后置快速基线

- `example/testexample.jl` 退出码 0；网格仍为 `n_theta=80`、`ne=1682`、`nT=1763`，总时间步数/CallModel calls 仍为 19。
- 科学输出按控制台记录精度与前置锚点一致：最终电压 `3.9438 V`、容量 `0.0833 Ah`、温度 `298.15–299.00 K`、`D_max=D_mean=0`、最大法向分离 `7.4202e-13 m`、断裂单元 0、环向应力 `-1.8433–4.1344 MPa`、切向剪应力 `-0.26184–0.21666 MPa`。
- CZM 从 `44.804 s` 降为 `6.943 s`，减少约 `84.50%`；每次调用平均从 `2358.126 ms` 降为 `365.423 ms`。
- 墙钟从 `55.976 s` 降为 `18.237 s`，减少约 `67.42%`。墙钟只作本机性能证据，不作为科学一致性判据。
- `example/testexample.jl` 不生成 `output/testexample/metrics.toml`；强制门档案是 `Simplify/baseline/testexample/metrics.toml`。本轮控制台的 mesh/result 各字段逐项匹配该档案记录精度；档案明确 `timing_exact=false`。

## 最终范围审查

- tracked 源码修改只有 `src/CouplingState.jl` 与 `src/CzmSolve.jl`；新增测试为 `test/test_czm_factorization_cache.jl`。
- `src/parameters/` 与 `src/SetParams.jl` 差异数为 0，确认正确实验参数未修改。
- scoped whitespace issue 数为 0；LF→CRLF 与用户级 Git ignore 权限提示是环境 warning，不影响检查退出状态。
- `.mimosa/`、`output/` 和其他用户未跟踪内容未纳入实现差异，也未清理。
- 方案 A 达成 spec：严格内容门、原生 `factorize`、basic 非几何单路径、异常与状态提交语义保持、科学快速门一致，并取得可测量提速。

## 问题与处理

| 问题 | 处理 |
|---|---|
| 初始记录称 `superpowers` 不可用 | 已安装 14 个官方技能并在本轮纠正记录；继续遵守 brainstorming 的设计批准门 |
| Git 无法读取 `C:\Users\user\.config\git\ignore` | 记录为环境 warning；不改全局配置，使用当前仓库证据继续 |
| 历史称 `lu(A)\\b == A\\b` 可位级复用，但 Julia 1.11.2 小矩阵实测不成立 | 不采用该论证；检查 `factorize`/生产 dispatch，并以真实 `K_bc` 的逐位实验作门 |
| 显式 `lu` 可能改变直接 `\\` 的算法选择 | 使用 `factorize(K_bc)` 保持直接求解 dispatch；真实生产矩阵仍须逐位 A/B 门 |
| 真实矩阵 A/B 探针首次 `include_string` 路径错误 | 保留原探针绝对文件名作为 include 上下文后重试；不创建临时源码文件 |

## 资源

- `docs/superpowers/plans/2026-04-20-czm-vectorized-solver-plan.md`
- `docs/superpowers/specs/2026-04-20-czm-vectorized-solver-design.md`
- `docs/superpowers/plans/2026-08-19-perf-optimization-spme-thermal-czm.md`
- `docs/superpowers/specs/2026-08-19-perf-optimization-spme-thermal-czm-design.md`
- `docs/planning-with-files/09_向量化CZM/`
- `docs/planning-with-files/01_CZM瓶颈/`

---

每两次搜索/查看后更新本文件。
