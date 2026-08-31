# 位移协调源码审查发现

## 当前已知

- 审查对象是当前工作树，而不是历史提交或技术文档中的目标架构。
- 需要分别判断四类关系：实体网格内部协调、cohesive 双侧位移跳跃、Φ 跨匝协调、SP–涂层接触/摩擦。
- 历史记忆记录过旧快照中 `phi_pairs` 只有拓扑、没有力学消费者；当前源码可能已变化，尚未验证。
- 当前源码已出现明确变化：`Jellyrollmodel.jl` 定义 `merge_phi_pairs` 并生成 `mesh_bonded`，`test_czm_phi_merge.jl` 将其描述为默认 Φ 完美粘结；这推翻了旧快照“没有机械消费者”作为当前事实的可能性，仍需检查后续装配是否确实使用 `mesh_bonded`。
- `Option.jl` 明确标注 `friction_mu` 为预留且当前无消费者；`MechState.contact::Nothing` 也被标注为 SP–涂层接触预留位。
- 当前工作树原本已有与本任务无关的 `src/Initialisation.jl`、一张输出 PNG、多个输出/工具/文档未跟踪项；审查不会改动或清理它们。

## 证据表

| 主题 | 源码/测试锚点 | 事实 | 结论状态 |
|---|---|---|---|
| 网格与 DOF | `src/Jellyrollmodel.jl:681-716`, `src/CzmMesh.jl:79-127` | bulk 初始共享界面节点；仅真实 PE–PCC/NE–NCC 面复制节点并重写外层单元连接 | 已确认拓扑语义 |
| cohesive 位移跳跃 | `src/CzmMesh.jl:96-175`, `src/czm.jl:66-101, 118-289` | cohesive 两侧坐标相同但节点号/DOF 不同；`B_global=[-N_bottom,+N_top]`，经局部 `R=[n;t]` 得到 top−bottom 跳跃，牵引以 `BᵀT` 成对回装 | 已确认弱式协调 |
| Φ 跨匝约束 | `src/Jellyrollmodel.jl:591-630, 726-748`, `src/CzmMesh.jl:35`, `src/CouplingState.jl:238` | `outer_n` 被重映射为对应 `inner_n`，删除 outer 节点并连续重编号；生产网格创建和 MechState 初始化均指向 `mesh_bonded` | 已确认拓扑消元为共享 DOF，待检查 cohesive 再分裂 |
| 边界条件 | `src/CzmBC.jl`, `src/CzmSolve.jl:59-122` | 内圈（由 `fix_inner` 控制）和外圈节点的 x/y DOF 置零；求解线性系统时用相对对角罚项，并在试探位移上显式回写 Dirichlet 值 | 待核对生产调用一致性 |
| 求解与状态提交 | `src/CzmSolve.jl:175-277, 280-667, 670-988`; `src/CouplingState.jl:390-475`; `src/Solve.jl:241-275` | 所有调度分支从 `ms.u_prev` 起步并共享总系统；仅收敛提交，主耦合入口对非收敛直接报错，不记录失败位移 | 已确认 |
| SP 接触/摩擦 | `src/Option.jl:55`, `src/CouplingState.jl:207`; 全局引用搜索 | `friction_mu` 仅在默认值测试出现，`contact` 仅为 `Nothing` 字段；无残差、切线或状态消费者 | 已确认未实现 |

## 暂不下结论的事项

- `phi_pairs` 的存在不等于它已进入机械残差或约束装配。
- cohesive 单元共享或重复节点的外观不等于位移连续；必须检查连接表、局部 DOF 提取及牵引残差。
- `friction_mu` 选项的存在不等于 Coulomb 接触已实现。

## 网格协调事实

- `build_czm_submesh` 先建立 `(n_layers+1)` 条螺旋节点线；相邻 Q4 层在公共径向边界直接引用同一条螺旋的相同节点号，因此未经界面分裂的位置天然共享两个平移 DOF。
- Φ 处理不是罚函数、Lagrange 乘子或求解后平均：`merge_phi_pairs` 在构网时把每个 outer 节点号替换为 inner 节点号，删除 outer 节点行，再把单元连接重编号。只要后续装配使用 `mesh_bonded`，该节点对的 `u_x/u_y` 就是同一个未知量，属于强式完美粘结。
- `.mesh` 保留未合并物理身份与 `phi_pairs`，`.mesh_bonded` 承担求解拓扑；`phi_keep` 记录 raw→bonded 保留行。此双网格设计必须避免后续误用未合并 `.mesh` 求解。
- `create_czm_mesh` 明确从 `mesh_bonded` 开始，仅识别 PE–PCC 与 NE–NCC 四类真实箔–涂层面；它为公共边节点建立 memoized 副本，并把外层 Q4 的连接替换成副本。这样两侧初始几何重合，但有独立 `u_x/u_y`，允许形成非零位移跳跃。
- SP–涂层以及同材细分边界未被选为 cohesive 面，仍沿用公共节点的强式连续；这不是接触语义，而是完美粘结语义。
- cohesive 单元节点顺序固定为 `[bottom_lo, bottom_hi, top_hi, top_lo]`，并另存一致切向顺序的 `nodes_bottom=[lo,hi]`、`nodes_top=[lo_copy,hi_copy]`，为局部位移跳跃计算提供双侧 DOF。

## cohesive 方程事实

- 局部标架中切向 `t` 沿 θ 递增边，法向 `n` 不依赖固定左法线，而由内层单元质心指向外层单元质心校正，避免把径向张开误判成压缩。
- cohesive DOF 顺序为底面两节点后接顶面两节点；`B_global` 对底面赋 `-N1/-N2`，对顶面赋 `+N2/+N1`，所以插值得到的是 `u_top-u_bottom`，再由 `R` 分解为法向 `δ_n` 和切向 `δ_t`。
- 分离进入双线性牵引–分离律与算法切线；内力 `B_localᵀ[T_n,T_t]` 回装到同一 8 DOF。由于 B 两侧符号相反，界面牵引在上下表面形成作用—反作用，提供弱式位移协调，而不是强制 `u_top=u_bottom`。
- 当界面刚度有限、软化或损伤增长时允许位移跳跃；只有未分裂的共享节点/Φ 合并节点是严格零跳跃。
- 线弹性 bulk 由 Q4 节点号生成 `2n-1/2n` 全局 DOF；几何非线性/J2 分支仍对同一 `bulk_element` 和相同全局 DOF 装配，仅改变单元残差/切线和材料状态，不改变协调拓扑。

## 总系统与边界事实

- 总系统在同一个 `ndof=2*nnode` 空间中相加：`K_total=K_bulk+K_coh`、`f_int=f_bulk+f_coh`，残差为外力与热化学等效力减总内力。因此 bulk 位移场与界面分离不是两个串联求解器，而是同一机械平衡方程中的共同未知量。
- 生产边界识别将外边界节点固定 x/y；`fix_inner=true` 时也固定内边界节点。该约束用于去除/限定整体运动和外形边界，并不替代层间共享 DOF 或 cohesive 牵引协调。
- `apply_bc_czm` 采用 `penalty=1e6*maximum(abs,diag(K))` 的相对罚法修改线性系统；求解器另外用 `apply_czm_dirichlet!` 把试探/更新后的受约束位移精确回写为目标值。需要继续检查每个 basic/load-substep/arc-length 分支是否都执行这两层处理。

## 求解与状态事实

- basic、载荷子步和旧 arc-length 分支都从上一收敛态 `copy(ms.u_prev)` 开始，使用相同 `assemble_coupled_system`；线搜索或弧长只改变增量求解/路径控制，不改变节点共享和 cohesive 跳跃定义。
- 这些分支都把受约束 DOF 的残差写成 `val-u[dof]`，在求解修正时调用 `apply_bc_czm`，并在接受试探位移后显式施加 Dirichlet 值。
- damage 在位移平衡迭代中冻结，整体/子步收敛后再按界面分组更新；这是一种交错的状态更新，不是额外的位移协调方程。
- `MechState` 的 `u_prev` 与 `damage_states` 只在结果 `converged=true` 时原位提交。basic 失败会返回起始位移；旧 arc/load-substep 失败可能在 `CZMResult.displacement` 暴露未完成载荷路径的试探/部分收敛位移，但不会污染 `ms.u_prev`。上层是否会误用失败结果需检查。
- 当前统一调度 `solve_czm_step` 根据 `iter_method` 选择 basic/load_substep/arc_length；若同时开启几何非线性，arc_length 改走严格的 `solve_czm_arc_geo_step`。该分支在参考态、每个弧长步和最终 λ=1 平衡中都维持相同 BC 与总系统，失败抛错并回滚子步，λ=1 验收后才提交。
- 生产耦合入口 `update_czm_damage!` 对所有输出做有限性检查，并对 `converged=false` 抛出异常；`Solve.jl` 只有在该调用成功返回后才保存 CZM 快照和计算宏观应力。因此旧 arc/load-substep 的失败结果不会在正常主流程中成为有效输出。
- `MechState.contact` 固定为 `Nothing`；塑性状态仅在机械平衡收敛后通过一次 `commit_plastic=true` 的一致装配提交。状态提交改变本构历史，不改变节点/界面协调拓扑。

## 现有验证覆盖

- `test/test_czm_phi_merge.jl` 验证物理 Φ 对坐标重合、outer 原节点被删除、inner 保留、bonded 节点数准确减少、连接索引连续，以及装配矩阵有限/近似对称。测试标题称“绑定对位移恒等”，但数值断言本质是拓扑唯一 DOF 与装配烟雾检查，没有保留两个独立 DOF 可直接相减；零相对位移来自节点合并的构造性保证。
- `test/unit_czm_bilinear.jl` 对复制后的 cohesive 顶面施加独立法向/切向位移，验证 `δ_n/δ_t`、解析双线性牵引、损伤起始、卸载/重载不可逆性以及接近断裂。这直接覆盖“允许跳跃并由牵引协调”的单元级行为。
- `Mechanical.thermal_diffusion_stress_2D` 的 czm-off 宏观力学路径直接在 `mesh_bonded` 上求解，不创建复制界面节点或 cohesive 单元。因此该路径的所有层间面（包括 PE–PCC/NE–NCC）都是共享 DOF 完美粘结；它不能表现 CZM 分离/损伤。
- `identify_boundary_nodes` 对实际传入的 `Mesh.nlen` 或扩展后 `CohesiveMesh.nnode` 逐节点做几何边界判定，所以生产 CZM BC 作用在真实求解节点空间；不是拿未合并 raw 网格索引套到 bonded/复制后的 DOF 上。

## 验证结果

- 2026-08-31 运行 `test/test_czm_phi_merge.jl`（Julia 1.11.2，单线程项目环境，`GKSwstype=100`）：退出码 0；三个 testset 分别 8/8、2/2、1/1 通过。
- 2026-08-31 运行 `test/unit_czm_bilinear.jl`（同环境）：退出码 0；Mode I 单调 90/90、卸载/重载 83/83、Mode II 切向 61/61，共 234/234 通过。Mode II 当前 `model1` 输出 `D_max=0`，该组测试验证切向分离与初始切向牵引，不应被表述为切向损伤已增长。

## 审查结论

当前源码的“位移协调”分为三个层次：

1. **强式拓扑协调**：普通层间公共边与 Φ 跨匝缝使用唯一节点号/唯一 x-y DOF，位移严格相等。
2. **弱式界面协调**：PE–PCC/NE–NCC 真实面复制节点形成双侧 DOF，允许法/切向跳跃；CZM 牵引以作用—反作用形式进入总残差，并以算法切线耦合两侧。
3. **全局边界约束**：内/外圈 Dirichlet 条件限制整体边界运动；它独立于层间协调机制。

未实现的第四层是 SP–涂层单边接触/摩擦。当前这类界面仍是共享节点完美粘结，无法开缝、闭合、滑移或发生 stick/slip；不得把 `phi_pairs`、`contact::Nothing` 或 `friction_mu` 默认值当作接触能力。
