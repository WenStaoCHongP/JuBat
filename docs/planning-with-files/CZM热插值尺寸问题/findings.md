# CZM 热插值尺寸问题发现

## 用户要求

- 当前 `Solve.jl`/`CouplingState.jl` 的严格失败修改正确。
- 处理新暴露的尺寸问题，并分析造成原因。
- 修改业务文件前仍需报告具体问题并等待批准。

## 已知现象

- `test/smoke_czm_redesign.jl`：插值矩阵第二维为 2524，热状态长度为 1322。
- `example/testexample.jl`：插值矩阵第二维为 3366，热状态长度为 1763。
- 两组矩阵列数均约为热单元节点数的两倍，提示插值矩阵可能按另一套加密或重复节点集合构造。
- `test/test_thermal_to_czm_interp.jl`、损伤映射与缓存测试均通过，说明局部函数契约不能覆盖端到端网格组合错误。

## 待收集证据

- `mesh_data.thermal2D`、`mesh_data.czm_submesh`、`case.mesh["thermal2D"]`、`case.czm_mesh` 的节点数与对象来源。
- `create_czm_mesh` 调用和 `build_thermal_to_czm_interp` 的参数顺序与构造时机。
- `compute_czm_strain_inputs` 使用的矩阵及温度向量来源。
- 缓存键是否包含热网格身份/节点数。

## 调用链证据

1. `setup_thermal2D_mesh` 在默认路径强制 `use_merged=true`，并将 `case.mesh["thermal2D"]` 设为 `mesh_data.thermal2D_merged`。
2. `create_czm_mesh(czm_submesh, thermal_mesh, param)` 使用传入的 `thermal_mesh` 调用 `build_thermal_to_czm_interp`；矩阵列数直接等于该网格的 `nlen`。
3. `example/testexample.jl:85` 与 `test/smoke_czm_redesign.jl:37` 在 setup 后仍传入 `mesh_data.thermal2D`（未合并网格），而主求解温度向量来自 `case.mesh["thermal2D"]`（合并网格）。
4. 正确的定向测试均传入 `case.mesh["thermal2D"]`，因此矩阵列数与温度状态一致并通过。
5. 多个网格敏感性示例存在与主例相同的错误调用，说明问题不是缓存陈旧，而是调用方同时持有两套热网格时选错对象。

## 当前根因假设

- 2026-07-21 的 v2 修改将默认活动热网格改为强制合并网格，但若干 CZM 示例仍沿用旧调用，把未合并热网格传给 `create_czm_mesh`。
- `create_czm_mesh` 没有能力知道哪一个网格是 `case` 的活动热网格，也没有在求解前校验 `size(thermal_to_czm, 2) == length(T_nodes)`；错误因此直到矩阵乘法才暴露，随后又曾被 `Solve.jl` 的 catch 吞掉。

## 运行时尺寸证据

| 算例网格 | 未合并节点 | 合并节点 | 活动节点 | 插值矩阵（错误输入） | 插值矩阵（活动网格） | 热单元数（两者） |
|---|---:|---:|---:|---|---|---:|
| `nθ=60, nθ_czm=20` | 2524 | 1322 | 1322 | `3969 × 2524` | `3969 × 1322` | 1261 |
| `nθ=80, nθ_czm=80` | 3366 | 1763 | 1763 | `15849 × 3366` | `15849 × 1763` | 1682 |

- 失败堆栈中的 2524/1322 与 3366/1763 精确等于上述未合并/活动节点数，不是近似关系。
- 合并前后热单元行数和编号范围不变；`thermal_elem_map` 的最大值分别为 1261、1682，均合法。因此无需重建单元映射，只需让节点插值矩阵绑定活动合并网格。
- `compute_czm_strain_inputs` 当前在 `M * T_nodes` 后才断言输出行数，缺少乘法前的列数不变量检查。

## 理论复核（进行中）

### `Theory/04_CZM.md` §3.0a–§3.3

- C-skip-thermal 明确规定：CZM 子网格不求解温度，节点温度必须由“粗热网格”的唯一活动温度场正向插值得到。
- `thermal_elem_map` 是 CZM 子网格单元到粗热单元的 O(1) 映射；它服务于单元温升和 SOC 映射，与节点插值矩阵是两条不同映射链。
- 理论上的粗热网格是热方程实际求解所在的网格，不是网格生成器中任意一份几何相近的候选网格。因此 `thermal_to_czm` 的列空间必须严格等于活动热状态空间。
- 当前使用未合并网格建矩阵、再乘合并网格温度向量，不仅尺寸错误，也违反“单一粗热温度场 → CZM 子网格”的理论接口。
- 合并前后单元编号一致支持继续使用 `thermal_elem_map`；但不能据此推导节点自由度也等价，因为合并操作改变了节点空间维数和连续性约束。

## 受影响调用点

- 当前耦合执行路径共发现 10 个文件、11 处错误调用：8 个公开示例、1 个 CZM 冒烟测试，以及 `tools/check_czm_methods_coupled.jl` 中 2 处。
- `docs/planning-with-files/损伤异常/` 下另有 6 个可执行诊断探针沿用同一错误调用；如果重跑也会失败，建议在同一机械修复批次同步纠正。
- 正常通过的单元测试统一传入 `case.mesh["thermal2D"]`，可作为修复范式。
- `example/coupled_czm_thermal_example.jl` 在 `setup_thermal2D_mesh` 之前构造 CZM 网格，不能只替换参数，还必须把 CZM 构造移动到活动热网格设置之后。
- `example/循环验证/czm_cycle_example.jl` 将单独构造的 `czm_mesh` 传给 `solve_cycling`，也必须绑定活动热网格。
- `tools/czm_baseline_probe.jl`、`tools/czm_convergence_diag.jl`、`tools/verify_czm_standalone.jl` 和 `tools/verify_czm_system.jl` 是纯 CZM standalone 路径，没有 `setup_thermal2D_mesh` 生成的活动热状态；它们使用未合并网格是合法的，不纳入修复。

## 技术决策

| 决策 | 理由 |
|---|---|
| 暂不修改业务代码 | 根因尚未确认，避免用裁剪或补零再次掩盖网格不变量错误 |
| 不新增绑定 `Case` 的薄包装入口 | 现有迁移文档和通过的测试已经规定正确调用为 `create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)`；新增包装会扩大 API 而不能修复旧三参数调用 |

## 候选修复结构

1. 将 10 个当前耦合文件和 6 个可执行历史探针中的错误第二参数改为 `case.mesh["thermal2D"]`；耦合示例同时把 CZM 网格构造移动到 setup 之后。
2. 在 `compute_czm_strain_inputs` 的矩阵乘法前增加严格尺寸不变量，错误信息同时报告矩阵列数、活动热网格节点数和温度向量长度；不做裁剪、补零或重建回退。
3. 在 `test/test_thermal_to_czm_interp.jl` 增加负向断言：把未合并网格构造的 CZM 网格绑定到合并热 case 时，必须在矩阵乘法前抛出 `DimensionMismatch`。
4. 用插值/应变定向测试、CZM 冒烟和 `example/testexample.jl` 验证。

## 问题记录

| 问题 | 处理 |
|---|---|
| 一条组合文档检索因部分模式无匹配返回退出码 1 | 已获得目标调用链；改用精确文件/字符串查询，不重复原组合命令 |
| 首次规划记录补丁引用了不存在的模板空表格 | 读取当前文件后改用准确锚点的小补丁 |
| 一条全目录文档组合检索有结果但因无完整匹配返回退出码 1 | 保存已返回的迁移文档证据，后续改用 `--glob '*.jl'` 精确枚举可执行脚本 |

## 拟修改文件（等待批准）

- 生产守卫：`src/CouplingState.jl`。
- 回归测试：`test/test_thermal_to_czm_interp.jl`、`test/smoke_czm_redesign.jl`。
- 当前示例/工具：`example/testexample.jl`、`example/coupled_czm_thermal_example.jl`、`example/循环验证/czm_cycle_example.jl`、`example/内聚力验证/czm_grid_convergence.jl`、两代 4/5 网格敏感性脚本、`tools/check_czm_methods_coupled.jl`。
- 可执行历史探针：`docs/planning-with-files/损伤异常/` 下 6 个 CZM 探针。
- 不修改 `CzmMesh.jl`、网格生成算法或 standalone 工具。

## 资源

- `src/Jellyrollmodel.jl`
- `src/CzmMesh.jl`
- `src/CouplingState.jl`
- `src/SetCase.jl`
- `test/smoke_czm_redesign.jl`
- `example/testexample.jl`

### `Theory/04_CZM.md` §3.6.0a–§3.6.3 复核补充

- 式 (3.95a) 明确定义 `M ∈ R^(N_czm_nodes × N_thermal_nodes)`，所以 `size(M, 2)` 必须等于当前热状态向量的节点自由度数；这不是可回退的输入差异，而是耦合接口不变量。
- 理论要求温度节点映射、单元温升映射和 SOC 映射各司其职：`T_czm_nodes = M*T_thermal_nodes`，而 `dT_czm[e]`、PE/NE 的 `Δsoc[e]` 通过 `thermal_elem_map[e]` 对父热单元直接索引。不能用“单元编号相同”替代节点空间一致性。
- “合并网格必须用于 CZM”不是理论结论。准确结论是：构造 `M` 时必须使用本次热方程实际求解所在的活动网格。当前配置的活动网格恰好是 merged；若以后恢复未合并热网格求解，正确输入也应随 `case.mesh["thermal2D"]` 改变。
- 式 (3.95b)–(3.95c) 把 `thermal_elem_map` 定义为每个 CZM 实体单元到一个有效父热单元的全映射。现实现若对越界映射 `if` 跳过并保留零温升/SOC，会制造无理论依据的零载荷，应改成统一边界校验后直接索引。
- 式 (3.95c) 还规定 PE/NE 的 `Δsoc` 是必需耦合输入；当前缺键时先生成零浓度再减 `cs0`，会得到 `Δsoc=-cs0` 的伪化学载荷，而不是零载荷。缺键或长度不符时应立即报错。
- 式 (3.95d)–(3.95f) 表明损伤反馈是 CZM 接口到活动粗热单元的 max-reduction，再进入 `k_n_eff(D)` 等系数；这进一步说明节点插值矩阵与损伤归并矩阵不能混为一条映射。

### `Theory/06_热源.md` §5.1–§5.6

- 热方程只求解一个活动温度场 `T(s,n,t)`；§5.5 的热应变由该场的局部 `ΔT = T - T_ref` 唯一生成。理论没有第二套“未合并节点温度”可供 CZM 使用，所以从未参与求解的候选节点空间构造 `M` 没有物理意义。
- 损伤对热场的反馈通过粗热单元上的 `D → k_n_eff(D)` 和 `D → Q_contact(D)` 完成；这与 §3.6.0b 的反向归并一致，也说明 merged/unmerged 的选择应由活动热离散决定，而不是由 CZM 调用方自行猜测。
- §5.4 明确标注旧 `x(r)` 公式是 A 方案残留、当前多 SPMe 架构应使用逐热单元体积平均热源；因此本次诊断不应借用 §5.4 的旧径向映射为节点插值提供依据。
- 本次排查时 Theory 曾混用 `N_lay` 与 `N_coh=2`；后续“界面术语统一”任务已改为 `N_face,repeat^coh=4`（真实面数），并与 2 种本构类型分开。该历史文档漂移与 `M` 的列数和活动热节点空间无关，不改变本次尺寸根因。
- §5.5 强调空间非均匀 `ΔT` 才是热应力驱动；任何尺寸或父映射错误后用零载荷回退都会抹除局部热应力，直接造成损伤结果失真，必须 fail-fast。

## 理论结合现代码后的复评结论

1. 原诊断的直接根因正确且证据闭合：失败调用用未合并候选网格构造 `M`，活动温度状态来自 merged 网格，因此矩阵定义域与状态空间不同。
2. 诊断措辞应从“CZM 必须使用合并网格”修正为“CZM 插值必须绑定 `case.mesh["thermal2D"]` 所代表的活动热网格”。merged 只是当前 `setup_thermal2D_mesh` 默认配置的结果，不是永久理论要求。
3. `T_czm_nodes` 当前只被 `compute_czm_strain_inputs` 返回，未被后续 CZM 装配读取；但式 (3.95a) 明确要求该映射，删除这次矩阵乘法只会隐藏错误调用，不能作为修复。
4. `Jellyrollmodel.jl` 中“`thermal_elem_map` 索引到未合并粗热单元”的注释过时。该映射实际依赖的是合并前后保持不变的热单元编号；节点是否合并与父单元编号是两个不变量。
5. 原候选方案只增加矩阵列数守卫仍不完整。同一生产入口还应严格校验：`length(T_nodes)==active_mesh.nlen`、`size(M)==(czm_nodes, active_mesh.nlen)`、`length(thermal_elem_map)==ne_czm`、所有父单元索引在 `1:n_thermal_elem` 内、SOC 键存在且选取当前列后长度等于 `n_thermal_elem`。
6. 这些校验通过后应直接索引，删除越界 `if` 跳过和 SOC 缺键默认值；不得裁剪、补零、重建或继续计算。
7. 调用方仍应统一传入活动网格；生产守卫是防止未来错误重新潜入，不替代调用点修正。

## `CouplingState.jl` 实施设计（已批准）

- `compute_czm_strain_inputs` 是尺寸错误首次进入运算的位置，也是式 (3.95a)–(3.95c) 三条输入链的共同入口，严格校验集中放在该函数中。
- 当前 `update_czm_damage!` 的有限值检查位于 `compute_czm_strain_inputs` 调用之后，无法保护矩阵乘法、温度索引或 SOC 索引；输入有限值与尺寸检查必须前移到共同入口。
- 保留 SOC 向量与历史矩阵两种已存在的表示；若为矩阵，仍选当前最后一列，但在选择后严格要求长度等于活动热单元数。
- 错误类型按语义区分：缺失必需键用 `KeyError`，尺寸契约违例用 `DimensionMismatch`，非法父单元值或非有限输入用 `ArgumentError`。
- 现有正向测试使用活动热网格和长度匹配的 SOC 向量，应在严格校验下继续通过；本阶段不改测试文件。
- 已实现的校验顺序保证 `M*T_nodes` 之前完成活动节点数、矩阵行列数、父单元映射、Q4 连接、SOC 表示/长度及有限值检查；通过后所有父单元均直接索引。
- 现有错误调用现在会在 `M*T_nodes` 前准确报告两套节点空间，不再依赖 Julia 矩阵乘法的通用错误，也不会进入任何热/化学应变计算。
- 正向测试确认材料分发科学结果未改变；矩阵探针确认历史 SOC 仍取最后一列。
- `test/smoke_czm_redesign.jl` 无需重排：其 `setup_thermal2D_mesh` 已先执行，单行改用 `case.mesh["thermal2D"]` 后完整仿真通过并越过原尺寸断点。
- `example/testexample.jl` 单行修复后高分辨率 CZM 实际执行 19 次并全部收敛；电化学/热学/损伤打印指标保持冻结值，只有此前被失败路径压成零的法向分离恢复为 `1.3527e-14 m`，从而改变自动缩放的结果图及 PNG 哈希。
## 下一调用方：`example/coupled_czm_thermal_example.jl`

- 第 104 行在 `setup_thermal2D_mesh` 之前，用候选网格 `mesh_data.thermal2D` 创建 `czm_mesh`；第 118 行之后才得到活动热网格 `case.mesh["thermal2D"]`。
- 该顺序会让 CZM 插值矩阵继续绑定未合并候选节点，而求解状态绑定活动合并节点，复现本轮已确认的尺寸来源错配。
- 修复必须同时调整顺序与参数：先完成 Case 的热网格安装，再用活动热网格创建 CZM；不能只增加尺寸回退或吞掉错误。
