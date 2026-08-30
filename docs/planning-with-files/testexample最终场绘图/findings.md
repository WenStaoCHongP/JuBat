# testexample 最终场绘图：发现与决策

## 用户要求

- 修改 `example/testexample.jl`。
- 删除其他绘图代码。
- 仅绘制最终时刻温度场、环向应力场和切向应力场分布图。
- 后续明确要求不用散点图，改为云图形式。
- 最新要求：应力云图必须使用力学/CZM bulk 子网格，而不是热学粗网格；温度图未要求改变。

## 已知项目约束

- `example/testexample.jl` 是强制行为基线入口，不能改变网格、步数或科学结果。
- 绘图输出路径应为 `output/testexample/`，并基于 `@__DIR__` 构造。
- 缺失或非法工程状态必须显式失败，不能补零或静默回退。
- 历史验证将 `example/testexample.jl` 与冻结 PNG 哈希一起作为门禁；本次用户明确要求改变绘图产物，因此旧 PNG 哈希预期失效，但求解与科学结果仍须保持不变。

## 待确认

- 脚本现有全部绘图调用及生成文件名。
- 温度、环向应力、切向应力最终帧的索引方式和单位恢复路径。
- 删除非目标绘图后是否有无用导入、变量或辅助函数可一并窄范围移除。

## 现状发现

- 当前脚本使用 `Plots`，仅生成 `output/testexample/testexample_results.png`。
- 当前图为三联时间历程：温度/电压、CZM 分离、CZM 损伤/容量；并非空间场分布。
- 用户要求的三个最终场需要从现有热/力学场后处理接口取得，不能从上述时间历程图直接裁剪得到。
- 当前脚本的 `opt.mechanicalmodel = "none"`，但 `opt.czm_enabled = true`；宏观力学/CZM 位移与应力可能仍由 CZM 路径计算，必须核实实际结果键。
- 温度节点时序的既有输出键为 `"thermal2D temperature at nodes [K]"`，并另有最终节点温度键；应优先使用有量纲最终场输出。
- 代码中“particle surface tangential stress”属于颗粒尺度扩散应力，不应未经确认当作果冻卷横截面的宏观环向应力场。
- “环向应力”和“切向应力”在普通极坐标语境可能同义；需从项目现有宏观应力命名判断用户所指的两个独立分量，避免凭名称臆造。
- `example/jellyroll_stress_displacement.jl` 已有可复用的 Q4 单元填色绘图函数和最终场重建范式。
- 该示例通过 `JuBat.thermal_diffusion_stress_2D` 得到宏观涂层尺度的 `σ_xx/σ_yy/σ_xy`，并用 `scale.E_coat` 恢复为 Pa；这与颗粒尺度应力明确不同。
- `thermal_diffusion_stress_2D` 目前不直接返回极坐标分量；若目标为宏观环向/切向分量，需对每个单元质心处的全局应力张量作明确坐标变换。
- 理论坐标中 `σ_ss` 是切向/环向正应力，`σ_sn` 是切向剪应力；为形成两个不重复的目标场，本次将输出明确命名为环向正应力 `σ_θθ` 与切向剪应力 `τ_rθ`。
- 全局到局部的张量旋转采用单元质心径向单位向量 `n=(cosθ,sinθ)` 与环向单位向量 `t=(-sinθ,cosθ)`：`σ_θθ=tᵀσt`，`τ_rθ=nᵀσt`。
- `NormaliseParam` 将 `cell.T0` 除以 `scale.T_ref`，并将热膨胀系数乘以 `scale.T_ref`；因此 `thermal_diffusion_stress_2D` 的 `T_nodes` 输入必须是无量纲温度。
- 现有 `jellyroll_stress_displacement.jl` 直接把 `[K]` 温度传入该函数，与其内部归一化契约不一致；本次不会复制这一单位错误，而是将最终节点温度除以 `scale.T_ref` 后重建应力。
- `Solve` 直接提供 `"thermal2D final temperature at nodes [K]"`；单元 SOC 历史在后处理中无量纲直传，可用最后一列重建最终宏观应力。
- 当前 `PostProcessing` 仅输出 CZM 位移、界面牵引/分离/损伤，并不输出 bulk 单元应力；因此目标宏观应力图采用项目已有的 `thermal_diffusion_stress_2D` 后处理路径，并明确属于涂层尺度 Q4 场，不声称是 CZM 界面牵引场。
- 热网格坐标按项目约定为 `x*=x/L`；绘图时将坐标乘 `scale.L`，轴标签才可写为米。
- 完整算例保持冻结科学指标：1682 个热单元、1763 个节点、19 步、最终电压 3.9438 V、最终容量 0.0833 Ah、温度 298.15–299.00 K、最大法向分离 1.5174e-12 m、零损伤。
- 首次图像核对显示三张图均完整可读；切向剪应力范围约 `-17.861–25.699 MPa`，小于环向应力 `30.760–94.340 MPa`，因此两个应力场应分别使用各自以零为中心的对称色标。
- `output/testexample/testexample_results.png` 是 2026-08-24 的冻结基线产物，本次不删除；新脚本中已不存在生成该图的调用。
- 项目记忆指向现有显式 8 层 Q4/CZM 拓扑及 `thermal_elem_map`；当前续改应现场核实字段所在类型、长度和索引范围后再使用。
- 现场源码确认 `CzmSubmesh.thermal_elem_map[e_mech]` 定义为每个力学 bulk 单元到父热单元的直接索引；同一映射已被 `compute_czm_strain_inputs` 用于温度与 SOC 耦合。
- `create_czm_mesh` 以 `mesh_bonded` 建立最终拓扑，之后为 cohesive 界面复制节点并重写 `CohesiveMesh.bulk_element`；最终力学求解几何应取 `CohesiveMesh.node + bulk_element`，不能误用 `bulk_mesh` 作为最终连接。
- `CohesiveMesh.czm_submesh` 保留原 `thermal_elem_map`，其单元行序与 `bulk_element` 一致；实施前还需用当前 nθ=80 算例数值核对长度和索引范围。
- 当前 nθ=80 数值探针确认：热单元 1682、最终力学 bulk 单元 13,456、力学节点 20,276、`thermal_elem_map` 长度 13,456、父索引范围 1–1682，完整性检查通过。
- 首次完整力学子网格云图运行 exit code 0，冻结数值指标保持一致；视觉检查确认 13,456 个子单元均被绘制，但逐单元边线过密，会削弱云图可读性。

## 后续发现（2026-08-29）：环向应力图不出现正负极“一拉一压”的原因

> **2026-08-29 参数更正注记**：本节分析基于当时 `PE.Omega = −7.28e-7`（未审 placeholder）。同日经用户决策更正为 `+7.88e-7`（嵌锂膨胀），并按“层分辨应力求解”任务实现层分辨应力（耦合在线导出 + 固体工具函数）后，涂层一拉一压（放电 NE 拉 / PE 压）已经出现；本节“抹平”机理分析对旧粗网格实现仍然成立，见 `docs/planning-with-files/层分辨应力求解/`。

用户提问：为什么 `final_hoop_stress_field.png` 观察不到正负极涂层一拉一压的应力分布（图为近似均匀拉应力）。本节为分析结论，未改动任何代码。

- 结论：这是当前宏观应力计算链路的构造性结果，不是绘图或坐标旋转错误；层间拉压对比在进入有限元求解前就已被抹平，现有模型在该图上不可能呈现 PE/NE 交替拉压。
- 原因一（特征应变求和抹平）：`thermal_diffusion_stress_2D` 在热粗网格上求解，每个热 Q4 单元横跨一个完整 8 层卷绕重复单元，单元仅有一个标量特征应变 `ε₀ = α_eff·ΔT + β_n·Δsoc_n + β_p·Δsoc_p`（`src/Mechanical.jl` `epsilon_0_elem` 装配处）；PE 与 NE 的膨胀贡献在求解前已加成为一个数。
- 原因二（单一材料）：整个域共用一组 `E_eff/ν_eff/α_eff`，全部取 `:PE_PCC` 界面参数占位；`src/Mechanical.jl` 该处注释明确"per-interface 化由 CZM 路径负责"，即分层刚度尚未进入 bulk 应力求解（`compute_czm_params_per_interface` 目前只服务 CZM 界面参数）。
- 原因三（力学子网格仅几何加密）：绘图把父热单元应力经 `thermal_elem_map` 直接继承到 13,456 个力学 bulk 单元（13,456 = 8 × 1682）；同一重复单元内 8 个分层 bulk 单元的应力张量完全相同，图上相邻 PE/NE 涂层必然同色。
- 原因四（载荷同号 + 全约束边界）：`src/parameters/Jellyroll.jl` 中 `NE.Omega = +3.1e-6`、`PE.Omega = -7.28e-7`；放电时 Δsoc_n<0（脱锂收缩）、Δsoc_p>0 且 Ω_p<0（嵌锂也收缩），两项贡献永远同号叠加，全场 ε₀ 同号（放电收缩、充电膨胀）。叠加内外圈节点全部 `fixed_xy` 全约束（`identify_boundary_nodes` + 罚刚度），受约束的整体收缩转化为近似均匀的环向拉应力（本次算例约 +30~94 MPa），即图上那片均匀拉应力；径向微弱梯度来自逐单元 SPMe 电流分布差异。
- 量纲链核实：`variables["thermal2D element soc_n/p"]` 存归一化化学计量比（`src/CallModel.jl` 对颗粒浓度取均值，`SetParams` 中 `cs0 = cs0/cs_max`）；`Omega` 归一化为 `Ω·cs_max`，故 `β·Δsoc` 量纲自洽，排除"单位错误导致看不到拉压"的假设。
- 该算例 `opt.mechanicalmodel = "none"`，颗粒尺度扩散应力（`Calstressdisp`，颗粒表面拉/中心压）也完全未参与计算。
- 当前均匀拉应力作为"单一材料抹平模型 + 全约束边界受约束收缩"的解本身自洽，不是 bug；它回答的不是层分辨力学问题。
- 若要看到 PE/NE 一拉一压（如放电时 NE 相对受拉、PE 相对受压），需实现层分辨 bulk 力学：在现有 `CzmSubmesh`/`thermal_elem_map` 拓扑上给 8 层各自分配刚度与各自特征应变（NE 涂层 β_n·Δsoc_n、PE 涂层 β_p·Δsoc_p、SP/集流体仅热应变）后求解平衡方程；即 `src/Mechanical.jl` 注释中推迟的 per-interface 化扩展到 bulk 应力场。
- 短期替代方案：可分别绘制分层特征应变场 `β_n·Δsoc_n` 与 `β_p·Δsoc_p` 的空间分布对比图，先行展示层间膨胀差，无需等层分辨求解器落地。

## 技术决策

| 决策 | 理由 |
|---|---|
| 仅删除绘图及其专属准备代码 | 保持求解和科学契约不变 |
| 复用现有 Q4 场填色逻辑 | 避免引入未经验证的新绘图库或插值 |
| 将第二个应力图明确标为切向剪应力 `τ_rθ` | 避免把与环向正应力同义的“切向应力”重复绘制 |
| 使用最终温度专用结果键与 SOC 历史最后一列 | 明确满足“最终时刻”，避免近邻时间索引 |
| 应力图标注为宏观涂层尺度 Q4 后处理 | 区分于颗粒应力和 CZM 界面牵引 |
| 两个应力图使用独立的零中心对称色标 | 保留正负号语义并提高较小切向剪应力场的可读性 |
| 使用实际 Q4 单元多边形按单元值直接填色 | 形成连续覆盖的云图，不使用散点或规则网格插值 |
| 最终节点温度取 Q4 四节点算术平均后绘制 | 与项目既有节点到单元温度映射契约一致 |
| 温度图保留热网格，应力图使用力学 bulk 子网格 | 精确响应用户最新的绘图网格要求 |
| 应力云图几何取 `case.czm_mesh.node` 与 `bulk_element` | 使用实际力学装配拓扑，包括 cohesive 界面节点复制后的连接 |
| 先映射全局应力分量，再在力学子单元质心旋转 | 让每个子单元使用自身的局部极坐标方向，而不是直接复制父单元旋转结果 |
| 映射后长度和父索引范围均显式校验 | 让拓扑不一致直接失败，不截断、不补值 |
| 保留子单元填色但隐藏 Q4 边线 | 仍严格使用力学子网格，同时避免 13,456 条密集网格边界遮盖云图 |

## 资源

- `example/testexample.jl`
- `Simplify/baseline.md`
- `Simplify/baseline/testexample/`
- `example/jellyroll_stress_displacement.jl`
- `src/Mechanical.jl`
- `src/CouplingState.jl`
- `src/CallModel.jl`
- `src/SetParams.jl`
- `src/parameters/Jellyroll.jl`

## 问题与处理

| 问题 | 处理 |
|---|---|
| `rg` 搜索中包含了不存在的 `param/` 目录并返回错误 | 已确认不影响其余目录结果；后续搜索只使用实际存在目录 |
| 搜索时误列不存在的 `src/NormaliseParam.jl`、`src/Case.jl` | 改用 `rg --files src` 定位实际文件后再核对 |
