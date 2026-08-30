# 堆芯塌陷力学建模：任务计划

## 目标

按 Theory v3.0 工况 C 契约补齐 JuBat 力学模块的非线性核心，最终实现堆芯塌陷（L4/C4）模拟能力。求解策略先采用 operator-splitting（Theory §6.6.5 允许），反馈闭环以同一时间步外层固定点实现。

## 关键建模决策（已由用户确认）

| 决策 | 内容 |
|---|---|
| 各向异性路线 | **不存在本征正交各向异性**：每个材料都是各向同性的，各向异性只由多种材料层叠加等效产生（切向并联/法向串联）；Theory 中的横观各向同性 5 常数表述属于理解错误，须修复理论文档（Batch 0'） |
| SP–电极界面 | 先按**完美接触（完美粘结）**假设：层内 SP 面共享节点天然粘结，跨匝 `phi_pairs` 绑定；单边接触+Coulomb 摩擦后置为 Batch 8 |
| 配置入口 | `czm_enabled` 为唯一主开关；几何非线性/J2 塑性/预应力/连续反馈默认关；spec v1.5 起 Φ 完美粘结并入 CZM 网格默认构造，不再有独立开关，结果必须携带近似声明 |
| 求解架构 | 先 operator-splitting；‖ΔD‖>0.05 等触发条件验证必要时再升级 monolithic |
| 新反馈 | 全部 opt-in；连续 R_contact(D)/A_eff(D) 与界面热阻恢复不改变默认基线 |
| 物理工况包络 | 用户确认 `Δsoc∈[-1,1]`、`|ΔT|≤20°C`；超出范围的载荷只作数值回归压力测试，不作为能力验收或物理阻塞项 |
| 文档 | 双份：`docs/superpowers/specs/`（规格）+ 本目录三件套；均需更新 |

## 当前阶段

Batch 5 实现与原收尾门禁已完成；2026-08-23 提交审查修复 R0–R6 complete。2026-08-24 已纠正 D10 的机械验证断点：夹具误固定 `Γ_in,free`、winding cohesive 法向反向，以及 geo 弧长 λ=0 预应力参考态只做一次 Newton/失败预测子步不回滚。按用户随后确认的 `Δsoc∈[-1,1]`、`|ΔT|≤20°C` 物理包络复算，16/16 边界/轴点全部收敛、`D=0`、最大起损比 0.438、最大 `Δ_core=0.002755`（约 5.29 μm）；`Δsoc=1.5/2.0` 的大损伤或极限点只保留为超范围压力测试，不再作为当前工况断点。专项与全量代码测试 34/34 通过。法向修正引起的冻结力学漂移随后已由用户接受，并纳入 `testexample` v4–v8 基线；该修复批次可提交。损伤影响有效电化学面积的 Batch 6A 尝试已按用户要求完整回滚，Batch 6/7 保持未开始。修复记录见 `../31_堆芯塌陷提交修复/task_plan.md`。

## 批次

### Batch 0'：Theory 各向异性表述修复
- [x] 修订 `Theory/03`：删除单层本征横观各向同性/正交各向异性 5 常数闭合，改为逐层各向同性 + 叠层等效（Voigt 并联/Reuss 串联）
- [x] 公式编号重排单调：(2.35) 逐层稳定性 → (2.36) 逐层各向同性矩阵 → (2.37) 等效参数集 → (2.38) 等效矩阵 → (2.39)/(2.40) 工况选择
- [x] 联动检查 Theory/02、04、05、06、07、08、09、10：均为叠层等效或通用符号，无需改动
- [x] grep 验证无本征各向异性残留；全套测试 26/27（唯一失败 `unit_czm_eigenstrain.jl` 为界面术语统一任务记录的既有缺口）
- **状态：** complete

### Batch 1：规格文档 + 基线冻结 + 非线性机械核边界
- [x] 规格文档已写入 `docs/superpowers/specs/2026-08-20-core-collapse-mechanics-design.md`，并经后续 v1.5 决策同步。
- [x] 开工前已登记 `test/` 现状；既有 eigenstrain 缺口在后续参数一致性修复后纳入全量绿色门禁。
- [x] 按 AGENTS 9.6 核对 `example/testexample.jl` 对 `Simplify/baseline/testexample/`（退出码/网格步数/metrics.toml/PNG SHA-256）。
- [x] `Option` 子选项与 `assemble_bulk_residual_tangent` 已实现；未实现/冲突组合保持显式失败。
- **验收门**：现有路径基线按记录精度一致；默认关闭的新增路径不改变 `testexample`；现行 standalone v3 快照不变。
- **状态：** implementation complete；现行门禁入口为 `tools/verify_czm_standalone.jl`，v3 与 `testexample` 已复核。

### Batch 2（C1）：von Kármán 残差 + 几何刚度 K_G
- 逐层各向同性 Q4（现有 `moduli_of` 路径不动）+ von Kármán 残差 + 一致 K_G；不叠加 CLT D
- 测试：有限差分切线、刚体零应变、patch、自由膨胀零应力、K_G→0 退化回冻结解
- **验收门（C1）**：塑性关、接触关、K_G→0 按记录精度退回冻结基线
- **状态：** implementation complete；完全 Green–Lagrange/几何刚度路径及回归锚通过。

### Batch 3（C2/L1）：PCC/NCC J2 塑性
- 参数集新增 `σ_y0`（Al 60 MPa / Cu 200 MPa）、`H→0+`（AGENTS 9.4 防御模式）；高斯点平面应力 J2 返回映射 + 一致切线 + `εp/κ` 提交—回滚；跨半循环持久化
- 测试：单轴屈服、卸载、KKT、塑性耗散非负、切线有限差分
- **验收门（C2/L1）**：强侧约束+完美圆仅预屈曲场；塑性关退回 C1
- **状态：** implementation complete；历史塑性应变只扣除一次，局部不收敛显式失败。

### Batch 4'（小批次）：Φ 跨匝完美粘结绑定
- Φ 完美粘结已并入 CZM 默认求解网格；`phi_pairs` 保留未合并力学网格真实节点对，`mesh_bonded`/`phi_keep` 单独承担自由度合并
- 测试：绑定对位移一致；基线不变
- **状态：** implementation complete；旧 5 参数 `CzmSubmesh` 构造器已恢复。

### Batch 5（C4-lite）：多圈状态 + 路径跟踪 + Δ_core
- 状态持久化挂 Case 层（仿 `czm_layout.u_prev`）：`u/几何`、`εp/κ`、cohesive D、`Δ_core`；**修改 `src/CycleSolver.jl` phase 交接**，参考构型来自持久化状态而非重建完美圆；失败步全回滚
- `solve_czm_arc_length_step` 扩展为全机械残差路径跟踪；分岔附近步长缩减与可诊断终止
- **验收门（C4-lite）**：完美粘结近似下多圈 Δ_core 与 D 联合增长；单循环全塌陷不算成功；默认基线一致
- **状态：** implementation complete / baseline accepted；BC、法向、弧长参考态与失败子步修复已通过专项及全量 34/34，物理包络 16/16 收敛且未起损。超范围高载停滞不再阻塞当前工况；法向引起的可解释漂移已进入 `testexample` v4–v8 基线。包络内仅出现微米级不圆度、无 `D` 联合增长，故仍不作塌陷能力宣称。

### Batch 6：损伤–电–热反馈闭环（opt-in）
- `czm_continuous_feedback=true`：`R_contact(D)/A_eff(D)` 进分支电压/BV 残差；旧面积权重代理保持默认
- 按 `cohesive_to_thermal` N:1 归约恢复界面热阻 + `Q_contact`；D=0 时热结果与现模型一致
- 同步外层固定点（Theory 07 策略 B）；不收敛报错、不提交部分状态
- **验收门**：D=0 退化一致；D↑ 方向正确
- **状态：** deferred / not implemented。2026-08-23 的 `A_eff(D)` 单机制代码与测试已完整回滚；理论先明确损伤对应的物理反应面积、正负极/界面作用域、物质与能量守恒以及与旧面积权重、完全失效和后续接触项的组合规则，再重新立项。

### Batch 7：系统验证与宣称门禁
- C1+C2 全量、C4-lite（带"完美粘结近似、SP 先皱机理未建模"声明）、耗散/守恒检查、网格/步长敏感性、Φ 完整性、运行时间冒烟预算
- 每批跑 `example/testexample.jl` 冻结门禁；任一指标变化即停止定位

### Batch 8（后置）：SP–涂层单边接触 + Coulomb 摩擦
- 替换 Batch 4' 完美粘结：层内 3 个 SP 邻接面 + 跨匝 `phi_pairs` 的 SP–PE 面；法向零抗拉罚/增广 Lagrange，切向 μ=0.10（0.05–0.4 敏感性）stick/slip
- **验收门**：C3+C5；解锁 C4 全量与"可模拟堆芯塌陷"完整宣称

## 约束

- 不修改 `md/`、`src/`、`test/` 于 Batch 0'（已完成：仅 Theory/03 文档修订）
- 规格文档未经用户评审不进入代码实现
- 每批遵守 AGENTS 9.6/9.7 基线门禁；不以"数值接近"放行
- 逐批提交到 `codex/src-physics-modularization`，沿用现有提交风格；每批更新本目录三件套与总索引
