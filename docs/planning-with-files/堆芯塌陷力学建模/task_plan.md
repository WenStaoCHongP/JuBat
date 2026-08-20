# 堆芯塌陷力学建模：任务计划

## 目标

按 Theory v3.0 工况 C 契约补齐 JuBat 力学模块的非线性核心，最终实现堆芯塌陷（L4/C4）模拟能力。求解策略先采用 operator-splitting（Theory §6.6.5 允许），反馈闭环以同一时间步外层固定点实现。

## 关键建模决策（已由用户确认）

| 决策 | 内容 |
|---|---|
| 各向异性路线 | **不存在本征正交各向异性**：每个材料都是各向同性的，各向异性只由多种材料层叠加等效产生（切向并联/法向串联）；Theory 中的横观各向同性 5 常数表述属于理解错误，须修复理论文档（Batch 0'） |
| SP–电极界面 | 先按**完美接触（完美粘结）**假设：层内 SP 面共享节点天然粘结，跨匝 `phi_pairs` 绑定；单边接触+Coulomb 摩擦后置为 Batch 8 |
| 配置入口 | `czm_enabled` 为唯一主开关；新增子选项（几何非线性/J2 塑性/Φ 绑定/连续反馈）全部默认关，`czm_enabled=true` 单独使用行为与现状完全一致 |
| 求解架构 | 先 operator-splitting；‖ΔD‖>0.05 等触发条件验证必要时再升级 monolithic |
| 新反馈 | 全部 opt-in；连续 R_contact(D)/A_eff(D) 与界面热阻恢复不改变默认基线 |
| 文档 | 双份：`docs/superpowers/specs/`（规格）+ 本目录三件套；均需更新 |

## 当前阶段

Batch 0' Complete — Batch 1 Pending（等待规格文档评审）

## 批次

### Batch 0'：Theory 各向异性表述修复
- [x] 修订 `Theory/03`：删除单层本征横观各向同性/正交各向异性 5 常数闭合，改为逐层各向同性 + 叠层等效（Voigt 并联/Reuss 串联）
- [x] 公式编号重排单调：(2.35) 逐层稳定性 → (2.36) 逐层各向同性矩阵 → (2.37) 等效参数集 → (2.38) 等效矩阵 → (2.39)/(2.40) 工况选择
- [x] 联动检查 Theory/02、04、05、06、07、08、09、10：均为叠层等效或通用符号，无需改动
- [x] grep 验证无本征各向异性残留；全套测试 26/27（唯一失败 `unit_czm_eigenstrain.jl` 为界面术语统一任务记录的既有缺口）
- **状态：** complete

### Batch 1：规格文档 + 基线冻结 + 非线性机械核边界
- [ ] 写规格文档（`docs/superpowers/specs/2026-08-17-core-collapse-mechanics-design.md`），含 Δ_core 离散定义（初定：内圈节点半径对最优拟合圆的最大偏差）；**须经用户评审后才开始实现**
- [ ] 开工前登记 `test/` 现状，把 `unit_czm_eigenstrain.jl` 既有失败标记为前置状态
- [ ] 按 AGENTS 9.6 核对 `example/testexample.jl` 对 `Simplify/baseline/testexample/`（退出码/网格步数/metrics.toml/PNG SHA-256）
- [ ] `src/Option.jl` 新增默认关子选项；`src/Czm.jl` 新增 `assemble_bulk_residual_tangent`：bulk/cohesive/contact 三槽位分列（contact 槽位留空）；位移/塑性/接触依赖切线不再沿用常量 `K_bulk` 缓存；弹性层材料切线可按单元缓存
- **验收门**：现有路径基线按记录精度一致；全开关关时新入口 ≡ `K_bulk*u`；`tools/czm_baseline_probe.jl` 三方法快照不变

### Batch 2（C1）：von Kármán 残差 + 几何刚度 K_G
- 逐层各向同性 Q4（现有 `moduli_of` 路径不动）+ von Kármán 残差 + 一致 K_G；不叠加 CLT D
- 测试：有限差分切线、刚体零应变、patch、自由膨胀零应力、K_G→0 退化回冻结解
- **验收门（C1）**：塑性关、接触关、K_G→0 按记录精度退回冻结基线

### Batch 3（C2/L1）：PCC/NCC J2 塑性
- 参数集新增 `σ_y0`（Al 60 MPa / Cu 200 MPa）、`H→0+`（AGENTS 9.4 防御模式）；高斯点平面应力 J2 返回映射 + 一致切线 + `εp/κ` 提交—回滚；跨半循环持久化
- 测试：单轴屈服、卸载、KKT、塑性耗散非负、切线有限差分
- **验收门（C2/L1）**：强侧约束+完美圆仅预屈曲场；塑性关退回 C1

### Batch 4'（小批次）：Φ 跨匝完美粘结绑定
- 工况 C 开关下将 `phi_pairs` 节点对绑定（合并重合节点或 MPC tie）：零穿透、零分离、零滑移；默认关不影响基线
- 测试：绑定对位移一致；基线不变

### Batch 5（C4-lite）：多圈状态 + 路径跟踪 + Δ_core
- 状态持久化挂 Case 层（仿 `czm_layout.u_prev`）：`u/几何`、`εp/κ`、cohesive D、`Δ_core`；**修改 `src/CycleSolver.jl` phase 交接**，参考构型来自持久化状态而非重建完美圆；失败步全回滚
- `solve_czm_arc_length_step` 扩展为全机械残差路径跟踪；分岔附近步长缩减与可诊断终止
- **验收门（C4-lite）**：完美粘结近似下多圈 Δ_core 与 D 联合增长；单循环全塌陷不算成功；默认基线一致

### Batch 6：损伤–电–热反馈闭环（opt-in）
- `czm_continuous_feedback=true`：`R_contact(D)/A_eff(D)` 进分支电压/BV 残差；旧面积权重代理保持默认
- 按 `cohesive_to_thermal` N:1 归约恢复界面热阻 + `Q_contact`；D=0 时热结果与现模型一致
- 同步外层固定点（Theory 07 策略 B）；不收敛报错、不提交部分状态
- **验收门**：D=0 退化一致；D↑ 方向正确

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
