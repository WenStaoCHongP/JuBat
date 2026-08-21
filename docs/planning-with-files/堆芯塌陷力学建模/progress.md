# 堆芯塌陷力学建模：进度记录

## 会话：2026-08-17

### 计划制定（brainstorming，architectural 路径）

- **状态：** complete
- 三个并行探索完成：界面术语统一任务契约、力学/CZM 代码现状、Theory 工况 C 验收要求（C1–C5/L1–L4）。
- 用户四项决策确认：全路线分批、czm 开关为主开关内部子选项、operator-splitting 先行、双文档。
- 计划自评审发现并修正：P0 本构缺参报错须限定工况 C 路径、连续反馈必须 opt-in、参数集扩展遗漏、CycleSolver 遗漏、记录精度判据。
- 用户两次修订：①各向异性=逐材料各向同性+叠加等效，修复理论错误；②SP–电极先完美接触，接触/摩擦后置 Batch 8。

### Batch 0'：Theory 各向异性表述修复

- **状态：** complete
- 修订仅涉及 `Theory/03_本构理论.md`（详见 `findings.md` 修订表）：§2.1 表格与注 2.1b、统一二维约化、注 2/注 2.3、§2.3.5 输出、§2.4.2 CLT Q̄、§2.4.4 整节重写（公式 2.35–2.40 重排单调）、§2.4.5、§2.6（含 ε_zz 各向同性闭合 (2.42a)/(2.47)）、§2.6.4 输出与自查门槛。
- 联动核查 Theory/00/02/04/05/06/07/08/09/10/CLAUDE：均无需改动。
- 建立本目录三件套并更新 `docs/planning-with-files/index.md`（29 目录、93 文件）。

## 会话：2026-08-20

### 理论验证 + FEM 可行性评审 + spec v1.1 修订

- **状态：** complete（spec 待用户评审）
- 通读 Theory 00/02/03/04/07/08 工况 C 相关段落 + 代码锚点核验（Czm/CzmSolve/CouplingState/Jellyrollmodel）。
- 发现 5 项理论缺陷（Δ_core 基线污染、(1.30a) 坐标系错配、R-EC-1 疑似双计、Q_contact j 语义、J2 表述），详见 findings。
- FEM 可行性确认：全部标准算法，风险在网格长宽比与 C4-lite 可达性，非算法本身。
- 用户 4 项决策（D8–D11）确认并写入 spec：位移基 Δ_core、完全 GL TL、探针后停止、预应力 opt-in Batch 2'。
- spec 更新：头部 v1.1 行、决策表 D8–D11、§3.2/3.3/3.4/3.5/3.6 修订、新增 §3.7、§4.1/§5/§7/§8/§9 联动、§10.1 评审记录；自评审补齐 §1↔§8 SP 热膨胀交叉引用。

### 第二轮独立复核（全 13 文件 + 代码交叉核验）

- **状态：** complete（待用户就 4 项存疑决策商讨）
- 4 个并行子代理分片复核（01/03、05/06、07/08、09-11），逐条重读原文与代码独立验证；子代理若干条经核查为表述不准或与独立判断不符，未采纳。
- **更正第一轮两处错误结论**：`Q_contact` 中 `j` 是界面电流密度（A/m²）而非分支总电流，且体积源需乘 `a_s`；`t_repeat` 代码实参为 373.6 μm 而非理论表的 284 μm。
- 新增确定发现 E1–E13，其中影响决策的三项：E1 理论几何参数表与代码不符（L_turn 3.5×、L_spiral 7.8×、匝数 2.2×，疑因把直径当半径）；E2 `(1−D)⁻²` 由「疑似双计」升级为**确定内部矛盾**（§3.6.1 几何推导 vs §4.5.2 材料劣化辩护互斥，支撑文献不可核查）；E13 每层厚度方向仅 1 个 Q4、薄层长宽比 29–38:1 ⟹ 单层起皱不可解、叠层椭圆化可解。
- 用户 4 项决策 D12–D15 确认并写入 spec v1.2：`(1−D)⁻²` 延后至 Batch 6 专项推导、网格探针后定能力边界、几何表全库重算 + 字段溯源、理论修订拆 Batch 0''/Batch 7。
- spec v1.2 更新：头部版本行、决策表 D12–D15、§3.5 数值更正、§3.6 `j`/`a_s` 更正与 D12 三候选、新增 §3.8 网格能力边界与探针、§4.1 新增探针脚本与 Theory 行、§7 新增 Batch 0''/2'' 行、§8 风险表四项重写、§9 批次顺序插入 0''/2'' 并列出 8+5 项清单、§10.2 复核记录。

## 验证记录

| 检查 | 预期 | 当前结果 | 状态 |
|---|---|---|---|
| grep `横观各向同性` | 仅注 2.1b 废除声明 | 符合 | 通过 |
| grep `正交各向异性` | 仅叠层等效矩阵形式描述 | 符合 | 通过 |
| grep `E_s^i/G_{sn}^i/ν_{sz}/ν_{nz}` | 0 | 0 | 通过 |
| 公式 tag 顺序 | 2.35→2.40 单调 | 符合 | 通过 |
| 尾随空白 | 0 | 0 | 通过 |
| 全套测试 | 与既有状态一致 | 26 通过/1 失败（`unit_czm_eigenstrain.jl` 既有缺口） | 通过 |
| Git 变更范围 | 仅 Theory/03 + 本目录 + 索引 | 符合 | 通过 |

## 错误日志

| 时间 | 错误 | 尝试 | 处理 |
|---|---:|---:|---|
| 2026-08-17 | ugrep 对混合中文+LaTeX 花括号模式报 invalid character class | 2 | 拆分简单模式分别检索 |
| 2026-08-17 | Git Bash 下 sed/grep 输出中文文件名乱码 | 多次 | 行号仍有效；内容核对改用 Read 工具 |

## 会话：2026-08-21

### 实现计划评审 + spec v1.3 / 计划 v1.1 修订（方案 B）

- **状态：** complete（Batch 0'' 待用户放行）
- 评审 Batch 0''+1 实现计划（`49ec452`）并对 HEAD `e117fd2` 逐项核实：3 处前提失效（探针被删/PNG 路径迁移/eigenstrain 已修复）+ 4 项次要缺口，详见 findings"实现计划评审"节；核对无误项（行号、签名、几何数值、基线冻结值）未改动。
- 用户决策**方案 B**：Batch 1 基线门禁改用 `tools/verify_czm_standalone.jl`；spec v1.3（§头部/§5/§7/§10.3）与计划 v1.1（Task 1 重写为复核、Task 2 KKT 扩清理、Task 3 计数订正、测试预期 22/22 全绿、PNG 路径更正、偏差登记第三项、修订记录 8 条）同步修订。
- **基线已冻结**：实测运行 `verify_czm_standalone.jl`（exit 0），快照写入 `baseline_czm_standalone.md`——网格 10946 节点/6728 bulk/3364 cohesive；8×3 收敛表（0.1–5.0 全 OK，10.0 三方法 FAIL 冻结）；Summary 三方法各 7/8、total_iter 16/794/2187。工具既有瑕疵两条登记 findings，不修改。
- 待提交变更：spec v1.3、计划 v1.1、`baseline_czm_standalone.md`、findings/progress/index（提交须经用户授权绕过 Mimosa 守卫）。

### 会话：2026-08-21（续）：Batch 0''+1 计划执行

#### Task 1：基线复核

- **状态：** complete
- 重跑 `verify_czm_standalone.jl`（HEAD `40757d7`，与冻结同环境），全部数据行与冻结快照**逐位一致**（含 10.0 水平 FAIL 与 stall 警告）——基线可复现，Batch 1 门禁有效。
- findings 的工具瑕疵登记在位（`:66`/`:134` 传参不一致、未合并 thermal2D）。

#### Task 2–4：Batch 0'' 理论修订（8 项阻塞矛盾）

- **状态：** complete
- ① `07` 同号 §6.4.5 去重（保留位移阈值代数更新版）+ `:573` KKT 引用清除；执行中另清理 5 处现行 KKT 表述（§1.1 表 ×2、§6.9 输入 3、§6.9.6 校核表、本章小结）。
- ② `K_uu` 统一三项定义（式 6.8 补 cohesive 项、`K_G(σ)→K_G(S)`；6.54 改复述）；`:1073` 引用错位转 Batch 7。
- ③ 新增 §6.10 柱面弧长法（Crisfield 式 6.88–6.91，λ 缩放本征应变增量）。
- ④ §2.6 后补 Φ 跨匝配对施加方式（节点合并/主从消元，非弱形式项）。
- ⑤ `02` 式 (1.30a) 后加 κ_ss 与 C⁰ Q4 不兼容实现注记（D9，含能力边界：单层起皱不可解）。
- ⑥ `01` 插入层编号双约定声明（5 材料类型 vs 8 物理层序）。
- ⑦ `A_eff` 全库统一无量纲（计划 3 处 + 执行发现 7 处；`09:62` R_contact 单位 Ω→Ω·m²；D12 冻结区未触碰）。
- ⑧ 几何参数全库按代码重算（21.03 匝 / 0.7727 m / 373.6 μm / γ∈[0.3485°,1.7744°]，02 表新增"代码字段"列；计划外同源值 6 组一并重算，见 findings）。
- **门禁**：`rg "46.6|0.132"` 0 命中；`rg "284"` 仅 `06:69` ρ≈2846 豁免；派生比值 0 残留；层编号声明/A_eff/κ_ss 三 grep 通过。
- **全套测试**：`test/runtests.jl` 22/22 通过（3m37.6s，含 eigenstrain 60/60），文档批次零回归。
- 提交：Task 2 `d28b3b4`、Task 3 `879440c`、Task 4（本批）。

#### Task 5–7：Batch 1（Option 子选项 + bulk 残差/切线统一入口 + 接线）

- **状态：** complete
- Task 5（`cd7be6e`）：`Option` 新增 6 个默认关子选项（`czm_geo_nonlinear`/`czm_winding_prestress`/`czm_j2_plasticity`/`czm_phi_bond`/`czm_continuous_feedback`/`czm_friction_mu`），TDD 三 testset 锁定默认值契约。
- Task 6（`016a8b9`）：`src/czm.jl` 新增 `assemble_bulk_residual_tangent`（线弹性槽位：`f=K_bulk*u`、`K_tan=K_bulk`；`geo_nl`/`plasticity`/`mech_state` 非默认即 `error`；`K_bulk_cached` 透传同一对象）。6/6 testset 通过，含与 `K_bulk*u` 逐位等价、缓存同一对象（`===`）、对称性、维度/槽位报错。
- Task 7：`assemble_coupled_system` 改走新入口（`assemble_coupled_system_full` 间接经它）。与 spec §4.1/§4.2 偏差三项已登记（Option 提前全加、`PlasticState`/`MechHistory` 延后、`K_bulk_cached` 签名扩充）。
- **三道门禁实测**：全套 `runtests.jl` **24/24**（4m09s，含两个新测试文件）；`verify_czm_standalone.jl` 快照与冻结基线**逐位一致**（8×3 表含 10.0 FAIL 条目与 Summary）；`testexample.jl` exit 0，网格 1682/1763、步数 19、电压 4.0367/3.9438/0.0929、容量 0.0833 Ah、温度 298.15–299.00 K、D 0.0000%、分离 1.2557e-14 m、断裂 0，全部与冻结表一致，PNG SHA-256 `4ba6207c…e932` 一致。
- **A/B 复核**：冻结表"CZM converged updates 19/19"不对应脚本打印行（日志实为 18 条 `converged=true` debug 行）；经 stash A/B 确认**接线前后同为 18 行、同一 PNG SHA**，属登记口径而非本批漂移；已在本节记录，不改冻结表（该行非脚本打印指标）。
- `Simplify/baseline.md` 批次表已追加 Batch 1 行（PASS）。

### 会话：2026-08-21（续二）：Batch 2（K_G/C1）执行

- **状态：** complete
- 计划评审通过（用户批准 D-B2-1/D-B2-2，计划 `2026-08-21-core-collapse-batch2-plan.md`，提交 `a4b58a4`）。
- Task 1（`7cc8b0d`）：`gl_element_residual_tangent`（完全 GL、SVK、标准初应力 K_G）+ `geo_nl` 槽位启用 + `eigenstrain` 关键字；5/5 testset（FD 切线、30° 刚体转动机器零内力、u=0 切线≡线性刚度 rtol 1e-12、自由膨胀零应力、K_G 压缩方向性）。
- Task 2（`aea1e85`）：mech_core 夹具改传 `case.param`（归一化网格对齐生产）；Batch 1 的"geo_nl 未实现须报错"断言随槽位实现演进为"缓存冲突契约"断言。
- Task 3（`b551dac`）：`geo_nl`/`eigenstrain` 贯穿 `assemble_coupled_system` → basic/load_substep（含 backtrack 与 newton 线搜索）→ `update_czm_damage!` 生产调用点；`arc_length+geo_nl` 显式 error（Batch 5）；C1 2/2（K_G→0 极限 geo≈线性 rtol 1e-4；geo off 与 Batch 1 逐位一致）。
- **门禁执行中发现并修复**：`b551dac` 的 4 空格缩进 replace_all 以子串方式误伤嵌套 arc 函数内同文本行（`eig_kwargs` 未定义 → 探针 arc 方法 exit 1）——**Gate B 门禁捕获**，回滚 `1bb0026`，重跑探针与冻结基线逐位一致。教训：对"缩进即作用域"的嵌套函数禁用短前缀 replace_all。
- **三道门禁全绿**：全套 **26/26**（4m21s）；`verify_czm_standalone` 快照逐位一致；testexample 全部冻结指标与 PNG SHA `4ba6207c…e932` 一致（geo_nl 默认关零漂移）。
- 执行偏差（计划小瑕疵，现场修正）：① FD 测试 dof=5000 超出 nθ=8 网格自由度 → 改动态取；② isapprox 自定义 norm 单参调用 MethodError → 改稠密比较；③ Batch 1 geo_nl 报错断言过期（上文 Task 2）。

### 会话：2026-08-21（续三）：Batch 3（J2 塑性/C2）执行

- **状态：** complete
- 计划 `2026-08-21-core-collapse-batch3-plan.md`（`b52c4ba`，用户指示直接执行）；参数联网核实（findings"Batch 3 参数核实"节：Al 15μm 硬态 UTS 150–290 MPa/屈服≈0.9UTS、ED Cu 108–441 MPa、薄硬态箔硬化近耗尽——spec 缺省 60/200/0 落点合理）。
- Task 1（`4fcb243`）：`src/CzmPlasticity.jl`——`PlasticState`、`foil_params_of`、4×4 Newton 平面应力 J2 一致返回映射（切线由收敛系统线性化精确导出）、物理箔参数（PCC 70GPa/0.33/60MPa、NCC 110GPa/0.34/200MPa、H=0，σ_czm 归一）；`unit_czm_j2.jl` 5/5（弹性不变/屈服面停留/卸载永久应变/一致切线 FD/组合回面）。执行修正：`J[4,1:3].=n'` 广播维错、`[I zeros]` 拼接错、单轴 σ11=σ_y 断言在平面应力关联流动（n=[1,−0.5,0]）下不成立改为不变量断言。
- Task 2（`495287c`）：GL 单元 plastic/commit_to、装配 mech_state 位置参数、求解器全链路、`CzmLayout.plastic_states` 持久化、收敛后 `commit_plastic!`；集成 4/4（C2-lite：Δsoc=−0.3 失配驱动下 PCC/NCC κ>0 其余层 0、回归锚、不可逆有界、非法组合报错）。执行修正三处：mech_state 误作关键字、提交元组缺嵌套、**fill(PlasticState(),ne,4) 可变对象别名共享**（经典陷阱，测试与生产初始化都已改独立构造）。
- **语义说明（D-B3-2 推论）**：同载荷重解存在有界的再预测漂移（κ 二次增长 ~1.6×），单调不可逆性保持；生产时序每步新 Δε* 下漂移为二阶小量。
- **三道门禁全绿**：全套 28/28（5m07s）；探针快照逐位；testexample 冻结指标与 PNG SHA 一致（塑性默认关零漂移）。

## 5 问重启检查

| 问题 | 回答 |
|---|---|
| 当前在哪里？ | Batch 3（J2/C2）完成：塑性全链路 + 持久化，三道门禁全绿（28/28） |
| 下一步去哪里？ | Batch 4'（Φ 完美粘结小批次）或 Batch 2'（预应力）立项出计划，经用户评审后执行 |
| 目标是什么？ | 补齐工况 C 力学非线性核心，最终实现堆芯塌陷模拟（L4/C4） |
| 已学到什么？ | 见 `findings.md` |
| 已做什么？ | 见本文件上述记录 |
