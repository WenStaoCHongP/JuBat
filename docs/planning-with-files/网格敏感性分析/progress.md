# Progress: 网格敏感性分析

## Session Log

### 2026-04-22 计划落盘

**已完成：**

- 创建了三份规划文件：`task_plan.md`、`findings.md`、`progress.md`
- 记录了电化学、热学和内聚力三类网格的分析范围和指标
- 记录了热 / 内聚网格按 cohesive 特征长度换算 `nθ` 区间的规则
- 明确了 lumped 热模型不直接提供空间温度梯度这一约束

**当前状态：**

- 仍处于计划阶段
- 尚未运行任何网格敏感性计算
- 尚未生成误差表或收敛图

### 2026-04-22 口径再确认

**已确认：**

- 热网格独立判据采用热穿透深度，$\delta_T \sim \sqrt{\alpha t_{\max}}$
- 内聚力特征长度采用 $l_c = G_c \cdot E / \sigma_{\max}^2$，其中 $E$ 取体材料杨氏模量
- CZM 指标同时保留载荷-位移曲线和界面牵引-分离曲线
- 热模型采用固定热源的热传导基准

**下一步：**

1. 确认 cohesive 特征长度的最终定义
2. 把热 / 内聚网格的 4 个角度划分点定下来
3. 按三个分析方向分别跑参考网格与粗网格对比
4. 将每轮结果回写到 `findings.md`

### 2026-04-22 评审与修订（rev 2）

**评审发现的关键问题（已全部解决）：**

1. **热穿透深度 $\delta_T$ 不适用**：$\delta_T \approx 194\ \text{mm} \gg R_{out} = 10.15\ \text{mm}$，无法导出有意义的 nθ 下界。替换为 Biot 数分析 + 系统加密。
2. **$E_{eff}$ 已确认**：`Mechanical.jl:181` 的厚度加权平均 $E_{eff} \approx 25.6\ \text{GPa}$，由此 $l_c \approx 96\ \mu\text{m}$。
3. **纯热路径已确认**：`example/热模块验证/thermal_verify.jl` 提供了 `opt.model = "thermal"` + 均匀热源的现成模式。
4. **载荷-位移曲线已澄清**：纯机械模型施加位移边界条件触发损伤。
5. **nθ 定义已确认**：每周的周向分辨率。

**已更新的文件：**

- `findings.md` §3: 完整的 Biot 数推导、$l_c$ 计算过程、nθ 定义
- `findings.md` §6.5: 口径确认更新
- `plans/2026-04-22-grid-sensitivity-analysis-plan.md`: Thermal Track 和 Cohesive Track 全面修订
- `specs/2026-04-22-grid-sensitivity-analysis-design.md`: §2.1-2.3 和 §3.2 更新

**当前状态：**

- 热网格 nθ 候选值已确定：$\{20, 40, 80, 160\}$
- CZM nθ 区间已确定：126–664（待取整为 4 个值）
- 尚未运行任何网格敏感性计算

**下一步：**

1. 将 CZM nθ 区间 126–664 三等分为 4 个具体值
2. 按三个分析方向分别跑参考网格与粗网格对比
3. 将每轮结果回写到 `findings.md`

### 2026-04-22 追加能量守恒检查

**新增内容：**

- 在 findings.md §7 中定义了全场耦合能量守恒方程
- 在 plan §4 中新增 Energy Conservation Check 章节
- 能量项：电功 $W_{elec}$、热能变化 $\Delta E_{th}$、弹性应变能 $\Delta E_{el}$、断裂能 $E_{frac}$、边界热损失 $Q_{loss}$、化学能变化 $\Delta E_{chem}$
- 简化备用方案：热-弹-断裂子系统检查（不依赖化学能定义）
- 验收标准：$\epsilon_R < 1\%$，无系统性增长趋势

**代码侧状态：**

- ✅ 电功：`cell voltage [V]` × `cell current [A]` 直接可用
- ⚠️ 热能：需新增热容矩阵加权积分
- ⚠️ 弹性应变能：`assemble_bulk_stiffness` 已有，需加 $\frac{1}{2}u^Ku$ 积分
- ⚠️ 断裂能：`damage_states` 和单元长度已有，需加 $G_c \cdot l_e \cdot D_e$ 求和
- ⚠️ 边界热损失：需提取边界节点温度和面积
- ⚠️ 化学能：`element soc_n/soc_p` 和 OCV 函数已有，需组装积分
- ✅ `czm traction normal/tangent` 和 `czm separation normal/tangent` 确认存在（`PostProcessing.jl:77-80`）

**当前状态：**

- 能量守恒方程已定义，实现路径已明确
- 需要新增后处理函数来计算各能量项

### 2026-04-22 脚本编写完成

**已创建 5 个脚本**（位于 `example/网格敏感性/`）：

| 脚本 | 文件 | 功能 |
|------|------|------|
| 1 | `1_cohesive_characteristic_length.jl` | 计算 E_eff、l_c、CZM nθ 区间、热 nθ 候选 |
| 2 | `2_electrochemical_mesh_sensitivity.jl` | 4 组 (Nn,Ns,Np) SPMe+lumped thermal 收敛分析 |
| 3 | `3_thermal_mesh_sensitivity.jl` | 4 组 nθ={20,40,80,160} Jellyroll 电-热耦合收敛分析 |
| 4 | `4_czm_mesh_sensitivity.jl` | 4 组 CZM nθ（由 l_c 自动确定）全耦合收敛分析 |
| 5 | `5_energy_conservation_check.jl` | 简化方案能量守恒检查：Q_gen = ΔE_th + Q_loss + E_frac + R |

**图片输出目录**：`output/mesh_sensitivity/`

**当前状态：**

- 5 个脚本全部编写完成
- 尚未运行
- 下一步：按顺序执行脚本，回写运行结果到 findings.md
