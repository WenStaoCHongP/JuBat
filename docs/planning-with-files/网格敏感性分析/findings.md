# 网格敏感性分析发现记录

**创建日期**: 2026-04-22

## 1. 任务约束

- 电化学网格候选值固定为 `(40, 20, 40)`、`(20, 10, 20)`、`(20, 5, 20)`、`(10, 5, 10)`。
- 热网格由 Biot 数分析 + 系统加密确定候选值 `nθ = {20, 40, 80, 160}`（每周周向分辨率），验证收敛性。
- 内聚力网格按 cohesive 特征长度 $l_c \approx 96\ \mu\text{m}$ 换算为 `nθ` 区间，再三等分得到 4 个点。
- 目标判断标准统一为：各候选网格相对其最精细网格的偏差小于 5%。

## 2. 代码侧可用输出

### 2.1 电化学 / 热输出

`PostProcessing` 已经能输出这些关键量：

- `cell voltage [V]`
- `temperature [K]`
- `thermal2D temperature at nodes [K]`
- `thermal2D Q_* [W/m3]`

这意味着：

- 电压曲线可以直接做网格对比
- 温度峰值可以直接从时间序列或节点场中提取
- 对 lumped 模型，梯度指标不取空间梯度，而是温度随时间的梯度 dT/dt，需要在后处理里从温度-时间序列额外计算

### 2.2 机械 / CZM 输出

机械与内聚力侧已经有可用指标：

- `diffusion stress vonMises`
- `thermal stress vonMises`
- `displacement x`
- `displacement y`
- `czm D_max`
- `czm D_mean`
- `czm n_fractured`
- `czm traction normal`
- `czm traction tangent`
- `czm separation normal`
- `czm separation tangent`

这些量足够支持：

- 应力峰值
- 载荷-位移曲线
- 起裂时刻
- 裂纹拓展速率

## 3. 热 / 内聚网格区间的换算口径

### 3.1 热网格——Biot 数分析 + 系统加密

**为何不用热穿透深度 $\delta_T$：**

各向异性导热参数（代码位置 `Jellyroll.jl:172-173`）：

$$\lambda_r = \frac{\sum d_i}{\sum d_i / \lambda_i} \approx 1.32\ \text{W/m/K} \quad\text{(径向/through-plane, 串联)}$$

$$\lambda_t = \frac{\sum d_i \lambda_i}{\sum d_i} \approx 25.3\ \text{W/m/K} \quad\text{(切向/in-plane, 并联)}$$

各向异性比 $\lambda_t / \lambda_r \approx 19\times$。

等效热扩散率：

$$\alpha_t = \frac{\lambda_t}{\rho c_p} \approx \frac{25.3}{2813 \times 860} \approx 1.04 \times 10^{-5}\ \text{m}^2/\text{s}$$

热穿透深度（$t_{\max} = 3600$ s）：

$$\delta_T = \sqrt{\alpha_t \times 3600} \approx 194\ \text{mm} \gg R_{out} = 10.15\ \text{mm}$$

$\delta_T$ 远大于域尺寸，无法从中推导出有意义的最小 nθ 下界。

**Biot 数分析（替代方案）：**

切向 Biot 数：

$$Bi_t = \frac{h \cdot R_{out}}{\lambda_t} = \frac{10 \times 0.01015}{25.3} \approx 0.004 \ll 0.1$$

$Bi_t = 0.004$ 说明**温度场在 $\theta$ 方向几乎轴对称**。因此：

1. 热精度对 nθ 极不敏感——即使 nθ = 20，热结果也可能与 nθ = 200 差异 < 1%
2. 真正需要 nθ 的原因是**热源空间分辨率**（per_element_spme 时），而非热传导本身
3. 径向温度梯度由螺旋几何的圈数决定（~22 圈），不由 nθ 决定

**最终方案：直接给出经验候选值，验证收敛**

$$n\theta \in \{20, 40, 80, 160\}$$

这些值跨越了 8 倍的分辨率范围，足以观察收敛趋势。运行时监控 T_max 和温度空间均匀性。

### 3.2 内聚力网格——特征长度换算

Cohesive 特征长度 $l_c = G_c \cdot E_{eff} / \sigma_{\max}^2$ 中各参数已确认：

| 参数 | 值 | 代码位置 |
|------|-----|---------|
| $G_c$ | 25.3 J/m² | `Jellyroll.jl` Cohesive 法向断裂能 |
| $\sigma_{\max}$ | 82 MPa | `Jellyroll.jl` Cohesive 法向峰值强度 |
| PE.E | 37.5 GPa | `Jellyroll.jl:28` |
| NE.E | 15.0 GPa | `Jellyroll.jl:56` |
| $E_{eff}$ | $\approx$ 25.6 GPa | `Mechanical.jl:181` 厚度加权平均 |

$$E_{eff} = \frac{NE.E \cdot t_{NE} + PE.E \cdot t_{PE}}{t_{NE} + t_{PE}} = \frac{15 \times 85.2 + 37.5 \times 75.6}{85.2 + 75.6} \approx 25.6\ \text{GPa}$$

$$l_c = \frac{25.3 \times 25.6 \times 10^9}{(82 \times 10^6)^2} \approx 9.6 \times 10^{-5}\ \text{m} \approx 96\ \mu\text{m}$$

换算为 nθ 区间（要求单元弧长 < $l_c$）：

- 外圈约束：$n\theta > 2\pi R_{out} / l_c \approx 2\pi \times 10.15 / 0.096 \approx 664$
- 内圈约束：$n\theta > 2\pi R_{in} / l_c \approx 2\pi \times 1.92 / 0.096 \approx 126$

将该区间三等分可得 4 个 nθ 水平（需结合实际计算成本取整）。

### 3.3 nθ 定义确认

`nθ` 是**每周的周向分辨率**（per-revolution circumferential resolution），由 `jellyroll_collector_seed_mesh(param; nθ=ntheta)` 控制。参考 `thermal_verify.jl:80` 中的 ring mesh 定义。

总单元数 $\approx n\theta \times N_{turns}$，其中 $N_{turns} \approx 22$。

### 3.4 待确认项

- CZM 的 nθ 区间取整后是否需要强制保持严格递增
- 热网格与内聚力网格是否需要最终落到同一组几何节点（建议各自独立即可）

## 4. 指标定义建议

- 电压变化：用全时域最大相对偏差，或在容量坐标上重采样后的最大相对偏差
- 温度峰值：用 `max(T(t))` 对比参考解
- 最大温度梯度：对 lumped 模型取温度随时间的梯度峰值 dT/dt；若是分布式热场，则另行按空间梯度定义
- 应力峰值：用 Von Mises 峰值对比参考解
- 载荷-位移曲线：保留常规定义的全局响应曲线，用最大相对偏差或曲线面积差对比参考解
- 界面牵引-分离曲线：保留 traction-separation 响应曲线，用最大相对偏差或曲线面积差对比参考解
- 起裂时刻：用首次出现失稳或 `D_max` 达到阈值的时刻
- 裂纹拓展速率：用 `n_fractured(t)` 的增长斜率或前后差分

## 5. 风险点

- lumped 模型不直接给出空间梯度，必须提前说明梯度指标的来源；这里的梯度指标统一指温度随时间的梯度 dT/dt
- 若参考值接近 0，相对误差会失真，需要改成绝对误差
- 载荷-位移曲线和起裂时刻可能受步长控制影响，建议使用统一的时间步策略

## 6.5 讨论后的口径确认

- 热网格独立判据：Biot 数分析表明温度场近乎轴对称（$Bi_t \approx 0.004$），直接用经验候选值 $n\theta \in \{20, 40, 80, 160\}$，验证收敛即可
- 内聚力特征长度：$l_c = G_c \cdot E_{eff} / \sigma_{\max}^2 \approx 96\ \mu\text{m}$，其中 $E_{eff} \approx 25.6\ \text{GPa}$ 取自 `Mechanical.jl:181` 的厚度加权平均
- 纯热基准：采用均匀体积热源，参考 `example/热模块验证/thermal_verify.jl` 的实现模式（`opt.model = "thermal"`，`q_func = (r,theta,t) -> q0`）
- CZM 载荷-位移曲线：纯机械模型施加位移边界条件触发损伤演化，不依赖热-化学耦合
- nθ 定义：每周的周向分辨率

## 7. 全场耦合能量守恒检查

### 7.1 能量平衡方程

对全耦合系统（电化学 + 热 + 弹性 + CZM），基于热力学第一定律：

$$\frac{dE_{total}}{dt} = P_{elec}(t) - Q_{loss}(t)$$

其中 $E_{total} = E_{thermal} + E_{elastic} + E_{chem}$ 为系统总储能，$P_{elec} = V \cdot I$ 为电功率输入，$Q_{loss}$ 为边界热损失。

展开各项后，残余量定义为：

$$R(t) = W_{elec}(t) - Q_{loss}(t) - \Delta E_{thermal}(t) - \Delta E_{elastic}(t) - E_{fracture}(t) - \Delta E_{chem}(t)$$

相对误差：

$$\epsilon_R(t) = \frac{|R(t)|}{|W_{elec}(t)|}$$

对于能量守恒的离散格式，$\epsilon_R$ 应随网格加密趋于零。

### 7.2 各项能量定义与计算方式

| 能量项 | 定义 | 计算方式 | 代码可用性 |
|--------|------|----------|-----------|
| 电功 $W_{elec}$ | $\int_0^t V(\tau) \cdot I(\tau) \, d\tau$ | 对 `cell voltage [V]` 和 `cell current [A]` 时间序列做梯形积分 | ✅ 直接可用 |
| 热能变化 $\Delta E_{th}$ | $\int_V \rho c_p [T(t) - T_0] \, dV$ | 对节点温度场做 $\sum_n M_n [T_n(t) - T_n(0)]$（热容矩阵加权） | ⚠️ 需新增：提取热容矩阵 $M$ 和节点温度差 |
| 弹性应变能 $\Delta E_{el}$ | $\frac{1}{2} u^T K_{bulk} u \Big\vert_0^t$ | 从 CZM 位移场和体刚度矩阵计算 | ⚠️ 需新增：`assemble_bulk_stiffness` 已有 (`czm.jl:329`)，需加积分 |
| 断裂能 $E_{frac}$ | $\sum_{e \in \text{cohesive}} G_c \cdot l_e \cdot D_e(t)$ | 对 cohesive 单元的损伤值加权求和 | ⚠️ 需新增：`D_e` 已存于 `damage_states`，单元长度 `l_e` 已有 |
| 边界热损失 $Q_{loss}$ | $\int_0^t \int_\Gamma h[T(\tau) - T_{amb}] \, d\Gamma \, d\tau$ | 从边界节点温度和换热系数积分 | ⚠️ 需新增：需提取边界节点温度和面积 |
| 化学能变化 $\Delta E_{chem}$ | $\int_V [U_{OCV}(SOC(t)) - U_{OCV}(SOC_0)] \cdot F \cdot c_{s,max} \, dV$ | 从各单元 SOC 变化和 OCV 曲线计算 | ⚠️ 需新增：`element soc_n/soc_p` 和 OCV 函数已有 |

### 7.3 简化方案

若完整化学能计算过于复杂，可先退化为**热-弹-断裂**子系统的能量守恒检查：

$$R_{sub}(t) = Q_{generated}(t) - \Delta E_{th}(t) - Q_{loss}(t) - \Delta E_{el}(t) - E_{frac}(t)$$

其中 $Q_{generated} = \int_0^t \sum_e q_e(\tau) A_e \, d\tau$ 直接从 `total heat source [W]` 积分得到（已可用）。

该子系统的优势：不依赖化学能定义，仅需确认热源→热能+热损失+弹性能+断裂能的平衡。

### 7.4 误差来源

- 时间积分误差（梯形法的截断误差）
- 空间离散误差（网格精度对能量积分精度的影响）
- 非线性迭代残差（CZM Newton-Raphson 未完全收敛时，能量不平衡）
- 归一化与量纲转换的累积舍入

### 7.5 验收标准

- $\epsilon_R(t) < 1\%$ 在整个仿真时长内
- 残余 $R(t)$ 不应有系统性增长趋势（否则表示能量泄漏）
- 在最细网格上应达到最小残差，粗网格上残差应更大

---

## 6. 计划评审记录（2026-04-22）

### 6.1 严重问题（必须修复）

#### P1：热网格由 CZM 特征长度决定——物理依据不足

当前计划将热网格的 `nθ` 区间由 cohesive 特征长度 $l_c$ 推导。这是**物理上不成立的**：

- **热分辨率需求**由热梯度决定（傅里叶数 $\text{Fo} = \alpha \Delta t / \Delta x^2$、Biot 数 $\text{Bi} = hL/k$）
- **CZM 分辨率需求**由过程区长度决定
- 两者是独立物理量，不存在因果关联

**建议**：热网格和 CZM 网格应有各自独立的分辨率判据。热网格可基于：
1. 经验法则（先粗后细，观察温度场收敛）
2. 热穿透深度 $\delta_T \sim \sqrt{\alpha t}$
3. 或直接给定 4 个 nθ 水平（如 60, 80, 120, 200）

#### P2：特征长度公式中的 $E$ 定义不清

$l_c = G_c \cdot E / \sigma_{\max}^2$ 是 Irwin 型过程区长度估计。但：
- 公式中的 $E$ 应为**体材料杨氏模量**，不是 cohesive 惩罚刚度 $K_n$
- 代码中 $K_n = 2.4 \times 10^{15}$ Pa/m 是刚度（不是模量），单位不同
- 需要从参数中找到或确认体材料的 $E$ 值

若代入法向参数做量级估算（假设需要找到体材料 E）：
- $G_c = 25.3$ J/m²，$\sigma_{\max} = 82$ MPa
- 若 $E \sim 100$ GPa（典型电极材料量级），则 $l_c \approx 100\times10^9 \times 25.3 / (82\times10^6)^2 \approx 0.38$ mm
- 若 $E \sim 1$ GPa，则 $l_c \approx 3.8$ μm

量级差三个数量级，**必须明确 $E$ 的来源和数值**。

#### P3：CZM 指标与电池仿真场景不匹配

计划列出的指标 "载荷-位移曲线" 和 "裂纹拓展速率" 与电池仿真场景有偏差：

- 电池 CZM 没有外部机械载荷，而是由**扩散应力 + 热应力**驱动
- 传统的 load-displacement curve 在此处无意义
- 应改为 **separation-time 曲线** 或 **D_max(t) 演化曲线**

**建议修正后的 CZM 指标**：
| 原指标 | 建议替换为 | 数据来源 |
|--------|-----------|----------|
| 应力峰值 | 保留 vonMises 峰值 | `diffusion stress vonMises` / `thermal stress vonMises` |
| 载荷-位移曲线 | separation-time 曲线 | `czm δ_max_n [m]` |
| 起裂时刻 | D_max 首达阈值时刻 | `czm D_max` |
| 裂纹拓展速率 | n_fractured(t) 斜率 | `czm n_fractured` |

### 6.2 重要问题（建议修复）

#### P4：电化学网格未包含颗粒分辨率

当前候选值只变化电解液网格 `(Nn, Ns, Np)`，未变动颗粒网格 `(Nrn, Nrp)`：
- 颗粒扩散方程的精度同样受网格影响
- 建议明确说明 Nrn/Nrp 是否固定（如固定为默认 10），还是纳入变化

#### P5："热模型基准" 已澄清，但热源冻结方式仍需确认

计划第 2.2 节现已明确为 "固定热源的热传导基准"，但实现时仍需确认：
- 热源如何施加？是均匀假设、还是从耦合算例中提取后固定？
- 若热源固定，不同 nθ 下热传导方程的网格敏感性才有意义
- 若热源仍然耦合 SPMe，则不是 "纯热" 而是 "电-热耦合"

#### P6：缺少时间步控制策略

不同网格方案对比时，时间离散误差会干扰空间离散误差的判断：
- 应固定 `dt` 或固定自适应策略参数
- 建议在执行前统一设定 `opt.dt` 为同一区间

#### P7：部分 CZM 输出变量未经代码确认

计划 §2.2 列出的以下变量**未在代码中确认存在**：
- `czm traction normal` — 待确认
- `czm traction tangent` — 待确认
- `czm separation tangent` — 待确认

已确认存在的 CZM 输出：
- `czm D_max`、`czm D_mean`、`czm n_fractured` ✅
- `czm δ_max_n [m]`、`czm δ_mean_n [m]` ✅

### 6.3 改进建议（可选）

#### S1：gsorder 应声明为固定参数

`gsorder` 影响积分精度，与 `nθ` 存在交互。应在计划中明确声明 gsorder 固定（如 2）。

#### S2：增加计算成本跟踪

网格敏感性分析应同时记录每组的运行时间和内存占用，用于绘制精度-成本权衡曲线。

#### S3：5% 阈值可按指标分化

| 指标 | 建议阈值 | 理由 |
|------|---------|------|
| 电压 | 1% | 电池模型精度要求高 |
| 温度峰值 | 2-3% | 工程可接受 |
| CZM 损伤 | 5% | 高度非线性，收敛慢 |

#### S4：热网格与 CZM 网格不一定要共用 nθ

热模型和 CZM 模型的网格分辨率需求不同：
- 热模型可能只需要较粗网格就能收敛
- CZM 需要 finer mesh 来捕捉过程区
- 强制共用会要么浪费热计算资源，要么 CZM 精度不足

### 6.4 评审总结

| 类别 | 项目数 | 结论 |
|------|--------|------|
| 严重问题 | 3 | P1 物理依据不足、P2 参数不明、P3 指标错配 |
| 重要问题 | 4 | P4-P7 需补充说明 |
| 改进建议 | 4 | S1-S4 可提升计划质量 |
| 总体评价 | — | 计划框架清晰、交付物定义完整，但核心物理动机需修正后再执行 |

## 8. 统计指标体系升级（2026-04-29）

### 8.1 问题

原有三个 Track 的所有对比指标都基于单点值（snapshot），对局部波动敏感，无法反映整条曲线或空间场的整体收敛质量。

### 8.2 解决方案

将所有指标统一替换为基于 RMSPE 的统计量：

- **电化学 Track**: V(t) RMSPE、T(t) RMSPE、dT/dt(t) RMSPE
- **热学 Track**: T_max(t) RMSPE、T_range(t) RMSPE、空间场 RMSPE（时间平均）
- **CZM Track**: D_max(t) RMSPE、n_frac(t) RMSPE、δ_max_n(t) RMSPE、牵引-分离面积偏差
- **能量守恒**: 保留 ε_R(t) 瞬时值，新增归一化 RMS 残余 ε_R,rms

### 8.3 关键设计决策

- 误差公式：RMSPE（相对均方根百分比误差），带零点保护（threshold = 1e-3 * max|y_ref|）
- 验收阈值：统一 5%
- 时间对齐：手写线性插值 `align_to_ref`，不依赖外部包
- 牵引-分离面积偏差：选取最终时刻 D 值最大的单元进行对比
- 能量残余：归一化 RMS（非 RMSPE，避免除零）

### 8.4 旧指标放弃理由

- 角变化收敛：Bi_t ≈ 0.004 导致角变化极小，RMSPE 零点保护大量跳过
- 应力峰值：空间分布不均匀，单点意义有限
- 损伤起始时间：事件时间，RMSPE 不适用
- 载荷-位移曲线：需纯机械位移 BC，与电池实际驱动不符

### 8.5 参考文件

- Spec: `docs/superpowers/specs/2026-04-29-grid-sensitivity-statistical-metrics-design.md`
- Plan: `docs/superpowers/plans/2026-04-29-statistical-metrics-implementation-plan.md`
