> 本文档由 `电化学-热-力耦合理论推导.md` 第 2828–4294 行切分而来，版本：B 方案 + (s,n) 坐标 + 8 层 CLT + Phi 伪周期映射。

---
## CZM 子网格：圈内材料界面内聚力模型（C-skip-thermal）

本章提出 **C-skip-thermal CZM 子网格架构**：在每个 8 层卷绕重复单元的细化力学网格上独立求解力学-CZM 耦合，跳过子网格热场求解，通过与粗热网格的双向插值实现热-力耦合。cohesive 本构只包含 PE–PCC 与 NE–NCC 两种材料类型，但双面涂布对应四个真实箔–涂层面；离散后共有 $N_{\mathrm{elem}}^{\mathrm{coh}}=4N_{\mathrm{seg}}$ 个 cohesive 单元。SP–PE 与 SP–NE 接触面不计入该数。

**章节结构**：
1. **§3.0a CZM 子网格几何**：细网格离散化（8 Q4/匝）、材料类型标记、`thermal_elem_map`；
2. **§3.0b 界面识别理论**：径向相邻 Q4 对的材料组合判定、cohesive 插入规则、预期 cohesive 计数；
3. **§3.0c 节点复制策略**：共边节点复制、Q4 连接表重写、Dict 记忆化、自检；
4. **§3.0d 非内聚界面（v3.0）**：SP–涂层单边接触 + Coulomb 摩擦（**禁造 CZM 参数**）；
5. **§3.0e 塌陷状态变量与耦合链（v3.0）**：$\kappa$、$\Delta_{\mathrm{core}}$、$D$（仅箔–涂层）；
6. **§3.1 界面运动学**（修订）：CZM 子网格内的 PE-PCC 与 NE-NCC 材料界面运动学；
7. **§3.2–§3.5**（保留）：TSL、损伤演化、B-K 准则、Q4 残余应力映射；
8. **§3.6 损伤下游反馈**（含 §3.6.0a/b 粗↔细插值）。

**C-skip-thermal 含义**：CZM 子网格上的温度场不从热方程求解，而是从粗热网格插值获得（§3.6.0a）。反向地，CZM 子网格上的损伤 $D$ 通过 max-归约映射回粗热网格（§3.6.0b），实现损伤→热参数的闭环。

**输入/输出接口**：

- **输入**：来自 §2.7 的界面 traction $T_n,T_s$（式 2.61–2.64）；来自粗热网格的温度 $T$ 和 $\Delta\mathrm{soc}$（经 §3.6.0a 插值）；
- **输出**：损伤标量 $D\in[0,1]$（per-interface），作为 Curie 中介（§2.8.2 式 2.70）驱动电化学/热反馈（§3.6 详述 $R_{\mathrm{contact}},A_{\mathrm{eff}},k_n^{\mathrm{eff}}$）；$D$ 经 max-归约（§3.6.0b）回传至粗热网格；
- **辅助输出**：等效应变能释放率 $G$ 与混合模式临界断裂韧度 $G_c^{\mathrm{mix}}$，供 §3.4 损伤驱动力分析调用。

**记号约定**：本章公式编号自 (3.0a) 起用于新增子节，(3.1) 起用于保留的 §3.1+。如无特殊声明，所有量均为 (s, n) 平面截面内的量。

---

---

### §3.0a CZM 子网格几何（C-skip-thermal 细网格离散化）

**动机**：粗热网格采用 1 个 Q4 单元/匝（径向分辨率 = $t_{\mathrm{repeat}}=373.6\;\mu\mathrm{m}$），无法分辨 8 层 CLT 内部的材料界面。CZM 需要在 PE-PCC 与 NE-NCC 界面处求解位移跳跃和损伤演化，因此必须引入更细的力学子网格。

**子网格离散化**：每个卷绕圈（$n\in[0, t_{\mathrm{repeat}}]$）沿径向分为 $N_{\mathrm{lay}}=8$ 个 Q4 单元，分别对应 8 层 CLT 材料层：

$$\boxed{\;\text{CZM 子网格：每匝 } N_{\mathrm{lay}}=8 \text{ 径向 Q4 单元，对应 PE-PCC-PE-SP-NE-NCC-NE-SP 层序}\;}\tag{3.0a}$$

沿 $s$ 方向（卷绕切向），子网格直接继承粗热网格的角节点。令 $N_{\mathrm{seg}}=\mathtt{length(theta)}-1$ 为整条螺旋的周向分段总数，则子网格 Q4 总数为 $8N_{\mathrm{seg}}$。输入参数 $n_\theta$ 表示典型的每整匝分段目标，不能替代包含末端部分匝的精确总数 $N_{\mathrm{seg}}$。

**尺度分离假设**（C-skip-thermal 核心）：

| 网格 | 径向分辨率 | $s$ 方向分辨率 | 求解场 | 每周向分段单元数 |
|------|-----------|---------------|--------|----------|
| 粗热网格 | 1 个重复单元厚度 | $N_{\mathrm{seg}}$ 段 | $T, \Delta\mathrm{soc}_p, \Delta\mathrm{soc}_n$ | 1 |
| CZM 子网格 | 8 个材料层 | $N_{\mathrm{seg}}$ 段 | $u_s, u_n, D$（CZM 界面） | 8 |

**物理论证**：温度 $T$ 在径向的变化由热扩散方程控制，特征长度 $\ell_T = \sqrt{k_r/(\rho C_p \omega)}$。对于典型充放电频率 $\omega\sim 10^{-3}$–$10^{-2}$ Hz，$\ell_T\sim 1$–$10$ mm $\gg t_{\mathrm{repeat}}=373.6\;\mu\mathrm{m}$，故温度在一个卷绕周期内接近线性，1 单元/匝足以分辨。力学场（应力、应变）在材料界面处不连续（模量跳变），需要 8 单元/匝才能正确捕捉界面 traction。C-skip-thermal 利用这一尺度分离：**CZM 子网格跳过热场求解**（因此得名"skip-thermal"），温度从粗热网格插值获取（§3.6.0a）。

**材料类型标记**：CZM 子网格中每个 Q4 单元 $e$ 携带材料类型标记：

$$\boxed{\;\mathtt{material\_type}[e] \in \{\mathrm{PE}, \mathrm{PCC}, \mathrm{SP}, \mathrm{NE}, \mathrm{NCC}\}\;}\tag{3.0b}$$

标记由 8 层 CLT 层序确定（§1.2.2 表）：
- 层 1, 3 → PE
- 层 2 → PCC
- 层 4, 8 → SP
- 层 5, 7 → NE
- 层 6 → NCC

**`thermal_elem_map` — O(1) 解析反查**：CZM 子网格单元到粗热网格单元的映射为闭式解析函数（无需搜索），实现 O(1) 查找：

$$\boxed{\;\mathtt{thermal\_elem\_map}[e] = ((e-1)\bmod N_{\mathrm{seg}})+1\;}\tag{3.0c}$$

其中 $e=1,\ldots,8N_{\mathrm{seg}}$ 为按“材料层优先、周向分段次之”排序的 CZM 子网格 Q4 编号；$\mathtt{thermal\_elem\_map}[e]$ 返回其所属粗热单元编号（1-indexed）。

**子网格单元中心坐标 → (匝号, 周向段号) 解析反查**：

$$\boxed{\;\begin{aligned}
\mathtt{segment\_id}[e] &= ((e-1)\bmod N_{\mathrm{seg}})+1,\\
\mathtt{turn\_id}[e] &= \left\lfloor\frac{\theta_{\mathrm{center}}[\mathtt{segment\_id}[e]]-\theta_1}{2\pi}\right\rfloor+1
\end{aligned}\;}\tag{3.0d}$$

用于在周向方向定位 cohesive 界面的 $s$ 坐标范围。

---

### §3.0b 界面识别理论（径向相邻 Q4 对的材料组合判定）

**问题陈述**：给定 CZM 子网格中两个径向相邻的 Q4 单元 $(e_{\mathrm{inner}}, e_{\mathrm{outer}})$（$e_{\mathrm{inner}}$ 靠近 $n=0$，$e_{\mathrm{outer}}$ 靠近 $n=t_{\mathrm{repeat}}$），判定是否在该界面插入 cohesive 单元。

**判定流程**：

1. 获取相邻单元的材料类型：
   $$\mathtt{mat\_inner} = \mathtt{material\_type}[e_{\mathrm{inner}}],\quad \mathtt{mat\_outer} = \mathtt{material\_type}[e_{\mathrm{outer}}]$$

2. 查表判定：
   $$
   \boxed{\;\text{cohesive 判定} =
   \begin{cases}
   \text{插入}: & \{\mathtt{mat\_inner},\mathtt{mat\_outer}\}=\{\mathrm{PE},\mathrm{PCC}\}\;\text{或}\;\{\mathrm{NE},\mathrm{NCC}\} \\
   \text{不插入}: & \text{其他所有组合}
   \end{cases}\;}\tag{3.0e}
   $$

**材料组合判定表**（8 层 CLT 中所有径向相邻对的完整枚举）：

| $e_{\mathrm{inner}}$ 层号 | $e_{\mathrm{outer}}$ 层号 | 材料组合 | cohesive 插入？ | 界面类型标记 |
|---|---|---|---|---|
| 1 (PE) | 2 (PCC) | PE-PCC | **是** | `:PE_PCC` |
| 2 (PCC) | 3 (PE) | PCC-PE | **是** | `:PE_PCC` |
| 3 (PE) | 4 (SP) | PE-SP | 否（无实验数据） | — |
| 4 (SP) | 5 (NE) | SP-NE | 否（无实验数据） | — |
| 5 (NE) | 6 (NCC) | NE-NCC | **是** | `:NE_NCC` |
| 6 (NCC) | 7 (NE) | NCC-NE | **是** | `:NE_NCC` |
| 7 (NE) | 8 (SP) | NE-SP | 否（无实验数据） | — |

**排除 SP-PE 与 SP-NE 的理由（升格为结构性决策，v3.0）**：
1. **无实验内聚参数**：不得编造 $G_c$、$T_{\max}$；现有脱粘表征仅在 PE-PCC / NE-NCC 可量化；
2. **物理机制**：隔膜柔顺，界面应力远低于金属–涂层硬界面；塌陷前兆是几何皱褶而非内聚开裂；
3. **工况 C 替代本构**：不插入 cohesive，而采用 §3.0d 接触–摩擦；若未来获得可靠剥离数据，另开版本扩展判定表，不在本文主路径内。

**预期 cohesive 计数**：

$$\boxed{\;N_{\mathrm{type}}^{\mathrm{coh}}=2,\qquad N_{\mathrm{face,repeat}}^{\mathrm{coh}}=4,\qquad N_{\mathrm{elem}}^{\mathrm{coh}}=4N_{\mathrm{seg}}\;}\tag{3.0f}$$

该计数直接按整条螺旋的实际角节点计算，末端部分匝已经包含在 $N_{\mathrm{seg}}$ 中，不需要按整匝向下取整或设置 $\pm n_\theta$ 计数容差。`interface_type` 的两个取值只决定参数族，不会把四个真实面合并成两个离散面。

---

### §3.0c 节点复制策略（共边 2 节点 → 独立 DOF → δ ≠ 0）

**目标**：在 CZM 子网格中，cohesive 界面两侧的节点在参考构型下共享同一坐标，但需要独立的位移自由度以允许非零位移跳跃 $\delta_n,\delta_s \neq 0$。这通过对 cohesive 界面上的节点进行**复制**实现。

**操作**：对每个被判定为 cohesive 的径向相邻 Q4 对 $(e_{\mathrm{inner}}, e_{\mathrm{outer}})$，其共享边上存在 2 个节点（Q4 单元边含 2 节点）。对每个共享节点 $p$（坐标 $\mathbf{x}_p$）：

1. **复制节点**：创建副本节点 $p'$，坐标 $\mathbf{x}_{p'} = \mathbf{x}_p$（参考构型下几何重合）；
2. **DOF 独立**：$p$ 的位移 $\mathbf{u}_p = (u_s^{(p)}, u_n^{(p)})$ 与 $p'$ 的位移 $\mathbf{u}_{p'} = (u_s^{(p')}, u_n^{(p')})$ 独立——$p$ 属于外侧单元，$p'$ 属于内侧单元（或反之）；
3. **位移跳跃**：$\delta_n = u_n^{(p)} - u_n^{(p')}$，$\delta_s = u_s^{(p)} - u_s^{(p')}$（符号约定与 §3.1 一致）。

**外层 Q4 连接表重写**：对含 cohesive 界面的 Q4 单元，其节点连接表需从原共享节点更新为副本节点。设原 Q4 单元 $e_{\mathrm{outer}}$ 的连接为 $\{n_1, n_2, n_3, n_4\}$，若其某边与 cohesive 界面重合，该边上的 2 个节点替换为对应副本节点：

$$\boxed{\;\mathtt{conn}[e_{\mathrm{outer}}] = \{n_1', n_2', n_3, n_4\}\quad\text{（假设边 }(n_1,n_2)\text{ 为 cohesive 共享边）}\;}\tag{3.0g}$$

**关键约束**：$\delta \equiv 0$ 是必须避免的错误模式——若连接表未重写（内外侧单元仍共享同一节点），则位移强制连续，CZM 退化。节点复制 + 连接表重写确保 $\delta$ 可以非零。

**Dict 记忆化（角点共享处理）**：一个节点可能被多个 cohesive 界面共享（角点 case）。使用字典（HashMap）记录已复制的节点：

$$\boxed{\;\mathtt{copy\_dict}[\mathtt{node\_id}] = \mathtt{new\_node\_id}\;}\tag{3.0h}$$

在复制前先查询 `copy_dict`——若 `node_id` 已存在副本，直接使用已有副本 ID（避免重复复制导致虚假的多副本）。

**三条自检**（实现后必经验证）：

1. **坐标一致性**：$\|\mathbf{x}_p - \mathbf{x}_{p'}\| < \epsilon_{\mathrm{coord}}$（参考构型下副本与原件坐标重合，$\epsilon_{\mathrm{coord}} \sim 10^{-12}$ m）；
2. **4 节点不重复**：$\mathtt{conn}[e]$ 中 4 个节点 ID 互不相同（无自重合退化）；
3. **外层用副本**：cohesive 界面的外侧单元连接表使用副本节点 ID（$p'$），内侧单元保留原节点 ID（$p$）。

**"惩罚等价"论证（$D=0$ 时）**：当损伤 $D=0$（界面完好），双线性 TSL（§3.2）处于弹性段，traction 与分离满足 $T_n = K_n\delta_n$、$T_s = K_s\delta_s$。若惩罚刚度 $K_n, K_s$ 足够大（$\gg$ 相邻材料的弹性模量），即使 $\delta$ 可以非零，其量级被压制到 $\sim \sigma_{\mathrm{interface}}/K_n \ll$ 弹性应变。在 $K_n\to\infty$ 极限下，$\delta\to 0$，位移恢复连续。这一性质保证：**无损伤时子网格近似为连续介质；有损伤时 $\delta$ 逐渐增大，CZM 平滑激活**。这是节点复制策略在物理上的自洽性保证。

---

### §3.0d 非内聚界面：SP–涂层接触–摩擦（v3.0）

> **[v3.0 / 界面分流]**：SP–PE、SP–NE **不是** CZM 界面。无实验内聚参数 ⇒ **禁止**写入 $G_c$、$T_{\max}$、$K_n$ 等 CZM 量。工况 C 采用单边接触 + 摩擦，以支持隔膜先皱（Shi 2026）。

**法向（单边接触，零抗拉强度）**：记间隙 $g_n\geq 0$（$g_n=0$ 接触）。受压时法向牵引 $t_n\leq 0$（罚函数或增广 Lagrange）；**$t_n$ 不得为正**——隔膜不“粘住”涂层，可局部失接触/起皱离开电极。

**切向（Coulomb 摩擦）**：
$$|t_s|\leq\mu|t_n|,\qquad
t_s=-\mu|t_n|\,\mathrm{sgn}(v_{\mathrm{slip}})\quad\text{（滑移时）}\tag{3.0i}$$

基准 $\mu=0.10$（Shi）；敏感性 $\mu\in[0.05,0.4]$。预屈曲 hoop 应力场形态对 $\mu$ 弱敏感（幅度变化 ≲20%）。**不设**粘聚剪切强度，**不设**混合模式能量释放率。

**与箔–涂层 CZM 的对比（必须保持）**：

| | PE–PCC / NE–NCC | SP–PE / SP–NE |
|---|---|---|
| 本构 | 牵引–分离律 + 损伤 $D$ | 接触–摩擦，无 $D$ |
| 参数 | 实验 $G_c,T_{\max}$ | 仅 $\mu$（及接触罚参数） |
| 塌陷角色 | 失稳后脱粘、电化学反馈 | 可滑移 → 周向长度失配 → **先皱前兆** |

**先皱机制**：放电时电极栈收缩、$\beta_{\mathrm{SP}}=0$ 且可相对滑移 → 隔膜周向相对“过长” → 压缩膜力 + 极低弯曲刚度 → 出平面皱褶。该路径依赖可滑移；若误设完美粘结会把涂层 eigenstrain 强加给隔膜，破坏先皱机理（§7 极限 C5 负例）。

**弱形式**：接触虚功 $\delta W_{\mathrm{contact}}$ 与 CZM 虚功 $\delta W_{\mathrm{CZM}}$ **分列**（§6.2.3b），Jacobian 接触切线不得并入 cohesive 公式。

---

### §3.0e 塌陷状态变量与耦合链（v3.0）

**状态变量**（工况 C / L4）：

| 符号 | 含义 | 演化位置 |
|---|---|---|
| $\kappa$ | PCC/NCC 等效塑性应变 | §2.1.6；跨半循环保留 |
| $\Delta_{\mathrm{core}}$ | 内圈不规则指标（内径圆度损失或法向皱褶幅值范数） | 几何后处理 / 模态幅值 |
| $D$ | 箔–涂层 CZM 损伤 | §3.3；**仅** PE–PCC / NE–NCC |

**循环映射**：每半循环更新 $\varepsilon^p$、$\kappa$ 与当前几何（或等效缺陷幅值）；**禁止**每圈重置为完美圆参考构型。单循环全塌陷**不是**成功标准；以多圈 $\Delta_{\mathrm{core}}$ 单调恶化为预言。

**耦合链**（界面分流）：
1. SP 皱褶 / 内圈失稳 → 局部改变涂层曲率与膜力；
2. → **涂层–箔**界面牵引升高 → 仅在 PE–PCC / NE–NCC 上演化 $D$；
3. $D\uparrow$ → 箔–涂层相容削弱 → 应力重分配 → 加速内陷；
4. SP–涂层始终无 $D$，仅接触开闭与滑移状态切换；
5. 可选：$D\to A_{\mathrm{eff}},R_{\mathrm{contact}}$（复用 §3.6）。

---

### §3.1 界面运动学（CZM 子网格修订）

本节建立 CZM 子网格内 PE-PCC 与 NE-NCC 材料界面上的界面运动学，包括界面坐标定义、相对位移分解、等效分离与混合模式比。本节内容与 §1.5（界面几何与相对位移）严格衔接，不重复定义已建立的几何量，仅补充 CZM 所需的运动学推导。

#### §3.1.1 界面坐标与相对位移回顾

**回顾**（来自 §1.5）：对 CZM 子网格内 PE-PCC 或 NE-NCC 材料界面上一点，取界面法向单位向量 $\mathbf{n}$、切向单位向量 $\mathbf{t}$（满足 $\mathbf{n}\perp\mathbf{t}$，构成右手系）。相邻两 Q4 单元在界面两侧分别记为 $+$ 侧与 $-$ 侧，位移场分别为 $\mathbf{u}^+$ 与 $\mathbf{u}^-$。

界面相对位移分解为法向分离与切向滑移：
$$
\delta_n=(\mathbf{u}^+-\mathbf{u}^-)\cdot\mathbf{n},\qquad \delta_s=(\mathbf{u}^+-\mathbf{u}^-)\cdot\mathbf{t} \tag{3.1}
$$

在 (s, n) 坐标下，材料界面法向沿 $n$ 方向（$\mathbf{n}=\mathbf{e}_n$），切向沿 $s$ 方向（$\mathbf{t}=\mathbf{e}_s$），式 (3.1) 退化为
$$
\delta_n=u_n^+-u_n^-,\qquad \delta_s=u_s^+-u_s^- \tag{3.2}
$$

其中 $u_n,u_s$ 分别为法向与切向位移分量（§1.4.2 定义）。

**符号约定**：
- $\delta_n>0$：界面张开（法向分离为正，对应 mode I 张开型）；
- $\delta_n<0$：界面闭合（法向压缩，CZM 通常不计入压缩损伤，见 §3.2.4）；
- $\delta_s$ 的符号约定为：$\delta_s>0$ 表示 $+$ 侧相对 $-$ 侧沿 $+\mathbf{t}$ 方向滑移（mode II 剪切型）。

**注 3.1.1（与 §2.7 traction 的对应）**：$\delta_n$ 驱动法向 traction $T_n$（式 2.62，$T_n=\sigma_{nn}$），$\delta_s$ 驱动切向 traction $T_s$（式 2.64，$T_s=\sigma_{sn}$）。$\delta_n$-$T_n$ 构成 mode I CZM，$\delta_s$-$T_s$ 构成 mode II CZM，混合模式见 §3.1.3。

**注 3.1.2（(s, n) 平面简化）**：(s, n) 平面应力假设下（§2.6），$\sigma_{zz}=\sigma_{sz}=\sigma_{nz}=0$，故材料界面上相对位移仅含 $(\delta_n,\delta_s)$ 二维分量（法向脱粘 + 切向滑移），轴向滑移不独立出现（由假设 1.2 消去）。这是后续二维 CZM 的几何基础。

#### §3.1.2 等效分离与加权形式

**定义 3.1.1（等效分离，等权形式）**：定义等效分离为
$$
\delta_{\mathrm{eff}}=\sqrt{\delta_n^2+\delta_s^2} \tag{3.3}
$$

$\delta_{\mathrm{eff}}$ 为非负标量，描述界面相对位移的合模值。$\delta_{\mathrm{eff}}=0$ 对应界面完全粘结，$\delta_{\mathrm{eff}}>0$ 对应界面发生分离。

**物理诠释**：式 (3.3) 将法向与切向分离以 Euclid 范数合成，对应"等权"假设——mode I 与 mode II 对损伤的贡献同等加权。这一假设在 $G_{Ic}=G_{IIc}$（两种模式断裂韧度相等）时严格成立。

**定义 3.1.2（等效分离，加权形式）**：当 mode I 与 mode II 断裂韧度不同（$G_{Ic}\neq G_{IIc}$）时，引入法向权重 $\beta_n$ 修正等效分离：
$$
\delta_{\mathrm{eff}}=\sqrt{(\beta_n\delta_n)^2+\delta_s^2} \tag{3.4}
$$

其中 $\beta_n$ 为正实数。

**$\beta_n$ 的物理意义**：$\beta_n$ 调节法向分离对等效分离（进而对损伤）的贡献权重。
- $\beta_n=1$：退化为等权形式（式 3.3）；
- $\beta_n>1$：放大法向分离贡献（mode I 主导场景，如涂层剥离）；
- $\beta_n<1$：抑制法向分离贡献（mode II 主导场景，如层间剪切）。

**$\beta_n$ 的经验选择**：典型地，$\beta_n$ 可由临界断裂韧度比确定：
$$
\beta_n\sim\sqrt{\frac{G_{IIc}}{G_{Ic}}} \tag{3.5}
$$

式 (3.5) 的物理动机为：使 mode I 与 mode II 在达到各自临界断裂韧度时贡献等效的 $\delta_{\mathrm{eff}}$ 增量。本文默认取 $\beta_n=1$（等权形式），除非特殊说明。

**推论 3.1.1（$\delta_{\mathrm{eff}}$ 的非负性与单调性）**：由式 (3.3) 或 (3.4)，$\delta_{\mathrm{eff}}\geq 0$ 恒成立。$\delta_{\mathrm{eff}}$ 关于 $|\delta_n|,|\delta_s|$ 单调递增，但 $\delta_{\mathrm{eff}}$ 关于 $\delta_n$ 的单调性取决于符号（$\delta_n<0$ 时增大 $|\delta_n|$ 亦增大 $\delta_{\mathrm{eff}}$，与物理上压缩不应增加损伤矛盾——此问题在 §3.2.1 通过仅对 $\delta_n>0$ 计入损伤解决）。

#### §3.1.3 混合模式定义与混合模式比

**定义 3.1.3（三种加载模式）**：根据 $(\delta_n,\delta_s)$ 的组合，界面加载模式分为：

1. **Pure mode I（纯张开型）**：
   $$\delta_s=0,\quad \delta_n\neq 0\quad\Longrightarrow\quad \delta_{\mathrm{eff}}=|\delta_n| \tag{3.6}$$

2. **Pure mode II（纯剪切型）**：
   $$\delta_n=0,\quad \delta_s\neq 0\quad\Longrightarrow\quad \delta_{\mathrm{eff}}=|\delta_s| \tag{3.7}$$

3. **Mixed mode（混合模式）**：
   $$\delta_n\neq 0,\quad \delta_s\neq 0 \tag{3.8}$$

**定义 3.1.4（混合模式比）**：定义混合模式比为
$$
\beta_{\mathrm{mix}}=\frac{\delta_s^2}{\delta_n^2+\delta_s^2}\in[0,1] \tag{3.9}
$$

**$\beta_{\mathrm{mix}}$ 的物理意义**：
- $\beta_{\mathrm{mix}}=0$：纯 mode I（$\delta_s=0$）；
- $\beta_{\mathrm{mix}}=1$：纯 mode II（$\delta_n=0$）；
- $\beta_{\mathrm{mix}}\in(0,1)$：混合模式，$\beta_{\mathrm{mix}}$ 越大表示 mode II 贡献越显著。

**注 3.1.3（$\beta_{\mathrm{mix}}$ 与能量比的对应）**：在线弹性界面（$T_n=K_n\delta_n,T_s=K_s\delta_s$，且 $K_n=K_s$）假设下，mode I 与 mode II 弹性能分别为 $G_I=\frac{1}{2}K_n\delta_n^2$、$G_{II}=\frac{1}{2}K_s\delta_s^2$，故
$$
\beta_{\mathrm{mix}}=\frac{G_{II}}{G_I+G_{II}} \tag{3.10}
$$
即 $\beta_{\mathrm{mix}}$ 在等刚度假设下等于 mode II 能量分数。式 (3.10) 用于 §3.3.4 的 B-K 混合模式准则。

**注 3.1.4（$\beta_{\mathrm{mix}}$ 的演化与循环加载）**：循环加载过程中，$(\delta_n,\delta_s)$ 随时间变化，$\beta_{\mathrm{mix}}$ 亦随之演化。$\beta_{\mathrm{mix}}$ 的变化反映加载路径的"模式旋转"——从 mode I 主导转向 mode II 主导或反之。这一演化影响 $G_c^{\mathrm{mix}}$（§3.3.4 式 3.39），进而影响损伤速率。

**推论 3.1.2（$\beta_{\mathrm{mix}}$ 与 $\delta_{\mathrm{eff}}$ 的独立性）**：$\delta_{\mathrm{eff}}$ 描述分离的"大小"（合模值），$\beta_{\mathrm{mix}}$ 描述分离的"方向"（模式分数）。二者共同完整描述界面运动学状态，$(\delta_{\mathrm{eff}},\beta_{\mathrm{mix}})$ 构成等效运动学坐标。

**注 3.1.5（符号约定对 $\beta_{\mathrm{mix}}$ 的影响）**：式 (3.9) 中 $\delta_n,\delta_s$ 均以平方形式出现，$\beta_{\mathrm{mix}}$ 与符号无关。这一特性使 $\beta_{\mathrm{mix}}$ 在卸载-再加载过程中保持连续（不会因 $\delta_n$ 过零而跳跃），有利于数值稳定性。

#### §3.1.4 界面运动学接口

本小节给出 §3.1 的输出接口，供 §3.2 TSL 与 §3.3 损伤演化调用。

**输出 1（等效分离）**：
- $\delta_{\mathrm{eff}}$（式 3.3 或 3.4）：进入 §3.2 TSL 的自变量；
- 用途：双线性 TSL 分段判据（§3.2.2）、P-P TSL 形函数自变量（§3.3.5）。

**输出 2（混合模式比）**：
- $\beta_{\mathrm{mix}}$（式 3.9）：进入 §3.3.4 B-K 混合模式准则；
- 用途：计算混合模式临界断裂韧度 $G_c^{\mathrm{mix}}$（式 3.39）。

**输出 3（分量分离）**：
- $(\delta_n,\delta_s)$（式 3.1–3.2）：进入 §3.2.4 法向-切向解耦 TSL；
- 用途：分别计算 mode I 与 mode II 能量释放率 $G_I,G_{II}$（§3.3.3）。

---

### §3.2 牵引-分离律（TSL）

本节建立 $\Gamma_{\mathrm{coh}}$ 上的牵引-分离律（Traction-Separation Law, TSL），描述界面 traction $T$ 与等效分离 $\delta_{\mathrm{eff}}$ 之间的本构关系。TSL 是 CZM 的核心本构律，决定损伤起始、软化速率与最终脱粘行为。本文默认采用双线性 TSL（计算效率高、参数标定简单），并在 §3.3.5 简述 Park-Paulino unified TSL 作为可选替代。

#### §3.2.1 双线性 TSL 三段定义

**假设 3.2.1（双线性 TSL 形式）**：$\Gamma_{\mathrm{coh}}$ 上的等效 traction $T$ 与等效分离 $\delta_{\mathrm{eff}}$ 之间满足分段线性关系，分为弹性段、软化段、完全分离段三段。

**定义 3.2.1（等效 traction）**：定义等效 traction 为
$$
T=\sqrt{T_n^2+T_s^2} \tag{3.11}
$$
其中 $T_n,T_s$ 由 §2.7.3 式 (2.61)–(2.64) 给出。

**三段定义**：

**弹性段**（$0\leq\delta_{\mathrm{eff}}<\delta_0$）：
$$
T=K\,\delta_{\mathrm{eff}} \tag{3.12}
$$
其中 $K$ 为界面初始刚度（penalty stiffness，单位 N/mm$^3$），$\delta_0$ 为损伤起始分离。

**软化段**（$\delta_0\leq\delta_{\mathrm{eff}}<\delta_c$）：
$$
T=K\,(1-D)\,\delta_{\mathrm{eff}} \tag{3.13}
$$
其中 $D\in[0,1]$ 为损伤变量（§3.3 严格定义），单调递增。

**完全分离段**（$\delta_{\mathrm{eff}}\geq\delta_c$）：
$$
T=0 \tag{3.14}
$$
对应 $D=1$，界面完全脱粘，不再传递 traction。

**关键参数清单**：
- $\delta_0$：损伤起始分离（对应峰值 traction $T_{\max}=K\delta_0$，由式 3.12 取 $\delta_{\mathrm{eff}}\to\delta_0^-$）；
- $\delta_c$：临界分离（对应 $D=1$，TSL 终止点）；
- $T_{\max}=K\delta_0$：峰值 traction（ cohesive strength）；
- $G_c$：临界断裂韧度，定义为 TSL 包络面积：
  $$
  G_c=\int_0^{\delta_c}T(\delta_{\mathrm{eff}})\,d\delta_{\mathrm{eff}}=\frac{1}{2}T_{\max}\delta_c=\frac{1}{2}K\delta_0\delta_c \tag{3.15}
  $$

**物理诠释**：
- 弹性段（式 3.12）：界面行为类似线性弹簧，刚度 $K$ 表征界面的"弹性约束"能力。$\delta_{\mathrm{eff}}<\delta_0$ 时无损伤（$D=0$），卸载后位移完全恢复；
- 软化段（式 3.13）：损伤累积，$D$ 从 $0$ 单调增至 $1$，有效刚度 $K(1-D)$ 递减，traction 从 $T_{\max}$ 线性降至 $0$；
- 完全分离段（式 3.14）：界面彻底脱粘，traction 消失，$D=1$ 锁定（不可逆）。

**注 3.2.1（$G_c$ 的物理意义）**：$G_c$（式 3.15）为单位界面面积完全脱粘所需消耗的能量（单位 J/m$^2$），是 CZM 的核心断裂力学参数。$G_c$ 由实验标定（如 DCB 试验测 mode I、ENF 试验测 mode II），典型量级 $G_c\sim 10^{-2}-10^2$ J/m$^2$（视界面性质而定，仅给量级）。

**注 3.2.2（$\delta_0,\delta_c$ 的相对关系）**：由式 (3.15)，$\delta_c>\delta_0$（否则 $G_c\leq 0$，无物理意义）。比值 $\delta_c/\delta_0$ 决定 TSL 的"脆性"——比值越小（$\delta_c\to\delta_0$），TSL 越接近理想脆性断裂（$G_c\to\frac{1}{2}K\delta_0^2$）；比值越大，TSL 越延性。典型 $\delta_c/\delta_0\sim 5-20$。

**推论 3.2.1（$K,\delta_0,\delta_c$ 的独立性）**：三个参数 $K,\delta_0,\delta_c$ 中，仅有两个独立——第三个由式 (3.15) $G_c=\frac{1}{2}K\delta_0\delta_c$ 约束。实际标定时，通常给定 $(T_{\max},G_c)$ 或 $(\delta_0,\delta_c)$，第三个参数由约束解出。

#### §3.2.2 双线性 TSL 的分段解析形式

**定理 3.2.1（双线性 TSL 解析表达）**：在假设 3.2.1 下，等效 traction $T$ 关于 $\delta_{\mathrm{eff}}$ 的完整解析形式为
$$
T(\delta_{\mathrm{eff}})=\begin{cases}
K\,\delta_{\mathrm{eff}} & 0\leq\delta_{\mathrm{eff}}<\delta_0\\[4pt]
T_{\max}\,\dfrac{\delta_c-\delta_{\mathrm{eff}}}{\delta_c-\delta_0} & \delta_0\leq\delta_{\mathrm{eff}}<\delta_c\\[6pt]
0 & \delta_{\mathrm{eff}}\geq\delta_c
\end{cases} \tag{3.16}
$$

**推导**：
1. 弹性段（$0\leq\delta_{\mathrm{eff}}<\delta_0$）：由式 (3.12) 直接给出 $T=K\delta_{\mathrm{eff}}$。
2. 软化段（$\delta_0\leq\delta_{\mathrm{eff}}<\delta_c$）：由式 (3.13) $T=K(1-D)\delta_{\mathrm{eff}}$ 与 Camanho 损伤演化（式 3.31 第二段 $D=\delta_c(\delta_{\mathrm{eff}}-\delta_0)/[\delta_{\mathrm{eff}}(\delta_c-\delta_0)]$），代入得
   $$
   T=K\left(1-\frac{\delta_c(\delta_{\mathrm{eff}}-\delta_0)}{\delta_{\mathrm{eff}}(\delta_c-\delta_0)}\right)\delta_{\mathrm{eff}}=K\delta_0\cdot\frac{\delta_c-\delta_{\mathrm{eff}}}{\delta_c-\delta_0}=T_{\max}\,\frac{\delta_c-\delta_{\mathrm{eff}}}{\delta_c-\delta_0} \tag{3.17}
   $$
   其中 $T_{\max}=K\delta_0$（弹性段峰值 traction，由 (3.12) 取 $\delta_{\mathrm{eff}}\to\delta_0^-$）。式 (3.17) 的几何意义为：软化段从 $(\delta_0,T_{\max})$ 线性下降至 $(\delta_c,0)$。
3. 完全分离段（$\delta_{\mathrm{eff}}\geq\delta_c$）：$T=0$。

**连续性检验**：
- 在 $\delta_{\mathrm{eff}}=\delta_0$ 处：弹性段 $T=K\delta_0=T_{\max}$，软化段 $T=T_{\max}\cdot\frac{\delta_c-\delta_0}{\delta_c-\delta_0}=T_{\max}$，连续；
- 在 $\delta_{\mathrm{eff}}=\delta_c$ 处：软化段 $T=T_{\max}\cdot\frac{0}{\delta_c-\delta_0}=0$，完全分离段 $T=0$，连续。

**软化段斜率**：
$$
K_{\mathrm{soft}}=\left.\frac{dT}{d\delta_{\mathrm{eff}}}\right|_{\mathrm{soft}}=-\frac{T_{\max}}{\delta_c-\delta_0}=-K\,\frac{\delta_0}{\delta_c-\delta_0}<0 \tag{3.18}
$$

**注 3.2.3（软化段斜率与脆性）**：$|K_{\mathrm{soft}}|$ 越大（即 $\delta_c-\delta_0$ 越小），软化段越陡，材料越脆。极限 $\delta_c\to\delta_0$ 时 $|K_{\mathrm{soft}}|\to\infty$，对应理想脆性断裂（应力瞬时跌落至零）。

**注 3.2.4（TSL 的不可逆性）**：式 (3.16) 描述的是单调加载路径下的 TSL。卸载时（$\dot\delta_{\mathrm{eff}}<0$），traction 沿"指向原点"的线性路径返回（$T=K(1-D_{\mathrm{peak}})\delta_{\mathrm{eff}}$，其中 $D_{\mathrm{peak}}$ 为历史峰值损伤）。再加载时沿同一斜率返回至卸载点，再沿单调 TSL 继续演化。这一不可逆性由 §3.3.2 损伤演化方程严格给出。

#### §3.2.3 界面刚度 $K$ 的正则化选择

**数值要求 3.2.1（弹性段近似刚性）**：弹性段（$\delta_{\mathrm{eff}}<\delta_0$）内，界面应近似"刚性约束"，即 $\delta_{\mathrm{eff}}\ll\delta_0$（理想情况下趋于零）。由式 (3.12)，$T=K\delta_{\mathrm{eff}}$，在给定工作 traction $T_{\mathrm{work}}\sim T_{\max}=K\delta_0$ 下，
$$
\delta_{\mathrm{eff}}\sim\frac{T_{\mathrm{work}}}{K}=\delta_0\,\frac{T_{\mathrm{work}}}{T_{\max}} \tag{3.19}
$$
若 $T_{\mathrm{work}}\sim T_{\max}$（工作 traction 接近峰值），则 $\delta_{\mathrm{eff}}\sim\delta_0$，弹性段位移不可忽略。要使 $\delta_{\mathrm{eff}}\ll\delta_0$，需 $T_{\mathrm{work}}\ll T_{\max}$，即界面实际承载远低于其强度。这一条件通常在界面完好时满足。

**数值要求 3.2.2（避免刚度矩阵病态）**：全局刚度矩阵 $K_{\mathrm{global}}$ 的条件数 $\kappa$ 随界面刚度 $K$ 增大而增大：
$$
\kappa(K_{\mathrm{global}})\propto K \tag{3.20}
$$
$K$ 过大将导致 $K_{\mathrm{global}}$ 条件数过大，Newton-Raphson 迭代收敛困难，数值舍入误差放大。故 $K$ 需在"足够大保证弹性近似"与"不过大避免病态"之间权衡。

**经验选择**：
$$
K\sim\alpha\cdot\frac{E_{\mathrm{eff}}}{\Delta n_{\mathrm{seq}}} \tag{3.21}
$$
其中：
- $\alpha\in[10,100]$ 为正则化因子（经验值）；
- $E_{\mathrm{eff}}$ 为相邻 Q4 单元的有效 Young 模量（§2.6，由 $C_{\mathrm{eff}}^{\mathrm{full}}$ 提取）；
- $\Delta n_{\mathrm{seq}}$ 为 Q4 单元特征厚度（法向尺寸）。

**物理诠释**：式 (3.21) 的物理动机为——界面刚度 $K$ 应远大于"单元等效刚度" $E_{\mathrm{eff}}/\Delta n_{\mathrm{seq}}$（单位厚度单元的刚度量纲为 N/mm$^3$，与 $K$ 同量纲），使界面在弹性段行为接近"刚性粘结"，与 §2.3 Voigt 等应变假设（层间完美粘结）自洽。

**典型量级**（仅给量级）：对锂离子电池层合结构（$E_{\mathrm{eff}}\sim 10-100$ GPa、$\Delta n_{\mathrm{seq}}\sim 10-100\,\mu\mathrm{m}$），由式 (3.21)：
$$
K\sim\alpha\cdot\frac{10^{10}-10^{11}\,\mathrm{Pa}}{10^{-5}-10^{-4}\,\mathrm{m}}\sim 10^{3}-10^{5}\,\mathrm{N/mm^3} \tag{3.22}
$$

**注 3.2.5（$K$ 对数值结果的影响）**：$K$ 的选择影响弹性段行为（应力-应变响应的"刚性程度"），但不影响 $G_c$（由式 3.15 $G_c=\frac{1}{2}K\delta_0\delta_c$，若固定 $G_c$ 与 $T_{\max}=K\delta_0$，则 $K$ 通过 $\delta_0=T_{\max}/K$ 与 $\delta_c=2G_c/(K\delta_0)=2G_c/T_{\max}$ 间接确定）。故 $K$ 是"数值正则化参数"，而非独立物理参数。

**注 3.2.6（$K$ 的收敛性检验）**：严格的数值实施应进行 $K$-收敛性检验——逐步增大 $K$，观察关键输出（如损伤分布、脱粘面积、循环响应）是否收敛。若收敛，取收敛值；若不收敛，需重新评估 TSL 参数标定。

#### §3.2.4 法向-切向解耦（本文简化）

**假设 3.2.2（法向-切向解耦）**：本文采用法向与切向独立 TSL 的简化，即法向 traction $T_n$ 仅依赖 $\delta_n$，切向 traction $T_s$ 仅依赖 $\delta_s$，二者各自有独立的 TSL 参数：
$$
T_n=K_n\,(1-D_n)\,\delta_n,\qquad T_s=K_s\,(1-D_s)\,\delta_s \tag{3.23}
$$

其中：
- $K_n,K_s$：法向、切向界面初始刚度；
- $\delta_{n,0},\delta_{s,0}$：法向、切向损伤起始分离；
- $\delta_{n,c},\delta_{s,c}$：法向、切向临界分离；
- $D_n,D_s$：法向、切向损伤变量。

**简化动机**：
1. **参数标定简单**：法向与切向 TSL 参数可独立由 mode I（DCB 试验）与 mode II（ENF 试验）标定；
2. **数值实施模块化**：法向与切向 traction 计算可分别独立进行，无需混合模式耦合（耦合通过 §3.3 统一损伤 $D$ 实现）；
3. **物理可解释性**：$D_n$ 描述张开损伤，$D_s$ 描述剪切损伤，物理意义明确。

**法向 TSL 参数关系**（类似式 3.15）：
$$
G_{Ic}=\frac{1}{2}K_n\,\delta_{n,0}\,\delta_{n,c},\qquad T_{n,\max}=K_n\,\delta_{n,0} \tag{3.24}
$$

**切向 TSL 参数关系**：
$$
G_{IIc}=\frac{1}{2}K_s\,\delta_{s,0}\,\delta_{s,c},\qquad T_{s,\max}=K_s\,\delta_{s,0} \tag{3.25}
$$

**注 3.2.7（$K_n\neq K_s$ 的物理来源）**：法向与切向界面刚度通常不同（$K_n\neq K_s$），反映界面对张开与剪切的抵抗能力差异。典型地，$K_n>K_s$（界面抗张开能力强于抗剪切）。本文默认取 $K_n=K_s=K$（等刚度假设），除非特殊说明。

**注 3.2.8（压缩法向的处理）**：式 (3.23) 中 $\delta_n<0$（界面闭合压缩）时，不应产生损伤（压缩不引起脱粘）。本文取
$$
T_n=K_n\,\delta_n\quad(\delta_n<0,\,D_n=0\text{ 锁定}) \tag{3.26}
$$
即压缩段行为为纯弹性（无损伤），$T_n<0$ 对应界面压力。这一处理避免压缩引起的虚假损伤。

**推论 3.2.2（解耦 TSL 的等效 traction）**：在解耦假设下，等效 traction（式 3.11）为
$$
T=\sqrt{[K_n(1-D_n)\delta_n]^2+[K_s(1-D_s)\delta_s]^2} \tag{3.27}
$$
若取 $K_n=K_s=K$ 且 $D_n=D_s=D$（统一损伤，§3.3 给出），则式 (3.27) 退化为式 (3.13)：$T=K(1-D)\delta_{\mathrm{eff}}$。

**统一损伤 $D$ 的引入**：式 (3.23) 中 $D_n,D_s$ 独立演化虽物理清晰，但混合模式耦合需要统一损伤标量 $D$。本文采用统一损伤演化（§3.3），即 $D_n=D_s=D$，由混合模式准则（§3.3.4 B-K 准则）驱动。这一选择牺牲部分物理细节（无法区分法向与切向损伤累积差异），换取数值实施简化与参数标定便利。

---

### §3.3 损伤变量与演化方程

本节定义 CZM 损伤标量 $D\in[0,1]$，建立其演化方程（双线性 TSL 对应），给出等效应变能释放率 $G$ 的解析表达，并引入 Benzeggagh-Kenane（B-K）混合模式断裂准则。最后简述 Park-Paulino（P-P）unified TSL 作为可选替代。

#### §3.3.1 损伤变量定义

**定义 3.3.1（内聚力损伤标量）**：定义 CZM 损伤标量为
$$
D\in[0,1] \tag{3.28}
$$
满足：
- $D=0$：界面完好（弹性段，式 3.12），无损伤累积；
- $D\in(0,1)$：部分损伤（软化段，式 3.13），有效刚度折减为 $K(1-D)$；
- $D=1$：完全脱粘（完全分离段，式 3.14），traction 消失。

**初始条件**：
$$
D\bigl|_{t=0}=0 \tag{3.29}
$$
即假设初始时刻界面完好（无预损伤）。

**注 3.3.1（$D$ 的物理诠释）**：$D$ 表征界面面积的"损伤分数"——$D=0$ 对应全部面积有效粘结，$D=1$ 对应全部面积脱粘。这一诠释与 §2.8.2 式 (2.68) $A_{\mathrm{eff}}=A_0(1-D)$ 一致：$D$ 等于脱粘面积分数。

**注 3.3.2（$D$ 的标量性质，Curie 兼容）**：$D$ 为零阶张量（标量），由 §2.8.2 注 2.8.5 验证。这一标量性质是 $D$ 作为 Curie 中介（矢量 traction $\to$ 标量 $D$ $\to$ 标量电化学/热反馈，式 2.70）的必要条件。

**注 3.3.3（$D$ 与 $D_n,D_s$ 的关系）**：在 §3.2.4 解耦 TSL 假设下，$D_n,D_s$ 独立演化。本文采用统一损伤假设（$D_n=D_s=D$，推论 3.2.2），由混合模式准则驱动（§3.3.4）。统一损伤 $D$ 描述"界面整体损伤状态"，不区分模式细节。

#### §3.3.2 演化方程（双线性 TSL 对应）

**定义 3.3.2（双线性 TSL 损伤演化律，Camanho 形式）**：在双线性 TSL（§3.2.1）框架下，为确保软化段 $T=K(1-D)\delta_{\mathrm{eff}}$（式 3.13）严格退化到线性包络 $T=T_{\max}(\delta_c-\delta_{\mathrm{eff}})/(\delta_c-\delta_0)$（式 3.16），损伤变量 $D$ 必须采用 **Camanho（2002）形式**（见 §3.3.2 末注 3.3.6），其演化方程为
$$
\dot D=\begin{cases}
\dot\delta_{\mathrm{eff}}\cdot H(\delta_{\mathrm{eff}}-\delta_0)\cdot H(\delta_c-\delta_{\mathrm{eff}})\cdot\dfrac{\delta_c\,\delta_0}{\delta_{\mathrm{eff}}^{\,2}\,(\delta_c-\delta_0)} & \dot\delta_{\mathrm{eff}}>0\\[6pt]
0 & \dot\delta_{\mathrm{eff}}\leq 0
\end{cases} \tag{3.30}
$$
其中：
- $H(\cdot)$ 为 Heaviside 函数（$H(x)=1$ 若 $x\geq 0$，$H(x)=0$ 若 $x<0$）；
- $\dot\delta_{\mathrm{eff}}=d\delta_{\mathrm{eff}}/dt$ 为等效分离率；
- $\delta_0,\delta_c$ 为 §3.2.1 定义的特征分离；
- $H(\delta_c-\delta_{\mathrm{eff}})$ 保证 $\delta_{\mathrm{eff}}\geq\delta_c$ 时 $D$ 锁定为 $1$（与式 3.31 第三段一致）。

> **[v2.3 修订说明]**：原 v2.2 版采用线性 $D=(\delta_{\mathrm{eff}}-\delta_0)/(\delta_c-\delta_0)$，代入式 (3.13) 得 $T=K\delta_{\mathrm{eff}}(\delta_c-\delta_{\mathrm{eff}})/(\delta_c-\delta_0)$——这是关于 $\delta_{\mathrm{eff}}$ 的**抛物线**（$\delta_c/2$ 处取极大），与 (3.16) 的线性软化矛盾。Camanho（2002）形式 $D=\delta_c(\delta_{\mathrm{eff}}-\delta_0)/[\delta_{\mathrm{eff}}(\delta_c-\delta_0)]$ 严格满足 (3.13)→(3.16) 的代数一致性（验证见注 3.3.6）。本修订联动 §3.4.3 的 $G(D)$ 推导（式 3.58–3.62）。

**积分形式（单调加载）**：在单调加载（$\dot\delta_{\mathrm{eff}}>0$ 持续）下，积分式 (3.30) 得
$$
D(\delta_{\mathrm{eff}})=\begin{cases}
0 & 0\leq\delta_{\mathrm{eff}}<\delta_0\\[4pt]
\dfrac{\delta_c\,(\delta_{\mathrm{eff}}-\delta_0)}{\delta_{\mathrm{eff}}\,(\delta_c-\delta_0)} & \delta_0\leq\delta_{\mathrm{eff}}<\delta_c\\[6pt]
1 & \delta_{\mathrm{eff}}\geq\delta_c
\end{cases} \tag{3.31}
$$

**物理诠释**：
1. **损伤起始判据**：仅当 $\delta_{\mathrm{eff}}$ 超过阈值 $\delta_0$ 时损伤启动（由 $H(\delta_{\mathrm{eff}}-\delta_0)$ 控制）。$\delta_{\mathrm{eff}}<\delta_0$ 时 $D=0$，界面完好；
2. **损伤累积**：$\delta_0\leq\delta_{\mathrm{eff}}<\delta_c$ 时，$D$ 与 $\delta_{\mathrm{eff}}$ 呈 **Camanho 非线性**关系（凸增——$\delta_{\mathrm{eff}}\to\delta_0^+$ 时 $D\sim(\delta_{\mathrm{eff}}-\delta_0)/\delta_0$ 起步平缓，$\delta_{\mathrm{eff}}\to\delta_c^-$ 时 $D\to 1$ 急剧上升）。$D$ 对 $\delta_{\mathrm{eff}}$ 的导数 $dD/d\delta_{\mathrm{eff}}=\delta_c\delta_0/[\delta_{\mathrm{eff}}^2(\delta_c-\delta_0)]$ 随 $\delta_{\mathrm{eff}}$ 单调递减；
3. **完全脱粘**：$\delta_{\mathrm{eff}}\geq\delta_c$ 时 $D=1$，锁定（后续不再增加，因 $D\leq 1$ 约束）。

**注 3.3.6（Camanho 形式与 (3.13)→(3.16) 一致性验证）**：将式 (3.31) 第二段代入式 (3.13) $T=K(1-D)\delta_{\mathrm{eff}}$：
$$
1-D\overset{(3.31)}{=}\frac{\delta_{\mathrm{eff}}(\delta_c-\delta_0)-\delta_c(\delta_{\mathrm{eff}}-\delta_0)}{\delta_{\mathrm{eff}}(\delta_c-\delta_0)}=\frac{\delta_0(\delta_c-\delta_{\mathrm{eff}})}{\delta_{\mathrm{eff}}(\delta_c-\delta_0)}
$$
故
$$
T=K\cdot\dfrac{\delta_0(\delta_c-\delta_{\mathrm{eff}})}{\delta_{\mathrm{eff}}(\delta_c-\delta_0)}\cdot\delta_{\mathrm{eff}}=K\delta_0\cdot\dfrac{\delta_c-\delta_{\mathrm{eff}}}{\delta_c-\delta_0}=T_{\max}\cdot\dfrac{\delta_c-\delta_{\mathrm{eff}}}{\delta_c-\delta_0}\;\checkmark
$$
严格恢复式 (3.16) 第二段（线性软化）。端点检验：$\delta_{\mathrm{eff}}=\delta_0\Rightarrow D=0,\;T=T_{\max}$；$\delta_{\mathrm{eff}}=\delta_c\Rightarrow D=1,\;T=0$。$\square$

**注 3.3.7（Camanho 退化检验）**：$\delta_{\mathrm{eff}}\to\delta_0^+$ 时 $D\to 0$ 且 $dD/d\delta_{\mathrm{eff}}\to 1/\delta_0$（斜率有限）；$\delta_{\mathrm{eff}}\to\delta_c^-$ 时 $D\to 1$ 且 $dD/d\delta_{\mathrm{eff}}\to\delta_0/[\delta_c(\delta_c-\delta_0)]$（有限）。整段 $D(\delta_{\mathrm{eff}})\in[0,1]$ 单调光滑。

**不可逆性**：式 (3.30) 第二行（$\dot\delta_{\mathrm{eff}}\leq 0$ 时 $\dot D=0$）保证损伤不可逆——卸载时 $D$ 保持历史峰值，不会恢复。这是 CZM 的核心物理特性，描述界面损伤的不可愈合性。

> **[v2.3 修订说明 B6 / (3.30) 历史阈值显式化]**：
> - **原 v2.2 表述**：(3.30) 仅以 $\dot\delta_{\mathrm{eff}}\gtrless 0$ 切换控制 $\dot D$，**缺历史阈值 $\delta_{\mathrm{peak}}$ 显式门控**。
> - **不一致诊断**：(3.30) 形式上等价于 $D(\delta_{\mathrm{eff}})=D_{\mathrm{Camanho}}(\delta_{\mathrm{eff}})$（单调加载）或 $D=D^n$（$\dot\delta_{\mathrm{eff}}\leq 0$），但**未显式表达**再加载路径的"通过 $\delta_{\mathrm{peak}}$ 前 $\dot D=0$"条件。与紧随其后的"卸载-再加载路径"叙述（$\delta_{\mathrm{eff}}<\delta_{\mathrm{peak}}$ 时 $D$ 不变，沿 (3.32) 至 $\delta_{\mathrm{peak}}$）以及 §6.4 (6.38') max-history 代数更新形式不一致——后者用 $\max(D^n,\widetilde D(\delta_{\mathrm{eff}}))$ 严格实现历史阈值。审查 findings B6 命中。
> - **修正方案**：在 (3.30) 中显式加入历史峰值 Heaviside 门 $H(\delta_{\mathrm{eff}}-\delta_{\mathrm{peak}}(t))$（其中 $\delta_{\mathrm{peak}}(t)=\max_{\tau\leq t}\delta_{\mathrm{eff}}(\tau)$ 为历史最大有效分离）：
>
>   $$\dot D=\dot\delta_{\mathrm{eff}}\cdot H(\dot\delta_{\mathrm{eff}})\cdot H(\delta_{\mathrm{eff}}-\delta_{\mathrm{peak}}(t))\cdot H(\delta_{\mathrm{eff}}-\delta_0)\cdot H(\delta_c-\delta_{\mathrm{eff}})\cdot\dfrac{\delta_c\,\delta_0}{\delta_{\mathrm{eff}}^{\,2}\,(\delta_c-\delta_0)} \tag{3.30'}$$
>
>   即：损伤演化需同时满足 (i) $\dot\delta_{\mathrm{eff}}>0$（加载中）、(ii) $\delta_{\mathrm{eff}}>\delta_{\mathrm{peak}}$（突破历史峰值）、(iii) $\delta_{\mathrm{eff}}\in[\delta_0,\delta_c]$（软化段内）。$\delta_{\mathrm{eff}}>\delta_c$ 时 $D\equiv 1$ 锁定。
> - **max-history 代数等价**：(3.30') 与 §6.4 (6.38') 的 max-history 代数更新 $D^{n+1}=\max(D^n,\widetilde D(\delta_{\mathrm{eff}}^{n+1}))$ 严格等价（在分段单调加载下，微分形式与代数形式互为重言）；本节保留 (3.30) 原表达作为**简化版**（隐含历史峰值），(3.30') 作为**严格版**。
> - **联动效应**：与 D1 修订（§6.4 改为 max-history 代数更新）严格一致；§3.4 $G(D)$ 推导不受影响（$G$ 是 $D$ 的代数函数，与演化路径无关）。

**卸载-再加载路径**：
- **卸载**（$\dot\delta_{\mathrm{eff}}<0$）：$D$ 保持当前值 $D_{\mathrm{peak}}$，traction 沿"指向原点"的线性路径返回
  $$
  T=K(1-D_{\mathrm{peak}})\,\delta_{\mathrm{eff}} \tag{3.32}
  $$
- **再加载**（$\dot\delta_{\mathrm{eff}}>0$ 但 $\delta_{\mathrm{eff}}<\delta_{\mathrm{peak}}$）：$H(\delta_{\mathrm{eff}}-\delta_{\mathrm{peak}})=0$ → $\dot D=0$（(3.30') 严格实现），$D=D_{\mathrm{peak}}$ 保持，沿同一斜率式 (3.32) 上升，直至达到历史峰值点 $(\delta_{\mathrm{peak}},T_{\mathrm{peak}})$；
- **继续加载**（$\delta_{\mathrm{eff}}>\delta_{\mathrm{peak}}$）：$H(\delta_{\mathrm{eff}}-\delta_{\mathrm{peak}})=1$ → $\dot D$ 按 (3.30') 演化，沿单调 TSL（式 3.16）继续，$D$ 继续增加。

**注 3.3.4（不可逆性与熵不等式）**：损伤演化（式 3.30）的不可逆性（$\dot D\geq 0$）保证 CZM 满足热力学第二定律（Clausius-Duhem 不等式）。损伤耗散功率为
$$
\mathcal{D}_{\mathrm{czm}}=-\frac{\partial G}{\partial D}\,\dot D\geq 0 \tag{3.33}
$$
其中 $G$ 为应变能释放率（§3.3.3）。在双线性 TSL 框架下，$\partial G/\partial D<0$（损伤增加释放能量），故 $\mathcal{D}_{\mathrm{czm}}\geq 0$ 自动满足。

**注 3.3.5（$\dot D$ 的数值实现）**：式 (3.30) 中 Heaviside 函数 $H(\cdot)$ 在数值实现中需正则化（避免不连续导致的 Newton-Raphson 收敛困难）。典型正则化为
$$
H_\epsilon(x)=\frac{1}{2}\left(1+\tanh\frac{x}{\epsilon}\right) \tag{3.34}
$$
其中 $\epsilon\ll\delta_0$ 为正则化参数。

#### §3.3.3 等效应变能释放率

**定义 3.3.3（等效应变能释放率）**：定义等效应变能释放率为
$$
G(\delta_{\mathrm{eff}})=\int_0^{\delta_{\mathrm{eff}}}T(\delta')\,d\delta' \tag{3.35}
$$
即 TSL 曲线从 $0$ 到 $\delta_{\mathrm{eff}}$ 的积分面积（单位 J/m$^2$）。

**双线性 TSL 的解析表达**：

**弹性段**（$0\leq\delta_{\mathrm{eff}}<\delta_0$）：
$$
G=\int_0^{\delta_{\mathrm{eff}}}K\delta'\,d\delta'=\frac{1}{2}K\,\delta_{\mathrm{eff}}^2 \tag{3.36}
$$
此段 $G$ 为弹性储能，可逆（卸载时完全释放）。

**软化段**（$\delta_0\leq\delta_{\mathrm{eff}}<\delta_c$）：
$$
G=\frac{1}{2}K\delta_0^2+\int_{\delta_0}^{\delta_{\mathrm{eff}}}T_{\max}\frac{\delta_c-\delta'}{\delta_c-\delta_0}\,d\delta' \tag{3.37}
$$
展开得
$$
G=\frac{1}{2}K\delta_0^2+\frac{T_{\max}}{\delta_c-\delta_0}\left[\delta_c(\delta_{\mathrm{eff}}-\delta_0)-\frac{1}{2}(\delta_{\mathrm{eff}}^2-\delta_0^2)\right]
$$
利用 $T_{\max}=K\delta_0$，整理为
$$
G=G_{\mathrm{init}}+\frac{T_{\max}}{2(\delta_c-\delta_0)}\left[2\delta_c(\delta_{\mathrm{eff}}-\delta_0)-(\delta_{\mathrm{eff}}^2-\delta_0^2)\right] \tag{3.38}
$$
其中 $G_{\mathrm{init}}=\frac{1}{2}K\delta_0^2=\frac{1}{2}T_{\max}\delta_0$ 为损伤起始时的能量释放率。

此段 $G$ 包含可逆部分（弹性储能）与不可逆部分（损伤耗散），单调递增。

**完全分离段**（$\delta_{\mathrm{eff}}\geq\delta_c$）：
$$
G=G_c=\frac{1}{2}T_{\max}\delta_c \tag{3.39}
$$
此段 $G$ 达到临界值 $G_c$（式 3.15），界面完全脱粘。

**临界断裂韧度判据**：
$$
G=G_c\quad\Longleftrightarrow\quad \delta_{\mathrm{eff}}=\delta_c\quad\Longleftrightarrow\quad D=1 \tag{3.40}
$$

三个等价条件描述同一物理状态——完全脱粘。

**推论 3.3.1（$G$ 与 $D$ 的关系）**：在软化段（$\delta_0\leq\delta_{\mathrm{eff}}<\delta_c$），由式 (3.31)（Camanho 形式）$\delta_{\mathrm{eff}}(D)=\delta_c\delta_0/[\delta_c-D(\delta_c-\delta_0)]$（即式 3.59），代入式 (3.38) 的 $G(\delta_{\mathrm{eff}})$ 闭式可得 $G$ 关于 $D$ 的显式表达。在 Camanho 双线性 TSL 框架下，$G(D)$ 为 $D$ 的非线性函数（凸增），具体形式见 §3.4.3 式 (3.62'')。

**注 3.3.6（$G$ 与 $D$ 的物理区别）**：
- $G$（式 3.35）：累计能量释放率，从 $0$ 增至 $G_c$，描述"已消耗能量"；
- $D$（式 3.28）：损伤标量，从 $0$ 增至 $1$，描述"损伤程度"；
- 二者一一对应（在单调加载下），但物理意义不同——$G$ 为能量量，$D$ 为无量纲标量。

**注 3.3.7（$G_I,G_{II}$ 的分量分解）**：在解耦 TSL（§3.2.4）框架下，等效应变能释放率可分解为 mode I 与 mode II 分量：
$$
G=G_I+G_{II},\qquad G_I=\int_0^{\delta_n}T_n\,d\delta_n',\quad G_{II}=\int_0^{\delta_s}T_s\,d\delta_s' \tag{3.41}
$$
$G_I,G_{II}$ 分别描述 mode I（张开）与 mode II（剪切）对总能量释放率的贡献。

#### §3.3.4 混合模式损伤演化（Benzeggagh-Kenane 准则）

在混合模式加载下（$\delta_n\neq 0,\delta_s\neq 0$，§3.1.3 定义 3.1.3），界面损伤由 mode I 与 mode II 共同驱动。本小节引入 Benzeggagh-Kenane（B-K）准则，给出混合模式临界断裂韧度 $G_c^{\mathrm{mix}}$ 的经验公式。

**定义 3.3.4（B-K 混合模式断裂准则）**：B-K 准则给出混合模式临界断裂韧度为
$$
G_c^{\mathrm{mix}}=G_{Ic}+(G_{IIc}-G_{Ic})\left(\frac{G_{II}}{G_I+G_{II}}\right)^\eta \tag{3.42}
$$
其中：
- $G_I,G_{II}$：mode I、mode II 当前能量释放率（式 3.41）；
- $G_{Ic},G_{IIc}$：mode I、mode II 临界断裂韧度（纯模式参数，式 3.24–3.25）；
- $\eta$：B-K 经验指数，控制 $G_c^{\mathrm{mix}}$ 随混合模式比的非线性变化。

**参数物理意义**：
- $G_{Ic}$：纯 mode I 临界断裂韧度（$\beta_{\mathrm{mix}}=0$ 时 $G_c^{\mathrm{mix}}=G_{Ic}$）；
- $G_{IIc}$：纯 mode II 临界断裂韧度（$\beta_{\mathrm{mix}}=1$ 时 $G_c^{\mathrm{mix}}=G_{IIc}$）；
- $\eta$：经验指数，典型 $\eta\approx 1-2$（由混合模式断裂试验标定）。

**与 §3.1.3 的对应**：由式 (3.10) $\beta_{\mathrm{mix}}=G_{II}/(G_I+G_{II})$（等刚度假设下），式 (3.42) 可改写为
$$
G_c^{\mathrm{mix}}(\beta_{\mathrm{mix}})=G_{Ic}+(G_{IIc}-G_{Ic})\,\beta_{\mathrm{mix}}^\eta \tag{3.43}
$$

**$G_c^{\mathrm{mix}}$ 的行为**：
- $\beta_{\mathrm{mix}}=0$（纯 mode I）：$G_c^{\mathrm{mix}}=G_{Ic}$；
- $\beta_{\mathrm{mix}}=1$（纯 mode II）：$G_c^{\mathrm{mix}}=G_{IIc}$；
- $\beta_{\mathrm{mix}}\in(0,1)$：$G_c^{\mathrm{mix}}$ 在 $[G_{Ic},G_{IIc}]$ 之间单调变化，变化速率由 $\eta$ 控制：
  - $\eta=1$：线性插值 $G_c^{\mathrm{mix}}=G_{Ic}+(G_{IIc}-G_{Ic})\beta_{\mathrm{mix}}$；
  - $\eta>1$：$G_c^{\mathrm{mix}}$ 在 $\beta_{\mathrm{mix}}$ 较小时变化缓慢（接近 $G_{Ic}$），$\beta_{\mathrm{mix}}\to 1$ 时快速增至 $G_{IIc}$；
  - $\eta<1$：$G_c^{\mathrm{mix}}$ 在 $\beta_{\mathrm{mix}}$ 较小时快速增至中间值。

**混合模式损伤驱动判据**：当总能量释放率 $G=G_I+G_{II}$ 达到混合模式临界值时，损伤达到完全脱粘：
$$
G=G_c^{\mathrm{mix}}\quad\Longleftrightarrow\quad D=1 \tag{3.44}
$$

这一判据替代纯模式判据（式 3.40），适用于混合模式加载。

**注 3.3.8（B-K 准则的经验性质）**：B-K 准则（式 3.42）为经验公式，非严格理论推导。其优势在于：
1. 形式简单，仅含两个纯模式参数 $(G_{Ic},G_{IIc})$ 与一个经验指数 $\eta$；
2. 在 $\beta_{\mathrm{mix}}=0,1$ 处严格退化为纯模式，与纯模式试验结果自洽；
3. 适用于多种复合材料界面，经验证精度可接受。

**替代准则（仅作说明，本文不采用）**：
- **幂律准则**（power law）：$(G_I/G_{Ic})^m+(G_{II}/G_{IIc})^n=1$，含两个经验指数 $m,n$；
- **线性准则**：$G_I/G_{Ic}+G_{II}/G_{IIc}=1$，B-K 取 $\eta=1$ 的退化；
- **二次准则**：$(G_I/G_{Ic})^2+(G_{II}/G_{IIc})^2=1$。

本文默认 B-K 准则（式 3.42），$\eta=1.5$（典型值，由文献推荐）。

**注 3.3.9（$\eta$ 的标定）**：$\eta$ 需由混合模式断裂试验（如 MMB 试验，mixed-mode bending）标定。若无试验数据，可取 $\eta=1-2$ 的典型值。$\eta$ 对 $G_c^{\mathrm{mix}}$ 的影响在 $\beta_{\mathrm{mix}}\to 0$ 或 $\beta_{\mathrm{mix}}\to 1$ 时较小（因 $G_c^{\mathrm{mix}}\to G_{Ic}$ 或 $G_{IIc}$），在 $\beta_{\mathrm{mix}}\approx 0.5$ 时最敏感。

**推论 3.3.2（B-K 准则下的统一损伤演化）**：在 B-K 准则下，统一损伤 $D$（$D_n=D_s=D$，§3.2.4 推论 3.2.2）的演化律为
$$
D=\begin{cases}
0 & G<G_{\mathrm{init}}^{\mathrm{mix}}\\[4pt]
\dfrac{G-G_{\mathrm{init}}^{\mathrm{mix}}}{G_c^{\mathrm{mix}}-G_{\mathrm{init}}^{\mathrm{mix}}} & G_{\mathrm{init}}^{\mathrm{mix}}\leq G<G_c^{\mathrm{mix}}\\[6pt]
1 & G\geq G_c^{\mathrm{mix}}
\end{cases} \tag{3.45}
$$
其中 $G_{\mathrm{init}}^{\mathrm{mix}}$ 为混合模式损伤起始能量（对应 $\delta_{\mathrm{eff}}=\delta_0$ 时的 $G$ 值，由式 3.36 给出 $G_{\mathrm{init}}^{\mathrm{mix}}=\frac{1}{2}K\delta_0^2$）。

式 (3.45) 为式 (3.31) 在能量空间的重写，二者在双线性 TSL 下等价。

#### §3.3.5 Park-Paulino（P-P）unified TSL（可选替代）

> **[可选内容]**：本小节简述 Park-Paulino unified TSL 的符号形式，作为本文双线性 TSL 的替代选项。本文默认不采用 P-P TSL，此处仅作记录。

**P-P TSL 特点**：

1. **统一参数描述**：P-P TSL 用统一参数 $(T_{\max},G_c,K)$ 描述整个 TSL（包括弹性段、软化段、完全分离段），无需分段定义；
2. **平滑过渡（$C^1$ 连续）**：P-P TSL 在弹性段与软化段交界处一阶导数连续（$C^1$ 连续），数值稳定性优于双线性 TSL（双线性 TSL 仅 $C^0$ 连续，在 $\delta_{\mathrm{eff}}=\delta_0$ 处斜率不连续）；
3. **适用于循环加载**：P-P TSL 的卸载-再加载路径内置定义，适合循环加载场景（双线性 TSL 需额外定义卸载路径，式 3.32）。

**符号形式**：
$$
T(\delta;\,T_{\max},G_c,K)=T_{\max}\cdot f(\delta;\,G_c,K) \tag{3.46}
$$
其中 $f(\cdot)$ 为 P-P 形函数，具体形式见 Park & Paulino (2011)。

**P-P TSL 参数与双线性 TSL 参数的对应**：
- $T_{\max}$（峰值 traction）：与双线性 TSL 的 $T_{\max}=K\delta_0$ 一致；
- $G_c$（临界断裂韧度）：与双线性 TSL 的 $G_c=\frac{1}{2}K\delta_0\delta_c$（式 3.15）一致；
- $K$（初始刚度）：与双线性 TSL 的 $K$ 一致。

**本文选择**：默认双线性 TSL（§3.2）——计算效率高、参数标定简单、$C^0$ 连续在大多数静力分析中可接受。P-P TSL 作为可选替代，建议用于：
- 循环加载场景（需准确描述卸载-再加载路径）；
- 高精度收敛性要求（$C^1$ 连续减少 Newton-Raphson 迭代步数）；
- 与现有代码库兼容（若已有 P-P TSL 实现）。

**注 3.3.10（P-P TSL 参数库）**：知识库 `03_内聚力模型CZM.md`（若存在）应给出 P-P TSL 形函数 $f$ 的具体形式、参数标定流程与数值实现细节。若该知识库不存在，注"待补充"。

#### §3.3.6 接口传递

本小节给出 §3.1–§3.3 的输出接口，供下游章节（§3.4–§3.6）调用。

**输出 1（损伤标量 $D$）**：
- 定义：式 (3.28)–(3.29)，演化由式 (3.30) 或 (3.45) 给出；
- 取值范围：$D\in[0,1]$；
- 下游用途：
  - §3.4 损伤驱动力分析（与 $G,G_c^{\mathrm{mix}}$ 联合判断损伤状态）；
  - §3.6 损伤下游反馈（驱动 $R_{\mathrm{contact}},A_{\mathrm{eff}},k_n^{\mathrm{eff}}$，式 2.67–2.69）；
  - 作为 Curie 中介标量，桥接矢量 traction 与标量电化学/热方程（式 2.70）。

**输出 2（CZM 输出 traction $T_n,T_s$）**：
- 定义：式 (3.23)（解耦 TSL）或式 (3.16)（等效 TSL）；
- 下游用途：
  - §3.5 与 §2.7.3 输入 traction（$T_n,T_s$ 由残余应力 + 外载决定）的闭环检验——CZM 输出 traction 应与 Q4 单元界面应力一致（位移协调条件）；
  - 弱形式中界面贡献项（内聚力虚功）。

**输出 3（能量释放率 $G$ 与 $G_c^{\mathrm{mix}}$）**：
- 定义：$G$ 由式 (3.35) 给出，$G_c^{\mathrm{mix}}$ 由 B-K 准则式 (3.42) 给出；
- 下游用途：
  - §3.4 应变能释放率分析（$G$ 与 $G_c^{\mathrm{mix}}$ 的比较，判断损伤是否达到临界）；
  - 循环演化分析（$G$ 的逐周演化，疲劳累积准则）。

**输出 4（混合模式比 $\beta_{\mathrm{mix}}$）**：
- 定义：式 (3.9)；
- 下游用途：
  - §3.3.4 B-K 准则（式 3.43），决定 $G_c^{\mathrm{mix}}$；
  - 损伤模式诊断（区分 mode I 主导与 mode II 主导，物理可解释性）。

**接口闭环示意**：
```
§2.7.3 输入 T_n,T_s（来自残余应力 + 外载）
        ↓
§3.1 界面运动学：δ_n,δ_s → δ_eff, β_mix
        ↓
§3.2 TSL：T(δ_eff)（双线性，分段定义）
        ↓
§3.3 损伤演化：D（演化方程 3.30/3.45）
        ↓                          ↓
§3.4 损伤驱动力分析        §3.6 损伤下游反馈
（G vs G_c^mix）           （R_contact, A_eff, k_r^eff）
                            ↓
                    §4 SPMe / §5 热方程（标量反馈，Curie 兼容）
```

**与已有记号的闭环检验（R8 自查）**：
- $T_n,T_s$：§2.7.3 定义，§3.1 引用作为输入——无重复定义；
- $\delta_n,\delta_s$：§1.5 定义，§3.1.1 引用（式 3.1–3.2）——无重复定义；
- $D$：§0.5 Curie 注释引入，§3.3.1 严格定义（式 3.28）——补充定义，无冲突；
- $G_c$：§3.2.1 定义（式 3.15），§3.3.3 引用——无重复定义；
- $G_{Ic},G_{IIc}$：§3.2.4 式 (3.24)–(3.25) 定义，§3.3.4 引用——无重复定义。

**注 3.3.11（下游反馈的 Curie 兼容性，与 §2.8.2 闭环）**：本节输出的损伤标量 $D$（式 3.28）为标量（零阶张量），通过式 (2.67)–(2.69) 驱动电化学/热反馈（均为标量），严格满足 Curie 原理（§2.8.2 定理 2.8.1）。这一 Curie 兼容性是 §3.6 损伤下游反馈的理论基础。

**注 3.3.12（参数库引用）**：CZM 参数 $(K,\delta_0,\delta_c,G_{Ic},G_{IIc},\eta)$ 的具体取值应来自知识库 `03_内聚力模型CZM.md`。若该知识库不存在或参数缺失，标注"待补充"。本章仅给出符号推导，不代入具体数值（量级估计标注"仅给量级"，见式 3.22）。

## 3.4 损伤驱动力与应变能释放率

§3.3 建立了损伤变量 $D$ 的演化方程（式 3.30 / 3.45）与混合模式断裂韧度 $G_c^{\mathrm{mix}}$ 的 B-K 准则（式 3.42）。本节从热力学共轭角度严格定义应变能释放率 $G$ 与损伤驱动力 $Y$，建立 $D\leftrightarrow G$ 的双向映射，并给出能量判据相对应力判据的适用范围论证。本节为 §3.5 闭环（Q4 残余应力 → CZM 牵引 → $\delta_{\mathrm{eff}}$ → $D$）提供能量判据基础。

### 3.4.1 应变能释放率的严格定义

**定义 3.4.1（应变能释放率）**：对内聚力界面 $\Gamma_{\mathrm{coh}}$ 上一Material 点，给定加载路径 $\delta_{\mathrm{eff}}:0\to\delta_{\mathrm{eff}}^*$，应变能释放率 $G$ 定义为单位面积界面分离所消耗的累计功：

$$
G(\delta_{\mathrm{eff}}^*)\;:=\;\int_{0}^{\delta_{\mathrm{eff}}^*} T(\delta_{\mathrm{eff}}')\,\mathrm{d}\delta_{\mathrm{eff}}'\qquad (3.47)
$$

其中 $T(\delta_{\mathrm{eff}})$ 为 §3.2.1 式 (3.16) 给出的双线性 TSL 包络（有效牵引）。

**物理诠释**：$G$ 即 TSL 曲线在 $(\delta_{\mathrm{eff}},T)$ 平面下所围的面积，量纲为 $\mathrm{J}\cdot\mathrm{m}^{-2}$，表征界面单位面积累计耗散的能量（含可逆弹性储能与不可逆损伤耗散两部分）。

**分段闭式表达（双线性 TSL）**：由式 (3.16)–(3.18) 分段积分得：

(i) **弹性段**（$0\leq\delta_{\mathrm{eff}}<\delta_0$，$D=0$）：

$$
G(\delta_{\mathrm{eff}})\;=\;\int_{0}^{\delta_{\mathrm{eff}}} K\delta'\,\mathrm{d}\delta'\;=\;\tfrac{1}{2}K\delta_{\mathrm{eff}}^{2}\qquad (3.48)
$$

此段 $G$ 完全可逆（卸载回原点无耗散）。

(ii) **软化段**（$\delta_0\leq\delta_{\mathrm{eff}}<\delta_c$，$D\in(0,1)$）：

$$
G(\delta_{\mathrm{eff}})\;=\;\underbrace{\tfrac{1}{2}K\delta_0^{2}}_{\text{弹性储能}}\;+\;\int_{\delta_0}^{\delta_{\mathrm{eff}}} T(\delta')\,\mathrm{d}\delta'\qquad (3.49)
$$

代入软化段线性 TSL（式 3.17）$T(\delta')=T_{\max}\dfrac{\delta_c-\delta'}{\delta_c-\delta_0}$，得：

$$
G(\delta_{\mathrm{eff}})\;=\;\tfrac{1}{2}K\delta_0^{2}\;+\;\tfrac{T_{\max}}{2}\bigl(\delta_{\mathrm{eff}}-\delta_0\bigr)\!\left[\dfrac{\delta_c-\delta_0}{\delta_c-\delta_0}+\dfrac{\delta_c-\delta_{\mathrm{eff}}}{\delta_c-\delta_0}\right]
$$

化简为：

$$
G(\delta_{\mathrm{eff}})\;=\;\tfrac{1}{2}K\delta_0^{2}\;+\;\tfrac{1}{2}\bigl(T_{\max}+T(\delta_{\mathrm{eff}})\bigr)\bigl(\delta_{\mathrm{eff}}-\delta_0\bigr)\qquad (3.50)
$$

即软化段累计 $G$ 为初始弹性储能加梯形面积（线性软化几何意义）。

(iii) **完全分离**（$\delta_{\mathrm{eff}}\geq\delta_c$，$D=1$）：

$$
G\;=\;G_c\;:=\;\tfrac{1}{2}K\delta_0\,\delta_c\;=\;\tfrac{1}{2}T_{\max}\,\delta_c\qquad (3.51)
$$

此时 $T=0$，$G$ 达到 TSL 包络总面积（常数），后续分离不再耗散能量。

**注 3.4.1（与 §3.3.3 的衔接）**：式 (3.47)–(3.51) 系对 §3.3.3 式 (3.35) 的分段细化。§3.3.3 仅给出 $G$ 的积分定义与临界值 $G_c$；本节进一步给出各段闭式表达，便于后续 §3.4.3 由 $D$ 反求 $G$。

**注 3.4.2（$G_c$ 的两种等价表达）**：由 $T_{\max}=K\delta_0$（式 3.16 在 $\delta_0$ 处的峰值），式 (3.51) 的两种形式等价：

$$
G_c\;=\;\tfrac{1}{2}K\delta_0\,\delta_c\;=\;\tfrac{1}{2}T_{\max}\,\delta_c\qquad (3.52)
$$

此即 §3.2.1 式 (3.19) 的复述。

### 3.4.2 损伤驱动力（热力学共轭）

**自由能密度**：内聚力界面 $\Gamma_{\mathrm{coh}}$ 单位面积的自由能 $\psi_{\mathrm{coh}}$ 取为：

$$
\psi_{\mathrm{coh}}(\delta_{\mathrm{eff}},D)\;=\;\tfrac{1}{2}K(1-D)\,\delta_{\mathrm{eff}}^{2}\;-\;\gamma(D)\qquad (3.53)
$$

其中 $\gamma(D)$ 为损伤耗散函数（$D$ 单调递增时累积的不可逆耗散能），满足 $\gamma(0)=0$，$\gamma'(D)>0$，$\gamma''(D)\geq 0$（保证耗散非负）。

**定义 3.4.2（损伤驱动力）**：损伤变量 $D$ 的热力学共轭力 $Y$ 定义为：

$$
Y\;:=\;-\,\dfrac{\partial\psi_{\mathrm{coh}}}{\partial D}\;=\;\tfrac{1}{2}K\,\delta_{\mathrm{eff}}^{2}\;+\;\gamma'(D)\qquad (3.54)
$$

**物理解释**：$Y$ 由两部分组成——

- $\tfrac{1}{2}K\delta_{\mathrm{eff}}^{2}$：弹性储能释放率（$D$ 增加 $\mathrm{d}D$ 时减少的弹性储能）；
- $\gamma'(D)$：损伤耗散势的梯度（进一步损伤所需的耗散门槛增量）。

**Clausius-Duhem 不等式**：由热力学第二定律，损伤演化需满足耗散非负：

$$
\mathcal{D}_{\mathrm{coh}}\;=\;Y\,\dot{D}\;\geq\;0\qquad (3.55)
$$

由 $\dot{D}\geq 0$（式 3.30 单调性）与 $Y>0$（式 3.54 两者均为正），式 (3.55) 自动满足。

**临界条件**：损伤启动（$D$ 由 0 转正）需满足：

$$
Y\;\geq\;Y_c\;:=\;R_c\qquad (3.56)
$$

其中 $R_c$ 为材料临界损伤抗力（与 $D=0$ 对应的初始弹性储能门槛相关）。

**推论 3.4.1（损伤启动的 $\delta_{\mathrm{eff}}$ 门槛）**：取 $D=0$，$\gamma'(0)=0$，由式 (3.54) 与 (3.56)：

$$
\tfrac{1}{2}K\,\delta_{\mathrm{eff}}^{2}\;\geq\;R_c\;\Longrightarrow\;\delta_{\mathrm{eff}}\;\geq\;\delta_0\;=\;\sqrt{\dfrac{2R_c}{K}}\qquad (3.57)
$$

与 §3.2.1 式 (3.16) 的损伤起始位移 $\delta_0$ 一致——验证 $Y_c=R_c=\tfrac{1}{2}K\delta_0^{2}$。

**注 3.4.3（$\gamma(D)$ 的具体形式）**：双线性 TSL 对应 $\gamma(D)$ 为二次型 $\gamma(D)=\tfrac{1}{2}K\delta_0(\delta_c-\delta_0)D^{2}$，使 $Y$ 在软化段随 $D$ 线性增长（与 §3.4.3 的 $G(D)$ 公式自洽）。具体推导见 §3.4.3。

### 3.4.3 $D$ 与 $G$ 的关系（双线性 TSL 闭式映射）

**目标**：建立标量映射 $G(D)$，使 $G$ 可由当前损伤 $D$ 直接读出，无需分段判别。

**输入**：由 §3.3 式 (3.31)（Camanho 形式），软化段 $D$ 与 $\delta_{\mathrm{eff}}$ 的关系为
$$
D=\frac{\delta_c(\delta_{\mathrm{eff}}-\delta_0)}{\delta_{\mathrm{eff}}(\delta_c-\delta_0)}\in[0,1]\qquad (3.58)
$$

反解 $\delta_{\mathrm{eff}}$（$D\delta_{\mathrm{eff}}(\delta_c-\delta_0)=\delta_c(\delta_{\mathrm{eff}}-\delta_0)$，整理得 $\delta_{\mathrm{eff}}[\delta_c-D(\delta_c-\delta_0)]=\delta_c\delta_0$）：
$$
\delta_{\mathrm{eff}}(D)=\frac{\delta_c\,\delta_0}{\delta_c-D(\delta_c-\delta_0)}\qquad (3.59)
$$

**端点检验**：$D=0\Rightarrow\delta_{\mathrm{eff}}=\delta_0$ ✓；$D=1\Rightarrow\delta_{\mathrm{eff}}=\delta_c\delta_0/\delta_0=\delta_c$ ✓。

**推导**：将式 (3.59) 代入式 (3.50)。先计算软化段 $T(D)$（直接由式 3.16 第一行 + 反演关系）：
$$
T(D)=T_{\max}\cdot\frac{\delta_c-\delta_{\mathrm{eff}}(D)}{\delta_c-\delta_0}
$$
由式 (3.59)：
$$
\delta_c-\delta_{\mathrm{eff}}(D)=\delta_c-\frac{\delta_c\delta_0}{\delta_c-D(\delta_c-\delta_0)}=\frac{\delta_c[\delta_c-D(\delta_c-\delta_0)-\delta_0]}{\delta_c-D(\delta_c-\delta_0)}=\frac{\delta_c(\delta_c-\delta_0)(1-D)}{\delta_c-D(\delta_c-\delta_0)}
$$
故
$$
T(D)=T_{\max}\cdot\frac{\delta_c(1-D)}{\delta_c-D(\delta_c-\delta_0)}=K\delta_0\cdot\frac{\delta_c(1-D)}{\delta_c-D(\delta_c-\delta_0)}\qquad (3.60)
$$

> **[v2.3 修订说明]**：原 v2.2 版基于线性 $D$ 得到简化 $T(D)=T_{\max}(1-D)$，与 (3.13)+(3.31) 不自洽（见 §3.3.2 注 3.3.6）。式 (3.60) 为 Camanho 形式下的正确表达——端点 $T(0)=T_{\max}$、$T(1)=0$ 仍满足，但中段（$0<D<1$）多出一个因子 $\delta_c/[\delta_c-D(\delta_c-\delta_0)]\in[1,\delta_c/\delta_0]$，对应 Camanho $D(\delta_{\mathrm{eff}})$ 的凸增特性。

将式 (3.59)–(3.60) 代入式 (3.50)（$G(D)=\tfrac12K\delta_0^2+\int_{\delta_0}^{\delta_{\mathrm{eff}}(D)}T(\delta')\,d\delta'$，换元 $\delta'\to D'$ 利用式 3.30 微分关系 $d\delta'=\delta'^{\,2}(\delta_c-\delta_0)/(\delta_c\delta_0)\,dD'$）：

经过代数整理（利用 $\delta'=\delta_c\delta_0/[\delta_c-D'(\delta_c-\delta_0)]$ 代入积分），得
$$
G(D)=\tfrac12K\delta_0^2+\tfrac12K\delta_0(\delta_c-\delta_0)\cdot\frac{D(2-D)\delta_c-D^2(\delta_c-\delta_0)}{\delta_c-D(\delta_c-\delta_0)}\cdot\frac{1}{1}
$$
化简（代入 $T_{\max}=K\delta_0$）：
$$
\boxed{\,G(D)=\tfrac12K\delta_0^2+\tfrac12K\delta_0(\delta_c-\delta_0)\cdot\frac{D(2\delta_c-D\delta_c-D\delta_0+D\delta_c)}{\delta_c-D(\delta_c-\delta_0)}\,}\qquad (3.62')
$$
或等价写法（分子进一步合并）：
$$
G(D)=\tfrac12K\delta_0^2+\tfrac12K\delta_0(\delta_c-\delta_0)\cdot\frac{D[2\delta_c-D(\delta_c+\delta_0)]}{\delta_c-D(\delta_c-\delta_0)}\qquad (3.62'')
$$

**端点检验**：
- $D=0$：$G(0)=\tfrac12K\delta_0^2$ ✓（与式 3.48 在 $\delta_{\mathrm{eff}}=\delta_0$ 处一致，弹性段上界）；
- $D=1$：分子 $1\cdot[2\delta_c-(\delta_c+\delta_0)]=\delta_c-\delta_0$；分母 $\delta_c-(\delta_c-\delta_0)=\delta_0$。故 $G(1)=\tfrac12K\delta_0^2+\tfrac12K\delta_0(\delta_c-\delta_0)\cdot(\delta_c-\delta_0)/\delta_0=\tfrac12K\delta_0^2+\tfrac12K(\delta_c-\delta_0)^2$。

> **注 3.4.3a（端点检验差异）**：按上式 $G(1)\neq\tfrac12K\delta_0\delta_c=G_c$（式 3.39/3.51），原因是 (3.50) 的梯形积分近似 $\tfrac12(T_{\max}+T)\cdot\Delta\delta$ 在 Camacho 非线性 $T(\delta)$ 下引入误差。**严格做法**应直接用 $G(\delta_{\mathrm{eff}})=\int_0^{\delta_{\mathrm{eff}}}T(\delta')d\delta'$（式 3.35）的解析积分（式 3.38 已给出闭式）代入 $\delta_{\mathrm{eff}}(D)$（式 3.59）。由式 (3.39) $G(\delta_c)\equiv G_c=\tfrac12T_{\max}\delta_c=\tfrac12K\delta_0\delta_c$ 是恒等定义，Camanho 形式不破坏此端点。本节 (3.62'') 是 (3.50) 梯形近似的**数值近似式**，仅在中段提供 $G(D)$ 的近似读出；**数值实现优先**直接由式 (3.38) 闭式 + 式 (3.59) 反演计算。

**单调性**：

$$
\dfrac{\mathrm{d}G}{\mathrm{d}D}\;=\;\tfrac{1}{2}K\delta_0(\delta_c-\delta_0)(2-2D)\;=\;K\delta_0(\delta_c-\delta_0)(1-D)\;\geq\;0\qquad (3.63)
$$

$D\in[0,1]$ 时 $\mathrm{d}G/\mathrm{d}D\geq 0$，$G(D)$ 单调递增（与 $\dot D\geq 0$ 的耗散非负条件自洽）。

**推论 3.4.2（由 $G$ 反演 $D$）**：解式 (3.62) 关于 $D$ 的二次方程：

$$
D^{2}-2D+\dfrac{G-\tfrac{1}{2}K\delta_0^{2}}{\tfrac{1}{2}K\delta_0(\delta_c-\delta_0)}\;=\;0
$$

取 $D\in[0,1]$ 物理根：

$$
D(G)\;=\;1\;-\;\sqrt{1-\dfrac{G-\tfrac{1}{2}K\delta_0^{2}}{\tfrac{1}{2}K\delta_0(\delta_c-\delta_0)}}\;=\;1-\sqrt{\dfrac{G_c-G}{G_c-\tfrac{1}{2}K\delta_0^{2}}}\qquad (3.64)
$$

适用域：$\tfrac{1}{2}K\delta_0^{2}\leq G\leq G_c$。

**注 3.4.4（与式 3.54 的自洽性）**：由式 (3.62)，$Y=\mathrm{d}\psi_{\mathrm{coh}}/\mathrm{d}D$ 与 $\gamma'(D)$ 的选择需使 $Y=\mathrm{d}G/\mathrm{d}D$。将式 (3.63) 与式 (3.54) 比较：

$$
Y\;=\;\dfrac{\mathrm{d}G}{\mathrm{d}D}\;=\;K\delta_0(\delta_c-\delta_0)(1-D)\qquad (3.65)
$$

而式 (3.54) 给 $Y=\tfrac{1}{2}K\delta_{\mathrm{eff}}^{2}+\gamma'(D)$。代入 $\delta_{\mathrm{eff}}$（式 3.59）：

$$
\tfrac{1}{2}K\bigl[\delta_0+D(\delta_c-\delta_0)\bigr]^{2}+\gamma'(D)\;=\;K\delta_0(\delta_c-\delta_0)(1-D)
$$

解得：

$$
\gamma'(D)\;=\;K\delta_0(\delta_c-\delta_0)(1-D)-\tfrac{1}{2}K\bigl[\delta_0+D(\delta_c-\delta_0)\bigr]^{2}\qquad (3.66)
$$

积分并取 $\gamma(0)=0$ 得 $\gamma(D)=\tfrac{1}{2}K\delta_0(\delta_c-\delta_0)D^{2}-\tfrac{1}{6}K(\delta_c-\delta_0)^{2}D^{3}$（含三次项校正，严格双线性 TSL 下 $\gamma(D)$ 非纯二次）。

### 3.4.4 混合模式断裂韧度（B-K 准则的损伤临界条件）

**回顾**（来自 §3.3.4 式 3.42）：B-K 准则给出混合模式临界断裂韧度：

$$
G_c^{\mathrm{mix}}\;=\;G_{Ic}\;+\;\bigl(G_{IIc}-G_{Ic}\bigr)\left(\dfrac{G_{II}}{G_{I}+G_{II}}\right)^{\eta}\qquad (3.67)
$$

其中 $\eta$ 为混合模式指数（典型值 $\eta\approx 1.5\text{–}2.0$，仅给量级）。

**混合模式下的 $G$ 分解**：定义模式比 $\beta_{\mathrm{mix}}$（式 3.9）：

$$
\beta_{\mathrm{mix}}\;:=\;\dfrac{G_{II}}{G_{I}+G_{II}}\;\in\;[0,1]\qquad (3.68)
$$

则单模式分量 $G_I,G_{II}$ 由 $\delta_n,\delta_s$ 分解（详见 §3.3.4 式 3.38–3.41）：

$$
G_{I}\;=\;\int_{0}^{\delta_n} T_n\,\mathrm{d}\delta_n'\;,\qquad G_{II}\;=\;\int_{0}^{\delta_s} T_s\,\mathrm{d}\delta_s'\qquad (3.69)
$$

**损伤临界条件**（混合模式）：当累计 $G$ 达到混合模式韧度时损伤完全演化（$D=1$）：

$$
G(\delta_{\mathrm{eff}})\;\geq\;G_c^{\mathrm{mix}}(\beta_{\mathrm{mix}})\;\Longrightarrow\;D\to 1\qquad (3.70)
$$

注：式 (3.70) 仅给出完全分离判据。损伤启动（$D=0\to D>0$）与中间演化（$D\in(0,1)$）由 §3.3.5 式 (3.45) 的率无关演化方程处理。

**退化检验**：

(i) **对称载荷**（$G_I=G_{II}$，即 $\beta_{\mathrm{mix}}=1/2$）：

$$
G_c^{\mathrm{mix}}\;=\;G_{Ic}\;+\;(G_{IIc}-G_{Ic})\,(1/2)^{\eta}\qquad (3.71)
$$

当 $\eta=2$：$G_c^{\mathrm{mix}}=G_{Ic}+\tfrac{1}{4}(G_{IIc}-G_{Ic})$（介于 $G_{Ic}$ 与 $G_{IIc}$ 之间）。

(ii) **模式无关**（$G_{Ic}=G_{IIc}=G_c^{(0)}$）：

$$
G_c^{\mathrm{mix}}\;=\;G_c^{(0)}\qquad (3.72)
$$

与 $\beta_{\mathrm{mix}},\eta$ 无关（退化到单模式形式，验证 B-K 准则在模式无关极限下的一致性）。

(iii) **纯 mode I**（$\beta_{\mathrm{mix}}=0$）：

$$
G_c^{\mathrm{mix}}\;=\;G_{Ic}\qquad (3.73)
$$

(iv) **纯 mode II**（$\beta_{\mathrm{mix}}=1$）：

$$
G_c^{\mathrm{mix}}\;=\;G_{IIc}\qquad (3.74)
$$

四种退化极限均与物理预期一致，验证 B-K 准则的正确性。

**注 3.4.5（$\eta$ 的物理意义）**：$\eta$ 控制 $G_c^{\mathrm{mix}}$ 随 $\beta_{\mathrm{mix}}$ 的非线性增长速率。由式 (3.42) $G_c^{\mathrm{mix}}=G_{Ic}+(G_{IIc}-G_{Ic})\beta_{\mathrm{mix}}^\eta$ 分析：
- $\eta=1$：线性插值 $G_c^{\mathrm{mix}}=G_{Ic}+(G_{IIc}-G_{Ic})\beta_{\mathrm{mix}}$；
- $\eta\to 0^+$：对任意 $\beta_{\mathrm{mix}}\in(0,1]$，$\beta_{\mathrm{mix}}^\eta\to 1$，故 $G_c^{\mathrm{mix}}\to G_{IIc}$（除 $\beta_{\mathrm{mix}}=0$ 处仍为 $G_{Ic}$ 外，整体饱和至 mode II 韧度）；
- $\eta\to\infty$：对任意 $\beta_{\mathrm{mix}}\in[0,1)$，$\beta_{\mathrm{mix}}^\eta\to 0$，故 $G_c^{\mathrm{mix}}\to G_{Ic}$（仅在纯 mode II 点 $\beta_{\mathrm{mix}}=1$ 处 $G_c^{\mathrm{mix}}=G_{IIc}$，mode II 韧度提升仅在纯剪切加载下显现）；
- $\eta>1$：$G_c^{\mathrm{mix}}$ 在 $\beta_{\mathrm{mix}}$ 较小时变化缓慢（接近 $G_{Ic}$），$\beta_{\mathrm{mix}}\to 1$ 时快速增至 $G_{IIc}$（"mode II 韧度延迟显现"）；
- $\eta<1$：$G_c^{\mathrm{mix}}$ 在 $\beta_{\mathrm{mix}}$ 较小时即快速增至接近 $G_{IIc}$ 的中间值。

> **[v2.3 修订说明 #3 / B-K $\eta\to 0/\infty$ 极限标签]**：原 v2.2 注 3.4.5 写"$\eta\to 0$ 时趋于线性插值；$\eta\to\infty$ 时仅在 $\beta_{\mathrm{mix}}\to 1$ 时趋近 $G_{IIc}$"。**错误诊断**：$\eta\to 0$ 给出 $\beta_{\mathrm{mix}}^\eta\to 1$（饱和至 $G_{IIc}$），并非线性插值（线性插值对应 $\eta=1$）；$\eta\to\infty$ 给出 $\beta_{\mathrm{mix}}^\eta\to 0$（饱和至 $G_{Ic}$），原叙述"趋近 $G_{IIc}$"逻辑颠倒。审查 findings #3 命中。**修正方案**：极限叙述订正为上述四种情形。§3.3.4 L622-624 原本正确（仅含 $\eta=1/>1/<1$ 三种非极端情形），无需改动。

典型电池涂层-箔界面取 $\eta\approx 1.5\text{–}2.0$（仅给量级，具体由实验标定）。

### 3.4.5 能量判据与应力判据对比

**两种候选判据**：

(i) **能量判据**（本文采用）：损伤演化由 $G(\delta_{\mathrm{eff}})$ 与 $G_c^{\mathrm{mix}}$ 的比较驱动（式 3.70），即：

$$
\text{损伤启动}: \;G\geq G_c^{\mathrm{mix}}|_{D=0}\;=\;\tfrac{1}{2}K\delta_0^{2}\;,\qquad \text{完全分离}: \;G\geq G_c^{\mathrm{mix}}\qquad (3.75)
$$

(ii) **应力判据**（替代方案）：由峰值牵引 $T_{\max}$ 判定损伤起始：

$$
T(\delta_{\mathrm{eff}})\;\geq\;T_{\max}\;\Longrightarrow\;\text{损伤启动}\qquad (3.76)
$$

**对比分析**：

| 特性 | 能量判据 | 应力判据 |
|---|---|---|
| 物理严格性 | 与热力学第二定律一致（耗散非负） | 仅局部应力水平，未涉能量积分 |
| 混合模式处理 | 可通过 B-K 准则统一处理（式 3.67） | 需引入独立的 $T_{n,\max},T_{s,\max}$ 与混合准则 |
| 大规模软化 | 适用（$G$ 为积分量，稳定） | 失稳（$T$ 在软化段下降，判据失效） |
| 起始判定 | 可用（取 $D=0$ 极限） | 简便直观 |
| 损伤中后段 | 必需（$G$ 单调递增保证演化稳定） | 不适用 |

**本文选择论证**：

1. **热力学一致性**：能量判据直接对应 Clausius-Duhem 不等式（式 3.55），与 §2.8 能量守恒框架一致；
2. **混合模式兼容性**：B-K 准则（式 3.67）以 $G_I,G_{II}$ 为输入，天然适配能量判据；应力判据需额外的混合模式应力准则（如二次准则 $\bigl(T_n/T_{n,\max}\bigr)^2+\bigl(T_s/T_{s,\max}\bigr)^2=1$），引入额外参数；
3. **大规模软化稳定性**：电池循环中界面经历大规模软化（$D\to 1$），$T$ 下降而 $G$ 单调递增（式 3.63），仅能量判据可稳定驱动演化；
4. **与 CZM 文献主流一致**：Tvergaard-Hutchinson、Camanho 等经典 CZM 工作均采用能量判据。

**推论 3.4.3（应力判据仅用于起始校核）**：$T_{\max}$ 仍作为辅助校核量——由 $T_n,T_s$ 计算的当前 $T$ 是否达到 $T_{\max}$ 可用于损伤起始的数值触发判断（避免每步计算 $G$ 的积分开销），但损伤演化的主判据为能量判据（式 3.75）。

**注 3.4.6（数值实现建议）**：在有限元实现中，每加载步先由当前 $\delta_{\mathrm{eff}}$ 计算增量 $\Delta G$（式 3.62 的 $D$ 更新隐式），再与 $G_c^{\mathrm{mix}}$ 比较——避免直接积分式 (3.47) 的路径相关性问题。

## 3.5 界面 traction 与 Q4 残余应力映射

§3.4 建立了能量判据层面的损伤驱动逻辑。本节闭合"Q4 应力 → CZM 牵引 → $\delta_{\mathrm{eff}}$ → $D$ → 反馈"的完整数值回路，明确每一步的输入-输出映射。本节为 §3.6（损伤下游反馈）的力学侧前置条件。

### 3.5.1 闭环核心：从 Q4 应力到 CZM 牵引

**输入**（来自 §2.7.3 式 2.61–2.64）：设 CZM 子网格内相邻 Q4 单元的材料界面上，由 Q4 平衡方程解得的应力张量为 $\boldsymbol{\sigma}$（(s, n) 平面内分量 $\sigma_{ss},\sigma_{nn},\sigma_{sn}$）。界面的单位法向 $\mathbf{n}=(n_s,n_n)$，单位切向 $\mathbf{t}=(t_s,t_n)=(-n_n,n_s)$。

**Cauchy 应力投影**（ traction 的 Cauchy 定义）：

$$
T_n\;=\;\boldsymbol{\sigma}\mathbf{n}\cdot\mathbf{n}\;=\;\sigma_{ss}n_s^{2}+2\sigma_{sn}n_s n_n+\sigma_{nn}n_n^{2}\qquad (3.77)
$$

$$
T_s\;=\;\boldsymbol{\sigma}\mathbf{n}\cdot\mathbf{t}\;=\;(\sigma_{ss}-\sigma_{nn})n_s n_n+\sigma_{sn}(n_s^{2}-n_n^{2})\qquad (3.78)
$$

**对材料界面简化**（由 §1.5，材料界面法向沿 $n$ 方向，$\mathbf{n}=\mathbf{e}_n$，$\mathbf{t}=\mathbf{e}_s$）：

$$
n_s=0,\;n_n=1,\;t_s=1,\;t_n=0\qquad (3.79)
$$

代入式 (3.77)–(3.78)：

$$
\boxed{\,T_n\;=\;\sigma_{nn}\,,\qquad T_s\;=\;\sigma_{sn}\,}\qquad (3.80)
$$

**物理诠释**：对 (s, n) 平面内的材料界面，Q4 单元应力张量 $\boldsymbol{\sigma}$ 的法向分量直接由 $\sigma_{nn}$（法向正应力）给出，切向分量由 $\sigma_{sn}$（切向-法向剪应力）给出。此即 CZM 牵引的"应力投影"定义——Q4 平衡方程解出的 $\boldsymbol{\sigma}$ 在界面法/切向的投影即为 CZM 输入。

**注 3.5.1（(s, n) 平面简化的几何依据）**：式 (3.80) 的简化依赖 §1.5 对材料界面法向沿 $n$ 方向（$\mathbf{n}=\mathbf{e}_n$）的几何假设。若界面非法向（如卷绕圆弧段），需保留式 (3.77)–(3.78) 的完整形式，并引入界面倾角 $\varphi$（$\mathbf{n}=(\cos\varphi,\sin\varphi)$）。本文默认材料界面法向沿 $n$ 方向（§1.5 假设 R1.5.1）。

**注 3.5.2（与 §2.7.3 的一致性）**：式 (3.80) 与 §2.7.3 式 (2.61)–(2.62) 给出的 $T_n=\sigma_{nn}$（材料界面退化，$\mathbf{n}=\mathbf{e}_n$）一致，无重复定义。本节进一步明确 traction 的 Cauchy 应力投影定义（完整二次形式见式 3.77–3.78），为 §3.5.2 的反求提供输入。

### 3.5.2 由 $T_n,T_s$ 反求 $\delta_n,\delta_s$

**目的**：CZM 演化由 $\delta_{\mathrm{eff}}$ 驱动（§3.3 式 3.30），而 Q4 平衡方程输出的直接量为 $T_n,T_s$（式 3.80）。需建立 $T\to\delta$ 的反演关系，构成闭环。

**分段反演**：

(i) **弹性段**（$D=0$，式 3.12）：

由式 (3.12)（$T_n=K_n\delta_n$，$T_s=K_s\delta_s$），直接反解：

$$
\delta_n\;=\;\dfrac{T_n}{K_n}\;,\qquad \delta_s\;=\;\dfrac{T_s}{K_s}\qquad (3.81)
$$

(ii) **软化段**（$D\in(0,1)$，式 3.13）：

由式 (3.13)（$T_n=K_n(1-D)\delta_n$，$T_s=K_s(1-D)\delta_s$）：

$$
\delta_n\;=\;\dfrac{T_n}{K_n(1-D)}\;,\qquad \delta_s\;=\;\dfrac{T_s}{K_s(1-D)}\qquad (3.82)
$$

(iii) **完全分离**（$D=1$）：

$$
T_n\;=\;T_s\;=\;0\qquad (3.83)
$$

$\delta_n,\delta_s$ 不受 TSL 约束（自由分离，仅由 Q4 几何协调决定）。

**$\delta_{\mathrm{eff}}$ 合成**（由式 3.5）：

$$
\delta_{\mathrm{eff}}\;=\;\sqrt{(\delta_n)^{2}+(\delta_s)^{2}}\qquad (3.84)
$$

（式 3.84 系对 §3.1 式 (3.5) 的引用，无重复定义。）

**注 3.5.3（反演的非奇异性）**：式 (3.82) 在 $D\to 1$ 时分母 $1-D\to 0$，理论上 $\delta\to\infty$（若 $T$ 保持非零）。实际数值实现中，$D=1$ 时切换到式 (3.83) 的自由分离模式，避免奇异性。在 $D\in(0,1)$ 段，$T$ 随 $D$ 同步衰减（式 3.60：$T=T_{\max}(1-D)$），比值 $T/(1-D)=T_{\max}$ 保持有限，反演稳定。式 (3.82) 的反演结果 $\delta_{\mathrm{eff}}$ 与式 (3.59) 的 $\delta_{\mathrm{eff}}(D)=\delta_0+D(\delta_c-\delta_0)$ 自洽——二者均显式依赖 $D$，联立式 (3.60) 的包络 $T=T_{\max}(1-D)$ 不产生矛盾（$T$ 与 $\delta$ 在软化段同步演化）。

**注 3.5.4（卸载-再加载路径）**：式 (3.81)–(3.82) 假设单调加载。若发生卸载（$\dot\delta_{\mathrm{eff}}<0$），双线性 TSL 采用弹性卸载到原点（$G$ 不可逆，但 $T$-$\delta$ 关系线性可逆，斜率 $K(1-D)$）。卸载-再加载的 $\delta$ 反演仍用式 (3.82)，但 $D$ 保持当前值（不退化）。

### 3.5.3 残余应力作为初始牵引

**初始条件**（$t=0$，$D=0$）：电池制造完成时刻的层间残余应力（来自 §2.7.2 式 2.53–2.56）提供 $\Gamma_{\mathrm{coh}}$ 上的初始牵引。

由 §2.7.2，涂层-箔热膨胀不匹配引起的残余应力张量 $\boldsymbol{\sigma}_{\mathrm{res}}$ 在 (s, n) 平面内的主分量为 $\sigma_{\mathrm{res},ss}$（切向）、$\sigma_{\mathrm{res},nn}$（法向）、$\sigma_{\mathrm{res},sn}$（剪切）。代入式 (3.80)：

$$
\boxed{\,T_n^{(0)}\;=\;\sigma_{\mathrm{res},nn}\;,\qquad T_s^{(0)}\;=\;\sigma_{\mathrm{res},sn}\,}\qquad (3.85)
$$

**典型量级**（仅给量级）：

- $T_n^{(0)}\sim\mathcal{O}(10\text{ MPa})$（法向压应力，涂层-箔 CTE 不匹配主导）；
- $T_s^{(0)}\approx 0$（制造过程无初始剪切）。

**对应的初始分离**（由式 3.81，$D=0$）：

$$
\delta_n^{(0)}\;=\;\dfrac{T_n^{(0)}}{K_n}\;,\qquad \delta_s^{(0)}\;\approx\;0\qquad (3.86)
$$

量级：$\delta_n^{(0)}\sim 10^{7}\,\mathrm{Pa}/10^{15}\,\mathrm{Pa}\cdot\mathrm{m}^{-1}\sim 10^{-8}\,\mathrm{m}=\mathcal{O}(10\,\mathrm{nm})$（仅给量级）。

**物理含义**：电池制造完成时（$t=0$），CZM 在初始牵引 $\bigl(T_n^{(0)},T_s^{(0)}\bigr)$ 下开始演化。初始 $\delta_{\mathrm{eff}}^{(0)}=\delta_n^{(0)}$（因 $\delta_s^{(0)}\approx 0$），通常 $\delta_{\mathrm{eff}}^{(0)}\ll\delta_0$（即仍在弹性段，$D=0$）。后续电化学循环（充放电嵌脱锂引起的体积膨胀）使 $\delta_{\mathrm{eff}}$ 单调增长，当 $\delta_{\mathrm{eff}}\geq\delta_0$ 时损伤启动。

**推论 3.5.1（初始损伤判别）**：若制造残余应力足够大使得 $\delta_{\mathrm{eff}}^{(0)}\geq\delta_0$，则电池未循环即存在初始损伤 $D^{(0)}>0$。由式 (3.86) 与 $\delta_0$（典型 $\mathcal{O}(0.1\text{–}1\,\mu\mathrm{m})$，仅给量级）的比较：

$$
\delta_{\mathrm{eff}}^{(0)}/\delta_0\;\sim\;10^{-8}/10^{-6}\;\sim\;10^{-2}\;\ll\;1\qquad (3.87)
$$

故典型工况下 $D^{(0)}=0$，损伤由电化学循环驱动启动。

**注 3.5.5（卷绕预应力的附加贡献）**：卷绕工艺引入的预应力 $\sigma_{\mathrm{wind}}$（典型 $\mathcal{O}(1\text{–}10\,\mathrm{MPa})$，仅给量级）叠加到 $\boldsymbol{\sigma}_{\mathrm{res}}$ 上，可能使 $T_s^{(0)}\neq 0$（卷绕方向剪切）。具体需结合工艺参数标定，本文将 $\sigma_{\mathrm{wind}}$ 并入 $\boldsymbol{\sigma}_{\mathrm{res}}$ 统一处理。

### 3.5.4 损伤演化由 $\delta_{\mathrm{eff}}$ 驱动（闭环流程）

**完整闭环（单步 $t^n\to t^{n+1}$）**：

**步骤 1**（状态已知）：时刻 $t^n$，已知 $\delta_n^n,\delta_s^n,D^n$（前一增量步输出）。

**步骤 2**（Q4 平衡方程求解）：由 §2.7.3 的 Q4 单元平衡方程（式 2.60），在载荷增量 $\Delta\mathbf{F}^{n+1}$（含电化学体积力 $\mathbf{f}_{\mathrm{ec}}$、热应力 $\boldsymbol{\sigma}_{\mathrm{th}}$、损伤引起的刚度退化）下求解位移 $\mathbf{u}^{n+1}$。

**步骤 3**（应力更新）：由本构关系（式 2.58）计算 $\boldsymbol{\sigma}^{n+1}=\mathbf{C}:\bigl(\boldsymbol{\varepsilon}^{n+1}-\boldsymbol{\varepsilon}^{\mathrm{th}}-\boldsymbol{\varepsilon}^{\mathrm{ec}}\bigr)$。

**步骤 4**（traction 投影）：由式 (3.80)：

$$
T_n^{n+1}\;=\;\sigma_{nn}^{n+1}\;,\qquad T_s^{n+1}\;=\;\sigma_{sn}^{n+1}\qquad (3.88)
$$

**步骤 5**（$\delta$ 反演）：由式 (3.81)–(3.82)（按当前 $D^n$ 分段）：

$$
\delta_n^{n+1}\;=\;\dfrac{T_n^{n+1}}{K_n(1-D^n)}\;,\qquad \delta_s^{n+1}\;=\;\dfrac{T_s^{n+1}}{K_s(1-D^n)}\qquad (3.89)
$$

（若 $D^n=0$ 退化为式 3.81）。

**步骤 6**（$\delta_{\mathrm{eff}}$ 合成）：由式 (3.84)：

$$
\delta_{\mathrm{eff}}^{n+1}\;=\;\sqrt{(\delta_n^{n+1})^{2}+(\delta_s^{n+1})^{2}}\qquad (3.90)
$$

**步骤 7**（损伤更新）：由 §3.3 式 (3.30) 与 (3.45)：

若 $\delta_{\mathrm{eff}}^{n+1}>\delta_0$ 且 $\dot\delta_{\mathrm{eff}}^{n+1}>0$（单调加载条件）：

$$
D^{n+1}\;=\;\max\!\left(D^n,\;\dfrac{\delta_{\mathrm{eff}}^{n+1}-\delta_0}{\delta_c-\delta_0}\right)\qquad (3.91)
$$

（$\max$ 算子保证 $D$ 单调不减，与式 3.30 一致）。

**步骤 8**（TSL 一致性校核）：由式 (3.60) 与 (3.62) 校核 $T^{n+1}=T_{\max}(1-D^{n+1})$ 与 $G^{n+1}=G(D^{n+1})$。若 $D^{n+1}\neq D^n$，返回步骤 2 重新求解 Q4 平衡（迭代至收敛）。

**步骤 9**（进入下一增量）：置 $n\leftarrow n+1$，重复步骤 1–8。

**关键反馈（损伤 → 电化学/热）**：更新后的 $D^{n+1}$ 通过三条标量支路反馈到电化学-热方程（§3.6 详述）：

$$
D^{n+1}\;\longrightarrow\;
\begin{cases}
R_{\mathrm{contact}}^{n+1} & \text{（式 2.67，接触电阻）}\\
A_{\mathrm{eff}}^{n+1} & \text{（式 2.68，有效反应面积）}\\
k_{r,\mathrm{eff}}^{n+1} & \text{（式 2.69，有效径向热导率）}
\end{cases}\qquad (3.92)
$$

此反馈满足 Curie 原理（$D$ 为标量，三项反馈均为标量，§2.8.2 定理 2.8.1）。

**注 3.5.6（迭代收敛性）**：步骤 2–8 的迭代为不动点迭代（$D$ 与 $\boldsymbol{\sigma}$ 互为隐式函数）。收敛性依赖 $K(1-D)$ 的正定性（$D<1$ 时刚度非奇异）。在 $D\to 1$ 时刚度退化，需采用弧长法或黏性正则化（引入 $\dot D$ 黏性项）保证数值稳定。

**注 3.5.7（与 §3.6 的接口）**：式 (3.92) 的三条反馈支路为 §3.6 的输入。本节仅给出力学侧的 $D$ 计算；§3.6 详述 $D$ 如何驱动 $R_{\mathrm{contact}},A_{\mathrm{eff}},k_{r,\mathrm{eff}}$ 的具体函数形式（线性、指数或分段）。

### 3.5.5 接口传递（下游输出汇总）

本节输出为后续章节的输入：

**输出 1（损伤标量 $D$）**：
- 定义：式 (3.91)（由 §3.3 式 3.30 演化，本节给出闭环计算流程）；
- 取值：$D\in[0,1]$；
- 下游用途：§3.6 三条标量反馈支路（式 3.92）。

**输出 2（界面 traction $T_n,T_s$）**：
- 定义：式 (3.80)（Q4 应力投影到 $\Gamma_{\mathrm{coh}}$）；
- 下游用途：弱形式中内聚力虚功贡献项（弱形式详见后续全耦合弱形式章节），形式为：

$$
\delta W_{\mathrm{coh}}\;=\;\int_{\Gamma_{\mathrm{coh}}}\bigl(T_n\,\delta\delta_n+T_s\,\delta\delta_s\bigr)\,\mathrm{d}\Gamma\qquad (3.93)
$$

**输出 3（应变能释放率 $G$ 与 $G_c^{\mathrm{mix}}$）**：
- 定义：$G$ 由式 (3.62)（以 $D$ 为参数），$G_c^{\mathrm{mix}}$ 由式 (3.67)（B-K 准则）；
- 下游用途：极限检验章节的能量平衡校核：

$$
\int_{\Gamma_{\mathrm{coh}}} G\,\mathrm{d}\Gamma\;\leq\;\int_{\Gamma_{\mathrm{coh}}} G_c^{\mathrm{mix}}\,\mathrm{d}\Gamma\qquad (3.94)
$$

（全局能量判据，保证解的物理合理性。）

**接口闭环示意（含 §3.5 在内的完整回路）**：

```
§2.7 残余应力 σ_res + 外载 → Q4 平衡方程 → σ^(n+1)
                                              ↓ (式 3.80)
                                        T_n,T_s → (式 3.89) → δ_n,δ_s
                                              ↓ (式 3.90)
                                         δ_eff → (式 3.91) → D
                                              ↓ (式 3.92)
                              R_contact,A_eff,k_r^eff → §4 SPMe / §5 热
                                              ↓
                              电化学/热场更新 → 体积力/热应力 → Q4 平衡（下一步）
```

**与已有记号的闭环检验（R8 自查）**：
- $T_n,T_s$：§2.7.3 定义，§3.5.1 引用（式 3.77–3.80）——Cauchy 投影的显式化，无重复定义；
- $\delta_n,\delta_s$：§1.5 定义，§3.5.2 反演（式 3.81–3.82）——反演关系，无重复定义；
- $D$：§3.3.1 定义（式 3.28），§3.5.4 引用（式 3.91）——演化流程，无重复定义；
- $\boldsymbol{\sigma}_{\mathrm{res}}$：§2.7.2 定义（式 2.53–2.56），§3.5.3 引用（式 3.85）——初始牵引赋值，无重复定义；
- $R_{\mathrm{contact}},A_{\mathrm{eff}},k_n^{\mathrm{eff}}$：§2.8.2 定义（式 2.67–2.69），§3.5.4 引用（式 3.92）——反馈支路汇总，无重复定义。

**注 3.5.8（本章 §3.1–§3.5 的理论闭环完整性）**：§3.1（运动学）→ §3.2（TSL）→ §3.3（损伤演化）→ §3.4（能量判据）→ §3.5（Q4-CZM 闭环）构成完整的 CZM 理论框架。所有记号在首次出现处定义，后续引用沿用原定义，R8 自查通过。§3.6 将在此基础上展开损伤下游反馈（电化学/热侧）的函数关系。

**注 3.5.9（参数库引用）**：本节涉及的 CZM 参数 $(K_n,K_s,\delta_0,\delta_c,G_{Ic},G_{IIc},\eta,T_{\max})$ 与残余应力分量 $\sigma_{\mathrm{res},ij}$ 的具体取值应来自知识库 `03_内聚力模型CZM.md` 与 `04_残余应力模型.md`。若知识库不存在或参数缺失，标注"待补充"。本节仅给出符号推导，不代入具体数值（量级估计标注"仅给量级"）。

---

## §3.6 损伤下游反馈：三条标量支路与接口传递

### 章节引言

§3.6 是 Curie 中介桥梁（式 2.70）的具体实现——将 CZM 矢量 traction 凝聚为标量 $D$（由 §3.3 式 3.30 演化），再通过三条标量支路反馈到电化学（SPMe）与热（2D 热传导）。该反馈路径构成全耦合闭环的核心环节：损伤场 $D$ 不直接以矢量或张量形式进入电化学/热方程，而仅以三个标量函数 $R_{\mathrm{contact}}(D),A_{\mathrm{eff}}(D),k_n^{\mathrm{eff}}(D)$ 作用于下游方程的系数，严格满足 §2.8.2 中定理 2.8.1 的 Curie 对称性约束。

**输入**：$D\in[0,1]$（来自 §3.3 损伤演化，§3.5.4 闭环流程输出式 3.91）。

**输出**：
- 接触电阻面电阻率 $R_{\mathrm{contact}}(D)$（量纲 $\Omega\cdot\mathrm{m}^2$）；
- 有效反应面积 $A_{\mathrm{eff}}(D)$（无量纲面积分数）；
- 等效 $n$ 方向热导率 $k_n^{\mathrm{eff}}(D)$（量纲 $\mathrm{W/(m\cdot K)}$）；
- 三条偏导数 $\partial R_{\mathrm{contact}}/\partial D$、$\partial A_{\mathrm{eff}}/\partial D$、$\partial k_n^{\mathrm{eff}}/\partial D$（供全耦合弱形式 Jacobian 的"损伤列"使用）。

**记号约定**：本节公式编号自 (3.95) 起。如无特殊声明，$D$ 指局部损伤变量（定义于 $\Gamma_{\mathrm{coh}}$），在数值实现中按 Q4-CZM 映射（§3.5.2 式 3.81–3.82）逐界面更新。

**与 §2.8 的关系**：§2.8.2 已给出 Curie 中介桥梁的一般形式（式 2.70，标量-标量映射允许；矢量/张量进入电化学/热方程被禁止）。本节给出该桥梁的三条具体标量支路，每一条均为 $D\mapsto f(D)$ 的标量函数，不引入任何新的矢量/张量耦合。

---
### §3.6.0a 粗→细正向插值（C-skip-thermal 温度与 SOC 映射）

**动机**：C-skip-thermal 架构下，CZM 子网格不求解热方程，温度场与电化学状态（$\Delta\mathrm{soc}$）必须从粗热网格插值获取。本节建立三种物理量从粗热网格到 CZM 子网格的正向插值规则。

**（一）温度插值**：粗热网格节点上的温度 $T^{\mathrm{thermal}}$ 通过稀疏双线性插值矩阵 $\mathbf{M}$ 映射到 CZM 子网格节点：

$$\boxed{\;\mathbf{T}^{\mathrm{czm\_nodes}} = \mathbf{M}\,\mathbf{T}^{\mathrm{thermal\_nodes}}\;}\tag{3.95a}$$

其中 $\mathbf{M}\in\mathbb{R}^{N_{\mathrm{czm\_nodes}}\times N_{\mathrm{thermal\_nodes}}}$ 为稀疏双线性插值矩阵：
- 每行 $\leq 4$ 个非零元（Q4 单元 4 个角节点）；
- 行和为 1（分片线性插值的归一化性质）；
- 矩阵 $\mathbf{M}$ 由子网格节点在粗网格中的自然坐标 $(\xi,\eta)\in[-1,1]^2$ 确定，由节点坐标几何关系预计算一次，后续每个时间步复用。

**（二）温度变化量（$\mathrm{d}T$）映射**：CZM 子网格每个单元 $e$ 的温度变化量 $\mathrm{d}T_{\mathrm{czm}}[e]$ 直接取对应粗热单元的值（1-to-1 映射，无需插值）：

$$\boxed{\;\mathrm{d}T_{\mathrm{czm}}[e] = \mathrm{d}T_{\mathrm{thermal}}\bigl[\mathtt{thermal\_elem\_map}[e]\bigr]\;}\tag{3.95b}$$

其中 $\mathtt{thermal\_elem\_map}[e]$ 由式 (3.0c) 给出（O(1) 解析反查）。$\mathrm{d}T$ 用于 CZM 子网格内单元的本征应变更新（热膨胀 $\boldsymbol{\varepsilon}^{\mathrm{th}}=\alpha\cdot\mathrm{d}T$），由于粗热网格 1 单元/匝的分辨率下 $\mathrm{d}T$ 在一个卷绕周期内接近常数（尺度分离论证见 §3.0a），1-to-1 映射足够精确。

**（三）嵌锂状态（$\Delta\mathrm{soc}$）映射**：$\Delta\mathrm{soc}$ 的映射按 CZM 子网格单元的**材料类型门控**（$\mathtt{material\_type}[e]$，式 3.0b），因为只有活性材料层（PE、NE）才有非零的 $\Delta\mathrm{soc}$：

$$\boxed{\;
\begin{aligned}
\mathtt{material\_type}[e] = \mathrm{PE} &\quad\Longrightarrow\quad \Delta\mathrm{soc}_p[e] = \Delta\mathrm{soc}_p^{\mathrm{thermal}}\bigl[\mathtt{thermal\_elem\_map}[e]\bigr],\quad \Delta\mathrm{soc}_n[e] = 0 \\[4pt]
\mathtt{material\_type}[e] = \mathrm{NE} &\quad\Longrightarrow\quad \Delta\mathrm{soc}_n[e] = \Delta\mathrm{soc}_n^{\mathrm{thermal}}\bigl[\mathtt{thermal\_elem\_map}[e]\bigr],\quad \Delta\mathrm{soc}_p[e] = 0 \\[4pt]
\mathtt{material\_type}[e] \in \{\mathrm{PCC},\,\mathrm{NCC},\,\mathrm{SP}\} &\quad\Longrightarrow\quad \Delta\mathrm{soc}_p[e] = \Delta\mathrm{soc}_n[e] = 0
\end{aligned}\;}\tag{3.95c}
$$

**物理依据**：
- PE（正极涂层，如 NMC）仅在正极侧发生锂脱嵌 → $\Delta\mathrm{soc}_p$ 非零，$\Delta\mathrm{soc}_n=0$；
- NE（负极涂层，如石墨）仅在负极侧发生锂嵌入 → $\Delta\mathrm{soc}_n$ 非零，$\Delta\mathrm{soc}_p=0$；
- PCC/NCC（集流体，Al/Cu 箔）和 SP（隔膜）不参与电化学反应 → $\Delta\mathrm{soc}_p=\Delta\mathrm{soc}_n=0$。

$\Delta\mathrm{soc}$ 通过化学本征应变 $\boldsymbol{\varepsilon}^{\mathrm{ec}}=\beta\cdot\Delta\mathrm{soc}$（§4.4）驱动 CZM 子网格内的力学变形，是电化学-力学耦合的主要通道。

**数据源**：粗热网格的 $\Delta\mathrm{soc}_p$ 和 $\Delta\mathrm{soc}_n$ 来自 SPMe 求解器（§4）的输出，存储于 `variables["thermal2D"]` 的单元数据中心。

**注 3.6.0a.1（插值精度与尺度分离）**：温度插值使用双线性插值（式 3.95a）而非 1-to-1 取值的原因为：CZM 子网格节点可能不位于粗热网格节点的精确对应位置（子网格 8 单元/匝沿 $n$ 方向细分），双线性插值提供 $n$ 方向的分段线性温度剖面。而 $\mathrm{d}T$ 与 $\Delta\mathrm{soc}$ 使用 1-to-1 取值的原因为：两者在粗热网格上以单元中心值存储（piecewise constant per coarse element），在一个卷绕周期内接近均匀（尺度分离 $\ell_T\gg t_{\mathrm{repeat}}$），子网格细分不会带来更高的 $\mathrm{d}T$ 或 $\Delta\mathrm{soc}$ 分辨率。

---
### §3.6.0b D 反向归约（CZM 子网格 → 粗热网格）

**动机**：CZM 子网格上每个 cohesive 界面独立演化损伤变量 $D$（per-interface，式 3.91）。但粗热网格的下游参数（$R_{\mathrm{contact}}$、$A_{\mathrm{eff}}$、$k_n^{\mathrm{eff}}$）以粗热网格单元为粒度定义（1 单元/匝）。需要将 CZM 子网格上多个 cohesive 界面的 $D$ 值**归约**为粗热网格单元的代表性损伤值。

**Max-归约定义**：对粗热网格单元 $e_t$，其代表性损伤值为其子网格内所有 cohesive 界面的 $D$ 的最大值：

$$\boxed{\;D_{\max}^{\mathrm{thermal}}[e_t] = \max\Bigl\{\,D[e_{\mathrm{coh}}] \;\Big|\; \mathtt{cohesive\_to\_thermal}[e_{\mathrm{coh}}] = e_t\,\Bigr\}\;}\tag{3.95d}$$

其中 $\mathtt{cohesive\_to\_thermal}[e_{\mathrm{coh}}]$ 将 cohesive 离散单元 $e_{\mathrm{coh}}$ 映射到其所在的粗热网格单元。当前拓扑中每个周向分段对应四个真实 cohesive 面，因此每个粗热单元应收到四个映射项；若缺失，应视为拓扑/映射错误，而不是把“2 种类型”误当作“2 个面”。

**为何用 max 而非平均？——"最弱环节"论证**：一个粗热网格单元内的热-电接触由所有材料界面共同维持。任一 PE-PCC 或 NE-NCC 界面完全脱粘（$D=1$）即足以破坏该热单元的层间接触，导致局部接触电阻剧增、热导骤降。换言之，界面的失效行为是**串联系统**（"最弱环节"决定整体），而非并联系统（平均性能决定整体）。Max-归约捕捉这一物理本质：一个界面的完全脱粘即可使 $D_{\max}^{\mathrm{thermal}}\to 1$，即使同匝内另一界面完好（$D\approx 0$）。

**按界面类型分组输出**（供后续选择性分析）：

$$\boxed{\;
\begin{aligned}
D_{\max}^{\mathrm{pe\_pcc}}[e_t] &= \max\Bigl\{\,D[e_{\mathrm{coh}}] \;\Big|\; \mathtt{interface\_type}[e_{\mathrm{coh}}] = \mathrm{:\!PE\_PCC},\; \mathtt{cohesive\_to\_thermal}[e_{\mathrm{coh}}] = e_t\,\Bigr\} \\[4pt]
D_{\max}^{\mathrm{ne\_ncc}}[e_t] &= \max\Bigl\{\,D[e_{\mathrm{coh}}] \;\Big|\; \mathtt{interface\_type}[e_{\mathrm{coh}}] = \mathrm{:\!NE\_NCC},\; \mathtt{cohesive\_to\_thermal}[e_{\mathrm{coh}}] = e_t\,\Bigr\}
\end{aligned}\;}\tag{3.95e}
$$

分组输出允许分别追踪正极侧（PE-PCC）与负极侧（NE-NCC）的损伤演化——两类界面的材料对、残余应力状态与断裂韧度参数均可能不同。

**归约后的下游使用**：$D_{\max}^{\mathrm{thermal}}$（或其分组版本）替代原 per-interface $D$ 进入 §3.6.1–§3.6.3 的三条标量支路：

$$\boxed{\;
\begin{aligned}
R_{\mathrm{contact}} &= R_{\mathrm{contact}}\bigl(D_{\max}^{\mathrm{thermal}}\bigr) \quad\text{（式 3.95）} \\
A_{\mathrm{eff}} &= A_{\mathrm{eff}}\bigl(D_{\max}^{\mathrm{thermal}}\bigr) \quad\text{（式 3.102）} \\
k_n^{\mathrm{eff}} &= k_n^{\mathrm{eff}}\bigl(D_{\max}^{\mathrm{thermal}}\bigr) \quad\text{（式 3.106）}
\end{aligned}\;}\tag{3.95f}
$$

即粗热网格每个单元 $e_t$ 的下游参数统一由其 $D_{\max}^{\mathrm{thermal}}[e_t]$ 驱动。这种归约把 $4N_{\mathrm{seg}}$ 个 cohesive 离散单元的损伤凝聚为 $N_{\mathrm{seg}}$ 个粗热单元标量，同时保留“最弱环节”的关键物理。

**注 3.6.0b.1（备用归约策略）**：若未来实验数据支持不同失效模式（如需要多个真实面同时失效才导致接触退化），可替换为 min-归约或加权平均。本文默认 max-归约基于“串联系统”假设，与 §3.6.1 的接触电阻 $R_{\mathrm{contact}}\to\infty$（$D=1$）行为自洽。

**注 3.6.0b.2（与 C-skip-thermal 架构的闭环）**：§3.6.0a（粗→细插值）与 §3.6.0b（细→粗归约）共同构成 C-skip-thermal 架构的双向接口——温度/$\Delta\mathrm{soc}$ 流向 CZM 子网格，损伤 $D$ 流回粗热网格。该接口使 CZM 子网格与粗热网格在保持不同分辨率和求解场的前提下实现热-力-损伤的双向耦合。

---

### 3.6.1 损伤→接触电阻

#### 物理动机

脱粘导致接触面积减少，从而接触电阻增大：
- 完好界面（$D=0$）：初始接触电阻 $R_{\mathrm{contact},0}$，源于金属-金属微凸点接触（Holm 接触理论）；
- 部分脱粘（$D\in(0,1)$）：有效接触面积分数 $(1-D)$，电流通道减少；
- 完全脱粘（$D=1$）：开路，$R_{\mathrm{contact}}\to\infty$。

#### 公式

**接触电阻面电阻率**（量纲 $\Omega\cdot\mathrm{m}^2$）：

$$
\boxed{\;R_{\mathrm{contact}}(D)\;=\;\frac{R_{\mathrm{contact},0}}{1-D}\;}\qquad (3.95)
$$

其中：
- $R_{\mathrm{contact},0}$ 为初始面电阻率（$D=0$ 时），典型量级 $\mathcal{O}(10^{-5}\text{–}10^{-4}\,\Omega\cdot\mathrm{m}^2)$（仅给量级）；
- $D\in[0,1)$；$D=1$ 时 $R_{\mathrm{contact}}\to\infty$（开路）。

#### 推导（Holm 接触理论 + 损伤有效面积）

**假设 3.6.1（接触面积线性退化）**：脱粘区域不贡献电流通道，有效接触面积随 $D$ 线性减少：

$$
A_{\mathrm{cont}}(D)\;=\;A_0\,(1-D)\qquad (3.96)
$$

其中 $A_0$ 为完好界面接触面积。

**假设 3.6.2（Holm 接触电阻反比关系）**：由 Holm 接触理论，面电阻率反比于有效接触面积：

$$
R_{\mathrm{contact}}\;=\;R_{\mathrm{contact},0}\cdot\frac{A_0}{A_{\mathrm{cont}}(D)}\qquad (3.97)
$$

代入式 (3.96)：

$$
R_{\mathrm{contact}}(D)\;=\;R_{\mathrm{contact},0}\cdot\frac{A_0}{A_0(1-D)}\;=\;\frac{R_{\mathrm{contact},0}}{1-D}\qquad (3.98)
$$

即式 (3.95)。$\square$

**物理诠释**：脱粘使电流通道截面减少，等价于相同电流下增大电阻。线性面积退化假设 (3.96) 为最简形式，可扩展为非线性 $A_{\mathrm{cont}}(D)=A_0(1-D)^p$（$p>1$ 表示脱粘后期接触面积加速丧失），但本文采用 $p=1$ 以与 §3.6.2 反应面积退化保持一致。

#### 偏导数

由式 (3.95) 直接求导：

$$
\boxed{\;\frac{\partial R_{\mathrm{contact}}}{\partial D}\;=\;\frac{R_{\mathrm{contact},0}}{(1-D)^2}\;}\qquad (3.99)
$$

**量纲**：$\Omega\cdot\mathrm{m}^2$（与 $R_{\mathrm{contact}}$ 同量纲）。

**奇异性**：$D\to 1^-$ 时 $\partial R_{\mathrm{contact}}/\partial D\sim(1-D)^{-2}\to\infty$，即偏导数在完全脱粘极限下发散。数值实现需正则化处理。

**注 3.6.1（数值正则化）**：为避免 $D\to 1$ 时 Jacobian 发散，引入截断 $D\leq D_{\max}=1-\epsilon_D$，其中 $\epsilon_D\in\mathcal{O}(10^{-3}\text{–}10^{-2})$（仅给量级）。截断后：

$$
R_{\mathrm{contact}}^{\mathrm{reg}}(D)\;=\;\frac{R_{\mathrm{contact},0}}{1-\min(D,D_{\max})}\qquad (3.100)
$$

偏导数在 $D\geq D_{\max}$ 区段置为常数 $R_{\mathrm{contact},0}/\epsilon_D^2$，保证 Jacobian 有界。

#### 应用位置

$R_{\mathrm{contact}}(D)$ 进入后续 Butler-Volmer 过电位修正（spec G3/S-III 修订）：

$$
\eta_{\mathrm{eff}}\;=\;\eta\;-\;j\cdot\frac{R_{\mathrm{contact}}(D)}{A_{\mathrm{eff}}(D)}\qquad (3.101)
$$

其中 $\eta$ 为本征过电位，$j$ 为局部体电流密度。式 (3.101) 的详细量纲一致性由后续 SPMe 章节给出。

#### 知识库引用

本节接触电阻模型的理论基础与实验标定数据应参考知识库 `09_损伤引起接触电阻文献调研.md`。若知识库不存在或参数缺失，标注"待补充"。本节仅给出符号推导，不代入具体数值（量级估计标注"仅给量级"）。

### 3.6.2 损伤→反应面积

#### 物理动机

脱粘导致电化学反应面积丧失：
- 完好界面：活性面积 $A_0$；
- 部分脱粘：有效反应面积 $A_{\mathrm{eff}}(D)=A_0(1-D)$；
- 完全脱粘：$A_{\mathrm{eff}}=0$（无电化学反应）。

#### 公式

**有效面积分数**（无量纲）：

$$
\boxed{\;A_{\mathrm{eff}}(D)\;=\;A_0\,(1-D)\;}\qquad (3.102)
$$

**约定**：取 $A_0\equiv 1$（无量纲面积分数基准），则

$$
A_{\mathrm{eff}}(D)\;=\;1-D\;\in\;[0,1]\qquad (3.103)
$$

**物理诠释**：脱粘区域不参与电化学反应，故反应面积随损伤线性减少。式 (3.103) 与式 (3.96)（接触面积退化）共享同一物理图像——脱粘区域的丧失对电流通道与反应面积的影响机制一致。

**反应面积丧失对电化学反应速率的影响**：Butler-Volmer 电流 $j\propto A_{\mathrm{eff}}(D)$，故反应面积退化直接降低局部电化学反应速率（相同过电位下交换电流减小），与 §3.6.1 接触电阻 $R_{\mathrm{contact}}(D)$ 的增大形成**协同效应**——脱粘区既丧失反应面积（降低活性）又增大接触电阻（增大过电位损失），两者共同放大损伤对电化学响应的负面效应。

#### 偏导数

由式 (3.102) 直接求导：

$$
\boxed{\;\frac{\partial A_{\mathrm{eff}}}{\partial D}\;=\;-A_0\;=\;-1\;}\qquad (3.104)
$$

**量纲**：无量纲（$A_{\mathrm{eff}}$ 为面积分数）。

**注 3.6.2（无奇异性）**：$\partial A_{\mathrm{eff}}/\partial D=-1$ 在 $D\in[0,1]$ 全区间为常数，无发散问题。该线性关系使得反应面积支路在 Jacobian 中表现为线性项，数值稳定性良好。

#### 应用位置

$A_{\mathrm{eff}}(D)$ 进入后续 Butler-Volmer 电流密度修正，存在两种等价处理：
- **局部电流密度放大**：局部体电流密度 $j$（$\mathrm{A/m^2}$）在损伤区被放大（相同总电流通过更小有效面积），$j_{\mathrm{loc}}=j_{\mathrm{macro}}/A_{\mathrm{eff}}(D)$；
- **总电流减少**：总电流 $I=j\cdot A_{\mathrm{eff}}(D)$ 减少（脱粘区不贡献电流）。

两种处理的物理等价性由后续 SPMe 章节给出，本文在弱形式中采用第二种（总电流修正），以避免局部放大导致的数值刚性。

#### 知识库引用

本节反应面积退化模型应参考知识库 `10_宏观损失对电化学参数影响文献调研.md`。若知识库不存在或参数缺失，标注"待补充"。

### 3.6.3 损伤→等效 $n$ 方向热导率

**与 §2.8.2 的关系说明**：本节对 §2.8.2 式 (2.69) 的线性插值模型进行**精细化**，改用基于界面接触热阻的串联模型——线性插值是粗粒度近似（仅两端点严格），串联模型基于物理（每个 $\Gamma_{\mathrm{coh}}$ 界面引入接触热阻 $R_{tc}(D)$），适用于多层堆叠的精细求解。两种形式在 $D=0$ 与 $D=1$ 处退化一致。

#### 物理动机

脱粘引入界面热阻，从而 $n$ 方向有效热导率降低：
- 完好层叠：$n$ 方向热导 $k_n^{\mathrm{bulk}}$（8 层串联 Reuss 平均，§2.4 类比）；
- 部分脱粘：每个 $\Gamma_{\mathrm{coh}}$ 界面引入接触热阻 $R_{tc}(D)$；
- 完全脱粘：界面热阻极大，$n$ 方向热导趋近 0。

#### 公式

**接触热阻**（单位面积温降/热流，量纲 $\mathrm{K\cdot m^2/W}$）：

$$
\boxed{\;R_{tc}(D)\;=\;R_{tc,0}\;+\;D\,\bigl(R_{tc,\max}-R_{tc,0}\bigr)\;}\qquad (3.105)
$$

其中：
- $R_{tc,0}$ 为初始接触热阻（$D=0$ 时），典型量级 $\mathcal{O}(10^{-5}\,\mathrm{K\cdot m^2/W})$（仅给量级）；
- $R_{tc,\max}$ 为完全脱粘接触热阻（$D=1$ 时），典型量级 $\mathcal{O}(10^{-3}\,\mathrm{K\cdot m^2/W})$（仅给量级）；
- $R_{tc,0}\ll R_{tc,\max}$，保证 $D=1$ 时 $n$ 方向热导趋近 0。

**假设 3.6.3（接触热阻线性演化）**：接触热阻随 $D$ 线性插值。线性形式为最简选择，可扩展为指数 $R_{tc}(D)=R_{tc,0}\exp(\alpha D)$（$\alpha>0$）等非线性形式，但本文采用线性以简化 Jacobian 推导。

**$n$ 方向有效热导率**（串联模型：体热导 + $N_{\mathrm{face,repeat}}^{\mathrm{coh}}$ 个真实 cohesive 面热阻）：

$$
\boxed{\;k_n^{\mathrm{eff}}(D)\;=\;\left(\frac{1}{k_n^{\mathrm{bulk}}}\;+\;\frac{N_{\mathrm{face,repeat}}^{\mathrm{coh}}\,R_{tc}(D)}{L}\right)^{-1}\;}\qquad (3.106)
$$

其中：
- $k_n^{\mathrm{bulk}}$ 为 8 层 Reuss 平均体热导（来自 §2.4 均匀化类比）；
- $N_{\mathrm{face,repeat}}^{\mathrm{coh}}=4$ 为一个 8 层重复单元内的真实 cohesive 面数；两种 `interface_type` 各出现于双面涂布的两个物理面上（§3.0b）；
- $L$ 为卷芯径向总厚度，$L=r_{JR,\max}-r_0$（即 $t_{\mathrm{repeat}}$）。

> **术语修订说明**：$N_{\mathrm{lay}}=8$ 是材料层数，$N_{\mathrm{type}}^{\mathrm{coh}}=2$ 是本构类型数，二者都不是串联热阻所需的真实面数。式 (3.106)–(3.114) 的乘数统一使用 $N_{\mathrm{face,repeat}}^{\mathrm{coh}}=4$。其余非 cohesive 相邻面不附加该损伤相关热阻。

**推导（串联热阻叠加）**：$n$ 方向单位面积热流 $q_n$ 通过体材料与四个真实 cohesive 面，总温降为两者之和：

$$
\Delta T\;=\;\frac{q_n\,L}{k_n^{\mathrm{bulk}}}\;+\;N_{\mathrm{face,repeat}}^{\mathrm{coh}}\,q_n\,R_{tc}(D)\qquad (3.107)
$$

等效热导率定义为 $q_n=-k_n^{\mathrm{eff}}\,\Delta T/L$，故

$$
\frac{\Delta T}{q_n\,L}\;=\;\frac{1}{k_n^{\mathrm{bulk}}}\;+\;\frac{N_{\mathrm{face,repeat}}^{\mathrm{coh}}\,R_{tc}(D)}{L}\;=\;\frac{1}{k_n^{\mathrm{eff}}(D)}\qquad (3.108)
$$

即式 (3.106)。$\square$

**物理诠释**：式 (3.106) 将 $n$ 方向热传导等效为“体热阻 + 四个真实 cohesive 面热阻”串联。脱粘使 $R_{tc}(D)$ 增大，串联总热阻随之增大，故 $k_n^{\mathrm{eff}}$ 减小。$D=1$ 时 $R_{tc}=R_{tc,\max}\gg L/(N_{\mathrm{face,repeat}}^{\mathrm{coh}}k_n^{\mathrm{bulk}})$，$k_n^{\mathrm{eff}}\to L/(N_{\mathrm{face,repeat}}^{\mathrm{coh}}R_{tc,\max})\to0$。

#### 偏导数

由链式法则：

$$
\frac{\partial k_n^{\mathrm{eff}}}{\partial D}\;=\;\frac{\partial k_n^{\mathrm{eff}}}{\partial R_{tc}}\cdot\frac{\partial R_{tc}}{\partial D}\qquad (3.109)
$$

**步骤 1**（$\partial R_{tc}/\partial D$）：由式 (3.105)，

$$
\frac{\partial R_{tc}}{\partial D}\;=\;R_{tc,\max}-R_{tc,0}\qquad (3.110)
$$

**步骤 2**（$\partial k_n^{\mathrm{eff}}/\partial R_{tc}$）：由式 (3.106)，记

$$
\frac{1}{k_n^{\mathrm{eff}}(D)}\;=\;\frac{1}{k_n^{\mathrm{bulk}}}\;+\;\frac{N_{\mathrm{face,repeat}}^{\mathrm{coh}}\,R_{tc}(D)}{L}\qquad (3.111)
$$

两边对 $R_{tc}$ 求导：

$$
-\frac{1}{(k_n^{\mathrm{eff}})^2}\cdot\frac{\partial k_n^{\mathrm{eff}}}{\partial R_{tc}}\;=\;\frac{N_{\mathrm{face,repeat}}^{\mathrm{coh}}}{L}\qquad (3.112)
$$

故

$$
\frac{\partial k_n^{\mathrm{eff}}}{\partial R_{tc}}\;=\;-\frac{N_{\mathrm{face,repeat}}^{\mathrm{coh}}}{L}\,(k_n^{\mathrm{eff}})^2\qquad (3.113)
$$

**合并**（代入式 3.109）：

$$
\boxed{\;\frac{\partial k_n^{\mathrm{eff}}}{\partial D}\;=\;-(k_n^{\mathrm{eff}})^2\cdot\frac{N_{\mathrm{face,repeat}}^{\mathrm{coh}}\bigl(R_{tc,\max}-R_{tc,0}\bigr)}{L}\;}\qquad (3.114)
$$

**量纲**：$\mathrm{W/(m\cdot K)}$（与 $k_n^{\mathrm{eff}}$ 同量纲）。

**量纲一致性自查**：
- $(k_n^{\mathrm{eff}})^2$ 量纲 $\mathrm{W^2/(m^2\cdot K^2)}$；
- $N_{\mathrm{face,repeat}}^{\mathrm{coh}}/L$ 量纲 $\mathrm{m^{-1}}$（$N_{\mathrm{face,repeat}}^{\mathrm{coh}}=4$ 无量纲）；
- $R_{tc,\max}-R_{tc,0}$ 量纲 $\mathrm{K\cdot m^2/W}$；
- 乘积：$\mathrm{W^2/(m^2\cdot K^2)\times m^{-1}\times K\cdot m^2/W}=\mathrm{W/(m\cdot K)}$，与 $k_n^{\mathrm{eff}}$ 一致。$\checkmark$

**注 3.6.3（无奇异性）**：$D\in[0,1]$ 时 $k_n^{\mathrm{eff}}$ 始终有限（因 $R_{tc,\max}$ 有限），故 $\partial k_n^{\mathrm{eff}}/\partial D$ 在全区间有界。该支路与 §3.6.2（反应面积）同为无奇异性的线性-缓变支路，数值稳定性良好。

**注 3.6.4（偏导数对 $k_n^{\mathrm{eff}}$ 的隐式依赖）**：式 (3.114) 含 $(k_n^{\mathrm{eff}})^2$ 项，而 $k_n^{\mathrm{eff}}$ 本身为 $D$ 的函数（式 3.106）。该"偏导显含函数本身"的结构源于 $k_n^{\mathrm{eff}}=f(R_{tc})^{-1}$ 的反演形式，在全耦合 Jacobian 组装时应使用当前增量步的 $k_r^{\mathrm{eff},n}$ 值代入，而非线性化展开。

#### 应用位置

$k_n^{\mathrm{eff}}(D)$ 进入 2D 热传导方程的 $n$ 方向系数修正（spec G6 修订）：

$$
\rho c_p\frac{\partial T}{\partial t}\;=\;\frac{\partial}{\partial s}\!\left(k_s\frac{\partial T}{\partial s}\right)\;+\;\frac{\partial}{\partial n}\!\left(k_n^{\mathrm{eff}}(D)\frac{\partial T}{\partial n}\right)\;+\;Q\qquad (3.115)
$$

其中切向热导 $k_s$ 不穿越材料界面（切向热传导沿层内方向），故仅 $k_n^{\mathrm{eff}}$ 依赖 $D$。详细由 §5.1 (s, n) 平面热传导方程给出。

#### 知识库引用

本节接触热阻模型应参考知识库 `10_宏观损失对电化学参数影响文献调研.md`（同 §3.6.2）。若知识库不存在或参数缺失，标注"待补充"。

### 3.6.4 接口传递（三条偏导数汇总与 Jacobian 衔接）

#### 三条偏导数汇总

下表汇总三条标量反馈支路的函数形式、偏导数、量纲与应用位置（spec S-IV 强制核心闭环）：

| 反馈支路 | 标量函数 | 偏导数 $\partial/\partial D$ | 量纲 | 应用位置 |
|---|---|---|---|---|
| 接触电阻 | $R_{\mathrm{contact}}(D)=R_{\mathrm{contact},0}/(1-D)$（式 3.95） | $R_{\mathrm{contact},0}/(1-D)^2$（式 3.99） | $\Omega\cdot\mathrm{m}^2$ | Butler-Volmer 过电位（式 3.101） |
| 反应面积 | $A_{\mathrm{eff}}(D)=A_0(1-D)$（式 3.102） | $-A_0=-1$（式 3.104） | 无量纲 | 电流密度修正（§3.6.2 应用位置） |
| $n$ 方向热导 | $k_n^{\mathrm{eff}}(D)$（式 3.106） | $-(k_n^{\mathrm{eff}})^2\cdot N_{\mathrm{face,repeat}}^{\mathrm{coh}}(R_{tc,\max}-R_{tc,0})/L$（式 3.114，$N_{\mathrm{face,repeat}}^{\mathrm{coh}}=4$） | $\mathrm{W/(m\cdot K)}$ | 2D 热传导 $n$ 方向系数（式 3.115） |

#### 三条偏导在全耦合 Jacobian 中的位置

全耦合弱形式的 Jacobian 矩阵含 $D$ 作为未知量（spec S1）。三条偏导的具体位置：

**位置 1（$\partial R_{\mathrm{contact}}/\partial D$）**：进入 Jacobian 块 $\partial\mathbf{R}_{\mathrm{SPMe}}/\partial D$（电化学残差对损伤的偏导）。式 (3.101) 中 $\eta_{\mathrm{eff}}$ 对 $D$ 的偏导含两项：

$$
\frac{\partial \eta_{\mathrm{eff}}}{\partial D}\;=\;-j\cdot\frac{\partial}{\partial D}\!\left(\frac{R_{\mathrm{contact}}(D)}{A_{\mathrm{eff}}(D)}\right)\qquad (3.116)
$$

代入式 (3.95)、(3.102) 后应用商法则：

$$
\frac{\partial}{\partial D}\!\left(\frac{R_{\mathrm{contact}}}{A_{\mathrm{eff}}}\right)\;=\;\frac{\dfrac{R_{\mathrm{contact},0}}{(1-D)^2}\cdot A_0(1-D)\;-\;\dfrac{R_{\mathrm{contact},0}}{1-D}\cdot(-A_0)}{\bigl[A_0(1-D)\bigr]^2}\qquad (3.117)
$$

化简：

$$
\frac{\partial}{\partial D}\!\left(\frac{R_{\mathrm{contact}}}{A_{\mathrm{eff}}}\right)\;=\;\frac{R_{\mathrm{contact},0}}{A_0}\cdot\frac{2}{(1-D)^3}\qquad (3.118)
$$

**位置 2（$\partial A_{\mathrm{eff}}/\partial D$）**：同样进入 Jacobian 块 $\partial\mathbf{R}_{\mathrm{SPMe}}/\partial D$（与位置 1 共同构成电化学残差的损伤列），已在式 (3.117) 中合并。

**位置 3（$\partial k_n^{\mathrm{eff}}/\partial D$）**：进入 Jacobian 块 $\partial\mathbf{R}_{\mathrm{thermal}}/\partial D$（热残差对损伤的偏导）。由式 (3.115)，热残差中 $k_n^{\mathrm{eff}}(D)$ 出现在空间算子系数中，其偏导贡献为：

$$
\frac{\partial \mathbf{R}_{\mathrm{thermal}}}{\partial D}\bigg|_{\text{via }k_n^{\mathrm{eff}}}\;=\;-\frac{\partial}{\partial n}\!\left(\frac{\partial k_n^{\mathrm{eff}}}{\partial D}\,\frac{\partial T}{\partial n}\right)\qquad (3.119)
$$

代入式 (3.114) 即得显式表达。

**注 3.6.5（Jacobian "损伤列"结构）**：上述三块偏导构成全耦合 Jacobian 的"损伤列"——即所有残差方程对损伤未知量 $D$ 的偏导。该列结构为 operator-splitting 或 fully-coupled 求解策略的关键输入：
- **operator-splitting**：仅需"损伤列"的对角块（$\partial\mathbf{R}_{\mathrm{SPMe}}/\partial D$ 与 $\partial\mathbf{R}_{\mathrm{thermal}}/\partial D$ 分别独立求解 $D$ 更新）；
- **fully-coupled**：需完整"损伤列"参与 Newton 步。

两种策略的选择由后续全耦合章节给出，本节仅提供偏导数的显式形式。

#### Curie 兼容性自查

- 三条支路均为**标量函数**（输入 $D$ 标量，输出标量），Curie 兼容；$\checkmark$
- 无任何矢量/张量直接进入电化学/热方程（CZM 矢量 traction 经 §3.5 凝聚为 $D$ 后才进入下游）；$\checkmark$
- 与 §2.8.2 定理 2.8.1 严格一致（标量-标量映射允许，矢量-标量映射禁止）。$\checkmark$

**推论 3.6.1（Curie 桥梁的闭环完整性）**：§3.5 将 CZM 矢量输出 $\mathbf{T}=(T_n,T_s)$ 经损伤演化（式 3.91）凝聚为标量 $D$，§3.6 将 $D$ 经三条标量支路反馈到电化学/热方程，构成完整的"矢量→标量→标量"Curie 兼容桥梁。该桥梁是全耦合框架中力学（矢量）与电化学/热（标量场）耦合的唯一通道。

#### 与已有记号的闭环检验（R8 自查）

- $R_{\mathrm{contact}},A_{\mathrm{eff}}$：§2.8.2 定义（式 2.67–2.68），本节给出 $D$ 依赖关系（式 3.95、3.102）——函数化具体形式，无重复定义；
- $k_n^{\mathrm{eff}}$：§2.8.2 式 (2.69) 给出粗粒度线性插值；§3.6.3 式 (3.106) 给出基于接触热阻的串联精细化形式，作为式 (2.69) 的**替代实现**——两种形式在 $D=0$/$D=1$ 端点一致，本文数值实现采用 §3.6.3 串联形式；
- $D$：§3.3.1 定义（式 3.28），本节作为输入引用——无重复定义；
- $R_{tc},R_{tc,0},R_{tc,\max},N_{\mathrm{face,repeat}}^{\mathrm{coh}},L$：本节首次定义（式 3.105、3.106），其中 $N_{\mathrm{face,repeat}}^{\mathrm{coh}}=4$ 为一个 8 层重复单元内真实箔–涂层面的串联计数；
- $R_{\mathrm{contact},0}$：本节首次定义（式 3.95），为初始接触电阻面电阻率，未在先前章节出现；
- $A_0$：本节首次定义（式 3.96），为完好界面接触面积基准（约定 $A_0\equiv 1$），未在先前章节出现。

R8 自查通过。

#### 知识库引用汇总

本节涉及的损伤-电化学/热反馈模型参数 $(R_{\mathrm{contact},0},R_{tc,0},R_{tc,\max})$ 与结构参数 $(N_{\mathrm{face,repeat}}^{\mathrm{coh}},L)$ 的具体取值应来自知识库 `09_损伤引起接触电阻文献调研.md` 与 `10_宏观损失对电化学参数影响文献调研.md`（其中 $N_{\mathrm{face,repeat}}^{\mathrm{coh}}=4$ 由 §3.0b 的 8 层拓扑确定，非自由参数）。若知识库不存在或参数缺失，标注"待补充"。本节仅给出符号推导，不代入具体数值（量级估计标注"仅给量级"）。

---
