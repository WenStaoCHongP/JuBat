# 颗粒尺度 vs 宏观尺度力学无量纲化对比

> **目的**：澄清 JuBat 中两个力学尺度（颗粒 / 极片-界面）在无量纲化策略上的本质差异，避免跨尺度混用导致的数值失真。
> **关联**：`md/15_颗粒与极片模量区分.md`、`src/SetParams.jl`、`src/Mechanical.jl`、`src/CouplingState.jl`
>
> **⚠ 更新（2026-07-22 重设计 v2 已实施）**：本文 §3.1 中 `δ_czm = L` 的描述已过时。
> δ_czm 现锚定为 `2·G_c_pe_pcc/σ_max_pe_pcc`（断裂能定义的临界分离，锚定界面 δ_c* ≡ 1），
> 装配层引入换算因子 Λ = L/δ_czm。同时修复了体/界面应力空间混装（`moduli_of` 双重再缩放）
> 与输出还原尺度误用（分离 ×r0→×δ_czm，牵引 ×E_n/E_p→×σ_czm）。
> 见 `docs/planning-with-files/力学模块修改/宏观力学模块无量纲化重设计.md`。

---

## 1. 两个尺度的物理对象

| 维度 | 颗粒尺度 | 宏观尺度 |
|------|---------|---------|
| 物理对象 | 单个活性物质颗粒（AM particle） | 极片涂层 + 集流体 + 隔膜 + 层间内聚力单元 |
| 几何尺度 | 颗粒半径 `rs` ~ 5–10 μm | 螺旋叠合厚度 `L` ~ 100–200 μm；CZM 界面厚度 → 0 |
| 应力来源 | 锂浓度梯度（颗粒内扩散）引起的本征应变 | 热膨胀 `α·ΔT` + 宏观浓度应变 `β·Δsoc` + 界面分离 `δ` |
| 典型模量 | `PE.E`/`NE.E` ~ 10–100 GPa | `PE.E_coat`/`NE.E_coat` ~ 0.1–1 GPa；`SP.E`/`PCC.E`/`NCC.E` ~ 1–100 GPa |
| 典型应力水平 | 1–100 MPa（扩散应力） | 10–200 MPa（CZM 峰值 `σ_max` ~ 80 MPa） |
| 求解入口 | `Calstressdisp`（`src/Mechanical.jl:112`） | `thermal_diffusion_stress_2D` + `assemble_coupled_system`（`src/Mechanical.jl:165`, `src/czm.jl`） |

---

## 2. 颗粒尺度的无量纲化

### 2.1 参考量（在 `ChooseCell` 中设定）

定义于 `src/SetParams.jl:309-310`：

$$
E_{\mathrm{scale},p} \;=\; c_{s,\max,p}\,R\,T_{\mathrm{ref}}, \qquad
E_{\mathrm{scale},n} \;=\; c_{s,\max,n}\,R\,T_{\mathrm{ref}}
$$

量纲验证：`[mol/m³] · [J/(mol·K)] · [K] = [J/m³] = [Pa]`，与弹性模量同量纲。物理意义即**理想溶液的化学压参考**，与扩散应力的生成机制同源。

### 2.2 字段归一化（在 `NormaliseParam` 中执行）

`src/SetParams.jl:364-369, 382-396`：

$$
E_p^{*} = \frac{E_p}{E_{\mathrm{scale},p}}, \qquad
E_n^{*} = \frac{E_n}{E_{\mathrm{scale},n}}
$$

$$
c_s^{*} = \frac{c_s}{c_{s,\max}}, \qquad
\Omega^{*} = \Omega \cdot c_{s,\max}
$$

注意 `Omega`（偏摩尔体积 `[m³/mol]`）**乘以** `cs_max` 而非除，使其无量纲化为 O(1) 数（典型值 10⁻⁶ m³/mol × 3×10⁴ mol/m³ ≈ 0.03）。

### 2.3 颗粒扩散应力的无量纲形式

原始物理公式（球形颗粒）：

$$
\sigma_r(0) = \frac{2\,\Omega\,E}{9(1-\nu)}\,(\bar{c}_s - c_s(0)), \qquad
\sigma_\theta(r_s) = \frac{\Omega\,E}{3(1-\nu)}\,(\bar{c}_s - c_s(r_s))
$$

代入归一化变量（`src/Mechanical.jl:134-135` 实际计算式）：

$$
\sigma_r^{*}(0) \;=\; \frac{2\,\Omega^{*}\,E_p^{*}}{9(1-\nu)}\,(\bar{c}_s^{*} - c_s^{*}(0))
\;=\; \frac{\sigma_r(0)}{E_{\mathrm{scale},p}}
$$

$$
\sigma_\theta^{*}(r_s) \;=\; \frac{\Omega^{*}\,E_p^{*}}{3(1-\nu)}\,(\bar{c}_s^{*} - c_s^{*}(r_s))
\;=\; \frac{\sigma_\theta(r_s)}{E_{\mathrm{scale},p}}
$$

**关键性质**：颗粒尺度下 **E 与 σ 共享同一参考量** `E_{\mathrm{scale}}`。归一化后 `E_p*`、`σ_r*`、`σ_θ*` 都是 O(1) 数，刚度矩阵良态。

---

## 3. 宏观尺度的无量纲化

### 3.1 参考量（在 `ChooseCell` 中设定）

定义于 `src/SetParams.jl:311-339`：

$$
E_{\mathrm{coat}} \;=\; \frac{2\,t_{PE}\,E_{PE}^{\mathrm{coat}} + 2\,t_{NE}\,E_{NE}^{\mathrm{coat}} + 2\,t_{SP}\,E_{SP} + t_{PCC}\,E_{PCC} + t_{NCC}\,E_{NCC}}{2\,t_{PE} + 2\,t_{NE} + 2\,t_{SP} + t_{PCC} + t_{NCC}}
$$

即"全叠合厚度加权平均"模量，量纲 `[Pa]`。物理意义是**电极叠合层的等效宏观刚度参考**。

$$
\sigma_{\mathrm{czm}} \;=\; \sigma_{\max}^{PE\text{-}PCC}, \qquad
\delta_{\mathrm{czm}} \;=\; L, \qquad
G_{\mathrm{czm}} \;=\; \sigma_{\mathrm{czm}}\,\delta_{\mathrm{czm}}, \qquad
K_{\mathrm{czm}} \;=\; \frac{\sigma_{\mathrm{czm}}}{\delta_{\mathrm{czm}}}
$$

锚点选取：用 PE-PCC 界面的法向强度 `σ_max_pe_pcc` 作为 CZM 应力归一化参考（自锚点策略，`σ_max_pe_pcc* ≡ 1`）。

### 3.2 字段归一化（在 `NormaliseParam` 中执行）

**体模量类**（`src/SetParams.jl:365, 392, 415-443`）：

$$
E_{coat,p}^{*} = \frac{E_{coat,p}}{E_{\mathrm{coat}}}, \quad
E_{coat,n}^{*} = \frac{E_{coat,n}}{E_{\mathrm{coat}}}, \quad
E_{SP}^{*} = \frac{E_{SP}}{E_{\mathrm{coat}}}, \quad
E_{PCC}^{*} = \frac{E_{PCC}}{E_{\mathrm{coat}}}, \quad
E_{NCC}^{*} = \frac{E_{NCC}}{E_{\mathrm{coat}}}
$$

**CZM 本构类**（`src/SetParams.jl:481-501`），以 PE-PCC 为例（NE-NCC 同构）：

$$
\sigma_{\max}^{*} = \frac{\sigma_{\max}}{\sigma_{\mathrm{czm}}}, \quad
\delta_{0}^{*} = \frac{\delta_0}{\delta_{\mathrm{czm}}}, \quad
\delta_c^{*} = \frac{\delta_c}{\delta_{\mathrm{czm}}}, \quad
G_c^{*} = \frac{G_c}{G_{\mathrm{czm}}}, \quad
K_n^{*} = \frac{K_n}{K_{\mathrm{czm}}}
$$

Mode II（切向）同样 5 字段：`τ_max*`、`δ_0_t*`、`δ_c_t*`、`G_c_t*`、`K_t*`。两界面合计 20 个 CZM 本构字段。

### 3.3 派生量：CZM 内部 `E_eff`（"双重再缩放"）

`src/CouplingState.jl:306-307`：

$$
E_{\mathrm{eff}}^{CZM} \;=\; E_{coat,p}^{*} \;\times\; \frac{E_{\mathrm{coat}}}{\sigma_{\mathrm{czm}}}
\;=\; \frac{E_{coat,p}}{\sigma_{\mathrm{czm}}}
$$

这一步把 **E_coat-normalized** 的涂层模量转换到 **σ_czm-normalized** 的 CZM 内部空间，供 `assemble_bulk_stiffness` 与 CZM 牵引-分离律在同一参考下耦合。

数值示例（Jellyroll 参数集）：
- `E_coat_p` ≈ 5×10⁸ Pa
- `σ_czm` ≈ 8×10⁷ Pa
- `E_eff^{CZM}` ≈ 6.25（无量纲，O(1)）

---

## 4. 核心差异

| 维度 | 颗粒尺度 | 宏观尺度 |
|------|---------|---------|
| E 参考量 $E_{\mathrm{scale}}$ | $c_{s,\max}\,R\,T_{\mathrm{ref}}$（化学压参考） | $E_{\mathrm{coat}}$（厚度加权宏观模量） |
| σ 参考量 $\sigma_{\mathrm{scale}}$ | **= $E_{\mathrm{scale}}$**（与 E 同源） | **= $\sigma_{\max}^{PE\text{-}PCC}$**（界面强度独立锚点） |
| E 与 σ 是否同参考 | ✅ 是 | ❌ 否（解耦） |
| 几何尺度耦合 | 通过 $r_s$ 隐式（颗粒球坐标） | 显式 `δ_czm = L`，`K_czm = σ_czm / L` |
| 派生量 | 无（颗粒公式中 E 直接驱动应力） | 有：`E_eff^{CZM}` 需双重再缩放 |
| 输出应力还原 | $\sigma_{\mathrm{dim}} = \sigma^{*} \cdot E_{\mathrm{scale}}$ | $\sigma_{\mathrm{dim}} = \sigma^{*} \cdot \sigma_{\mathrm{czm}}$ |

### 4.1 为什么颗粒尺度 E 与 σ 共享参考

扩散应力物理上由化学势梯度驱动，本构为 $\sigma \propto E \cdot \Omega \cdot \Delta c$。用 $c_{s,\max}\,R\,T$ 作为参考时：

$$
\sigma^{*} \;\propto\; \frac{E \cdot \Omega \cdot \Delta c}{c_{s,\max}\,R\,T}
\;=\; E^{*} \cdot \Omega^{*} \cdot \Delta c^{*}
$$

三项归一化因子自然抵消（$c_{s,\max}$ 被 $\Omega^{*}=\Omega\,c_{s,\max}$ 与 $\Delta c^{*}=\Delta c/c_{s,\max}$ 内部吸收），结果是无量纲 O(1) 数。

### 4.2 为什么宏观尺度 E 与 σ 解耦

宏观尺度有两个**独立测量**的物理量：
1. **涂层模量** $E_{coat}$：由拉伸实验测得，反映极片整体刚度
2. **界面强度** $\sigma_{\max}$：由脱层实验测得，反映涂层-集流体结合力

二者没有物理绑定关系（$\sigma_{\max}/E_{coat}$ 可以是 0.01 到 1 之间任意值）。若强行用 $E_{coat}$ 作为 $\sigma$ 参考量：
- 对硬界面（$\sigma_{\max} \sim E_{coat}$）：CZM 本构在归一化空间 O(1)，正常
- 对软界面（$\sigma_{\max} \ll E_{coat}$）：$\sigma_{\max}^{*} \ll 1$，CZM 本构几乎全部落在弹性区，损伤变量不演化，数值上"看不见脱层"

**自锚点策略**（$\sigma_{\mathrm{czm}} \equiv \sigma_{\max}^{PE\text{-}PCC}$）保证 CZM 应力-分离曲线在归一化空间峰值恒为 1，与材料硬度无关，本构数值稳定。代价是引入双重再缩放：体刚度项的 $E_{coat}^{*}$ 必须乘 $E_{\mathrm{coat}}/\sigma_{\mathrm{czm}}$ 才能与 CZM 牵引同空间装配。

---

## 5. 跨尺度混用的后果

### 5.1 若 CZM 用 $E_{\mathrm{coat}}$ 归一化 σ

对 Jellyroll（$E_{coat} \approx 8\times 10^{10}$ Pa，$\sigma_{\max} \approx 8\times 10^{7}$ Pa）：

$$
\sigma_{\max}^{*} \;=\; \frac{\sigma_{\max}}{E_{\mathrm{coat}}} \;\approx\; 10^{-3}
$$

CZM 损伤演化方程中 `σ/σ_max` 在加载初期远小于 1，损伤变量 D 长期停在 0，**脱层永不发生**。同时刚度矩阵 `K_n* = K_n / (E_{coat}/L) ≈ 10^{-3}`，条件数恶化。

### 5.2 若颗粒扩散用 $\sigma_{\mathrm{czm}}$ 归一化

颗粒扩散应力 $\sigma_{\mathrm{diff}} \sim 10$ MPa，$\sigma_{\mathrm{czm}} \sim 80$ MPa：

$$
\sigma_{\mathrm{diff}}^{*} \;\approx\; 0.1
$$

而 $E_p^{*}$ 若仍按 $E_p / E_{\mathrm{scale},p}$（典型 O(1)），则扩散应力公式 $\sigma^{*} \propto E^{*}\,\Omega^{*}\,\Delta c^{*}$ 输出与"σ 参考为 $\sigma_{\mathrm{czm}}$"不一致，导致**输出量纲错位**：保存到 `variables` 的应力数值再乘 `σ_czm` 还原时，与实际应力差一个 $E_{\mathrm{scale},p}/\sigma_{\mathrm{czm}} \approx 10^{-2}$ 因子。

### 5.3 防御机制

| 位置 | 检查 |
|------|------|
| `ChooseCell`（`SetParams.jl:312-314`） | `PE.E_coat == 0` 触发 `@warn`，`scale.E_coat = 0` |
| `ChooseCell`（`SetParams.jl:333-334`） | `σ_max_pe_pcc == 0` 触发 `@warn`，CZM 归一化将产生 `Inf/NaN` |
| `compute_czm_params_per_interface`（`CouplingState.jl:297-303`） | 7 个 `@assert` 拦截零/负模量、零 `σ_czm`、零 `E_coat` |
| `thermal_diffusion_stress_2D`（`Mechanical.jl:167`） | 入口 `@assert PE.E_coat > 0 && NE.E_coat > 0` |

---

## 6. 字段速查表

| `Scale` 字段 | 表达式 | 量纲 | 服务对象 |
|------|------|------|------|
| `E_p` | `PE.cs_max · R · T_ref` | Pa | 颗粒扩散应力（正极） |
| `E_n` | `NE.cs_max · R · T_ref` | Pa | 颗粒扩散应力（负极） |
| `E_coat` | 全叠合厚度加权平均 | Pa | 宏观体模量归一化 |
| `σ_czm` | `cohesive.σ_max_pe_pcc` | Pa | CZM 应力归一化（自锚点） |
| `δ_czm` | `scale.L` | m | CZM 分离位移归一化 |
| `G_czm` | `σ_czm · δ_czm` | J/m² | CZM 断裂能归一化 |
| `K_czm` | `σ_czm / δ_czm` | Pa/m | CZM 刚度归一化 |

---

## 7. 总结

**颗粒尺度**走"化学压参考"路线：用 $c_{s,\max}\,R\,T$ 同时归一化 $E$ 和 $\sigma$，扩散应力公式 $\sigma^{*} \propto E^{*}\,\Omega^{*}\,\Delta c^{*}$ 在归一化空间自然闭合，无派生量。

**宏观尺度**走"界面强度自锚点"路线：$E_{coat}$ 与 $\sigma_{\max}$ 解耦归一化，各自 O(1)；CZM 装配时通过 $E_{\mathrm{eff}}^{CZM} = E_{coat,p}^{*} \cdot E_{\mathrm{coat}}/\sigma_{\mathrm{czm}}$ 把体刚度项重缩放到 CZM 应力空间，保证切向-法向-体刚度三者同参考。

两套尺度的参考量之间**不存在简单的换算关系**，代码中通过 `Scale` 结构体的不同字段隔离，并在入口处用 `@assert` 防御跨尺度误用。
