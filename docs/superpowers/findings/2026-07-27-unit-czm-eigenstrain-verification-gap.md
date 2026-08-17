# 单元级 CZM 条带验证：本征应变工况未通过原因与参数调研议程

> **日期**：2026-07-27  
> **关联**：`docs/superpowers/specs/2026-07-23-unit-czm-strip-verification-design.md`  
> **脚本**：`test/unit_czm_bilinear.jl`、`test/unit_czm_eigenstrain.jl`  
> **参数**：`src/parameters/Jellyroll.jl`  
> **后续**：文献调查合适的极片模量 / 界面 CZM 参数（见 §5）

---

## 1. 验证范围与“未通过”指什么

| 脚本 | 目标 | 结果 |
|------|------|------|
| `unit_czm_bilinear.jl` | 位移 BC 驱动，验证双线性 T–δ（Mode I / 卸载 / Mode II） | **通过**（本构与装配路径正常） |
| `unit_czm_eigenstrain.jl` | 合成 \(\Delta T/\Delta\mathrm{soc}\) + 固定端，弹性段/损伤段与物理量级一致 | **在“实际量级参数”下不能自然进入损伤** |

因此“单元级验证未通过”主要指：**在贴近实际的涂层模量 + 实际 SOC 变化 + 原 Jellyroll CZM 强度下，本征应变驱动无法使界面达到起裂**；并非网格生成或双线性本构实现错误。

为演示损伤，脚本曾对 CZM 做临时 retune（\(\sigma_{\max}\times 0.1\)，\(\delta_c/\delta_0=50\)），这属于**验证用降强度**，不能当作生产参数已自洽。

---

## 2. 当前问题的力学抽象（一句话）

**平面应力多层弹性条带（8 层 Q4）+ 4 个双线性 COH2D4 界面；热/化学本征应变驱动；集流体与顶边固支、全体 \(u_x=0\)。**  
问：在 \(\Delta T=10\,\mathrm{K}\)、\(\Delta\mathrm{soc}=0.8\) 下，NE–NCC / PE–PCC 能否出现 Mode I 张开损伤（\(\delta_n>\delta_0\Rightarrow D>0\)）。

### 2.1 几何（物理）

- 参考长度 \(L \approx 1.73\times 10^{-4}\,\mathrm{m}\)
- 条带 \(W\approx H\approx 3.74\times 10^{-4}\,\mathrm{m}\)
- 层序：PE → PCC → PE → SP → NE → NCC → NE → SP  
  厚度约：75.6 / 16 / 75.6 / 12 / 85.2 / 12 / 85.2 / 12 µm

### 2.2 体材料（当前 Jellyroll）

| 材料 | \(E\) | \(\nu\) |
|------|-------|---------|
| PE / NE 涂层 | **500 MPa** | 0.30 / 0.25 |
| SP / PCC / NCC | 500 MPa | 0.30 |

本征应变（生产契约）：

\[
\varepsilon_0=\alpha_{\mathrm{eff}}\Delta T+\beta_n\Delta\mathrm{soc}_n+\beta_p\Delta\mathrm{soc}_p
\]

在 \(\Delta T=10\,\mathrm{K}\)、\(\Delta\mathrm{soc}=0.8\) 时量级约为：

| 层 | \(\varepsilon_T\) | \(\varepsilon_{\mathrm{chem}}\) | \(\varepsilon_0\) |
|----|-------------------|----------------------------------|-------------------|
| PE | \(\sim 1.5\times 10^{-4}\) | \(\sim -1.2\%\) | \(\sim -1.2\%\)（收缩） |
| NE | \(\sim 1.5\times 10^{-4}\) | \(\sim +2.7\%\) | \(\sim +2.8\%\)（膨胀） |
| 其它 | \(\sim 1.5\times 10^{-4}\) | 0 | \(\sim 0.015\%\) |

### 2.3 界面 CZM（Jellyroll 原始 vs 脚本 retune）

| 量 | PE–PCC 原始 | NE–NCC 原始 | 脚本 retune（×0.1, \(\delta_c/\delta_0=50\)） |
|----|-------------|-------------|-----------------------------------------------|
| \(\sigma_{\max}\) | 82 MPa | 92 MPa | ≈ 8.2 / 9.2 MPa |
| \(K_n\) | \(2.4\times 10^{17}\) Pa/m | \(1.2\times 10^{17}\) Pa/m | 由 \(\sigma_{\max}/\delta_0\) 重算 |
| \(G_c\) | 25.3 J/m² | 6.2 J/m² | 随 \(\sigma_{\max}\) 同比缩小（δc 名义保留策略下） |

### 2.4 边界条件（本征应变脚本当前设定）

1. 全体节点 \(u_x=0\)
2. PCC、NCC 全部节点 \(u_y=0\)
3. 顶边 \(u_y=0\)
4. 底边 \(u_y\) 自由  

准静态斜坡；起裂后可 hold 载荷以改善收敛。

---

## 3. 未通过的根因（定量）

### 3.1 应力天花板：涂层撑不住界面强度

固定端本征应变问题中，界面可达到的牵引力受涂层刚度与本征应变量级限制：

\[
\sigma \lesssim E_{\mathrm{coat}}\cdot|\varepsilon_0|
\approx 500\,\mathrm{MPa}\times(0.03\text{–}0.04)
\approx 15\text{–}20\,\mathrm{MPa}
\]

而原始起裂强度 \(\sigma_{\max}\approx 82\text{–}92\,\mathrm{MPa}\)。  

**结论**：在 \(E_{\mathrm{coat}}=500\,\mathrm{MPa}\) + 实际 \(\Delta\mathrm{soc}=0.8\) 下，**弹性应力上限低于 CZM 峰值**，无论怎样调整固支组合，\(\delta/\delta_0\) 都到不了 1。

实测扫描（原始 \(\sigma_{\max}\)，\(\Delta\mathrm{soc}=0.8\)）：

| BC | \(\max|\delta|/\delta_0\)（量级） | \(D\) |
|----|-----------------------------------|-------|
| PCC+NCC+顶固支 | \(\lesssim 0.18\) | 0 |
| 底+顶固支 | \(\lesssim 0.06\) | 0 |
| 强冷却 \(\Delta T=-100\sim-300\,\mathrm{K}\) | 仍 \(\ll 1\) | 0 |

### 3.2 一维串联近似（解释“只改 BC 无效”）

\[
\frac{\delta}{\delta_0}
\sim
\frac{\bigl|\sum_i h_i\varepsilon_{0,i}\bigr|}
{\sigma_{\max}\sum_i h_i/E_i+\delta_0}
\]

- \(E\) 变小 → 分母变大 → \(\delta/\delta_0\) 下降（软涂层吸走失配）。  
- 仅降低 \(K_n\)（增大 \(\delta_0\)）对软涂层往往**更难**起裂。  
- 有效杠杆：**提高 \(E_{\mathrm{coat}}\)**、**降低 \(\sigma_{\max}\)**、或 **加大 \(|\varepsilon_0|\)**（超物理 SOC / 大幅冷却）。

### 3.3 其它曾干扰验证的因素（次要）

| 因素 | 现象 | 处理 |
|------|------|------|
| \(\delta_c/\delta_0\sim 10^3\)（原 \(K_n\) 极大） | 一过起裂 Newton 易发散 | 设定 \(\delta_c/\delta_0\sim 20\text{–}50\) 或线搜索 + hold |
| 损伤仅在收敛后 commit | \(\delta>\delta_0\) 但步未收敛 → 打印 \(D=0\) | 减小步长 / hold；勿在未收敛迭代中错误 commit |
| \(\delta_n<0\) | 压缩罚接触，不进 Mode I 损伤 | 需净张开（收缩约束或合适失配符号） |
| 打印误用同一 \(\delta_0\) | PE/NE 界面阈值不同 | 按 `interface_type` 分别归一 |

### 3.4 参数矛盾一览

| 组合 | 能否自然损伤（\(\Delta\mathrm{soc}=0.8\)） |
|------|---------------------------------------------|
| \(E_{\mathrm{coat}}\sim\mathrm{GPa}\) + 原 \(\sigma_{\max}\) | 较易（曾验证） |
| \(E_{\mathrm{coat}}=500\,\mathrm{MPa}\) + 原 \(\sigma_{\max}\) | **否**（应力天花板） |
| \(E_{\mathrm{coat}}=500\,\mathrm{MPa}\) + \(\sigma_{\max}\times 0.1\) | **是**（脚本 retune，非文献标定） |
| 仅改 BC、不改 \(E\)/\(\sigma_{\max}\) | **否** |

---

## 4. 已确认“通过”的部分（避免误判）

1. **网格路径（方案 C）**：`create_unit_czm_strip` → 生产 `create_czm_mesh`，4×COH、界面类型、副本节点断言通过。  
2. **双线性本构**：位移驱动全曲线、卸载、Mode II 测试通过。  
3. **本征应变装配管线**：`assemble_thermal_chemical_load` 与自由涂层解析（\(\varepsilon_x=0\) 时 \(u=-\varepsilon_0(1+\nu)h\)）在弹性、合适 BC 下一致。  
4. **符号约定**：\(\delta_n>0\) 张开可损伤；\(\delta_n<0\) 压缩罚刚度。

失败点集中在：**宏观模量–界面强度–本征应变量级三者未自洽**，而不是 CZM 公式实现本身。

---

## 5. 后续文献调查议程（寻找合适参数）

目标：使 **\(E_{\mathrm{coat}}\cdot|\varepsilon_0|\) 与 \(\sigma_{\max}\)、\(G_c\)、\(K_n\)** 落在同一可起裂、可收敛的物理窗口，并尽量有文献/实验依据。

### 5.1 极片（coating）有效模量

| 调研问题 | 备注 |
|----------|------|
| 干燥态 / 电解液浸润态 \(E_{\mathrm{coat}}\) 典型范围 | 干态常更高；湿态可降至百 MPa 级 |
| 正极 NMC / LFP vs 负极石墨–硅 涂层 | 与 Jellyroll 化学体系对齐 |
| 是否用厚度加权等效模量代替均一 \(E_{\mathrm{coat}}\) | 与 `compute_effective_coating_modulus` 对照 |
| 泊松比 \(\nu_{\mathrm{coat}}\) 对平面应力本征应变的影响 | 验证中 \(\varepsilon_x=0\) 时用 \((1+\nu)\) |

**验收线索**：文献值应使 \(E|\varepsilon_0|\) 至少与目标 \(\sigma_{\max}\) 同量级，或明确说明界面先于体破坏。

### 5.2 涂层–集流体界面 CZM

| 调研问题 | 备注 |
|----------|------|
| Mode I \(\sigma_{\max}\)、\(G_c\)（剥离 / 划痕 / DCB 等） | 正极箔 vs 负极箔可能差一个数量级 |
| 初始刚度 \(K_n\) 或 \(\delta_0=\sigma_{\max}/K_n\) | 避免 \(\delta_0\) 亚纳米 + \(\delta_c/\delta_0\sim 10^3\) 的病态 |
| 湿态 / 循环后界面退化 | 若 \(\sigma_{\max}\) 降至 ~10 MPa，则与 500 MPa 涂层更匹配 |
| Mode II / 混合模式是否需要独立标定 | 当前 model1 以法向为主 |

**验收线索**：推荐一组 \((\sigma_{\max},G_c,K_n)\) 使 \(\delta_c/\delta_0\sim O(10\text{–}100)\)，且在 \(\Delta\mathrm{soc}\sim 0.5\text{–}1\)、代表性 BC 下 \(\delta/\delta_0\) 可过 1。

### 5.3 化学膨胀系数 \(\Omega\) / \(\beta\)

| 调研问题 | 备注 |
|----------|------|
| 宏观极片体积应变 vs SOC（非颗粒 \(\Omega\)） | 代码用 \(\beta=\Omega/3\) 与颗粒场同源，需确认宏观适用性 |
| PE 负膨胀、NE 正膨胀的幅度 | 决定失配符号与哪一侧界面先张 |
| 是否应对 PE/NE 使用不同 \(\alpha\) | 现统一 \(\alpha_{\mathrm{eff}}\) |

### 5.4 边界条件与工况（文献中的力学设定）

| 调研问题 | 备注 |
|----------|------|
| 卷绕电芯局部“集流体近似刚支”是否常见 | 对应 PCC/NCC 固支假设 |
| 自由膨胀 vs 外壳约束 | 影响张开/压缩 |
| 充放电半循环哪一侧更易脱粘 | 指导 \(\Delta\mathrm{soc}\) 符号与观察界面 |

### 5.5 建议的文献检索关键词（英文）

- `electrode coating Young's modulus electrolyte soaked`
- `current collector delamination cohesive zone` / `peel strength anode cathode`
- `traction-separation law battery electrode interface`
- `volume expansion graphite silicon electrode SOC`
- `jellyroll mechanical stress interlayer debonding`

### 5.6 调研产出物（建议）

1. 参数表：\(E_{\mathrm{coat}},\nu,\sigma_{\max},G_c,K_n\)（或 \(\delta_0,\delta_c\)）+ 文献出处 + 干/湿态。  
2. 自洽检查：计算 \(E|\varepsilon_0|/\sigma_{\max}\)；目标 \(\gtrsim 1\)（可起裂）或明确“界面永不先破坏”。  
3. 回填 `Jellyroll.jl` 后重跑：  
   - `julia --project=. test/unit_czm_bilinear.jl`  
   - `julia --project=. test/unit_czm_eigenstrain.jl`（去掉临时 \(\sigma_{\max}\times 0.1\) retune）

---

## 6. 短期建议（在文献到位前）

1. **文档层面**：将“本征应变脚本通过”定义为两档——  
   - A：装配/弹性解析一致（不要求损伤）；  
   - B：损伤可达（需参数自洽）。  
2. **代码层面**：保留 retune 开关并打印 \(E|\varepsilon_0|\) 与 \(\sigma_{\max}\) 比值，避免静默降强度。  
3. **参数层面**：优先查湿态 \(E_{\mathrm{coat}}\) 与箔–涂层剥离 \(\sigma_{\max}/G_c\)；二者同时改比单独改 BC 有效。

---

## 7. 相关文件

| 路径 | 作用 |
|------|------|
| `test/unit_czm_eigenstrain.jl` | 本征应变 + BC +（当前）CZM retune |
| `test/unit_czm_bilinear.jl` | 位移驱动本构验证（已通过） |
| `src/CzmUnitMesh.jl` | 单元条带网格 |
| `src/parameters/Jellyroll.jl` | \(E_{\mathrm{coat}}\)、CZM 原始参数 |
| `docs/superpowers/specs/2026-07-23-unit-czm-strip-verification-design.md` | 原规格（弹性段 + 固定端） |
| `docs/superpowers/findings/2026-07-22-nondimensionalization-particle-vs-macro.md` | 颗粒/宏观尺度与 \(\sigma_{\mathrm{czm}}\) |

---

## 8. 结论摘要

1. **双线性本构与网格路径已通过**；未通过的是“实际量级参数下本征应变致损伤”自洽性。  
2. **根因**：\(E_{\mathrm{coat}}=500\,\mathrm{MPa}\) 时 \(E|\varepsilon_0|\ll\sigma_{\max}^{\mathrm{原始}}\)，存在应力天花板；只改边界无法越过。  
3. **后续**：通过文献确定匹配的涂层模量与界面 \((\sigma_{\max},G_c,K_n)\)，再去掉脚本临时降强度，重跑单元级验证。
