# 参数与单位换算（无量纲 ↔ SI）

本文档汇总 JuBat 中 SPMe 与二维分布热模型耦合时常用物理量的无量纲与 SI 单位换算规则，便于与理论公式对照与结果后处理。

> 注：无量纲化由 `SetParams.NormaliseParam` 完成。带 “*” 表示无量纲量（例如 T*）。下文“常量/尺度”取自 `param_dim.scale`。

## 常用尺度定义（来自 param_dim.scale）
- 参考温度：`T_ref` [K]
- 长度尺度：`L` [m]（电极堆叠厚度）；颗粒半径尺度：`r0` [m]；`a0 = 1/r0` [1/m]
- 时间尺度：`t0 = 3600` [s]
- 法拉第常数：`F` [C/mol]；气体常数：`R` [J/(mol·K)]
- 电势尺度：`phi = (R T_ref / F)` [V]
- 典型电流：`I_typ` [A]；电池面积：`A_cell` [m²]（由 `param_dim.cell.area`）
- 电流密度尺度：`j_scale = I_typ / (a0 · L · A_cell)` [A/m²]
- 电解液电导率尺度：`kappa = L·I_typ / (phi·A_cell)` [S/m]
- 固相电导率尺度：`sig = L·I_typ / (phi·A_cell)` [S/m]
- 电解液扩散系数尺度：`De = L² / te` [m²/s]（其中 `te = F·ce0·A_cell·L / I_typ`）
- 颗粒扩散尺度：`Ds_p = r0² / ts_p`、`Ds_n = r0² / ts_n` [m²/s]（`ts_* = F·cs_max·A_cell·L / I_typ`）
- 反应速率常数尺度：`k_p = j_scale / cs_max / sqrt(ce0)`、`k_n = j_scale / cs_max / sqrt(ce0)`

## 快速换算表（SI ← 无量纲）

| 物理量 | 符号 | 无量纲变量 | SI 单位 | 换算（SI ← 无量纲） | 备注 |
|---|---|---|---|---|---|
| 温度 | T | T* | K | T = T* · T_ref | `Variables["temperature"]` |
| 时间 | t | t* | s | t = t* · t0 |  |
| 长度（厚度/坐标） | x, t_n, t_p, t_sp | x*, t_n*, … | m | x = x* · L | 电极/隔膜厚度同理 |
| 颗粒半径 | r_s | r_s* | m | r_s = r_s* · r0 |  |
| 比表面积 | a_s | a_s* | 1/m | a_s = a_s* · a0 | a0 = 1/r0 |
| 固相锂浓度 | c_s | c_s* | mol/m³ | c_s = c_s* · cs_max | 表面/平均/体内同理 |
| 电解液浓度 | c_e | c_e* | mol/m³ | c_e = c_e* · ce0 |  |
| 界面电流密度 | j | j* | A/m² | j = j* · j_scale | `Variables["… interfacial current density"]` |
| 交换电流密度 | j0 | j0* | A/m² | j0 = j0* · j_scale | 由 k、c_surf、c_e 计算 |
| 过电位/电势/电压 | η, φ, V | η*, φ*, V* | V | 量 = 量* · phi | `phi = R T_ref / F` |
| OCP 熵变 | dU/dT | (dU/dT)* | V/K | (dU/dT) = (dU/dT)* · phi / T_ref | 代码中以 `dUdT*` 存储 |
| 固相扩散系数 | D_s | D_s* | m²/s | 对 NE: D_s = D_s* · Ds_n；对 PE: D_s = D_s* · Ds_p | `ElectrodeDiffusion` |
| 电解液扩散系数 | D_e | D_e* | m²/s | D_e = D_e* · De | `ElectrolyteDiffusion` |
| 固相电导率 | σ | σ* | S/m | σ = σ* · sig | `ElectrodePotential` |
| 电解液电导率 | κ | κ* | S/m | κ = κ* · kappa | `ElectrolytePotential` |
| 反应速率常数 | k | k* | A·m²·(mol⁻¹·⁵) | 对 NE/PE: k = k* · k_n / k_p | 按电极使用对应尺度 |
| 对流换热系数 | h | h* | W/(m²·K) | h = h* · (phi·I_typ) / (A_cell·T_ref) | 见 `NormaliseParam` |

> 说明：表中“对 NE/PE 使用对应尺度”表示在 `NormaliseParam` 中 NE/PE 分别用 `Ds_n/Ds_p` 与 `k_n/k_p` 归一化。

## 热源（用于 2D 热）与单位
耦合到二维热方程的体热源 q（`variables["heat_source_fields"]`）在理论 SI 形式为：

- 反应热：q_rxn = a_s · j · η  [W/m³]
- 可逆热：q_rev = a_s · j · T · (dU/dT)  [W/m³]
- 欧姆热（P2D 场）：q_ohm,s = σ_eff |∇φ_s|²，q_ohm,e = κ_eff |∇φ_e|²  [W/m³]
- SPMe 近似层平均（均流）：
  - 固相：P_s,layer ≈ I² (t_layer/σ_eff)/3，q_ohm,s = P_s,layer / t_layer
  - 电解液：P_e,NE ≈ I² (t_n/κ_ne)/3，P_e,SP ≈ I² (t_sp/κ_sp)，P_e,PE ≈ I² (t_p/κ_pe)/3；q_ohm,e = P_e,layer / t_layer

在无量纲实现中，热源也为无量纲；如需与 SI 对照，请先按上表将 j、η、T、(dU/dT)、σ、κ 等还原为 SI 后代入公式，q 的单位即为 W/m³。

## 备注与注意
- 代码注释曾标注“热参数未完全无量纲化”，但当前实现对 `h`、`cell.heat_Q`、`T_amb/T0` 等已有缩放；若做全 SI 热计算，建议直接使用 `param_dim` 的热物性重建热模型，并用 SI 公式计算 q。
- `U` 与 `dUdT` 在 `NormaliseParam` 中被重载为无量纲函数：`U* = U/phi`，`(dU/dT)* = (dU/dT)/phi · T_ref`。还原时按表中公式即可。

---

## 本示例（SPMe + Jellyroll 2D 强耦合）专用归一化与差异
对应示例：`example/jellyroll_strong_coupling_debug.jl`

本例启用如下特性：
- 并联分流 + 公共端电压求解：`opt.parallel_solve_V = true`
- 单元级 SPMe 强耦合：`opt.coupling_mode = "strong"`，`opt.per_element_spme = true`
- 集流体播种几何（collector-seeded）层权重：`opt.collector_seeded = true`
- 热模型使用二维分布式与“方案B”无量纲：`opt.thermalmodel = "distributed2D"`，`opt.units_thermal = "nd"`

与常规定义的主要区别：
- 电化学电流不再视为单一均匀支路，而是按热网格元素进行并联分流；通过一个“公共端电压”`Vc*`与各支路（元素）方程联立求解。
- 元素级几何采用“层权重”（layer_weights）从 NE/SP/PE/PCC/NCC 映射到元素，有效参数在元素上做带权平均；串联传输量（如厚向欧姆电阻）按串联规则聚合。
- 二维热方程以“面热源”驱动（厚度方向被折算），使用面积尺度进行热通量/换热系数的无量纲化（方案B）。

### 1) 电化学归一化（本例特有）
记：`A_cell` 为整电芯可参与反应的总面积；`A_e` 为热网格元素 e 的在位面积；`A_e* = A_e / A_cell`。

1. 总电流与并联约束（无量纲）：
   - 处方总电流（示例中 `opt.Current` 返回的量）为 `I* = I / I_typ`。
   - 并联分流满足约束：
     $$\sum_e A_e^*\, I_e^* = I^*$$
     其中 `I_e*` 为元素 e 的无量纲支路电流。

2. 公共端电压与支路关系（无量纲）：
   - 端电压无量纲：`V* = V / phi`；示例中通过 `Vc*` 与每个元素内的 SPMe 关系式联立，满足各支路端电压一致：
     $$ V_c^* = U_n^* - U_p^* + (\eta_n^* - \eta_p^*) + I_e^*\, R_{\text{tot},e}^* $$
     其中 `U_*^* = U_*/phi`，`η_*^* = η_*/phi`。

3. 元素总电阻的无量纲化与构造：
   - SI 串联厚向总电阻：
     $$ R_{\text{tot},e} = \frac{t_{\text{PCC}}}{\sigma_{\text{PCC}} A_e} + \frac{t_n}{\sigma_{s,n}^{\text{eff}} A_e} + \left(\frac{t_n}{\kappa_{\text{ne}} A_e} + \frac{t_{\text{sp}}}{\kappa_{\text{sp}} A_e} + \frac{t_p}{\kappa_{\text{pe}} A_e}\right) + \frac{t_p}{\sigma_{s,p}^{\text{eff}} A_e} + \frac{t_{\text{NCC}}}{\sigma_{\text{NCC}} A_e} $$
   - 无量纲电阻（使得 `Δφ* = I* R*` 成立）：
     $$ R_{\text{tot},e}^* = R_{\text{tot},e}\, \frac{I_{\text{typ}}}{\phi} $$
   - 若某材料/层在元素 e 的占比为 `f_{k,e}`（来自 layer_weights，且 `\sum_k f_{k,e}=1`），则对“并联/体平均”类有效量（如 `a_s`、`j0` 等）做加权：
     $$ X_e^* = \sum_k f_{k,e} X_k^* $$
     对“串联”类（如厚向导电/导离子）用上式的串联公式（基于厚度 `t_k` 与 `A_e`）。

4. 元素级反应量与热源所需的无量纲变量：
   - `j*`、`η*`、`U*`、`(dU/dT)*`、`T*` 均与“通用换算表”一致（见上文）。
   - 本例“单元级 SPMe”直接在元素 e 上计算其本地 `j_e^*`、`η_e^*` 等（温度、参数亦可随 e 变化）。

### 2) 热归一化（本例特有，方案B）
本例二维热模型使用面热源驱动，先由电化学得到各层体热源 `q_{\text{layer},e}`，再按层厚度折算为面热源并在元素上汇总。

1. 面热源与尺度：
   - 选取面热通量尺度：
     $$ s_0 = \frac{\phi\, I_{\text{typ}}}{A_{\text{cell}}} \quad [\,\text{W}/\text{m}^2\,] $$
   - 元素总面热源的无量纲化：
     $$ s_e^* = \frac{s_e}{s_0} = \sum_{\text{layer}} q_{\text{layer},e}^*\, t_{\text{layer}}^* $$
     其中 `q^* = q / q_0`，`t^* = t / L`，`q_0 = s_0 / L = (\phi I_{\text{typ}})/(A_{\text{cell}} L)`。

2. 各层体热源（与“通用”一致，先用无量纲变量计算，再按上式折算到面热源）：
   - 反应热：`q_rxn^* \propto a_s^*\, j^*\, \eta^*`
   - 可逆热：`q_rev^* \propto a_s^*\, j^*\, T^*\, (dU/dT)^*`
   - 欧姆热（SPMe 近似层平均）：
     - 固相：`q_{ohm,s}^* \propto (I^*)^2 / \sigma_{\text{eff}}^*`
     - 电解液：`q_{ohm,e}^* \propto (I^*)^2 / \kappa_{\text{eff}}^*`
     实际实现可按上文“SPMe 近似层平均（均流）”给出的层功率关系先求 SI，再转回无量纲。

3. 对流换热（边界条件）与方案B对应无量纲：
   - `h* = h / h_0`，其中 `h_0 = (\phi I_{\text{typ}})/(A_{\text{cell}} T_{\text{ref}})`（见上文“快速换算表”），确保能量方程中面换热项与 `s_e^*` 同量纲。

4. SI 重建速查：
   - 面热源：`s_e [W/m^2] = s_e^* · (\phi I_{\text{typ}} / A_{\text{cell}})`
   - 体热源：`q_e [W/m^3] = q_e^* · (\phi I_{\text{typ}} / (A_{\text{cell}} L))`
   - 公共端电压：`V_c [V] = V_c^* · \phi`
   - 元素电流：`I_e [A] = I_e^* · I_{\text{typ}}`

### 3) 小结与实用提示
- 若仅做对比与后处理，建议在电化学侧优先保持无量纲变量（`j*`、`η*` 等），在热侧汇总为 `s_e^*` 后再按需要还原为 SI 面热源。
- “collector-seeded” 仅影响 `f_{k,e}`（层权重）与 `A_e`（几何），归一化尺度与公式不变；当网格裁剪为环带/螺旋时，`\sum_e A_e = A_{\text{cell}}` 成立，使 `\sum_e A_e^* = 1`，并联约束更为直观。
