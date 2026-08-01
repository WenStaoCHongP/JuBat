# SetParams.jl

- **源文件**: `src/SetParams.jl`
- **行数**: 522 行
- **函数/struct 计数**: 10 个 struct（`Electrode`、`Separator`、`CurrentCollector`、`Electrolyte`、`Cell`、`Tab`、`Binder`、`Cohesive`、`Scale`、`Params`）；2 个函数（`ChooseCell`、`NormaliseParam`）
- **职责**: 定义全部有量纲参数 struct，提供 `ChooseCell` 工厂与 `NormaliseParam` 无量纲化入口
- **相关技术文档**: `md/01_参数定义与归一化.md`、`md/10_参数传递与模块架构.md`、`md/15_颗粒与极片模量区分.md`

## 数据结构

> 顶部 L1-L36 为共享字段命名约定文档注释（NE/PE/SP/EL/NCC/PCC/Cell/Tab/BD 域前缀 + theta/thickness/lambda/Ds 等变量含义）。

### `Electrode`（L37-L71，`@with_kw mutable struct`）

单个电极（正极或负极）的物理参数。同时承担「颗粒」与「极片涂层」两个尺度的属性。

| 字段 | 类型 | 默认值 | 含义/单位 | 归一化意义 |
|------|------|--------|-----------|------------|
| `theta_100` | Float64 | 0 | 100% SOC 化学计量比 | passthrough |
| `theta_0` | Float64 | 0 | 0% SOC 化学计量比 | passthrough |
| `thickness` | Float64 | 0 | 电极厚度 [m] | ÷ scale.L |
| `lambda` | Float64 | 0 | 导热率 [W/(m·K)] | ÷ scale.lambda |
| `Ds` | Float64 | 0 | 颗粒扩散系数 [m²/s] | ÷ scale.Ds_p / scale.Ds_n |
| `rho` | Float64 | 0 | 密度 [kg/m³] | ÷ scale.rho |
| `heat_Q` | Float64 | 0 | 比热容 [J/(kg·K)] | × ρ·L³·T_ref/(t0·P_ref) |
| `eps` | Float64 | 0 | 孔隙率 | passthrough |
| `eps_fi` | Float64 | 0 | 填料体积分数 | passthrough |
| `eps_s` | Float64 | 0 | 活性物质体积分数（= 1 − eps − eps_fi，由 `ChooseCell` 计算） | passthrough |
| `brugg` | Float64 | 0 | Bruggeman 指数 | passthrough |
| `k` | Float64 | 0 | 反应速率常数 [A·m^2.5/mol^1.5] | ÷ scale.k_p / scale.k_n |
| `cs_max` | Float64 | 0 | 颗粒最大锂浓度 [mol/m³] | passthrough |
| `cs0` | Float64 | 0 | 初始锂浓度 [mol/m³] | ÷ cs_max |
| `rs` | Float64 | 0 | 颗粒半径 [m] | ÷ scale.r0 |
| `as` | Float64 | 0 | 比表面积 [1/m] | ÷ scale.a0 |
| `sig` | Float64 | 0 | 电导率 [S/m] | ÷ scale.sig |
| `E` | Float64 | 0 | 颗粒弹性模量 [Pa]（用于颗粒扩散应力） | ÷ scale.E_p / scale.E_n |
| `nu` | Float64 | 0 | 颗粒泊松比 | passthrough |
| `alphaT` | Float64 | 0 | 热膨胀系数 [1/K] | × T_ref |
| `Omega` | Float64 | 0 | 偏摩尔体积 [m³/mol] | × cs_max |
| `E_coat` | Float64 | 0 | 极片涂层宏观弹性模量 [Pa] | ÷ scale.E_coat |
| `nu_coat` | Float64 | 0 | 极片涂层泊松比 | passthrough |
| `Eac_D` | Float64 | 0 | 扩散 Arrhenius 活化能 [J/mol] | ÷ (R·T_ref) |
| `Eac_k` | Float64 | 0 | 反应 Arrhenius 活化能 [J/mol] | ÷ (R·T_ref) |
| `alpha` | Float64 | 0 | Alpha 因子（传递系数相关） | passthrough |
| `U` | Function | x->0.0 | 开路电位 OCP [V] | ÷ scale.phi |
| `dUdT` | Function | x->0.0 | OCP 熵变 [V/K] | ÷ phi × T_ref |
| `M_d`, `K_d`, `M_p`, `K_p` | SparseMatrixCSC | spzeros(0,0) | 颗粒质量/刚度矩阵（运行期填） | — |

### `Separator`（L73-L84，`@with_kw mutable struct`）

隔膜参数。

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `thickness` | Float64 | 0 | 厚度 [m] |
| `lambda` | Float64 | 0 | 导热率 [W/(m·K)] |
| `rho` | Float64 | 0 | 密度 [kg/m³] |
| `heat_Q` | Float64 | 0 | 比热容 [J/(kg·K)] |
| `eps`, `eps_fi`, `brugg` | Float64 | 0 | 孔隙率/填料/Bruggeman |
| `E` | Float64 | 0 | 弹性模量 [Pa] |
| `nu` | Float64 | 0 | 泊松比 |
| `alphaT` | Float64 | 0 | 热膨胀系数 [1/K] |

### `CurrentCollector`（L86-L95，`@with_kw mutable struct`）

集流体参数（PCC/NCC 共用）。

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `thickness` | Float64 | 0 | 厚度 [m] |
| `lambda` | Float64 | 0 | 导热率 [W/(m·K)] |
| `rho` | Float64 | 0 | 密度 [kg/m³] |
| `heat_Q` | Float64 | 0 | 比热容 [J/(kg·K)] |
| `sig` | Float64 | 0 | 电导率 [S/m] |
| `E`, `nu`, `alphaT` | Float64 | 0 | 弹性/泊松比/热膨胀 |

### `Electrolyte`（L97-L107，`@with_kw mutable struct`）

电解液参数（多为函数形式，依赖浓度与温度）。

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `De` | Function | (x,y=0)->0 | 扩散系数 [m²/s] |
| `kappa` | Function | (x,y=0)->0 | 离子电导率 [S/m] |
| `dlnf_dlnc` | Function | (x,y=0)->0 | 活度系数导数 |
| `rho`, `heat_Q` | Float64 | 0 | 密度/比热容 |
| `tplus` | Float64 | 0 | 迁移数 |
| `ce0` | Float64 | 0 | 初始锂浓度 [mol/m³] |
| `Eac_D`, `Eac_k` | Float64 | 0 | Arrhenius 活化能 |

### `Cell`（L109-L139，`@with_kw mutable struct`）

电池单体级几何与热学参数。

| 字段 | 类型 | 默认值 | 含义/单位 |
|------|------|--------|-----------|
| `length`, `width`, `wrapper` | Float64 | 0 | 几何尺寸 [m] |
| `I1C` | Float64 | 0 | 1C 电流 [A] |
| `no_layers` | Int32 | 0 | 卷绕层数 |
| `capacity` | Float64 | 0 | 容量 [Ah] |
| `cooling_surface`, `area` | Float64 | 0 | 冷却面积 / 截面积 [m²] |
| `v_h`, `v_l` | Float64 | 0 | 上下限电压 [V] |
| `volume` | Float64 | 0 | 体积 [m³] |
| `rho`, `mass` | Float64 | 0 | 密度/质量 |
| `alphaT`, `heat_Q`, `h` | Float64 | 0 | 热膨胀/比热容/换热系数 |
| `T0` | Float64 | 298 | 初始温度 [K] |
| `T_amb` | Float64 | 0 | 环境温度 [K] |
| `Rin`, `Rout` | Float64 | 0 | Jellyroll 内/外半径 [m] |
| `layer` | Float64 | 0 | 单层叠合厚度 [m] |
| `lambda_r`, `lambda_t` | Float64 | 0 | 径向/切向等效导热率 [W/(m·K)] |
| `Nr_th`, `Nθ_th` | Int | 0 | Jellyroll 顶视热网格划分 |
| `n_windings` | Int | 0 | 卷绕数 |

### `Tab`（L141-L148，`@with_kw mutable struct`）

极耳参数。

| 字段 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `length`, `width`, `area` | Float64 | 0 | 几何 [m / m / m²] |
| `h` | Float64 | 0 | 换热系数 |
| `theta_pos`, `theta_neg` | Vector{Float64} | Float64[] | 极耳化学计量比 |

### `Binder`（L152-L154，`@with_kw mutable struct`）

粘结剂参数（仅密度）。

| 字段 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `rho` | Float64 | 0 | 密度 [kg/m³] |

### `Cohesive`（L156-L193，`@with_kw mutable struct`）

内聚力模型参数，按 PE-PCC 界面与 NE-NCC 界面分组（法向 `n`、切向 `t` 各 5 参数：σ_max、K、δ_0、G_c、δ_c）。

| 字段组 | 字段 | 含义/单位 |
|--------|------|-----------|
| PE-PCC 法向 | `σ_max_pe_pcc`, `K_n_pe_pcc`, `δ_0_pe_pcc`, `G_c_pe_pcc`, `δ_c_pe_pcc` | 最大牵引 [Pa]、刚度 [Pa/m]、损伤起始位移 [m]、断裂能 [J/m²]、临界位移 [m] |
| PE-PCC 切向 | `τ_max_pe_pcc`, `K_t_pe_pcc`, `δ_0_pe_pcc_t`, `G_c_pe_pcc_t`, `δ_c_pe_pcc_t` | 同上，Mode II |
| NE-NCC 法向 | `σ_max_ne_ncc`, `K_n_ne_ncc`, `δ_0_ne_ncc`, `G_c_ne_ncc`, `δ_c_ne_ncc` | 同上 |
| NE-NCC 切向 | `τ_max_ne_ncc`, `K_t_ne_ncc`, `δ_0_ne_ncc_t`, `G_c_ne_ncc_t`, `δ_c_ne_ncc_t` | 同上 |
| 共用 | `czm_model`, `eta` | 模型选择 / BK 指数 |
| 界面热阻 | `h_c0`=1e7, `k_air`=0.026, `lambda_m`=70e-9, `beta`=1.0, `threshold`=70e-9 | 间隙导热模型参数 |
| 粘性 | `tau_visc`=0.0 | 物理松弛时间 [s] |

### `Scale`（L194-L233，`@with_kw mutable struct`）

归一化参考尺度集合（由 `ChooseCell` 填充）。

| 字段组 | 字段 | 含义 |
|--------|------|------|
| 长度 | `L`=1e-6, `r0`=1e-6, `a0`=1/r0 | 长度/颗粒半径/比表面积尺度 |
| 时间 | `t0`=3600 | 统一时间尺度 [s] |
| 热力学 | `T_ref`=298, `F`, `R` | 参考温度 / 法拉第常数 / 气体常数 |
| 电化学 | `j`, `Ds_p`, `Ds_n`, `ts_p`, `ts_n`, `te`, `De`, `phi`, `sig`, `kappa` | 电流密度/扩散/时间/电位/电导率尺度 |
| 浓度 | `cp_max`, `cn_max`, `ce` | 最大浓度尺度 |
| 反应 | `k_p`, `k_n` | 反应速率尺度 |
| 电流 | `I_typ`, `R_cell` | 典型电流 / 电池电阻 |
| 力学（颗粒） | `E_n`, `E_p` | 颗粒扩散应力尺度 |
| 力学（极片） | `E_coat` | 极片宏观模量参考 [Pa] |
| 热 | `rho`, `P_ref`, `lambda`, `q`, `h` | 统一能量尺度热参数 |
| CZM | `σ_czm`, `δ_czm`, `G_czm`, `K_czm` | CZM 归一化锚点（PE-PCC 界面） |

### `Params`（L235-L247，`@with_kw mutable struct`）

顶层参数容器，聚合各子 struct。

| 字段 | 类型 | 含义 |
|------|------|------|
| `PE`, `NE` | Electrode | 正/负极 |
| `EL` | Electrolyte | 电解液 |
| `SP` | Separator | 隔膜 |
| `cell` | Cell | 单体 |
| `PCC`, `NCC` | CurrentCollector | 集 流 体 |
| `tab` | Tab | 极耳 |
| `binder` | Binder | 粘结剂 |
| `scale` | Scale | 归一化尺度 |
| `cohesive` | Cohesive | CZM 参数（默认空 Cohesive()） |

## 函数清单

### `ChooseCell`（L249-L348）

```julia
function ChooseCell(CellType::String="LG M50")
```

**职责**: 按型号名称加载参数文件并计算 `scale` 中全部参考尺度（电化学 + 热学 + 力学 + CZM）。

**关键逻辑**:
- L257-L268：根据 `CellType`（"LG M50" / "Northrop" / "Enertech" / "Jellyroll" / "Ring"）`include` 对应参数文件
- L269-L272：补全 `eps_s = 1 - eps - eps_fi` 与 `as = 3·eps_s/rs`
- L273-L291：若 `cell.rho` / `cell.heat_Q` 近 0（阈值 1e-8），用厚度加权回退值
- L293-L310：填充电化学 `scale` 尺度（I_typ、L、j、ts_p/n、te、Ds_p/n、De、phi、sig、kappa、cp/cn_max、E_n/p）
- L312-L324：极片模量尺度 `E_coat` —— 缺失时 `@warn` 并保持 0；否则按全叠合厚度加权
- L325-L332：热学尺度 `P_ref = phi·I_typ`，`lambda`、`q`、`h`（Biot 数）
- L333-L346：CZM 锚点 —— 缺失 cohesive 数据时 `@warn`；`δ_czm` 锚点取 `2·G_c/σ_max`，缺失时回退 `scale.L`

**跨文件依赖**: `parameters/LGM50.jl`、`parameters/Northrop.jl`、`parameters/Enertech.jl`、`parameters/Jellyroll.jl`、`parameters/Ring.jl`

### `NormaliseParam`（L350-L521）

```julia
function NormaliseParam(param_dim::Params)
```

**职责**: 接收有量纲 `Params`，返回深拷贝后的无量纲 `Params`（`scale` 子结构保持有量纲）。

**关键逻辑**:
- L357：`deepcopy(param_dim)` 后逐字段归一化
- L359-L385：正极字段归一化（theta passthrough、cs0÷cs_max、thickness÷L、Ds÷Ds_p、k÷k_p、sig÷sig、E÷E_p、E_coat÷E_coat、U÷phi、dUdT÷phi×T_ref、lambda÷lambda、heat_Q 用统一能量公式）
- L386-L411：负极对称（÷Ds_n、÷E_n、÷k_n）
- L413-L428：隔膜（含 `scale.E_coat > 0` 守卫）
- L430-L455：PCC/NCC（含 `scale.E_coat > 0` 守卫）
- L457-L461：电解液（用 `Base.invokelatest` 包装避免 world-age 问题）
- L463-L484：cell 与 tab 几何归一化
- L486-L518：cohesive per-interface 归一化（÷σ_czm / ÷δ_czm / ÷G_czm / ÷K_czm）+ 界面热阻参数（h_c0·L/lambda 等）

**跨文件依赖**: `Params`、`Scale`（同文件）

## 省略项

- L36 注释行 `# param_dim.Tab.width = ...`、L149-L150 注释行：已注释的代码片段，不计入条目

### [DEBUG]

无

> 上述 `@warn`（L313、L334）是参数完整性校验的结构化警告，按规范不算 [DEBUG]。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L339 | `# cohesive 数据缺失时回退到 scale.L（保持非 CZM 参数集的旧行为，Λ = 1）。if ... > 0 ... else scale.δ_czm = scale.L end`（L339-L344） | 注释含「回退/旧行为」；非 CZM 参数集的兜底，CZM 启用时由前置 `@warn` 拦截，风险低 |
| L312 | `if PE.E_coat == 0 \|\| NE.E_coat == 0 ... @warn ... scale.E_coat = 0`（L312-L324） | 缺失字段时 `scale.E_coat` 保持 0；下游 `NormaliseParam` 中用 `E_coat > 0` 守卫跳过力学归一化。属「兜底」行为，但伴随 `@warn`，风险低 |
| L333 | `if σ_max_pe_pcc == 0 \|\| G_c_pe_pcc == 0 @warn ... end`（L333-L335） | 同上，CZM 锚点缺失的兜底（不直接吞错，伴随 `@warn`） |

> 上述均为带 `@warn` 的显式回退，非静默吞错，按规范风险列说明但不算严重 PLACEHOLDER。

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L273 | `if abs(rho) < 1e-8 ... 厚度加权公式 ... else ... end` + `if abs(heat_Q) < 1e-8 ... 厚度加权公式 ... else ... end`（L273-L291） | 连续两段相似守卫 + 长公式（heat_Q 表达式 > 100 字符），可抽取 `layer_weighted_average(fields, layers)` 辅助函数 |
| L359 | 正极/负极归一化块（各 ~25 行）高度对称（L359-L411） | 重复结构可参数化为 `_normalise_electrode!(param, param_dim, which::Symbol)` 减半代码 |
