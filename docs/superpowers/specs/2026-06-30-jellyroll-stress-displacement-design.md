# Jellyroll 电池二维应力/位移场后处理脚本设计

- 日期：2026-06-30
- 类型：新增示例脚本（零侵入）
- 关联文件：`example/testexample.jl`、`src/Mechanical.jl`
- 关联函数：`JuBat.thermal_diffusion_stress_2D(case, variables)`

---

## 1. 背景与动机

`src/Mechanical.jl:165` 中已实现 `thermal_diffusion_stress_2D`：基于二维 Q4 平面应力 FEM，以
"热膨胀应变 `α_eff · ΔT`" 与 "锂化膨胀应变 `β_n · Δsoc_n + β_p · Δsoc_p`" 作为初始应变，
求解并输出单元应力 `σ_xx / σ_yy / σ_xy / σ_vm [Pa]` 与节点位移 `U_x / U_y [m]`。

但该函数：

1. 未被 `Solve` 流程自动调用（`mechanicalmodel = "full"` 仅触发颗粒尺度 `Calstressdisp`，不影响宏观 2D 场）；
2. 现有 `example/mechanical_example.jl` 是 LG M50 颗粒尺度验证，并非 Jellyroll 2D 宏观场脚本。

因此需要一个面向 **Jellyroll 果冻卷电池** 的示例脚本：单次放电后，在用户指定的时间节点
输出宏观应力场与位移场云图，**不启用 CZM**，仅传统固体力学部分。

## 2. 目标与非目标

**目标：**

- 新建 `example/jellyroll_stress_displacement.jl`，单文件、零侵入、不修改 `src/`。
- 复用 `testexample.jl` 的求解器配置（Jellyroll、1C、3600 s、nθ=360、per_element_spme、热耦合开启）。
- 与 testexample 的差异：`opt.czm_enabled = false`（testexample 中为 true）；
  `mechanicalmodel` 两脚本均为 `"none"`（testexample 已是颗粒应力关闭状态）。
- 对用户在脚本顶部显式声明的若干物理时间节点，调用 `thermal_diffusion_stress_2D` 计算场，
  绘制"所有节点拼成一张大图"。
- 控制台打印每个节点的 `σ_vm` 峰值、`|U|` 峰值，作为数量级自检。

**非目标：**

- 不修改核心求解器（`Solve`、`Variables.jl`、`PostProcessing.jl` 等）。
- 不导出 CSV/VTK 等场数据文件（仅 PNG）。
- 不编写正式单元测试（脚本属于 `example/`）。
- 不实现颗粒尺度应力（`Calstressdisp`）相关绘图。

## 3. 应力源

复用 `thermal_diffusion_stress_2D` 现有行为：**热 + 扩散耦合**作为初始应变共同驱动应力场。
不在脚本侧屏蔽任何一项；不引入纯热应力 / 纯扩散应力分支。

## 4. 求解器配置

```julia
param_dim = JuBat.ChooseCell("Jellyroll")
param_dim.cell.v_l = 2.5
param_dim.cell.v_h = 4.2

opt = JuBat.Option()
opt.Current       = x -> 5.0          # 1C
opt.model         = "SPMe"
opt.Nn = 10; opt.Ns = 5; opt.Np = 10
opt.Nrn = 10; opt.Nrp = 10
opt.gsorder      = 2
opt.dimension    = 1
opt.mechanicalmodel = "none"          # 不启用颗粒尺度应力耦合

opt.time = [0.0, 3600.0]
opt.dt   = [0.5, 10]
opt.dtType    = "auto"
opt.jacobi    = "update"
opt.solveType = "Crank-Nicolson"

opt.thermal_enabled = true
opt.thermalmodel    = "distributed2D"
opt.thermal_dim     = "2D"
opt.cool_method     = "surface"
opt.per_element_spme = true

opt.czm_enabled = false                # 关键差异：关闭 CZM
```

网格：`mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=360, gsorder=2)`，
`case = JuBat.setup_thermal2D_mesh(case, mesh_data)`。
因 `czm_enabled = false`，不创建 `case.czm_mesh`。

> **E_coat 前置保障**：`thermal_diffusion_stress_2D` 入口有 `@assert PE.E_coat > 0 && NE.E_coat > 0`
> （`Mechanical.jl:167`），Jellyroll 参数集 `PE.E_coat = 4.6e8`、`NE.E_coat = 4.8e8` 满足，无需脚本侧处理。

## 5. 时间节点选择

脚本顶部声明：

```julia
const PLOT_TIMES_S = [600, 1800, 3600]   # 用户可改
```

选择算法：

1. 取 `t_s = result["time [s]"]`。
2. 对每个 `t_target` 计算 `ti = argmin(abs.(t_s .- t_target))`。
3. 对超出 `[t_s[1], t_s[end]]` 的目标，`@warn` 并钳到端点。
4. 去重 + 按 `ti` 升序排序，保证子图列顺序与时间递增一致。

## 6. 单节点场重建（核心算法）

对每个时间索引 `ti`：

```julia
variables_ti = Dict{String, Union{Array{Float64},Float64}}(
    "T_nodes"                 => result["thermal2D temperature at nodes [K]"][:, ti],
    "thermal2D element soc_n" => result["thermal2D element soc_n"][:, ti],
    "thermal2D element soc_p" => result["thermal2D element soc_p"][:, ti],
)
vars_out = JuBat.thermal_diffusion_stress_2D(case, variables_ti)
```

从 `vars_out` 读取：

- `diffusion stress xx / yy / xy / vonMises`（单元常数，`ne` 长度）
- `displacement x / y`（节点值，`nlen` 长度）

### 6.1 单位约定（已核实）

> 本函数内部 `param = case.param`（`NormaliseParam` 后的归一化空间）；以下 `param.X` 均指归一化后
> 数值，`param_dim.X` 仅出现于 `case.param_dim`。

- **温度**：`thermal_diffusion_stress_2D` 内部 `T0 = param.cell.T0`（K）、`Tref = param.scale.T_ref`，
  直接做 `dT = T_elem - T0`，因此期望**有量纲 K**。`result["thermal2D temperature at nodes [K]"]`
  在 PostProcessing 中 `*= scale.T_ref`，正好是 K。一致。
- **SOC（已确认一致）**：`thermal_diffusion_stress_2D` 内部 `soc_ref_n = param.NE.cs0`，注意是
  `case.param` 即 `NormaliseParam` 后的归一化值 `= param_dim.NE.cs0 / param_dim.NE.cs_max`
  （见 `SetParams.jl` 中 `NormaliseParam` 对 `cs0` 的归一化）。`result["thermal2D element soc_n"]`
  来源于 `CallModel.jl` 的 `mean(csn_data)`，同样是 cs_max 归一化值（求解器全程在归一化空间，
  PostProcessing 对 `soc_n/p` 字段未做任何缩放）。**两者同量纲，脚本侧无需任何缩放，不要乘除
  `scale.cn_max`**。实现时仅需打印 `param.NE.cs0` 与字段范围做一次数量级 sanity check
  （预期都在 [0, 1] 区间）。
- **输出**：函数内部已对位移乘 `param.scale.L`、应力使用 `param.scale` 还原后写出 `[Pa]` 与 `[m]`，
  脚本侧无需再次缩放。

## 7. 绘图

- 依赖：`Plots.jl`（不引入新依赖）。
- 自写 Q4 cell plotting：用 `Plots.Shape` 把每个单元四角连成多边形填色；单元常数场（应力）
  以四节点平均值或单元中心色块呈现；节点场（位移）先插值到单元中心或用节点填色四边形。
- **布局**：一张大图，行=场分量，列=时间节点（共 5 行 × N 列）：
  - 行 1：`σ_xx [MPa]`
  - 行 2：`σ_yy [MPa]`
  - 行 3：`σ_xy [MPa]`
  - 行 4：`σ_vm [MPa]`
  - 行 5：`|U| [µm]`，并叠加放大变形网格（放大系数 `DEF_SCALE`）
- `DEF_SCALE` 默认自适应：`DEF_SCALE = 0.05 * L_ref / max(|U|_global)`，使最大变形 ≈ 5% 特征长度；
  也可在脚本顶部手动覆写为常数。
- 每子图标题：`t = 600 s` 等；轴标签统一 `x [m]`、`y [m]`。
- 输出：`output/jellyroll_stress_displacement.png`，`size=(400*N, 2000)` 左右，DPI 默认。

## 8. 错误处理

- 量纲断言：`@assert length(T_nodes) == mesh_th.nlen`、`length(soc_n) == ne`。
- 若 `haskey(result, "thermal2D temperature at nodes [K]") == false` → `error(...)` 提示用户检查
  `opt.thermalmodel == "distributed2D"`。
- 时间节点超出范围 → `@warn` 钳到端点（见 §5）。
- `thermal_diffusion_stress_2D` 内部已有 `K_mech \ F_mech` 的 try/catch，失败时返回零位移并 `@warn`；
  脚本无需额外处理，但会在数量级自检时发现并提示。

## 9. 数量级自检（非正式）

首个时间节点跑通后，控制台 `@printf`：

- `σ_vm_max [MPa]`（预期 MPa–GPa 量级，与极片涂层模量 `E_coat ~ 5e8 Pa` 协调）
- `|U|_max [µm]`（预期 µm–mm 量级）
- 对应的单元索引 / 节点索引
- 顺带打印 `param.NE.cs0` 与 `result["thermal2D element soc_n"][:, ti]` 的 min/max，作为 §6.1
  归一化一致性的 sanity check（预期都在 [0, 1] 区间）。

若量级显著异常（应力 < 1 kPa 或 > 1 TPa），排查网格/参数；若 SOC 范围不在 [0, 1]，再回头
排查归一化路径。

## 10. 输出文件

| 文件 | 说明 |
|------|------|
| `example/jellyroll_stress_displacement.jl` | 新建脚本 |
| `output/jellyroll_stress_displacement.png` | 大图（所有时间节点 × 所有场分量） |

## 11. 不做的事（明确边界）

- 不修改任何 `src/` 文件。
- 不增加新的 `Option` 字段。
- 不修改 `thermal_diffusion_stress_2D` 的签名或行为。
- 不导出场数据为 CSV/VTK。
- 不编写正式单元测试。
- 不在 `Solve` 内自动调用宏观力学函数。
