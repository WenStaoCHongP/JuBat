# AGENTS.md

本文档为 Codex 提供 JuBat 项目中 **Jellyroll 果冻卷电池 SPMe-二维分布式热-CZM 耦合模型**的开发指导。

---

## 1. 项目概述

JuBat 是基于 Julia 的电池建模框架，采用二阶有限元方法 (FEM) 求解方程。本项目专注于 **Jellyroll（果冻卷）型号电池**的多物理场耦合仿真：

| 物理场 | 模型 | 技术文档 |
|--------|------|----------|
| 电化学 | SPMe (Single Particle Model with electrolyte) | `md/04_电化学模型_SPMe.md` |
| 热学 | 二维分布式热模型 (各向异性导热) | `md/05_热模型_二维分布式.md` |
| 力学 | 内聚力模型 (CZM) | `md/06_内聚力模型_CZM.md` |
| 耦合 | 界面热阻模型 | `md/07_界面热阻模型.md` |

**参考文献**: W. Ai, Y. Liu, Improving the convergence rate of Newman's battery model using 2nd order finite element method, J. Energy Storage. 67 (2023) 107512.

---

## 2. 技术文档索引

详细技术文档位于 `md/` 目录，按层次组织：

### 第一层：参数与基础 (01-03)
| 编号 | 文档 | 内容 |
|------|------|------|
| 01 | `01_参数定义与归一化.md` | 物理参数、电化学/热学/机械归一化、Jellyroll专用参数 |
| 02 | `02_几何与网格.md` | 阿基米德螺旋线、collector-seeded网格、COH2D4单元、CohesiveMesh结构 |
| 03 | `03_边界条件.md` | 侧面/极耳冷却、界面热阻边界条件 |

### 第二层：模型实现 (04-07)
| 编号 | 文档 | 内容 |
|------|------|------|
| 04 | `04_电化学模型_SPMe.md` | 颗粒扩散、电解液守恒、Butler-Volmer、机械耦合 |
| 05 | `05_热模型_二维分布式.md` | 能量方程、分层热源、各向异性导热、极坐标FVM |
| 06 | `06_内聚力模型_CZM.md` | 双线性牵引-分离律、CZMResult结构、热-化学载荷 |
| 07 | `07_界面热阻模型.md` | 间隙导热系数、损伤耦合、使用建议 |

### 第三层：算法与求解 (08-10)
| 编号 | 文档 | 内容 |
|------|------|------|
| 08 | `08_逐单元算法.md` | 多SPMe并行架构、状态向量设计、分层热源计算 |
| 09 | `09_分流求解器.md` | Newton-Raphson分流、截止电压检测、CZM失效处理 |
| 10 | `10_参数传递与模块架构.md` | Case/variables结构、CycleSolver、耦合数据流 |

### 第四层：验证方案 (11-13)
| 编号 | 文档 | 内容 |
|------|------|------|
| 11 | `11_电化学验证方案.md` | SPMe验证、验收标准 |
| 12 | `12_热模型验证方案.md` | 圆环精确解、FVM验证 |
| 13 | `13_耦合验证方案.md` | 电-热-CZM全耦合验证 |

---

## 3. 快速开始

### 3.1 模块加载

```julia
include("src/JuBat.jl")
using .JuBat
```

### 3.2 基础仿真流程

```julia
# 1. 选择 Jellyroll 电池参数
param_dim = JuBat.ChooseCell("Jellyroll")

# 2. 配置选项
opt = JuBat.Option()
opt.model = "SPMe"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.per_element_spme = true
opt.czm_enabled = true
opt.mechanicalmodel = "full"

# 3. 创建案例
case = JuBat.SetCase(param_dim, opt)

# 4. 创建热网格 (Jellyroll 专用)
mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)

# 5. 求解
result = JuBat.Solve(case)
```

### 3.3 循环仿真

```julia
cycle_opt = JuBat.CycleOption(
    n_cycles = 50,
    I_charge = 5.0,
    I_discharge = 5.0,
    t_charge = 3600,
    t_discharge = 3600,
    V_upper = 4.2,
    V_lower = 2.5,
    SOC_init = 0.05
)
result = JuBat.solve_cycling(case, cycle_opt)
```

---

## 4. 核心架构

### 4.1 多 SPMe 并行架构

当 `opt.per_element_spme = true` 时，每个热单元拥有独立的 SPMe 模型：

```
状态向量结构:
yt_global = [yt_chem[1]; yt_chem[2]; ...; yt_chem[ne]; T_nodes]

每个单元状态:
yt_chem[e] = [cn_surf[1:Nrn]; cp_surf[1:Nrp]; ce[1:Nel]]

矩阵结构:
M_global = blockdiag(M_elems..., MT)
```

详见 `md/08_逐单元算法.md`。

### 4.2 Jellyroll 几何

采用阿基米德螺旋线描述：`r(θ) = a + bθ`

- 内半径 `a = R_in`
- 螺旋增长率 `b = t_repeat / (2π)`
- 层序: PE → PCC → PE → SP → NE → NCC → NE → SP

详见 `md/02_几何与网格.md`。

### 4.3 耦合数据流

```
SPMe (电化学)
    │ 传递内热源 q_e (反应热+可逆热+欧姆热)
    ▼
Thermal2D (热模型)
    │ 传递单元均温 T_e
    ▼
SPMe (温度影响反应速率、电导率)
    │
    │ 传递扩散应力 + 热应力
    ▼
CZM (内聚力模型)
    │ 传递损伤状态 D
    ▼
┌─────────────────────────────┐
│ 分流求解器 (影响电流分布)    │
│ 界面热阻模型 (影响层间导热)  │
└─────────────────────────────┘
```

详见 `md/10_参数传递与模块架构.md`。

---

## 5. 关键配置选项

### 5.1 Option 结构 (`src/Option.jl`)

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `model` | "SPM" | 模型类型，Jellyroll 用 "SPMe" |
| `Nrp`, `Nrn` | 10 | 正/负极颗粒网格点数 |
| `Np`, `Ns`, `Nn` | 10 | 电极/隔膜网格点数 |
| `gsorder` | 2 | 高斯积分阶数 |
| `Current` | x->5.0 | 电流函数 I(t) |
| `time` | [0, 3600] | 仿真时间范围 |
| `dt` | [1e-6, 10] | 自适应时间步范围 |
| `solveType` | "Crank-Nicolson" | 时间离散格式 |

### 5.2 热学选项

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `thermal_enabled` | false | 启用热耦合 |
| `thermalmodel` | "none" | 热模型类型，用 "distributed2D" 或 "ring2D_polar" |
| `cool_method` | "surface" | 冷却方式: "surface" 或 "tab" |
| `per_element_spme` | false | 启用逐单元 SPMe |

### 5.3 机械/CZM 选项

| 字段 | 默认值 | 说明 |
|------|--------|------|
| `mechanicalmodel` | "none" | 机械模型，用 "full" |
| `czm_enabled` | false | 启用 CZM |
| `czm_update_interval` | 1 | 损伤更新间隔 |
| `czm_iter_method` | "basic" | CZM 迭代方法 |
| `czm_area_loss_enabled` | false | 启用渐进式有效面积损失（D > threshold 时缩减面积） |
| `czm_area_loss_threshold` | 0.83 | 面积开始缩减的 D 阈值 |

---

## 6. 输出变量

### 6.1 基础变量

| 键 | 类型 | 说明 |
|----|------|------|
| `time [s]` | Vector | 时间序列 |
| `cell voltage [V]` | Vector | 端电压 |
| `temperature` | Float64 | 体积平均温度 (K) |

### 6.2 热学变量

| 键 | 说明 |
|----|------|
| `thermal2D temperature [K]` | 节点温度场 |
| `thermal2D element current` | 各单元电流 |
| `thermal2D element area` | 各单元面积 |
| `heat_source_fields` | 各单元热源 |

### 6.3 分层热源 (W/m³)

| 键 | 说明 |
|----|------|
| `thermal2D Q_rxn_NE [W/m3]` | 负极反应热 |
| `thermal2D Q_rev_NE [W/m3]` | 负极可逆热 |
| `thermal2D Q_SP [W/m3]` | 隔膜欧姆热 |
| `thermal2D Q_rxn_PE [W/m3]` | 正极反应热 |
| `thermal2D Q_PCC/NCC [W/m3]` | 集流体欧姆热 |

### 6.4 CZM 变量

| 键 | 说明 |
|----|------|
| `czm D_max` | 最大损伤值 |
| `czm D_mean` | 平均损伤值 |
| `czm n_fractured` | 完全断裂数量 |
| `czm δ_max_n [m]` / `czm δ_mean_n [m]` | 最大/平均法向分离 |

注意：`soh` 不是单次求解 result 字典的键，而是 `CyclingResult.soh::Vector{Float64}` 字段（`cycle_summary.csv` 中同名列）。

---

## 7. 代码文件组织

### 7.1 核心求解

| 文件 | 功能 |
|------|------|
| `src/Solve.jl` | 主求解器（时间推进、纯热分支） |
| `src/CallModel.jl` | 模型调度：CallModel、CallModel_MultiSPMe |
| `src/SPMe.jl` | SPMe 模型、SPMe_element |
| `src/CycleSolver.jl` | 循环求解器、PhaseResult/CycleResult |

### 7.2 热模型

| 文件 | 功能 |
|------|------|
| `src/ThermalDistributed.jl` | 二维分布式热模型、边界条件 |
| `src/ThermalPolar2D.jl` | 极坐标 FVM 求解器 |

### 7.3 CZM

| 文件 | 功能 |
|------|------|
| `src/CzmMesh.jl` | CZM 界面单元拓扑与 `create_czm_mesh` |
| `src/czm.jl` | CZM 损伤状态、材料映射与系统装配 |
| `src/CzmBC.jl` | CZM 边界节点识别与边界条件施加 |
| `src/CzmSolve.jl` | Newton-Raphson CZM 求解 |
| `src/Mechanical.jl` | 应力计算 |

### 7.4 Jellyroll

| 文件 | 功能 |
|------|------|
| `src/Jellyrollmodel.jl` | 螺旋几何、collector-seeded网格 |
| `src/parameters/Jellyroll.jl` | Jellyroll 参数集 |

### 7.5 并行求解

| 文件 | 功能 |
|------|------|
| `src/Parallelsolution.jl` | 分流求解器、solve_branch_currents_newton |

---

## 8. 调试与验证

### 8.1 调试选项

```julia
opt.debug_coupling = true
opt.debug_log_path = "output/debug.log"
```

### 8.2 验证脚本

热模型验证脚本位于 `example/热模块验证/`：

| 文件 | 说明 | 状态 |
|------|------|------|
| `thermal_verify.jl` | 圆环精确解验证 | ✓ 已更新 |
| `thermal_error_source_analysis.jl` | 热误差来源分析 | ✓ 已修正几何尺度 (2026-03-15) |
| `thermal_equivalent_lumped_compare.jl` | 2D等效集总量对比 | ✓ 已修正几何尺度 (2026-03-15) |

**重要说明**：
- 验证脚本已更新以使用统一归一化方案
- 网格坐标已归一化（`x* = x / L`），从网格积分得到的 `A_elem` 是无量纲面积
- 计算物理功率时必须乘以 `scale.L^2` 转换为物理面积
- 详见 `docs/plans/thermal_verify_normalization_findings.md`

### 8.3 常见问题

| 问题 | 解决方案 |
|------|----------|
| 收敛困难 | 减小 `dt_max`，检查参数归一化 |
| 温度异常 | 检查热源计算、边界条件 |
| CZM 不收敛 | 降低 `czm_update_interval`，检查刚度参数 |
| 截止电压误触发 | 检查分流求解器预因子 |
| 功率计算量级异常 | 检查是否正确使用 `scale.L^2` 转换无量纲面积 |
| 热源归一化错误 | 确保使用 `scale.q = P_ref / L³` 进行归一化 |

---

## 9. 开发约定

### 9.1 无量纲化（统一归一化方案）

> **重要更新 (2026-03-12)**: JuBat 采用**统一时间尺度**和**统一能量尺度**策略。

内部计算使用无量纲变量，参数在 `NormaliseParam` 中归一化：

**统一参考尺度**：
- **时间尺度**: `t_0 = 3600` s（电化学与热模型统一）
- **功率尺度**: `P_ref = φ_ref × I_typ`（统一能量尺度）
- **长度尺度**: `L = t_PE + t_SP + t_NE`（电极堆叠厚度）
- **温度尺度**: `T_ref = 298` K
- **电流尺度**: `I_typ`（1C 电流）

**热参数归一化**（基于能量守恒原则）：
- **密度**: `ρ* = ρ / ρ_ref`
- **比热容**: `c* = c × ρ_ref × L³ × T_ref / (t_0 × P_ref)`
- **体积热容**: `(ρc)* = ρc × L³ × T_ref / (t_0 × P_ref)`
- **导热率**: `k* = k × L × T_ref / P_ref`
- **热源**: `q* = q × L³ / P_ref`（即 `q* = q / scale.q`）

**关键特征**：
- 电化学和热模型使用相同时间尺度，无需额外时间缩放因子
- 所有热参数归一化基于统一的电功率参考 `P_ref`
- 网格坐标已归一化（`x* = x / L`），因此网格积分得到的面积是无量纲的

**结果还原**：
结果需通过 `param_dim.scale` 中的参考值还原为物理单位。

详见 `md/01_参数定义与归一化.md`。

### 9.2 状态管理

- 循环仿真中，`final_state` 作为下一相位的 `initial_state`
- 损伤状态跨周期累积
- 可通过 `reset_T_each_cycle = true` 重置温度

### 9.3 网格分辨率

- `nθ` 同时控制热网格与分层力学/CZM 网格的周向分辨率 (典型值 80-200)；力学网格必须直接继承热网格的实际角节点，不得另设 `nθ_czm` 或独立角度数组
- CZM 子网格用 `czm_enabled=true` 启用；`thermal_elem_map` 必须由共享角段的父子拓扑直接构造
- `CzmSubmesh.phi_pairs` 保存 `(outer_node_at_θ, inner_node_at_θ+2π)` 力学匝间配对，供后续接触/约束装配使用
- cohesive 计数统一为：2 种本构/材料类型（`:PE_PCC`、`:NE_NCC`），每个 8 层卷绕重复单元 4 个真实箔–涂层面，完整离散单元数 `4 * (length(theta)-1)`；`phi_pairs` 不参与该计数
- `gsorder` 控制积分精度 (2-3)
- 平衡精度与计算成本

### 9.4 弹性模量字段层级（颗粒 vs 极片）

JuBat 区分两个尺度的弹性模量，**不可混用**：

| 字段 | 物理对象 | 典型值 | 用于 |
|------|----------|--------|------|
| `PE.E` / `NE.E` | 活性物质**颗粒**（particle） | ~1e10 Pa | 颗粒扩散应力（`Calstressdisp`） |
| `PE.E_coat` / `NE.E_coat` | **极片**涂层（coating） | ~5e8 Pa | CZM 应变驱动、二维宏观应力 |
| `SP.E` / `PCC.E` / `NCC.E` | 隔膜/集流体连续层 | — | 同极片，参与厚度加权 |
| `scale.E_p` / `scale.E_n` | 颗粒扩散应力归一化尺度 | cs_max·R·T_ref | 颗粒应力无量纲化 |
| `scale.E_coat` | 极片宏观模量参考（全叠合厚度加权） | — | CZM/宏观应力无量纲化 |

**CZM/二维宏观应力统一入口**：`compute_czm_params_per_interface(case)`（`src/CouplingState.jl:302`）→ 返回 `CzmParamCache`，按界面（`:PE_PCC` / `:NE_NCC`）提供 `CzmInterfaceParams`。**E_eff 按界面直接取涂层模量（PE-PCC 用 `PE.E_coat`、NE-NCC 用 `NE.E_coat`），不做全叠合厚度加权**；全叠合厚度加权仅保留在参考尺度 `scale.E_coat` 的定义中。

**防御**：缺失 `E_coat` 时 `ChooseCell` 触发 `@warn`，`compute_czm_params_per_interface`（`src/CouplingState.jl:307-313`）与 `thermal_diffusion_stress_2D` 入口（`src/Mechanical.jl:166-167`）处 `@assert` 拦截。

详见 `md/15_颗粒与极片模量区分.md`。

### 9.5 planning-with-files 文件存放约定

- 使用 `planning-with-files` 技能生成的 `task_plan.md`、`findings.md`、`progress.md`，统一保存到 `docs/planning-with-files/<中文任务名>/`。
- `<中文任务名>` 采用与 `docs/planning-with-files/` 现有子目录一致的简短中文主题命名，例如 `代码简化计划评审`。
- 不得将这三个文件直接创建在项目根目录；若误建，应迁移到对应任务子目录，并避免覆盖已有任务记录。
- `docs/planning-with-files/index.md` 是规划文件总索引；新增、迁移或删除任务目录/文件后，须同步更新任务说明、时间、Git 修改次数与跟踪状态。

### 9.6 代码简化基线约定

- `example/testexample.jl` 是 `src/` 代码简化的强制行为基线，基线档案位于 `Simplify/baseline/testexample/`，总入口为 `Simplify/baseline.md`。
- 每个简化批次前后均须以 Julia 1.11.2、单线程、`GKSwstype=100` 和 `--startup-file=no` 运行同一入口与参数。
- 退出码、网格/步数、`metrics.toml` 中的科学结果（按记录精度）以及结果 PNG SHA-256 必须与基线一致；运行耗时仅供参考，不作严格判定。
- 任一强制指标不一致时，停止后续简化并定位或回退该批次，不得以“数值接近”代替基线一致。

### 9.7 严格契约判断的后续简化计划

- 当前为暴露静默回退而加入的类型、维度、有限性和收敛判断应先保留，例如 `Variable_update!` 对历史矩阵与当前变量形状的显式检查。
- 待完整正常路径测试与强制行为基线全部通过，并确认实际工程流程不会触发这些错误判断后，可在独立的代码简化批次中删除已证明冗余的显式分支，以保持库实现简洁。
- 删除判断时必须依赖直接赋值、索引或求解操作的原生失败语义；不得重新引入补零、保留旧值、截断数组、取首元素或继续计算等静默回退。
- 简化前后仍须遵守第 9.6 节基线约定；科学结果、网格/步数和 PNG SHA-256 任一变化都表示该判断不可删除。

### 9.8 SP–涂层接触与摩擦后续开发建议

- 当前力学网格只建立 Φ 匝间节点配对，尚未装配 SP–PE、SP–NE 的单边法向接触、相对滑移与 Coulomb 摩擦；不得把共享节点的完美粘结结果解释为已实现接触。
- 后续实现应复用 `CzmSubmesh.phi_pairs` 及内部 SP–涂层界面双侧节点，单独建立接触状态、间隙、法向罚/增广 Lagrange 与切向 stick-slip 方程，不得并入 COH2D4 牵引–分离本构。
- SP–涂层没有可靠 $G_c$、$T_{\max}$ 时不得伪造 CZM 参数；基准摩擦系数及敏感性范围应由实验或明确文献给出。
- 接触扩展必须新增：零穿透、零抗拉、摩擦锥、stick/slip 切换、Φ 配对完整性及 `example/testexample.jl` 强制基线测试。

---

## 10. 示例文件

| 文件 | 说明 |
|------|------|
| `example/minimal_example.jl` | 基础仿真 |
| `example/SPMe_Thermal_example.jl` | 电化学-热耦合 |
| `example/czm_cycle_example.jl` | CZM 循环仿真 |
| `example/testexample.jl` | 全耦合仿真 |
| `example/jellyroll_stress_displacement.jl` | 二维应力/位移场后处理（无 CZM） |
