# parameters/Ring.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 66
- 路径: `src/parameters/Ring.jl`

### 主要参数块

| 参数块 | 类型 | 说明 |
|--------|------|------|
| `PE` | `Electrode` | 正极参数（与Jellyroll一致厚度） |
| `NE` | `Electrode` | 负极参数（与Jellyroll一致厚度） |
| `EL` | `Electrolyte` | 电解液参数（仅ce0） |
| `SP` | `Separator` | 隔膜参数（仅厚度） |
| `PCC` | `CurrentCollector` | 正极集流体（空） |
| `NCC` | `CurrentCollector` | 负极集流体（空） |
| `tab` | `Tab` | 极耳（空） |
| `cohesive` | `Cohesive` | CZM参数（仅界面热阻相关） |
| `cell` | `Cell` | 电池级别参数 |
| `scale` | `Scale` | 无量纲化尺度参数 |
| `param_dim` | `Params` | 完整参数集合 |

### 关键参数值

**电极参数**（仅设置厚度和力学参数，不含电化学）：
- PE厚度: 75.6um, E=50GPa, nu=0.30, alphaT=1e-5/K
- NE厚度: 85.2um, E=20GPa, nu=0.28, alphaT=3e-6/K
- SP厚度: 12um

**电池几何**（与Jellyroll一致）：
- Rin = 1.92mm, Rout = 10.15mm
- width = 65mm, length = 1.58
- volume = pi*(Rout^2 - Rin^2)*width

**等效热学参数**（直接给出，非分层计算）：
- lambda_r = 1.318 W/m/K
- lambda_t = 25.26 W/m/K
- rho = 2813 kg/m3
- heat_Q = 860 J/kg/K

**界面热阻参数**（与Jellyroll一致）：
- h_c0 = 1e7, k_air = 0.026, lambda_m = 70nm, beta = 1.0

## 功能描述

本文件定义了圆环（Ring）型号的电池参数集，是Jellyroll参数的简化版本。设计原则：

1. **尺度一致性**：使用与Jellyroll相同的电极厚度（PE=75.6um, NE=85.2um, SP=12um），确保 `scale.L = PE + NE + SP` 一致，无量纲化后热源量纲无需桥接。

2. **热学等效**：直接给出等效热学参数（lambda_r, lambda_t, rho, heat_Q），而非像Jellyroll那样通过分层参数计算。这使得Ring模型可用于快速热分析。

3. **最小化参数**：仅设置热/力学分析所需的最少参数，不包含完整电化学参数（Ds, k, OCV等）。

4. **界面热阻**：保留了界面热阻参数（h_c0, k_air等），使得Ring模型也可测试界面热阻对热传导的影响。

## 依赖关系

### 该文件依赖
- `src/SetParams.jl` — 所有参数类型定义

### 哪些文件调用该文件
- `src/ChooseCell.jl` — `ChooseCell("Ring")` 时include此文件
- `src/ring.jl` — 使用 param.cell.Rin, param.cell.Rout 生成圆环网格
- `src/ThermalDistributed.jl` — 使用热学参数进行热分析

## 耦合分析

本文件是 **纯热分析/简化验证** 的参数集：

- **与distributed2D热模型耦合**：提供等效热学参数用于热有限元计算。由于直接给出lambda_r和lambda_t，无需 `thermal_anisotropic_conductivity_2d` 的分层计算。

- **与multi-SPMe耦合**：不直接参与。Ring参数集缺少完整电化学参数，无法驱动SPMe求解器。若需电-热耦合，应使用Jellyroll参数集。

- **与CZM耦合**：不直接参与。Ring参数集的cohesive参数仅包含界面热阻相关参数，缺少完整CZM参数（sigma_max_n, K_n, G_c_n等）。

- **与Jellyroll的关系**：Ring是Jellyroll的简化版本，使用相同的几何尺寸和尺度参数。适用于：
  - 热模型验证（圆环解析解对比）
  - 快速热分析（无需考虑卷绕层界面）
  - 热模型参数敏感性分析
