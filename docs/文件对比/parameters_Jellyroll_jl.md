# parameters/Jellyroll.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 180
- 路径: `src/parameters/Jellyroll.jl`

### 主要参数块

| 参数块 | 类型 | 说明 |
|--------|------|------|
| `PE` | `Electrode` | 正极参数（Chen2020 NMC811） |
| `NE` | `Electrode` | 负极参数（Chen2020 Graphite） |
| `EL` | `Electrolyte` | 电解液参数（LiPF6 in EC:EMC） |
| `SP` | `Separator` | 隔膜参数 |
| `PCC` | `CurrentCollector` | 正极集流体参数（Al） |
| `NCC` | `CurrentCollector` | 负极集流体参数（Cu） |
| `tab` | `Tab` | 极耳参数 |
| `cell` | `Cell` | 电池级别参数（几何、热、电压范围） |
| `cohesive` | `Cohesive` | CZM内聚力参数 |
| `scale` | `Scale` | 无量纲化尺度参数（由NormaliseParam填充） |
| `param_dim` | `Params` | 完整参数集合 |

### 关键参数值

**电极参数**：
- PE厚度: 75.6um, rs: 5.22um, cs_max: 63104 mol/m3
- NE厚度: 85.2um, rs: 5.86um, cs_max: 33133 mol/m3
- SP厚度: 12um

**力学参数**：
- PE: E=37.5GPa, nu=0.20, alphaT=1e-5/K, Omega=-7.28e-7 m3/mol
- NE: E=15.0GPa, nu=0.28, alphaT=3e-6/K, Omega=3.1e-6 m3/mol

**CZM参数**：
- sigma_max_n = 82 MPa（法向最大牵引力）
- K_n = 2.4e17 Pa/m（初始刚度）
- G_c_n = 25.3 J/m2（法向断裂能）
- tau_max_t = 0.15 MPa（切向最大牵引力）
- eta = 1.45（BK准则指数）

**极耳参数**：
- theta_pos = [15pi], theta_neg = [44pi]（极耳角度位置）
- width = 4mm, length = 7.425mm

**电池几何**：
- Rin = 1.92mm, Rout = 10.15mm（半径）
- layer = 2*(PE+NE+SP) + PCC + NCC = 2*(75.6+85.2+12) + 16 + 12 = 373.2um
- width = 65mm（轴向高度）

## 功能描述

本文件定义了Jellyroll（果冻卷）型号电池的完整物理参数集，基于Chen2020（LG M50）数据。参数涵盖：

1. **电化学参数**：颗粒扩散系数（Ds）、反应速率常数（k）、最大锂浓度（cs_max）、化学计量范围（theta_0/theta_100）、开路电位函数（U/OCV）、电解液扩散率和电导率。

2. **热学参数**：各层导热率（lambda）、密度（rho）、比热容（heat_Q）。由这些参数可计算：
   - 径向等效导热 `lambda_r = sum(t_i) / sum(t_i/lambda_i)` = 1.318 W/m/K
   - 周向等效导热 `lambda_t = sum(t_i*lambda_i) / sum(t_i)` = 25.26 W/m/K

3. **力学参数**：弹性模量（E）、泊松比（nu）、热膨胀系数（alphaT）、部分摩尔体积（Omega）。用于扩散应力和热应力计算。

4. **CZM参数**：法向/切向的双线性牵引-分离律参数，混合模式BK准则指数。

5. **界面热阻参数**：h_c0、k_air、lambda_m、beta，用于损伤影响下的界面导热计算。

6. **几何参数**：螺旋结构由 Rin, Rout, layer 定义。极耳位置通过角度 theta_pos/theta_neg 指定。

## 依赖关系

### 该文件依赖
- `src/SetParams.jl` — 所有参数类型定义（Electrode, Electrolyte, Separator, CurrentCollector, Tab, Cell, Cohesive, Scale, Params, Binder）

### 哪些文件调用该文件
- `src/ChooseCell.jl` — `ChooseCell("Jellyroll")` 时include此文件
- `src/Jellyrollmodel.jl` — 使用 param.cell.Rin, param.cell.Rout, param.cell.layer 等几何参数
- `src/Materialmatrix.jl` — 使用 param.NE.lambda, param.PE.rho 等热学参数
- `src/Mechanical.jl` — 使用 param.NE.E, param.PE.Omega 等力学参数
- `src/CycleSolver.jl` — 使用 param_dim.NE.cs_max 等维度参数

## 耦合分析

本文件是 **全耦合模型** 的物理参数基础：

- **与multi-SPMe耦合**：电极电化学参数（Ds, k, cs_max, U, theta范围）直接决定SPMe模型的动力学行为。
- **与distributed2D热模型耦合**：热学参数（lambda, rho, heat_Q）用于计算各向异性导热和体积热容。层厚度权重决定了分层热源的空间分布。
- **与CZM耦合**：cohesive参数定义了界面脱粘行为。力学参数（E, nu, alphaT, Omega）决定了热-化学载荷到应力/损伤的转化关系。
- **与Jellyroll几何耦合**：cell.Rin, cell.Rout, cell.layer 定义了螺旋几何的基本参数。tab.theta_pos/theta_neg 定义了极耳位置。
