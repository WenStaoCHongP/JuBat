# SPMe 模型欧姆热计算修改计划


**创建日期**: 2026-03-24
**状态**: pending
**关联文档**: `md/05_热模型_二维分布式.md` §2.2.1


---

## 1. 目标


将 `Thermal.jl` 第 12 行 `Q_ohm = 0` 修改为 PyBaMM SPMe 模型的欧姆热计算方式，使集总热模型的欧姆热与 PyBaMM 对齐。


---

## 2. 问题分析


### 2.1 当前代码问题


**文件**: `src/Thermal.jl` 第 6-14 行


```julia
if case.opt.model == "SPM" || case.opt.model == "SPMe"
    I_app = variables["cell current"]
    eta_n = variables["negative electrode overpotential"][1]
    eta_p = variables["positive electrode overpotential"][end]
    csn_surf = variables["negative particle surface lithium concentration"][1]
    csp_surf = variables["positive particle surface lithium concentration"][end]
    Q_ohm = 0  # ❌ SPMe 模型中欧姆热被忽略
    Q_rxn = abs(I_app * (eta_p - eta_n))
    Q_rev = abs(I_app) * T * (param.PE.dUdT(csp_surf) - param.NE.dUdT(csn_surf))
```


**问题**：
- SPMe 模型中 `Q_ohm = 0`，完全忽略欧姆热
- 导致热模型与 PyBaMM 验证结果差异较大


### 2.2 PyBaMM 的欧姆热计算方法（电势差积分）


> **重要**：PyBaMM 使用的是**电势差积分**方法，而不是简化公式 $I^2/(3\sigma)$。


PyBaMM SPMe 模型的欧姆热计算步骤：


**Step 1: 从 SPMe 解析解重构电势场**


SPMe 虽然不求解电势 PDE，但通过解析公式构造电势分布：


**固相电势降**：
$$\Delta\phi_s = \frac{I}{3}\left(\frac{L_n}{\sigma_n} + \frac{L_p}{\sigma_p}\right)$$


**电解液电势降**：
$$\Delta\phi_e = 2T(1-t^+)\frac{\Delta c_e}{c_{e,0}} - I \cdot R_{el} - \Delta\phi_s$$


其中电解液电阻：
$$R_{el} = \frac{L_n}{3\kappa_n^{eff}} + \frac{L_{sp}}{\kappa_{sp}^{eff}} + \frac{L_p}{3\kappa_p^{eff}}$$


**Step 2: 计算欧姆热（功率积分）**


PyBaMM 对各区域进行积分计算欧姆热：


**固相欧姆热**（负极和正极）：
$$Q_{ohm,s,k} = \int_{V_k} \sigma_k^{eff} (\nabla\phi_s)^2 dV = I \cdot \Delta\phi_{s,k}$$


**电解液相欧姆热**：
$$Q_{ohm,e} = \int_{V_e} \kappa^{eff} (\nabla\phi_e)^2 dV = I \cdot \Delta\phi_e$$


**总欧姆热**：
$$Q_{ohm} = I \cdot (\Delta\phi_s + \Delta\phi_e)$$


**关键区别**：
- PyBaMM 用 **$I \cdot \Delta\phi$**（电势差积分）
- JuBat 简化公式用 **$I^2/(3\sigma)$**（电阻功率）


两者在**线性电势分布**时数学等价，但 PyBaMM 的方法能捕捉：
- 电解液浓度变化对电势梯度的影响
- 反应电流分布不均匀的影响


**数值示例（1.0C 放电）**：

| 时间 (s) | PyBaMM Q_ohm (W) | 说明 |
|----------|-----------------|------|
| 0 | 0.135 | 初始状态 |
| 500 | 0.311 | 电解液浓度变化 |
| 999.5 | 0.313 | 接近稳态 |
| 3569.5 | 0.316 | 放电结束 |

> **数据来源**: `output/thermal_equivalent_lumped_compare.csv`


### 2.3 JuBat 已有的相关计算


在 `SPMe.jl` 第 165-172 行，已有电势差计算：


```julia
# 固相电势降
dphi_S = I_app / 3 * (param.NE.thickness / param.NE.sig + param.PE.thickness / param.PE.sig)

# 电解液有效电导率
kappa_ne = param.EL.kappa(param.EL.ce0, T) * param.NE.eps ^ param.NE.brugg
kappa_pe = param.EL.kappa(param.EL.ce0, T) * param.PE.eps ^ param.PE.brugg
kappa_sp = param.EL.kappa(param.EL.ce0, T) * param.SP.eps ^ param.SP.brugg

# 电解液电阻
R_EL = param.NE.thickness / kappa_ne / 3.0 + param.SP.thickness / kappa_sp + param.PE.thickness / kappa_pe / 3.0

# 电解液电势降
dphi_e = 2.0 * T * (1 - param.EL.tplus) * (csp_av - csn_av)/param.EL.ce0 .- I_app * R_EL .- dphi_S
```


**可复用**：`dphi_S` 和 `R_EL` 可用于计算欧姆热。


---
## 3. 修改方案
### 3.1 方案 A：直接功率公式（推荐）
**优点**：简单直接，与 PyBaMM 完全一致
**缺点**：假设线性电势分布
**公式**：
$$Q_{ohm} = I^2 \cdot \left(\frac{L_n}{3\sigma_n^{eff}} + \frac{L_p}{3\sigma_p^{eff}} + R_{el}\right)$$
其中：
- $\sigma_n^{eff} = \sigma_n \cdot \epsilon_{s,n}$ （负极有效电导率）
- $\sigma_p^{eff} = \sigma_p \cdot \epsilon_{s,p}$ （正极有效电导率）
- $R_{el}$ 为电解液电阻（如上定义）
### 3.2 方案 B：电势差积分
**优点**：更接近物理本质
**缺点**：需要额外计算
**公式**：
$$Q_{ohm} = I \cdot \Delta\phi_{ohm}$$
其中 $\Delta\phi_{ohm}$ 是总欧姆电势降（不含反应过电势和开路电势）。
### 3.3 推荐方案
**采用方案 A**，原因：
1. 与 PyBaMM 实现一致
2. 公式简单明确
3. 可直接复用 `SPMe.jl` 中已有的参数


---
## 4. 详细修改计划
### Phase 1: 添加欧姆热计算函数 [pending]
**新文件**: `src/ThermalOhmic.jl`（可选）或在 `Thermal.jl` 中添加
**函数定义**：
```julia
"""
    compute_ohmic_heat_spme(param, I_app, T)
计算 SPMe 模型的欧姆热。
# 参数
- `param`: 归一化参数结构
- `I_app`: 无量纲电流
- `T`: 无量纲温度
# 返回
- `Q_ohm`: 无量纲欧姆热功率
# 公式
Q_ohm = I² × R_ohm_total
其中 R_ohm_total = R_s_n + R_s_p + R_el
- R_s_n = L_n / (3 × σ_n^eff)
- R_s_p = L_p / (3 × σ_p^eff)
- R_el = L_n/(3×κ_n^eff) + L_sp/κ_sp^eff + L_p/(3×κ_p^eff)
"""
function compute_ohmic_heat_spme(param, I_app, T)
    # 1. 固相有效电导率
    sig_n_eff = param.NE.sig * param.NE.eps_s  # 负极
    sig_p_eff = param.PE.sig * param.PE.eps_s  # 正极
    
    # 2. 电解液有效电导率（温度相关）
    kappa_ne = param.EL.kappa(param.EL.ce0, T) * param.NE.eps ^ param.NE.brugg
    kappa_pe = param.EL.kappa(param.EL.ce0, T) * param.PE.eps ^ param.PE.brugg
    kappa_sp = param.EL.kappa(param.EL.ce0, T) * param.SP.eps ^ param.SP.brugg
    
    # 3. 各部分电阻（无量纲）
    R_s_n = param.NE.thickness / (3.0 * sig_n_eff)   # 负极固相电阻
    R_s_p = param.PE.thickness / (3.0 * sig_p_eff)   # 正极固相电阻
    R_el = (param.NE.thickness / kappa_ne / 3.0 
          + param.SP.thickness / kappa_sp 
          + param.PE.thickness / kappa_pe / 3.0)     # 电解液电阻
    
    # 4. 总欧姆电阻
    R_ohm_total = R_s_n + R_s_p + R_el
    
    # 5. 欧姆热功率
    Q_ohm = I_app^2 * R_ohm_total
    
    return Q_ohm
end
```
### Phase 2: 修改 Thermal.jl [pending]
**文件**: `src/Thermal.jl` 第 6-14 行
**修改前**：
```julia
if case.opt.model == "SPM" || case.opt.model == "SPMe"
    I_app = variables["cell current"]
    eta_n = variables["negative electrode overpotential"][1]
    eta_p = variables["positive electrode overpotential"][end]
    csn_surf = variables["negative particle surface lithium concentration"][1]
    csp_surf = variables["positive particle surface lithium concentration"][end]
    Q_ohm = 0  # ❌ 修改此处
    Q_rxn = abs(I_app * (eta_p - eta_n))
    Q_rev = abs(I_app) * T * (param.PE.dUdT(csp_surf) - param.NE.dUdT(csn_surf))
```
**修改后**：
```julia
if case.opt.model == "SPM" || case.opt.model == "SPMe"
    I_app = variables["cell current"]
    eta_n = variables["negative electrode overpotential"][1]
    eta_p = variables["positive electrode overpotential"][end]
    csn_surf = variables["negative particle surface lithium concentration"][1]
    csp_surf = variables["positive particle surface lithium concentration"][end]
    
    # 计算 SPMe 模型的欧姆热（PyBaMM 方法）
    Q_ohm = compute_ohmic_heat_spme(param, I_app, T)
    
    Q_rxn = abs(I_app * (eta_p - eta_n))
    Q_rev = abs(I_app) * T * (param.PE.dUdT(csp_surf) - param.NE.dUdT(csn_surf))
```
### Phase 3: 验证修改 [pending]
**验证脚本**: 新建 `example/热模块验证/ohmic_heat_verification.jl`
```julia
# 验证欧姆热计算
include("src/JuBat.jl")
using .JuBat
# 1. 创建案例
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
opt.model = "SPMe"
opt.thermal_enabled = true
opt.thermalmodel = "lumped"
case = JuBat.SetCase(param_dim, opt)
# 2. 求解
result = JuBat.Solve(case)
# 3. 与 PyBaMM 数据对比
pybamm_data = CSV.read("src/data/pybamm_SPMe_LGM50_1.0C.csv", DataFrame)
# 4. 绘制对比图
# ... 欧姆热随时间变化对比
```
**验收标准**：

| 时间 (s) | PyBaMM Q_ohm (W) | JuBat Q_ohm (W) | 误差 |
|----------|-----------------|-----------------|------|
| 0 | 0.135 | < 0.14 | < 5% |
| 500 | 0.311 | < 0.32 | < 5% |
| 999.5 | 0.313 | < 0.32 | < 5% |

> **数据来源**: `output/thermal_equivalent_lumped_compare.csv`

### Phase 4: 文档更新 [pending]
**需要更新的文档**：
1. `md/05_热模型_二维分布式.md`：
   - 更新 §2.2.1 中的 JuBat 公式
   - 添加实现说明
2. `docs/thermal_verify/ohmic_heat_modification_plan.md`：
   - 更新状态为 completed
---
## 5. 公式推导记录
### 5.1 固相欧姆热
**假设**：电流均匀分布，电势线性分布
**负极**：
$$Q_{ohm,s,n} = \int_0^{L_n} \sigma_n^{eff} \left(\frac{d\phi_s}{dx}\right)^2 A dx$$
线性分布：$\phi_s(x) = \phi_s(0) + \frac{\Delta\phi_s}{L_n} x$
$$\frac{d\phi_s}{dx} = \frac{\Delta\phi_s}{L_n} = \frac{I}{\sigma_n^{eff} A}$$
$$Q_{ohm,s,n} = \sigma_n^{eff} \cdot \left(\frac{I}{\sigma_n^{eff} A}\right)^2 \cdot A \cdot L_n = \frac{I^2 L_n}{\sigma_n^{eff} A}$$
**考虑有限元积分因子**：$\times \frac{1}{3}$
$$Q_{ohm,s,n} = \frac{I^2 L_n}{3 \sigma_n^{eff} A}$$
**无量纲化**：
$$Q_{ohm,s,n}^* = \frac{(I^*)^2 L_n^*}{3 \sigma_n^{eff,*}}$$
### 5.2 电解液相欧姆热
**类似推导**：
$$Q_{ohm,e} = I^2 \cdot R_{el}$$
$$R_{el} = \frac{L_n}{3\kappa_n^{eff}} + \frac{L_{sp}}{\kappa_{sp}^{eff}} + \frac{L_p}{3\kappa_p^{eff}}$$
**注意**：隔膜层没有因子 1/3，因为隔膜中电势梯度是均匀的（无反应电流）。
---
## 6. 注意事项
### 6.1 无量纲处理
JuBat 内部使用无量纲变量，需要确保：
- `I_app` 是无量纲电流（已归一化）
- `param.NE.thickness` 是无量纲厚度
- `param.NE.sig` 是无量纲电导率
- 返回的 `Q_ohm` 是无量纲功率
### 6.2 与分布式模型的一致性
修改后，`Thermal.jl` 的集总欧姆热应与 `ThermalDistributed.jl` 中简化公式一致：
```julia
# ThermalDistributed.jl 中的公式
Q_ohm_s_NE = I_local^2 / (3.0 * sig_n_eff)
Q_ohm_e_NE = I_local^2 / (3.0 * kappa_ne)
```
### 6.3 SPM 模型
SPM 模型不考虑电解液，因此：
- 如果 `case.opt.model == "SPM"`，`R_el = 0`
- 或者 SPM 模型保持 `Q_ohm = 0`
---
## 7. 相关文件
| 文件 | 说明 | 状态 |
|------|------|------|
| `src/Thermal.jl` | 需要修改的文件 | pending |
| `src/SPMe.jl` | 参考电势计算 | 无需修改 |
| `src/ThermalDistributed.jl` | 分布式欧姆热 | 无需修改 |
| `md/05_热模型_二维分布式.md` | 技术文档 | 需要更新 |
| `src/data/pybamm_SPMe_LGM50_*.csv` | 验证数据 | 无需修改 |
---
## 8. 执行进度
- [ ] Phase 1: 添加 `compute_ohmic_heat_spme` 函数
- [ ] Phase 2: 修改 `Thermal.jl` 第 12 行
- [ ] Phase 3: 创建验证脚本并运行
- [ ] Phase 4: 更新技术文档
---
## 9. 参考资料
1. PyBaMM SPMe 模型文档: https://docs.pybamm.org/en/stable/source/examples/notebooks/models/SPMe.html
2. `md/05_热模型_二维分布式.md` §2.2.1 欧姆热计算方法对比
3. `src/SPMe.jl` 第 165-172 行电势计算
