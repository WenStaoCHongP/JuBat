# Thermal.jl

- **源文件**: `src/Thermal.jl`
- **行数**: 81 行
- **函数/struct 计数**: 1 个独立函数；0 个 struct
- **职责**: 集总参数（lumped）热模型，根据电化学模型类型（SPM / SPMe / 完整 P2D）计算电池单体总产热（反应热 + 可逆热 + 欧姆热）与对外散热，返回质量矩阵与净热源
- **相关技术文档**: `md/05_热模型_二维分布式.md`、`md/07_界面热阻模型.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `ThermalLumped(case::Case, variables::Dict{String, Union{Array{Float64},Float64}})` — L1-L80

集总热模型入口，返回 `(MT, Q_in - q)`，其中 `MT = m·c·heat_Q` 为 1×1 质量矩阵、`Q_in - q` 为净热源（产热 − 对流散热）。

- 根据 `case.opt.model` 分三个分支计算产热分量 `Q_rxn`/`Q_ohm`/`Q_rev`：
  - **SPM 分支** L6-L14：仅反应热 + 可逆热，`Q_ohm = 0`（SPM 无电解液描述）
  - **SPMe 分支** L15-L42：基于解析公式计算电解液电势降 `dphi_e`，进而 `Q_ohm = -I·dphi_e`；含 Bruggeman 有效电导率 `κ_eff = κ(ce,T)·ε^brugg`
  - **默认（P2D）分支** L43-L75：通过高斯点积分 `IntV` 计算各相（固相 n/p、电解液 e）欧姆热；电解液欧姆热 `Q_ohm_e` 含浓度梯度耦合项
- 净热源：`Q_in - q`，其中散热 `q = h·A_cool·(T − T_amb)`
- 写入：`variables["thermal lumped internal heat"] = [Q_in]`
- 跨文件依赖：`IntV`（积分工具），参数访问 `param.cell.mass/heat_Q/cooling_surface/T_amb/h`、`param.PE/NE/SP/EL`

## 省略项

无。

### [DEBUG]

无。

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L12 | `Q_ohm = 0`（SPM 分支） | SPM 模型无电解液物理，强制将欧姆热置零会低估总产热；如需精确热源应改用 SPMe 或 P2D 分支。属模型近似而非临时占位 |

### [COMPLEX-CHECK]

无。注：L73 单行电解液欧姆热表达式 ~175 字符，但属单赋值语句（非条件链），不触发 ≥3 `&&` 或 ≥3 层嵌套规则。
