# ElectrodePotential.jl

- **源文件**: `src/ElectrodePotential.jl`
- **行数**: 19 行
- **函数/struct 计数**: 1 个独立函数；0 个 struct
- **职责**: 装配电极固相电势（电荷守恒）方程的质量矩阵 M（零）与刚度矩阵 K，用于 P2D 模型
- **相关技术文档**: `md/04_电化学模型_SPMe.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `ElectrodePotential(electrode::Electrode, mesh::Mesh, mlen, T)` — L1-L19

装配 `M du/dt = Ku + F` 的 M、K。

- M：`spzeros(mlen, mlen)`（电势方程无时间项）
- K 系数：`sig_eff = sig·eps_s`（有效电子电导率，含固相体积分数校正），`coeff = sig_eff·weight·detJ`
- 默认参数：`T=1.0`（当前未在公式中使用，预留温度依赖接口）
- 跨文件依赖：`Assemble`

## 省略项

无。

### [DEBUG]

无。

### [PLACEHOLDER]

无。

### [COMPLEX-CHECK]

无。

---
