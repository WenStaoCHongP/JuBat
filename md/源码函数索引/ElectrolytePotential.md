# ElectrolytePotential.jl

- **源文件**: `src/ElectrolytePotential.jl`
- **行数**: 23 行
- **函数/struct 计数**: 1 个独立函数；0 个 struct
- **职责**: 装配电解液液相电势（电荷守恒）方程的质量矩阵 M（零）与刚度矩阵 K，按三段取有效离子电导率，用于 P2D 模型
- **相关技术文档**: `md/04_电化学模型_SPMe.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `ElectrolytePotential(param::Params, mesh::Mesh, mlen, variables)` — L1-L23

装配 `M du/dt = Ku + F` 的 M、K。

- M：`spzeros(mlen, mlen)`（电势方程无时间项）
- 从 `variables` 取三段高斯点浓度 `ce_{n,sp,p}_gs` 与温度 `T`
- K 系数：`kappa(ce,T)·eps^brugg`（有效离子电导率，Bruggeman 校正），分段拼接后 `coeff = -[...]·weight·detJ`
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
