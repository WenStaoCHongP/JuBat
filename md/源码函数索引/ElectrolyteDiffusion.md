# ElectrolyteDiffusion.jl

- **源文件**: `src/ElectrolyteDiffusion.jl`
- **行数**: 30 行
- **函数/struct 计数**: 1 个独立函数；0 个 struct
- **职责**: 装配电解液 Li⁺ 浓度扩散方程的质量矩阵 M 与刚度矩阵 K，按负极/隔膜/正极分段取孔隙率与有效扩散系数
- **相关技术文档**: `md/04_电化学模型_SPMe.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `ElectrolyteDiffusion(param::Params, mesh::Mesh, mlen, variables)` — L1-L30

装配 `M du/dt = Ku + F` 的 M、K。

- 从 `variables` 取三段（负极/隔膜/正极）高斯点浓度 `ce_{n,sp,p}_gs` 与温度 `T`
- M 系数：分段孔隙率 `eps`，`coeff = [eps_ne; eps_sp; eps_pe]·weight·detJ`
- K 系数：`De_eff = De(ce,T)·eps^brugg·Arrhenius(Eac_D,T)`，分段拼接后 `coeff = -De_eff·weight·detJ`
- 跨文件依赖：`Assemble`、`Arrhenius`

## 省略项

无。

### [DEBUG]

无。

### [PLACEHOLDER]

无。

### [COMPLEX-CHECK]

无。

---
