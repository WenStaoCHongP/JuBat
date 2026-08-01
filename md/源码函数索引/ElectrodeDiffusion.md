# ElectrodeDiffusion.jl

- **源文件**: `src/ElectrodeDiffusion.jl`
- **行数**: 18 行
- **函数/struct 计数**: 1 个独立函数；0 个 struct
- **职责**: 装配电极颗粒（球形径向）扩散方程的质量矩阵 M 与刚度矩阵 K，支持力-化学耦合（应力影响扩散系数）
- **相关技术文档**: `md/04_电化学模型_SPMe.md`

## 数据结构

本文件无独立 struct 定义。

## 函数清单

### `ElectrodeDiffusion(electrode::Electrode, mesh::Mesh, mlen, c, theta, T)` — L1-L18

装配 `M du/dt = Ku + F` 的 M、K（球形坐标，含 r² 权重）。

- M 系数：`coeff = x²·weight·detJ`（球坐标体积权重）
- K 系数：`Ds_eff = Ds·Arrhenius(Eac_D,T)·(1 + θ·c)`，`coeff = -Ds_eff·x²·weight·detJ`
- 应力耦合：`theta` 为应力耦合系数（来自 `Mechanicaloutput`），缺省 0
- 默认参数：`c=zeros, theta=0.0, T=1.0`
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
