# Assemble.jl

- **源文件**: `src/Assemble.jl`
- **行数**: 41 行
- **函数/struct 计数**: 2 个函数（无 struct）
- **职责**: 有限元全局矩阵与载荷向量装配——稀疏系统矩阵装配（`Assemble`）、1D 载荷向量装配（`Assemble1D`）
- **相关技术文档**: `md/05_热模型_二维分布式.md`、`md/04_电化学模型_SPMe.md`

## 数据结构

本文件无独立 struct 定义。返回 `SparseMatrixCSC`（系统矩阵）和 `Vector{Float64}`（载荷向量）。

## 函数清单

### `Assemble(Vi, Vj, Ni, Nj, coeff, mlen1, mlen2=mlen1) -> SparseMatrixCSC` — L1-L26

装配系统矩阵 K(Vi, Vj) = ∫ Ni·Nj·coeff·weight·detJ。

- 由 `Ni`、`Nj` 的形状推断 Gauss 点数 `gslen` 与每维节点数 `gslen1`、`gslen2`（L9-L11）
- 预分配 COO 三元组数组 `KI`、`KJ`、`KV`（L12-L14），总长度 `gslen * gslen1 * gslen2`
- 双层 for 循环（L16-L23）逐 (i, j) 对填充：行索引 `Vi[:,i]`、列索引 `Vj[:,j]`、值 `Ni[:,i] .* Nj[:,j] .* coeff`
- 末尾 `sparse(KI, KJ, KV, mlen1, mlen2)` 构造稀疏矩阵（L24）

### `Assemble1D(Vi, Ni, coeff, mlen1) -> Vector{Float64}` — L28-L41

装配载荷向量 F(Vi) = ∫ Ni·coeff·weight·detJ。

- 由 `Ni` 第一维推断 Gauss 点数 `gslen1`（L35）
- 初始化 `F = zeros(mlen1)`（L36）
- 单层 for 循环（L37-L39）：`F[Vi[i,:]] .+= Ni[i,:] .* coeff[i]`（注意是原位累加 `+=`）

## 省略项

无。所有 function 均有独立条目。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试 `@info` / `@warn`。

### [PLACEHOLDER]

无。本文件无 TODO/FIXME/占位/magic-number fallback。

### [COMPLEX-CHECK]

无。本文件无 ≥3 的 `&&` 链、≥3 层嵌套 if 或长条件表达式。
