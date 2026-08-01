# ring.jl

- **源文件**: `src/ring.jl`
- **行数**: 71 行
- **函数/struct 计数**: 1 个独立函数
- **职责**: 圆环（极坐标）网格生成器——按等角度 + 等径向分割构造 Q4 网格，返回网格 + 内外边界节点索引
- **相关技术文档**: `md/12_热模型验证方案.md`

## 数据结构

无独立 struct 定义。本文件构造并返回 `Mesh`（定义于 `src/SetMesh.jl`）。

## 函数清单

### `ring_mesh(param; ntheta, nr, phase, gsorder) -> NamedTuple` — L7-L70

生成圆环 Q4 网格（极坐标规则网格）。

- 参数：`param.cell.Rin` / `Rout` 控制半径，`ntheta`（周向分段，默认 40），`nr`（径向分段，默认 20），`phase`（角度相位，默认 0），`gsorder`（高斯阶，默认 2）
- L9-L14：参数校验（`Rout > Rin`、`ntheta >= 3`、`nr >= 1`），失败抛 `error`
- L16-L30：节点生成——`r` 等距分布于 `[Rin, Rout]`，`theta` 等距分布于 `[0, 2π)`（不闭合，最后一点由 `it_next` 回绕到 index 1）；索引函数 `idx(ir, it) = (ir-1)*ntheta + it`
- L32-L48：单元连接——每个 Q4 单元四节点 `(ir, it) → (ir+1, it) → (ir+1, it+1) → (ir, it+1)`，`it == ntheta` 时 `it_next = 1` 实现周向环绕
- L51-L61：方向校正——用 Shoelace 面积 `0.5·Σ(x_i·y_{i+1} - x_{i+1}·y_i)` 检测方向，`area < 0` 时交换节点 2/4 保证逆时针（detJ > 0）
- L63-L64：调 `GetGS` 构造高斯点，封装为 `Mesh("Q4", ...)`
- L66-L69：返回 NamedTuple `(mesh, inner_nodes, outer_nodes, r, theta, nr, ntheta)`；`theta` 返回值附加首点 +2π 用于绘图闭合
- 跨文件依赖：`GetGS`、`Mesh`（`src/SetMesh.jl`）、`ShapeFunction2D`（经 `GetGS` 调用）

## 省略项

无。本文件唯一函数已独立列出。

### [DEBUG]

无。本文件无 `println` / `@show` / 调试用途的 `@info` / `@warn`。

### [PLACEHOLDER]

无。本文件无 `TODO` / `FIXME` / 兜底值 / 静默 try-catch；初始值 `zeros(...)` 是合理的零初始化。

### [COMPLEX-CHECK]

无。本文件条件表达式均为单条件（`Rout > Rin`、`ntheta >= 3`、`nr >= 1`、`it == ntheta`、`area < 0.0`），无 ≥3 条件的 `&&` 链或 ≥3 层嵌套。
