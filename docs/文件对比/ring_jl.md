# ring.jl

## 文件状态: 新增 (Parameters_Design分支)

## 文件概况
- 行数: 83
- 路径: `src/ring.jl`

### 主要函数/方法列表

| 函数签名 | 行号 | 说明 |
|----------|------|------|
| `ring_mesh(param; ntheta, nr, phase, gsorder)` | L7 | 生成圆环区域Q4网格，返回mesh + 边界节点 + 网格参数 |

## 功能描述

本文件实现了圆环（Ring）区域的规则Q4网格生成器，用于圆柱电池热模型的简化几何建模。主要功能包括：

1. **网格生成**（`ring_mesh`）：
   - 在 `[Rin, Rout]` x `[0, 2pi]` 区域生成结构化Q4网格
   - 使用 `nr+1` 个径向节点和 `ntheta` 个周向节点
   - 节点按极坐标排列后转换为笛卡尔坐标
   - 最后一个周向节点不与第一个闭合（开环），简化边界条件处理

2. **单元方向保证**：
   - 通过计算有向面积检查单元方向
   - 负面积时交换节点2和4，确保逆时针方向（detJ > 0）

3. **返回值**：
   - `mesh`: Mesh对象（Q4网格）
   - `inner_nodes`: 内壁节点索引
   - `outer_nodes`: 外壁节点索引
   - `r`: 径向坐标数组
   - `theta`: 周向坐标数组（含闭合点 theta[1] + 2pi）
   - `nr`, `ntheta`: 网格分辨率参数

与 `Jellyrollmodel.jl` 的螺旋网格不同，ring.jl 生成的是简单的同心圆环网格，不包含卷绕层界面信息。适用于验证热模型（如圆环精确解对比）或作为简化几何使用。

## 依赖关系

### 该文件依赖
- `src/SetMesh.jl` — `Mesh`结构体、`GetGS`高斯积分函数（通过条件include）
- `LinearAlgebra` — Julia标准库

### 哪些文件调用该文件
- `src/JuBat.jl` — `include("ring.jl")`（L30）
- `example/热模块验证/thermal_verify.jl` — 使用 `ring_mesh` 生成验证用圆环网格
- 简化圆柱模型场景 — 使用Ring参数集 + ring_mesh 替代完整的Jellyroll螺旋网格

## 耦合分析

本文件是distributed2D热模型的**简化几何提供者**，主要用于：

- **与distributed2D热模型耦合**：生成的圆环网格可直接用作 `case.mesh["thermal2D"]`。适用于不需要考虑卷绕层界面的简化热分析场景。

- **验证用途**：由于圆环区域存在解析解（如稳态径向导热），ring_mesh 可用于验证FEM热求解器的正确性。参见 `example/热模块验证/thermal_verify.jl`。

- **与Jellyroll参数一致**：Ring参数集（`parameters/Ring.jl`）使用与Jellyroll相同的电极厚度和尺度参数，确保 `scale.L` 一致，热源量纲无需桥接。

- **不参与CZM耦合**：由于圆环网格不包含界面节点对信息，无法创建CZM网格。因此ring.jl仅用于纯热分析或电-热耦合（不含力学）。
