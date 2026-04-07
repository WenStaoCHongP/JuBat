# SetCase.jl 优化方案

> 日期: 2026-04-01
> 文件: `src/SetCase.jl`
> 状态: 修改 (106→111 行)
> main 分支行数: 104

---

## 1. main 分支现状

### 1.1 Case 结构体

```julia
mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Union{Array{Int64}, Int64}}
end
```

5 个字段，职责清晰：
- `param_dim`: 有量纲参数
- `param`: 无量纲参数
- `opt`: 求解选项
- `mesh`: 字符串键网格字典
- `index`: DOF 索引字典

### 1.2 SetCase() 函数 (106 行)

流程：`NormaliseParam` → 按 model 类型建 mesh/index → `Case()` 构造

支持的 model 路径：
- SPM / SPMe → 粒子网格 + (SPMe: 电解液网格)
- P2D / sP2D → 粒子网格 + 电解液 + 电极网格
- lumped 热 → index 追加 temperature

---

## 2. 当前分支变更 (Parameters_Design)

### 2.1 变更清单

| 位置 | 变更 |
|------|------|
| 行 41 | `"sP2D"` 条件移除 → 只留 `"P2D"` |
| 行 76-82 | 新增 `distributed2D` thermalmodel 分支 |
| 行 83 | `"sP2D"` 条件移除 → 只留 `"P2D"` |
| 行 92-93 | 新增 `multi_spme_layout` 字段初始化 |
| 行 98-105 | Case 新增第 6 字段 |

### 2.2 distributed2D 分支

```julia
elseif opt.thermalmodel == "distributed2D"
    index["temperature"] = [v0 + 1]
```

仅注册一个 temperature 索引（代表值），实际温度 DOF 在 multi_spme_layout 中管理。

---

## 3. 优化方案

### 3.1 约束

- **SPM/SPMe/P2D 的 mesh/index 构建逻辑完全不动**
- 仅修改 Parameters_Design 新增的部分

### 3.2 Case 结构体：3 字段替代 1 个 Dict

```julia
# 旧:
mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Union{Array{Int64}, Int64}}
    multi_spme_layout::Dict{String,Any}
end

# 新:
mutable struct Case
    param_dim::Params
    param::Params
    opt::Option
    mesh::Dict{String, Mesh}
    index::Dict{String, Union{Array{Int64}, Int64}}
    # ---- 类型化替代 ----
    layout::Union{Nothing, MultiSPMeLayout}    # 布局索引
    geometry::Union{Nothing, MeshGeometry}     # 几何拓扑
    czm_mesh::Union{Nothing, CohesiveMesh}     # CZM 网格
end

# 便捷构造器（main 分支 5 参数路径）
Case(pd, p, o, m, i) = Case(pd, p, o, m, i, nothing, nothing, nothing)
```

### 3.3 SetCase() 函数改动

```julia
# 旧 (行 93):
case = Case(param_dim, param, opt, mesh, index, Dict{String,Any}())

# 新:
case = Case(param_dim, param, opt, mesh, index)
# 使用 5 参数便捷构造器，layout/geometry/czm_mesh 自动为 nothing
```

distributed2D 分支不变：

```julia
elseif opt.thermalmodel == "distributed2D"
    index["temperature"] = [v0 + 1]
```

### 3.4 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| Case 字段数 | 6 (含 Dict{String,Any}) | 8 (3 个 Union{Nothing, T}) |
| 类型安全 | 无 (Any) | 有 (struct 字段) |
| main 分支影响 | N/A | 零（便捷构造器兼容） |
| 拼写错误检测 | 运行时 | 编译时 |

---

## 4. 不动的部分（明确列出）

以下代码**完全不动**：
- `NormaliseParam(param_dim)` 调用
- SPM 分支 (行 10-22)
- SPMe 分支 (行 23-40) 的 mesh/index 构建
- P2D 分支 (行 41-74) 的 mesh/index 构建
- lumped thermal 分支 (行 76-78)
- P2D 电位 index (行 83-91)
