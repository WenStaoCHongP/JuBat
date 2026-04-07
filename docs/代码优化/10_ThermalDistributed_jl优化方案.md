# ThermalDistributed.jl 优化方案

> 日期: 2026-04-01 (修订)
> 文件: `src/ThermalDistributed.jl`
> 状态: 新增 (425 行)

---

## 1. 现状分析

### 1.1 函数清单

| 函数 | 行范围 | 行数 | 职责 |
|------|--------|------|------|
| `ThermalDistributed2D` | 1-47 | 47 | 2D FEM 热刚度/容量矩阵装配 |
| `apply_convection_bc` | 49-104 | 56 | Newton 对流边界条件 |
| `apply_cool_method` | 107-185 | 79 | 冷却方式 (none/tab/surface) |
| `ThermalDistributed2D_BC` | 187-219 | 33 | BC 总入口 (含 CZM 界面热阻) |
| `ThermalDistributed2D_Ring` | 221-270 | 50 | 极坐标 FEM 热模型 |
| `ThermalRing2D_BC` | 272-279 | 8 | 极坐标 BC |
| `compute_heat_sources` | 281-391 | 111 | **11 层热源逐一计算** |
| `compute_heat_sources_with_czm` | 393-425 | 33 | CZM 活跃单元过滤 |

### 1.2 核心问题

1. **`compute_heat_sources` 是最大的单函数 (111 行)**：11 层热源逐一计算 + 层权重 + 单位转换
2. **`compute_heat_sources_with_czm` 是简单包装**：仅调用 `compute_heat_sources` 后用 `active_mask` 过滤
3. **热源计算与 FEM 装配职责混杂**：ThermalDistributed.jl 应只负责 FEM 装配和 BC

---

## 2. 优化方案

### 2.1 约束

- **不新增函数**
- `compute_heat_sources` 保留在 ThermalDistributed.jl 内（不迁出）
- `compute_heat_sources_with_czm` 保留
- 仅做类型替换和命名修正

### 2.2 `compute_heat_sources` 内部修正

#### 关键：保留层权重 `fks` 和单位转换

当前代码正确使用了：
- `fks = jellyroll_element_properties(case, ...)` 获取层面积权重
- `q_rxn_ne[e] = fks[e,1] * Q_rxn_NE` 应用层权重
- `q_total[e] * (scale.L^3 / case.param_dim.cell.volume)` 单位转换

**这些逻辑必须完整保留，不能简化。** 方案仅做：

1. `case.multi_spme_layout["key"]` → `case.layout.key` 类型替换
2. `jellyroll_element_properties` 直接调用改为通过 `case.geometry` 获取 `fks`
3. 函数重命名: `compute_heat_sources` → `compute_heat_sources`（保持不变，已经是 snake_case）

### 2.3 层权重获取方式修正

```julia
# 旧: 直接调用 jellyroll_element_properties
fks = jellyroll_element_properties(case, ...)

# 新: 从预计算的 geometry 获取
fks = case.geometry.layer_weights
```

这消除了 `ThermalDistributed` 对 `Jellyrollmodel` 的直接耦合泄漏（00 文档标记的高严重度问题）。

### 2.4 `compute_heat_sources_with_czm` 保留

```julia
# 保持原样，仅做类型替换
function compute_heat_sources_with_czm(case, variables, variables_elems, I_e, Te_prev, areas, czm_mesh, mesh_th)
    q_total, q_fields = compute_heat_sources(case, variables, variables_elems, I_e, Te_prev, areas, mesh_th)
    # ... active_mask 过滤逻辑保持 ...
end
```

### 2.5 FEM 装配函数不做修改

| 函数 | 操作 |
|------|------|
| `ThermalDistributed2D` | 不动 |
| `apply_convection_bc` | 不动 |
| `apply_cool_method` | 不动 |
| `ThermalDistributed2D_BC` | 不动 |
| `ThermalDistributed2D_Ring` | 不动 |
| `ThermalRing2D_BC` | 不动 |

---

## 3. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| ThermalDistributed.jl 行数 | 425 | ~415 (仅类型替换) |
| 热源计算独立性 | 耦合在 FEM 文件中 | 保留（仅消除 Jellyrollmodel 耦合） |
| 层权重正确性 | 正确 | 正确（保持不变） |
| 单位转换 | 正确 | 正确（保持不变） |
| 新增函数 | 0 | 0 |
| 新增文件 | 0 | 0（不创建 ThermalHeatSource.jl） |
