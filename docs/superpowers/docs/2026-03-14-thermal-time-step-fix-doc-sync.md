# 热模型时间步修复与文档同步

**日期**: 2026-03-14  
**原因**: `Solve.jl` 中热模型分支的时间步处理方式修改，需要同步更新相关技术文档

---

## 修复概述

### 核心变化

**文件**: `src/Solve.jl` 第38-42行

```julia
if case.opt.model == "thermal"
    t0 = case.opt.time[1]/case.param.scale.t0      # 无量纲初始时间
    t_end = case.opt.time[end]/case.param.scale.t0   # 无量纲终止时间
    dt = case.opt.dt[1]/case.param.scale.t0        # 无量纲时间步长
    vars = case.multi_spme_layout["thermal_variables"]
    update_fn = case.multi_spme_layout["thermal_update_fn"]
    # ... 后续求解使用无量纲时间
end
```

### 关键特性

1. **自动时间归一化**：
   - 输入的物理时间 (`opt.time`, `opt.dt`) 自动除以 `param.scale.t0` 转换为无量纲时间
   - 求解器内部所有时间变量均为无量纲

2. **与电化学模型一致性**：
   - 热模型专用分支 (`model == "thermal"`) 的时间处理方式与电化学-热耦合分支 (`model == "SPMe"`) 保持一致
   - 统一使用 `scale.t0` 进行时间归一化

3. **结果时间转换**：
   - 求解结果的 `ring_sol.time` 为无量纲时间
   - 需要乘以 `scale.t0` 还原为物理时间

---

## 验证脚本修改

### 文件: `example/热模块验证/jellyroll_vs_ring_thermal_compare.jl`

#### 修改1: 热源更新函数 (第249-254行)

**修改前**:
```julia
t_vec = collect(t)
function update_heat_source(t_phys, vars)
    q_at_t = [interp_linear(t_vec, q_ring_nd[e, :], [t_phys])[1] for e in 1:ne]
    vars["heat_source_fields"] = q_at_t
end
```

**修改后**:
```julia
t_vec_nd = collect(t) ./ param_ring.scale.t0  # 转换为无量纲时间
function update_heat_source(t_nd, vars)
    q_at_t = [interp_linear(t_vec_nd, q_ring_nd[e, :], [t_nd])[1] for e in 1:ne]
    vars["heat_source_fields"] = q_at_t
end
```

**原因**: `update_heat_source` 现在接收无量纲时间参数（由 `Solve.jl` 自动归一化）

#### 修改2: 求解后时间转换 (第271-272行)

**修改前**:
```julia
ring_sol = JuBat.Solve(case_ring)
A_elem_ring = compute_element_areas(mesh_ring)
T_vol_ring = [area_weighted_mean(compute_element_mean(mesh_ring, ring_sol.T_hist[:, k] .* T_ref), A_elem_ring)
              for k in 1:length(ring_sol.time)]
```

**修改后**:
```julia
ring_sol = JuBat.Solve(case_ring)

# ring_sol.time 现在是无量纲时间，需要转换回物理时间
t_ring_phys = ring_sol.time .* param_ring.scale.t0
A_elem_ring = compute_element_areas(mesh_ring)
T_vol_ring = [area_weighted_mean(compute_element_mean(mesh_ring, ring_sol.T_hist[:, k] .* T_ref), A_elem_ring)
              for k in 1:length(t_ring_phys)]
```

**原因**: `ring_sol.time` 现在返回无量纲时间，需要转换回物理时间以保持一致性

#### 修改3: 结果返回 (第287行)

**修改前**:
```julia
return (; t=ring_sol.time, T_vol=T_vol_ring, T_hist=ring_sol.T_hist, mesh=mesh_ring,
```

**修改后**:
```julia
return (; t=t_ring_phys, T_vol=T_vol_ring, T_hist=ring_sol.T_hist, mesh=mesh_ring,
```

**原因**: 返回物理时间而不是无量纲时间，保持与 Jellyroll 结果的一致性

#### 修改4: 调试日志更新 (第265行, 第283行)

更新调试日志以记录无量纲时间和物理时间的转换：
- 添加 `dt_nd = opt.dt[1]/param_ring.scale.t0`
- 添加 `t_vec_nd_first` 和 `t_vec_physical_first` 记录
- 添加 `t_ring_nd_first` 和 `t_ring_phys_first` 记录

---

## 技术文档同步

### 文档1: `md/05_热模型_二维分布式.md`

**位置**: 第112-120行 (时间离散章节)

**修改内容**:
- 更新时间离散公式，移除独立的 `t_{th}` 尺度
- 明确统一时间尺度 `t_0 = 3600` s
- 添加实现细节代码示例
- 添加验证脚本使用说明

**关键更新**:
```julia
### 3.2 时间离散

**后退欧拉格式**：

$$
\left( \frac{M_T}{\Delta t^*} + K_T \right) T^{n+1} = \frac{M_T}{\Delta t^*} T^n + F_T + F_{boundary}
$$

**统一时间尺度**（2026-03-14 更新）：
- 使用统一时间尺度 \(t_0 = 3600\) s
- 热模型时间步长 \(\Delta t^* = \Delta t / t_0\)（无量纲）
- 与电化学模型使用相同时间尺度，无需额外缩放因子
```

### 文档2: `md/01_参数定义与归一化.md`

**位置**: 第542-558行 (时间步处理章节)

**修改内容**:
- 添加热模型专用分支的说明
- 添加 `Solve.jl` 实现代码示例
- 添加验证脚本使用说明
- 扩展优势列表，添加接口统一性

**关键更新**:
```julia
### 5.3 时间步处理

**热模型专用分支（`Solve.jl`）**：

当 `opt.model == "thermal"` 时，求解器会自动进行时间归一化：

```julia
# src/Solve.jl 第38-42行
if case.opt.model == "thermal"
    t0 = case.opt.time[1]/case.param.scale.t0      # 无量纲初始时间
    t_end = case.opt.time[end]/case.param.scale.t0   # 无量纲终止时间
    dt = case.opt.dt[1]/case.param.scale.t0        # 无量纲时间步长
    # ... 后续使用无量纲时间求解
end
```

**验证脚本使用说明**：
```julia
# 输入时间需要转换为无量纲
t_vec_nd = collect(t_phys) ./ param.scale.t0  # 物理时间 → 无量纲时间

# 热源更新函数接收无量纲时间
function update_heat_source(t_nd, vars)
    # t_nd 是无量纲时间
    vars["heat_source_fields"] = interpolate_heat_source(t_nd)
end

# 求解结果的时间需转换回物理时间
ring_sol = JuBat.Solve(case)
t_phys = ring_sol.time .* param.scale.t0  # 无量纲时间 → 物理时间
```
```

### 文档3: `md/10_参数传递与模块架构.md`

**位置**: 第450行后 (新增章节)

**修改内容**:
- 新增第6章：求解器分支说明
- 详细说明热模型专用分支的特性
- 提供验证脚本使用示例
- 对比热模型分支与电化学-热耦合分支的区别

**关键更新**:
```julia
## 6. 求解器分支说明

### 6.1 热模型专用分支 (`Solve.jl`)

当 `opt.model == "thermal"` 时，求解器会进入热模型专用分支，该分支有特殊的时间处理方式：

**关键特性**：

1. **自动时间归一化**：
   - 输入的物理时间自动除以 `param.scale.t0` 转换为无量纲时间
   - 求解器内部所有时间变量均为无量纲

2. **热源更新函数接口**：
   - `thermal_update_fn` 接收无量纲时间参数
   - 验证脚本需要先将物理时间转换为无量纲时间

3. **结果时间转换**：
   - 求解结果的 `ring_sol.time` 为无量纲时间
   - 需要乘以 `scale.t0` 还原为物理时间
```

### 文档4: `md/12_热模型验证方案.md`

**位置**: 第97-123行 (测试脚本章节)

**修改内容**:
- 添加第5.2节：使用热模型专用分支
- 提供完整的验证脚本示例
- 说明时间尺度转换的关键要点

**关键更新**:
```julia
### 5.2 使用热模型专用分支

当使用 `model == "thermal"` 进行验证时，需要注意时间尺度的转换：

```julia
# example/热模块验证/jellyroll_vs_ring_thermal_compare.jl

opt.model = "thermal"  # 使用热模型专用分支
opt.thermalmodel = "ring2D_polar"
opt.time = [0.0, 60.0]  # 物理时间 [s]
opt.dt = [0.2, 2.0]      # 物理时间步 [s]

# 预先将物理时间转换为无量纲时间
t_vec_nd = collect(t_phys) ./ param_ring.scale.t0

function update_heat_source(t_nd, vars)
    # t_nd 是无量纲时间（Solve.jl 会自动除以 t0）
    q_at_t = [interp_linear(t_vec_nd, q_ring_nd[e, :], [t_nd])[1] for e in 1:ne]
    vars["heat_source_fields"] = q_at_t
end
```
```

---

## 修改总结

### 核心原则

1. **统一时间尺度**：热模型专用分支与电化学-热耦合分支使用相同的时间归一化方式
2. **自动转换**：`Solve.jl` 自动进行时间归一化，验证脚本无需手动处理
3. **一致性**：结果时间变量统一为无量纲，输出时需转换回物理时间

### 影响范围

| 组件 | 修改类型 | 影响 |
|------|----------|------|
| `src/Solve.jl` | 核心修复 | 热模型分支时间处理方式 |
| `example/热模块验证/jellyroll_vs_ring_thermal_compare.jl` | 同步修改 | 验证脚本适配新时间处理 |
| `md/05_热模型_二维分布式.md` | 文档更新 | 时间离散章节 |
| `md/01_参数定义与归一化.md` | 文档更新 | 时间步处理章节 |
| `md/10_参数传递与模块架构.md` | 文档更新 | 新增求解器分支说明 |
| `md/12_热模型验证方案.md` | 文档更新 | 测试脚本章节 |

### 验收标准

- [x] `Solve.jl` 热模型分支正确进行时间归一化
- [x] 验证脚本适配新的时间处理方式
- [x] 技术文档同步更新，反映最新实现
- [x] 验证脚本运行正常，温度曲线合理
- [x] Jellyroll 与 Ring 模型温度演化一致性良好

---

## 参考文档

- `md/01_参数定义与归一化.md` - 时间尺度和参数归一化
- `md/05_热模型_二维分布式.md` - 热模型时间离散
- `md/10_参数传递与模块架构.md` - 求解器分支说明
- `md/12_热模型验证方案.md` - 热模型验证方法
- `example/热模块验证/jellyroll_vs_ring_thermal_compare.jl` - 验证脚本
