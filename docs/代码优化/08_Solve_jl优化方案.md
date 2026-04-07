# Solve.jl 优化方案

> 日期: 2026-04-01 (修订)
> 文件: `src/Solve.jl`
> 状态: 修改 (155→717 行, +362%)
> main 分支行数: 155

---

## 1. main 分支现状

### 1.1 函数清单

| 函数 | 行数 | 职责 |
|------|------|------|
| `Solve(case)` | ~100 | 纯时间步进器 |
| `CallModel(case, yt, t; jacobi)` | ~25 | 薄调度器 (SPM/SPMe/P2D/sP2D) |
| `RecordMatrix!(case, M, K)` | 10 | 稀疏矩阵记录 |
| `ErrorEstimation(case, y_old, y_new, coeff)` | 20 | dt 自适应 |

### 1.2 核心设计模式

```
Solve: 纯时间步进器
  ├── 初始化 y0, dt, theta
  ├── while t <= t_end
  │     CallModel → (M, K, F, variables, y_phi)
  │     时间离散: y_new = Mt \ (Kt*y_old + Ft)
  │     ErrorEstimation → dt自适应
  │     Variable_update!
  │     电压截止 → break
  └── PostProcessing
```

**main 分支 Solve 约 100 行，架构清晰。**

---

## 2. 当前分支膨胀根因

| 新增功能 | 塞入位置 | 行数 | 应在位置 |
|----------|----------|------|----------|
| 文件日志重定向 | Solve 开头 | ~35 | 保留在 Solve（仅精简） |
| MultiSPMe 初始化分支 | Solve 中部 | ~40 | 保留在 Solve |
| distributed2D 热场初始化 | Solve 中部 | ~40 | 保留在 Solve |
| 单元截止电压检测 | Solve 循环内 | ~60 | 保留在 Solve |
| 物理单位后处理 | Solve 尾部 | ~60 | 保留在 Solve |
| `CallModel_MultiSPMe` | Solve 下半部 | 148 | **迁出到 CallModel.jl** |
| `CallModel` 扩展 | Solve 下半部 | 109 | **迁出到 CallModel.jl** |

**核心策略**: 将 `CallModel` 和 `CallModel_MultiSPMe` 搬到 `CallModel.jl`，Solve.jl 仅保留步进器 + 初始化 + 后处理。所有逻辑保持内联，不提取子函数。

---

## 3. 优化方案

### 3.1 拆分策略

```
Solve.jl (717 → ~480 行)           迁出到 CallModel.jl:
──────────────────────────────     ──────────────────
Solve()           ~430行(保留)
RecordMatrix!      ~10行(保留)
ErrorEstimation    ~20行(保留)
                                   CallModel_MultiSPMe  148行(整体搬家)
                                   CallModel            109行(整体搬家)
```

### 3.2 Solve() 保持原结构

Solve() 主函数不拆分子函数。所有初始化、时间步进、截止检测、后处理逻辑保持内联。

仅做以下修改：
1. `case.multi_spme_layout["key"]` → `case.layout.key` 类型替换
2. `case.param.cell.T0` 直接修改 bug → 用局部变量
3. `czm_mesh` 作用域问题 → 在函数开头声明为 `nothing`

### 3.3 Bug 修复（在 Solve 内部）

#### 3.3.1 `czm_mesh` 作用域 (行 531 附近)

```julia
# 旧: czm_mesh 可能在某些分支未定义
# 新: 在 Solve() 函数开头声明
czm_mesh = nothing
czm_params = nothing

# 后续分支中赋值:
if case.opt.czm_enabled && case.czm_mesh !== nothing
    czm_mesh = case.czm_mesh
    czm_params = ...
end
```

#### 3.3.2 `case.param.cell.T0` 副作用 (行 210 附近)

```julia
# 旧:
case.param.cell.T0 = Tm  # 直接修改参数！

# 新:
T_ref_current = Tm  # 使用局部变量
# 后续所有引用 case.param.cell.T0 的地方改为 T_ref_current
```

---

## 4. CallModel.jl (新文件)

### 4.1 内容

从 Solve.jl 整体搬出 `CallModel_MultiSPMe`（148行）和 `CallModel`（109行），**不做任何逻辑拆分**，仅做：

1. 函数重命名: `CallModel_MultiSPMe` → `call_model_multi_spme`, `CallModel` → `call_model`
2. 类型替换: `case.multi_spme_layout["key"]` → `case.layout.key`
3. 修复热源计算中对 `jellyroll_element_properties` 的直接调用 → 通过 `case.geometry` 传递

### 4.2 热源→FT 映射（关键修复）

当前 `call_model_multi_spme` 中 FT 映射逻辑已存在但可能不完整。检查并确保：

```julia
# 热源 (逐单元) 映射到 F 向量 (逐节点)
# 需要确认映射逻辑正确，不能简单 FT .+= q_fields
```

---

## 5. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| Solve.jl 行数 | 717 | ~480 |
| CallModel.jl 行数 | 0 (在 Solve 内) | ~260 (整体搬家) |
| `Solve()` 主函数行数 | ~430 | ~430（仅修 bug + 类型替换） |
| Bug 修复 | 0 | 2 (czm_mesh 作用域 + T0 副作用) |
| 新增函数 | 0 | 0（仅搬家 + 改名） |
| main 分支 SPM/SPMe 标准路径 | N/A | 零改动 |
| `RecordMatrix!` | 不动 | 不动 |
| `ErrorEstimation` | 不动 | 不动 |
