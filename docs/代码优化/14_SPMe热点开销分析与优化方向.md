# SPMe 热点开销分析与优化方向

> 日期: 2026-04-08
> 基于: §6.2 实测数据 + 全源码审查 + 计时代码覆盖分析

---

## 1. 数据基线

### 1.1 实测计时（testexample, czm_enabled=false, nθ=80）

| 模块 | 报告占比 | 修正后估计占比 | 说明 |
|------|----------|--------------|------|
| SPMe 求解 | 91.59% | **~92-93%** | 代表性状态 SPMe_variables 未计量，被低估 ~1% |
| 热分布式模型 | 7.85% | **~7.3-7.5%** | 辅助变量循环归入 thermal，被高估 ~0.3% |
| 分流求解器 | 0.56% | ~0.55% | 基本准确 |
| CZM 模型 | 0% | 0% | 未启用 |

### 1.2 计时盲区（详见 §6.2 计时覆盖分析）

`CallModel_MultiSPMe` 中以下操作**未被任何模块计时覆盖**:

| 盲区 | 位置 | 估计占比 | 性质 |
|------|------|----------|------|
| 代表性状态 SPMe_variables | L55-57 | ~1% | SPMe 工作被漏计 |
| StandardVariables(case, 1) | L15 | ~0.1% | Dict + ~50 数组创建 |
| 状态提取 ne × copy | L21-24 | ~0.3% | 向量拷贝 |
| 面积计算循环 | L27-34 | ~0.1% | Gauss 点累加 |
| blockdiag 稀疏拼接 | L100-102 | ~0.5-1% | ne 个稀疏矩阵拼合 |
| 全局 blockdiag | L140-142 | ~0.2-0.5% | 最终矩阵拼接 |

**结论**: 分母仅含已计量部分，所有报告比率偏高 ~2-4 个百分点。SPMe 的真实占比高于报告值。

---

## 2. SPMe 单步开销构成

### 2.1 调用链

每个时间步的 SPMe 工作量：

```
CallModel_MultiSPMe (每步调用一次)
│
├─ [未计量] SPMe_variables(case, mean(yt_chem), ...)    ×1   ← 代表性状态
│
├─ [已计量] Threads.@threads for e in 1:ne              ×ne
│   └─ SPMe_element(case, yt_chem[e], t, e; ...)
│       ├─ SPMe_variables(case, yt_e, t; I_app, T_e)        ← 创建 StandardVariables Dict
│       │   ├─ StandardVariables(case, 1)                      ← ~50 数组分配
│       │   ├─ 状态提取 (for i in keys(case.index))            ← Dict 键遍历 + 拷贝
│       │   ├─ Gauss 点插值 (3 个区域)                          ← 向量化运算
│       │   ├─ Butler-Volmer (j0, eta)                         ← Arrhenius + asinh
│       │   ├─ 电解液电位 (kappa, R_EL, dphi_e)                ← 闭包调用 + IntV 积分
│       │   └─ OCV + dUdT + cell voltage                       ← 闭包调用
│       │
│       ├─ Mechanicaloutput(case, variables_e)                 ← Calstressdisp ×2
│       │   └─ Calstressdisp(NE) + Calstressdisp(PE)          ← 应力耦合
│       │
│       ├─ ElectrodeDiffusion(NE, ...)                         ← Assemble ×1 (M矩阵)
│       │   └─ Assemble(Vi, Vj, Ni, Ni, coeff, mlen)          ← 分配 KI/KJ/KV
│       ├─ ElectrodeDiffusion(PE, ...)                         ← Assemble ×1 (M矩阵)
│       │   └─ Assemble(...)
│       ├─ ElectrodeDiffusion(NE, ...)                         ← Assemble ×1 (K矩阵)
│       │   └─ Assemble(...)
│       ├─ ElectrodeDiffusion(PE, ...)                         ← Assemble ×1 (K矩阵)
│       │   └─ Assemble(...)
│       │
│       ├─ ElectrolyteDiffusion(param, mesh_el, ...)           ← Assemble ×2 (M + K)
│       │   └─ Assemble(...) ×2
│       │
│       ├─ SPMe_BC(case, variables_e)                          ← Assemble1D ×1
│       │   └─ Assemble1D(...)
│       │
│       └─ blockdiag(M_np, M_pp, M_el)                        ← 稀疏块对角
│
├─ [已计量-thermal] 辅助变量提取 (for e in 1:ne)              ← 10+ Dict 访问/单元
├─ [已计量-thermal] compute_heat_sources(...)
├─ [已计量-thermal] ThermalDistributed2D(...)
└─ [已计量-thermal] ThermalDistributed2D_BC(...)
```

### 2.2 每步操作量统计（ne=80 典型配置）

| 操作类别 | 每步次数 | 单次分配量 | 总分配次数 |
|----------|----------|-----------|-----------|
| `StandardVariables(case, 1)` | 81 | ~50 个数组 + 1 个 Dict | ~4050 个数组 |
| `Assemble` (粒子扩散) | 320 (80×4) | 3 个数组 (KI/KJ/KV) | 960 个数组 |
| `Assemble` (电解液扩散) | 160 (80×2) | 3 个数组 | 480 个数组 |
| `Assemble1D` (BC 源项) | 80 (80×1) | 1 个向量 | 80 个向量 |
| `sparse()` 转换 | 160 (80×2) | 1 个稀疏矩阵 | 160 个矩阵 |
| `vec()` 拷贝 | 80 (80×1) | 1 个向量 | 80 个向量 |
| `Dict` 创建 (SPMe_variables 返回) | 81 | 1 个 Dict + ~20 键 | 81 个 Dict |
| `blockdiag` (单元级) | 80 (80×1) | 2 个稀疏矩阵 | 160 个矩阵 |
| **合计** | | | **~5051+ 次堆分配** |

---

## 3. 开销来源分析（按严重度排序）

### 3.1 来源 A: StandardVariables 重复分配（预计占比 15-25%）

**位置**: `SPMe_variables` → `StandardVariables(case, 1)`（L122）

```julia
function SPMe_variables(case, yt, t; ...)
    variables = StandardVariables(case, 1)   # ← 每次调用创建全新 Dict
    ...
end
```

`StandardVariables`（`Variables.jl:1-141`）为 multi-SPMe 模式创建 **~50 个预分配数组**，包括：
- 基础电化学: ~14 个数组
- SPMe 电解液: ~7 个数组
- 热分布式: ~30 个数组（`thermal2D q_rxn_ne` 等，但 SPMe_element 内部不需要这些热数组）
- CZM: ~2 个数组

**问题**:
1. 每步 81 次调用 × ~50 数组 = **~4050 次堆分配**，触发 Julia GC
2. 热分布式数组（~30 个）在 SPMe_element 内部**完全无用**，白白分配
3. SPMe 模式下 `Nn=Np=1`，很多 `zeros(Float64, 1, 1)` 数组只含 1 个元素

**影响**: 分配 + GC 压力。在 Julia 中，GC 暂停会阻塞所有线程。

### 3.2 来源 B: Assemble 数组分配（预计占比 10-15%）

**位置**: `ElectrodeDiffusion.jl`, `ElectrolyteDiffusion.jl` → `Assemble()`

`Assemble`（`Assemble.jl:1-26`）每次调用分配 3 个临时数组：

```julia
function Assemble(Vi, Vj, Ni, Nj, coeff, mlen1, mlen2)
    KI = zeros(Int64, gslen * gslen1 * gslen2)    # ← 分配
    KJ = deepcopy(KI)                              # ← 再分配
    KV = zeros(Float64, gslen * gslen1 * gslen2)   # ← 分配
    ...
    K = sparse(KI, KJ, KV, mlen1, mlen2)          # ← 分配稀疏矩阵
    return K
end
```

每个 `SPMe_element` 调用 6 次 `Assemble` + 1 次 `Assemble1D`：
- 粒子扩散: 4 次 Assemble（NE M矩阵、NE K矩阵、PE M矩阵、PE K矩阵）
- 电解液扩散: 2 次 Assemble（M矩阵、K矩阵）
- BC 源项: 1 次 Assemble1D

每步总计: **80 × 6 = 480 次 Assemble**，产生 **~1440 个临时数组**。

**问题**: 粒子网格和电解液网格在所有单元之间**完全相同**（`case.mesh["negative particle"]` 对所有 e 都一样），因此 KI/KJ 的内容对每个单元也相同，仅 KV 不同。这是巨大的冗余。

### 3.3 来源 C: Dict 操作开销（预计占比 5-10%）

**位置**: `SPMe_variables`, `SPMe_element`, `Mechanicaloutput`

整个 SPMe 管线使用 `Dict{String, Union{Array{Float64}, Float64}}` 传递变量。每次调用涉及：

- `StandardVariables`: 创建 Dict + 插入 ~50 个键值对
- `SPMe_variables` L139-141: `for i in var_list; variables[i] = yt[case.index[i]]; end` — 遍历索引 Dict
- `SPMe_variables` L180-196: 写入 ~20 个结果
- `SPMe_BC` L97-98: 读取 2 个键
- `ElectrolyteDiffusion` L14-16: 读取 3 个键
- `Mechanicaloutput` L6-11: 读取 6 个键, L17-29: 写入 13 个键
- `CallModel.jl` L91-95: 每个单元将 `variables_e` 存入 `variables_elems[e]`
- `CallModel.jl` L117-127: 从 `variables_elems[e]` 读取 10+ 个键

**每步 Dict 操作估计**: 81 × (50 创建 + 20 写入 + 10 读取) ≈ **6480+ 次 Dict 操作**

**问题**: Julia 的 `Dict` 是基于哈希表的，每次查找需要计算哈希值 + 比较字符串。在每步 ne=80 × 每步多次访问的模式下，这不是零成本。

### 3.4 来源 D: 闭包函数调用（预计占比 3-8%）

**位置**: `SPMe_variables` L166-168, L177-178

```julia
kappa_ne_gs = param.EL.kappa(ce_n_gs, T) * ...   # ← 闭包调用
u_n = param.NE.U(cn_surf) .+ ...                   # ← OCV 闭包调用
... .* param.NE.dUdT(cn_surf)                      # ← dUdT 闭包调用
```

这些是通过 `Base.invokelatest` 包装的闭包（解决 Julia world-age 问题），每次调用有额外分发开销。在 multi-SPMe 模式下每步调用 81 × ~6 次 ≈ **486+ 次闭包调用**。

### 3.5 来源 E: 稀疏矩阵操作（预计占比 3-5%）

**位置**: `CallModel.jl` L92-94, L100-102, L140-142

```julia
# 每个单元内部
M_elems[e] = sparse(M_e)          # blockdiag(SparseMatrix) → sparse
K_elems[e] = sparse(K_e)

# 单元完成后
M_chem = blockdiag(M_elems...)    # 80 个稀疏矩阵 blockdiag
K_chem = blockdiag(K_elems...)

# 全局拼装
M = blockdiag(M_chem, sparse(MT)) # 最终 blockdiag
K = blockdiag(K_chem, sparse(KT))
```

`blockdiag` 对稀疏矩阵需要分配新的合并数组。每步总计:
- 80 × 2 = 160 次 `sparse()` 转换
- 1 + 1 = 2 次 80-矩阵 `blockdiag`
- 1 + 1 = 2 次全局 `blockdiag`

### 3.6 来源 F: 代表性状态的冗余计算（预计占比 ~1%）

**位置**: `CallModel.jl` L55-63

```julia
yt_representative = mean(yt_chem)     # ne 个向量求均值 → 新向量
T_rep = mean(Te_prev)                 # ne 个标量求均值
vars_rep = SPMe_variables(case, yt_representative, t; I_app=I_total, T_e=T_rep)
```

这调用了一次**完整的 SPMe_variables**（含 StandardVariables 分配），仅为了给分流求解器提供初始电压估计。分流求解器本身仅占 0.56%，代表状态计算与其不成比例。

---

## 4. 开销占比估计

| 来源 | 估计占 SPMe 总时间 | 占总步时间 | 可优化性 |
|------|-------------------|-----------|---------|
| A: StandardVariables 分配 | 15-25% | **14-23%** | 高 — 预分配复用 |
| B: Assemble 临时数组 | 10-15% | **9-14%** | 中 — 预分配 KI/KJ |
| C: Dict 操作 | 5-10% | **5-9%** | 低 — 需改架构 |
| D: 闭包调用 | 3-8% | **3-7%** | 低 — Julia 语言限制 |
| E: 稀疏矩阵操作 | 3-5% | **3-5%** | 中 — 减少转换 |
| F: 代表性状态冗余 | ~1% | **~1%** | 高 — 可简化 |
| 物理计算 (BV/OCV/扩散) | 剩余 ~40-60% | **不可压缩** | — |

---

## 5. 优化方向

### 方向 1: StandardVariables 预分配复用（收益: 14-23%）

**核心思路**: 不再每次 `SPMe_element` 调用都创建新 Dict，而是预分配一份（或每线程一份），每次调用仅清零重用。

**方案 A — 对象池**（推荐，低风险）:
```julia
# CallModel_MultiSPMe 入口，每次调用创建一次
nthreads = Threads.nthreads()
var_pool = [StandardVariables(case, 1) for _ in 1:nthreads]

# 并行循环中
Threads.@threads for e in 1:ne
    tid = Threads.threadid()
    vars_e = var_pool[tid]
    # 清零所有数组
    for (k, v) in vars_e
        isa(v, Array{Float64}) && fill!(v, 0.0)
    end
    # ... 使用 vars_e ...
end
```

优点: 不修改 `SPMe_variables` 签名，兼容性好。
缺点: 需要确保 `SPMe_variables` 是覆盖写入（而非追加），需验证。

**方案 B — 精简 element-only 变量集**:
创建 `create_element_workspace(case)` 只包含 SPMe_element 实际需要的 ~20 个键（排除 ~30 个 thermal2D 键），减少无用分配。

优点: 减少 ~60% 的数组分配量。
缺点: 需修改 `SPMe_variables` 内部逻辑或新增变体。

### 方向 2: Assemble 预分配 KI/KJ（收益: 5-9%）

**核心思路**: 粒子网格的 KI/KJ 对所有单元相同，只预计算一次；KV 每次不同但可预分配缓冲区复用。

**观察**: `ElectrodeDiffusion` 对所有单元使用**相同的网格** (`case.mesh["negative particle"]`)，因此:
- KI（行索引）和 KJ（列索引）对所有单元**完全相同**
- 仅 coeff 和 KV 随单元变化

**方案**: 在 `SPMe_element` 外部预计算 KI/KJ 模板，内部只更新 KV：

```julia
# 预计算一次（Initialisation 或首次调用时）
KI_np_template, KJ_np_template = precompute_assemble_indices(mesh_np, mlen_np)
KI_pp_template, KJ_pp_template = precompute_assemble_indices(mesh_pp, mlen_pp)
KI_el_template, KJ_el_template = precompute_assemble_indices(mesh_el, mlen_el)

# SPMe_element 内部
KV_buf = KV_buffer[tid]  # 预分配
KV_buf .= Ds_eff .* mesh.gs.x.^2 .* mesh.gs.weight .* mesh.gs.detJ
K_np = sparse(KI_np_template, KJ_np_template, KV_buf, mlen, mlen)
```

优点: 消除 480 次 `zeros(Int64, ...)` 和 480 次 `deepcopy(KI)`。
缺点: 需新增 `precompute_assemble_indices` 和 `Assemble!` 变体。

### 方向 3: 减少稀疏矩阵转换（收益: 3-5%）

**核心思路**: `SPMe_element` 内部的 `blockdiag(M_np, M_pp, M_el)` 产生 DenseMatrix，外层再 `sparse()`。如果直接在稀疏格式工作，可消除中间步骤。

**观察**: `ElectrodeDiffusion` 和 `ElectrolyteDiffusion` 的输出已经是 `SparseMatrixCSC`（通过 `Assemble` 的 `sparse(KI,KJ,KV,...)` 构建），`blockdiag` 对稀疏矩阵也输出稀疏矩阵。因此 `sparse(M_e)` 实际上是 `sparse(SparseMatrixCSC)` — 一个冗余的类型转换。

**方案**: 移除 `CallModel.jl` L92-93 的 `sparse()` 调用：
```julia
# OLD:
M_elems[e] = sparse(M_e)
K_elems[e] = sparse(K_e)

# NEW: M_e, K_e 已经是 SparseMatrixCSC，直接存入
M_elems[e] = M_e
K_elems[e] = K_e
```

### 方向 4: 简化代表性状态计算（收益: ~1%）

**核心思路**: 用已求得的 `variables_elems[1]` 或均值直接构造 representative variables，避免完整 `SPMe_variables` 调用。

**方案**: 在并行循环完成后，直接从 `variables_elems` 合并：
```julia
# OLD:
yt_representative = mean(yt_chem)
vars_rep = SPMe_variables(case, yt_representative, t; I_app=I_total, T_e=T_rep)

# NEW: 直接用第一个单元结果，或用分流求解的 Vc
for (k, v) in variables_elems[1]
    variables[k] = v
end
```

需验证: 分流求解器 `compute_prefactors` 依赖 `variables` 中的浓度和 OCV 值，使用 element 1 的值是否满足精度要求。

### 方向 5: 减少闭包调用开销（收益: 3-7%，长期方向）

**核心思路**: `Base.invokelatest` 包装的 `U`、`dUdT`、`kappa`、`De` 闭包有分发开销。对固定参数文件，可缓存为直接函数引用。

**方案**: 在 `NormaliseParam` 阶段将 `Base.invokelatest(fn, args...)` 替换为 `invokelatest` 仅在首次调用时执行，后续缓存结果函数。

这是 Julia 语言的限制，优化空间有限，长期可考虑参数文件编译为原生函数。

### 方向 6: 减少 Dict 操作（收益: 5-9%，长期方向）

**核心思路**: 将 `Dict{String, ...}` 替换为类型化结构体，消除哈希查找和字符串比较。

**方案**: 新增 `SPMeElementVariables` struct:
```julia
struct SPMeElementVariables
    cn_surf::Float64
    cp_surf::Float64
    eta_n::Float64
    eta_p::Float64
    j0_n::Float64
    j0_p::Float64
    V_cell::Float64
    j_n::Float64
    j_p::Float64
    # ... ~20 个字段
end
```

这是**最大的架构变更**，需要修改 `SPMe_variables`、`SPMe_element`、`SPMe_BC`、`ElectrolyteDiffusion`、`Mechanicaloutput`、`compute_heat_sources` 等多个文件的接口。建议作为**独立重构阶段**，在上述低风险优化完成后再考虑。

---

## 6. 优化路线图

```
阶段 1 — 即时见效（1-2 天）                  预期收益
─────────────────────────────────────────────────────────
├─ deepcopy → copy (Solve.jl + Assemble.jl)     ~3%
├─ 缓存单元面积 (MultiSPMeLayout.areas)         ~0.5%
└─ @views 替代 extract_element_state             ~1%
                                                 小计: ~4.5%

阶段 2 — 预分配优化（3-5 天）                  预期收益
─────────────────────────────────────────────────────────
├─ StandardVariables 对象池 (方向 1 方案 A)       ~15-23%
├─ 移除冗余 sparse() 转换 (方向 3)               ~3-5%
└─ 简化代表性状态 (方向 4)                       ~1%
                                                 小计: ~19-29%

阶段 3 — Assemble 预分配（2-3 天）             预期收益
─────────────────────────────────────────────────────────
├─ 预计算 KI/KJ 索引模板 (方向 2)                ~5-9%
└─ 预分配 KV 缓冲区                              ~2-3%
                                                 小计: ~7-12%

阶段 4 — 架构重构（长期，可选）                预期收益
─────────────────────────────────────────────────────────
├─ Dict → 类型化 Struct (方向 6)                  ~5-9%
└─ 闭包优化 (方向 5)                             ~3-7%
                                                 小计: ~8-16%

总计潜在收益: ~40-60% SPMe 时间减少
换算为总步时间: ~37-55% wall-clock 减少
```

---

## 7. 与现有优化计划的关系

| 本文档方向 | 对应 `2026-04-07-simulation-speedup.md` Task | 状态 |
|-----------|----------------------------------------------|------|
| 方向 1 (StandardVariables 复用) | Task 7 (SPMe 工作区) | 已规划，高风险 |
| 方向 2 (Assemble 预分配) | Task 8 (预分配 Assemble) | 已规划，P3 |
| 方向 3 (sparse 转换) | 未覆盖 | **新增** |
| 方向 4 (代表性状态) | ~~Task 5~~ (已否决简化版) | 可选 |
| 方向 5 (闭包优化) | 未覆盖 | 长期方向 |
| 方向 6 (Dict → Struct) | 未覆盖 | 长期方向 |

### 新增发现（未在 speedup 计划中）

1. **`sparse()` 冗余转换**（方向 3）: `blockdiag` 对 `SparseMatrixCSC` 输入已返回 `SparseMatrixCSC`，外层 `sparse()` 是无操作。这行代码可直接删除。
2. **StandardVariables 分配了 ~30 个无用 thermal 数组**: `StandardVariables` 为 `distributed2D` 模式创建了 30 个热相关数组（`thermal2D q_rxn_ne` 等），但 `SPMe_element` 内部从不使用这些键。仅这一项就浪费了 ~60% 的 StandardVariables 分配量。
3. **Assemble 的 KI/KJ 在同网格单元间完全重复**: 80 个单元共享相同的粒子网格和电解液网格，因此 KI/KJ 索引数组重复计算 80 次。这是最大的"重复计算"来源。
