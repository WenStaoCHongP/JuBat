# SPMe-分布式热-CZM 仿真提速优化 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在浮点结果 bit 级不变的前提下，通过缓存、预分配、消除重复计算，降低 `example/testexample.jl`（SPMe-分布式热-CZM 全耦合基线场景）的仿真墙钟时间。

**Architecture:** 测量驱动分批优化。批次 0（Task 1-2）建立本机基线锚点并用探针脚本验证两个关键事实：(a) 电化学+热全局矩阵 M/K 是否跨步常量；(b) 各模块实测耗时占比。后续任务（Task 3-7）各带执行门（gate），由批次 0 写入 `findings.md` 的结论裁决是否执行与执行顺序。每个优化批次独立 commit、独立基线验收（四判据 bit 一致）。

**Tech Stack:** Julia 1.11.2（`D:\Julia-1.11.2\bin\julia.exe`）、SparseArrays/SuiteSparse UMFPACK、Profile 标准库。不引入任何新依赖。

**Spec:** `docs/superpowers/specs/2026-08-19-perf-optimization-spme-thermal-czm-design.md`（本计划从该 spec 出发，执行者须先读 spec）

## Global Constraints

- 基线运行环境（每个基线运行都必须用）：
  ```bash
  cd /c/Users/19303/Desktop/github/JuBat
  GKSwstype=100 JULIA_NUM_THREADS=1 /d/Julia-1.11.2/bin/julia.exe --startup-file=no example/testexample.jl
  ```
- 基线四判据（每批改动前后各跑一次本机基线，全部一致才通过）：退出码 0；网格/步数（控制台输出的单元数 ne、节点数 nT、`timing CallModel calls`）；`metrics.toml` 记录的科学结果（对照 `Simplify/baseline/testexample/metrics.toml`，按其记录精度）；本机锚点 PNG SHA-256（Task 1 建立；若 Task 1 判定 SHA 本机不稳定，则改为前三判据）。
- **数值约束（最高优先级）**：任何改动不得改变浮点运算的数值与顺序。只允许：缓存已证明跨步不变的对象、预分配替换新建（数值路径逐位一致）、跳过结果未被使用的计算。禁止：预聚合多个累加步（会改变求和顺序）、换求解器、改收敛容差。
- 每个任务单独 commit；四判据任一不一致 → `git revert` 该批 commit，在 findings.md 记录后跳过该手段。
- 规划文件目录：`docs/planning-with-files/29_仿真提速-SPMe热CZM/`（findings.md / task_plan.md / progress.md 三件套按 AGENTS.md 9.5 管理；`task_plan.md` 由执行者从本计划摘抄任务清单生成）。
- 所有 findings、commit message 用中文或中英混合，与仓库既有风格一致。
- 任务执行门：Task 3-7 开工前必须读 `docs/planning-with-files/29_仿真提速-SPMe热CZM/findings.md` 中 Task 2 写入的门裁决。Task 3-7 之间除标注的依赖外相互独立，**实际执行顺序按 findings.md 实测占比从高到低调整**（依赖：Task 6 依赖 Task 4 已合并；Task 7 依赖全部前置任务结束）。

---

### Task 1: 批次 0a — 本机基线锚点与 PNG SHA 稳定性验证

**Files:**
- Create: `docs/planning-with-files/29_仿真提速-SPMe热CZM/findings.md`
- Create: `docs/planning-with-files/29_仿真提速-SPMe热CZM/task_plan.md`
- Create: `docs/planning-with-files/29_仿真提速-SPMe热CZM/progress.md`
- Modify: `docs/planning-with-files/index.md`

**Interfaces:**
- Consumes: 无（首任务）
- Produces: `锚点-1/锚点-2` 两次基线运行记录（墙钟、四模块 timing 分解、PNG SHA 清单、退出码），写入 findings.md §锚点；PNG SHA 是否本机稳定的结论（findings.md §判据裁决）

- [ ] **Step 1: 运行基线第一次并记录**

```bash
cd /c/Users/19303/Desktop/github/JuBat
mkdir -p docs/planning-with-files/仿真提速-SPMe热CZM
{ time GKSwstype=100 JULIA_NUM_THREADS=1 /d/Julia-1.11.2/bin/julia.exe --startup-file=no example/testexample.jl ; } > /tmp/perf_anchor1.log 2>&1
echo "exit=$?"
grep -E "总单元数|总节点数|timing|SPMe 求解|分流求解器|热分布式模型|CZM 模型|CallModel" /tmp/perf_anchor1.log
sha256sum output/*.png 2>/dev/null
```

把以下内容记入 findings.md §锚点-1：`time` 的 real 值（墙钟）、四模块耗时分解行（SPMe 求解/分流求解器/热分布式模型/CZM 模型的 total 与 ratio%）、`timing CallModel calls` 值、全部 output/*.png 的 SHA-256 清单、退出码。

- [ ] **Step 2: 运行基线第二次并对比 SHA 稳定性**

```bash
{ time GKSwstype=100 JULIA_NUM_THREADS=1 /d/Julia-1.11.2/bin/julia.exe --startup-file=no example/testexample.jl ; } > /tmp/perf_anchor2.log 2>&1
sha256sum output/*.png 2>/dev/null
diff <(sha256sum output/*.png) <(grep -oE '^[0-9a-f]{64} ' /dev/null) 2>/dev/null; # 手动对比锚点-1与锚点-2的SHA清单
```

findings.md §判据裁决 写入结论：两次运行 PNG SHA 清单完全一致 → 后续批次以锚点-2 的 SHA 清单为对比锚点；不一致 → 后续批次放弃 PNG 判据，仅用退出码+网格/步数+metrics.toml 三判据（并在 findings 注明）。

- [ ] **Step 3: 初始化规划三件套并提交**

`task_plan.md` 内容：从本计划 Global Constraints 之后的任务清单抄录 Task 1-8 标题与门条件。`progress.md` 内容：`- 2026-08-19 Task 1 完成：锚点建立，SHA 稳定性=<结论>`。findings.md 结构：

```markdown
# 仿真提速-SPMe热CZM findings

## §锚点-1 / §锚点-2
（运行记录）

## §判据裁决
PNG SHA 本机稳定：<是/否>；后续四判据采用：<四判据/三判据>

## §实测占比（Task 2 填写）
## §事实核查（Task 2 填写）
## §批次裁决（Task 2 填写）
## §批次记录（各优化任务追加）
```

在 `docs/planning-with-files/index.md` 按其表格格式追加一行任务记录（说明、时间、Git 修改次数 0、跟踪状态 进行中）。

```bash
git add docs/planning-with-files/仿真提速-SPMe热CZM docs/planning-with-files/index.md
git commit -m "perf(批次0a): 建立本机基线锚点并验证PNG SHA稳定性"
```

---

### Task 2: 批次 0b — 矩阵常量性探针 + CZM 残留核对 + 批次裁决

**Files:**
- Create: `docs/planning-with-files/29_仿真提速-SPMe热CZM/probe_matrix_constancy.jl`
- Create: `docs/planning-with-files/29_仿真提速-SPMe热CZM/probe_equivalence.jl`
- Modify: `docs/planning-with-files/29_仿真提速-SPMe热CZM/findings.md`

**Interfaces:**
- Consumes: Task 1 的 findings.md 结构
- Produces: findings.md §实测占比（四模块 ratio%）、§事实核查（`M/K 跨步常量 = true/false`、CZM 单次更新迭代数）、§批次裁决（Task 3-7 每个任务 开/关 + 执行顺序）。后续所有任务的门从这里读。

- [ ] **Step 1: 写矩阵常量性探针脚本**

`probe_matrix_constancy.jl`（case 构造块逐行复制自 `example/testexample.jl:29-86`，保持完全一致的参数）：

```julia
# 探针：验证 CallModel 返回的全局 M/K 是否跨状态常量；测量 CZM 单次更新迭代数
include(joinpath(@__DIR__, "../../src/JuBat.jl"))
using .JuBat

param_dim = JuBat.ChooseCell("Jellyroll")
param_dim.cell.v_l = 2.5
param_dim.cell.v_h = 4.2
opt = JuBat.Option()
opt.Current = x -> 5.0
opt.model = "SPMe"
opt.Nn = 10; opt.Ns = 5; opt.Np = 10
opt.Nrn = 10; opt.Nrp = 10
opt.gsorder = 2
opt.dimension = 1
opt.mechanicalmodel = "none"
opt.time = [0.0, 60]
opt.dt = [0.5, 10]
opt.dtType = "auto"
opt.jacobi = "update"
opt.solveType = "Crank-Nicolson"
opt.thermal_enabled = true
opt.thermalmodel = "distributed2D"
opt.thermal_dim = "2D"
opt.cool_method = "surface"
opt.per_element_spme = true
opt.czm_enabled = true
opt.czm_fix_inner = false
opt.czm_iter_method = "basic"
opt.czm_load_steps = 10
opt.czm_tol = 1e-3

case = JuBat.SetCase(param_dim, opt)
mesh_data = JuBat.jellyroll_collector_seed_mesh(case.param; nθ=80, czm_enabled=true, gsorder=2)
case = JuBat.setup_thermal2D_mesh(case, mesh_data)
case.czm_mesh = JuBat.create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)

t = 0.0
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
M1, K1, F1, vars1, _ = JuBat.CallModel(case, y0, t; jacobi="update")

# 状态扰动：浓度类 DOF ±10%、温度 DOF +1%（DOF 总长 = ne*n_chem + nT）
y0b = copy(y0)
y0b .*= 1.05
M2, K2, F2, vars2, _ = JuBat.CallModel(case, y0b, t; jacobi="update")

println("=== 矩阵常量性 ===")
println("M identical: ", M1 == M2)
println("K identical: ", K1 == K2)
println("F identical (预期 false): ", F1 == F2)

println("=== CZM 单次更新 ===")
T_nodes = JuBat.get_thermal_dofs(y0, case.layout)
res = JuBat.update_czm_damage!(case, vars1, T_nodes)
println("iterations = ", res.iterations, ", converged = ", res.converged)
```

注意：`update_czm_damage!` 与 `get_thermal_dofs` 若未从 JuBat 导出，用 `JuBat.<名>` 限定访问即可（模块内函数均可限定访问）；若 `ModelInitialisation_MultiSPMe` 名字报错，先 `grep -rn "function ModelInitialisation" src/` 确认实际名。

- [ ] **Step 2: 运行探针并记录事实 3a/3b**

```bash
GKSwstype=100 JULIA_NUM_THREADS=1 /d/Julia-1.11.2/bin/julia.exe --startup-file=no docs/planning-with-files/29_仿真提速-SPMe热CZM/probe_matrix_constancy.jl 2>&1 | tee /tmp/probe.log
```

把 `M identical` / `K identical` / CZM `iterations` 写入 findings.md §事实核查。
- 若 `M identical && K identical` 为 true → 事实 3a 成立（矩阵常量），Task 3/4/6 门全开。
- 若 false → 找出变化块（在探针里追加打印 `M1 .!= M2` 的非零位置分布，确定是化学块还是热块），Task 3 只缓存相等的块，Task 4/6 关闭或降级。

- [ ] **Step 3: 写分解缓存等价性实验（Task 6 的前置门）**

`probe_equivalence.jl`（无占位，可直接运行）：

```julia
# 等价性 1：稀疏方阵 lu(A)\b 与 A\b 是否逐位一致（决定 Task 6 分解因子缓存可行性）
# 等价性 2：同 pattern 稀疏矩阵手写同序加法与 SparseArrays + 是否逐位一致（决定 Task 7 D1 可行性）
using SparseArrays, LinearAlgebra, Random
Random.seed!(20260819)

ok1 = true
for trial in 1:20
    n = rand(200:800)
    A = sprand(n, n, 0.005); A = A + A' + 2.0I   # 对称正定近似，模拟 M - θKdt
    b = randn(n)
    if (lu(A) \ b) != (A \ b); ok1 = false; end
end
println("lu(A)\\b == A\\b : ", ok1)

ok2 = true
for trial in 1:20
    n = rand(400:1200)
    Kb = sprand(n, n, 0.004)                     # 模拟 K_bulk
    Kc = sprand(n, n, 0.0008)                    # 模拟 K_coh（更稀疏）
    Kt = Kb + Kc                                 # 并集 pattern，只算一次
    rv = rowvals(Kt)
    manual = copy(nonzeros(Kt))
    for col in 1:size(Kt, 2)                     # 对 Kt 每个存储位置恰一次 a+b 加法
        for k in Kt.colptr[col]:(Kt.colptr[col+1]-1)
            i = rv[k]
            manual[k] = Kb[i, col] + Kc[i, col]
        end
    end
    Kman = SparseMatrixCSC(size(Kt,1), size(Kt,2), copy(Kt.colptr), copy(rv), manual)
    if Kman != Kt; ok2 = false; end
end
println("手写同序加法 == SparseArrays + : ", ok2)
```

（说明：等价性 2 中每个存储位置恰一次 `a+b` 加法、无重排，`Kb[i,col]`/`Kc[i,col]` 返回该位置存储值或 0.0——与 SparseArrays `+` 的合并语义逐位相同的假设由本实验验证。）

```bash
GKSwstype=100 JULIA_NUM_THREADS=1 /d/Julia-1.11.2/bin/julia.exe --startup-file=no docs/planning-with-files/29_仿真提速-SPMe热CZM/probe_equivalence.jl
```

两个 `true/false` 写入 findings.md §事实核查。

- [ ] **Step 4: CZM 向量化计划残留核对**

读 `docs/superpowers/plans/2026-04-20-czm-vectorized-solver-plan.md` 的任务勾选状态、`docs/planning-with-files/09_向量化CZM/progress.md`（若存在），并与当前代码对照（已知已落地：`CZMAssemblyCache` K_bulk 缓存、`assemble_czm_system` 的 K_coh pattern 复用+`fill!` 清零、`geom_cache`、`mul!` 无分配；已知残留：`src/czm.jl:569` `f_int_bulk = K_bulk * u` 与 `:575` `f_int_total = f_int_bulk + f_int_coh` 每次分配、`:572` `K_total = K_bulk + K_coh` 每次分配、`:251` `K_coh[dofs[a],dofs[b]] +=` 标量稀疏索引累加、`src/CzmSolve.jl:113` 线搜索全量组装）。残留清单写入 findings.md §事实核查，Task 7 只做残留且不与既有未完成项冲突的部分。

- [ ] **Step 5: 汇总占比与批次裁决写入 findings.md**

从 Task 1 锚点运行日志抄录四模块 ratio% 到 §实测占比。
（可选交叉确认：spec §2 允许 Profile 火焰图佐证占比。注意 `testexample.jl:56` 已设 `debug_coupling=true` 且 timing 累计本就在 result 键中——timing 开销已是基线固有口径，占比与基线一致，故 Profile 降为可选；若跳过，在 findings §实测占比 注明"占比采用基线 debug 固有口径，Profile 未做"。）
§批次裁决按以下规则逐任务写 开/关 与执行顺序（占比 ≥10% 才开）：
- Task 3（热矩阵缓存）：门 = 3a 中热矩阵块相等（`K identical` 含 KT_BC；若全局 false 但热块相等仍可开）
- Task 4（全局矩阵缓存）：门 = 3a 为 true
- Task 5（单元循环分配消除）：门 = SPMe 求解 ratio ≥10%
- Task 6（按 dt 分解缓存）：门 = 3a 为 true 且等价性 1 为 true；依赖 Task 4 已合并
- Task 7（CZM 装配优化）：门 = CZM ratio ≥10% 且 Step 4 残留核对完成；D1 手写加法部分另需等价性 2 为 true
- Task 8（收尾）：恒开

同时写入按占比排序后的执行顺序（例：`Task 7 → Task 3 → Task 4 → Task 6 → Task 5`）与总墙钟目标（各开门任务对应模块的 timing 合计 × 保守折减系数 0.5，即"预期总提速"值）。

- [ ] **Step 6: 提交**

```bash
git add docs/planning-with-files/仿真提速-SPMe热CZM
git commit -m "perf(批次0b): 矩阵常量性探针、等价性实验、CZM残留核对与批次裁决"
```

---

### Task 3: [门: findings §批次裁决] A 组 — 热矩阵不变量缓存

**Files:**
- Create: `src/PerfCache.jl`
- Modify: `src/ThermalDistributed.jl:1-47`（`ThermalDistributed2D` 加缓存分支）
- Modify: `src/ThermalDistributed.jl:49-158`（`apply_convection_bc`/`apply_cool_method` 加 `apply_to_K::Bool=true` kwarg）
- Modify: `src/ThermalDistributed.jl:160-196`（`ThermalDistributed2D_BC` 加缓存分支）
- Modify: `src/JuBat.jl`（include PerfCache.jl）
- Modify: `src/SetCase.jl:100-117`（Case 加 `perf_cache` 字段）

**Interfaces:**
- Consumes: Task 2 §事实核查 3a（KT 含 BC 后跨步常量）
- Produces:
  - `mutable struct PerfCache`（src/PerfCache.jl），含字段 `thermal::Union{Nothing, ThermalMatrixCache}`、`spme_loop::Union{Nothing, SpmeLoopCache}`、`mat_by_dt::Dict{Float64,Any}`、`factor_by_dt::Dict{Float64,Any}`（后两个供 Task 6 用，本任务只建空 Dict）
  - `mutable struct ThermalMatrixCache`，字段：`MT::SparseMatrixCSC{Float64,Int64}`、`KT::SparseMatrixCSC{Float64,Int64}`、`KT_bc::SparseMatrixCSC{Float64,Int64}`、`mesh_id::UInt`
  - `function ensure_perf_cache!(case::Case)::PerfCache`（懒创建 `case.perf_cache`）
  - `case.perf_cache::Union{Nothing, PerfCache}` 新字段（5 参数兼容构造器同步加 nothing）

- [ ] **Step 1: 建 PerfCache.jl 与 Case 字段**

`src/PerfCache.jl`：

```julia
"""
    PerfCache

SPMe-分布式热-CZM 支线性能缓存（spec 2026-08-19 §5）。
所有缓存对象均已通过探针验证为跨时间步不变（bit 级一致前提），
按 objectid 失效判定，挂在 Case.perf_cache 上。
"""
mutable struct ThermalMatrixCache
    MT::SparseMatrixCSC{Float64,Int64}
    KT::SparseMatrixCSC{Float64,Int64}
    KT_bc::SparseMatrixCSC{Float64,Int64}
    mesh_id::UInt
end

mutable struct SpmeLoopCache   # Task 5 使用，本任务先定义空壳
    ws_pool::Union{Nothing, Vector{Dict{String,Union{Array{Float64},Float64}}}}
    var_list::Union{Nothing, Vector{String}}
    variables_elems::Union{Nothing, Vector{Dict{String,Union{Array{Float64},Float64}}}}
end
SpmeLoopCache() = SpmeLoopCache(nothing, nothing, nothing)

mutable struct PerfCache
    thermal::Union{Nothing, ThermalMatrixCache}
    spme_loop::SpmeLoopCache
    mat_by_dt::Dict{Float64,Any}
    factor_by_dt::Dict{Float64,Any}
end
PerfCache() = PerfCache(nothing, SpmeLoopCache(), Dict{Float64,Any}(), Dict{Float64,Any}())

function ensure_perf_cache!(case::Case)
    if case.perf_cache === nothing
        case.perf_cache = PerfCache()
    end
    return case.perf_cache
end
```

`src/SetCase.jl`：`mutable struct Case` 末尾加一行 `perf_cache::Union{Nothing, PerfCache}`；5 参数兼容构造器改为：

```julia
function Case(param_dim, param, opt, mesh, index)
    Case(param_dim, param, opt, mesh, index, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end
```

然后 `grep -rn "Case(" src/ example/ --include="*.jl"` 核对所有 `Case(` 位置构造调用点（预期只有 `SetCase.jl` 主路径与 5 参数兼容构造器），逐点同步补 `nothing`。`src/JuBat.jl` 在 include 列表（`SetCase.jl` include 之后）加 `include("PerfCache.jl")`。

**验证**（手动一次性，不入库）：

```bash
GKSwstype=100 JULIA_NUM_THREADS=1 /d/Julia-1.11.2/bin/julia.exe --startup-file=no -e 'include("src/JuBat.jl"); using .JuBat; c = JuBat.SetCase(JuBat.ChooseCell("Jellyroll"), JuBat.Option()); println(c.perf_cache === nothing ? "字段OK" : "ERROR")'
```

- [ ] **Step 2: BC 函数加 apply_to_K 开关**

`apply_convection_bc` 与 `apply_cool_method` 签名加 `apply_to_K::Bool=true`，循环体内每个 K 写入行改为 `apply_to_K && (K[a, a] += ke11)` 形式（F 写入行不动）。例（apply_convection_bc）：

```julia
function apply_convection_bc(K, F, mesh, is_outer, case; edge_cache=nothing, apply_to_K::Bool=true)
    ...
        apply_to_K && (K[a, a] += ke11; K[a, b] += ke12)
        apply_to_K && (K[b, a] += ke12; K[b, b] += ke22)
        F[a] += fe1; F[b] += fe2
    ...
```

`apply_cool_method` 的 surface 分支同理（`K[ni, nj] -= wt * Ni_g * Nj_g` 包进 `apply_to_K &&`），tab 分支同理。默认 `true` 时数值路径与原代码逐位一致。

- [ ] **Step 3: ThermalDistributed2D 与 ThermalDistributed2D_BC 加缓存分支**

`ThermalDistributed2D` 开头加缓存分支（矩阵命中则跳过 5 次 `Assemble` 与 4 次稀疏加法；FT 每步照原样算）：

```julia
function ThermalDistributed2D(case::Case, variables::Dict{String,<:Any})
    mesh = case.mesh["thermal2D"]
    param = case.param
    nnode = mesh.nlen
    pc = ensure_perf_cache!(case)
    if pc.thermal !== nothing && pc.thermal.mesh_id == objectid(mesh)
        MT, KT = pc.thermal.MT, pc.thermal.KT
    else
        # ……原有 MT/KT 组装代码不动（:17-37）……
        pc.thermal = ThermalMatrixCache(MT, KT, KT, objectid(mesh))  # KT_bc 先占位为 KT
    end
    # ========== 载荷向量 ==========（原有 :39-44 不动，FT 每步新算）
    ...
    return MT, KT, FT
end
```

`ThermalDistributed2D_BC` 加缓存分支（KT_bc 命中直接返回；F 每步走原路径但 K 写入关闭）：

```julia
function ThermalDistributed2D_BC(KT, FT, case::Case, t::Float64)
    mesh = case.mesh["thermal2D"]
    pc = ensure_perf_cache!(case)
    if pc.thermal !== nothing && pc.thermal.mesh_id == objectid(mesh) && pc.thermal.KT_bc !== pc.thermal.KT
        K = pc.thermal.KT_bc
        F = copy(FT)
        F = apply_convection_bc(K, F, mesh, nothing, case; edge_cache=case.geometry.boundary_edges, apply_to_K=false)[2]
        F = apply_cool_method(K, F, mesh, case; apply_to_K=false)[2]
        return K, F
    end
    # 首次：原路径完整跑一遍（copy + 全部 BC），随后把结果缓存为 KT_bc
    K = copy(KT)
    F = copy(FT)
    K, F = apply_convection_bc(K, F, mesh, nothing, case; edge_cache=case.geometry.boundary_edges)
    K, F = apply_cool_method(K, F, mesh, case)
    pc.thermal === nothing || (pc.thermal.KT_bc = copy(K))
    return K, F
end
```

注意：`ThermalDistributed2D_BC` 的 F 路径（copy(FT) + 逐边/逐高斯点 `+=`）与原实现保持同一边顺序与高斯顺序 → 逐位一致。KT_bc 占位判断 `KT_bc !== KT` 用于区分"还没算过 BC"（构建时占位为 KT 本身）。

**等价性自检**（一次性脚本，可放 /tmp）：

```julia
include(joinpath(@__DIR__, "../src/JuBat.jl")); using .JuBat
# ……（probe_matrix_constancy.jl 的 case 构造块）……
y0 = JuBat.ModelInitialisation_MultiSPMe(case)
_, _, FT1, vars1, _ = JuBat.CallModel(case, y0, 0.0; jacobi="update")   # 首次：原路径建缓存
Kc1 = copy(case.perf_cache.thermal.KT_bc)
_, _, FT2, vars2, _ = JuBat.CallModel(case, y0b, 0.0; jacobi="update")  # 二次：走缓存
Kc2 = case.perf_cache.thermal.KT_bc
println("KT_bc 复用不变: ", Kc1 == Kc2)
```

- [ ] **Step 4: 基线验收（四判据）**

```bash
{ time GKSwstype=100 JULIA_NUM_THREADS=1 /d/Julia-1.11.2/bin/julia.exe --startup-file=no example/testexample.jl ; } > /tmp/perf_A.log 2>&1
```

对照锚点：退出码、网格/步数、metrics、PNG SHA（或三判据）。任一不一致 → `git revert`，findings §批次记录写"失败+现象"并结束本任务。一致 → 记录"热分布式模型"timing 前后对比到 §批次记录。

- [ ] **Step 5: 提交**

```bash
git add src/PerfCache.jl src/ThermalDistributed.jl src/JuBat.jl src/SetCase.jl
git commit -m "perf(热): 热矩阵不变量缓存，跳过每步5次Assemble与BC重复装配（bit级一致）"
```

---

### Task 4: [门: findings §批次裁决] C1 — 全局矩阵与 blockdiag 结果缓存

**Files:**
- Modify: `src/PerfCache.jl`（PerfCache 加 `global_MK::Union{Nothing, Tuple{SparseMatrixCSC{Float64,Int64},SparseMatrixCSC{Float64,Int64}}}` 字段，构造器补 nothing）
- Modify: `src/CallModel.jl:135-138, 175-186`（`CallModel_MultiSPMe` 第 5/7/8 步加缓存分支）

**Interfaces:**
- Consumes: Task 2 §事实核查（全局 M/K 跨步常量 = true）；Task 3 的 `ensure_perf_cache!`
- Produces: `case.perf_cache.global_MK`（缓存的 `(M, K)` 元组）；`CallModel_MultiSPMe` 签名不变，返回的 M/K 在缓存命中时为同一对象（调用方 `Solve.jl` 只读不写，`K_old = copy(K_new)` 已有 copy 保护）

- [ ] **Step 1: CallModel_MultiSPMe 装配段加缓存分支**

`src/CallModel.jl` 第 5 步（:135-138）与第 7/8 步（:175-186）改为（Task 3 已让 `ThermalDistributed2D`/`ThermalDistributed2D_BC` 内部带矩阵缓存，本任务的 `global_MK` 只负责两级 blockdiag 结果缓存；FT 与 `F_T_bc` 每步必须新算，不能跳过）：

```julia
    pc = ensure_perf_cache!(case)

    # 5) 装配电化学全局矩阵（缓存命中时跳过；单元循环仍照跑，其深化留 Task 5 可选 B4）
    if pc.global_MK === nothing
        M_chem = blockdiag(M_elems...)
        K_chem = blockdiag(K_elems...)
    end
    F_chem = vcat(F_elems...)
    ...（第 6 步热源计算 :140-173 不动）...
    # 7) 装配热学矩阵（矩阵走 Task 3 缓存，FT/F_T_bc 每步新算）
    t_thermal_ns = time_ns()
    t_czm_model_s = 0.0
    ...（原 :142-157 czm/热源分支不动）...
    MT, KT, FT = ThermalDistributed2D(case, variables)
    t_ratio = 1.0
    MT = t_ratio == 1.0 ? MT : MT .* t_ratio   # .*1.0 为 IEEE 恒等，跳过分配与逐元素乘
    KT_bc, F_T_bc = ThermalDistributed2D_BC(KT, FT, case, t)
    t_thermal_s = (time_ns() - t_thermal_ns) * 1e-9

    # 8) 全局拼装（缓存命中直接复用；调用方 Solve.jl 对 M/K 只读，K_old = copy(K_new) 已有保护）
    if pc.global_MK === nothing
        M = blockdiag(M_chem, MT)
        K = blockdiag(K_chem, KT_bc)
        pc.global_MK = (M, K)
    else
        M, K = pc.global_MK
    end
    F = [F_chem; F_T_bc]
```

- [ ] **Step 2: 基线验收（四判据）**

同 Task 3 Step 4 命令，日志存 `/tmp/perf_C1.log`。一致 → 记录 timing 到 §批次记录（关注"热分布式模型"与总墙钟）。

- [ ] **Step 3: 提交**

```bash
git add src/PerfCache.jl src/CallModel.jl
git commit -m "perf(全局): 缓存常量全局M/K，跳过每步两级blockdiag与MT缩放（bit级一致）"
```

---

### Task 5: [门: findings §批次裁决] B 组 — 单元循环分配消除

**Files:**
- Modify: `src/PerfCache.jl`（`SpmeLoopCache` 已有 ws_pool/var_list 字段）
- Modify: `src/CallModel.jl:118-131`（ws_pool/var_list 走缓存）
- Modify: `src/CallModel.jl:242-264`（新增 `copy_element_results!` 原位变体并替换 :131 调用点）
- Modify: `src/SPMe.jl:118-126`（var_list 参数化）

**Interfaces:**
- Consumes: Task 3 的 `SpmeLoopCache`（字段 `ws_pool`、`var_list`、`variables_elems`）
- Produces:
  - `copy_element_results!(dst::Dict{String,Union{Array{Float64},Float64}}, vars_e)` —— 原位变体，键集合与原 `copy_element_results` 完全一致的 17 键
  - `SPMe_variables!(ws, case, yt, t; I_app, T_e, var_list)` —— 新 kwarg `var_list::Union{Nothing,Vector{String}}=nothing`，传入时跳过 `collect(keys(case.index))` 与 filter

- [ ] **Step 1: SPMe_variables! 加 var_list kwarg**

`src/SPMe.jl:101` 签名加 `var_list::Union{Nothing,Vector{String}}=nothing`；:118-121 改为：

```julia
    if var_list === nothing
        var_list = collect(keys(case.index))
        if T_e !== nothing
            var_list = filter(k -> k != "temperature", var_list)
        end
    end
    for i in var_list
        if haskey(ws, i)
            ws[i] = yt[case.index[i]]
        end
    end
```

（var_list 只控制赋值遍历顺序，赋值互不依赖 → 顺序无关数值；缓存的 var_list 来自同一 `case.index` 的同一构造顺序。）

- [ ] **Step 2: 新增 copy_element_results! 原位变体**

`src/CallModel.jl` 在原函数后新增（键列表逐字对应原函数 `copy_element_results`（:242-264）的 17 键，勿增删）：

```julia
"""
    copy_element_results!(dst, vars_e)

copy_element_results 的原位变体：把 17 个结果键的引用写入复用的 dst Dict。
值引用每次来自当步 workspace 的新数组/标量，跨步覆盖引用安全。
"""
function copy_element_results!(dst::Dict{String,Union{Array{Float64},Float64}}, vars_e)
    for k in (
        "negative electrode overpotential", "positive electrode overpotential", "cell voltage",
        "negative electrode exchange current density", "positive electrode exchange current density",
        "negative electrode interfacial current density", "positive electrode interfacial current density",
        "negative electrode open circuit potential", "positive electrode open circuit potential",
        "temperature",
        "electrolyte lithium concentration at negative electrode Gauss point",
        "electrolyte lithium concentration at positive electrode Gauss point",
        "electrolyte lithium concentration at separator Gauss point",
        "negative particle surface lithium concentration", "positive particle surface lithium concentration",
        "negative particle lithium concentration", "positive particle lithium concentration",
    )
        dst[k] = vars_e[k]
    end
    return dst
end
```

- [ ] **Step 3: 单元循环接缓存（ws_pool、variables_elems、var_list 跨步复用）**

`src/CallModel.jl:112-133` 改为（`variables_elems` 首步创建 ne 个 Dict 后跨步复用——spec B2 的预分配复用本意；值引用每步被 `copy_element_results!` 覆盖为当步新数组，下游 `compute_heat_sources` 与 ：160-172 当步读完即弃，无跨步别名）：

```julia
    pc = ensure_perf_cache!(case)
    slc = pc.spme_loop
    if slc.ws_pool === nothing
        slc.ws_pool = [create_element_workspace(case) for _ in 1:Threads.nthreads()]
    end
    if slc.variables_elems === nothing
        slc.variables_elems = [Dict{String,Union{Array{Float64},Float64}}() for _ in 1:ne]
    end
    if slc.var_list === nothing
        slc.var_list = filter(k -> k != "temperature", collect(keys(case.index)))
    end
    ws_pool = slc.ws_pool
    variables_elems = slc.variables_elems

    M_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    K_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    F_elems = Vector{Vector{Float64}}(undef, ne)

    t_spme_ns = time_ns()
    Threads.@threads for e in 1:ne
        tid = Threads.threadid()
        ws_e = ws_pool[tid]
        M_e, K_e, F_e, vars_e = SPMe_element(case, yt_chem[e], t, e; I_e=I_e[e], T_e=Te_prev[e], jacobi=jacobi, workspace=ws_e, var_list=slc.var_list)
        M_elems[e] = M_e
        K_elems[e] = K_e
        F_elems[e] = vec(F_e)
        copy_element_results!(variables_elems[e], vars_e)
    end
    t_spme_s = (time_ns() - t_spme_ns) * 1e-9
```

（`ws_pool` 每步重建与 `var_list` 的 keys 收集+filter 均已消除；`SPMe_element` 需加 kwarg `var_list=nothing` 并透传给 `SPMe_variables!`，`src/SPMe.jl:37-93`。）

- [ ] **Step 4: 基线验收（四判据）+ 提交**

命令同 Task 3 Step 4，日志 `/tmp/perf_B.log`，记录"SPMe 求解"timing 前后对比。

```bash
git add src/PerfCache.jl src/CallModel.jl src/SPMe.jl
git commit -m "perf(SPMe循环): ws_pool/var_list跨步缓存，copy_element_results原位化（bit级一致）"
```

- [ ] **Step 5（可选，门: 本批后 SPMe ratio 仍 ≥10%）: B4 — 单元矩阵跳过**

`SPMe_element` 加 kwarg `with_matrices::Bool=true`；`false` 时 `ElectrodeDiffusion`/`ElectrolyteDiffusion` 组装段跳过、`M_e/K_e` 返回缓存引用（`PerfCache` 加 `elem_matrices::Union{Nothing,Tuple{Vector,SparseMatrixCSC,SparseMatrixCSC}}`——首步收集每单元 `M_e/K_e` 与 `blockdiag` 结果）。门触发才做，做完同样四判据验收 + 单独 commit。

---

### Task 6: [门: findings §批次裁决；依赖 Task 4 已合并] C2 — 按 dt 缓存线性分解因子

**Files:**
- Modify: `src/Solve.jl:211-215`（主循环求解段）

**Interfaces:**
- Consumes: Task 4 的 `pc.global_MK`（M/K 常量 → `Mt/Kt` 只随 dt 变）；Task 2 等价性 1 结论 `lu(A)\b == A\b`；Task 3 的 `ensure_perf_cache!`、`pc.mat_by_dt`、`pc.factor_by_dt`
- Produces: 无新接口（Solve 内部优化）

- [ ] **Step 1: 主循环求解段改造**

`src/Solve.jl:211-215` 改为：

```julia
        pc = ensure_perf_cache!(case)
        if !haskey(pc.mat_by_dt, dt)
            pc.mat_by_dt[dt] = convert(SparseMatrixCSC{Float64,Int}, M_new - theta * K_new * dt)
        end
        Mt = pc.mat_by_dt[dt]
        Kt = (1 - theta) * K_old * dt + M_new
        Ft = theta * F_new * dt + (1 - theta) * F_old * dt
        rhs = Kt * y_old[vc] + Ft
        if haskey(pc.factor_by_dt, dt)
            y_c = pc.factor_by_dt[dt] \ rhs
        else
            F_lu = lu(Mt; check=true)
            pc.factor_by_dt[dt] = F_lu
            y_c = F_lu \ rhs
        end
        y_new = vcat(y_c, y_phi)
```

要点：`Mt` 缓存的是同一矩阵对象（同 dt 下 `M_new - theta*K_new*dt` 的浮点结果逐位可再现，缓存省去重复构造与分配）；`Kt`/`Ft`/`rhs` 构造式逐字保留原顺序；`lu(Mt)\rhs` 与原 `Mt \ rhs` 的等价性由 Task 2 等价性 1 背书。若等价性 1 为 false（门已关）则本任务整体不执行。
注意 dt 缓存 key 用原浮点值精确匹配（dt 由 `dt/2`、`dt*2`、`abs(RunTime[vt]-t)` 产生，相同来源产生相同位模式；不同位模式的 dt 各占一条缓存条目，条目数上限为不同 dt 值个数，本场景为个位数到十位数）。
**边界**：`ModelInitialisation` 后的首次求解（Solve.jl:180 `y_c = (M_old - K_old*dt_init)\...`）不动——dt_init=1e-8 一次性。

- [ ] **Step 2: 基线验收（四判据）+ 提交**

命令同 Task 3 Step 4，日志 `/tmp/perf_C2.log`。本批预期是最大单项收益（每步省一次 UMFPACK 符号+数值分解），记录总墙钟与四模块 timing。

```bash
git add src/Solve.jl
git commit -m "perf(求解): 按dt缓存Mt矩阵与LU分解因子，连续同dt步复用ldiv（bit级一致）"
```

---

### Task 7: [门: findings §批次裁决] D 组 — CZM 装配分配消除

**Files:**
- Modify: `src/czm.jl:546-578`（`assemble_coupled_system` K_total/f_int 预分配 + `assemble_K` 开关）
- Modify: `src/czm.jl:88-131`（`assemble_czm_system` nzval 直索引；workspace 加 `cohesive_nzidx`）
- Modify: `src/CzmSolve.jl:107-116`（线搜索传 `assemble_K=false`）
- Modify: `src/czm.jl`（`CZMAssemblyWorkspace` 加 `K_total_buf`、`f_int_bulk_buf`、`f_int_total_buf`、`cohesive_nzidx` 字段——struct 定义处按现有构造函数同步）

**Interfaces:**
- Consumes: Task 2 等价性 2 结论（手写同序加法 `==` SparseArrays `+`）；Task 2 残留核对清单
- Produces: `assemble_coupled_system(...; assemble_K::Bool=true)` 新 kwarg（`false` 时首个返回值为 `nothing`，仅供只取 f_int 的调用方使用）；其余签名不变

- [ ] **Step 1: workspace 加缓冲字段**

`CZMAssemblyWorkspace`（`src/czm.jl` 中 struct 定义与构造函数）加字段并初始化为 nothing：

```julia
    K_total_buf::Union{Nothing, SparseMatrixCSC{Float64,Int64}}
    K_total_mapK::Union{Nothing, Vector{Int}}      # K_total 每存储位置 → K_bulk nzval 下标（0=无）
    K_total_mapC::Union{Nothing, Vector{Int}}      # K_total 每存储位置 → K_coh nzval 下标（0=无）
    f_int_bulk_buf::Union{Nothing, Vector{Float64}}
    f_int_total_buf::Union{Nothing, Vector{Float64}}
    cohesive_nzidx::Union{Nothing, Matrix{Int}}    # n_coh × 64：每单元 64 项的 K_coh nzval 下标
```

- [ ] **Step 2: assemble_czm_system 组装行改 nzval 直索引**

先在 `src/czm.jl` 加列内二分查找辅助函数（文件内私有）：

```julia
# 列内二分查找稀疏矩阵 (i,j) 的 nzval 下标；不存在返回 0（pattern 固定时仅供初始化使用）
function _nz_index(K::SparseMatrixCSC, i::Int, j::Int)
    lo, hi = K.colptr[j], K.colptr[j+1] - 1
    while lo <= hi
        mid = (lo + hi) >>> 1
        if K.rowval[mid] == i; return mid
        elseif K.rowval[mid] < i; lo = mid + 1
        else; hi = mid - 1; end
    end
    return 0
end
```

`src/czm.jl:101-128` pattern 首建块内、`sparse(I_pat, J_pat, V_pat, ndof, ndof)` 与 `ws.K_coh = K_coh` 之后，同步构建 `cohesive_nzidx`（pattern 与下方组装循环的 a/b 顺序固定一一对应）：

```julia
        nzidx = zeros(Int, n_coh, 64)
        for i in 1:n_coh
            if geom_cache !== nothing
                dofs_ = geom_cache[i].dofs
            else
                elem = czm_mesh.cohesive_elements[i]
                n1, n2 = elem.nodes_bottom; n4, n3 = elem.nodes_top
                dofs_ = [2*n1-1,2*n1,2*n2-1,2*n2,2*n3-1,2*n3,2*n4-1,2*n4]
            end
            c = 0
            for a in 1:8, b in 1:8
                c += 1
                nzidx[i, c] = _nz_index(K_coh, dofs_[a], dofs_[b])
            end
        end
        ws.cohesive_nzidx = nzidx
```

组装循环 `src/czm.jl:248-253` 改为（`K_coh[dofs[a],dofs[b]] += v` 的 getindex/setindex 语义 = 读当前存储值→加→写回，与 nzval 直加逐位一致）：

```julia
        if ws.cohesive_nzidx !== nothing
            c = 0
            for a in 1:8
                ws.f_int_coh[dofs[a]] += ws.f_int_e[a]
                for b in 1:8
                    c += 1
                    K_coh.nzval[ws.cohesive_nzidx[i, c]] += ws.K_e[a, b]
                end
            end
        else
            # 旧 workspace 回退原写法
            for a in 1:8
                ws.f_int_coh[dofs[a]] += ws.f_int_e[a]
                for b in 1:8
                    K_coh[dofs[a], dofs[b]] += ws.K_e[a, b]
                end
            end
        end
```

- [ ] **Step 3: assemble_coupled_system 预分配 + assemble_K 开关**

`src/czm.jl:546-578` 改为：

```julia
    # 固体内力（预分配 mul!，与 K_bulk * u 同一稀疏 matvec 路径）
    if ws !== nothing
        if ws.f_int_bulk_buf === nothing || length(ws.f_int_bulk_buf) != ndof
            ws.f_int_bulk_buf = zeros(Float64, ndof)
        end
        mul!(ws.f_int_bulk_buf, K_bulk, u)
        f_int_bulk = ws.f_int_bulk_buf
    else
        f_int_bulk = K_bulk * u
    end

    # 总内力（copyto! + .+= 与向量 + 逐位一致）
    ... f_int_total = f_int_bulk + f_int_coh 的预分配等价实现：
    if ws !== nothing
        ws.f_int_total_buf === nothing && (ws.f_int_total_buf = zeros(Float64, ndof))
        copyto!(ws.f_int_total_buf, f_int_bulk)
        ws.f_int_total_buf .+= f_int_coh
        f_int_total = ws.f_int_total_buf
    else
        f_int_total = f_int_bulk + f_int_coh
    end

    if assemble_K
        if ws !== nothing && ws.K_total_buf !== nothing
            # 预分配同序加法：每存储位置恰一次 K_bulk 值 + K_coh 值（等价性 2 已验证）
            nzT = nonzeros(ws.K_total_buf)
            for k in 1:length(nzT)
                nzT[k] = (ws.K_total_mapK[k] > 0 ? K_bulk.nzval[ws.K_total_mapK[k]] : 0.0) +
                         (ws.K_total_mapC[k] > 0 ? K_coh.nzval[ws.K_total_mapC[k]] : 0.0)
            end
            K_total = ws.K_total_buf
        else
            K_total = K_bulk + K_coh
            if ws !== nothing
                ws.K_total_buf = copy(K_total)
                mapK = zeros(Int, nnz(K_total)); mapC = zeros(Int, nnz(K_total))
                for col in 1:size(K_total,2)
                    for k in K_total.colptr[col]:(K_total.colptr[col+1]-1)
                        mapK[k] = _nz_index(K_bulk, K_total.rowval[k], col)
                        mapC[k] = _nz_index(K_coh,  K_total.rowval[k], col)
                    end
                end
                ws.K_total_mapK = mapK; ws.K_total_mapC = mapC
            end
        end
    else
        K_total = nothing
    end
    return K_total, f_int_total, separations, tractions
```

map 构建已并入上面的 else 分支（首次 `+` 后一次性构建，之后走预分配直加路径）。

签名加 `assemble_K::Bool=true`。**注意**：`K_coh` 在 ws 中是共享对象，同一次 `update_czm_damage!` 的多次调用间 `assemble_czm_system` 每次都 fill! 重算——本任务的 K_total 直加每次都用当次 K_coh 值，语义不变。
**范围注记**：spec D3（`clone_czm_mesh_with_damage` 每步克隆消除）不进本任务；若 Task 8 Step 2 判定 CZM 收益不及预期，D3 连同更深的装配项一起走增补机制。

- [ ] **Step 4: 线搜索路径跳过 K_total 组装**

`src/CzmSolve.jl:113` 改为：

```julia
        _, f_int_trial, _, _ = assemble_coupled_system(czm_mesh, u_trial, param_cache; damage_states=damage_states, K_bulk_cached=K_bulk_cached, geom_cache=geom_cache, ws=ws, visc_beta=visc_beta, assemble_K=false)
```

（线搜索只用 f_int_trial；跳过 K_total 组装不改变任何被使用的数值。）

- [ ] **Step 5: CZM 单步等价性自检 + 基线验收（四判据）+ 提交**

单步自检（一次性脚本）：case 构造同探针；固定 `u`、`damage_states` 副本，改造前后各调一次 `solve_czm_step`（或直接跑基线对比）——**以基线四判据为准**（CZM 损伤场进 PNG 与 metrics），单步脚本可省略，直接四判据。

```bash
{ time GKSwstype=100 JULIA_NUM_THREADS=1 /d/Julia-1.11.2/bin/julia.exe --startup-file=no example/testexample.jl ; } > /tmp/perf_D.log 2>&1
```

记录"CZM 模型"timing 前后对比到 §批次记录。

```bash
git add src/czm.jl src/CzmSolve.jl
git commit -m "perf(CZM): K_total/f_int预装配缓冲、nzval直索引、线搜索跳过K_total（bit级一致）"
```

---

### Task 8: 收尾 — 汇总、目标确认、文档同步

**Files:**
- Modify: `docs/planning-with-files/29_仿真提速-SPMe热CZM/findings.md`（§批次记录补总表）
- Modify: `docs/planning-with-files/29_仿真提速-SPMe热CZM/progress.md`
- Modify: `docs/planning-with-files/index.md`（任务收尾状态）
- Modify: `docs/superpowers/specs/2026-08-19-perf-optimization-spme-thermal-czm-design.md`（状态行改为"已实施"，附结论一行）

**Interfaces:**
- Consumes: 全部前置任务的 §批次记录
- Produces: 汇总报告（findings.md §批次记录总表：每批 timing 前后、四判据结果、总墙钟前后、相对 Task 2 §批次裁决预期目标的达成情况）

- [ ] **Step 1: 汇总表**

findings.md §批次记录追加总表：| 批次 | 模块 timing 前 | 后 | 总墙钟前 | 后 | 判据 |。对比 §批次裁决 的预期目标值，写达成/未达成与原因（未开门的任务标"门关未做"）。

- [ ] **Step 2: 增补机制说明（如需要）**

若某开门批次完成后其模块 ratio 仍 ≥10%（收益不及预期），在 findings 写明"建议增补批次"及具体入手点（如 Task 5 Step 5 B4、或 CZM 需另行细化的装配项），**不在本计划内实施**——增补须回到 spec 补充后另立批次。

- [ ] **Step 3: 文档同步 + 提交**

```bash
git add docs/planning-with-files/仿真提速-SPMe热CZM docs/planning-with-files/index.md docs/superpowers/specs/2026-08-19-perf-optimization-spme-thermal-czm-design.md
git commit -m "perf(收尾): 提速批次汇总与规划文档收尾"
```

---

## 执行门速查（从 findings.md §批次裁决 读）

| 任务 | 门条件 | 依赖 |
|------|--------|------|
| Task 3 热矩阵缓存 | 热矩阵块跨步常量 | Task 2 |
| Task 4 全局 M/K 缓存 | 3a 全局 true | Task 2、Task 3 结构（PerfCache） |
| Task 5 单元循环分配消除 | SPMe ratio ≥10% | Task 3（SpmeLoopCache） |
| Task 6 按 dt 分解缓存 | 3a true 且等价性 1 true | Task 4 已合并 |
| Task 7 CZM 装配优化 | CZM ratio ≥10%、残留核对完成；D1 部分另需等价性 2 true | Task 2 |
| Task 8 收尾 | 恒开 | 全部前置结束（含门关跳过的确认记录） |
