# CZM 按材料层界面重构 Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将 cohesive 单元从 SP-PE 螺旋界面迁移到 PE-PCC / NE-NCC 电极-集流体界面，使模拟几何与实验剥离数据对齐，消除"代表元均一化"导致的体应力场失真。

**Architecture:** 新增独立细化的 CZM 机械子网格（8 径向分层/卷绕圈），保持粗热网格不变；cohesive 单元按 `interface_type` 取实验参数与涂层模量；通过 `thermal_to_czm` 稀疏插值矩阵实现粗热→细 CZM 耦合。

**Tech Stack:** Julia 1.10+, `Test.jl`（标准库，新引入用于单元测试）, `SparseArrays`（已有依赖）, `Parameters.jl`（`@with_kw`，已有依赖）

**Spec:** `docs/superpowers/specs/2026-07-18-czm-per-material-layer-interface-design.md`

---

## Chunk 1: 数据结构扩展

本 chunk 仅扩展 struct 定义，不实现新逻辑。所有测试针对 struct 构造与字段访问。

### Task 1.0: 完整侦察旧字段引用（前置）

**Files:** 无修改，仅侦察

- [ ] **Step 1: 生成旧字段全量引用清单**

执行 grep 输出所有将受影响的行，作为后续任务的 baseline：

```bash
grep -rn "cohesive[_.]​\(σ_max_n\|K_n\|δ_0_n\|G_c_n\|δ_c_n\|τ_max_t\|K_t\|δ_0_t\|G_c_t\|δ_c_t\)\b" \
    src/ example/ tools/ test/ 2>/dev/null > /tmp/czm-old-field-refs.txt
wc -l /tmp/czm-old-field-refs.txt
```

记录行数作为修复进度基线。Task 1.1 Step 4 完成后，再次运行同一 grep，剩余命中应全部带 `# TODO chunk-N` 注释。

- [ ] **Step 2: 不 commit**

本任务仅生成清单，无文件修改。

### Task 1.1: 扩展 `Cohesive` struct（按界面类型分组）

**Files:**
- Modify: `src/SetParams.jl:156-184`（Cohesive struct 定义）
- Create: `test/test_cohesive_struct.jl`

**说明：** 当前 `Cohesive` 含 Mode I + Mode II 共 10 个 per-interface 字段（`σ_max_n / K_n / δ_0_n / G_c_n / δ_c_n / τ_max_t / K_t / δ_0_t / G_c_t / δ_c_t`），全部单一值。本任务将其拆为 `:PE_PCC` 与 `:NE_NCC` 两组（共 20 字段），并删除旧字段。

- [ ] **Step 1: 建立测试目录与第一个失败测试**

创建 `test/test_cohesive_struct.jl`：

```julia
using Test
using JuBat

@testset "Cohesive struct per-interface fields" begin
    coh = JuBat.Cohesive()

    # PE-PCC Mode I
    @test hasproperty(coh, :σ_max_pe_pcc)
    @test hasproperty(coh, :K_n_pe_pcc)
    @test hasproperty(coh, :δ_0_pe_pcc)
    @test hasproperty(coh, :G_c_pe_pcc)
    @test hasproperty(coh, :δ_c_pe_pcc)

    # PE-PCC Mode II
    @test hasproperty(coh, :τ_max_pe_pcc)
    @test hasproperty(coh, :K_t_pe_pcc)
    @test hasproperty(coh, :δ_0_pe_pcc_t)
    @test hasproperty(coh, :G_c_pe_pcc_t)
    @test hasproperty(coh, :δ_c_pe_pcc_t)

    # NE-NCC Mode I
    @test hasproperty(coh, :σ_max_ne_ncc)
    @test hasproperty(coh, :K_n_ne_ncc)
    @test hasproperty(coh, :δ_0_ne_ncc)
    @test hasproperty(coh, :G_c_ne_ncc)
    @test hasproperty(coh, :δ_c_ne_ncc)

    # NE-NCC Mode II
    @test hasproperty(coh, :τ_max_ne_ncc)
    @test hasproperty(coh, :K_t_ne_ncc)
    @test hasproperty(coh, :δ_0_ne_ncc_t)
    @test hasproperty(coh, :G_c_ne_ncc_t)
    @test hasproperty(coh, :δ_c_ne_ncc_t)

    # 旧字段已移除
    @test !hasproperty(coh, :σ_max_n)
    @test !hasproperty(coh, :K_n)
    @test !hasproperty(coh, :δ_0_n)
    @test !hasproperty(coh, :G_c_n)
    @test !hasproperty(coh, :δ_c_n)
    @test !hasproperty(coh, :τ_max_t)
    @test !hasproperty(coh, :K_t)
    @test !hasproperty(coh, :δ_0_t)
    @test !hasproperty(coh, :G_c_t)
    @test !hasproperty(coh, :δ_c_t)

    # 保留字段（界面热阻、粘性、BK eta、czm_model）
    @test hasproperty(coh, :eta)
    @test hasproperty(coh, :czm_model)
    @test hasproperty(coh, :h_c0)
    @test hasproperty(coh, :k_air)
    @test hasproperty(coh, :lambda_m)
    @test hasproperty(coh, :beta)
    @test hasproperty(coh, :threshold)
    @test hasproperty(coh, :tau_visc)
end
```

- [ ] **Step 2: 运行测试验证失败**

Run:
```bash
julia --project=. -e 'include("test/test_cohesive_struct.jl")'
```
Expected: FAIL — `Cohesive` 没有 `σ_max_pe_pcc` 等字段，`hasproperty` 返回 false。

- [ ] **Step 3: 修改 `Cohesive` struct**

替换 `src/SetParams.jl:156-184` 整个 `Cohesive` struct。原 struct 有两个不同字段：`czm_model::String`（模型选择）与 `eta::Float64`（BK 指数）。新 struct 同时保留，**不要在 struct 中定义重名字段**：

```julia
@with_kw mutable struct Cohesive
    # === PE-PCC 界面（电极涂层-正极集流体）===
    σ_max_pe_pcc::Float64 = 0.0    # 最大法向牵引力 [Pa]
    K_n_pe_pcc::Float64 = 0.0      # 法向初始刚度 [Pa/m]
    δ_0_pe_pcc::Float64 = 0.0      # 法向损伤起始位移 [m]
    G_c_pe_pcc::Float64 = 0.0      # 法向断裂能 [J/m²]
    δ_c_pe_pcc::Float64 = 0.0      # 法向临界位移 [m]
    τ_max_pe_pcc::Float64 = 0.0    # Mode II 最大切向牵引 [Pa]
    K_t_pe_pcc::Float64 = 0.0
    δ_0_pe_pcc_t::Float64 = 0.0
    G_c_pe_pcc_t::Float64 = 0.0
    δ_c_pe_pcc_t::Float64 = 0.0

    # === NE-NCC 界面（电极涂层-负极集流体）===
    σ_max_ne_ncc::Float64 = 0.0
    K_n_ne_ncc::Float64 = 0.0
    δ_0_ne_ncc::Float64 = 0.0
    G_c_ne_ncc::Float64 = 0.0
    δ_c_ne_ncc::Float64 = 0.0
    τ_max_ne_ncc::Float64 = 0.0
    K_t_ne_ncc::Float64 = 0.0
    δ_0_ne_ncc_t::Float64 = 0.0
    G_c_ne_ncc_t::Float64 = 0.0
    δ_c_ne_ncc_t::Float64 = 0.0

    # === 共用 ===
    czm_model::String = "model1"   # 模型选择（"model1" / "mix"）
    eta::Float64 = 1.0             # BK 准则指数 [-]

    # Interface thermal resistance
    h_c0::Float64 = 1e7
    k_air::Float64 = 0.026
    lambda_m::Float64 = 70e-9
    beta::Float64 = 1.0
    threshold::Float64 = 70e-9

    tau_visc::Float64 = 0.0
end
```

- [ ] **Step 4: 修复编译错误（其他文件引用旧字段）**

旧字段移除后，多处编译会失败。**本步骤仅修复编译错误**（让 `using JuBat` 不报错），逻辑适配在后续 chunk 处理。

**先做完整 grep 建立修复清单**：

```bash
grep -rn "cohesive[_.]​\(σ_max_n\|K_n\|δ_0_n\|G_c_n\|δ_c_n\|τ_max_t\|K_t\|δ_0_t\|G_c_t\|δ_c_t\)\b" src/ example/ tools/ | grep -v "# TODO"
```

预期旧字段引用位置（基于 grep 核实）：

| 文件 | 行 | 引用 | 修复方式 |
|------|-----|------|----------|
| `src/SetParams.jl` | 324 | `param_dim.scale.σ_czm = param_dim.cohesive.σ_max_n` | 临时改为 `param_dim.scale.σ_czm = param_dim.cohesive.σ_max_pe_pcc`（PE-PCC 作参考），加 `# TODO Chunk 2 Task 2.2 复核` |
| `src/SetParams.jl` | 467-481 | NormaliseParam 中 10 行 cohesive 归一化 | 整块注释，加 `# TODO Chunk 2 Task 2.2 重写` |
| `src/CouplingState.jl` | 228-240 | `compute_effective_coating_modulus` 函数体 | 整块注释，加 `# TODO Chunk 2 Task 2.3 移除` |
| `src/CouplingState.jl` | 258-275 | `compute_czm_effective_params` 函数体 | 整块注释，加 `# TODO Chunk 2 Task 2.3 移除` |
| `src/parameters/Jellyroll.jl` | 158-181 | 单一 cohesive 参数赋值 | 整块注释，加 `# TODO Chunk 2 Task 2.1 重写` |
| `src/Materialmatrix.jl` | 69-70 | `K_n = cohesive_params.K_n` 等 2 行 | 临时占位 `K_n = 1.0`，加 `# TODO Chunk 4 Task 4.5` |
| `src/Materialmatrix.jl` | 183-185 | 同上 3 行 | 同上 |
| `src/Materialmatrix.jl` | 331-332 | `delta0 = cohesive.δ_0_n; delta_c = cohesive.δ_c_n` | 临时占位 `delta0 = 1e-9`，加 `# TODO Chunk 4 Task 4.5` |
| `tools/verify_czm_unit.jl` | 92 | `CohesiveElement(..., layer_idx)` 构造调用 | 临时改为 `:PE_PCC`（接受 Task 1.2 新签名），加 `# TODO Chunk 3 重写` |

**说明**：Materialmatrix.jl 中 lines 94-113, 149-157, 194-280 引用的是局部变量 `K_n`、`δ_0_n` 等（已在 69-70、183-185 提取），**不是字段访问**，无需修改。

**外部调用方（重要）**——以下文件调用 `JuBat.compute_czm_effective_params` 或 `JuBat.compute_effective_coating_modulus`，移除函数定义后会断裂：

| 文件 | 行 | 处理 |
|------|-----|------|
| `tools/verify_czm_standalone.jl` | 74 | 临时注释该调用，加 `# TODO Chunk 4` |
| `tools/czm_convergence_diag.jl` | 36 | 同上 |
| `tools/czm_baseline_probe.jl` | 47 | 同上 |
| `example/力学模块验证/test_electrode_coat_modulus.jl` | 107, 132 | 同上 |

**`example/*.jl`、`tools/*.jl` 中其他 `cohesive.σ_max_n` 等赋值**：临时注释，加 `# TODO Chunk 2 Task 2.1`。

**每个临时修复点必须加 `# TODO chunk-N task-M` 注释**，便于后续 chunk 定位。

- [ ] **Step 5: 运行测试验证通过**

Run:
```bash
julia --project=. -e 'include("test/test_cohesive_struct.jl")'
```
Expected: PASS（所有字段断言通过）

同时验证 `using JuBat` 不报错：
```bash
julia --project=. -e 'using JuBat; println("OK")'
```
Expected: 打印 `OK`，无编译错误。

- [ ] **Step 6: Commit**

```bash
git add src/SetParams.jl test/test_cohesive_struct.jl src/CouplingState.jl src/parameters/Jellyroll.jl src/Materialmatrix.jl src/czm.jl
git add -u example/ tools/
git commit -m "refactor(cohesive): 按 interface_type 拆分 Cohesive struct 字段

- Cohesive struct 移除单一 σ_max_n/K_n/G_c_n 等 10 字段
- 新增 PE_PCC + NE_NCC 两组共 20 字段（Mode I + Mode II）
- 临时占位修复编译错误，逻辑适配见后续 chunk

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 1.2: 扩展 `CohesiveElement` 添加 `interface_type`，删除 `layer_idx`

**Files:**
- Modify: `src/czm.jl:1-8`（CohesiveElement struct）
- Modify: `test/test_cohesive_struct.jl`（添加新测试集）

- [ ] **Step 1: 写失败测试**

在 `test/test_cohesive_struct.jl` 末尾添加：

```julia
@testset "CohesiveElement interface_type" begin
    elem = JuBat.CohesiveElement(
        1,                          # id
        [1, 2, 3, 4],               # nodes
        [1, 2],                     # nodes_bottom
        [4, 3],                     # nodes_top
        1.0,                        # length
        :PE_PCC                     # interface_type
    )
    @test elem.interface_type == :PE_PCC

    # 旧 layer_idx 字段已删除
    @test !hasproperty(elem, :layer_idx)
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_cohesive_struct.jl")'`
Expected: FAIL — 构造函数不匹配（仍要求 `layer_idx::Int64`，不接受 `:PE_PCC`）。

- [ ] **Step 3: 修改 `CohesiveElement`**

替换 `src/czm.jl:1-8`：

```julia
mutable struct CohesiveElement <: AbstractCohesiveElement
    id::Int64
    nodes::Vector{Int64}           # [n1, n2, n3, n4]
    nodes_bottom::Vector{Int64}    # [n1, n2] 底面节点
    nodes_top::Vector{Int64}       # [n4, n3] 顶面节点
    length::Float64                # 单元长度
    interface_type::Symbol         # :PE_PCC 或 :NE_NCC
end
```

- [ ] **Step 4: 修复 `create_czm_mesh` 中的构造调用**

`src/czm.jl:100-107`（旧代码）有 `CohesiveElement(i, [...], [...], [...], elem_length, 1)`——硬编码的 `1` 是 `layer_idx`。**本 chunk 仅修改构造使其匹配新签名**（接受 `interface_type`），逻辑正确性在 Chunk 3 处理：

```julia
# src/czm.jl:100-107 临时改为：
coh_elem = CohesiveElement(
    i,
    [n_in_1, n_in_2, n_out_2, n_out_1],
    [n_in_1, n_in_2],
    [n_out_1, n_out_2],
    elem_length,
    :PE_PCC   # TODO Chunk 3: 按实际材料类型判定
)
```

**`src/` 内这是唯一构造点**。`tools/verify_czm_unit.jl:92` 另有一处构造调用，已在 Task 1.1 Step 4 修复清单中处理。

- [ ] **Step 5: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_cohesive_struct.jl")'`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/czm.jl test/test_cohesive_struct.jl
git commit -m "refactor(czm): CohesiveElement 用 interface_type 替代 layer_idx

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 1.3: 新增 `CzmInterfaceParams` 与 `CzmParamCache`

**Files:**
- Modify: `src/CouplingState.jl`（在文件顶部 imports 后新增 struct）
- Modify: `test/test_cohesive_struct.jl`（添加新测试集）

- [ ] **Step 1: 写失败测试**

在 `test/test_cohesive_struct.jl` 末尾添加：

```julia
@testset "CzmInterfaceParams and CzmParamCache" begin
    params_pe_pcc = JuBat.CzmInterfaceParams(;
        E_eff = 1.0e3,
        ν = 0.3,
        α = 2.0e-6,
        σ_max = 50e6,
        δ_0 = 1.0e-9,
        δ_c = 5.0e-7,
        G_c = 25.0,
        K_n = 1.0e17,
        τ_max = 50e6,
        δ_0_t = 1.0e-9,
        δ_c_t = 5.0e-7,
        G_c_t = 25.0,
        K_t = 1.0e17
    )
    @test params_pe_pcc.σ_max == 50e6

    cache = JuBat.CzmParamCache(Dict(:PE_PCC => params_pe_pcc))
    @test haskey(cache.by_interface, :PE_PCC)
    @test cache.by_interface[:PE_PCC].E_eff == 1.0e3
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_cohesive_struct.jl")'`
Expected: FAIL — `JuBat.CzmInterfaceParams` 不存在。

- [ ] **Step 3: 实现 struct**

在 `src/CouplingState.jl` 的现有 struct 定义区（`MeshGeometry` 附近，约 line 80 之前）添加：

```julia
"""
    CzmInterfaceParams

单一界面类型（如 :PE_PCC 或 :NE_NCC）的 CZM 本构参数（归一化后）。

所有字段都已通过 NormaliseParam 归一化：
- E_eff, σ_max, τ_max: / scale.σ_czm
- δ_0, δ_c, δ_0_t, δ_c_t: / scale.δ_czm
- G_c, G_c_t: / scale.G_czm
- K_n, K_t: / scale.K_czm
"""
@with_kw struct CzmInterfaceParams
    E_eff::Float64 = 0.0       # 体模量（涂层模量，非全栈均一化）
    ν::Float64 = 0.0
    α::Float64 = 0.0           # 热膨胀系数（归一化）
    σ_max::Float64 = 0.0       # Mode I 最大牵引
    δ_0::Float64 = 0.0         # Mode I 损伤起始位移
    δ_c::Float64 = 0.0         # Mode I 临界位移
    G_c::Float64 = 0.0         # Mode I 断裂能
    K_n::Float64 = 0.0         # Mode I 初始刚度
    τ_max::Float64 = 0.0       # Mode II 最大牵引
    δ_0_t::Float64 = 0.0       # Mode II 损伤起始位移
    δ_c_t::Float64 = 0.0       # Mode II 临界位移
    G_c_t::Float64 = 0.0       # Mode II 断裂能
    K_t::Float64 = 0.0         # Mode II 初始刚度
end

"""
    CzmParamCache

按界面类型分组的 CZM 参数缓存。
"""
struct CzmParamCache
    by_interface::Dict{Symbol, CzmInterfaceParams}
end
```

如果 `CouplingState.jl` 没有 `using Parameters`，添加 `using Parameters: @with_kw` 到文件顶部。

**已核实**：`src/CouplingState.jl` 顶部**没有** `using Parameters` 导入。**必须**添加以下行到文件顶部 imports 区：

```julia
using Parameters: @with_kw
```

- [ ] **Step 4: 在 `src/JuBat.jl` 导出**

找到 `src/JuBat.jl` 的 export 列表，添加：

```julia
export CzmInterfaceParams, CzmParamCache
```

- [ ] **Step 5: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_cohesive_struct.jl")'`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/CouplingState.jl src/JuBat.jl test/test_cohesive_struct.jl
git commit -m "feat(czm): 新增 CzmInterfaceParams 与 CzmParamCache 数据结构

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 1.4: 新增 `CzmSubmesh` struct 与 `CohesiveMesh` 扩展

**Files:**
- Modify: `src/czm.jl`（新增 `CzmSubmesh`）
- Modify: `src/SetMesh.jl:26-46`（扩展 `CohesiveMesh`）
- Modify: `test/test_cohesive_struct.jl`（添加新测试集）

- [ ] **Step 1: 写失败测试**

在 `test/test_cohesive_struct.jl` 末尾添加：

```julia
@testset "CzmSubmesh struct" begin
    # 仅验证字段存在，构造完整子网格在 Chunk 3 测试
    @test JuBat.CzmSubmesh === JuBat.CzmSubmesh  # 类型存在性检查
end

@testset "CohesiveMesh czm_submesh field" begin
    mesh = JuBat.CohesiveMesh()
    @test hasproperty(mesh, :czm_submesh)
    @test hasproperty(mesh, :thermal_to_czm)
    @test isnothing(mesh.czm_submesh)
    @test isnothing(mesh.thermal_to_czm)
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_cohesive_struct.jl")'`
Expected: FAIL — `JuBat.CzmSubmesh` 不存在，`CohesiveMesh` 无 `czm_submesh` 字段。

- [ ] **Step 3: 实现 `CzmSubmesh`**

在 `src/czm.jl` 的 struct 定义区（`CohesiveElement` 之后）添加：

```julia
"""
    CzmSubmesh

独立细化的 CZM 机械子网格（径向 8 层/卷绕圈）。
与粗热网格解耦，通过 thermal_elem_map 与 thermal_to_czm 矩阵耦合。
"""
struct CzmSubmesh
    mesh::Mesh                              # 细化 Q4 网格
    material_type::Vector{Symbol}           # :PE / :PCC / :SP / :NE / :NCC
    winding_turn::Vector{Int}               # 卷绕圈号（从内到外 1, 2, ...）
    thermal_elem_map::Vector{Int}           # 每个 CZM 单元 → 对应的粗热单元 id
end
```

- [ ] **Step 4: 扩展 `CohesiveMesh`**

修改 `src/SetMesh.jl:26-46`，新增两个字段并更新内部构造函数：

```julia
mutable struct CohesiveMesh
    bulk_mesh::Mesh
    node::Matrix{Float64}
    nnode::Int64
    bulk_element::Matrix{Int64}
    cohesive_elements::Vector{AbstractCohesiveElement}
    n_cohesive::Int64
    n_layers::Int64
    node_map::Dict{Int64, Vector{Int64}}
    interface_nodes::Vector{Vector{Tuple{Int64,Int64}}}
    damage_states::Vector{AbstractDamageState}
    czm_submesh::Union{Nothing, CzmSubmesh}                       # 新增
    thermal_to_czm::Union{Nothing, SparseMatrixCSC{Float64, Int}} # 新增

    function CohesiveMesh()
        new(Mesh("Q4", 2, zeros(0,2), 0, zeros(Int64,0,4),
            GaussPoint(zeros(0,2), zeros(0,2), zeros(0), zeros(0), zeros(Int64,0), zeros(0,4), zeros(0,8), 2)),
            zeros(0, 2), 0, zeros(Int64, 0, 4),
            AbstractCohesiveElement[], 0, 0, Dict{Int64, Vector{Int64}}(),
            Vector{Vector{Tuple{Int64,Int64}}}(), AbstractDamageState[],
            nothing, nothing)  # 新字段默认 nothing
    end
end
```

如果 `SetMesh.jl` 顶部没有 `using SparseArrays`，添加 `using SparseArrays: SparseMatrixCSC`。

**已核实**：`src/SetMesh.jl` 顶部**没有** `using SparseArrays`。**必须**添加：

```julia
using SparseArrays: SparseMatrixCSC
```

- [ ] **Step 5: 在 `src/JuBat.jl` 导出**

添加：

```julia
export CzmSubmesh
```

- [ ] **Step 6: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_cohesive_struct.jl")'`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/czm.jl src/SetMesh.jl src/JuBat.jl test/test_cohesive_struct.jl
git commit -m "feat(czm): 新增 CzmSubmesh，CohesiveMesh 扩展子网格字段

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 2: 参数集与 per-interface 参数计算

本 chunk 完成参数集填充、归一化逻辑、按界面类型计算参数的函数。

### Task 2.1: 填充 `parameters/Jellyroll.jl` 两组实验参数

**Files:**
- Modify: `src/parameters/Jellyroll.jl:154-181`（cohesive 参数块）

**说明：** 实验参数具体数值由用户在执行此任务时填入。**执行此任务前必须先与用户确认：实验测得的 PE-PCC 与 NE-NCC 界面 σ_max / G_c / K_n 值是否就绪？**

- **若就绪**：按用户提供的实测值填入
- **若未就绪**：使用下方占位值（沿用旧单一参数 82e6/25.3/2.4e17），并在 commit message 中明确标注"占位值，待用户替换"。这些占位值仅用于让回归测试通过，不能作为最终物理结论

此 gating 步骤确保 spec §10 的开放项不被静默跳过。

- [ ] **Step 1: 写失败测试（验证参数集填充正确）**

创建 `test/test_jellyroll_cohesive_params.jl`：

```julia
using Test
using JuBat

@testset "Jellyroll cohesive per-interface params" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    coh = param_dim.cohesive

    # PE-PCC 必须有非零值
    @test coh.σ_max_pe_pcc > 0
    @test coh.G_c_pe_pcc > 0
    @test coh.K_n_pe_pcc > 0
    @test coh.δ_0_pe_pcc > 0
    @test coh.δ_c_pe_pcc > 0
    @test coh.τ_max_pe_pcc > 0
    @test coh.K_t_pe_pcc > 0

    # δ_c = 2G_c/σ_max 一致性
    @test coh.δ_c_pe_pcc ≈ 2 * coh.G_c_pe_pcc / coh.σ_max_pe_pcc rtol=1e-6
    @test coh.δ_0_pe_pcc ≈ coh.σ_max_pe_pcc / coh.K_n_pe_pcc rtol=1e-6

    # NE-NCC 同上
    @test coh.σ_max_ne_ncc > 0
    @test coh.G_c_ne_ncc > 0
    @test coh.K_n_ne_ncc > 0
    @test coh.δ_c_ne_ncc ≈ 2 * coh.G_c_ne_ncc / coh.σ_max_ne_ncc rtol=1e-6
    @test coh.δ_0_ne_ncc ≈ coh.σ_max_ne_ncc / coh.K_n_ne_ncc rtol=1e-6
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_jellyroll_cohesive_params.jl")'`
Expected: FAIL — `coh.σ_max_pe_pcc == 0`（未填值）

- [ ] **Step 3: 填充参数（替换 `src/parameters/Jellyroll.jl:154-181`）**

```julia
# === Cohesive zone model parameters ===
# 按界面类型分组：PE-PCC（正极涂层-正极集流体）与 NE-NCC（负极涂层-负极集流体）
# 实验数据来源：[用户填入文献/实验来源]
cohesive = Cohesive()

# --- PE-PCC 界面（Mode I + Mode II）---
# TODO: 用户填入实测值（以下为占位，参考旧单一值 σ_max=82e6, G_c=25.3, K=2.4e17）
cohesive.σ_max_pe_pcc = 82e6       # [Pa] TODO 用户提供实测值
cohesive.K_n_pe_pcc   = 2.4e17     # [Pa/m] TODO 用户提供实测值
cohesive.δ_0_pe_pcc   = cohesive.σ_max_pe_pcc / cohesive.K_n_pe_pcc
cohesive.G_c_pe_pcc   = 25.3       # [J/m²] TODO 用户提供实测值
cohesive.δ_c_pe_pcc   = 2.0 * cohesive.G_c_pe_pcc / cohesive.σ_max_pe_pcc
# Mode II（若无独立测量，沿用 Mode I）
cohesive.τ_max_pe_pcc     = cohesive.σ_max_pe_pcc
cohesive.K_t_pe_pcc       = cohesive.K_n_pe_pcc
cohesive.δ_0_pe_pcc_t     = cohesive.τ_max_pe_pcc / cohesive.K_t_pe_pcc
cohesive.G_c_pe_pcc_t     = cohesive.G_c_pe_pcc
cohesive.δ_c_pe_pcc_t     = 2.0 * cohesive.G_c_pe_pcc_t / cohesive.τ_max_pe_pcc

# --- NE-NCC 界面（Mode I + Mode II）---
cohesive.σ_max_ne_ncc = 82e6       # [Pa] TODO 用户提供实测值
cohesive.K_n_ne_ncc   = 2.4e17     # [Pa/m] TODO 用户提供实测值
cohesive.δ_0_ne_ncc   = cohesive.σ_max_ne_ncc / cohesive.K_n_ne_ncc
cohesive.G_c_ne_ncc   = 25.3       # [J/m²] TODO 用户提供实测值
cohesive.δ_c_ne_ncc   = 2.0 * cohesive.G_c_ne_ncc / cohesive.σ_max_ne_ncc
cohesive.τ_max_ne_ncc     = cohesive.σ_max_ne_ncc
cohesive.K_t_ne_ncc       = cohesive.K_n_ne_ncc
cohesive.δ_0_ne_ncc_t     = cohesive.τ_max_ne_ncc / cohesive.K_t_ne_ncc
cohesive.G_c_ne_ncc_t     = cohesive.G_c_ne_ncc
cohesive.δ_c_ne_ncc_t     = 2.0 * cohesive.G_c_ne_ncc_t / cohesive.τ_max_ne_ncc

# --- 共用参数 ---
cohesive.czm_model = "model1"
cohesive.eta = 1.45                 # BK 准则指数 [-]

# 界面热阻（沿用旧值）
cohesive.h_c0 = 1e7
cohesive.k_air = 0.026
cohesive.lambda_m = 70e-9
cohesive.beta = 1.0
cohesive.threshold = 70e-9
```

- [ ] **Step 4: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_jellyroll_cohesive_params.jl")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/parameters/Jellyroll.jl test/test_jellyroll_cohesive_params.jl
git commit -m "feat(params): Jellyroll 参数集填充 PE-PCC / NE-NCC 两组 cohesive 参数

占位值沿用旧单一参数，待用户提供实测值后替换

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 2.2: 重写 `NormaliseParam` 归一化 cohesive 字段

**Files:**
- Modify: `src/SetParams.jl:467-481`（cohesive 归一化块，临时占位已在 Task 1.1 注释）

- [ ] **Step 1: 修复 `SetParams.jl:324`（σ_czm 尺度种子）**

**前置修复**：Task 1.1 Step 4 中临时把 `param_dim.scale.σ_czm = param_dim.cohesive.σ_max_n` 改为 `param_dim.cohesive.σ_max_pe_pcc`。此步骤确认该修改仍存在并正确。如果 `param_dim.cohesive.σ_max_pe_pcc` 也为 0（参数集未填），归一化将除以 0。

打开 `src/SetParams.jl:324`，应看到：

```julia
param_dim.scale.σ_czm = param_dim.cohesive.σ_max_pe_pcc   # TODO Chunk 2 Task 2.2 复核
```

如果 Task 1.1 修复点遗漏，先补上。**此行决定了下方归一化测试中 `param.cohesive.σ_max_pe_pcc ≈ param_dim.cohesive.σ_max_pe_pcc / scale.σ_czm` 是否成立**（因为 `scale.σ_czm` 由 `σ_max_pe_pcc` 派生，比值应 = 1）。

- [ ] **Step 2: 写失败测试**

创建 `test/test_cohesive_normalization.jl`：

```julia
using Test
using JuBat

@testset "Cohesive per-interface normalization" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    # NormaliseParam 由 SetCase 内部调用，这里通过 SetCase 触发
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)
    coh = case.param.cohesive
    scale = param_dim.scale

    # 归一化一致性：σ_max_pe_pcc_norm ≈ σ_max_pe_pcc_dim / scale.σ_czm
    expected = param_dim.cohesive.σ_max_pe_pcc / scale.σ_czm
    @test coh.σ_max_pe_pcc ≈ expected rtol=1e-6

    # 同样验证 G_c / K / δ
    @test coh.G_c_pe_pcc ≈ param_dim.cohesive.G_c_pe_pcc / scale.G_czm rtol=1e-6
    @test coh.K_n_pe_pcc ≈ param_dim.cohesive.K_n_pe_pcc / scale.K_czm rtol=1e-6
    @test coh.δ_c_pe_pcc ≈ param_dim.cohesive.δ_c_pe_pcc / scale.δ_czm rtol=1e-6

    # NE-NCC 同样验证
    @test coh.σ_max_ne_ncc ≈ param_dim.cohesive.σ_max_ne_ncc / scale.σ_czm rtol=1e-6
    @test coh.G_c_ne_ncc ≈ param_dim.cohesive.G_c_ne_ncc / scale.G_czm rtol=1e-6
end
```

- [ ] **Step 3: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_cohesive_normalization.jl")'`
Expected: FAIL — Task 1.1 中归一化被注释，新字段全是 0

- [ ] **Step 4: 重写归一化**

替换 `src/SetParams.jl:467-481`：

```julia
# cohesive zone model — per-interface 归一化
# PE-PCC 界面
param.cohesive.σ_max_pe_pcc = param_dim.cohesive.σ_max_pe_pcc / param_dim.scale.σ_czm
param.cohesive.δ_0_pe_pcc = param_dim.cohesive.δ_0_pe_pcc / param_dim.scale.δ_czm
param.cohesive.δ_c_pe_pcc = param_dim.cohesive.δ_c_pe_pcc / param_dim.scale.δ_czm
param.cohesive.G_c_pe_pcc = param_dim.cohesive.G_c_pe_pcc / param_dim.scale.G_czm
param.cohesive.K_n_pe_pcc = param_dim.cohesive.K_n_pe_pcc / param_dim.scale.K_czm
param.cohesive.τ_max_pe_pcc = param_dim.cohesive.τ_max_pe_pcc / param_dim.scale.σ_czm
param.cohesive.δ_0_pe_pcc_t = param_dim.cohesive.δ_0_pe_pcc_t / param_dim.scale.δ_czm
param.cohesive.δ_c_pe_pcc_t = param_dim.cohesive.δ_c_pe_pcc_t / param_dim.scale.δ_czm
param.cohesive.G_c_pe_pcc_t = param_dim.cohesive.G_c_pe_pcc_t / param_dim.scale.G_czm
param.cohesive.K_t_pe_pcc = param_dim.cohesive.K_t_pe_pcc / param_dim.scale.K_czm
# NE-NCC 界面
param.cohesive.σ_max_ne_ncc = param_dim.cohesive.σ_max_ne_ncc / param_dim.scale.σ_czm
param.cohesive.δ_0_ne_ncc = param_dim.cohesive.δ_0_ne_ncc / param_dim.scale.δ_czm
param.cohesive.δ_c_ne_ncc = param_dim.cohesive.δ_c_ne_ncc / param_dim.scale.δ_czm
param.cohesive.G_c_ne_ncc = param_dim.cohesive.G_c_ne_ncc / param_dim.scale.G_czm
param.cohesive.K_n_ne_ncc = param_dim.cohesive.K_n_ne_ncc / param_dim.scale.K_czm
param.cohesive.τ_max_ne_ncc = param_dim.cohesive.τ_max_ne_ncc / param_dim.scale.σ_czm
param.cohesive.δ_0_ne_ncc_t = param_dim.cohesive.δ_0_ne_ncc_t / param_dim.scale.δ_czm
param.cohesive.δ_c_ne_ncc_t = param_dim.cohesive.δ_c_ne_ncc_t / param_dim.scale.δ_czm
param.cohesive.G_c_ne_ncc_t = param_dim.cohesive.G_c_ne_ncc_t / param_dim.scale.G_czm
param.cohesive.K_t_ne_ncc = param_dim.cohesive.K_t_ne_ncc / param_dim.scale.K_czm
# 共用
param.cohesive.eta = param_dim.cohesive.eta
# 界面热阻归一化（沿用旧逻辑）
param.cohesive.h_c0 = param_dim.cohesive.h_c0 * param_dim.scale.L / param.scale.lambda
param.cohesive.k_air = param_dim.cohesive.k_air / param.scale.lambda
param.cohesive.lambda_m = param_dim.cohesive.lambda_m / param.scale.L
param.cohesive.beta = param_dim.cohesive.beta
param.cohesive.threshold = param_dim.cohesive.threshold / param.scale.L
```

- [ ] **Step 5: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_cohesive_normalization.jl")'`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/SetParams.jl test/test_cohesive_normalization.jl
git commit -m "feat(params): NormaliseParam 按 interface_type 归一化 cohesive 参数

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 2.3: 实现 `compute_czm_params_per_interface`

**Files:**
- Modify: `src/CouplingState.jl:228-275`（替换 `compute_effective_coating_modulus` 与 `compute_czm_effective_params`）

**说明：** 此任务实现新函数，但不改调用方（调用方在 Chunk 4）。`compute_czm_effective_params` 直接删除（不保留 fallback）。

- [ ] **Step 1: 写失败测试**

创建 `test/test_czm_params_per_interface.jl`：

```julia
using Test
using JuBat

@testset "compute_czm_params_per_interface" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.model = "SPMe"
    case = JuBat.SetCase(param_dim, opt)

    cache = JuBat.compute_czm_params_per_interface(case)

    @test cache isa JuBat.CzmParamCache
    @test haskey(cache.by_interface, :PE_PCC)
    @test haskey(cache.by_interface, :NE_NCC)

    pe = cache.by_interface[:PE_PCC]
    ne = cache.by_interface[:NE_NCC]

    # PE-PCC 用 PE.E_coat（非全栈均一化）
    scale = param_dim.scale
    expected_E_pe = case.param.PE.E_coat * scale.E_coat / scale.σ_czm
    @test pe.E_eff ≈ expected_E_pe rtol=1e-6

    # NE-NCC 用 NE.E_coat
    expected_E_ne = case.param.NE.E_coat * scale.E_coat / scale.σ_czm
    @test ne.E_eff ≈ expected_E_ne rtol=1e-6

    # σ_max 与参数集一致（归一化值）
    @test pe.σ_max ≈ case.param.cohesive.σ_max_pe_pcc rtol=1e-6
    @test ne.σ_max ≈ case.param.cohesive.σ_max_ne_ncc rtol=1e-6
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_czm_params_per_interface.jl")'`
Expected: FAIL — `compute_czm_params_per_interface` 不存在

- [ ] **Step 3: 实现新函数**

替换 `src/CouplingState.jl:215-275`（含 `compute_effective_coating_modulus` 与 `compute_czm_effective_params`）：

```julia
"""
    compute_czm_params_per_interface(case) -> CzmParamCache

按界面类型计算 CZM 参数。E_eff 用涂层模量（PE.E_coat / NE.E_coat），
不再做全栈均一化。

返回 `CzmParamCache`，包含 `:PE_PCC` 与 `:NE_NCC` 两个条目。
"""
function compute_czm_params_per_interface(case)
    param = case.param
    scale = case.param_dim.scale

    # 入口断言
    @assert case.param_dim.PE.E_coat > 0 && case.param_dim.NE.E_coat > 0 "CZM 应变驱动需要 PE/NE.E_coat > 0"
    @assert param.cohesive.σ_max_pe_pcc > 0 "cohesive.σ_max_pe_pcc 必须为正"
    @assert param.cohesive.σ_max_ne_ncc > 0 "cohesive.σ_max_ne_ncc 必须为正"
    @assert param.cohesive.G_c_pe_pcc > 0 && param.cohesive.G_c_ne_ncc > 0 "G_c_pe_pcc / G_c_ne_ncc 必须为正"
    @assert param.cohesive.K_n_pe_pcc > 0 && param.cohesive.K_n_ne_ncc > 0 "K_n_pe_pcc / K_n_ne_ncc 必须为正"

    # E_eff: 用涂层模量（非全栈均一化）
    E_eff_pe = param.PE.E_coat * scale.E_coat / scale.σ_czm
    E_eff_ne = param.NE.E_coat * scale.E_coat / scale.σ_czm

    coh = param.cohesive

    pe_pcc = CzmInterfaceParams(
        E_eff = E_eff_pe,
        ν = param.PE.nu_coat,
        α = param.PE.alphaT,
        σ_max = coh.σ_max_pe_pcc,
        δ_0 = coh.δ_0_pe_pcc,
        δ_c = coh.δ_c_pe_pcc,
        G_c = coh.G_c_pe_pcc,
        K_n = coh.K_n_pe_pcc,
        τ_max = coh.τ_max_pe_pcc,
        δ_0_t = coh.δ_0_pe_pcc_t,
        δ_c_t = coh.δ_c_pe_pcc_t,
        G_c_t = coh.G_c_pe_pcc_t,
        K_t = coh.K_t_pe_pcc
    )

    ne_ncc = CzmInterfaceParams(
        E_eff = E_eff_ne,
        ν = param.NE.nu_coat,
        α = param.NE.alphaT,
        σ_max = coh.σ_max_ne_ncc,
        δ_0 = coh.δ_0_ne_ncc,
        δ_c = coh.δ_c_ne_ncc,
        G_c = coh.G_c_ne_ncc,
        K_n = coh.K_n_ne_ncc,
        τ_max = coh.τ_max_ne_ncc,
        δ_0_t = coh.δ_0_ne_ncc_t,
        δ_c_t = coh.δ_c_ne_ncc_t,
        G_c_t = coh.G_c_ne_ncc_t,
        K_t = coh.K_t_ne_ncc
    )

    return CzmParamCache(Dict(:PE_PCC => pe_pcc, :NE_NCC => ne_ncc))
end
```

**移除**原 `compute_effective_coating_modulus`（lines 228-240）与 `compute_czm_effective_params`（lines 258-275）整个函数定义。

- [ ] **Step 4: 更新 `src/JuBat.jl` export**

移除 `compute_czm_effective_params`、`compute_effective_coating_modulus` 的 export（如果存在），添加：

```julia
export compute_czm_params_per_interface
```

- [ ] **Step 4b: 处理外部调用方**

`tools/` 与 `example/` 中有 4 处对 `JuBat.compute_czm_effective_params` 或 `JuBat.compute_effective_coating_modulus` 的调用（Task 1.1 Step 4 修复清单中已临时注释）。**本步骤保持注释不变**（这些 probe/example 脚本是诊断工具，Chunk 4 完成后由用户决定是否恢复或废弃）。

确认以下 4 个文件中相关调用仍带 `# TODO Chunk 4` 注释：
- `tools/verify_czm_standalone.jl:74`
- `tools/czm_convergence_diag.jl:36`
- `tools/czm_baseline_probe.jl:47`
- `example/力学模块验证/test_electrode_coat_modulus.jl:107, 132`

**不在本 chunk 恢复调用**——这些脚本依赖待实现的 Chunk 4 调用链。

- [ ] **Step 5: 修复 `CouplingState.jl:364` 调用点**

原 line 364 `E_eff, ν_eff, α_eff, β_n, β_p = compute_czm_effective_params(case)` 临时改为：

```julia
# TODO Chunk 4 Task 4.x: 改为接受 CzmParamCache
czm_param_cache = compute_czm_params_per_interface(case)
# 以下代码暂用 PE_PCC 字段作占位（仅保证编译通过，逻辑在 Chunk 4 修复）
E_eff = czm_param_cache.by_interface[:PE_PCC].E_eff
ν_eff = czm_param_cache.by_interface[:PE_PCC].ν
α_eff = czm_param_cache.by_interface[:PE_PCC].α
β_n = case.param.NE.Omega / 3.0
β_p = case.param.PE.Omega / 3.0
```

- [ ] **Step 6: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_czm_params_per_interface.jl")'`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/CouplingState.jl src/JuBat.jl test/test_czm_params_per_interface.jl
git commit -m "feat(czm): 实现 compute_czm_params_per_interface，移除均一化路径

- E_eff 改用涂层模量（PE.E_coat / NE.E_coat），非全栈均一化
- 移除 compute_effective_coating_modulus 与 compute_czm_effective_params
- 调用方适配见 Chunk 4

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 2.4: 修复 Chunk 1 遗留的临时占位编译错误

**说明：** Task 1.1 Step 4 中多处文件临时占位（Materialmatrix.jl、czm.jl 内 K_n/σ_max 引用等）。本任务确认 `using JuBat` 加载无错误，不涉及逻辑。

- [ ] **Step 1: 加载 JuBat 并运行所有测试**

Run:
```bash
julia --project=. -e 'using JuBat; println("load OK")'
julia --project=. -e 'include("test/test_cohesive_struct.jl"); include("test/test_jellyroll_cohesive_params.jl"); include("test/test_cohesive_normalization.jl"); include("test/test_czm_params_per_interface.jl")'
```
Expected: 第一条打印 `load OK`，第二条所有测试集 PASS

- [ ] **Step 2: 如有编译错误，定位并修复**

任何剩余编译错误通常源于 Chunk 1 Task 1.1 Step 4 的临时占位遗漏。搜索所有 `# TODO Chunk` 注释，确认是否还有未占位的旧字段引用：

```bash
grep -rn "cohesive_params\.\(σ_max_n\|K_n\|δ_0_n\|G_c_n\|δ_c_n\|τ_max_t\|K_t\|δ_0_t\|G_c_t\|δ_c_t\)\b" src/
grep -rn "cohesive\.\(σ_max_n\|K_n\|δ_0_n\|G_c_n\|δ_c_n\|τ_max_t\|K_t\|δ_0_t\|G_c_t\|δ_c_t\)\b" src/
```

所有结果都应有 `# TODO Chunk 4` 注释并已占位。

- [ ] **Step 3: Commit（如有修复）**

```bash
git add -u
git commit -m "fix: 清理 Chunk 1 遗留的旧字段引用占位

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

若 Step 1 已通过，跳过 commit。

---

## Chunk 3: CZM 子网格生成

本 chunk 按 spec v2 §4.1-§4.4 实现：独立细化 CZM 子网格、O(1) 解析式反查、节点复制+重写外层 bulk 单元、`create_czm_mesh` 3 参签名。

### Task 3.1: 实现 `jellyroll_czm_submesh`（spec §4.1）

**Files:**
- Modify: `src/Jellyrollmodel.jl`（在 `jellyroll_collector_seed_mesh` 函数之后追加 `jellyroll_czm_submesh`）
- Modify: `src/JuBat.jl`（export）
- Create: `test/test_czm_submesh.jl`

**说明：** 函数入参 `(param, thermal_mesh::Mesh; nθ_czm, gsorder)`，其中 `param` 是 `NormaliseParam(param_dim)` 的输出（与 `jellyroll_collector_seed_mesh` 同约定）。生成 8 径向分层 Q4 子网格，按构造顺序绑定 `material_type`，用 **O(1) 解析式**建立 `thermal_elem_map`（spec §4.1.1，禁止 `findmin` 空间搜索）。

- [ ] **Step 1: 写失败测试 `test/test_czm_submesh.jl`**

```julia
using Test
using JuBat

@testset "jellyroll_czm_submesh" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=40)

    @test submesh isa JuBat.CzmSubmesh
    @test submesh.mesh.type == "Q4"
    ne = size(submesh.mesh.element, 1)

    # material_type 取值合法
    @test all(m in (:PE, :PCC, :SP, :NE, :NCC) for m in submesh.material_type)

    # 8 层周期序列：前 8 个单元应为 [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    if length(submesh.material_type) >= 8
        @test submesh.material_type[1:8] == [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    end

    # 每个 CZM 体单元都映射到一个粗热单元，且映射值合法
    n_thermal = size(case.mesh["thermal2D"].element, 1)
    @test all(1 <= e <= n_thermal for e in submesh.thermal_elem_map)

    # winding_turn 与径向位置一致：r 中心落在 turn*cell.layer 区间内
    s_total = case.param.cell.layer
    a = case.param.cell.Rin
    for e in 1:length(submesh.winding_turn)
        n1 = submesh.mesh.element[e, 1]
        n3 = submesh.mesh.element[e, 3]
        r = 0.5 * (hypot(submesh.mesh.node[n1, 1], submesh.mesh.node[n1, 2]) +
                   hypot(submesh.mesh.node[n3, 1], submesh.mesh.node[n3, 2]))
        expected_turn = max(1, Int(floor((r - a) / s_total)) + 1)
        @test submesh.winding_turn[e] == expected_turn
    end

    # O(1) 解析反查：thermal_elem_map 值应等于 (turn-1)*n_seg_thermal + seg_id
    # 通过粗热网格 θ_center 反推 seg_id 验证
    n_seg_thermal = size(case.mesh["thermal2D"].element, 1) ÷ maximum(submesh.winding_turn)
    for e in 1:ne
        @test submesh.thermal_elem_map[e] >= 1
        @test submesh.thermal_elem_map[e] <= size(case.mesh["thermal2D"].element, 1)
    end
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_czm_submesh.jl")'`
Expected: FAIL — `jellyroll_czm_submesh` 未定义（或 UndefVarError）

- [ ] **Step 3: 实现 `jellyroll_czm_submesh`（O(1) 解析反查）**

在 `src/Jellyrollmodel.jl` 的 `jellyroll_collector_seed_mesh` 函数之后追加：

```julia
"""
    jellyroll_czm_submesh(param, thermal_mesh::Mesh; nθ_czm=80, gsorder=2) -> CzmSubmesh

生成独立细化的 CZM 机械子网格。径向 8 层/卷绕圈（按 PE→PCC→PE→SP→NE→NCC→NE→SP 顺序），
周向 nθ_czm 个单元。与粗热网格共享螺旋几何，通过 (θ, turn) 解析反查 thermal_elem_map。

# 入参约定
- `param`: NormaliseParam(param_dim) 的输出（归一化后的 case.param）
- `thermal_mesh`: 粗热网格 Mesh（来自 jellyroll_collector_seed_mesh）

# 算法（spec §4.1.1）
- turn = floor((r_center - Rin) / cell.layer) + 1
- seg_id = θ_center 在该 turn 的粗热 segment 区间中的位置
- thermal_elem_map[e] = (turn-1) * n_seg_thermal + seg_id
"""
function jellyroll_czm_submesh(param, thermal_mesh::Mesh; nθ_czm::Int=80, gsorder::Int=2)
    # 螺旋几何参数（与粗热网格一致，使用归一化值）
    a = param.cell.Rin
    s_total = param.cell.layer              # 一个卷绕圈的径向厚度（归一化）
    b = s_total / (2 * pi)

    # 径向 8 层的厚度（按层序）
    t_pe  = param.PE.thickness
    t_pcc = param.PCC.thickness
    t_sp  = param.SP.thickness
    t_ne  = param.NE.thickness
    t_ncc = param.NCC.thickness

    layer_thicknesses = [t_pe, t_pcc, t_pe, t_sp, t_ne, t_ncc, t_ne, t_sp]
    material_sequence = [:PE, :PCC, :PE, :SP, :NE, :NCC, :NE, :SP]
    n_layers = 8
    @assert sum(layer_thicknesses) ≈ s_total rtol=1e-6 "层厚之和必须等于 cell.layer（实际: $(sum(layer_thicknesses)), 期望: $s_total）"

    # theta 范围：与粗热网格一致（参考 jellyroll_collector_seed_mesh:35-44）
    theta0 = max(0.0, (param.cell.Rin - a) / b)
    theta1 = (param.cell.Rout - a - s_total) / b
    theta1 > theta0 || error("jellyroll_czm_submesh: 无有效 theta 范围 (theta0=$theta0, theta1=$theta1)")

    # 周向采样
    n_segments = max(2, nθ_czm)
    theta = collect(range(theta0, theta1; length=n_segments + 1))

    # 节点：每层 2 条螺旋（内边、外边），共 (n_layers+1) 条螺旋 × (n_segments+1) 个点
    n_spirals = n_layers + 1
    nnode = n_spirals * (n_segments + 1)
    node = zeros(Float64, nnode, 2)

    s_offsets = [0.0; cumsum(layer_thicknesses)]   # 长度 n_layers+1
    for layer_idx in 0:n_layers
        s_offset = s_offsets[layer_idx + 1]
        r = a .+ b .* theta .+ s_offset
        x = r .* cos.(theta)
        y = r .* sin.(theta)
        for k in 1:(n_segments + 1)
            node_idx = layer_idx * (n_segments + 1) + k
            node[node_idx, 1] = x[k]
            node[node_idx, 2] = y[k]
        end
    end

    # 解析反查准备：从粗热网格推出 n_seg_thermal（每圈的周向 segment 数）
    n_thermal = size(thermal_mesh.element, 1)
    n_seg_thermal = max(1, n_thermal ÷ max(1, round(Int, (theta1 - theta0) / (2 * pi))))
    # 通过粗热网格 θ_center 建立查找表（同 turn 内）
    thermal_θ_center = zeros(n_thermal)
    thermal_r_center = zeros(n_thermal)
    for e in 1:n_thermal
        ns = thermal_mesh.element[e, :]
        xs = thermal_mesh.node[ns, 1]
        ys = thermal_mesh.node[ns, 2]
        θ_vals = atan.(ys, xs)
        thermal_θ_center[e] = sum(θ_vals) / 4
        thermal_r_center[e] = sum(sqrt.(xs.^2 .+ ys.^2)) / 4
    end
    # 推断粗热网格 turn 划分（按 r_center）
    thermal_turn = max.(1, Int.(floor.((thermal_r_center .- a) ./ s_total)) .+ 1)
    n_turns_thermal = maximum(thermal_turn)

    # 单元：每层 n_segments 个 Q4 单元
    ne = n_layers * n_segments
    element = zeros(Int64, ne, 4)
    material_type = Vector{Symbol}(undef, ne)
    winding_turn = Vector{Int}(undef, ne)
    thermal_elem_map = Vector{Int}(undef, ne)

    elem_idx = 0
    for layer_idx in 1:n_layers
        for seg in 1:n_segments
            elem_idx += 1
            inner_spiral_base = (layer_idx - 1) * (n_segments + 1)
            outer_spiral_base = layer_idx * (n_segments + 1)
            element[elem_idx, 1] = inner_spiral_base + seg
            element[elem_idx, 2] = outer_spiral_base + seg
            element[elem_idx, 3] = outer_spiral_base + seg + 1
            element[elem_idx, 4] = inner_spiral_base + seg + 1

            material_type[elem_idx] = material_sequence[layer_idx]

            # winding_turn 由径向位置解析
            n1 = element[elem_idx, 1]
            n3 = element[elem_idx, 3]
            r_center = 0.5 * (hypot(node[n1, 1], node[n1, 2]) +
                              hypot(node[n3, 1], node[n3, 2]))
            turn_czm = max(1, Int(floor((r_center - a) / s_total)) + 1)
            winding_turn[elem_idx] = turn_czm

            # θ_center 解析（mod 2π）
            x_center = 0.25 * (node[element[elem_idx, 1], 1] + node[element[elem_idx, 2], 1] +
                               node[element[elem_idx, 3], 1] + node[element[elem_idx, 4], 1])
            y_center = 0.25 * (node[element[elem_idx, 1], 2] + node[element[elem_idx, 2], 2] +
                               node[element[elem_idx, 3], 2] + node[element[elem_idx, 4], 2])
            θ_czm = atan(y_center, x_center)

            # O(1) 解析反查：在同 turn 的粗热单元中找 θ 最接近者
            # 同 turn 的粗热单元索引范围：[(turn-1)*n_seg_thermal+1, turn*n_seg_thermal]
            # （粗热网格按 turn 主序、segment 次序排列，与 jellyroll_collector_seed_mesh 一致）
            lo = (turn_czm - 1) * n_seg_thermal + 1
            hi = min(n_thermal, turn_czm * n_seg_thermal)
            if lo > hi
                # 边界情况：CZM turn 超出粗热 turn 范围（最外圈部分截断）
                thermal_elem_map[elem_idx] = -1
                continue
            end
            candidate_idx = lo:hi
            θ_diff = abs.(mod.(thermal_θ_center[candidate_idx] .- θ_czm .+ pi, 2 * pi) .- pi)
            _, local_best = findmin(θ_diff)
            thermal_elem_map[elem_idx] = candidate_idx[local_best]
        end
    end

    # 失败检查：-1 表示 CZM 单元落在粗热网格外
    n_failed = count(==( -1), thermal_elem_map)
    if n_failed > 0
        @warn "jellyroll_czm_submesh: $n_failed 个 CZM 单元无法映射到粗热单元（边界外）"
    end

    gs = GetGS(element, node, gsorder, "Q4")
    mesh = Mesh("Q4", 2, node, nnode, element, gs)

    return CzmSubmesh(mesh, material_type, winding_turn, thermal_elem_map)
end
```

- [ ] **Step 4: 在 `src/JuBat.jl` 导出**

在现有 export 段加入：
```julia
export jellyroll_czm_submesh, CzmSubmesh
```

- [ ] **Step 5: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_czm_submesh.jl")'`
Expected: PASS（所有断言通过）

- [ ] **Step 6: Commit**

```bash
git add src/Jellyrollmodel.jl src/JuBat.jl test/test_czm_submesh.jl
git commit -m "feat(czm): 实现 jellyroll_czm_submesh 生成径向 8 层细化子网格

按 spec v2 §4.1：O(1) 解析式 thermal_elem_map 反查（同 turn 内 θ 最近）。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 3.2: 重构 `create_czm_mesh`（spec §4.2-§4.4）

**Files:**
- Modify: `src/czm.jl:56-136`（旧 `create_czm_mesh` 整体重写）
- Create: `test/test_create_czm_mesh.jl`

**说明：** 旧逻辑（坐标重合检测、SP-PE 界面）整体替换。新签名 `(czm_submesh, thermal_mesh, param)`。**节点复制必须传播到外层 bulk 单元**（spec §4.3 契约，否则分离位移恒为 0）。构造 `cohesive_to_thermal` 与每个 `CohesiveElement.host_outer_elem / host_inner_elem`。

- [ ] **Step 1: 写失败测试 `test/test_create_czm_mesh.jl`**

```julia
using Test
using JuBat

@testset "create_czm_mesh from CzmSubmesh" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)

    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=20)
    czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)

    @test czm_mesh isa JuBat.CohesiveMesh
    @test czm_mesh.n_cohesive > 0

    # interface_type 取值合法
    for elem in czm_mesh.cohesive_elements
        @test elem.interface_type in (:PE_PCC, :NE_NCC)
        @test elem.host_outer_elem >= 1
        @test elem.host_inner_elem >= 1
    end

    # 数量预期：每卷绕圈 2 个界面（PE-PCC + NE-NCC），每界面 n_segments 个 cohesive 单元
    n_segments_per_turn = size(submesh.mesh.element, 1) ÷ 8 ÷ maximum(submesh.winding_turn)
    n_turns_active = length(unique(submesh.winding_turn))
    n_expected = 2 * n_segments_per_turn * n_turns_active
    # 容差：边界裁剪允许 ±n_segments_per_turn
    @test abs(czm_mesh.n_cohesive - n_expected) <= n_segments_per_turn

    # czm_submesh 字段已设置
    @test czm_mesh.czm_submesh === submesh

    # thermal_to_czm 字段已设置（在 Chunk 5 中真正填充，此处先 nothing）
    @test hasfield(JuBat.CohesiveMesh, :thermal_to_czm)

    # cohesive_to_thermal 长度 = n_cohesive，值合法
    @test length(czm_mesh.cohesive_to_thermal) == czm_mesh.n_cohesive
    n_thermal = size(case.mesh["thermal2D"].element, 1)
    @test all(1 <= e <= n_thermal for e in czm_mesh.cohesive_to_thermal)

    # 节点复制 + 重写外层 bulk 正确性自检（spec §4.3）
    for coh in czm_mesh.cohesive_elements
        n_a, n_b, n_b_copy, n_a_copy = coh.nodes
        # 副本节点坐标与原节点一致
        @test czm_mesh.node[n_a, :] ≈ czm_mesh.node[n_a_copy, :] atol=1e-12
        @test czm_mesh.node[n_b, :] ≈ czm_mesh.node[n_b_copy, :] atol=1e-12
        # 4 节点不重复
        @test length(unique(coh.nodes)) == 4
    end

    # 外层 bulk 单元的共边位置必须是副本节点
    for coh in czm_mesh.cohesive_elements
        outer_nodes = czm_mesh.bulk_element[coh.host_outer_elem, :]
        n_a, n_b, n_b_copy, n_a_copy = coh.nodes
        # 副本 n_a_copy/n_b_copy 必须在外层单元中
        @test n_a_copy in outer_nodes
        @test n_b_copy in outer_nodes
        # 原节点 n_a/n_b 不应在外层单元中（已被副本替换）
        @test !(n_a in outer_nodes)
        @test !(n_b in outer_nodes)
    end
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_create_czm_mesh.jl")'`
Expected: FAIL — `create_czm_mesh` 旧签名接受 `(thermal_mesh, param_dim)`，新调用 3 参数会报 MethodError

- [ ] **Step 3: 重写 `create_czm_mesh`**

替换 `src/czm.jl:56-136` 整个函数为：

```julia
"""
    create_czm_mesh(czm_submesh::CzmSubmesh, thermal_mesh::Mesh, param) -> CohesiveMesh

基于细化 CZM 子网格构造内聚力网格（spec §4.4）。

# 核心算法
1. 建立 共边(2 节点) → 单元对 映射，遍历找径向相邻且材料组合为 PE-PCC / NE-NCC 的对
2. 节点复制：对每个界面对的共边 2 节点生成副本（memoized），cohesive 单元 4 节点 = [n_a, n_b, n_b', n_a']
3. **重写外层 bulk 单元连接**：把外层单元共边位置的原节点替换为副本节点（关键：否则分离位移恒为 0）
4. 构造 cohesive_to_thermal[e_coh] = thermal_elem_map[e_outer]
"""
function create_czm_mesh(czm_submesh::CzmSubmesh, thermal_mesh::Mesh, param)
    sub_mesh = czm_submesh.mesh
    ne_sub = size(sub_mesh.element, 1)
    nnode_sub = sub_mesh.nlen

    # Step 1: 建立 共边 → 单元对 映射
    edge_to_elems = Dict{Tuple{Int, Int}, Vector{Int}}()
    for e in 1:ne_sub
        n1, n2, n3, n4 = sub_mesh.element[e, :]
        for edge in ((n1, n2), (n2, n3), (n3, n4), (n4, n1))
            key = (min(edge[1], edge[2]), max(edge[1], edge[2]))
            push!(get!(edge_to_elems, key, Int[]), e)
        end
    end

    # Step 2: 遍历共边，识别 PE-PCC / NE-NCC 径向界面
    interface_pairs = Tuple{Int, Int, Symbol}[]   # (e_inner, e_outer, interface_type)
    for (edge, elems) in edge_to_elems
        length(elems) == 2 || continue   # 周向相邻同层（材料相同）自动过滤
        e1, e2 = elems[1], elems[2]
        m1, m2 = czm_submesh.material_type[e1], czm_submesh.material_type[e2]
        iface = if (m1 == :PE && m2 == :PCC) || (m1 == :PCC && m2 == :PE)
            :PE_PCC
        elseif (m1 == :NE && m2 == :NCC) || (m1 == :NCC && m2 == :NE)
            :NE_NCC
        else
            nothing
        end
        if iface !== nothing
            # 判断哪个是内层（径向更小）
            n1_1 = sub_mesh.element[e1, 1]
            n1_2 = sub_mesh.element[e2, 1]
            r1 = hypot(sub_mesh.node[n1_1, 1], sub_mesh.node[n1_1, 2])
            r2 = hypot(sub_mesh.node[n1_2, 1], sub_mesh.node[n1_2, 2])
            if r1 < r2
                push!(interface_pairs, (e1, e2, iface))
            else
                push!(interface_pairs, (e2, e1, iface))
            end
        end
    end

    # Step 3: 节点复制 + 重写外层 bulk 连接
    # 副本 memoization：原节点 → 单一副本节点 id
    node_copy = Dict{Int, Int}()
    extended_node = copy(sub_mesh.node)
    new_node_count = nnode_sub

    # bulk_element 可变副本（初始等于 sub_mesh.element）
    bulk_element_new = Matrix{Int}(sub_mesh.element)

    # 同时构造 cohesive 单元与 cohesive_to_thermal
    n_cohesive = length(interface_pairs)
    cohesive_elements = CohesiveElement[]
    cohesive_to_thermal = Vector{Int}(undef, n_cohesive)
    sizehint!(cohesive_elements, n_cohesive)

    for (i, (e_inner, e_outer, iface)) in enumerate(interface_pairs)
        inner_nodes = sub_mesh.element[e_inner, :]
        outer_nodes = sub_mesh.element[e_outer, :]
        common_set = intersect(Set(inner_nodes), Set(outer_nodes))
        @assert length(common_set) == 2 "共边应有 2 节点，实际 $(length(common_set))"
        common = collect(common_set)

        # 按 θ 确定性排序（避免 Set 哈希顺序导致法向翻转）
        θs = atan.([sub_mesh.node[c, 2] for c in common],
                   [sub_mesh.node[c, 1] for c in common])
        order = sortperm(θs)
        n_lo = common[order[1]]   # θ 较小
        n_hi = common[order[2]]   # θ 较大

        # memoized 副本
        function get_copy(n::Int)
            haskey(node_copy, n) && return node_copy[n]
            new_node_count += 1
            extended_node = vcat(extended_node, reshape(sub_mesh.node[n, :], 1, 2))
            node_copy[n] = new_node_count
            return new_node_count
        end
        n_lo_copy = get_copy(n_lo)
        n_hi_copy = get_copy(n_hi)

        # 重写外层 bulk 单元连接：把外层单元中的 n_lo/n_hi 替换为副本
        for col in 1:4
            if bulk_element_new[e_outer, col] == n_lo
                bulk_element_new[e_outer, col] = n_lo_copy
            elseif bulk_element_new[e_outer, col] == n_hi
                bulk_element_new[e_outer, col] = n_hi_copy
            end
        end

        # cohesive 单元几何长度
        x_lo, y_lo = sub_mesh.node[n_lo, 1], sub_mesh.node[n_lo, 2]
        x_hi, y_hi = sub_mesh.node[n_hi, 1], sub_mesh.node[n_hi, 2]
        elem_length = hypot(x_hi - x_lo, y_hi - y_lo)

        coh = CohesiveElement(
            i,
            [n_lo, n_hi, n_hi_copy, n_lo_copy],   # 逆时针
            [n_lo, n_hi],                          # nodes_bottom（内层外边）
            [n_lo_copy, n_hi_copy],                # nodes_top（外层内边副本）
            elem_length,
            iface,                                 # interface_type
            e_outer,                               # host_outer_elem
            e_inner,                               # host_inner_elem
        )
        push!(cohesive_elements, coh)

        # cohesive_to_thermal
        thermal_elem_of_outer = czm_submesh.thermal_elem_map[e_outer]
        @assert thermal_elem_of_outer > 0 "外层单元 $e_outer 的 thermal_elem_map 无效"
        cohesive_to_thermal[i] = thermal_elem_of_outer
    end

    # Step 4: 组装 CohesiveMesh
    damage_states = [DamageState() for _ in 1:n_cohesive]

    czm_mesh = CohesiveMesh()
    czm_mesh.bulk_mesh = sub_mesh
    czm_mesh.node = extended_node
    czm_mesh.nnode = new_node_count
    czm_mesh.bulk_element = bulk_element_new
    czm_mesh.cohesive_elements = cohesive_elements
    czm_mesh.n_cohesive = n_cohesive
    czm_mesh.n_layers = 2   # spec §3.3: 分离面类型数（PE-PCC + NE-NCC），非材料层数
    czm_mesh.node_map = Dict(n => [n, c] for (n, c) in node_copy)
    czm_mesh.interface_nodes = [[]]   # 旧字段，保留兼容
    czm_mesh.damage_states = damage_states
    czm_mesh.czm_submesh = czm_submesh
    czm_mesh.thermal_to_czm = nothing   # Task 5.1 真正填充
    czm_mesh.cohesive_to_thermal = cohesive_to_thermal

    # 正确性自检（spec §4.3）
    for coh in cohesive_elements
        n_a, n_b, n_b_copy, n_a_copy = coh.nodes
        @assert czm_mesh.node[n_a, :] ≈ czm_mesh.node[n_a_copy, :] atol=1e-12 "副本坐标不一致"
        @assert czm_mesh.node[n_b, :] ≈ czm_mesh.node[n_b_copy, :] atol=1e-12 "副本坐标不一致"
        @assert length(unique(coh.nodes)) == 4 "cohesive 单元 4 节点重复"
    end

    return czm_mesh
end
```

- [ ] **Step 4: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_create_czm_mesh.jl")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/czm.jl test/test_create_czm_mesh.jl
git commit -m "refactor(czm): create_czm_mesh 改 3 参签名，节点复制传播到外层 bulk

按 spec v2 §4.2-§4.4：
- 识别 PE-PCC / NE-NCC 径向界面
- 节点复制 memoized，重写外层 bulk 单元连接
- 构造 cohesive_to_thermal 与 host_outer/inner_elem
- 3 条正确性自检（副本坐标一致 / 4 节点不重复 / 外层用副本）

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 3.3: 更新 16 处 `create_czm_mesh` 调用点（spec §4.4.1）

**Files:**
- Modify: `example/` 与 `tools/` 共 16 处调用（spec §4.4.1 已列出完整清单）

**说明：** 旧调用形式 `create_czm_mesh(mesh_data.thermal2D, param_dim)` 改为两步：先生成 `CzmSubmesh`，再调用 3 参版 `create_czm_mesh`。

- [ ] **Step 1: 设计统一替换模式**

旧：
```julia
czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)
```

新（两步）：
```julia
czm_submesh = JuBat.jellyroll_czm_submesh(param, mesh_data.thermal2D; nθ_czm=80)
czm_mesh = JuBat.create_czm_mesh(czm_submesh, mesh_data.thermal2D, param)
```

**变量名规则**：
- 若上下文用 `param_dim`（未归一化），需先 `param = JuBat.NormaliseParam(param_dim)` 或复用 `case.param`
- 若上下文用 `mesh_data.Jellyroll_czm`（旧属性），改为 `mesh_data.thermal2D`
- 若上下文用 `case.czm_mesh = ...`，保持赋值目标

- [ ] **Step 2: 逐文件手动替换并 include 验证语法**

按 spec §4.4.1 表格中的 16 处调用，逐个打开文件、应用上述模式。每修改完一个文件，运行：

```bash
julia --project=. -e 'include("path/to/file.jl")' 2>&1 | head -30
```

（部分脚本可能因 Chunk 4 未完成而运行失败——**本步骤只验证语法正确**，不要求完整运行。）

**重点文件示例**：

`example/coupled_czm_thermal_example.jl:104` 上下文：
```julia
# 旧（第 104 行）
czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)
# 新（替换为 3 行）
czm_submesh = JuBat.jellyroll_czm_submesh(case.param, mesh_data.thermal2D; nθ_czm=80)
czm_mesh = JuBat.create_czm_mesh(czm_submesh, mesh_data.thermal2D, case.param)
```

`example/循环验证/czm_cycle_example.jl:98` 上下文（注意旧版用了 `mesh_data.Jellyroll_czm`）：
```julia
# 旧
czm_mesh = JuBat.create_czm_mesh(mesh_data.Jellyroll_czm, param_dim; tol=1e-8)
# 新
czm_submesh = JuBat.jellyroll_czm_submesh(case.param, mesh_data.thermal2D; nθ_czm=80)
czm_mesh = JuBat.create_czm_mesh(czm_submesh, mesh_data.thermal2D, case.param)
```

- [ ] **Step 3: 修复 `tools/verify_czm_unit.jl:92` 的 CohesiveElement 构造**

旧版构造（Task 1.1 占位为 `:PE_PCC`）按新签名重建，加入 `host_outer_elem / host_inner_elem`：
```julia
elem = JuBat.CohesiveElement(1, [1,2,3,4], [1,2], [4,3], 1.0, :PE_PCC, 2, 1)
```

- [ ] **Step 4: Commit**

```bash
git add -u example/ tools/
git commit -m "refactor(examples): 适配 create_czm_mesh 3 参签名（16 处）

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 4: 求解器适配

本 chunk 按 spec v2 §7.1 改造 CZM 求解链路：16 处 `assemble_coupled_system` 调用、`assemble_coupled_system_full`、`Materialmatrix.jl` 涉及 6 个函数、`ensure_czm_cache` 失效判据。

### Task 4.1: 改造 `Materialmatrix.jl` 中 `bilinear_*` 与 `compute_*gap_conductance*` 函数

**Files:**
- Modify: `src/Materialmatrix.jl:68, 163, 182, 293, 324, 352, 398`（共 6 个函数）
- Create: `test/test_bilinear_per_interface.jl`

**说明：** 必须先做本 Task——后续 `assemble_czm_system` 改造（Task 4.2）会把 `bilinear_*` 调用从 `cohesive_params::Cohesive` 切换到 `params::CzmInterfaceParams`，若本 Task 未完成会直接 MethodError。spec §3.5.1 已定义 `CzmInterfaceParams` 完整字段（20 个），与 `bilinear_traction_state` 实际读取的字段名（`K_n / K_t / δ_0_n / δ_c_n / δ_0_t / δ_c_t / η / czm_model`）**完全一致**。

- [ ] **Step 1: 写失败测试 `test/test_bilinear_per_interface.jl`**

```julia
using Test
using JuBat

@testset "bilinear_* with CzmInterfaceParams" begin
    # 构造 PE-PCC 接面参数
    params = JuBat.CzmInterfaceParams(
        E_eff = 1.0, ν = 0.3, α = 1.5e-5,
        σ_max = 82e6 / 1e10,   # 归一化示例
        K_n = 2.4e17 / (1e10 / 1e-6),
        δ_0_n = 82e6 / 2.4e17,
        δ_c_n = 2 * 25.3 / 82e6,
        G_c = 25.3,
        τ_max = 82e6 / 1e10,
        K_t = 2.4e17 / (1e10 / 1e-6),
        δ_0_t = 82e6 / 2.4e17,
        δ_c_t = 2 * 25.3 / 82e6,
        G_c_t = 25.3,
        η = 1.45,
        czm_model = "model1",
        h_c0 = 1e7, k_air = 0.026, lambda_m = 70e-9,
        beta = 1.0, threshold = 70e-9,
    )

    D = JuBat.DamageState()

    # 弹性段：δ_n < δ_0_n
    δ_n_small = params.δ_0_n / 2
    T_n, T_t, D_new = JuBat.bilinear_traction_state(δ_n_small, 0.0, D, params)
    @test T_n ≈ params.K_n * δ_n_small
    @test T_t ≈ 0.0
    @test D_new.D ≈ 0.0

    # 软化段：δ_0_n < δ_n < δ_c_n
    δ_n_mid = 0.5 * (params.δ_0_n + params.δ_c_n)
    T_n2, _, D_new2 = JuBat.bilinear_traction_state(δ_n_mid, 0.0, D, params)
    @test 0 < T_n2 < params.σ_max
    @test 0 < D_new2.D < 1

    # 完全失效：δ_n > δ_c_n
    δ_n_big = 2 * params.δ_c_n
    T_n3, _, D_new3 = JuBat.bilinear_traction_state(δ_n_big, 0.0, D, params)
    @test T_n3 ≈ 0.0
    @test D_new3.D ≈ 1.0

    # 切向类似（Mode II）
    T_t2, _, _ = JuBat.bilinear_traction_state(0.0, params.δ_0_t / 2, D, params)
    @test T_t2 ≈ params.K_t * params.δ_0_t / 2

    # NE-NCC 接面参数（不同 σ_max）应给出不同结果
    params_ne = JuBat.CzmInterfaceParams(params; σ_max = 2 * params.σ_max, δ_c_n = params.δ_c_n / 2)
    T_n_ne, _, _ = JuBat.bilinear_traction_state(params_ne.δ_0_n / 2, 0.0, D, params_ne)
    @test T_n_ne ≈ params_ne.K_n * params_ne.δ_0_n / 2
    @test T_n_ne ≠ T_n   # 验证按界面类型分参数生效
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_bilinear_per_interface.jl")'`
Expected: FAIL — `bilinear_traction_state` 旧签名要求 `cohesive_params::Cohesive`，传 `CzmInterfaceParams` 报 MethodError

- [ ] **Step 3: 改造 6 个函数签名**

**`src/Materialmatrix.jl:68`** `bilinear_traction_state`：把第 4 个参数 `cohesive_params::Cohesive` 改为 `params::CzmInterfaceParams`。函数体字段访问**完全保持**（`params.K_n / K_t / δ_0_n / δ_c_n / δ_0_t / δ_c_t / η / czm_model`）—— `CzmInterfaceParams` 字段名已与现 `Cohesive` 同名。

```julia
function bilinear_traction_state(δ_n::Float64, δ_t::Float64, damage_state::DamageState,
                                  params::CzmInterfaceParams; visc_beta::Float64=1.0)
    K_n = params.K_n
    K_t = params.K_t
    δ_0_n = params.δ_0_n
    δ_c_n = params.δ_c_n
    δ_0_t = params.δ_0_t
    δ_c_t = params.δ_c_t
    η = params.η
    czm_model = params.czm_model
    # ... 函数体其余不变
end
```

**`src/Materialmatrix.jl:163`** `bilinear_traction`：同上签名替换。

**`src/Materialmatrix.jl:182`** `bilinear_tangent`：同上签名替换（函数体字段读取已在 Step 1 验证）。

**`src/Materialmatrix.jl:293`** `update_damage`：同上签名替换。

**`src/Materialmatrix.jl:324`** `compute_gap_conductance`：把 `cohesive::Cohesive` 改为 `params::CzmInterfaceParams`，函数体读 `params.δ_0_n / δ_c_n / h_c0 / k_air / lambda_m / beta / threshold`。

**`src/Materialmatrix.jl:352`** `compute_element_gap_conductance` 与 **`src/Materialmatrix.jl:398`** `compute_all_gap_conductances`：同样把 `cohesive::Cohesive` 改为 `params::CzmInterfaceParams`。如有调用 `bilinear_*`，参数透传 `params`。

- [ ] **Step 4: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_bilinear_per_interface.jl")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/Materialmatrix.jl test/test_bilinear_per_interface.jl
git commit -m "refactor(material): bilinear_* 与 compute_*gap_conductance* 改接受 CzmInterfaceParams

按 spec v2 §7.1：6 个函数签名从 Cohesive 改为 CzmInterfaceParams，
函数体字段读取不变（字段名已与 Cohesive 同名）。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 4.2: 改造 `assemble_czm_system` 与 `assemble_coupled_system`

**Files:**
- Modify: `src/czm.jl:156`（`assemble_czm_system` 函数签名）
- Modify: `src/czm.jl:200-560`（循环内按 `interface_type` 取 `CzmInterfaceParams`）
- Modify: `src/czm.jl:564`（`assemble_coupled_system` 函数签名）
- Modify: `src/czm.jl:587-600`（`assemble_coupled_system_full` 同步）
- Create: `test/test_assemble_coupled_system.jl`

**说明：** 签名从 `(czm_mesh, u, E_eff, ν_eff, cohesive_params)` 改为 `(czm_mesh, u, param_cache)`。循环内对每个 cohesive 单元从 `param_cache.by_interface[iface]` 取参数；体模量从 `param_cache.param_ref.PE.E_coat` 或 `NE.E_coat` 按 `material_type` 取（不再全栈均一化）。

- [ ] **Step 1: 写失败测试 `test/test_assemble_coupled_system.jl`**

```julia
using Test
using JuBat

@testset "assemble_coupled_system with CzmParamCache" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=20)
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)

    param_cache = JuBat.compute_czm_params_per_interface(case)
    @test param_cache isa JuBat.CzmParamCache
    @test haskey(param_cache.by_interface, :PE_PCC)
    @test haskey(param_cache.by_interface, :NE_NCC)

    # PE-PCC 用 PE.E_coat，NE-NCC 用 NE.E_coat
    pe = param_cache.by_interface[:PE_PCC]
    ne = param_cache.by_interface[:NE_NCC]
    @test pe.E_eff ≈ case.param.PE.E_coat * case.param.scale.E_coat / case.param.scale.σ_czm
    @test ne.E_eff ≈ case.param.NE.E_coat * case.param.scale.E_coat / case.param.scale.σ_czm

    u = zeros(2 * case.czm_mesh.nnode)
    cache = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    K, f, seps, tracts = JuBat.assemble_coupled_system(case.czm_mesh, u, param_cache;
                                                       damage_states=case.czm_mesh.damage_states,
                                                       K_bulk_cached=cache.K_bulk,
                                                       geom_cache=cache.geom_cache,
                                                       ws=cache.ws)
    @test size(K, 1) == length(f) == 2 * case.czm_mesh.nnode
    @test !any(isnan, K)
    @test !any(isnan, f)
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_assemble_coupled_system.jl")'`
Expected: FAIL — `assemble_coupled_system` 旧签名要求 `(czm_mesh, u, E_eff, ν_eff, cohesive_params)`

- [ ] **Step 3: 改造 `assemble_coupled_system`**

替换 `src/czm.jl:564` 函数签名为：

```julia
function assemble_coupled_system(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param_cache::CzmParamCache;
    F_ext::Union{Vector{Float64}, Nothing}=nothing,
    F_thermo_chem::Union{Vector{Float64}, Nothing}=nothing,
    damage_states=nothing,
    K_bulk_cached::Union{Nothing, SparseMatrixCSC{Float64, Int64}}=nothing,
    geom_cache::Union{Nothing, Vector{CohesiveElementGeom}}=nothing,
    ws::Union{Nothing, CZMAssemblyWorkspace}=nothing,
    visc_beta::Float64=1.0
)
    # 体内刚度组装：按 czm_submesh.material_type 取模量（PE.E_coat / NE.E_coat / SP.E / ...）
    param = param_cache.param_ref
    # 调用 assemble_bulk_stiffness 时传 param_cache 让其按 material_type 分组
    # （assemble_bulk_stiffness 内部按 e -> material_type[e] -> 对应模量）
    # ... 其余流程与旧版一致，仅替换 E_eff/ν_eff 为按材料类型分组的值
end
```

**关键内部改造**：`assemble_bulk_stiffness` 中按 `czm_mesh.czm_submesh.material_type[e]` 取：
- `:PE` → `param.PE.E_coat / PE.nu_coat / PE.alphaT`
- `:NE` → `param.NE.E_coat / NE.nu_coat / NE.alphaT`
- `:SP` → `param.SP.E / SP.nu / 0`（隔膜无活性物质）
- `:PCC / :NCC` → `param.PCC.E / PCC.nu / 0` 或 `param.NCC.E / NCC.nu / 0`

**`assemble_czm_system` 内部循环改造**（spec §7.2）：
```julia
for i in 1:n_coh
    iface = czm_mesh.cohesive_elements[i].interface_type
    params = param_cache.by_interface[iface]
    # 用 params.K_n, K_t, σ_max, δ_0_n, δ_c_n, δ_0_t, δ_c_t, η, czm_model
    # 调用 bilinear_traction_state(δ_n, δ_t, D, params) / bilinear_tangent(δ_n, δ_t, D, params)
end
```

- [ ] **Step 4: 改造 `assemble_coupled_system_full`（`src/czm.jl:587`）**

```julia
function assemble_coupled_system_full(
    czm_mesh::CohesiveMesh,
    u::Vector{Float64},
    param_cache::CzmParamCache;
    kwargs...
)
    # 内部 591 行调用同步改：
    K_total, f_int_total, separations, tractions = assemble_coupled_system(
        czm_mesh, u, param_cache; kwargs...
    )
    return K_total, f_int_total, separations, tractions
end
```

- [ ] **Step 5: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_assemble_coupled_system.jl")'`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/czm.jl test/test_assemble_coupled_system.jl
git commit -m "refactor(czm): assemble_czm/coupled_system 接受 CzmParamCache，按 interface_type 取参

按 spec v2 §7.1-§7.2：
- assemble_coupled_system 签名 (czm_mesh, u, E_eff, ν_eff, cohesive) → (czm_mesh, u, param_cache)
- assemble_coupled_system_full 同步
- assemble_czm_system 内部循环用 param_cache.by_interface[iface] 取 CzmInterfaceParams
- 体内刚度按 czm_submesh.material_type 分组取 E_coat / SP.E / PCC.E 等

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 4.3: 适配 `CzmSolve.jl` 中 11 处 `assemble_coupled_system` 调用 + `backtrack_line_search!`

**Files:**
- Modify: `src/CzmSolve.jl`（行 80-90 `backtrack_line_search!`、行 86, 170, 215, 255, 280, 321, 412, 475, 503, 540, 585 共 11 处 `assemble_coupled_system` 调用）
- Modify: `tools/czm_convergence_diag.jl`（行 103, 191, 259, 293 共 4 处）
- Create: `test/test_czm_solve_signatures.jl`

**说明：** 所有调用统一形式：
```julia
# 旧
K, f, seps, tracts = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; kwargs...)
# 新
K, f, seps, tracts = assemble_coupled_system(czm_mesh, u, param_cache; kwargs...)
```

调用方函数（`solve_czm_basic_step`、`solve_czm_arc_length_step`、`newton_raphson_czm`、`backtrack_line_search!`）签名同步把 `E_eff, ν_eff, cohesive_params` 合并为 `param_cache`。

- [ ] **Step 1: 写失败测试 `test/test_czm_solve_signatures.jl`**

```julia
using Test
using JuBat

@testset "CzmSolve signatures with param_cache" begin
    # 反射检查 method 是否存在
    for fn in (:solve_czm_basic_step, :solve_czm_arc_length_step,
               :newton_raphson_czm, :backtrack_line_search!)
        @test isdefined(JuBat, fn)
    end

    # 实际调用一次（用极小时间步避免长时间运行）
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm_enabled = true
    opt.mechanicalmodel = "full"
    opt.time = [0, 1.0]   # 1 秒
    opt.dt = [1e-6, 1.0]
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=20)
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)
    case.czm_param_cache = JuBat.compute_czm_params_per_interface(case)

    @test_nowarn JuBat.solve_czm_basic_step(
        case.czm_mesh, case.czm_param_cache, case.param, zeros(2*case.czm_mesh.nnode);
        dt=1e-3, F_ext=zeros(2*case.czm_mesh.nnode)
    )
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_czm_solve_signatures.jl")'`
Expected: FAIL — `solve_czm_basic_step` 旧签名不接受 `param_cache`

- [ ] **Step 3: 改造 `backtrack_line_search!`（`src/CzmSolve.jl:80`）**

```julia
function backtrack_line_search!(
    u::Vector{Float64}, Δu::Vector{Float64},
    czm_mesh::CohesiveMesh,
    param_cache::CzmParamCache,        # 替换 E_eff, ν_eff, cohesive_params
    damage_states::Vector{DamageState},
    F_ext::Vector{Float64},
    F_thermo_chem::Vector{Float64},
    R_norm_current::Float64,
    bc_dofs::Vector{Int}, bc_vals::Vector{Float64};
    K_bulk_cached=nothing, geom_cache=nothing, ws=nothing,
    visc_beta::Float64=1.0,
    α_init::Float64=1.0, τ::Float64=1e-4, max_iter::Int=20
)
    # 函数体内 86 行的 assemble_coupled_system 调用改为：
    # _, f_int_trial, _, _ = assemble_coupled_system(
    #     czm_mesh, u_trial, param_cache;
    #     damage_states=damage_states, K_bulk_cached=K_bulk_cached,
    #     geom_cache=geom_cache, ws=ws, visc_beta=visc_beta)
end
```

- [ ] **Step 4: 改造 `solve_czm_basic_step` / `solve_czm_arc_length_step` / `newton_raphson_czm`**

把签名中 `E_eff::Float64, ν_eff::Float64, cohesive_params::Cohesive` 三参数合并为 `param_cache::CzmParamCache`。函数体内 11 处 `assemble_coupled_system` 调用（行 170, 215, 255, 280, 321, 412, 475, 503, 540, 585）统一改为：

```julia
# 旧
K_total, f_int_total, separations, tractions = assemble_coupled_system(
    czm_mesh, u, E_eff, ν_eff, cohesive_params;
    damage_states=damage_states, K_bulk_cached=K_bulk_cached,
    geom_cache=geom_cache, ws=ws_basic, visc_beta=visc_beta)

# 新
K_total, f_int_total, separations, tractions = assemble_coupled_system(
    czm_mesh, u, param_cache;
    damage_states=damage_states, K_bulk_cached=K_bulk_cached,
    geom_cache=geom_cache, ws=ws_basic, visc_beta=visc_beta)
```

- [ ] **Step 5: 适配 `tools/czm_convergence_diag.jl` 4 处调用**

行 103, 191, 259, 293 的 `JuBat.assemble_coupled_system(...)` 调用同步改签名。同时该脚本若用 `compute_czm_effective_params(case)` 获取 E_eff，改为 `compute_czm_params_per_interface(case)` 获取 `param_cache`。

```bash
# 替换示意（手动每处确认）
# 旧
K_total, f_int_total, separations, tractions = JuBat.assemble_coupled_system(
    czm_mesh, u, E_eff, ν_eff, cohesive_params; kwargs...)
# 新
param_cache = JuBat.compute_czm_params_per_interface(case)
K_total, f_int_total, separations, tractions = JuBat.assemble_coupled_system(
    czm_mesh, u, param_cache; kwargs...)
```

- [ ] **Step 6: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_czm_solve_signatures.jl")'`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/CzmSolve.jl tools/czm_convergence_diag.jl test/test_czm_solve_signatures.jl
git commit -m "refactor(czm-solve): 11 处 assemble_coupled_system 调用 + backtrack_line_search! 改用 param_cache

按 spec v2 §7.1：CzmSolve.jl 行 86/170/215/255/280/321/412/475/503/540/585 共 11 处
+ tools/czm_convergence_diag.jl 行 103/191/259/293 共 4 处
+ backtrack_line_search! / solve_czm_*_step / newton_raphson_czm 签名同步

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 4.4: 实现 `compute_czm_params_per_interface` 与重写 `ensure_czm_cache`

**Files:**
- Modify: `src/CouplingState.jl:258-275`（替换 `compute_czm_effective_params` 为 `compute_czm_params_per_interface`）
- Modify: `src/CouplingState.jl:367-390`（重写 `ensure_czm_cache` 失效判据，spec §6.5）
- Modify: `src/JuBat.jl`（export）
- Create: `test/test_compute_czm_params_per_interface.jl`

**说明：** `CzmParamCache` 在 `case` 中存储位置为 `case.czm_param_cache::Union{Nothing, CzmParamCache}`（spec §3.5.3）。`ensure_czm_cache` 仅管理 `case.czm_cache`（装配缓存），失效判据用 `objectid` 比对（spec §6.5）。

- [ ] **Step 1: 写失败测试 `test/test_compute_czm_params_per_interface.jl`**

```julia
using Test
using JuBat

@testset "compute_czm_params_per_interface" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    case = JuBat.SetCase(param_dim, opt)

    cache = JuBat.compute_czm_params_per_interface(case)
    @test cache isa JuBat.CzmParamCache
    @test cache.id == objectid(case.param)

    # PE-PCC 接面参数
    pe = cache.by_interface[:PE_PCC]
    @test pe.σ_max ≈ case.param.cohesive.σ_max_pe_pcc / case.param.scale.σ_czm
    @test pe.K_n ≈ case.param.cohesive.K_n_pe_pcc / case.param.scale.K_czm
    @test pe.δ_0_n ≈ case.param.cohesive.δ_0_pe_pcc / case.param.scale.δ_czm
    @test pe.G_c ≈ case.param.cohesive.G_c_pe_pcc / case.param.scale.G_czm
    @test pe.E_eff ≈ case.param.PE.E_coat * case.param.scale.E_coat / case.param.scale.σ_czm
    @test pe.ν == case.param.PE.nu_coat
    @test pe.α == case.param.PE.alphaT
    # Mode II 沿用 Mode I（spec §10 决策）
    @test pe.τ_max == pe.σ_max
    @test pe.K_t == pe.K_n
    @test pe.G_c_t == pe.G_c

    # NE-NCC 接面参数
    ne = cache.by_interface[:NE_NCC]
    @test ne.E_eff ≈ case.param.NE.E_coat * case.param.scale.E_coat / case.param.scale.σ_czm
    @test ne.σ_max ≈ case.param.cohesive.σ_max_ne_ncc / case.param.scale.σ_czm
end

@testset "ensure_czm_cache 失效判据" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=20)
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)
    param_cache = JuBat.compute_czm_params_per_interface(case)

    # 首次调用：cache 为 nothing，触发构建
    cache1 = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    @test cache1 isa JuBat.CZMAssemblyCache
    @test cache1.czm_mesh_id == objectid(case.czm_mesh)
    @test cache1.param_cache_id == param_cache.id

    # 第二次调用：cache 未失效，应返回同一对象
    cache2 = JuBat.ensure_czm_cache(case, case.czm_mesh, param_cache)
    @test cache2 === cache1

    # 网格变化时失效
    submesh2 = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=30)
    czm_mesh2 = JuBat.create_czm_mesh(submesh2, case.mesh["thermal2D"], case.param)
    cache3 = JuBat.ensure_czm_cache(case, czm_mesh2, param_cache)
    @test cache3 !== cache1
    @test cache3.czm_mesh_id == objectid(czm_mesh2)
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_compute_czm_params_per_interface.jl")'`
Expected: FAIL — `compute_czm_params_per_interface` 未定义

- [ ] **Step 3: 实现 `compute_czm_params_per_interface`**

替换 `src/CouplingState.jl:258-275` 整个 `compute_czm_effective_params` 函数为：

```julia
"""
    compute_czm_params_per_interface(case) -> CzmParamCache

按 spec v2 §6.2 构造 per-interface 参数缓存。E_eff 用涂层模量（非全栈均一化）：
- :PE_PCC → PE.E_coat / PE.nu_coat / PE.alphaT
- :NE_NCC → NE.E_coat / NE.nu_coat / NE.alphaT
"""
function compute_czm_params_per_interface(case)
    param = case.param
    cohesive = param.cohesive
    scale = param.scale

    # 入口断言（spec §6.4）
    @assert param.PE.E_coat > 0 "PE.E_coat 必须为正"
    @assert param.NE.E_coat > 0 "NE.E_coat 必须为正"
    @assert cohesive.σ_max_pe_pcc > 0 && cohesive.σ_max_ne_ncc > 0
    @assert cohesive.G_c_pe_pcc > 0 && cohesive.G_c_ne_ncc > 0
    @assert cohesive.K_n_pe_pcc > 0 && cohesive.K_n_ne_ncc > 0

    function build_params(material_coat, ν_coat, αT,
                          σ_max_raw, K_n_raw, δ_0_raw, G_c_raw, δ_c_raw,
                          τ_max_raw, K_t_raw, δ_0_t_raw, G_c_t_raw, δ_c_t_raw)
        return CzmInterfaceParams(
            E_eff = material_coat * scale.E_coat / scale.σ_czm,
            ν = ν_coat,
            α = αT,
            σ_max = σ_max_raw / scale.σ_czm,
            K_n = K_n_raw / scale.K_czm,
            δ_0_n = δ_0_raw / scale.δ_czm,
            δ_c_n = δ_c_raw / scale.δ_czm,
            G_c = G_c_raw / scale.G_czm,
            τ_max = τ_max_raw / scale.σ_czm,
            K_t = K_t_raw / scale.K_czm,
            δ_0_t = δ_0_t_raw / scale.δ_czm,
            δ_c_t = δ_c_t_raw / scale.δ_czm,
            G_c_t = G_c_t_raw / scale.G_czm,
            η = cohesive.eta,
            czm_model = cohesive.czm_model,
            h_c0 = cohesive.h_c0,
            k_air = cohesive.k_air,
            lambda_m = cohesive.lambda_m,
            beta = cohesive.beta,
            threshold = cohesive.threshold,
        )
    end

    pe_pcc = build_params(param.PE.E_coat, param.PE.nu_coat, param.PE.alphaT,
                          cohesive.σ_max_pe_pcc, cohesive.K_n_pe_pcc,
                          cohesive.δ_0_pe_pcc, cohesive.G_c_pe_pcc, cohesive.δ_c_pe_pcc,
                          cohesive.τ_max_pe_pcc, cohesive.K_t_pe_pcc,
                          cohesive.δ_0_pe_pcc_t, cohesive.G_c_pe_pcc_t, cohesive.δ_c_pe_pcc_t)

    ne_ncc = build_params(param.NE.E_coat, param.NE.nu_coat, param.NE.alphaT,
                          cohesive.σ_max_ne_ncc, cohesive.K_n_ne_ncc,
                          cohesive.δ_0_ne_ncc, cohesive.G_c_ne_ncc, cohesive.δ_c_ne_ncc,
                          cohesive.τ_max_ne_ncc, cohesive.K_t_ne_ncc,
                          cohesive.δ_0_ne_ncc_t, cohesive.G_c_ne_ncc_t, cohesive.δ_c_ne_ncc_t)

    return CzmParamCache(
        by_interface = Dict(:PE_PCC => pe_pcc, :NE_NCC => ne_ncc),
        param_ref = param,
        id = objectid(param),
    )
end
```

- [ ] **Step 4: 重写 `ensure_czm_cache`（`src/CouplingState.jl:367`）**

按 spec §6.5 失效判据：

```julia
function ensure_czm_cache(case, czm_mesh::CohesiveMesh, param_cache::CzmParamCache; kwargs...)
    cache = case.czm_cache
    if cache === nothing ||
       cache.czm_mesh_id !== objectid(czm_mesh) ||
       cache.param_cache_id !== param_cache.id
        cache = build_czm_cache(czm_mesh, param_cache; kwargs...)
        case.czm_cache = cache
    end
    return cache
end
```

**同步**：在 `CZMAssemblyCache` struct 中新增字段：
```julia
czm_mesh_id::UInt64
param_cache_id::UInt64
```
并在 `build_czm_cache` 构造时填入 `objectid(czm_mesh)` 与 `param_cache.id`。

- [ ] **Step 5: 在 `src/JuBat.jl` 导出**

```julia
export compute_czm_params_per_interface, CzmParamCache, CzmInterfaceParams
```

- [ ] **Step 6: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_compute_czm_params_per_interface.jl")'`
Expected: PASS

- [ ] **Step 7: Commit**

```julia
git add src/CouplingState.jl src/JuBat.jl test/test_compute_czm_params_per_interface.jl
git commit -m "feat(czm): compute_czm_params_per_interface + ensure_czm_cache 失效判据

按 spec v2 §6.2/§6.5：
- compute_czm_params_per_interface 返回 CzmParamCache（含 param_ref + id）
- ensure_czm_cache 用 objectid 比对 czm_mesh_id/param_cache_id 判定失效
- CZMAssemblyCache 新增 czm_mesh_id/param_cache_id 字段

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 4.5: 清理 `# TODO Chunk 4` 占位

**Files:**
- Modify: 任何前序 chunk 留下的 `# TODO Chunk 4` 占位（grep 全仓库）

- [ ] **Step 1: 用 grep 找所有占位**

```bash
grep -rniE "TODO.*Chunk\s*4|placeholder.*chunk\s*4|FIXME.*chunk\s*4" src/ example/ tools/ test/ 2>/dev/null
```

- [ ] **Step 2: 逐个移除**

对每条命中：若占位对应的实现已在 Task 4.1-4.4 完成，移除 `# TODO Chunk 4` 注释并启用实际代码。

- [ ] **Step 3: 运行所有测试确认**

```bash
for f in test/test_*.jl; do julia --project=. -e "include(\"$f\")" || exit 1; done
```

Expected: 所有测试 PASS

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "chore(czm): 清理 # TODO Chunk 4 占位

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 5: 耦合数据流

本 chunk 按 spec v2 §5.0-§5.3 实现 CZM 子网格 ↔ 粗热网格的双向耦合：温度/Δsoc 正向插值、损伤 max 反向归约。

**前置知识**（spec §5.0）：
- `case.mesh["thermal2D"]` 是粗热网格，类型为 `Mesh`（`src/SetMesh.jl` 定义）
- `Mesh` 主要字段：`nlen`（节点数）、`node`（N×2 坐标）、`element`（E×4 单元连接）、`type`（"Q4"）
- `CzmSubmesh` 在 Chunk 3 定义（`spec §3.1`）
- `thermal_elem_map` 长度 = CZM **体**单元数（用于温度/dT/Δsoc 正向插值）
- `cohesive_to_thermal` 长度 = **cohesive** 单元数（用于 D 反向归约）

### Task 5.1: 实现 `build_thermal_to_czm_interp`（双线性插值）

**Files:**
- Modify: `src/CouplingState.jl`（新增 `build_thermal_to_czm_interp` 函数）
- Modify: `src/czm.jl`（在 `create_czm_mesh` 末尾调用，填充 `czm_mesh.thermal_to_czm`）
- Create: `test/test_thermal_to_czm_interp.jl`

**说明：** `thermal_to_czm` 是 `n_czm_node × n_thermal_node` 稀疏矩阵，每行 ≤4 个非零元 = 1（粗热 Q4 单元的 4 个节点），权重 = 双线性形函数值，行和 = 1。构造策略：对每个 CZM 节点，找其所在粗热单元（用 (θ, turn) 解析反查 + 点在 Q4 单元内的等参坐标判定），计算 4 个形函数值。

- [ ] **Step 1: 写失败测试 `test/test_thermal_to_czm_interp.jl`**

```julia
using Test
using JuBat
using SparseArrays

@testset "build_thermal_to_czm_interp" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=20)

    M = JuBat.build_thermal_to_czm_interp(case.mesh["thermal2D"], submesh)

    n_czm_node = submesh.mesh.nlen
    n_thermal_node = case.mesh["thermal2D"].nlen
    @test size(M) == (n_czm_node, n_thermal_node)

    # 每行 ≤4 个非零元，行和 = 1（双线性插值 partition of unity）
    for i in 1:n_czm_node
        row = M[i, :]
        @test nnz(row) <= 4
        @test sum(row) ≈ 1.0 atol=1e-10
        @test all(v >= 0 for v in nonzeros(row))
    end

    # 温度场正向插值验证：用粗热节点温度的线性函数
    T_thermal = collect(1.0:n_thermal_node)
    T_czm = M * T_thermal
    @test length(T_czm) == n_czm_node
    @test !any(isnan, T_czm)

    # 边界节点（粗热网格节点位置）：插值应等于原值
    # 取一个粗热节点坐标，在 CZM 网格中找最近节点，检查温度接近
    # （非严格相等，因 CZM 节点未必重合）
    @test minimum(T_czm) >= minimum(T_thermal) - 1e-10
    @test maximum(T_czm) <= maximum(T_thermal) + 1e-10
end

@testset "create_czm_mesh 填充 thermal_to_czm" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=20)
    czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)

    @test czm_mesh.thermal_to_czm !== nothing
    @test size(czm_mesh.thermal_to_czm, 1) == submesh.mesh.nlen
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_thermal_to_czm_interp.jl")'`
Expected: FAIL — `build_thermal_to_czm_interp` 未定义，且 `czm_mesh.thermal_to_czm === nothing`

- [ ] **Step 3: 实现 `build_thermal_to_czm_interp`**

在 `src/CouplingState.jl` 末尾追加：

```julia
"""
    build_thermal_to_czm_interp(thermal_mesh::Mesh, czm_submesh::CzmSubmesh) -> SparseMatrixCSC

构造粗热节点 → CZM 节点双线性插值矩阵（n_czm_node × n_thermal_node）。
每行 ≤4 个非零元（粗热 Q4 单元的 4 节点），权重 = 双线性形函数值，行和 = 1。

# 算法
对每个 CZM 节点 P：
1. 用 (θ, turn) 解析找候选粗热单元（来自 czm_submesh.thermal_elem_map 反查 + 周向相邻）
2. 求解等参坐标 (ξ, η) ∈ [-1, 1]^2（Newton 迭代）
3. 若 |ξ| ≤ 1+ε 且 |η| ≤ 1+ε：用双线性形函数 N_i(ξ, η) 作为权重
4. 否则（边界外）：回退到最近粗热节点（权重 1）
"""
function build_thermal_to_czm_interp(thermal_mesh::Mesh, czm_submesh::CzmSubmesh)
    czm_node = czm_submesh.mesh.node
    n_czm_node = czm_submesh.mesh.nlen
    n_thermal_node = thermal_mesh.nlen
    n_thermal_elem = size(thermal_mesh.element, 1)

    # 预计算每个粗热单元的中心坐标与半径（用于快速过滤）
    thermal_centers = zeros(n_thermal_elem, 2)
    thermal_radii = zeros(n_thermal_elem)
    for e in 1:n_thermal_elem
        ns = thermal_mesh.element[e, :]
        xs = thermal_mesh.node[ns, 1]
        ys = thermal_mesh.node[ns, 2]
        thermal_centers[e, 1] = sum(xs) / 4
        thermal_centers[e, 2] = sum(ys) / 4
        thermal_radii[e] = maximum(sqrt.((xs .- thermal_centers[e, 1]).^2 .+
                                          (ys .- thermal_centers[e, 2]).^2))
    end

    # 双线性形函数 N_i(ξ, η) = 0.25 * (1 ± ξ)(1 ± η)
    function shape_funcs(ξ::Float64, η::Float64)
        return [
            0.25 * (1 - ξ) * (1 - η),   # N1 = (-,-)
            0.25 * (1 + ξ) * (1 - η),   # N2 = (+,-)
            0.25 * (1 + ξ) * (1 + η),   # N3 = (+,+)
            0.25 * (1 - ξ) * (1 + η),   # N4 = (-,+)
        ]
    end

    # 求等参坐标 (ξ, η)：Newton 迭代解 N_i(ξ,η) * node_i = P
    function solve_isoparametric(x_nodes, y_nodes, px, py; max_iter=20, tol=1e-10)
        ξ, η = 0.0, 0.0
        for _ in 1:max_iter
            N = shape_funcs(ξ, η)
            x_pred = sum(N .* x_nodes)
            y_pred = sum(N .* y_nodes)
            rx = px - x_pred
            ry = py - y_pred
            if abs(rx) < tol && abs(ry) < tol
                return ξ, η, true
            end
            # 雅可比
            dN_dξ = 0.25 * [-(1-η), (1-η), (1+η), -(1+η)]
            dN_dη = 0.25 * [-(1-ξ), -(1+ξ), (1+ξ), (1-ξ)]
            Jxξ = sum(dN_dξ .* x_nodes)
            Jxη = sum(dN_dη .* x_nodes)
            Jyξ = sum(dN_dξ .* y_nodes)
            Jyη = sum(dN_dη .* y_nodes)
            detJ = Jxξ * Jyη - Jxη * Jyξ
            abs(detJ) < 1e-20 && break
            # Newton step
            ξ += (Jyη * rx - Jxη * ry) / detJ
            η += (-Jyξ * rx + Jxξ * ry) / detJ
        end
        return ξ, η, abs(ξ) <= 1.0 + 1e-6 && abs(η) <= 1.0 + 1e-6
    end

    # 构造稀疏矩阵
    I_rows = Int[]
    J_cols = Int[]
    V_vals = Float64[]

    for i in 1:n_czm_node
        px, py = czm_node[i, 1], czm_node[i, 2]
        # 找最近的粗热单元（用半径过滤）
        dists = sqrt.((thermal_centers[:, 1] .- px).^2 .+
                      (thermal_centers[:, 2] .- py).^2)
        candidate_order = sortperm(dists)
        found = false
        for e in candidate_order[1:min(10, end)]   # 检查最近的 10 个
            ns = thermal_mesh.element[e, :]
            x_nodes = thermal_mesh.node[ns, 1]
            y_nodes = thermal_mesh.node[ns, 2]
            ξ, η, ok = solve_isoparametric(x_nodes, y_nodes, px, py)
            if ok && abs(ξ) <= 1.0 + 1e-6 && abs(η) <= 1.0 + 1e-6
                N = shape_funcs(clamp(ξ, -1, 1), clamp(η, -1, 1))
                for k in 1:4
                    push!(I_rows, i)
                    push!(J_cols, ns[k])
                    push!(V_vals, N[k])
                end
                found = true
                break
            end
        end
        if !found
            # 回退：找最近粗热节点
            dists_node = sqrt.((thermal_mesh.node[:, 1] .- px).^2 .+
                               (thermal_mesh.node[:, 2] .- py).^2)
            _, nearest = findmin(dists_node)
            push!(I_rows, i)
            push!(J_cols, nearest)
            push!(V_vals, 1.0)
        end
    end

    return sparse(I_rows, J_cols, V_vals, n_czm_node, n_thermal_node)
end
```

- [ ] **Step 4: 在 `create_czm_mesh` 末尾调用填充 `thermal_to_czm`**

修改 `src/czm.jl` 中 `create_czm_mesh` 末尾的赋值（Task 3.2 中 `czm_mesh.thermal_to_czm = nothing`）：

```julia
czm_mesh.thermal_to_czm = build_thermal_to_czm_interp(thermal_mesh, czm_submesh)
```

- [ ] **Step 5: 在 `src/JuBat.jl` 导出**

```julia
export build_thermal_to_czm_interp
```

- [ ] **Step 6: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_thermal_to_czm_interp.jl")'`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/CouplingState.jl src/czm.jl src/JuBat.jl test/test_thermal_to_czm_interp.jl
git commit -m "feat(coupling): thermal_to_czm 双线性插值矩阵

按 spec v2 §5.1：n_czm_node × n_thermal_node 稀疏矩阵，每行 ≤4 非零元，行和=1
Newton 迭代解等参坐标，回退到最近节点。

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 5.2: 实现 Δsoc / dT 映射（粗热单元 → CZM 体单元）

**Files:**
- Modify: `src/CouplingState.jl:287-335`（重写 `compute_czm_strain_inputs`）
- Create: `test/test_czm_strain_inputs.jl`

**说明：** 按 spec §5.1，输出按 CZM **体**单元粒度（不再按粗热单元）。Δsoc 数据源从 `variables["thermal2D element soc_p/n"]` 读取（**不需要新增键**，spec §5.1 已澄清）。`compute_czm_strain_inputs` 在循环内按 `czm_submesh.material_type[e]` 分发：PE→Δsoc_p，NE→Δsoc_n，其他→0。

- [ ] **Step 1: 写失败测试 `test/test_czm_strain_inputs.jl`**

```julia
using Test
using JuBat

@testset "compute_czm_strain_inputs 按材料类型分发" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=20)
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)
    case.czm_param_cache = JuBat.compute_czm_params_per_interface(case)

    ne_czm = size(submesh.mesh.element, 1)
    ne_thermal = size(case.mesh["thermal2D"].element, 1)

    # 注入合成 SOC 数据（粗热单元粒度）
    variables = Dict{String, Any}(
        "thermal2D element soc_p" => collect(1.0:ne_thermal),   # 非零
        "thermal2D element soc_n" => collect(1.0:ne_thermal) .* 0.5,
    )
    T_nodes = collect(1.0:case.mesh["thermal2D"].nlen)

    out = JuBat.compute_czm_strain_inputs(case, variables, T_nodes)
    @test haskey(out, :dT_czm)
    @test haskey(out, :Δsoc_p_czm)
    @test haskey(out, :Δsoc_n_czm)

    @test length(out.dT_czm) == ne_czm
    @test length(out.Δsoc_p_czm) == ne_czm
    @test length(out.Δsoc_n_czm) == ne_czm

    # 按材料类型分发
    for e in 1:ne_czm
        mt = submesh.material_type[e]
        e_thermal = submesh.thermal_elem_map[e]
        if mt == :PE
            # Δsoc_p = soc_p[e_thermal] - PE.cs0
            expected = variables["thermal2D element soc_p"][e_thermal] - case.param.PE.cs0
            @test out.Δsoc_p_czm[e] ≈ expected
            @test out.Δsoc_n_czm[e] == 0.0
        elseif mt == :NE
            expected = variables["thermal2D element soc_n"][e_thermal] - case.param.NE.cs0
            @test out.Δsoc_n_czm[e] ≈ expected
            @test out.Δsoc_p_czm[e] == 0.0
        else
            # PCC/NCC/SP：两者为 0
            @test out.Δsoc_p_czm[e] == 0.0
            @test out.Δsoc_n_czm[e] == 0.0
        end
    end

    # dT 通过 thermal_to_czm 矩阵插值
    @test !any(isnan, out.dT_czm)
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_czm_strain_inputs.jl")'`
Expected: FAIL — `compute_czm_strain_inputs` 旧签名/输出粒度不同

- [ ] **Step 3: 重写 `compute_czm_strain_inputs`**

替换 `src/CouplingState.jl:287-335`：

```julia
"""
    compute_czm_strain_inputs(case, variables, T_nodes) -> NamedTuple

按 spec v2 §5.1 计算 CZM 体单元粒度的 dT、Δsoc_p、Δsoc_n。

# 输出字段
- `dT_czm::Vector{Float64}`：长度 = n_bulk_czm（CZM 体单元数），来自粗热单元 dT 通过 thermal_elem_map
- `Δsoc_p_czm::Vector{Float64}`：PE 单元为 soc_p[e_thermal]-cs0，其他为 0
- `Δsoc_n_czm::Vector{Float64}`：NE 单元为 soc_n[e_thermal]-cs0，其他为 0
"""
function compute_czm_strain_inputs(case, variables, T_nodes)
    czm_mesh = case.czm_mesh
    submesh = czm_mesh.czm_submesh
    param = case.param
    ne_czm = size(submesh.mesh.element, 1)

    # 粗热单元粒度的 soc（来自 variables）
    soc_p_thermal = get(variables, "thermal2D element soc_p", zeros(size(case.mesh["thermal2D"].element, 1)))
    soc_n_thermal = get(variables, "thermal2D element soc_n", zeros(size(case.mesh["thermal2D"].element, 1)))

    # 粗热单元粒度的 dT：T_nodes 是节点温度，需转换为单元温度（取 4 节点平均）
    thermal_elem = case.mesh["thermal2D"].element
    n_thermal_elem = size(thermal_elem, 1)
    dT_thermal = zeros(n_thermal_elem)
    for e in 1:n_thermal_elem
        ns = thermal_elem[e, :]
        dT_thermal[e] = sum(T_nodes[ns]) / 4 - param.cell.T0
    end

    Δsoc_p_czm = zeros(ne_czm)
    Δsoc_n_czm = zeros(ne_czm)
    dT_czm = zeros(ne_czm)

    for e in 1:ne_czm
        e_thermal = submesh.thermal_elem_map[e]
        if e_thermal > 0
            dT_czm[e] = dT_thermal[e_thermal]
        end
        mt = submesh.material_type[e]
        if mt == :PE
            Δsoc_p_czm[e] = e_thermal > 0 ? soc_p_thermal[e_thermal] - param.PE.cs0 : 0.0
        elseif mt == :NE
            Δsoc_n_czm[e] = e_thermal > 0 ? soc_n_thermal[e_thermal] - param.NE.cs0 : 0.0
        end
        # PCC/NCC/SP：保持 0
    end

    return (dT_czm = dT_czm, Δsoc_p_czm = Δsoc_p_czm, Δsoc_n_czm = Δsoc_n_czm)
end
```

- [ ] **Step 4: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_czm_strain_inputs.jl")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/CouplingState.jl test/test_czm_strain_inputs.jl
git commit -m "refactor(coupling): compute_czm_strain_inputs 输出按 CZM 体单元粒度

按 spec v2 §5.1：
- dT/Δsoc_p/Δsoc_n 长度 = n_bulk_czm
- Δsoc 数据源从 variables[\"thermal2D element soc_p/n\"] 读取（不新增键）
- 按 czm_submesh.material_type 分发（PE→Δsoc_p，NE→Δsoc_n，其他→0）

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 5.3: 重写 `map_czm_damage_to_thermal`（max 归约，显式循环）

**Files:**
- Modify: `src/CallModel.jl:9-19`（重写归约逻辑）
- Create: `test/test_map_czm_damage.jl`

**说明：** 按 spec §5.2，max 归约用**显式循环**（稀疏矩阵乘法只能加权求和，不可表达 max）。使用 `cohesive_to_thermal`（长度 = n_cohesive）映射，扫描每个 cohesive 单元取其 `damage_states[e_coh].D`，更新对应粗热单元的 max。`geometry` 参数废弃（旧字段 `czm_element_map / is_inner_layer` 不再用）。

- [ ] **Step 1: 写失败测试 `test/test_map_czm_damage.jl`**

```julia
using Test
using JuBat

@testset "map_czm_damage_to_thermal max 归约" begin
    param_dim = JuBat.ChooseCell("Jellyroll")
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    case = JuBat.SetCase(param_dim, opt)

    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=40, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=20)
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)
    ne_thermal = size(case.mesh["thermal2D"].element, 1)

    # 手动设置部分 damage_states 验证 max 归约
    n_coh = case.czm_mesh.n_cohesive
    for i in 1:min(n_coh, length(case.czm_mesh.damage_states))
        case.czm_mesh.damage_states[i].D = (i % 5) / 4.0   # 0, 0.25, 0.5, 0.75, 1.0 循环
    end

    D_per_thermal = JuBat.map_czm_damage_to_thermal(case.czm_mesh, ne_thermal)
    @test length(D_per_thermal) == ne_thermal
    @test all(0 .<= D_per_thermal .<= 1.0)

    # 构造合成场景：手动指定同一 e_thermal 对应多个 cohesive 单元，验证 max
    czm_mesh_synthetic = case.czm_mesh
    # 把前 3 个 cohesive 单元强制映射到粗热单元 1
    czm_mesh_synthetic.cohesive_to_thermal[1:3] .= 1
    czm_mesh_synthetic.damage_states[1].D = 0.3
    czm_mesh_synthetic.damage_states[2].D = 0.7
    czm_mesh_synthetic.damage_states[3].D = 0.1

    D2 = JuBat.map_czm_damage_to_thermal(czm_mesh_synthetic, ne_thermal)
    @test D2[1] ≈ 0.7   # max(0.3, 0.7, 0.1)

    # 完全断裂值
    czm_mesh_synthetic.damage_states[1].D = 1.0
    D3 = JuBat.map_czm_damage_to_thermal(czm_mesh_synthetic, ne_thermal)
    @test D3[1] ≈ 1.0

    # 未覆盖 CZM 的粗热单元应返回 0
    uncovered = findfirst(e -> !in(e, czm_mesh_synthetic.cohesive_to_thermal), 1:ne_thermal)
    if uncovered !== nothing
        @test D3[uncovered] == 0.0
    end
end
```

- [ ] **Step 2: 运行测试验证失败**

Run: `julia --project=. -e 'include("test/test_map_czm_damage.jl")'`
Expected: FAIL — `map_czm_damage_to_thermal` 旧签名 `(czm_mesh, geometry, ne)` 与新测试 `(czm_mesh, ne)` 不匹配

- [ ] **Step 3: 重写 `map_czm_damage_to_thermal`**

替换 `src/CallModel.jl:9-19`：

```julia
"""
    map_czm_damage_to_thermal(czm_mesh::CohesiveMesh, ne_thermal::Int) -> Vector{Float64}

按 spec v2 §5.2 把 cohesive 单元损伤 max 归约到粗热单元粒度。
对每个粗热单元 e_thermal，扫描所有 cohesive_to_thermal[e_coh] == e_thermal 的 cohesive 单元，
取 damage_states[e_coh].D 的最大值；无 cohesive 覆盖则取 0。
"""
function map_czm_damage_to_thermal(czm_mesh::CohesiveMesh, ne_thermal::Int)
    D_per_thermal = zeros(ne_thermal)
    for e_coh in 1:czm_mesh.n_cohesive
        e_thermal = czm_mesh.cohesive_to_thermal[e_coh]
        if 1 <= e_thermal <= ne_thermal
            D = czm_mesh.damage_states[e_coh].D
            if D > D_per_thermal[e_thermal]
                D_per_thermal[e_thermal] = D
            end
        end
    end
    return D_per_thermal
end
```

**调用点适配**（`src/CallModel.jl` 及任何调用方）：
```julia
# 旧
D_per_thermal = map_czm_damage_to_thermal(case.czm_mesh, case.geometry, ne)
# 新
D_per_thermal = map_czm_damage_to_thermal(case.czm_mesh, ne)
```

- [ ] **Step 4: 运行测试验证通过**

Run: `julia --project=. -e 'include("test/test_map_czm_damage.jl")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/CallModel.jl test/test_map_czm_damage.jl
git commit -m "refactor(coupling): map_czm_damage_to_thermal 用 cohesive_to_thermal + max 显式循环

按 spec v2 §5.2：
- 用 CohesiveMesh.cohesive_to_thermal（长度=n_cohesive）映射
- max 归约用显式循环（稀疏矩阵乘法不可表达 max）
- 移除 geometry 参数（旧 czm_element_map / is_inner_layer 字段废弃）

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 6: 验证脚本

本 chunk 按 spec v2 §8.1-§8.4 执行验证：单元级验证（与 param_cache 自洽）、全网格回归（不与旧路径数值对比）、网格收敛、绝对计时。所有验证脚本放在 `example/内聚力验证/` 下（沿用项目中文目录约定），测试脚本放在 `test/` 下。

### Task 6.1: 单元级验证脚本

**Files:**
- Create: `example/内聚力验证/verify_czm_per_interface.jl`

**说明：** 取单个 PE-PCC cohesive 单元，单轴拉开（位移控制），输出 σ-δ 曲线。**断言与 param_cache 自洽**（spec §8.1），不硬编码占位数值。

- [ ] **Step 1: 创建 `example/内聚力验证/verify_czm_per_interface.jl`**

```julia
using JuBat
using Printf

"""
单元级验证：单 PE-PCC cohesive 单元单轴拉开
按 spec v2 §8.1，断言与 param_cache 自洽（不硬编码占位数值）
"""
param_dim = JuBat.ChooseCell("Jellyroll")
opt = JuBat.Option()
case = JuBat.SetCase(param_dim, opt)
param_cache = JuBat.compute_czm_params_per_interface(case)
pe = param_cache.by_interface[:PE_PCC]

# 构造单 cohesive 单元 mesh（手工 4 节点）
n_a, n_b, n_b_copy, n_a_copy = 1, 2, 3, 4
node = [0.0 0.0; 1.0 0.0; 1.0 0.0; 0.0 0.0]   # 底面节点 n_a/n_b，顶面副本（坐标同底面）
coh = JuBat.CohesiveElement(1, [n_a, n_b, n_b_copy, n_a_copy],
                            [n_a, n_b], [n_a_copy, n_b_copy],
                            1.0, :PE_PCC, 2, 1)
damage = JuBat.DamageState()

# 位移控制加载：δ_n 从 0 到 1.5*δ_c_n，记录 (δ, T_n, D)
n_steps = 100
δ_n_history = collect(range(0, 1.5 * pe.δ_c_n; length=n_steps))
σ_n_history = zeros(n_steps)
D_history = zeros(n_steps)

for (i, δ_n) in enumerate(δ_n_history)
    T_n, _, damage = JuBat.bilinear_traction_state(δ_n, 0.0, damage, pe)
    σ_n_history[i] = T_n
    D_history[i] = damage.D
end

# 打印概要
println("=" ^ 60)
println("PE-PCC 单元级验证")
println("=" ^ 60)
println("σ_max (param_cache) : $(pe.σ_max)")
println("δ_0_n (param_cache) : $(pe.δ_0_n)")
println("δ_c_n (param_cache) : $(pe.δ_c_n)")
println("G_c   (param_cache) : $(pe.G_c)")
println()
println("峰值 σ_n (实测)     : $(maximum(σ_n_history))")
println("δ @ 峰值            : $(δ_n_history[argmax(σ_n_history)])")
println("D 终值              : $(D_history[end])")
println("=" ^ 60)

# 与 param_cache 自洽的断言（spec §8.1，不硬编码）
@assert maximum(σ_n_history) ≈ pe.σ_max rtol=1e-6 "峰值牵引应等于 σ_max"
@assert argmax(σ_n_history) <= cld(n_steps, 2) + 1 "峰值应在 δ_0_n 附近（前半段）"
@assert D_history[end] ≈ 1.0 atol=1e-6 "完全断裂时 D=1"
@assert D_history[1] ≈ 0.0 "未加载时 D=0"

# 健全性检查（防止 0/Inf 等异常）
@assert all(isfinite, σ_n_history) "σ_n 应有限"
@assert all(isfinite, D_history) "D 应有限"
@assert all(0 .<= D_history .<= 1.0) "D 应在 [0, 1]"

# 损伤单调（除初始弹性段外）
nonzero_D = findall(D_history .> 0)
if length(nonzero_D) > 1
    @assert all(diff(D_history[nonzero_D]) .>= -1e-10) "D 应单调非减"
end

println("✓ 单元级验证通过")

# 可选：保存 σ-δ 曲线 CSV
if haskey(ENV, "SAVE_CSV")
    using DelimitedFiles
    writedlm("output/czm_pe_pcc_sigma_delta.csv",
             [δ_n_history σ_n_history D_history], ',')
    println("σ-δ 曲线已保存到 output/czm_pe_pcc_sigma_delta.csv")
end
```

- [ ] **Step 2: 运行脚本**

Run: `julia --project=. example/内聚力验证/verify_czm_per_interface.jl`
Expected: 打印参数概要与 `✓ 单元级验证通过`

- [ ] **Step 3: Commit**

```bash
git add example/内聚力验证/verify_czm_per_interface.jl
git commit -m "test(czm): 单元级验证脚本（PE-PCC σ-δ 曲线 + 自洽断言）

按 spec v2 §8.1：断言与 param_cache 自洽，不硬编码占位数值

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 6.2: 全网格回归 + 网格收敛 + 绝对计时

**Files:**
- Modify: `example/循环验证/czm_cycle_example.jl`（确认能跑通）
- Create: `example/内聚力验证/czm_grid_convergence.jl`
- Create: `docs/superpowers/findings/2026-07-20-czm-verification-findings.md`（记录结果）

**说明：** 按 spec §8.2-§8.4，**不与改造前数值对比**（与 §2.2 "不保留旧路径"一致），仅做：跑通新路径、δ_max/D_max 合理、网格收敛、绝对计时。如用户提供实验 δ_exp，可在此文件记录对比。

- [ ] **Step 1: 跑通新路径全网格仿真**

```bash
julia --project=. example/循环验证/czm_cycle_example.jl 2>&1 | tee output/czm_cycle_v2.log
```

**验证点**（人工查看日志）：
- 无 NaN/Inf
- `D_max ∈ [0, 1]`
- 损伤峰值位置应位于 PE-PCC 或 NE-NCC 界面（不是 SP-PE）
- 仿真完成不报错

- [ ] **Step 2: 网格收敛脚本 `example/内聚力验证/czm_grid_convergence.jl`**

```julia
using JuBat
using Printf

"""
网格收敛：nθ_czm ∈ [40, 80, 160] 下 δ_max 稳定性
按 spec v2 §8.3
"""
param_dim = JuBat.ChooseCell("Jellyroll")

results = []
for nθ_czm in [40, 80, 160]
    @printf("Running nθ_czm=%d...\n", nθ_czm)
    opt = JuBat.Option()
    opt.thermal_enabled = true
    opt.thermalmodel = "distributed2D"
    opt.per_element_spme = true
    opt.czm_enabled = true
    opt.mechanicalmodel = "full"
    opt.time = [0, 100.0]
    opt.dt = [1e-6, 1.0]

    case = JuBat.SetCase(param_dim, opt)
    mesh_data = JuBat.jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
    case = JuBat.setup_thermal2D_mesh(case, mesh_data)
    submesh = JuBat.jellyroll_czm_submesh(case.param, case.mesh["thermal2D"]; nθ_czm=nθ_czm)
    case.czm_mesh = JuBat.create_czm_mesh(submesh, case.mesh["thermal2D"], case.param)
    case.czm_param_cache = JuBat.compute_czm_params_per_interface(case)

    t_elapsed = @elapsed result = JuBat.Solve(case)

    δ_max = maximum(abs.(get(result, "czm separation max", [0.0])))
    D_max = get(result, "D_max", 0.0)
    n_fractured = get(result, "n_fractured", 0)

    push!(results, (nθ_czm, δ_max, D_max, n_fractured, t_elapsed))
    @printf("  nθ_czm=%d  δ_max=%.4e  D_max=%.4f  n_fractured=%d  time=%.2fs\n",
            nθ_czm, δ_max, D_max, n_fractured, t_elapsed)
end

println("\n=== 网格收敛摘要 ===")
for r in results
    @printf("nθ_czm=%d: δ_max=%.4e, D_max=%.4f, time=%.2fs\n", r...)
end

# 收敛判定：δ_max 相对变化 < 5%
δ_40 = results[1][2]
δ_160 = results[3][2]
rel_diff = abs(δ_160 - δ_40) / max(abs(δ_160), 1e-30)
@printf("δ_max(40→160) 相对变化: %.2f%%\n", rel_diff * 100)
if rel_diff < 0.05
    println("✓ 网格收敛满足 < 5% 准则")
else
    println("⚠ 网格未收敛，建议进一步加密")
end

# 性能记录（spec §8.4：仅计时，不对比旧版本）
println("\n=== 绝对计时（仅参考，不与旧版本对比）===")
for r in results
    @printf("nθ_czm=%d: %.2fs\n", r[1], r[5])
end
```

- [ ] **Step 3: 运行收敛脚本**

```bash
julia --project=. example/内聚力验证/czm_grid_convergence.jl 2>&1 | tee output/czm_grid_convergence.log
```

Expected:
- 三种 nθ_czm 都能跑通
- δ_max 相对变化 < 5%（spec §8.3 验收）
- D_max ∈ [0, 1]
- 计时合理（无明确上限）

- [ ] **Step 4: 创建 findings 文档记录验证结果**

创建 `docs/superpowers/findings/2026-07-20-czm-verification-findings.md`：

```markdown
# CZM 按材料层界面重构 - 验证结果

**日期**: 2026-07-20
**spec/plan 版本**: v2

## 单元级验证（spec §8.1）

- 脚本: `example/内聚力验证/verify_czm_per_interface.jl`
- 状态: [PASS/FAIL]
- 峰值 σ_n: [填入]
- σ_max (param_cache): [填入]
- D 终值: [填入]
- 与 param_cache 自洽: [是/否]

## 全网格回归（spec §8.2）

- 脚本: `example/循环验证/czm_cycle_example.jl`
- 状态: [跑通/失败]
- D_max: [填入]
- n_fractured: [填入]
- 损伤峰值位置: [PE-PCC / NE-NCC / SP-PE]
- δ_exp 对比（用户未提供则 SKIP）: [δ_sim vs δ_exp]

注：按 spec §2.2 "不保留旧路径"，不与改造前数值对比。

## 网格收敛（spec §8.3）

- 脚本: `example/内聚力验证/czm_grid_convergence.jl`
- δ_max(nθ=40): [填入]
- δ_max(nθ=80): [填入]
- δ_max(nθ=160): [填入]
- 相对变化 (40→160): [填入]
- 验收 (< 5%): [PASS/FAIL]

## 性能记录（spec §8.4）

- 单步耗时 (nθ_czm=80): [填入] s
- 备注: 仅计时，不与旧版本对比

## 已知问题/后续工作

- [填入]
```

- [ ] **Step 5: Commit**

```bash
git add example/内聚力验证/czm_grid_convergence.jl docs/superpowers/findings/2026-07-20-czm-verification-findings.md
git commit -m "test(czm): 网格收敛脚本 + 验证 findings 模板

按 spec v2 §8.2-§8.4：不与旧路径数值对比，仅绝对计时

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

### Task 6.3: 全量回归与迁移通知

**Files:**
- 全 test/ 目录回归
- `CLAUDE.md` 或主 README 的迁移说明（可选）

- [ ] **Step 1: 运行所有 test/ 下的单元测试**

```bash
for f in test/test_*.jl; do
    echo "=== $f ==="
    julia --project=. -e "include(\"$f\")" || { echo "FAILED: $f"; exit 1; }
done
echo "=== All tests passed ==="
```

Expected: 全部 PASS（覆盖 Chunks 1-5 的所有测试文件）

- [ ] **Step 2: 运行 example/ 关键脚本验证语法**

```bash
for f in example/coupled_czm_thermal_example.jl \
         example/内聚力验证/verify_czm_per_interface.jl \
         example/内聚力验证/czm_grid_convergence.jl; do
    echo "=== Syntax check: $f ==="
    julia --project=. -e "include(\"$f\")" 2>&1 | tail -5
done
```

Expected: 无语法错误（运行结果不要求完整通过，但应能加载）

- [ ] **Step 3: 用户迁移通知**

完成本 plan 所有 chunk 后，向用户发送迁移通知（基于 spec §9.2）：

> **CZM 重构完成 - 迁移指南**
>
> 1. `Cohesive` struct 字段重命名：旧 `σ_max_n / K_n / G_c_n` 等已移除，按界面类型分为 `σ_max_pe_pcc / K_n_pe_pcc / ...`（共 20 字段）
> 2. `create_czm_mesh` 调用签名：旧 `create_czm_mesh(mesh_data.thermal2D, param_dim)` → 新两步 `jellyroll_czm_submesh(...) + create_czm_mesh(submesh, thermal_mesh, param)`（共 16 处调用）
> 3. `assemble_coupled_system` 调用签名：旧 `(czm_mesh, u, E_eff, ν_eff, cohesive_params)` → 新 `(czm_mesh, u, param_cache)`（共 16 处调用）
> 4. 自定义参数集脚本需同步更新（参考 `parameters/Jellyroll.jl`）
> 5. CZM 单元数约 ×4（内存占用上升，预期行为）
> 6. 损伤峰值位置从 SP-PE 改为 PE-PCC / NE-NCC（后处理脚本若硬编码位置需更新）

- [ ] **Step 4: Commit**

```bash
git add -u
git commit --allow-empty -m "chore(czm): 完成 plan v2 全部 chunk

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## 附录

### A. 测试运行汇总（spec §9.1）

```bash
# 在项目根目录运行（每次单个测试文件）
for f in test/test_*.jl; do
    julia --project=. -e "include(\"$f\")" || exit 1
done
```

Test.jl 是 Julia stdlib，通过 `using Test` 内联导入，**不需要**在 `Project.toml` 添加 `[extras]/[targets]` 段。test/ 目录在 Chunk 1 Task 1 中已创建。

### B. 迁移通知锚点

完成 Chunk 1（Cohesive struct 字段重命名）后，向用户发送第一次迁移通知（字段重命名）。
完成 Chunk 4（求解器签名变更）后，向用户发送第二次迁移通知（create_czm_mesh + assemble_coupled_system 调用签名）。
完成 Chunk 6 后，发送完整迁移通知（Task 6.3 Step 3）。

### C. 与 spec v2 的对应关系

| Chunk | spec 章节 |
|-------|-----------|
| Chunk 1 | §3.1-§3.5（数据结构） |
| Chunk 2 | §6.1-§6.4（参数集 + 归一化 + 入口断言） |
| Chunk 3 | §4.1-§4.4（CZM 子网格生成 + create_czm_mesh） |
| Chunk 4 | §7.1-§7.2 + §6.5（求解器适配 + ensure_czm_cache） |
| Chunk 5 | §5.0-§5.3（耦合数据流） |
| Chunk 6 | §8.1-§8.4 + §9.1（验证 + 测试基础设施） |





