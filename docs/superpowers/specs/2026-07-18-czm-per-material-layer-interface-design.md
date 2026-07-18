# CZM 内聚力网格按材料层界面重构设计

**日期**: 2026-07-18
**作者**: brainstorm with user
**分支**: czm-refactor
**状态**: 待审查

---

## 1. 背景与动机

### 1.1 问题现象

JuBat 项目当前 CZM 仿真得到的"分离位移过小"——模拟 δ_sim 与文献/实验剥离实验的 δ_exp 相比小 1–2 个量级。

### 1.2 根因定位（基于代码考古）

**事实 1**：`create_czm_mesh`（`src/czm.jl:56-136`）通过检测螺旋外圈/内圈节点坐标重合识别 cohesive 界面。由于螺旋层序为 `PE → PCC → PE → SP → NE → NCC → NE → SP`（首尾分别是 PE 与 SP），相邻卷绕圈相接时形成的是 **SP-PE 界面**——所有当前 cohesive 单元实际上都在 **SP-PE 界面**。

**事实 2**：用户可用的实验内聚力参数（σ_c / δ_c / G_c）**仅在电极-集流体界面**（PE-PCC、NE-NCC）通过剥离实验测得；SP-PE / SP-NE 界面**无实验数据**。

**事实 3**：所有当前 cohesive 单元的 `layer_idx` 在 `czm.jl:106` 被硬编码为 1；参数全部由 `compute_czm_effective_params`（`src/CouplingState.jl:258`）返回的**全栈均一化模量** `E_eff`（来自 `compute_effective_coating_modulus`，跨 8 层厚度加权）驱动。

**事实 4**：`param.cell.layer`（`src/parameters/Jellyroll.jl:128`）是 8 层完整卷绕圈的总厚度 `2(PE+NE+SP)+PCC+NCC`。`jellyroll_collector_seed_mesh`（`Jellyrollmodel.jl:17-75`）中每个 Q4 热单元的径向跨度 = `s_out - s_in = cell.layer` = 整个 8 层叠合。**热网格每圈只有 1 个径向单元**，PE-PCC / NE-NCC 界面位于热单元内部，无法通过相邻热单元对定位。

### 1.3 结论

模型几何与实验几何**根本不匹配**：模拟在 SP-PE 界面（无实验数据）跑 CZM，用户却期望对比 PE-PCC 界面（有实验数据）的剥离曲线。差异不仅是界面位置错位，还包括：

- **界面位置错位**（最根本）：cohesive 在 SP-PE，实验在 PE-PCC / NE-NCC
- **本构参数与界面不对应**：当前所有 cohesive 共用一组 (σ_max_n, K_n, δ_c_n, G_c_n, τ_max_t, K_t, δ_c_t, G_c_t, η)，无法反映不同界面的实验差异
- **体应力场失真**：均一化 E_eff（`compute_effective_coating_modulus` 全栈加权）影响 `assemble_coupled_system` 的体刚度（`czm.jl:564`），使 PE-PCC 界面附近的局部应力集中被抹平

**注意**：δ_c 公式 (`δ_c = 2G_c/σ_max`) 与 K_n (`K_n = σ_max/δ_0`) **不直接依赖 E_eff**——它们是直接设置的参数。所以"均一化→δ 偏小"的传递路径**不是**通过 K_0，而是通过体应力场失真。换言之，即使保持 K_n/σ_max 不变，只要把 cohesive 放到正确界面、并用正确的体模量（PE.E_coat / NE.E_coat），δ_sim 就可能恢复到实验量级。

### 1.4 路径选择

| 选项 | 决策 |
|------|------|
| A：保持螺旋界面位置 + 按 interface_type 分参数 | 否决：SP-PE 无实验数据 |
| B：热网格径向细化 8× + 全材料层界面 cohesive | 否决：热计算成本 ×8，且 SP 界面参数仍需猜测 |
| C'：细化热网格 + 仅 PE-PCC/NE-NCC cohesive | 否决：热计算成本仍 ×8 |
| **C-skip-thermal**：独立细化 CZM 子网格，粗热→细 CZM 插值耦合 | **采纳** |

---

## 2. 改造目标与边界

### 2.1 目标

1. 新增独立细化的 **CZM 机械子网格**（8 径向单元/卷绕圈），与现有粗热网格解耦
2. 在 CZM 子网格的 **PE-PCC** 与 **NE-NCC** 界面位置插入 cohesive 单元
3. 每个 cohesive 单元按 `interface_type` 取**实验测得**的 σ_c / δ_c / G_c 与对应涂层模量 E_coat（不再全栈均一化）
4. 模拟 δ_sim 与实验 δ_exp 可定量对比，落在 ±20% 范围

### 2.2 非目标（明确不做）

- 不改 `jellyroll_collector_seed_mesh` 的粗热网格生成逻辑（热模型粒度不变）
- 不改 SPMe 求解器（电化学仍按原热单元粒度）
- **不在 SP-PE / SP-NE 界面放 cohesive**（无实验数据，不臆造）
- 不改 Newton-Raphson 迭代算法本身（残差判据、线搜索策略不变）
- **不保留旧路径**作为 deprecated fallback（直接拆除）

**注**：NR 框架不变，但调用链上若干函数（`assemble_coupled_system`、`solve_czm_*_step`、`backtrack_line_search!`）的**签名**会变（参数传递从单一 E_eff 改为 `CzmParamCache`），这是数据流重构，不是迭代算法重构。

### 2.3 单向门声明

C-skip-thermal 完成后若发现 δ_sim 仍偏离 δ_exp，可升级到 C'（细化热网格）。本设计的数据结构（interface_type、CzmSubmesh、按界面类型分组的参数）在升级时全部可复用。

---

## 3. 数据结构变更

### 3.1 新增 `CzmSubmesh`（`src/czm.jl`）

```julia
struct CzmSubmesh
    mesh::Mesh                          # 细化 Q4 网格（8 径向 × nθ_czm 周向）
    material_type::Vector{Symbol}       # 每个单元：:PE / :PCC / :SP / :NE / :NCC
    winding_turn::Vector{Int}           # 每个单元所属卷绕圈号
    thermal_elem_map::Vector{Int}       # 每个 CZM 单元 → 对应的粗热单元 id（用于耦合）
end
```

### 3.2 `CohesiveElement` 扩展（`src/czm.jl:1-8`）

新增 `interface_type::Symbol` 字段（取值 `:PE_PCC` 或 `:NE_NCC`）。

`layer_idx::Int64` 字段：当前代码中无任何读取方（grep 确认仅在 `czm.jl:106` 构造函数处硬编码为 1）。**直接删除该字段**——既然没有任何消费方，保留只会留一个语义不清的"占位"。

### 3.3 `CohesiveMesh` 扩展（`src/SetMesh.jl:26-46`）

新增两个字段：
```julia
czm_submesh::Union{Nothing, CzmSubmesh}
thermal_to_czm::Union{Nothing, SparseMatrixCSC{Float64, Int}}  # 粗热→细 CZM 节点插值矩阵
```

### 3.4 `Cohesive` 参数 struct 扩展（`src/SetParams.jl:155-180`）

替换原单一参数组为按界面类型分组，**涵盖 Mode I、Mode II 与 BK 混合模式**（与现有 `parameters/Jellyroll.jl:159-174` 对应）：

```julia
# PE-PCC 界面（实验测得）
σ_max_pe_pcc::Float64    # Mode I 最大法向牵引 [Pa]
K_n_pe_pcc::Float64      # 法向初始刚度 [Pa/m]
δ_0_pe_pcc::Float64      # 法向损伤起始位移 [m]
G_c_pe_pcc::Float64      # Mode I 断裂能 [J/m²]
δ_c_pe_pcc::Float64      # 法向临界位移 [m]
τ_max_pe_pcc::Float64    # Mode II 最大切向牵引 [Pa]
K_t_pe_pcc::Float64      # 切向初始刚度 [Pa/m]
δ_0_pe_pcc_t::Float64    # 切向损伤起始位移 [m]
G_c_pe_pcc_t::Float64    # Mode II 断裂能 [J/m²]
δ_c_pe_pcc_t::Float64    # 切向临界位移 [m]

# NE-NCC 界面（实验测得）
# （同上结构，后缀 _ne_ncc）
σ_max_ne_ncc, K_n_ne_ncc, δ_0_ne_ncc, G_c_ne_ncc, δ_c_ne_ncc,
τ_max_ne_ncc, K_t_ne_ncc, δ_0_ne_ncc_t, G_c_ne_ncc_t, δ_c_ne_ncc_t
```

**BK 混合模式指数 η**：当前 `cohesive.eta` 全网格单一值。无证据表明 η 随界面类型显著变化（文献上 η 主要依赖材料混合模式响应而非界面方向），**保持单一** `cohesive.eta` 字段不变。

**移除**旧字段：`σ_max_n / K_n / δ_0_n / G_c_n / δ_c_n / τ_max_t / K_t / δ_0_t / G_c_t / δ_c_t`（不保留 fallback）。

### 3.5 新增 `CzmParamCache` 与 `CzmInterfaceParams`（`src/CouplingState.jl`）

```julia
struct CzmInterfaceParams
    E_eff::Float64      # 涂层模量（PE.E_coat 或 NE.E_coat），非全栈均一化
    ν::Float64
    α::Float64
    K0::Float64         # 初始刚度，由 E_eff / 涂层厚度推出
    σ_max::Float64
    δ_c::Float64
    G_c::Float64
end

struct CzmParamCache
    by_interface::Dict{Symbol, CzmInterfaceParams}
end
```

---

## 4. CZM 子网格生成与界面识别

### 4.1 新增 `jellyroll_czm_submesh`（`src/Jellyrollmodel.jl`）

```julia
function jellyroll_czm_submesh(param; nθ_czm::Int, gsorder::Int=2)
    # 1. 沿螺旋线生成 8 个径向分层的 Q4 单元/卷绕圈
    # 2. 径向分层边界（从内到外）：
    #    r_PE_inner, r_PE_outer = r_PE_inner + t_PE,
    #    r_PCC_inner, r_PCC_outer = r_PCC_inner + t_PCC,
    #    ...（按层序 PE → PCC → PE → SP → NE → NCC → NE → SP 累计）
    # 3. 每个单元的材料类型按生成顺序直接打标（不依赖坐标反推）
    # 4. thermal_elem_map：网格构造时一次性建立 CZM单元 → 粗热单元 映射
end
```

**关键设计点**：
- **不依赖坐标重合检测**（避免 `create_czm_mesh:75-80` 那种 tol=1e-8 脆弱判断）
- 材料类型直接从构造顺序确定
- `thermal_elem_map` 在网格构造时**一次性预建立**（amortized once）；由于 CZM 子网格与粗热网格共享螺旋几何参数（`param.cell.Rin`、`param.cell.layer`），映射可通过 θ 区间与径向圈号直接解析判定，无需空间搜索。运行时查询是 O(1) Dict/Vector 访问。

### 4.2 界面识别

在 CZM 子网格内部，遍历径向相邻 Q4 单元对 `(e_inner, e_outer)`：

| 材料组合 | 界面类型 | 是否插 cohesive |
|----------|----------|------------------|
| `(PE, PCC)` 或 `(PCC, PE)` | `:PE_PCC` | ✓ |
| `(NE, NCC)` 或 `(NCC, NE)` | `:NE_NCC` | ✓ |
| 其他（PE-SP、SP-NE、PE-PE 等） | — | ✗ |

### 4.3 Cohesive 单元几何与节点复制策略

每个 cohesive 单元的 4 节点来自相邻两个 Q4 单元的共边：
- 底面 = 内层单元的外边（2 节点）
- 顶面 = 外层单元的内边（2 节点）

**节点复制策略**（新逻辑，正向描述）：
- CZM 子网格的 Q4 单元生成时，**径向相邻的两个 Q4 单元共享一条边**（2 个共节点）
- 在共边位置插入 cohesive 单元时，复制该边的 2 个节点，使内层单元与外层单元各自拥有独立节点
- cohesive 单元的 4 个节点 = 内层单元的外边 2 节点（底面）+ 外层单元的内边 2 节点（顶面）
- 这样保证分离位移发生在两节点对之间，与现有 `CohesiveElement` 几何约定一致

### 4.4 重构 `create_czm_mesh`

- **旧逻辑**（`src/czm.jl:56-136`，坐标重合检测 + 螺旋界面）**移除**
- 新签名：`create_czm_mesh(czm_submesh::CzmSubmesh, param)` —— 接受细化子网格作为输入，按 §4.2 规则识别 PE-PCC / NE-NCC 界面，按 §4.3 策略复制节点并构造 cohesive 单元
- **调用点定位**：需在实施时 grep `create_czm_mesh(` 找到所有调用方并适配。已知 `create_czm_mesh` 当前在 `src/CallModel.jl` 之外的位置被调用（具体调用点在 plan 阶段核实）

---

## 5. 耦合数据流（粗热 → 细 CZM）

### 5.1 输入映射（粗热 → 细 CZM）

| 来源 | 目标 | 方法 |
|------|------|------|
| `T_nodes`（粗热网格节点温度） | CZM 子网格节点温度 | 双线性插值，预生成 `thermal_to_czm` 稀疏矩阵 |
| `dT_elem`（粗热单元温升） | CZM 单元 `dT` | 通过 `thermal_elem_map` 取对应粗单元值 |
| `Δsoc_n_elem`, `Δsoc_p_elem`（粗热单元粒度） | CZM 单元 Δsoc | 见下方契约 |

**Δsoc 映射契约**（按 CZM 子网格材料类型）：
- `material_type == :PE` 的 CZM 单元：`Δsoc_p = Δsoc_p_elem[thermal_elem_map[e]]`，`Δsoc_n = 0`
- `material_type == :NE` 的 CZM 单元：`Δsoc_n = Δsoc_n_elem[thermal_elem_map[e]]`，`Δsoc_p = 0`
- `material_type ∈ {:PCC, :NCC, :SP}` 的 CZM 单元：`Δsoc_p = Δsoc_n = 0`（这些材料不产生扩散应变，仅参与热应变）

### 5.2 输出反馈（CZM → 热/电）

**D_max_per_thermal_elem 归约规则**：
对每个粗热单元 `e_thermal`，扫描所有 `thermal_elem_map[e_czm] == e_thermal` 的 CZM 单元，取其 `damage_states[e_czm].D` 的**最大值**，赋给 `D_max_per_thermal_elem[e_thermal]`。若该粗热单元无对应 CZM 单元（例如不在 PE-PCC / NE-NCC 界面覆盖范围内），取 `D = 0`。

**输出变量粒度变更**：
- `D_max`、`D_mean`：改为按 interface_type 分组（`D_max_pe_pcc`、`D_max_ne_ncc`、`D_mean_pe_pcc`、`D_mean_ne_ncc`），同时保留全网格聚合值用于兼容
- `n_fractured`、`soh`：保持全网格聚合，不按 interface_type 拆分（聚合维度与现有 `variables` 字典约定一致，避免破坏 CSV 导出）

### 5.3 已知风险

粗热网格分辨率有限，一个粗热单元内的温度梯度被丢失。若 PE-PCC 界面附近存在显著温度梯度（例如集流体高热导率导致局部冷点），插值会平滑应力集中，可能让 δ_sim 偏小。

**对策**：先用 C-skip-thermal 验证整体假设；若 δ_sim 仍偏离 δ_exp > 20%，升级到 C''（细化热网格）。

---

## 6. 参数与归一化

### 6.1 参数集扩展（`src/parameters/Jellyroll.jl:154-180`）

在 `cohesive = Cohesive()` 块中填入两组实验值（待用户填入具体数值，单位：Pa, m, J/m²）：

```julia
# PE-PCC 界面（实验剥离测得）
cohesive.σ_max_pe_pcc = ...   # Pa
cohesive.K_n_pe_pcc   = ...   # Pa/m（保持与现有 K_n=2.4e17 量级一致或按实验调整）
cohesive.δ_0_pe_pcc   = cohesive.σ_max_pe_pcc / cohesive.K_n_pe_pcc
cohesive.G_c_pe_pcc   = ...   # J/m²
cohesive.δ_c_pe_pcc   = 2 * cohesive.G_c_pe_pcc / cohesive.σ_max_pe_pcc

# Mode II 参数（若实验未测，沿用现有 Mode II = Mode I 的简化）
cohesive.τ_max_pe_pcc     = cohesive.σ_max_pe_pcc
cohesive.K_t_pe_pcc       = cohesive.K_n_pe_pcc
cohesive.δ_0_pe_pcc_t     = cohesive.τ_max_pe_pcc / cohesive.K_t_pe_pcc
cohesive.G_c_pe_pcc_t     = cohesive.G_c_pe_pcc
cohesive.δ_c_pe_pcc_t     = 2 * cohesive.G_c_pe_pcc_t / cohesive.τ_max_pe_pcc

# NE-NCC 界面（同上结构）
# ...（待用户填入）
```

### 6.2 体模量按界面类型分别取值

替换 `compute_czm_effective_params`（`CouplingState.jl:258-275`）为：

```julia
function compute_czm_params_per_interface(case)
    # PE-PCC：用 PE 涂层模量（非全栈均一化）
    E_eff_pe_pcc = case.param.PE.E_coat * scale.E_coat / scale.σ_czm
    # NE-NCC：用 NE 涂层模量
    E_eff_ne_ncc = case.param.NE.E_coat * scale.E_coat / scale.σ_czm

    return CzmParamCache(Dict(
        :PE_PCC => CzmInterfaceParams(E_eff_pe_pcc, ν_pe, α_pe, ...),
        :NE_NCC => CzmInterfaceParams(E_eff_ne_ncc, ν_ne, α_ne, ...)
    ))
end
```

**关键点**：
- E_eff 不再跨 8 层加权（消除均一化导致体应力场失真）
- `ν_eff`、`α_eff` 也按涂层取值（PE.nu_coat / NE.nu_coat / PE.alphaT / NE.alphaT）
- **K_n 不变**（保持当前直接设置的惩罚刚度）；现有 `K_n = σ_max / δ_0` 模式保留，每个界面类型独立设置 K_n_pe_pcc / K_n_ne_ncc

### 6.3 归一化扩展（`src/SetParams.jl:469-479`）

`NormaliseParam` 中需为 20 个新字段（10/界面 × 2 界面）分别除以对应 scale。现有 scale 字段（`src/SetParams.jl:212, 220, 222, 327`）已定义：
- `scale.σ_czm`（牵引单位）—— 除 σ_max、τ_max 类
- `scale.G_czm`（断裂能单位）—— 除 G_c 类
- `scale.K_czm`（刚度单位，= σ_czm / δ_czm）—— 除 K_n、K_t 类
- `scale.δ_czm`（位移单位）—— 除 δ_0、δ_c 类

### 6.4 入口断言与迁移

`compute_czm_params_per_interface` 入口断言：
- `case.param.PE.E_coat > 0 && case.param.NE.E_coat > 0`
- `cohesive.σ_max_pe_pcc > 0 && cohesive.σ_max_ne_ncc > 0`
- `cohesive.G_c_pe_pcc > 0 && cohesive.G_c_ne_ncc > 0`
- `cohesive.K_n_pe_pcc > 0 && cohesive.K_n_ne_ncc > 0`

**迁移破坏性声明**：移除旧字段 `σ_max_n / K_n / G_c_n` 等后，现有 `parameters/Jellyroll.jl` 中相关行（159-171）需同步替换为按界面类型的新字段。任何引用旧字段的外部脚本（包括 `example/czm_*.jl`、`example/testexample.jl`）将**直接报错**——这是预期行为，强制用户迁移到按界面参数。

---

## 7. 求解器适配

### 7.1 受影响函数清单

| 函数 | 文件 | 改造内容 |
|------|------|----------|
| `compute_czm_effective_params` | `CouplingState.jl:258` | **重命名**为 `compute_czm_params_per_interface`，返回 `CzmParamCache` |
| `ensure_czm_cache` | `CouplingState.jl:367` | 接受 `CzmParamCache`，按 interface_type 分组缓存 |
| `assemble_czm_system` | `czm.jl:156` | 接受 `param_cache::CzmParamCache`，循环内按 `interface_type` 取参数 |
| `assemble_coupled_system` | `czm.jl:564` | 同上，签名从 `(czm_mesh, u, E_eff, ν_eff, cohesive_params)` 改为 `(czm_mesh, u, param_cache)` |
| `backtrack_line_search!` | `CzmSolve.jl:80` | 传递 `param_cache` 而非单一 E_eff |
| `solve_czm_basic_step`、`solve_czm_arc_step` 等 | `CzmSolve.jl:170, 215, 280, 321, 412, 503, 585` | 同上（共 7 处 `assemble_coupled_system` 调用需适配） |
| `compute_czm_strain_inputs` | `CouplingState.jl:287` | 输出按 CZM 子网格粒度（不再按粗热单元） |
| `Materialmatrix.jl` 中 `cohesive_params.K_n / K_t / δ_0_n / ...` | `Materialmatrix.jl:69-70, 183-185` 等 | 从单一结构取值改为按 `interface_type` 查 `param_cache.by_interface[iface]` |
| `map_czm_damage_to_thermal` | `CallModel.jl:9-19` | 适配按 interface_type 分组的 damage_states，重写归约规则（§5.2） |

### 7.2 `assemble_czm_system` 内部循环改造（伪代码）

```julia
for i in 1:n_coh
    iface = czm_mesh.cohesive_elements[i].interface_type
    params = param_cache.by_interface[iface]
    # 用 params.K0, params.σ_max, params.δ_c 计算单元刚度矩阵与内力
    # 其余流程不变
end
```

参数查表是 O(1) Dict 查找，性能影响可忽略。

### 7.3 损伤状态

`DamageState` 结构（`czm.jl:24-35`）不变。`damage_states` 向量维度 = cohesive 单元数（按新网格重建，比旧网格大约 ×4）。

---

## 8. 验证方案

### 8.1 单元级验证（先做）

新增 `example/内聚力验证/verify_czm_per_interface.jl`：
- 取单个 PE-PCC cohesive 单元，单轴拉开（位移控制）
- 输出 σ-δ 曲线，与实验剥离曲线对比
- **验收**：σ_max、δ_c、双线性形状三者定量匹配

### 8.2 全网格回归

- `example/czm_cycle_example.jl`：跑原循环仿真
- **基线说明**：旧路径跑在 SP-PE 界面（与实验几何不对应），改造前后 δ_max / D_max / n_fractured 的数值差异**不作为验收依据**——这是 apples-to-oranges 比较
- **有意义基线**：实验 δ_exp（PE-PCC / NE-NCC 剥离曲线）
- **验收**：新路径在 PE-PCC / NE-NCC 界面的 δ_sim 落在实验 δ_exp ± 20% 范围；同时新路径应显示损伤峰值**位于 PE-PCC / NE-NCC 界面**而非 SP-PE 界面

### 8.3 网格收敛

- 不同 `nθ_czm` ∈ [40, 80, 160] 下 δ_sim 稳定性
- **验收**：δ_sim 变化 < 5%；同时 δ_max 峰值位置（极坐标角度）应收敛到同一位置（误差 < 一个周向单元）

### 8.4 性能记录

- CZM 求解器单步耗时（改造前 vs 后）
- **验收**：单步耗时增加 < 30%（单元数约 ×4，参数查表 O(1)）

---

## 9. 实施步骤（粗略，详细 plan 后续生成）

1. 扩展 `Cohesive` struct（含 Mode I + Mode II 共 20 字段）与 `parameters/Jellyroll.jl`（参数集）
2. 实现 `CzmSubmesh` 与 `jellyroll_czm_submesh`
3. 重构 `create_czm_mesh`（基于 CzmSubmesh）
4. 实现 `CzmParamCache` + `compute_czm_params_per_interface`
5. 改造 `assemble_czm_system`、`assemble_coupled_system` 与 `CzmSolve.jl` 全部 7 处 `assemble_coupled_system` 调用
6. 改造 `Materialmatrix.jl`（按 interface_type 取 K_n / K_t / δ_0 等）
7. 实现粗热→细 CZM 插值（`thermal_to_czm` 矩阵）+ Δsoc 映射契约
8. 改造 `map_czm_damage_to_thermal`（按 §5.2 归约规则）
9. 单元级验证脚本
10. 全网格回归

### 9.1 用户可见的破坏性变更（迁移清单）

实施完成后，以下用户行为将受影响：

| 变更 | 用户需做什么 |
|------|-------------|
| `Cohesive` struct 字段重命名（旧 σ_max_n / K_n / ... 移除） | 更新自定义参数集脚本，使用新按界面类型字段 |
| `parameters/Jellyroll.jl` 中 cohesive 参数块整体替换 | 已由实施步骤 1 同步更新；自定义参数集用户需同步更新 |
| 任何 `example/*.jl` 引用旧 cohesive 字段 | 实施时同步更新所有 example 脚本 |
| CZM 单元数增加（约 ×4） | 内存占用上升；预期行为 |
| 损伤峰值位置改变（从 SP-PE → PE-PCC / NE-NCC） | 后处理脚本（含 `CsvExport`）若硬编码位置需更新 |

---

## 10. 待澄清/未决事项

- **实验参数具体数值**：用户需提供 PE-PCC / NE-NCC 的 σ_max、G_c、K_n 实验值（Mode II 若无独立测量，沿用 §6.1 的 Mode II = Mode I 简化）
- **nθ_czm 默认值**：建议 80（与现 nθ=80 默认一致），用户可在 `Option` 中调整。需在 Option struct 新增 `nθ_czm::Int` 字段
- **`create_czm_mesh` 现有调用点**：当前 spec 中调用链尚未完全核实（`CallModel.jl:12` 是 `map_czm_damage_to_thermal` 内部，非 `create_czm_mesh` 调用点）。实施时需 grep `create_czm_mesh(` 找到所有调用方并适配
- **Mode II 实验数据**：当前 `parameters/Jellyroll.jl:167-171` 的 Mode II 参数沿用 Mode I 值（`τ_max_t = σ_max_n`），新设计若保持这一简化，需在 spec 注释中明示

---

## 附录 A：相关文件清单

### A.1 必改文件
- `src/czm.jl` — CohesiveElement 扩展、create_czm_mesh 重构、CzmSubmesh 定义
- `src/SetMesh.jl` — CohesiveMesh 扩展
- `src/SetParams.jl` — Cohesive struct 扩展、NormaliseParam 新增字段归一化
- `src/parameters/Jellyroll.jl` — 填入两组实验参数
- `src/CouplingState.jl` — compute_czm_params_per_interface、ensure_czm_cache、CzmParamCache
- `src/CzmSolve.jl` — backtrack_line_search!、solve_czm_*_step 签名变更
- `src/Jellyrollmodel.jl` — 新增 jellyroll_czm_submesh
- `src/CallModel.jl` — 调用链适配

### A.2 新增文件
- `example/内聚力验证/verify_czm_per_interface.jl` — 单元级验证

### A.3 不变文件
- `src/SPMe.jl` — 电化学不变
- `src/ThermalDistributed.jl` — 热模型不变
- `src/jellyroll_collector_seed_mesh` — 粗热网格生成不变

---

## 附录 B：与 CLAUDE.md 约定的一致性

- ✅ 参数 struct 放在 `src/SetParams.jl`
- ✅ 状态/变量放在 `src/Variables.jl` 或最匹配的现有文件
- ✅ 函数合并优先于函数新增（移除 `compute_czm_effective_params`，合并为单一新函数）
- ✅ 不保留 deprecated 路径
- ✅ 不轻易新建文件（仅一个验证脚本）
