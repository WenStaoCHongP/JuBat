# CZM 内聚力网格按材料层界面重构设计

**日期**: 2026-07-18（2026-07-20 修订 v2）
**作者**: brainstorm with user
**分支**: czm-refactor
**状态**: 修订 v2（基于 plan reviewer 反馈补漏洞）

**v2 修订要点**（2026-07-20）：
- §3.2 `CohesiveElement` 新增 `host_outer_elem / host_inner_elem` 字段
- §3.3 `CohesiveMesh` 新增 `cohesive_to_thermal` 字段，澄清 `n_layers` 语义
- §3.5 `CzmInterfaceParams` 从 7 字段扩到 20 字段（含 Mode I/II/BK/热阻），新增 `param_ref / id`；§3.5.3 明确 `case.czm_param_cache` 存储位置
- §4.1 `jellyroll_czm_submesh` 签名加 `thermal_mesh`；§4.1.1 明确 `thermal_elem_map` 用 O(1) 解析式反查
- §4.2 明确每卷绕圈 cohesive 界面数 = 2（非 4）
- §4.3 节点复制策略加"重写外层 bulk 单元"契约与 3 条正确性自检
- §4.4 签名改为 3 参；§4.4.1 列出 16 处 create_czm_mesh 调用点
- §5.0 区分 `thermal_elem_map`（体单元维度）与 `cohesive_to_thermal`（cohesive 单元维度）
- §5.1 Δsoc 数据源回到 `variables["thermal2D element soc_p/n"]`
- §5.2 max 归约用显式循环（不可用矩阵乘法）
- §6.3 加 `SetParams.jl:324` 处理；§6.5 新增 `ensure_czm_cache` 失效判据
- §7.1 `assemble_coupled_system` 调用点 11→16（加 _full + 工具）
- §7.1 Materialmatrix.jl 涉及函数清单完整化
- §8.1 单元级验证断言改为"与 param_cache 自洽"（不硬编码占位数值）
- §8.2/8.4 与"不保留旧路径"调和，移除 git-stash 回滚步骤
- §9.1 测试基础设施约定（Test.jl stdlib、test/ 目录结构）
- §9.2 迁移清单加 create_czm_mesh 与 assemble_coupled_system 签名变更
- §10 移除已查清的 create_czm_mesh / assemble_coupled_system 未决项

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

### 2.4 分阶段启用策略（v2 修订：先禁用界面热阻）

**背景**：原设计在 CZM 启用时同时启用**界面热阻模型**（`compute_gap_conductance` 按损伤 D 与分离 δ_n 调整界面传热系数 `h_eff`；`setup_thermal2D_mesh` 自动选择**未合并网格**让界面热阻作为面接触热导进入二维热传导矩阵；见 `ThermalDistributed.jl:292-310`）。这会让 CZM 损伤场与温度场双向耦合，参数空间与收敛行为同时变化，难以独立验证 CZM 本构本身是否解决了 δ_sim 过小的问题。

**修订（2026-07-21）**：本 plan 实施时**先禁用界面热阻**，走传统二维热传导（合并网格，径向连续导热路径）。CZM 与热模型解耦为**单向耦合**——热→CZM（温度梯度作为 CZM 应力载荷），CZM→热（**不反馈**损伤对热导的影响）。

**实现方式（注释而非删除）**：保留现有界面热阻代码不动，仅在以下 3 处用注释禁用，便于后续 PR 取消注释恢复：

| # | 文件:行 | 禁用方式 | 说明 |
|---|---------|----------|------|
| 1 | `src/Jellyrollmodel.jl:531-534` | 把 `use_merged = !getfield(case_new.opt, :czm_enabled)` 注释掉，强制 `use_merged = true` | 即使 `czm_enabled=true` 也走合并网格，径向导热路径不被界面热阻截断 |
| 2 | `src/ThermalDistributed.jl:292-310` | 把整个 `if case.opt.czm_enabled ... end` 块（界面热阻 K 矩阵修改）注释掉 | CZM 损伤状态不再影响二维热传导矩阵 |
| 3 | `src/ThermalDistributed.jl:519` 附近 `get_active_elements` 调用 | 检查是否仍需调用；若仅用于热阻反馈，同步注释掉 | 断裂单元剔除是热阻模型的配套逻辑，禁用热阻时一并停用 |

每处注释统一格式：
```julia
# ============== [v2 修订 2026-07-21] 界面热阻暂禁用（先验证 CZM 本构）==============
# 原代码：xxxx
# =========================================================================================
```

**验收信号（解耦验证通过后可重新启用界面热阻）**：
- CZM 本身在单向耦合下使 δ_sim 落在 δ_exp ± 20%（§8.2 主验收）
- 全网格仿真无 NaN/Inf，D_max ∈ [0, 1]
- 重新启用后温度场不会跳变（差异 < 5%）

**非目标**：
- 不删除 `compute_gap_conductance` / `compute_all_gap_conductances` / `get_active_elements` / `effective_area_factor` 函数定义（保留备查，且后续 PR 会启用）
- 不删除 `CzmInterfaceParams` 中 `h_c0 / k_air / lambda_m / beta / threshold` 5 个热阻字段（保留，与 Cohesive 同名以便未来恢复）
- 不删除 `Cohesive` struct 中热阻字段（本次重构保留迁移）

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

新增 3 个字段：
```julia
interface_type::Symbol   # :PE_PCC 或 :NE_NCC
host_outer_elem::Int     # 外层 Q4 单元 id（在 czm_submesh.mesh.element 中的行号）
host_inner_elem::Int     # 内层 Q4 单元 id
```

`layer_idx::Int64` 字段：当前代码中无任何读取方（grep 确认仅在 `czm.jl:106` 构造函数处硬编码为 1）。**直接删除该字段**——既然没有任何消费方，保留只会留一个语义不清的"占位"。

`host_outer_elem` / `host_inner_elem` 在 `create_czm_mesh` 界面识别时填入，用于：
- 构造 `cohesive_to_thermal`（§5.0）
- 按界面类型分组归约（§5.2）


### 3.3 `CohesiveMesh` 扩展（`src/SetMesh.jl:26-46`）

新增字段：
```julia
czm_submesh::Union{Nothing, CzmSubmesh}               # 关联的细化子网格
thermal_to_czm::Union{Nothing, SparseMatrixCSC{Float64, Int}}  # 粗热节点 → 细 CZM 节点插值矩阵（n_czm_node × n_thermal_node，每行 ≤4 个非零元，行和=1）
cohesive_to_thermal::Vector{Int}                      # 每个 cohesive 单元 → 所属粗热单元 id（长度 = n_cohesive；用于 D 反向归约）
```

**`n_layers` 字段语义重定义**（`SetMesh.jl:35` 当前定义但语义不清）：原意"径向分离面数"，但旧代码硬编码为 2（与 8 材料层混淆）。本次重构明确 `n_layers = 2`（恒等：PE-PCC + NE-NCC 两个分离面类型），**不**表示材料层数（8）。新增 `czm_submesh.material_type` 才是 8 材料层语义的唯一来源。后续代码读取 `n_layers` 时须按"分离面类型数"理解。

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

#### 3.5.1 `CzmInterfaceParams` 字段定义

`CzmInterfaceParams` 必须涵盖 `Materialmatrix.jl:bilinear_traction_state`（行 68）与 `bilinear_tangent`（行 182）实际读取的**全部**字段，外加 `E_eff / ν / α` 用于体刚度与热化学载荷组装。**字段名与 `Cohesive` struct 保持一致**（避免 K_0 vs K_n 命名分裂）。

```julia
struct CzmInterfaceParams
    # ---- 体模量与热化学载荷（assemble_bulk_stiffness、assemble_thermal_chemical_load 用）----
    E_eff::Float64      # 涂层模量（PE.E_coat 或 NE.E_coat），非全栈均一化
    ν::Float64          # 涂层泊松比（PE.nu_coat 或 NE.nu_coat）
    α::Float64          # 涂层热膨胀系数（PE.alphaT 或 NE.alphaT）

    # ---- Mode I（法向）---- bilinear_* 与 compute_gap_conductance 用
    σ_max::Float64      # 最大法向牵引 [归一化]
    K_n::Float64        # 法向初始刚度 [归一化]
    δ_0_n::Float64      # 法向损伤起始位移 [归一化]
    δ_c_n::Float64      # 法向临界位移 [归一化]
    G_c::Float64        # Mode I 断裂能 [归一化]

    # ---- Mode II（切向）---- bilinear_* 用
    τ_max::Float64      # 最大切向牵引 [归一化]
    K_t::Float64        # 切向初始刚度 [归一化]
    δ_0_t::Float64      # 切向损伤起始位移 [归一化]
    δ_c_t::Float64      # 切向临界位移 [归一化]
    G_c_t::Float64      # Mode II 断裂能 [归一化]

    # ---- BK 混合模式 + 本构选择 ----
    η::Float64          # BK 准则指数（来自 Cohesive.eta，全网格单一值复制）
    czm_model::String   # 本构模型标识（来自 Cohesive.czm_model，"model1" 等）

    # ---- 界面热阻（compute_gap_conductance 用）----
    h_c0::Float64       # 完全接触界面传热系数
    k_air::Float64      # 空气导热系数
    lambda_m::Float64   # 界面微观粗糙度尺度
    beta::Float64       # 粗糙度指数
    threshold::Float64  # 间隙阈值
end
```

**与 `Cohesive` struct 的对应关系**：每个 `:PE_PCC` 接面的 `CzmInterfaceParams` 字段值来自 `Cohesive.σ_max_pe_pcc / K_n_pe_pcc / ...`，体模量来自 `param.PE.E_coat / PE.nu_coat / PE.alphaT`；`:NE_NCC` 接面对应 NE 字段。

**归一化说明**：`E_eff / σ_max / τ_max` 用 `scale.σ_czm` 归一化，`K_n / K_t` 用 `scale.K_czm`，`δ_0_n / δ_c_n / δ_0_t / δ_c_t` 用 `scale.δ_czm`，`G_c / G_c_t` 用 `scale.G_czm`，`h_c0 / k_air` 等热阻参数沿用现有归一化（`scale.K_czm` 与温度尺度复合，见 §6.3）。

#### 3.5.2 `CzmParamCache` 定义

```julia
struct CzmParamCache
    by_interface::Dict{Symbol, CzmInterfaceParams}   # 键 ∈ {:PE_PCC, :NE_NCC}
    param_ref::Params                                 # 用于 assemble_bulk_stiffness 读 PE/NE/E_coat 等
    id::UInt64                                        # objectid(param) 快速比对
end
```

#### 3.5.3 在 `case` 中的存储位置

新增 `case.czm_param_cache::Union{Nothing, CzmParamCache}` 字段（在 `Case` struct 中，与 `case.czm_cache::CZMAssemblyCache` 并列）：

| 字段 | 类型 | 职责 | 生命周期 |
|------|------|------|----------|
| `case.czm_param_cache` | `Union{Nothing, CzmParamCache}` | per-interface 材料参数（只读） | `SetCase` 后由用户/求解器一次性计算，整次仿真不变 |
| `case.czm_cache` | `CZMAssemblyCache` | 装配缓存（K_bulk、geom_cache、ws） | 每个 NR 迭代按需重建，由 `ensure_czm_cache` 管理 |

**约束**：`CzmParamCache` 一旦构造即为只读；`ensure_czm_cache` 仅管理 `czm_cache`，不重建 `czm_param_cache`。

---

## 4. CZM 子网格生成与界面识别

### 4.1 新增 `jellyroll_czm_submesh`（`src/Jellyrollmodel.jl`）

```julia
function jellyroll_czm_submesh(param, thermal_mesh::Mesh; nθ_czm::Int=80, gsorder::Int=2)
    # 1. 沿螺旋线生成 8 个径向分层的 Q4 单元/卷绕圈
    # 2. 径向分层边界（从内到外）：
    #    r_PE_inner, r_PE_outer = r_PE_inner + t_PE,
    #    r_PCC_inner, r_PCC_outer = r_PCC_inner + t_PCC,
    #    ...（按层序 PE → PCC → PE → SP → NE → NCC → NE → SP 累计）
    # 3. 每个单元的材料类型按生成顺序直接打标（不依赖坐标反推）
    # 4. thermal_elem_map：网格构造时一次性建立 CZM体单元 → 粗热单元 映射（O(1) 解析式，见下）
end
```

**关键设计点**：
- **入参**：`param`（归一化后的 `case.param`，与 `jellyroll_collector_seed_mesh` 同约定）+ `thermal_mesh`（粗热网格 `Mesh` 对象，用于解析反查）
- **不依赖坐标重合检测**（避免 `create_czm_mesh:75-80` 那种 tol=1e-8 脆弱判断）
- 材料类型直接从构造顺序确定

#### 4.1.1 `thermal_elem_map` 反查算法（O(1) 真正解析式 v3）

由于 CZM 子网格与粗热网格**共享螺旋几何参数**（`param.cell.Rin`、`param.cell.layer`），反查是真正 O(1) 解析式，**禁止使用 `findmin` 空间搜索**（O(n_thermal)，v2 审查发现旧实现仍含 `findmin`）：

```
前置：螺旋方程 r = a + b·θ_spiral + s_offset
       其中 a = Rin, b = cell.layer / (2π), s_offset 由 layer_idx 决定（构造时已知）

对每个 CZM 体单元 e_czm（已知 layer_idx 与 s_offset）：
    1. r_center = 0.5 * (r[node1] + r[node3])   # 用对角节点平均
    2. θ_spiral = (r_center - a - s_offset) / b  # 直接解析反解（无需 atan，无需 mod 2π）
    3. seg_global = clamp(floor(Int, (θ_spiral - theta0) / dθ_thermal) + 1, 1, n_thermal)
    4. thermal_elem_map[e_czm] = seg_global
```

**关键改进（v3 相对 v2）**：
- 不再用 `atan(y, x)` + `findmin` 在同 turn 候选中搜索——直接由 r 与 s_offset 反解全局 θ_spiral
- 不再需要从 `thermal_mesh.element/node` 推断 `n_seg_thermal`（不可靠）——改为强制要求 `nθ_thermal` kwarg
- 不再有"-1 映射失败"分支（clamp 保证总在 [1, n_thermal] 内）

**前置条件**：
- 调用方必须传入 `nθ_thermal`（粗热网格的周向分段数，即 `jellyroll_collector_seed_mesh(nθ=nθ_thermal)` 的同一变量）
- `n_thermal = size(thermal_mesh.element, 1)` 必须 `divrem(n_thermal, nθ_thermal)[2] == 0`（每圈单元数一致）
- `dθ_thermal = (theta1 - theta0) / n_thermal`

`thermal_mesh` 在函数内被**只读**使用一次以建立映射，之后运行时查询是 O(1) Vector 访问。

### 4.2 界面识别与 cohesive 单元数量预期

在 CZM 子网格内部，遍历径向相邻 Q4 单元对 `(e_inner, e_outer)`：

| 材料组合 | 界面类型 | 是否插 cohesive |
|----------|----------|------------------|
| `(PE, PCC)` 或 `(PCC, PE)` | `:PE_PCC` | ✓ |
| `(NE, NCC)` 或 `(NCC, NE)` | `:NE_NCC` | ✓ |
| 其他（PE-SP、SP-NE、PE-PE 等） | — | ✗ |

**每卷绕圈合格径向界面数 = 2**（一个 PE-PCC + 一个 NE-NCC）。整个网格预期 cohesive 单元数：

```
n_cohesive_expected = 2 × nθ_czm × n_turns_active
```

其中 `n_turns_active` 是 CZM 子网格覆盖的有效卷绕圈数（不含边界裁剪）。**测试断言必须用 `== n_cohesive_expected ± nθ_czm`**（±nθ_czm 容差来自最内/最外圈边界），不可用 `> 50% × ...` 这类过宽阈值。注意：旧 plan 中"4 × n_segments × n_turns"是**错误的**——一个卷绕圈的层序 `PE → PCC → PE → SP → NE → NCC → NE → SP` 内虽然有 7 个径向相邻对，但只有 2 对材料组合符合上表。

### 4.3 Cohesive 单元几何与节点复制策略

每个 cohesive 单元的 4 节点来自相邻两个 Q4 单元的共边：
- 底面 = 内层单元的外边（2 节点）
- 顶面 = 外层单元的内边（2 节点，为副本）

**节点复制策略**（新逻辑，正向描述）：
1. CZM 子网格的 Q4 单元生成时，**径向相邻的两个 Q4 单元共享一条边**（2 个共节点）
2. 在共边位置插入 cohesive 单元时，**复制该边的 2 个节点**（产生 2 个新节点，与原节点共享坐标）
3. **必须重写外层 Q4 单元的连接表**：把外层单元共边上的 2 个原节点替换为对应的副本节点。否则 cohesive 顶面节点与外层 bulk 单元共边节点不共享 DOF，分离位移恒为 0
4. 内层 Q4 单元的连接表**不变**（继续指向原节点）
5. cohesive 单元的 4 节点 = `[n_a, n_b, n_b_copy, n_a_copy]`（按逆时针，`n_a/n_b` 按 θ 递增排序）
6. 节点副本需**记忆化**：同一原节点若被多个 cohesive 共边引用（如角点处），只复制一次；用 `Dict{Int, Int}`（原节点 → 副本节点）实现

**正确性自检**（实现时必须在 `create_czm_mesh` 入口或测试中断言）：
- 对每个 cohesive 单元：`czm_mesh.node[n_a, :] ≈ czm_mesh.node[n_a_copy, :]` 且 `czm_mesh.node[n_b, :] ≈ czm_mesh.node[n_b_copy, :]`
- `length(unique(czm_mesh.cohesive_elements[i].nodes)) == 4`（4 节点不重复）
- 外层 bulk 单元的 4 节点中，共边位置必须是副本节点（非原节点）

**与现有 `CohesiveElement` 几何约定的一致性**：底/顶面节点顺序按 `czm.jl:1-8` docstring "顶面节点顺序与底面一致"。

### 4.4 重构 `create_czm_mesh`

- **旧逻辑**（`src/czm.jl:56-136`，坐标重合检测 + 螺旋界面）**移除**
- 新签名：`create_czm_mesh(czm_submesh::CzmSubmesh, thermal_mesh::Mesh, param) -> CohesiveMesh`
  - `czm_submesh`：细化子网格（含 `thermal_elem_map`，已预建）
  - `thermal_mesh`：粗热网格（用于建立 `thermal_to_czm` 与 `cohesive_to_thermal`）
  - `param`：归一化参数（用于 `Cohesive` 字段读取与 `scale` 信息）
- 实现职责：按 §4.2 识别 PE-PCC / NE-NCC 界面；按 §4.3 复制节点、重写外层 bulk 连接；构造 `CohesiveMesh` 含 `thermal_to_czm`、`cohesive_to_thermal` 字段
- **正确性自检**：构造完成后必须运行 §4.3 中列出的三条自检（throw AssertionError on failure）

#### 4.4.1 调用点清单（已核实，2026-07-20 grep）

src/ 内部仅函数定义本身，无调用。**外部调用共 16 处**（src + example + tools）：

| 文件 | 行号 | 当前调用形式 |
|------|------|--------------|
| `tools/verify_czm_system.jl` | 44 | `JuBat.create_czm_mesh(thermal_mesh, param_dim)` |
| `tools/verify_czm_standalone.jl` | 66, 129 | `JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)` |
| `tools/czm_convergence_diag.jl` | 26 | 同上 |
| `tools/czm_baseline_probe.jl` | 40, 71 | 同上 |
| `tools/check_czm_methods_coupled.jl` | 45, 123 | `JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)`（赋给 `case.czm_mesh`）|
| `tools/check_collector_mesh.jl` | 65 | `JuBat.create_czm_mesh(mesh, param)` |
| `example/coupled_czm_thermal_example.jl` | 104 | `JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)` |
| `example/testexample.jl` | 85 | 赋给 `case.czm_mesh` |
| `example/网格敏感性/5_energy_conservation_check.jl` | 69 | 同上 |
| `example/网格敏感性_v2/5_energy_conservation.jl` | 62 | 同上 |
| `example/网格敏感性/4_czm_mesh_sensitivity.jl` | 74 | 同上 |
| `example/网格敏感性_v2/4_czm_convergence.jl` | 74 | 同上 |
| `example/循环验证/czm_from_precomputed_example.jl` | 301 | `JuBat.create_czm_mesh(mesh_data.Jellyroll_czm, param_dim; tol=1e-8)` |
| `example/循环验证/czm_cycle_example.jl` | 98 | 同上 |

**统一替换模式**：

```julia
# 旧
czm_mesh = JuBat.create_czm_mesh(mesh_data.thermal2D, param_dim)

# 新（两步）
czm_submesh = JuBat.jellyroll_czm_submesh(param, mesh_data.thermal2D; nθ_czm=opt.nθ_czm)
czm_mesh = JuBat.create_czm_mesh(czm_submesh, mesh_data.thermal2D, param)
```

每处替换均需逐文件确认变量名（`param` vs `param_dim` vs `case.param`、`mesh_data.thermal2D` vs `case.mesh["thermal2D"]`）。

---

## 5. 耦合数据流（粗热 → 细 CZM）

### 5.0 两个映射字段的维度与职责

本设计区分**两条映射**，长度与索引各不相同，**不可混用**：

| 字段 | 所属结构 | 长度 | 索引语义 | 值语义 | 用途 |
|------|----------|------|----------|--------|------|
| `CzmSubmesh.thermal_elem_map` | 子网格 | `= size(mesh.element, 1)`（CZM **体**单元数） | `e_czm_bulk` ∈ [1, n_bulk] | 粗热单元 id | 温度、dT、Δsoc 的**正向**插值（粗热→细 CZM） |
| `CohesiveMesh.cohesive_to_thermal` | 内聚力网格 | `= n_cohesive`（cohesive 单元数） | `e_coh` ∈ [1, n_cohesive] | 粗热单元 id | 损伤 D 的**反向**归约（CZM→粗热） |

**构造关系**：`cohesive_to_thermal[e_coh] = thermal_elem_map[e_czm_outer]`，其中 `e_czm_outer` 是 cohesive 单元 `e_coh` 对应的外层 Q4 单元（在 `create_czm_mesh` 中界面识别时建立，存入 `CohesiveElement.host_outer_elem::Int` 字段）。

**新增 `CohesiveElement` 字段**（在 §3.2 已有 `interface_type` 基础上再加）：
```julia
host_outer_elem::Int   # 外层 Q4 单元 id（在 czm_submesh.mesh.element 中的行号）
host_inner_elem::Int   # 内层 Q4 单元 id
```
这两个字段在节点复制时填入，供 `cohesive_to_thermal` 构造与按界面类型分组归约使用。

### 5.1 输入映射（粗热 → 细 CZM）

| 来源 | 目标 | 方法 |
|------|------|------|
| `T_nodes`（粗热网格节点温度，长度 `n_thermal_node`） | CZM 体网格节点温度（长度 `n_czm_node`） | 双线性插值，预生成 `thermal_to_czm::SparseMatrixCSC`（`n_czm_node × n_thermal_node`，每行 ≤4 个非零元 = 1，行和 = 1）。每行非零元 = 该 CZM 节点所在粗热单元的 4 个节点，权重 = 双线性形函数值 |
| `dT_elem`（粗热单元温升，长度 `n_thermal_elem`） | CZM **体**单元 dT | `dT_czm_bulk[e] = dT_elem[thermal_elem_map[e]]`（1-to-1 取值） |
| `Δsoc_p_elem`, `Δsoc_n_elem`（粗热单元粒度） | CZM **体**单元 Δsoc | 按材料类型分发，见下方契约 |

**Δsoc 映射契约**（按 CZM 子网格材料类型，对 CZM **体**单元）：
- `material_type[e] == :PE`：`Δsoc_p[e] = Δsoc_p_elem[thermal_elem_map[e]]`，`Δsoc_n[e] = 0`
- `material_type[e] == :NE`：`Δsoc_n[e] = Δsoc_n_elem[thermal_elem_map[e]]`，`Δsoc_p[e] = 0`
- `material_type[e] ∈ {:PCC, :NCC, :SP}`：`Δsoc_p[e] = Δsoc_n[e] = 0`（这些材料不产生扩散应变，仅参与热应变）

**Δsoc 数据来源（在 `variables` 字典中）**：
- 现有键 `"thermal2D element soc_p"` / `"thermal2D element soc_n"`（粗热单元粒度的当前 soc）
- `compute_czm_strain_inputs` 计算 `Δsoc_X = soc_X - param.XE.cs0`（绝对差，与现有实现一致）
- **不需要**新增 `"Δsoc_p_thermal"` / `"Δsoc_n_thermal"` 键（旧 plan 误写）

### 5.2 输出反馈（CZM → 热/电）

**D 反向归约规则**：
对每个粗热单元 `e_thermal`，扫描所有 `cohesive_to_thermal[e_coh] == e_thermal` 的 **cohesive** 单元（注意：是 cohesive 单元，不是 CZM 体单元），取其 `damage_states[e_coh].D` 的**最大值**，赋给 `D_max_per_thermal_elem[e_thermal]`。若该粗热单元无对应 cohesive 单元，取 `D = 0`。

**实现形式**：max 归约不能用稀疏矩阵乘法表达（矩阵乘只能加权求和），用显式循环 + Dict 实现：

```julia
D_max_per_thermal_elem = zeros(n_thermal_elem)
for e_coh in 1:czm_mesh.n_cohesive
    e_thermal = czm_mesh.cohesive_to_thermal[e_coh]
    D = czm_mesh.damage_states[e_coh].D
    if D > D_max_per_thermal_elem[e_thermal]
        D_max_per_thermal_elem[e_thermal] = D
    end
end
```

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

`NormaliseParam` 中需为 20 个新 cohesive 字段（10/界面 × 2 界面）分别除以对应 scale。现有 scale 字段（`src/SetParams.jl:212, 220, 222, 327`）已定义：
- `scale.σ_czm`（牵引单位）—— 除 σ_max、τ_max 类
- `scale.G_czm`（断裂能单位）—— 除 G_c、G_c_t 类
- `scale.K_czm`（刚度单位，= σ_czm / δ_czm）—— 除 K_n、K_t 类
- `scale.δ_czm`（位移单位）—— 除 δ_0、δ_c 类

**新增归一化**（接口热阻参数）：`h_c0`（传热系数，单位 W/m²/K）、`k_air`（导热，W/m/K）、`lambda_m`、`threshold`（长度）、`beta`（无量纲）沿用 `parameters/Jellyroll.jl:177-181` 现有约定，**本次重构不修改热阻参数语义**（仍为单一值，与 `compute_gap_conductance` 现状一致）。

**新增 `SetParams.jl:324` 处理**：当前代码 `param_dim.scale.σ_czm = param_dim.cohesive.σ_max_n` 在移除 `σ_max_n` 后失效。改为：

```julia
# 用 PE-PCC 的 σ_max 作为 σ_czm 的参考（任意一个界面均可，仅用于无量纲化基准）
param_dim.scale.σ_czm = param_dim.cohesive.σ_max_pe_pcc
```

### 6.4 入口断言与迁移

`compute_czm_params_per_interface` 入口断言：
- `case.param.PE.E_coat > 0 && case.param.NE.E_coat > 0`
- `cohesive.σ_max_pe_pcc > 0 && cohesive.σ_max_ne_ncc > 0`
- `cohesive.G_c_pe_pcc > 0 && cohesive.G_c_ne_ncc > 0`
- `cohesive.K_n_pe_pcc > 0 && cohesive.K_n_ne_ncc > 0`

**迁移破坏性声明**：移除旧字段 `σ_max_n / K_n / G_c_n` 等后，现有 `parameters/Jellyroll.jl` 中相关行（159-171）需同步替换为按界面类型的新字段。任何引用旧字段的外部脚本（包括 `example/czm_*.jl`、`example/testexample.jl`）将**直接报错**——这是预期行为，强制用户迁移到按界面参数。

### 6.5 `ensure_czm_cache` 失效条件

`ensure_czm_cache` 仅管理 `case.czm_cache`（装配缓存），不重建 `case.czm_param_cache`。失效判据：

```julia
function ensure_czm_cache(case, czm_mesh, param_cache::CzmParamCache; ...)
    cache = case.czm_cache
    # 失效条件（任一触发即重建）：
    # 1. cache === nothing（首次调用）
    # 2. cache.czm_mesh_id !== objectid(czm_mesh)（网格变化）
    # 3. cache.param_cache_id !== param_cache.id（参数变化）
    if cache === nothing ||
       cache.czm_mesh_id !== objectid(czm_mesh) ||
       cache.param_cache_id !== param_cache.id
        cache = build_czm_cache(czm_mesh, param_cache; ...)
        case.czm_cache = cache
    end
    return cache
end
```

**约束**：`CZMAssemblyCache` 新增两个字段 `czm_mesh_id::UInt64` 与 `param_cache_id::UInt64`（在 `build_czm_cache` 构造时填入），用作快速比对。`czm_param_cache` 在 `SetCase` 后由 `compute_czm_params_per_interface(case)` 一次性构造并存入 `case.czm_param_cache`，整个仿真过程中**不变**。

---

## 7. 求解器适配

### 7.1 受影响函数清单

| 函数 | 文件 | 改造内容 |
|------|------|----------|
| `compute_czm_effective_params` | `CouplingState.jl:258` | **重命名**为 `compute_czm_params_per_interface`，返回 `CzmParamCache` |
| `ensure_czm_cache` | `CouplingState.jl:367` | 接受 `CzmParamCache`，失效判据见 §6.5 |
| `assemble_czm_system` | `czm.jl:156` | 接受 `param_cache::CzmParamCache`，循环内按 `interface_type` 取 `CzmInterfaceParams` |
| `assemble_coupled_system` | `czm.jl:564` | 签名从 `(czm_mesh, u, E_eff, ν_eff, cohesive_params)` 改为 `(czm_mesh, u, param_cache)` |
| `assemble_coupled_system_full` | `czm.jl:587` | 同步改造签名；内部 591 行调用 `assemble_coupled_system` 改为传 `param_cache` |
| `backtrack_line_search!` | `CzmSolve.jl:80` | 传递 `param_cache` 而非单一 E_eff |
| `solve_czm_basic_step` / `solve_czm_arc_length_step` / `newton_raphson_czm` | `CzmSolve.jl` | 全部 11 处 `assemble_coupled_system` 调用需适配（见 §7.1.1） |
| `compute_czm_strain_inputs` | `CouplingState.jl:287` | 输出按 CZM 子网格粒度（不再按粗热单元）；Δsoc 数据来源见 §5.1 |
| `bilinear_traction_state` / `bilinear_tangent` / `bilinear_traction` / `update_damage` | `Materialmatrix.jl:68, 163, 182, 293` | 签名从 `(δ, damage, cohesive_params::Cohesive)` 改为 `(δ, damage, params::CzmInterfaceParams)`，函数体读取 `params.K_n / K_t / δ_0_n / δ_c_n / δ_0_t / δ_c_t / η / czm_model` |
| `compute_gap_conductance` / `compute_element_gap_conductance` / `compute_all_gap_conductances` | `Materialmatrix.jl:324, 352, 398` | 接受 `params::CzmInterfaceParams`（用 `params.δ_0_n / δ_c_n / h_c0 / k_air / lambda_m / beta / threshold`），不再读 `cohesive::Cohesive`。**v2 §2.4**：签名重构保留，但调用点（`ThermalDistributed.jl:292-310`）已注释，本版本不实际调用 |
| `map_czm_damage_to_thermal` | `CallModel.jl:9-19` | 重写归约规则（§5.2）：用 `cohesive_to_thermal` 显式循环 + max |

#### 7.1.1 `assemble_coupled_system` 完整调用点清单（2026-07-20 grep）

**src/CzmSolve.jl 内部 11 处**：行 86, 170, 215, **255**, 280, 321, 412, **475**, 503, **540**, 585（粗体为旧 plan/spec 漏列的）

**src/czm.jl 内部 1 处**：行 591（`assemble_coupled_system_full` 内部调用）

**工具脚本 4 处**：`tools/czm_convergence_diag.jl` 行 103, 191, 259, 293

合计 **16 处**调用需同步改签名。所有调用统一形式：
```julia
# 旧
K, f, seps, tracts = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff, cohesive_params; kwargs...)
# 新
K, f, seps, tracts = assemble_coupled_system(czm_mesh, u, param_cache; kwargs...)
```

### 7.2 `assemble_czm_system` 内部循环改造（伪代码）

```julia
for i in 1:n_coh
    iface = czm_mesh.cohesive_elements[i].interface_type
    params = param_cache.by_interface[iface]
    # 用 params.K_n, params.K_t, params.σ_max, params.δ_0_n, params.δ_c_n, params.δ_0_t,
    #     params.δ_c_t, params.η, params.czm_model 计算单元刚度与内力
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
- **验收**（与参数集自洽，**不与 spec 中的占位数值硬编码比对**）：
  - `maximum(σ_n_history) ≈ param_cache.by_interface[:PE_PCC].σ_max`（rtol=1e-6）
  - `δ_at_peak ≈ param_cache.by_interface[:PE_PCC].δ_0_n`
  - `D_history[end] ≈ 1.0`（完全断裂）
  - 双线性形状（线性上升 → 线性软化 → 完全失效）肉眼可辨

### 8.2 全网格回归

- `example/czm_cycle_example.jl`：跑原循环仿真
- **基线说明**：旧路径跑在 SP-PE 界面（与实验几何不对应）。根据 §2.2 "不保留旧路径"，**不再跑改造前版本对比**——δ_max / D_max / n_fractured 的数值差异无意义（apples-to-oranges）
- **有意义基线**：实验 δ_exp（PE-PCC / NE-NCC 剥离曲线）
- **验收**：
  - 新路径在 PE-PCC / NE-NCC 界面的 δ_sim 落在实验 δ_exp ± 20% 范围（**用户未提供 δ_exp 时此条 SKIP**）
  - 损伤峰值**位于 PE-PCC / NE-NCC 界面**（不在 SP-PE）
  - `D_max ∈ [0, 1]`，无 NaN/Inf

**v2 修订（2026-07-21，§2.4）**：本阶段**界面热阻暂禁用**（走合并网格 + 传统二维热传导）。全网格回归只验证 CZM 本构 + 单向热耦合（热→CZM）能否解决 δ_sim 过小问题；不验证 CZM→热反馈。若 δ_sim 仍偏离 δ_exp > 20%，应优先排查 CZM 参数（σ_max / δ_c / E_coat），而非归因于热反馈缺失——后者留作后续 PR 的独立验证项。

### 8.3 网格收敛

- 不同 `nθ_czm` ∈ [40, 80, 160] 下 δ_sim 稳定性
- **验收**：δ_sim 变化 < 5%；同时 δ_max 峰值位置（极坐标角度）应收敛到同一位置（误差 < 一个周向单元）

### 8.4 性能记录（仅计时，不对比数值）

- CZM 求解器单步耗时：`@elapsed Solve(case)`
- **不与旧路径数值对比**（与 §2.2 "不保留旧路径"一致）。如需计时对比，用户须自行在新分支上跑旧版本——本 plan/spec 不内置 git-stash 回滚步骤
- **验收**：新版本单步耗时合理（无明确上限，记录用于后续优化参考）

---

## 9. 实施步骤（粗略，详细 plan 后续生成）

1. 扩展 `Cohesive` struct（含 Mode I + Mode II 共 20 字段）与 `parameters/Jellyroll.jl`（参数集）
2. 实现 `CzmSubmesh` 与 `jellyroll_czm_submesh`
3. 重构 `create_czm_mesh`（基于 CzmSubmesh，含节点复制+外层 bulk 重写）
4. 实现 `CzmParamCache` + `compute_czm_params_per_interface`
5. 改造 `assemble_czm_system`、`assemble_coupled_system`、`assemble_coupled_system_full` 与 `CzmSolve.jl` 全部 11 处 + `czm.jl:591` 共 12 处 `assemble_coupled_system` 调用，以及 `tools/czm_convergence_diag.jl` 4 处
6. 改造 `Materialmatrix.jl`（`bilinear_*`、`update_damage`、`compute_*gap_conductance*` 全部按 `CzmInterfaceParams` 取参）
7. 实现粗热→细 CZM 插值（`thermal_to_czm` 矩阵）+ Δsoc 映射契约（按 CZM 体单元粒度）
8. 改造 `map_czm_damage_to_thermal`（按 §5.2 归约规则，用 `cohesive_to_thermal`）
9. 单元级验证脚本
10. 全网格回归（不与旧路径对比数值）

### 9.1 测试基础设施约定

- **Test.jl 引入方式**：通过 `using Test` 内联导入。Test.jl 是 Julia stdlib，**不需要**在 `Project.toml` 添加 `[extras]` 或 `[targets]` 段（避免污染包依赖）
- **测试运行方式**：`julia --project=. -e 'include("test/test_xxx.jl")'`，不使用 `Pkg.test()`（后者需要 `[extras]` 配置且会激活整个 test 环境，过重）
- **test/ 目录结构**：在 Chunk 1 Task 1（数据结构扩展）首步创建 `test/` 目录；每个 Task 一个独立测试文件 `test/test_<feature>.jl`，文件内用 `@testset` 包裹
- **测试运行汇总命令**：`for f in test/test_*.jl; do julia --project=. -e "include(\"$f\")" || exit 1; done`

### 9.2 用户可见的破坏性变更（迁移清单）

实施完成后，以下用户行为将受影响：

| 变更 | 用户需做什么 |
|------|-------------|
| `Cohesive` struct 字段重命名（旧 σ_max_n / K_n / ... 移除） | 更新自定义参数集脚本，使用新按界面类型字段 |
| `parameters/Jellyroll.jl` 中 cohesive 参数块整体替换 | 已由实施步骤 1 同步更新；自定义参数集用户需同步更新 |
| 任何 `example/*.jl` 引用旧 cohesive 字段 | 实施时同步更新所有 example 脚本 |
| `create_czm_mesh` 调用签名变更（旧 2 参 → 新 3 参） | §4.4.1 列出的 16 处调用全部需迁移 |
| `assemble_coupled_system` 调用签名变更 | §7.1.1 列出的 16 处调用全部需迁移 |
| CZM 单元数增加（约 ×4） | 内存占用上升；预期行为 |
| 损伤峰值位置改变（从 SP-PE → PE-PCC / NE-NCC） | 后处理脚本（含 `CsvExport`）若硬编码位置需更新 |
| **界面热阻暂禁用**（spec §2.4 v2 修订） | 本版本即使 `czm_enabled=true` 也走合并网格 + 传统二维热传导（无 CZM→热反馈）。若用户依赖损伤影响温度场的旧行为，需知此功能暂停；恢复方式：取消 `ThermalDistributed.jl:292-310` 与 `Jellyrollmodel.jl:531-534` 的注释 |

---

## 10. 待澄清/未决事项

- **实验参数具体数值**：用户需提供 PE-PCC / NE-NCC 的 σ_max、G_c、K_n 实验值（Mode II 若无独立测量，沿用 §6.1 的 Mode II = Mode I 简化）
- **nθ_czm 默认值**：建议 80（与现 nθ=80 默认一致），用户可在 `Option` 中调整。需在 `Option` struct 新增 `nθ_czm::Int` 字段
- **`create_czm_mesh` 调用点**：✅ 已核实，见 §4.4.1（16 处）
- **`assemble_coupled_system` 调用点**：✅ 已核实，见 §7.1.1（16 处）
- **Mode II 实验数据**：当前 `parameters/Jellyroll.jl:167-171` 的 Mode II 参数沿用 Mode I 值（`τ_max_t = σ_max_n`），新设计保持这一简化，由 `compute_czm_params_per_interface` 内部复制到 `CzmInterfaceParams.τ_max / K_t / δ_0_t / δ_c_t / G_c_t`

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
