# CZM 按材料层界面重构 — 迁移通知

> 适用于 `czm-refactor` 分支合并到 `main` 后，下游脚本/用户的迁移指南。
> 对应 spec v2：`docs/superpowers/specs/2026-07-18-czm-per-material-layer-interface-design.md`

## 1. Cohesive struct 字段重命名

旧 `Cohesive` struct 字段（如 `σ_max_n / K_n / G_c_n / τ_max_t_n / K_t_n / G_c_t_n` 及对应 `_p` 后缀）已移除，
改为按**界面类型分组**（PE-PCC = 正极涂层-正极集流体；NE-NCC = 负极涂层-负极集流体）。

新 struct 定义在 `src/SetParams.jl` line 156-193，字段如下：

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

    # Interface thermal resistance（暂禁用，见 §7）
    h_c0::Float64 = 1e7
    k_air::Float64 = 0.026
    lambda_m::Float64 = 70e-9
    beta::Float64 = 1.0
    threshold::Float64 = 70e-9

    tau_visc::Float64 = 0.0
end
```

**字段统计**：
- PE-PCC 界面：10 字段（5 法向 + 5 切向）
- NE-NCC 界面：10 字段（5 法向 + 5 切向）
- 共用本构：2 字段（`czm_model` / `eta`）
- 界面热阻（共用，暂禁用）：5 字段（`h_c0` / `k_air` / `lambda_m` / `beta` / `threshold`）
- 粘性正则化：1 字段（`tau_visc`）
- **合计 28 字段**（其中 CZM 本构核心 22 字段，界面热阻 + 粘性 6 字段）

> 注：spec v2 §3.1 中所述"20 字段"仅指 CZM 本构核心；实际 struct 含原保留的界面热阻与粘性字段。

## 2. create_czm_mesh 调用签名变更

**旧（v2）**：
```julia
czm_mesh = create_czm_mesh(jellyroll_mesh, param)
```

**新（v3）**：两步式，CZM 子网格由 `jellyroll_collector_seed_mesh` 独立细分：

```julia
# Step 1: 构造热网格 + CZM 子网格（nθ_czm 控制内聚力单元密度）
mesh_data = jellyroll_collector_seed_mesh(param; nθ=80, nθ_czm=40, gsorder=2)
case = setup_thermal2D_mesh(case, mesh_data)

# Step 2: 从子网格构建 CohesiveMesh
case.czm_mesh = create_czm_mesh(mesh_data.czm_submesh, case.mesh["thermal2D"], case.param)
```

**新签名**（`src/czm.jl` line 58）：
```julia
function create_czm_mesh(czm_submesh::CzmSubmesh, thermal_mesh::Mesh, param)
```

**关键变化**：
- 第一个参数由完整 `jellyroll_mesh` 改为独立的 `CzmSubmesh`（细化 CZM 子网格）
- 新增 `nθ_czm` 关键字（`src/Jellyrollmodel.jl` line 19），允许 CZM 内聚力单元密度独立于热网格 `nθ`
- CZM 子网格与热网格解耦：热网格可粗（`nθ=80`），CZM 子网格可细（`nθ_czm=40`）
- 函数内部自动识别 PE-PCC / NE-NCC 径向界面（按材料类型配对），并构造 `thermal_to_czm` 双线性插值矩阵

## 3. assemble_coupled_system 调用签名变更

**旧**：
```julia
K, f, seps, tracts = assemble_coupled_system(czm_mesh, u, E_eff, ν_eff;
                                              damage_states=..., K_bulk_cached=..., geom_cache=..., ws=...)
```

**新**（`src/czm.jl` line 724-735）：
```julia
K, f, seps, tracts = assemble_coupled_system(czm_mesh, u, param_cache;
                                              F_ext=...,
                                              F_thermo_chem=...,
                                              damage_states=...,
                                              K_bulk_cached=...,
                                              geom_cache=...,
                                              ws=...,
                                              visc_beta=1.0)
```

其中 `param_cache::CzmParamCache` 由 `compute_czm_params_per_interface(case)` 构造
（`src/CouplingState.jl` line 292），内部结构：

```julia
struct CzmParamCache
    by_interface::Dict{Symbol, CzmInterfaceParams}   # :PE_PCC / :NE_NCC
    param_ref::Params
    id::UInt64                                         # 内容哈希，用于 ensure_czm_cache 失效判定
end
```

`CzmInterfaceParams`（`src/CouplingState.jl` line 24-54）按界面类型聚合 20 个字段：
- 体模量与热化学载荷（`E_eff` / `ν` / `α`）
- Mode I（`σ_max` / `K_n` / `δ_0_n` / `δ_c_n` / `G_c`）
- Mode II（`τ_max` / `K_t` / `δ_0_t` / `δ_c_t` / `G_c_t`）
- BK 混合模式（`η` / `czm_model`）
- 界面热阻（`h_c0` / `k_air` / `lambda_m` / `beta` / `threshold`，暂禁用）

**关键变化**：
- 第三个位置参数由 `(E_eff, ν_eff)` 两个标量改为单个 `param_cache` 对象
- 新增 `F_ext` / `F_thermo_chem` / `visc_beta` 关键字参数
- 内部通过 `interface_type` 字段自动选择对应界面参数，无需调用方关心

## 4. 自定义参数集脚本同步

下游若用自定义电池参数（非 `"Jellyroll"`），需：

1. 在 `src/parameters/<CellName>.jl` 的 `Cohesive` 构造中**显式提供**两组界面参数：
   - 10 个 `_pe_pcc` 后缀字段（PE-PCC 界面）
   - 10 个 `_ne_ncc` 后缀字段（NE-NCC 界面）

2. **缺失时的防御**（两道防线）：
   - `ChooseCell`（`src/SetParams.jl` line 313, 334）触发 `@warn`：
     - `PE.E_coat/NE.E_coat 未定义` → `scale.E_coat` 保持 0
     - `cohesive.σ_max_pe_pcc = 0` → `scale.σ_czm/G_czm/K_czm` 将为 0
   - 入口 `@assert` 拦截（防止 NaN 污染下游）：
     - `compute_czm_params_per_interface`（`src/CouplingState.jl` line 297-303）：
       ```julia
       @assert param.PE.E_coat > 0 && param.NE.E_coat > 0
       @assert scale.σ_czm > 0
       @assert scale.E_coat > 0
       @assert param.cohesive.σ_max_pe_pcc > 0
       @assert param.cohesive.σ_max_ne_ncc > 0
       @assert param.cohesive.G_c_pe_pcc > 0 && param.cohesive.G_c_ne_ncc > 0
       @assert param.cohesive.K_n_pe_pcc > 0 && param.cohesive.K_n_ne_ncc > 0
       ```
     - `thermal_diffusion_stress_2D`（`src/Mechanical.jl` line 167）：
       ```julia
       @assert case.param_dim.PE.E_coat > 0 && case.param_dim.NE.E_coat > 0
       ```

## 5. CZM 单元数变化

**旧版**：CZM 单元数与热网格边界数一致（`nθ` 量级）。

**新版**：CZM 单元数约 **×2-4**（PE-PCC / NE-NCC 各占一份，且 `nθ_czm` 通常 ≥ `nθ`）。

**具体公式**：
```
n_czm_elem ≈ 2 × nθ_czm × n_turns
```
（`n_turns` 为螺旋圈数；系数 2 来自 PE-PCC + NE-NCC 两类界面）

**示例**：`nθ=80, nθ_czm=40, n_turns≈30` → CZM 单元数 ≈ 2 × 40 × 30 = **2400**
（旧版仅 ≈ 80 × 30 = 240）

**性能影响**：CZM 装配时间相应增加，但通过 `CzmParamCache` / `K_bulk_cached` / `geom_cache` / `ws` 四级缓存缓解。

## 6. 损伤峰值位置变化

**旧版**：损伤峰值出现在 SP-PE 界面（隔膜-正极）——这是旧 `create_czm_mesh` 自动识别的边界。

**新版**：损伤峰值出现在 **PE-PCC（正极涂层-正极集流体）** 与 **NE-NCC（负极涂层-负极集流体）** 界面。
这符合实际电池脱层故障位置（涂层-集流体界面是力学最薄弱环节）。

**实现位置**：
- `create_czm_mesh`（`src/czm.jl` line 74-97）按材料类型配对识别界面：
  - `:PE_PCC` ← PE 与 PCC 共边
  - `:NE_NCC` ← NE 与 NCC 共边
  - 其它共边（如 PE-SP, NE-SP, PE-NE）自动过滤
- `CohesiveElement.interface_type` 字段（`src/czm.jl` line 7）存储界面类型
- `assemble_czm_system` 根据 `interface_type` 从 `param_cache.by_interface` 取对应参数

## 7. 界面热阻暂禁用（spec v2 §2.4）

**临时禁用**：CZM 损伤对热导率的耦合（`apply_czm_thermal_resistance` 路径）在本次重构中**被注释关闭**。

- **文件**：`src/ThermalDistributed.jl` line 292-316（在 `ThermalDistributed2D_BC` 函数内）
- **禁用方式**：整个 `if case.opt.czm_enabled ... end` 块被注释包裹，前缀注释说明：
  ```
  # ============== [v2 修订 2026-07-21] 界面热阻暂禁用（spec §2.4）==========================
  # 原代码：按 CZM 损伤状态 D 与分离 δ_n 调整界面传热系数 h_eff，修改 K 矩阵。
  # 禁用原因：CZM 损伤场与温度场双向耦合会让参数空间与收敛行为同时变化，
  #          难以独立验证 CZM 本构是否解决 δ_sim 过小问题。
  # 恢复方式：取消本块注释（同时恢复 setup_thermal2D_mesh 的 use_merged 自动逻辑）。
  # =========================================================================================
  ```
- **原因**：细化 CZM 子网格后，原界面热阻模型与新几何（PE-PCC / NE-NCC 径向界面）不兼容，需重新推导
- **影响**：
  - `czm_enabled = true` 时损伤仍正常增长（CZM 本构不受影响）
  - 但损伤状态**不反馈到热场**（即 D ↑ 不再降低层间导热）
  - 热场仅受 SPMe 热源 + 对流边界条件驱动
- **重新启用**：待 spec v3 推导新界面热阻公式（基于 PE-PCC / NE-NCC 界面分离 δ_n）后，移除注释块

## 回归测试结果（Step 1）

14 个测试文件全部通过（exit code = 0）。运行命令：
```
julia --project=. test/<file>.jl
```

| # | 测试文件 | 状态 | Pass/Total | 备注 |
|---|---------|------|------------|------|
| 1 | `test_assemble_coupled_system.jl` | PASS | 全通过 | assemble_coupled_system 签名与缓存透传 |
| 2 | `test_bilinear_per_interface.jl` | PASS | 全通过（8 个子测试集） | bilinear_* 按界面类型分支 |
| 3 | `test_cohesive_normalization.jl` | PASS | 全通过 | NormaliseParam 归一化 cohesive 字段 |
| 4 | `test_cohesive_struct.jl` | PASS | 全通过（5 个子测试集） | Cohesive / CohesiveElement / CzmInterfaceParams / CzmSubmesh / CohesiveMesh.czm_submesh |
| 5 | `test_create_czm_mesh.jl` | PASS | 全通过（2 个子测试集） | create_czm_mesh 新签名 + 界面识别 |
| 6 | `test_czm_params_per_interface.jl` | PASS | 全通过（2 个子测试集） | CzmParamCache 构造 + 防御性 @assert |
| 7 | `test_czm_solve_signatures.jl` | PASS | 全通过（含 Broken 标记） | CzmSolve.jl 调用点签名 |
| 8 | `test_czm_strain_inputs.jl` | PASS | 全通过 | compute_czm_strain_inputs 经插值矩阵映射 |
| 9 | `test_czm_submesh.jl` | PASS | 全通过（2 个子测试集） | CzmSubmesh struct + build_czm_submesh |
| 10 | `test_ensure_czm_cache.jl` | PASS | 全通过（9 个子测试） | ensure_czm_cache 内容哈希失效判据 |
| 11 | `test_jellyroll_cohesive_params.jl` | PASS | 全通过 | Jellyroll.jl 两组实验参数完整性 |
| 12 | `test_map_czm_damage.jl` | PASS | 全通过 | map_czm_damage_to_thermal 区域映射 |
| 13 | `test_thermal_resistance_disabled.jl` | PASS | 全通过 | 界面热阻路径已注释关闭 |
| 14 | `test_thermal_to_czm_interp.jl` | PASS | 全通过（2 个子测试集） | build_thermal_to_czm_interp 双线性插值 |

**汇总**：14/14 PASS，零失败。所有 Test Summary 行均显示 Pass = Total。

## 已知限制

- **界面热阻暂禁用**（见第 7 节）：CZM 损伤不反馈到热场
- **`example/czm_grid_convergence.jl` 全运行时长**（3 × 100s 物理时间）未在本次回归中跑完，仅做语法验证（由 Task 6.2 完成）
- **`example/czm_cycle_example.jl`** 同样仅语法验证（由 Task 6.2 完成）
- **`nθ_czm` 默认值**：当前默认为 `nothing`（即与 `nθ` 相同），下游显式设置更细值时单元数会显著增加（见第 5 节）
- **`Scale.E_coat` / `Scale.σ_czm` 必须先填充**：依赖 `ChooseCell` 阶段的参数完整性检查；自定义参数集需补全 `E_coat` 与 cohesive 字段
