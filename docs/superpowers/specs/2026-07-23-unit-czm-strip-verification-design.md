# 单元级 CZM 条带验证设计

日期：2026-07-23  
状态：已批准（brainstorming）  
参数源：`ChooseCell("Jellyroll")` → `SetCase` 归一化

## 1. 目标

实现**单元级内聚力网格生成**（1 个重复单元的 8 层叠层 = **8×Q4 + 4×COH2D4**），并产出两个独立验证脚本：

1. **位移边界条件**：不接入电化学与热源，验证双线性牵引-分离律（Mode I 全曲线、卸载/重加载、Mode II）。
2. **合成本征应变**：不跑真实 SPMe/分布式热，用合成 `ΔT`/`Δsoc` 斜坡驱动 `assemble_thermal_chemical_load`，在弹性段与 **1D 叠层解析解**定量比对界面分离位移。

完整电化学-热-CZM 耦合仍由现有 smoke（如 `test/smoke_czm_redesign.jl`）覆盖；本规格只验证条带几何上的本构与本征应变管线。

## 2. 范围

### 2.1 包含

- `src/CzmUnitMesh.jl`：`create_unit_czm_strip`
- 复用生产 `create_czm_mesh`（方案 C：生成 + 硬断言）
- `test/unit_czm_bilinear.jl`
- `test/unit_czm_eigenstrain.jl`
- 脚本内局部增量 Newton（自定义 `bc_dofs` / `bc_vals`）

### 2.2 不包含

- 修改生产 `create_czm_mesh` / `bilinear_*` / `solve_czm_step` 逻辑
- 真实 SPMe 步进、集总/分布式热求解、双向损伤反馈
- 环几何 `identify_bc_nodes_czm`、弧长法
- 损伤段解析闭式解、CI 强制挂接（可后续）

## 3. 架构

```
src/CzmUnitMesh.jl          ← create_unit_czm_strip(...)
        │ 构造 CzmSubmesh + 1 单元哑热网格
        ▼
   create_czm_mesh(...)     ← 生产路径（节点复制 / 外层重写）
        │ + 硬断言（4 coh、界面类型、副本坐标）
        ▼
   CohesiveMesh (条带)
        │
   ┌────┴────────────────────┐
   ▼                         ▼
test/unit_czm_bilinear.jl   test/unit_czm_eigenstrain.jl
  位移 BC 增量加载              合成 ΔT/Δsoc 斜坡
  Mode I + 卸载 + Mode II       assemble_*_full + 1D 解析解
```

**求解约定**：不调用 `solve_czm_step` / `identify_bc_nodes_czm`（假定环形几何）。脚本调用 `assemble_coupled_system` 或 `assemble_coupled_system_full` + `apply_bc_czm(; bc_dofs, bc_vals)`，自行管理 `DamageState` 与增量循环。

**参数约定**：`compute_czm_params_per_interface(case)`；模量经 `moduli_of`（σ_czm 空间）；位移↔分离用界面 `Λ`。

## 4. 网格生成器 `create_unit_czm_strip`

### 4.1 API

```julia
create_unit_czm_strip(param; width=nothing, y0=1.0, gsorder=2)
    -> (czm_mesh::CohesiveMesh, meta::NamedTuple)
```

文件：`src/CzmUnitMesh.jl`，在 `src/JuBat.jl` 中于 `czm.jl` 之后 `include`。

### 4.2 几何

- **形态**：平直矩形条带；叠层方向 = **+y**（法向恒为 +y）；宽度方向 = **x**。
- **层序**（自下而上，与 `Jellyrollmodel.jl` / `build_czm_submesh` 一致）：

  | 层 | 材料 | 厚度字段 |
  |----|------|----------|
  | 1 | PE | `param.PE.thickness` |
  | 2 | PCC | `param.PCC.thickness` |
  | 3 | PE | `param.PE.thickness` |
  | 4 | SP | `param.SP.thickness` |
  | 5 | NE | `param.NE.thickness` |
  | 6 | NCC | `param.NCC.thickness` |
  | 7 | NE | `param.NE.thickness` |
  | 8 | SP | `param.SP.thickness` |

- **宽度**：默认 `W = sum(layer_thicknesses)`（≈ `param.cell.layer`），可由 `width` 覆盖。
- **放置**：底边 `y = y0`（默认 `1.0`，无量纲，保证质心 `r > 0`）；`x ∈ [0, W]`。
- **节点编号**：沿 x 递增（左→右），满足 `create_czm_mesh` 用节点 id 排序共边的约定。

### 4.3 拓扑

- 9 条水平节点线 × 2 纵向节点 → **18 节点、8 个 Q4**（每层 1 个）。
- `CzmSubmesh.material_type` = 上表层序；`winding_turn` 全为 `1`。
- `thermal_elem_map` 全指向 `1`；配 **1 单元哑热网格**（仅满足 `create_czm_mesh` 签名与断言）。

### 4.4 生产路径 + 硬断言（方案 C）

1. 组装 `CzmSubmesh` + 哑热网格。
2. 调用生产 `create_czm_mesh(submesh, dummy_thermal, param)`。
3. 断言（失败即 `error`）：
   - `n_cohesive == 4`
   - `count(==(:PE_PCC), types) == 2` 且 `count(==(:NE_NCC), types) == 2`
   - 每个 cohesive：副本节点坐标与原节点重合（`atol=1e-12`）；4 节点互异
   - 外层 bulk 单元共边节点已替换为副本 id（分离位移可非零）

### 4.5 `meta` 字段

供脚本施加 BC，至少包含：

- `y_interfaces`：层界 y 坐标（含底/顶）
- `bottom_nodes` / `top_nodes`：底边、顶边节点索引
- `cohesive_ids`：长度 4 的 cohesive 单元索引
- `interface_types`：对应 `:PE_PCC` / `:NE_NCC`
- `layer_materials`：长度 8 的材料 Symbol

## 5. 脚本 1：`test/unit_czm_bilinear.jl`

### 5.1 局部求解器

增量 Newton：每载荷步设置 `bc_dofs`/`bc_vals`（罚法支持非零位移），组装耦合系统，更新损伤。底边 `u_x = u_y = 0`；顶边按工况施加位移；侧向 DOF 按工况约束，避免刚体/剪切污染。

位移与分离换算：`u_imposed = δ̃ / Λ`（`δ̃` 在分离空间，与 `δ_0*`/`δ_c*` 同参考）。

### 5.2 工况 A — Mode I 单调张开

- 顶边 `u_y` 从 `0` 增至 `> δ_c / Λ`。
- 记录界面平均 `(δ_n, T_n)`，与解析双线性律比对：
  - 弹性：`T = K_n · δ`（`δ < δ_0`）
  - 软化：线性降至 `(δ_c, 0)`
  - 断裂：`T ≈ 0`，`D ≈ 1`
- 容差：弹性段目标 `rtol ≈ 1e-3`；软化/断裂段可略松（多界面 + 体弹性）。

### 5.3 工况 B — 卸载 / 重加载

- 加载到 `δ_0 < δ < δ_c` → 卸至 0 → 再加载。
- 验收：卸载斜率 ≈ `K_n (1−D)`；重加载至历史 `δ_max` 前无新损伤；`D` 单调不减。

### 5.4 工况 C — Mode II 切向

- 顶边相对切向 `u_x`，`u_y = 0`（法向分离保持 ≈0）。
- 比对 `(δ_t, T_t)` 与切向双线性律。

### 5.5 输出

- `@testset` 硬断言；失败非零退出。
- 可选：`output/unit_czm_bilinear_Tdelta.png`。

## 6. 脚本 2：`test/unit_czm_eigenstrain.jl`

### 6.1 驱动（合成载荷）

不调用 SPMe / 热求解器。伪时间 `t ∈ [0, 1]` 上给定斜坡：

- `ΔT(t)`、`Δsoc_n(t)`、`Δsoc_p(t)`
- 幅值取 Jellyroll 量级，并**强制全程弹性**（`δ < δ_0`），以便与闭式解比对。
- **固定端约束下必须使 `Σ h_i ε₀,i < 0`（推荐冷却 `ΔT≤0`）**，以产生拉伸张开；正膨胀仅导致压缩接触，无法验证 Mode I 分离。

本征应变走生产契约 `assemble_thermal_chemical_load`（**不要**在脚本里另写按材料分支的 ε₀）：

```
ε₀[e] = α_eff * dT_elem[e] + β_n * Δsoc_n_elem[e] + β_p * Δsoc_p_elem[e]
```

- `α_eff`、`β_n = NE.Omega/3`、`β_p = PE.Omega/3`：与生产相同的标量入口（见 `CouplingState` / `Mechanical.jl`）。
- 按层填充长度-8 向量（与生产 `compute_czm_strain_inputs` 粒度一致）：
  - 全体层：`dT_elem[e] = ΔT(t)`
  - 仅 NE 层：`Δsoc_n_elem[e] = Δsoc_n(t)`，其余 0
  - 仅 PE 层：`Δsoc_p_elem[e] = Δsoc_p(t)`，其余 0
- 1D 解析解使用**同一** `ε₀[e]` 定义，保证比对对象一致。
- 不硬编码物理常数；从 `case.param` / `param_cache` 读取。

### 6.2 边界与求解

- 底边固定 + 顶边固定（模拟内外圈约束）→ 本征膨胀转化为界面张开。
- `assemble_coupled_system_full` + 自定义 BC Newton；`dT_elem` / `Δsoc_*_elem` 长度 = bulk 单元数（8）。

### 6.3 1D 解析解（仅弹性段）

将条带化为：8 段弹性杆 + 4 个法向弹簧（`K_n` 按界面类型）。

固定端约束示意：

```
Σ_i h_i (ε₀,i + σ_i / E_i*) + Σ_j δ_j = 0
界面：T_j = σ = K_n,j · δ_j（弹性）
```

平面应力下 `E*` 与 FEM `moduli_of` 一致。输出各界面 `δ_n^{anal}`、`T_n^{anal}`，与 FEM 比对；目标相对误差量级 `rtol ≈ 1e-2`（允许泊松/二维效应偏差）。

### 6.4 验收与输出

- `@testset`：4 界面分离相对误差；`D ≈ 0`；无 NaN。
- 可选：`output/unit_czm_eigenstrain_delta.png`（FEM vs 解析）。

### 6.5 与原需求标题的关系

原计划标题写「施加电化学 SPMe 和热集总模型」。本规格采用**合成本征应变斜坡**验证同一载荷管线（`assemble_thermal_chemical_load`），避免单元测试依赖完整电化学-热时间积分。端到端耦合仍由 smoke 覆盖。

## 7. 错误处理

| 情况 | 行为 |
|------|------|
| 网格硬断言失败 | 立即 `error`，附断言信息 |
| Newton 不收敛 | 测试失败；打印步号、残差、当前 `δ`/`D` |
| 解析比对超容差 | `@test` 失败；报告界面 id 与 FEM/解析值 |

## 8. 文件清单

| 路径 | 作用 |
|------|------|
| `src/CzmUnitMesh.jl` | `create_unit_czm_strip` |
| `src/JuBat.jl` | `include` 新文件 |
| `test/unit_czm_bilinear.jl` | 脚本 1 |
| `test/unit_czm_eigenstrain.jl` | 脚本 2 |
| `docs/superpowers/specs/2026-07-23-unit-czm-strip-verification-design.md` | 本规格 |
| `docs/planning-with-files/15_SPMe-热-内聚力耦合模型单元级验证/SPMe-热-内聚力耦合模型单元级验证.md` | 计划摘要（指向本规格） |

## 9. 验收清单

1. `create_unit_czm_strip` → 8 Q4 + 4 COH2D4，硬断言全过。
2. 脚本 1：Mode I 全曲线、卸载重载、Mode II 均过 `@testset`。
3. 脚本 2：弹性段分离与 1D 解析解在约定 `rtol` 内；`D≈0`。
4. 两脚本均可独立运行：`julia --project=. test/unit_czm_bilinear.jl` 与 `.../unit_czm_eigenstrain.jl`。
5. 不修改生产 `create_czm_mesh` / 双线性本构核心逻辑（仅新增与复用）。

## 10. 测试组织

- 位置：`test/`（可回归）；网格生成器在 `src/`（可复用）。
- 运行：独立脚本即可；不强制改全局 test runner。
- 出图：可选，写入 `output/`，失败不得依赖出图。
