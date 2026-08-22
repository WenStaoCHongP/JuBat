# 堆芯塌陷力学建模 Batch 2'（卷绕预应力 / opt-in）实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 实现 `czm_winding_prestress=true` 的卷绕预应力初始应力场 σ₀(r)（自平衡重分布经首个平衡求解），叠加进 GL 机械残差与 K_G，参数缺失显式 error，默认关零漂移；通过 spec §7 Batch 2' 验收门（默认基线逐位不变 + σ₀ 量级对照解析卷绕公式）。

**Architecture:** 两任务。Task 1 参数与 σ₀ 场：`Cell` 新增 `winding_T_ne/winding_T_pe`（Pa，0=未设置）、`NormaliseParam` `/σ_czm` 归一、`CzmPlasticity.jl` 旁新建 `winding_prestress_field(czm_mesh, param)`（逐单元全局系 σ₀，解析闭式）+ 单元测试（解析对照/量级/缺参）。Task 2 接入与门禁：`gl_element_residual_tangent` 逐单元 σ₀ 叠加（零时走原路径保逐位）、`prestress` kwarg 与 eigenstrain 同模式贯穿装配/求解器/生产调用点、`CzmLayout` 持久化 + 结果标记、集成测试（默认关逐位、开启收敛+K_G 方向性）、三道门禁与记录。

**Tech Stack:** Julia 1.11.2；`Test`/`LinearAlgebra`/`SparseArrays`。

**Spec:** `docs/superpowers/specs/2026-08-20-core-collapse-mechanics-design.md`（v1.3）§3.7（D11）、§4.1、§7 Batch 2' 行；参数依据：`宏观电池模型_正负极极片力学参数_文献整理.md` §10.4.1（T_w 负极侧 1–5 MPa、正极侧 0.5–2 MPa、卷芯预压 p0 0.2–1.0 MPa；工艺文献 + Snyder 2025 [14] + Anisotropic 2025 [14]）。

## Global Constraints

（同 Batch 2/3 计划：运行环境、强制行为基线 PNG SHA `4ba6207c…e932`、兼容性契约、AGENTS 9.7 错误处理、每 Task 一提交。）

### 设计决策（执行契约）

**D-B2'-1（σ₀ 场闭式模型冻结）**：按"卷入张力等应变分担 + 外层累积对数压力"两件套（薄层缠绕经典模型，与 §10.4.1 参数组一一对应）：
1. 卷入张力应变（每侧 web 等应变/Voigt 分担）：`ε_w,side = T_side / Ē_side`，`Ē_side = Σ_{k∈side}(E_k·t_k)/Σ t_k`（NE 侧 = {NE,NCC}，PE 侧 = {PE,PCC}，SP 取两侧均值；E_k 用 `moduli_of` 归一模量，t_k 归一厚度）。层 i 的卷入环向张力 `σ_t,i = E_i·ε_w,side(i)`。
2. 累积接触压力（外层缠绕压内层，经典 `p(r)=f·ln(R_end/r)`）：`f = (T_ne+T_pe)/t_repeat`（归一量），`R_end` 取网格最大节点半径；`σ_r0(r) = −p(r)`，`σ_θ0,i = σ_t,i − p(r̄_i)`（r̄_i 为单元质心半径）。
3. 逐单元变换到全局 (x,y)：`t̂=(−sinθ,cosθ)`、`n̂=(cosθ,sinθ)`，`σ_xx = σ_θ t̂x²+σ_r n̂x²` 等三分量（σ_xθ 项含 2t̂·n̂）。输出 `Vector{NTuple{3,Float64}}`（σ_czm 单位）。
**D-B2'-2（需 geo_nl）**：`czm_winding_prestress=true` 要求 `czm_geo_nonlinear=true`（K_G 消费初应力；线弹性常刚度槽位不承载初应力项，与 D-B3-1 同理）。`prestress` 与 `eigenstrain` 同为求解链 kwarg（`prestress=nothing` 缺省）。
**D-B2'-3（零值旁路保逐位）**：σ₀=(0,0,0) 时 `gl_element_residual_tangent` 走原路径（不执行 +0.0 加法），保证默认关时 GL 路径与 Batch 2/3 冻结输出逐位一致（规避 −0.0+0.0 类位级扰动）。
**D-B2'-4（首平衡重分布语义）**：σ₀ 不要求先自平衡——首个平衡求解（零本征应变步）按 BC 约束完成一致性重分布（spec §3.7 原文语义）；环向张力由外圈约束反力平衡（`czm_fix_inner=true` 默认内外圈固定）。

### 归一化契约

`Cell.winding_T_ne/pe`（Pa）在 `NormaliseParam` cell 段 `/σ_czm` 归一；`t_repeat` 用 `param.cell.layer`（归一值 2.162）；单元厚度 `t_k` 用 `param.{PE,NE,SP,PCC,NCC}.thickness`（归一）。

## File Structure

| 文件 | 动作 | 职责 |
|---|---|---|
| `src/SetParams.jl` | 修改 | `Cell` 两新字段（默认 0=未设置）；`NormaliseParam` cell 段两行 |
| `src/parameters/Jellyroll.jl` | 修改 | 显式 `cell.winding_T_ne = 3.0e6`（工艺典型中值）、`cell.winding_T_pe = 1.0e6`（§10.4.1） |
| `src/CzmPlasticity.jl` | 修改 | 新增 `winding_prestress_field`（紧邻塑性工具；同文件聚合 Batch 2'/3 机械工具） |
| `src/czm.jl` | 修改 | `gl_element_residual_tangent` σ₀ 叠加（零值旁路）；`assemble_bulk_residual_tangent` `prestress` kwarg（校验：与 plasticity 同关时 mech_state 规则不变；prestress 要求 geo_nl） |
| `src/CzmSolve.jl` | 修改 | `prestress=nothing` kwarg 贯穿（dispatch/basic/newton/backtrack，并入 eig_kwargs 元组） |
| `src/CouplingState.jl` | 修改 | `CzmLayout.winding_prestress::Union{Nothing,Vector{NTuple{3,Float64}}}`；`update_czm_damage!` 校验（需 geo_nl、T_ne/T_pe 至少一侧>0 否则 error 指明）、惰性计算并持久化 σ₀ 场、结果标记 |
| `src/Solve.jl` 或结果装配处 | 修改 | `winding_prestress` 开启时结果字典写入 `"winding prestress" = true`（与 `czm D_max` 键同站点；执行时 `rg -n "czm D_max" src/` 定位） |
| `test/test_czm_winding_prestress.jl` | 新建 | spec §7 Batch 2' 全部门禁 |

---

## Task 1: 参数与 σ₀ 场函数（TDD）

- [ ] **Step 1: 写失败测试** `test/test_czm_winding_prestress.jl`（testset ①②③）：

```julia
@testset "σ₀ 场解析对照（卷绕公式逐项）" begin
    # 构造 nθ=8 夹具（同 build_c1_fixture 模式），T_ne=3e6/T_pe=1e6（有量纲）经 NormaliseParam
    # 取 3 个不同半径的 PCC/NCC/NE 单元 e，手工按 D-B2'-1 公式计算期望 σ_θ0/σ_r0（归一单位），
    # 断言 winding_prestress_field 输出的全局分量经 (t̂,n̂) 旋转回 (θ,r) 后逐项 ≈（rtol 1e-12）
end
@testset "量级校验（§10.4.1 工艺区间）" begin
    # σ₀ 场环向应力绝对值 max ∈ [0.1, 10] MPa 等效（换算回有量纲：×σ_czm），卷芯 p(a) ∈ [0.1, 3] MPa
end
@testset "缺参/非法即 error（AGENTS 9.4/9.7）" begin
    # T_ne=T_pe=0 时 winding_prestress_field(czm_mesh, param) → ErrorException（指明 winding_T 未设置）
    # 负张力 → error（物理非法）
end
```

- [ ] **Step 2**: 运行确认 `UndefVarError: winding_prestress_field`。
- [ ] **Step 3**: 实现——`Cell` 字段 + `NormaliseParam`（cell 段，`param.cell.winding_T_ne = param_dim.cell.winding_T_ne / param.scale.σ_czm` 两行）+ Jellyroll.jl 显式值 + `winding_prestress_field`（完整 ~35 行：质心/半径、Ē_side、ε_w、p(r) 对数式、(θ,r)→(x,y) 旋转；参数校验前置）。测试体内给出手工期望值的计算代码（非魔法数）。
- [ ] **Step 4**: 3/3 通过 + `-e include` 干净 + 提交 `feat(czm): 卷绕预应力 σ₀ 场（等应变张力分担 + 对数累积压力，D-B2'-1）与参数`。

## Task 2: 接入、链路、持久化与集成测试（TDD）

- [ ] **Step 1: 追加失败测试**（testset ④⑤⑥）：

```julia
@testset "默认关逐位不变（回归锚）" begin
    # geo_nl=true、prestress 关：与 Batch 2 冻结 geo 解逐位（同载荷两次求解 == 断言，
    # 一次经新代码路径 prestress=nothing）
end
@testset "开启后首平衡收敛 + K_G 方向性" begin
    # 零本征应变 + prestress：求解收敛；u 场有限且非全零（重分布变形发生）；
    # 方向性：构造 n=2 椭圆化位移场 d = A·cos(2θ)·r̂（A=1e-3·R），
    # 断言 dot(d, K_prestress·d) < dot(d, K_0·d)（hoop 压缩降低椭圆化切线能量，
    # K 取 assemble_bulk_residual_tangent(u≈0) 两次带/不带 prestress）
end
@testset "求解器链路 + 结果标记" begin
    # solve_czm_step(...; geo_nl, eigenstrain, prestress=σ₀场) 收敛且与装配级一致；
    # update_czm_damage! 生产路径：opt 开 + T 缺 → error；开后 case.czm_layout.winding_prestress 非 nothing
end
```

- [ ] **Step 2**: 确认失败（prestress 关键字不存在）。
- [ ] **Step 3**: 实现——`gl_element_residual_tangent` 增 `σ0::NTuple{3,Float64}=(0.0,0.0,0.0)` 位置后参数：非零时 `S += σ0`（残差与 K_G 的 Ŝ 均用总应力；材料切线不变）；零值旁路（D-B2'-3）。`assemble_bulk_residual_tangent` 增 `prestress=nothing` kwarg（要求 geo_nl；逐单元切片传入）。`assemble_coupled_system`/`solve_czm_step`/basic/newton/backtrack：`prestress=nothing` 并入 eig_kwargs 元组与 backtrack 显式 kwargs（与 eigenstrain 完全同模式）。`CzmLayout` 增字段；`update_czm_damage!`：校验三元组（需 geo_nl、T 至少一侧>0、plasticity 组合合法）+ 惰性计算持久化 + 传参；结果标记键。
- [ ] **Step 4**: 全部 6 testset 通过；定向回归（geo_c1/geometric_stiffness/mech_core/unit_czm_j2/j2_integration/assemble_coupled_system）零失败。
- [ ] **Step 5**: 提交 `feat(czm): 卷绕预应力接入 GL 残差与 K_G 全链路（CzmLayout 持久化、结果标记、D-B2'-3 零值旁路）`。

## Task 3: 三道门禁与批次记录

- [ ] 全套 29/29；探针快照逐位（默认关路径 D-B2'-3 保证）；testexample 冻结指标 + PNG SHA 一致。
- [ ] `Simplify/baseline.md` 追加行；progress 追加 Batch 2' 小节（D-B2'-1..4 执行情况、量级校验实测值）；index 更新。
- [ ] 提交 `docs(baseline): Batch 2' 门禁记录`。

---

## 自评审记录

1. **spec 覆盖**：§3.7 全条目（opt-in/σ₀ 叠加残差与 K_G/缺参 error/首平衡重分布/结果标记/三验收）→ Task 1（场+量级）+ Task 2（接入+方向性+默认关逐位+标记）；§4.1 parameters 行 ✓；§7 Batch 2' 三测试 → testset ①②（解析+量级）、回归锚、④⑤ ✓。
2. **决策登记**：D-B2'-1（闭式模型——比 spec 点名的"Altmann 型"更简：等应变分担 + 对数压力均为薄层缠绕经典式，与 §10.4.1 参数组直接对接；差异在计划批准时冻结）；D-B2'-2（需 geo_nl）；D-B2'-3（零值旁路保逐位）；D-B2'-4（首平衡重分布）。
3. **类型一致性**：σ₀ 场 `Vector{NTuple{3,Float64}}` 在场函数/装配/布局三处一致；kwarg 名 `prestress` 与 eigenstrain 模式一致。
4. **无占位符**：testset ①③ 的期望值由测试内代码计算（公式同实现但独立编码），非魔法数；结果标记站点给出定位命令。
5. **风险**：首平衡收敛性（140× 刚度对比 + σ₀ 非自平衡）——load_substep 自适应子步兜底，测试用小 T（1 MPa 档）先验证收敛再上量级；若仍难收敛，缩小步长下限已内建（step/128）。
