# ThermalDistributed.jl 优化方案

> 日期: 2026-04-01
> 文件: `src/ThermalDistributed.jl`
> 状态: 新增 (425 行)

---

## 1. 现状分析

### 1.1 函数清单

| 函数 | 行范围 | 行数 | 职责 |
|------|--------|------|------|
| `ThermalDistributed2D` | 1-47 | 47 | 2D FEM 热刚度/容量矩阵装配 |
| `apply_convection_bc` | 49-104 | 56 | Newton 对流边界条件 |
| `apply_cool_method` | 107-185 | 79 | 冷却方式 (none/tab/surface) |
| `ThermalDistributed2D_BC` | 187-219 | 33 | BC 总入口 (含 CZM 界面热阻) |
| `ThermalDistributed2D_Ring` | 221-270 | 50 | 极坐标 FEM 热模型 |
| `ThermalRing2D_BC` | 272-279 | 8 | 极坐标 BC |
| `compute_heat_sources` | 281-391 | 111 | **11 层热源逐一计算** |
| `compute_heat_sources_with_czm` | 393-425 | 33 | CZM 活跃单元过滤 |

### 1.2 核心问题

1. **`compute_heat_sources` 是第二大的上帝函数 (111 行)**：11 层热源逐一计算 + 单位转换，逻辑冗长且与 FEM 装配无关
2. **`compute_heat_sources_with_czm` 是简单包装**：仅调用 `compute_heat_sources` 后用 `active_mask` 过滤，可内联
3. **热源计算与 FEM 装配职责混杂**：ThermalDistributed.jl 应只负责 FEM 装配和 BC，热源计算应独立

---

## 2. 优化方案

### 2.1 拆分 `compute_heat_sources` → ThermalHeatSource.jl

热源计算从 ThermalDistributed.jl 完全迁出到新文件 `src/ThermalHeatSource.jl`。

### 2.2 热源函数分层

```julia
# ThermalHeatSource.jl (新文件, ~180 行)

"""
计算所有 11 层逐单元热源 (无量纲)。

# 返回
- `q_total::Vector{Float64}`: ne 个单元的总热源
- `q_layers::NamedTuple`: 11 层热源 (各 ne 个值)
"""
function compute_heat_sources(case::Case, variables::Dict,
                               A_elem::Vector{Float64},
                               I_e::Vector{Float64},
                               Te_prev::Vector{Float64},
                               T_nodes::Vector{Float64})
    ne = length(I_e)
    param = case.param
    scale = case.param_dim.scale

    # 预提取共享量
    T_ref = scale.T_ref
    j_scale = scale.j
    phi_scale = scale.phi

    q_rxn_ne = zeros(ne)
    q_rev_ne = zeros(ne)
    q_ohm_s_ne = zeros(ne)
    q_ohm_e_ne = zeros(ne)
    q_sp = zeros(ne)
    q_rxn_pe = zeros(ne)
    q_rev_pe = zeros(ne)
    q_ohm_s_pe = zeros(ne)
    q_ohm_e_pe = zeros(ne)
    q_pcc = zeros(ne)
    q_ncc = zeros(ne)

    for e in 1:ne
        q_rxn_ne[e], q_rev_ne[e], q_ohm_s_ne[e], q_ohm_e_ne[e],
        q_sp[e],
        q_rxn_pe[e], q_rev_pe[e], q_ohm_s_pe[e], q_ohm_e_pe[e],
        q_pcc[e], q_ncc[e] = _compute_element_heat_source(
            e, variables, param, I_e[e], Te_prev[e], A_elem[e],
            T_ref, j_scale, phi_scale
        )
    end

    q_total = q_rxn_ne + q_rev_ne + q_ohm_s_ne + q_ohm_e_ne +
              q_sp + q_rxn_pe + q_rev_pe + q_ohm_s_pe + q_ohm_e_pe +
              q_pcc + q_ncc

    return q_total, (;
        q_rxn_ne, q_rev_ne, q_ohm_s_ne, q_ohm_e_ne,
        q_sp,
        q_rxn_pe, q_rev_pe, q_ohm_s_pe, q_ohm_e_pe,
        q_pcc, q_ncc
    )
end

"""单个单元的 11 层热源"""
function _compute_element_heat_source(e, variables, param, I_e, T_e, A_e,
                                       T_ref, j_scale, phi_scale)
    # 从 variables 中提取单元 e 的物理量
    eta_n = variables["thermal2D eta_n_e"][e]
    eta_p = variables["thermal2D eta_p_e"][e]
    dUdT_n = variables["thermal2D dUdT_n_e"][e]
    dUdT_p = variables["thermal2D dUdT_p_e"][e]
    sigma_eff = variables["thermal2D element conductivity"][e]

    # NE 反应热 + 可逆热
    q_rxn_ne = abs(I_e) * abs(eta_n) / A_e
    q_rev_ne = abs(I_e) * T_e * abs(dUdT_n) / A_e

    # NE 欧姆热 (固体 + 电解液)
    q_ohm_s_ne = I_e^2 * param.NE.Rs / (A_e * param.NE.thickness)
    q_ohm_e_ne = I_e^2 * param.EL.R_e_ne / (A_e * param.NE.thickness)

    # 隔膜欧姆热
    q_sp = I_e^2 * param.SP.R_sp / (A_e * param.SP.thickness)

    # PE 反应热 + 可逆热
    q_rxn_pe = abs(I_e) * abs(eta_p) / A_e
    q_rev_pe = abs(I_e) * T_e * abs(dUdT_p) / A_e

    # PE 欧姆热
    q_ohm_s_pe = I_e^2 * param.PE.Rs / (A_e * param.PE.thickness)
    q_ohm_e_pe = I_e^2 * param.EL.R_e_pe / (A_e * param.PE.thickness)

    # 集流体
    q_pcc = I_e^2 * param.cell.R_pcc / (A_e * param.cell.t_pcc)
    q_ncc = I_e^2 * param.cell.R_ncc / (A_e * param.cell.t_ncc)

    return q_rxn_ne, q_rev_ne, q_ohm_s_ne, q_ohm_e_ne,
           q_sp,
           q_rxn_pe, q_rev_pe, q_ohm_s_pe, q_ohm_e_pe,
           q_pcc, q_ncc
end

"""CZM 过滤版热源"""
function compute_heat_sources_with_czm(case, variables, A_elem, I_e, Te_prev, T_nodes)
    q_total, q_layers = compute_heat_sources(case, variables, A_elem, I_e, Te_prev, T_nodes)

    if case.czm_mesh !== nothing
        active = get_active_elements(case.czm_mesh, case.layout.ne)
        q_total .*= active
        # 过滤每层
        for k in keys(q_layers)
            getfield(q_layers, k) .*= active
        end
    end

    return q_total, q_layers
end
```

### 2.3 ThermalDistributed.jl 精简后

```
重构前 (425 行):                    重构后 (~220 行):
─────────────────────               ─────────────────────
ThermalDistributed2D     47行  →    保留                    47行
apply_convection_bc      56行  →    保留                    56行
apply_cool_method        79行  →    保留                    79行
ThermalDistributed2D_BC  33行  →    保留                    33行
ThermalDistributed2D_Ring 50行 →    保留                    50行
ThermalRing2D_BC          8行  →    保留                     8行
compute_heat_sources    111行  →    迁出 → ThermalHeatSource.jl
compute_heat_sources_... 33行  →    迁出 → ThermalHeatSource.jl
───────────────────────────────────────────────────────────
总计                    425行      ~275行 (FEM+BC) + ~180行 (新文件)
```

---

## 3. 不动的部分

| 函数 | 原因 |
|------|------|
| `ThermalDistributed2D` | FEM 装配核心，结构良好 |
| `apply_convection_bc` | BC 逻辑完整 |
| `apply_cool_method` | 冷却方式完整 |
| `ThermalDistributed2D_BC` | BC 入口（仅调用 `compute_heat_sources` 时需改 import） |
| `ThermalDistributed2D_Ring` | 极坐标路径独立 |
| `ThermalRing2D_BC` | 极坐标 BC |

---

## 4. 预期效果

| 指标 | 旧 | 新 |
|------|-----|-----|
| ThermalDistributed.jl 行数 | 425 | ~275 |
| ThermalHeatSource.jl 行数 | 0 | ~180 |
| 热源计算独立性 | 耦合在 FEM 文件中 | 完全独立 |
| `compute_heat_sources` 行数 | 111 | ~40 (入口) + ~50 (单元素) |
