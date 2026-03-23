# 热源计算代码重构设计

## 背景

Solve.jl 第633-744行与 ThermalDistributed.jl 中的热源计算函数存在功能重叠：
- Solve.jl 中逐单元计算11个分层热源分量
- ThermalDistributed.jl 中 `compute_element_heat_sources` 使用全局电化学变量
- 量纲转换复杂：无量纲 → 物理 → 无量纲
- 物理单位转换分散在求解代码中

## 设计目标

1. 消除代码冗余：合并 Solve.jl 与 ThermalDistributed.jl 中的热源计算逻辑
2. 简化量纲转换：计算过程保持无量纲，后处理统一转换
3. 统一接口：支持单 SPMe 和多 SPMe 模式

## 设计方案

### 1. Variables.jl 修改

将现有的带单位热源变量改为无量纲版本：

**修改位置**：第93-130行

**变量命名变更**：
| 原变量名 | 新变量名 |
|---------|---------|
| `thermal2D Q_rxn_NE [W/m3]` | `thermal2D q_rxn_ne` |
| `thermal2D Q_rev_NE [W/m3]` | `thermal2D q_rev_ne` |
| `thermal2D Q_ohm_s_NE [W/m3]` | `thermal2D q_ohm_s_ne` |
| `thermal2D Q_ohm_e_NE [W/m3]` | `thermal2D q_ohm_e_ne` |
| `thermal2D Q_SP [W/m3]` | `thermal2D q_sp` |
| `thermal2D Q_rxn_PE [W/m3]` | `thermal2D q_rxn_pe` |
| `thermal2D Q_rev_PE [W/m3]` | `thermal2D q_rev_pe` |
| `thermal2D Q_ohm_s_PE [W/m3]` | `thermal2D q_ohm_s_pe` |
| `thermal2D Q_ohm_e_PE [W/m3]` | `thermal2D q_ohm_e_pe` |
| `thermal2D Q_PCC [W/m3]` | `thermal2D q_pcc` |
| `thermal2D Q_NCC [W/m3]` | `thermal2D q_ncc` |

### 2. ThermalDistributed.jl 修改

**新增函数**：`compute_heat_sources`

```julia
function compute_heat_sources(case::Case, variables::Dict,
                              variables_elems::Union{Vector{<:Dict}, Nothing},
                              I_e::Vector{Float64}, T_e::Vector{Float64},
                              areas::Vector{Float64}; per_element_spme::Bool=false)
```

**功能**：
- 根据参数 `per_element_spme` 判断使用全局变量还是逐单元变量
- 计算各层热源分量（无量纲）
- 直接修改 `variables` 字典中的预分配数组
- 返回修改后的 `variables`

**删除函数**：
- `heatQ_Source`（被新函数替代）
- `compute_element_heat_sources`（被新函数替代）
- `heatQ_Source_with_czm`（功能合并到新函数，通过 `czm_mesh` 参数控制）

### 3. Solve.jl 修改

**修改位置**：第633-744行

**修改内容**：
- 删除约110行热源计算代码
- 替换为调用 `compute_heat_sources` 函数
- 保留辅助变量保存（eta_n_e, eta_p_e, dUdT_n_e, dUdT_p_e, soc_n, soc_p）

### 4. PostProcessing.jl 修改

**添加内容**：热源物理单位转换

```julia
if case.opt.thermalmodel == "distributed2D" && case.opt.model == "SPMe"
    q_scale = case.param_dim.scale.q

    result["thermal2D Q_rxn_NE [W/m3]"] = variables["thermal2D q_rxn_ne"][:, 1:v] .* q_scale
    result["thermal2D Q_rev_NE [W/m3]"] = variables["thermal2D q_rev_ne"][:, 1:v] .* q_scale
    # ... 其他热源分量 ...
end
```

## 数据流

```
Solve.jl                    ThermalDistributed.jl                PostProcessing.jl
    │                              │                              │
    ├─ variables_elems ────────► compute_heat_sources ──► variables
    │                              │                              │
    │                              ├─ 填充无量纲热源              │
    │                              │                              │
    └──────────────────────────────────────────────────────────► 转换为物理单位
                                                               │
                                                               ▼
                                                        result["Q_xxx [W/m3]"]
```

## 量纲转换简化

**重构前**（复杂）：
```
q_elem（无量纲）→ q_elem .* q_ec_scale（物理）→ ./ q_ref（无量纲）
```

**重构后**（简洁）：
```
q_elem（无量纲）→ PostProcessing 中 .* q_scale（物理单位输出）
```

## 影响范围

| 文件 | 改动量 | 风险 |
|------|--------|------|
| Variables.jl | 修改变量名 | 低（变量重命名） |
| ThermalDistributed.jl | 新增函数，删除旧函数 | 中（核心逻辑） |
| Solve.jl | 删除约110行，替换为函数调用 | 中（调用点） |
| PostProcessing.jl | 添加热源转换 | 低（新增代码） |

## 测试计划

1. 运行 `example/testexample.jl` 验证多 SPMe 模式热源计算
2. 运行 `example/SPMe_Thermal_example.jl` 验证单 SPMe 模式热源计算
3. 对比重构前后热源输出值（物理单位）应一致
4. 检查循环仿真中热源累积是否正确

## 回滚计划

如果出现问题，可以通过 git 回滚：
```bash
git checkout -- src/Variables.jl src/ThermalDistributed.jl src/Solve.jl src/PostProcessing.jl
```
