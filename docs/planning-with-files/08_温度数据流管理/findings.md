# Findings & Decisions

## 问题 1: 9 个死内存分配（预分配但从未写入）

以下变量在 `Variables.jl` 的 `StandardVariables` 中预分配，但在整个 `src/` 中从未被写入或读取：

| 变量键 | Variables.jl 行号 | 维度 | 估计浪费 |
|--------|-------------------|------|----------|
| `thermal2D temperature` | 97 | (ne, num) | ~0.6-1.6 MB |
| `thermal2D temperature  at nodes history` | 100 | (nT, num) | ~0.6-1.6 MB |
| `thermal2D displacement x` | 133 | (nT, num) | ~0.6-1.6 MB |
| `thermal2D displacement y` | 134 | (nT, num) | ~0.6-1.6 MB |
| `thermal2D element thermal stress` | 123 | (ne, num) | ~0.6-1.6 MB |
| `thermal2D element diffusion stress` | 124 | (ne, num) | ~0.6-1.6 MB |
| `thermal2D element total stress` | 125 | (ne, num) | ~0.6-1.6 MB |
| `thermal2D element diffusion strain` | 126 | (ne, num) | ~0.6-1.6 MB |
| `thermal2D element thermal strain` | 127 | (ne, num) | ~0.6-1.6 MB |

- **总计**: 9 个数组，~5-15 MB 无用内存/仿真
- **`thermal2D temperature  at nodes history`** 有双空格 bug（行 100）
- `thermal2D displacement x/y` 可能是为 CZM 位移预留的占位符，但实际 CZM 位移写入的是 `"czm displacement x/y"`

## 问题 2: 输出路径不一致

### PostProcessing.jl 还原的热变量（统一入口）
| 变量键 | 还原公式 | 行号 |
|--------|---------|------|
| `thermal2D q_rxn_ne` → `thermal2D Q_rxn_NE [W/m3]` | × scale.q | 52 |
| `thermal2D q_rev_ne` → `thermal2D Q_rev_NE [W/m3]` | × scale.q | 53 |
| `thermal2D q_ohm_s_ne` → `thermal2D Q_ohm_s_NE [W/m3]` | × scale.q | 54 |
| `thermal2D q_ohm_e_ne` → `thermal2D Q_ohm_e_NE [W/m3]` | × scale.q | 55 |
| `thermal2D q_sp` → `thermal2D Q_SP [W/m3]` | × scale.q | 56 |
| `thermal2D q_rxn_pe` → `thermal2D Q_rxn_PE [W/m3]` | × scale.q | 57 |
| `thermal2D q_rev_pe` → `thermal2D Q_rev_PE [W/m3]` | × scale.q | 58 |
| `thermal2D q_ohm_s_pe` → `thermal2D Q_ohm_s_PE [W/m3]` | × scale.q | 59 |
| `thermal2D q_ohm_e_pe` → `thermal2D Q_ohm_e_PE [W/m3]` | × scale.q | 60 |
| `thermal2D q_pcc` → `thermal2D Q_PCC [W/m3]` | × scale.q | 61 |
| `thermal2D q_ncc` → `thermal2D Q_NCC [W/m3]` | × scale.q | 62 |
| `thermal2D temperature at nodes` → `thermal2D temperature at nodes [K]` | × scale.T_ref | 63 |

### Solve.jl 直接输出（绕过 PostProcessing）
| 变量键 | 写入位置 | 还原公式 | 行号 |
|--------|---------|---------|------|
| `thermal2D element current` | Solve.jl:392-401 | 无量纲直传 | 392 |
| `thermal2D eta_n_e` | Solve.jl:392-401 | 无量纲直传 | 392 |
| `thermal2D eta_p_e` | Solve.jl:392-401 | 无量纲直传 | 392 |
| `thermal2D element soc_n` | Solve.jl:392-401 | 无量纲直传 | 392 |
| `thermal2D element soc_p` | Solve.jl:392-401 | 无量纲直传 | 392 |
| `thermal2D dUdT_n_e` | Solve.jl:392-401 | 无量纲直传 | 392 |
| `thermal2D dUdT_p_e` | Solve.jl:392-401 | 无量纲直传 | 392 |
| `thermal2D element voltages` | Solve.jl:392-401 | 无量纲直传 | 392 |
| `thermal2D element OCV` | Solve.jl:396 | 无量纲直传 | 396 |
| `thermal2D n_cutoff_elements` | Solve.jl:396 | 无量纲直传 | 396 |
| `thermal2D active_mask` | Solve.jl:396 | 无量纲直传 | 396 |
| `heat_source_fields` | Solve.jl:400 | 无量纲直传 | 400 |
| `total heat source` → `total heat source [W]` | Solve.jl:401 | × scale.L³ | 401 |

**问题**: 两条输出路径导致：
- 热源分量通过 PostProcessing 统一还原（× scale.q）
- 热源总量和其他变量通过 Solve.jl 直接输出（部分无量纲直传、部分手动 × scale.L³）
- 用户拿到结果后无法区分哪些已还原、哪些仍是归一化值

## 问题 3: 内部耦合变量（非 bug，但值得关注）

以下变量在 StandardVariables 中有完整历史预分配，但仅用于模块间耦合，历史数据不一定需要：

| 变量键 | 写入位置 | 读取位置 |
|--------|---------|---------|
| `thermal2D element area` | CallModel.jl:29,43,58 | 内部耦合 |
| `thermal2D element soc_n/soc_p` | CallModel.jl:136-137 | Mechanical.jl, CouplingState.jl, CycleData.jl |
| `heat_source_fields` | ThermalDistributed.jl:495 | ThermalPolar2D.jl, ThermalDistributed.jl |

这些变量确实有历史输出的价值（用于分析温度/SOC/热源的空间分布演化），暂不改动。

## 与 CZM 模块问题对比

| 问题类型 | CZM 模块 | 温度模块 |
|---------|---------|---------|
| 死内存分配 | 0 个 | **9 个**（更严重） |
| 键名不一致 | D_max vs max damage | 无（PostProcessing 用显示名，内部键一致） |
| 预分配缺失 | 6 个场变量未预分配 | 0 个（全部预分配，甚至过多） |
| 输出路径分散 | czm_output_to_variables + PostProcessing | Solve.jl 直接输出 + PostProcessing |
| Bug | 无 | 双空格键名 |

## Resources
- src/Variables.jl (变量定义, 行 88-135)
- src/PostProcessing.jl (还原, 行 49-64)
- src/Solve.jl (直接输出, 行 392-401)
- src/ThermalDistributed.jl (热源写入)
- src/CallModel.jl (element 变量写入)

---
*Update this file after every 2 view/browser/search operations*
