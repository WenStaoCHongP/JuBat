# 网格敏感性分析脚本评审

> 日期: 2026-04-22
> 分支: Parameters_Design (vs main)
> 评审范围: `example/网格敏感性/1_cohesive_characteristic_length.jl`, `example/网格敏感性/2_electrochemical_mesh_sensitivity.jl`, `example/网格敏感性/3_thermal_mesh_sensitivity.jl`, `example/网格敏感性/4_czm_mesh_sensitivity.jl`, `example/网格敏感性/5_energy_conservation_check.jl`
> 评审类型: 物理口径与输出一致性审查

---

## 1. 总体结论

这组脚本的结构是完整的，也通过了静态语法检查，但它们当前实现的物理口径和我们前面确认的网格敏感性方案并不一致。主要问题不是“能不能跑”，而是“跑出来的指标是否真的对应计划里的判据”。

**总体判断**: 目前这批脚本更像是一个探索性草稿，而不是可以直接作为网格敏感性结论依据的正式实现。

---

## 2. 严重问题

### P1：热网格脚本没有实现“固定热源基准”

`example/网格敏感性/3_thermal_mesh_sensitivity.jl` 声明的是“纯热模型（均匀体积热源 + 表面冷却）”，但实际代码仍然在跑 `SPMe`，并且设置了 `opt.per_element_spme = true`、`opt.thermalmodel = "distributed2D"`，还调用了 `setup_thermal2D_mesh`。

同时，脚本里定义了 `compute_q_uniform`，但没有任何调用。也就是说，当前热网格脚本并没有把热源冻结成固定热源场，因此不同 `nθ` 下的结果混杂了：

- SPMe 电化学耦合误差
- 热源分布误差
- 热网格离散误差

这会直接破坏热网格敏感性的可解释性。

### P2：cohesive 特征长度脚本和最终口径不一致

`example/网格敏感性/1_cohesive_characteristic_length.jl` 里用的是

$$
l_c = G_c \cdot E_{eff} / \sigma_{max}^2
$$

其中 `E_eff` 是按正/负极厚度加权得到的等效模量。这个做法和我们后来确认的口径不一致。当前方案要求 `E` 明确取体材料杨氏模量，而不是用一个厚度加权的等效模量去替代。

更关键的是，这个脚本还把热网格候选值直接写成了经验数组 `nθ_thermal = [20, 40, 80, 160]`，没有按热穿透深度来生成热网格候选区间。也就是说，热网格和 cohesive 网格在这里仍然被绑定到同一个脚本逻辑里，和最终方案冲突。

### P3：能量守恒检查中的 fracture energy 计算是失真的

`example/网格敏感性/5_energy_conservation_check.jl` 里，`E_frac` 依赖

- `result["czm damage [0-1]"]`
- `coh_lengths`

但当前代码中：

- `czm damage [0-1]` 并没有作为有效结果键导出，`PostProcessing` 的现有导出里也没有这个键
- `coh_lengths` 是通过 `elem.nodes[1]` 和 `elem.nodes[3]` 计算对角线长度得到的，不是 cohesive 单元的界面长度

因此，即使脚本能运行，`E_frac` 也会被系统性低估或直接保持为零，最终的能量残差 `R` 没有物理可信度。

---

## 3. 重要问题

### P4：CZM 网格脚本声明的指标和实际输出不一致

`example/网格敏感性/4_czm_mesh_sensitivity.jl` 的说明写着要对比 `D_max`、损伤演化速率、断裂单元数，并提到 `traction-separation` 曲线。

但实际脚本只稳定使用了这些输出：

- `czm D_max`
- `czm D_mean`
- `czm n_fractured`
- `czm δ_max_n [m]`

并没有真正把 traction-separation 曲线作为主比较对象，也没有把 load-displacement 曲线落成可复用的指标流程。结果就是“文档上承诺的曲线”与“代码实际产出的指标”不一致。

### P5：电化学网格脚本的热模型口径与注释不一致

`example/网格敏感性/2_electrochemical_mesh_sensitivity.jl` 的注释写的是“SPMe + lumped thermal”，但实现里设置的是 `opt.thermalmodel = "distributed2D"`，还构建了热网格。

如果目标是按照我们前面讨论过的 lumped 口径比较 `max |dT/dt|`，那这里的热模型配置就不对；如果目标是 distributed2D 热场，那脚本注释里 “lumped thermal” 的说法就应该删掉，避免后续误读。

### P6：热网格候选值仍然是硬编码经验值

`example/网格敏感性/3_thermal_mesh_sensitivity.jl` 中的 `THERMAL_Nθ = [20, 40, 80, 160]` 仍是硬编码。即便我们已经确认热网格应按热穿透深度独立定档，这个脚本也还没有把热穿透深度公式真正转成候选生成逻辑。

---

## 4. 建议修正顺序

1. 先把热网格脚本改成真正的固定热源基准，去掉 SPMe 耦合混杂。
2. 再把 cohesive 特征长度脚本拆成两条独立路径：热网格按热穿透深度，CZM 网格按 cohesive 长度。
3. 修改能量守恒检查，改用真实存在且可导出的 CZM 结果键，并使用 cohesive 单元自身长度。
4. 统一 CZM 脚本里的指标命名，把 traction-separation 或 load-displacement 明确成可输出、可对比的 surrogate。
5. 统一电化学脚本的热模型口径，避免 lumped / distributed2D 混写。

---

## 5. 可直接复用的有效输出键

从当前代码来看，下面这些键是已确认存在且可用于评审/后处理的：

- `cell voltage [V]`
- `temperature [K]`
- `thermal2D temperature at nodes [K]`
- `total heat source [W]`
- `czm D_max`
- `czm D_mean`
- `czm n_fractured`
- `czm δ_max_n [m]`
- `czm δ_mean_n [m]`
- `czm traction normal [Pa]`
- `czm traction tangent [Pa]`
- `czm separation normal [m]`
- `czm separation tangent [m]`

如果脚本想继续保留“界面损伤”类指标，就应该优先使用这些已确认导出的键，而不是依赖未导出的占位键。

---

## 6. 结论

这套脚本目前的最大问题不是语法，而是物理口径没有完全收敛：

- 热网格还没有真正独立出来
- cohesive 特征长度和热网格判据仍被混写
- 能量守恒检查的断裂能项不可信
- CZM 指标承诺和实际输出还没有对齐

**结论**: 这批脚本不建议直接作为正式评审结论依据，至少需要先修正热网格、CZM 指标和能量检查三条主线。
