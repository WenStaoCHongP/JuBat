# Option.jl

- **源文件**：`src/Option.jl`
- **行数**：94 行
- **函数/struct 计数**：3 个 struct（`CycleOption`、`CzmOptions`、`Option`）+ 1 个 enum（`PhaseType`）
- **职责**：定义单次求解、CZM 子配置与循环求解配置。
- **相关技术文档**：`md/10_参数传递与模块架构.md`、`md/09_分流求解器.md`、`md/06_内聚力模型_CZM.md`

## 数据结构

### `CycleOption` — L5-L25

循环数、四阶段时长、充放电电流、上下截止电压、初始 SOC、温度重置和循环时间步配置。

### `PhaseType` — L30-L34

`PHASE_CHARGE`、`PHASE_REST`、`PHASE_DISCHARGE`。

### `CzmOptions` — L35-L56

20 个 CZM 求解/耦合选项的唯一宿主：

| 字段 | 默认值 | 说明 |
|---|---:|---|
| `enabled` | false | 启用 CZM |
| `model` | `"model1"` | `model1` / `mix` |
| `update_interval` | 1 | 损伤更新时间步间隔 |
| `soh_threshold` | 0.8 | 循环 SOH 终止阈值 |
| `inner_exit_only` | true | 断裂时仅内圈退出电化学 |
| `fix_inner` | true | true 时内外圈均固定 |
| `iter_method` | `"basic"` | basic / load_substep / arc_length |
| `max_iter` | 100 | Newton 最大迭代数 |
| `tol` | 1e-4 | 收敛容差 |
| `load_steps` | 2 | 载荷子步数 |
| `arc_length_alpha` | 1.0 | 弧长系数 |
| `viscous_enabled` | false | 粘性正则化开关 |
| `viscous_tau` | 0.0 | 物理松弛时间 [s] |
| `area_loss_enabled` | false | 渐进有效面积损失 |
| `area_loss_threshold` | 0.83 | 面积缩减阈值 |
| `geo_nonlinear` | false | GL/TL 残差与初应力刚度 |
| `winding_prestress` | false | 卷绕预应力 |
| `j2_plasticity` | false | PCC/NCC J2 塑性 |
| `continuous_feedback` | false | 连续损伤–电–热反馈 |
| `friction_mu` | 0.10 | SP 摩擦预留值；当前无消费者 |

模型选择和粘性时间只从 `opt.czm.model/viscous_tau` 读取，不复制进材料参数。

### `Option` — L58-L94

主配置包含：

- 电化学/空间离散：`Np/Ns/Nn/Nrp/Nrn/model/meshType/gsorder/dimension`
- 时间推进与公共接口：`time/Current/coupleMethod/coupleOrder/y0/dt/dtType/dtThreshold/solveType/outputType/jacobi`
- 多物理场选择：`thermalmodel/mechanicalmodel/cite`
- 热模块：`thermal_enabled/thermal_dim/thermalmeshType/cool_method/collector_seeded/per_element_spme`
- CZM：`czm::CzmOptions = CzmOptions()`
- 调试：`debug_coupling/debug_log_path`

`mechanicalmodel` 保留顶层；所有原平铺 `czm_*` 字段已删除。公共 `opt.y0` 保留。

## 函数清单

本文件无独立函数。
