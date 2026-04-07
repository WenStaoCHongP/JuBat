# src 中温度相关变量梳理

## 说明
- 范围：仅统计 src 下与温度场、温度耦合、热源直接相关的变量定义与关键赋值。
- 归一化：JuBat 内部主要求解量多为无量纲温度（通常为 `T/T_ref`），后处理再转换为 K。
- 本文重点是“变量定义 + 位置 + 物理含义”，不展开全部调用链细节。

## 1. 核心温度状态变量（主状态/主输出）

| 变量名 | 首次/关键定义位置 | 物理含义 | 量纲/单位 | 主要用途 |
|---|---|---|---|---|
| `temperature` | `src/Variables.jl` 第 88 行；赋值见 `src/SPMe.jl` 第 195 行、`src/SPM.jl` 第 83 行、`src/P2D.jl` 第 283 行、`src/CallModel.jl` 第 148 行 | 电化学模型使用的电池特征温度（单温/代表温度） | 内部无量纲；输出时转 K | 进入动力学、传输、热源与电压计算 |
| `T_nodes` | `src/CallModel.jl` 第 47/62/149 行；初始化见 `src/Initialisation.jl` 第 38/105 行；提取函数见 `src/Initialisation.jl` 第 132 行 | 2D 热网格节点温度场 | 内部无量纲（注释中明确） | distributed2D 热自由度主变量，跨时间步传递 |
| `thermal2D temperature at nodes` | `src/Variables.jl` 第 98 行；更新见 `src/Solve.jl` 第 184/237 行 | 用于历史记录与输出的节点温度时序容器 | 内部无量纲；后处理转 K | 温度场结果存储与导出 |
| `thermal2D temperature` | `src/Variables.jl` 第 97 行；生成见 `src/Solve.jl` 第 424 行（由节点均值构造） | 单元尺度温度场（由节点场映射） | 输出时为 K | 单元级温度可视化与统计 |
| `thermal2D final temperature at nodes [K]` | `src/Solve.jl` 第 429 行 | 仿真结束时节点温度场（物理单位） | K | 最终状态输出 |
| `temperature [K]` | `src/PostProcessing.jl` 第 6 行 | 标量温度历史（由 `temperature` 还原） | K | 常规后处理输出 |

## 2. 温度相关初始化与状态布局变量

| 变量名 | 首次/关键定义位置 | 物理含义 | 量纲/单位 | 主要用途 |
|---|---|---|---|---|
| `T0_nodes` | `src/Initialisation.jl` 第 38、105 行 | 初始节点温度向量（全部置 `cell.T0`） | 与内部温度标度一致 | 初始化 distributed2D 热自由度 |
| `layout.thermal_range` | `src/CouplingState.jl` 第 15 行；使用见 `src/Initialisation.jl` 第 133/155 行 | 全局状态向量中热自由度索引区间 | 索引 | 从全局状态中切片/回填温度场 |
| `get_thermal_dofs(...)` | `src/Initialisation.jl` 第 132-133 行 | 从全局状态向量提取温度自由度 | 返回温度向量 | MultiSPMe 下温度场解耦访问 |
| `T_nodes_carry` | `src/Solve.jl` 第 187 行；`src/CycleData.jl` 第 71 行 | 时间推进或相位切换时携带的温度场缓存 | 内部无量纲 | 保持热场连续性，避免每步重建 |

## 3. 温度驱动的电化学耦合变量

| 变量名 | 首次/关键定义位置 | 物理含义 | 量纲/单位 | 主要用途 |
|---|---|---|---|---|
| `T`（局部变量） | 典型见 `src/SPMe.jl` 第 160-164 行、`src/Thermal.jl` 第 4 行、`src/ElectrolyteDiffusion.jl` 第 23 行、`src/ElectrolytePotential.jl` 第 17 行 | 当前参与本构/动力学计算的温度 | 内部无量纲 | Arrhenius、导电率/扩散率、过电势与电压计算 |
| `negative electrode temperature` | `src/Variables.jl` 第 25 行 | 负极温度变量槽位（预留/兼容） | 与内部温标一致 | 变量容器（当前主链多用 `temperature` 与 `T_nodes`） |
| `positive electrode temperature` | `src/Variables.jl` 第 26 行 | 正极温度变量槽位（预留/兼容） | 与内部温标一致 | 变量容器（当前主链多用 `temperature` 与 `T_nodes`） |
| `thermal2D dUdT_n_e` | `src/Variables.jl` 第 116 行；赋值见 `src/CallModel.jl` 第 122 行 | 负极开路电压温度导数（单元级） | OCV 对温度导数（归一化体系下） | 可逆热项与诊断输出 |
| `thermal2D dUdT_p_e` | `src/Variables.jl` 第 117 行；赋值见 `src/CallModel.jl` 第 123 行 | 正极开路电压温度导数（单元级） | 同上 | 可逆热项与诊断输出 |

## 4. 分布式热模型中的热源变量（与温度场直接耦合）

| 变量名 | 首次定义位置 | 关键赋值位置 | 物理含义 | 单位/量纲 |
|---|---|---|---|---|
| `heat_source_fields` | `src/Variables.jl` 第 101 行 | `src/ThermalDistributed.jl` 第 376 行 | 单元总热源场（进入热方程载荷） | 内部无量纲（后续可换算 W/m3） |
| `thermal2D q_rxn_ne` | `src/Variables.jl` 第 102 行 | `src/ThermalDistributed.jl` 第 377 行 | 负极反应热 | 同上 |
| `thermal2D q_rev_ne` | `src/Variables.jl` 第 103 行 | `src/ThermalDistributed.jl` 第 378 行 | 负极可逆热（熵热） | 同上 |
| `thermal2D q_ohm_s_ne` | `src/Variables.jl` 第 104 行 | `src/ThermalDistributed.jl` 第 379 行 | 负极固相欧姆热 | 同上 |
| `thermal2D q_ohm_e_ne` | `src/Variables.jl` 第 105 行 | `src/ThermalDistributed.jl` 第 380 行 | 负极电解液欧姆热 | 同上 |
| `thermal2D q_sp` | `src/Variables.jl` 第 106 行 | `src/ThermalDistributed.jl` 第 381 行 | 隔膜欧姆热 | 同上 |
| `thermal2D q_rxn_pe` | `src/Variables.jl` 第 107 行 | `src/ThermalDistributed.jl` 第 382 行 | 正极反应热 | 同上 |
| `thermal2D q_rev_pe` | `src/Variables.jl` 第 108 行 | `src/ThermalDistributed.jl` 第 383 行 | 正极可逆热（熵热） | 同上 |
| `thermal2D q_ohm_s_pe` | `src/Variables.jl` 第 109 行 | `src/ThermalDistributed.jl` 第 384 行 | 正极固相欧姆热 | 同上 |
| `thermal2D q_ohm_e_pe` | `src/Variables.jl` 第 110 行 | `src/ThermalDistributed.jl` 第 385 行 | 正极电解液欧姆热 | 同上 |
| `thermal2D q_pcc` | `src/Variables.jl` 第 111 行 | `src/ThermalDistributed.jl` 第 386 行 | 正极集流体欧姆热 | 同上 |
| `thermal2D q_ncc` | `src/Variables.jl` 第 112 行 | `src/ThermalDistributed.jl` 第 387 行 | 负极集流体欧姆热 | 同上 |
| `total heat source` | `src/Variables.jl` 第 132 行 | `src/ThermalDistributed.jl` 第 388 行 | 全电芯（面积加权）总热源 | 代码中作为总热功率量输出 |
| `thermal lumped internal heat` | `src/Variables.jl` 第 91 行；赋值见 `src/Thermal.jl` 第 77 行 | 集总热模型内部发热项 | 无量纲（后处理转 W/m3） | lumped 热模型源项 |

## 5. 温度-力学/CZM 耦合变量

| 变量名 | 首次/关键定义位置 | 物理含义 | 单位/量纲 | 主要用途 |
|---|---|---|---|---|
| `dT_elem` | `src/Mechanical.jl` 第 190 行；`src/CycleSolver.jl` 第 581 行 | 单元温升（相对参考温度） | K（或等效温差标度） | 计算热应变、驱动 CZM 热-化学载荷 |
| `thermal2D element thermal strain` | `src/Variables.jl` 第 127 行 | 单元热应变 | 无量纲应变 | 热-扩散耦合力学后处理 |
| `thermal2D element thermal stress` | `src/Variables.jl` 第 123 行 | 单元热应力 | 归一化应力量 | 力学场诊断与损伤分析 |
| `epsilon_0_elem = α_eff*dT + ...` | `src/Mechanical.jl` 第 242 行；`src/czm.jl` 第 363 行 | 热膨胀+浓度膨胀本征应变 | 应变 | 构造热-扩散耦合等效载荷 |

## 6. 输出与数据导出中的温度字段

| 输出字段/结构 | 位置 | 物理含义 |
|---|---|---|
| `result["thermal2D temperature at nodes [K]"]` | `src/Solve.jl` 第 413 行；`src/PostProcessing.jl` 第 63 行 | 节点温度时序（K） |
| `result["thermal2D temperature [K]"]` | `src/Solve.jl` 第 424 行 | 单元温度时序（K） |
| `result["thermal2D final temperature at nodes [K]"]` | `src/Solve.jl` 第 429 行 | 结束时节点温度（K） |
| `TimeStepData.T_nodes` | `src/CycleData.jl` 第 12 行 | 循环导出每个时间步的节点温度（K） |
| `TimeStepData.T_max / T_mean` | `src/CycleData.jl` 第 13-14 行 | 每步最大/平均温度（K） |

## 7. 简要结论
- 主温度状态有两套表达：
  - 标量 `temperature`（用于电化学参数与反应动力学耦合）。
  - 场变量 `T_nodes` / `thermal2D temperature at nodes`（用于 distributed2D 热传导）。
- 分布式热模型通过 `heat_source_fields` 和各分项 `thermal2D q_*` 驱动温度场演化。
- 力学/CZM 通过 `dT_elem` 将温度变化映射为热应变，从而影响损伤与应力演化。
