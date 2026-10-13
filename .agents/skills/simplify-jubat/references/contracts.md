# JuBat 的简化边界

依据：2026-10-14 读取的 `WenStaoCHongP/JuBat`，分支 `codex/src-physics-modularization`，提交 `59be8fd27e2a63dbc66e1678ccf4ee59e004202f`。
每次执行重读当前 AGENTS 与源码；以下路径不是冻结接口或完整 API 清单。

## 数据流与责任

| 链路 | 当前定位 | 简化必须保持的责任 |
|---|---|---|
| 模型调度与局部模型 | CallModel.jl、SPMe.jl | 逐单元化学状态，显式单元电流/温度输入，workspace 与返回变量的别名安全 |
| 时间推进 | Solve.jl、CycleSolver.jl、CycleData.jl | 外部状态长度、相位继承、实际步长、拒绝/接受时序、阶段截止 |
| 分流 | Parallelsolution.jl | 无量纲电流、有效面积权重、活动集合、各截止原因 |
| 热 FEM | ThermalDistributed.jl、Materialmatrix.jl、SetMesh.jl | 梯度映射、局部导热方向、积分权重、边界和源项符号 |
| 卷绕拓扑 | Jellyrollmodel.jl、CzmMesh.jl | 热/力共享角节点、Φ 配对、合并节点、父热单元映射 |
| 机械与界面 | CzmSolve.jl、czm.jl、CzmPlasticity.jl、Mechanical.jl、CzmBC.jl | 收敛提交、失败回滚、GP/单元口径、边界配置 |
| 导出 | Variables.jl、PostProcessing.jl、CsvExport.jl | 完整键名、单位、列形状、历史与末状态语义 |

## 必须保留的契约

- 量纲：内部按参考尺度归一化。时间按 `scale.t0`；电流按 `scale.I_typ`；二维积分面积按 `scale.L^2` 还原；结果变量按各自尺度还原。不要把字段名中“current”都解释成实际支路安培数。
- 分流：当前有效面积为 `A_e*(1-D_e)`，权重为有效面积归一化。约束是 `sum(w_e*I_e)=I_total`；物理总量贡献为 `w_e*I_e*I_typ`。不要改成 `sum(I_e)=I_total`。
- 面积失活和电压截止是不同原因；不要统一成一个无法追溯的布尔值。全面积失活带载失败，静置走明确诊断路径；事件字段须覆盖所有合法返回路径。
- 静置严格零支路是当前项目选择；不要因并联系统的一般物理直觉重新引入正负抵消电流。
- 粒子平均：节点算术平均与 r²/FE 体积平均不是等价简化。核实生产消费者实际用哪个键，禁止悄悄切换。径向粒子网格与电极厚度网格不同。
- 状态：`case.mech::MechState` 保存已提交位移、损伤、塑性等；材料参数归一化后冻结；mesh 保存拓扑及装配缓存。相位电压界限的临时设置与 finally 恢复是现有明确流程，不要与材料缓存参数混同。不要把演化状态迁回网格或通过缓存复用写坏历史。
- 圈边界和端点：同时核对 `outer_bc/fix_inner/fix_start/fix_end/endpoint_variant`；保持 CZM 与 bonded 域各自的节点身份和配置语义。不要用默认门覆盖合法非默认组合。
- 应变与材料：颗粒模量和极片模量属于不同尺度；热/化学本征应变按实际材料层消费。不要恢复跨层统一有效模量或统一电极膨胀。
- 数值历史：候选计数、接受计数、输出列号不同；观察点应使用实际 `h_used`，不能拿自适应后的下一步长替代。
- 接触：Φ 节点配对、完美粘结与真正 SP—涂层单边接触/摩擦不同；不要在简化中声称新增物理机制。
- 公开 API：保留默认值、方法派发、关键字、完整 result 键名、单位、返回类型及异常。修改公开接口属于另一个任务。

## 热模块特别检查

按实际矩阵布局推导，勿机械套用转置：

```julia
# 若局部梯度为 4×2，J 为 2×2：
J = transpose(dN_dxi) * x_e
dN_dx = dN_dxi * transpose(inv(J))
# 可讨论改成线性求解，但会改变浮点路径；不属于自动等价删除。
```

径/切向单位向量为 `e_r=(cosθ,sinθ)`、`e_t=(-sinθ,cosθ)` 时，
`k=λr*e_r*e_r' + λt*e_t*e_t'`，故 `k_xy=(λr-λt)*sinθ*cosθ`。
整体旋转同一网格及物理输入时，按对应节点编号比较离散算子应保持一致；真实螺旋接缝、极耳或非均匀输入并不要求同半径温度完全相同。

在上述创建时快照中，ShapeFunction2D 用了 J⁻¹ 而非 J⁻ᵀ，热张量 k_xy 使用了相反符号；Mechanical 的 bonded 装配另有局部重算梯度。源码审查和独立 Python 公式复现支持这两处问题，尚未用 Julia 全耦合验证。后续使用先确认是否已修复；不要将改正公式归类为“结果不变”的简化。

## 来源

- [项目规则快照](https://github.com/WenStaoCHongP/JuBat/blob/59be8fd27e2a63dbc66e1678ccf4ee59e004202f/AGENTS.md)
- [当前审查版本源码](https://github.com/WenStaoCHongP/JuBat/tree/59be8fd27e2a63dbc66e1678ccf4ee59e004202f/src)
