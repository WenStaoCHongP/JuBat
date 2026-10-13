# 简化的验证选择

## 先证明基线适用

读取当前 `AGENTS.md`、`Simplify/baseline.md` 及其档案。
创建时的快照登记 v11 / `testexample-20261007T152634+0800`，Julia 1.11.2，快速门脚本 60 s。
档案源码 HEAD 早于本技能所审分支，不可未经检查就当成当前分支必然通过的证明。

先核对实际脚本哈希、配置、依赖、线程、源码变化和用户认可的重冻结记录。若冻结档案不适用，记录不适用原因；仍可做父提交/候选提交同环境 A/B，但不得声称通过冻结基线，也不得自行替换档案。
保留已知错误与基线等价检查的区别：A/B 相同只证明行为保持，不能证明模型物理正确。

## 执行项目要求的门

在完整仓库根目录、已安装匹配 Julia 与依赖的环境执行：

```bash
OPENBLAS_NUM_THREADS=1 GKSwstype=100 julia -t1 --startup-file=no --project=. example/testexample.jl
```

Windows 使用当前基线提供的 PowerShell 配置，不复用旧机器的绝对路径。
运行前确认脚本仍为相应 60 s 配置，避免为了快门临时修改生产示例而未记录。
将 stdout/stderr 和运行身份保存到 `output/testexample/<批次名>/`。

按档案自身 schema 比较退出码、网格/节点/步数、电压、容量、温度、损伤、分离、应力及收敛统计；逐位项逐位比较，打印精度项按记录精度比较。耗时仅参考。
创建时 metrics 中的标准配置是 mix/geo/J2-on、面积反馈 OFF；这一个门不能覆盖 gp_basic、面积反馈 ON、非默认边界或长循环。

按 AGENTS §9.6 恢复受影响 KB 测试：
```bash
julia restore-from-kb.jl --what test
# 只有相关验证确实依赖 tools 时再恢复 tools。
```
KB 默认路径可能只存在于用户本机。恢复失败时报告缺少的验证，使用已存在的受影响验证，不伪造文件或 PASS。
恢复副本的回写和清理按当前 AGENTS 执行；删除前确认是本批创建且已同步的副本。

## 按实际变化选择专项

| 变化范围 | 最少覆盖 |
|---|---|
| 热算子、Q4 梯度或导热张量 | 常温/无源守恒；线性场梯度 patch；同节点编号的网格整体旋转；理想圆环均匀源轴对称性；生产 ThermalDistributed2D 路径 |
| 热边界或源项 | 冷却/产热总量、边界集合；surface/tab/none 的受影响分支；极耳组完整性与局部源注入 |
| 分流与面积反馈 | OFF/ON、部分/全部面积失活、带载/静置、混合电压截止、非法输入失败、加权电流守恒 |
| 状态切片或 workspace | 合法不同 Nrn/Nrp/Nn/Np 配置、SubArray 分派、输入不被修改、结果不被下次 workspace 覆盖 |
| 循环或时间步 | 相位连续性、温度重置、实际 h_used、拒绝/接受区分、截止原因、跨周期损伤/SOH |
| CZM 状态或装配 | gp_basic 及实际受影响方法；成功提交和不收敛回滚；GP 与单元归约；非默认 BC/端点组合 |
| 缓存或并行 | 参数冻结/失效条件、缓存首次与复用、线程 workspace 隔离、串行散射顺序；基线依然 -t1 |
| 导出或后处理 | 完整键/单位/列宽、实际流式输出、必要失败诊断；另跑 couple_example 三张最终 PNG 门 |

利用现有验证；只有缺少覆盖且关系到本次物理或状态风险时才增加针对性测试。
不要以对称矩形上的 Q4 测试代替斜/旋转单元；不要以热算子的对称矩阵、正特征值或 K*1≈0 单独证明方向正确。
不要把 `ThermalDistributed2D_Ring` 通过当成 `ThermalDistributed2D` 的导热张量通过；核对测试真实调用链。
旋转检查需同时旋转边界身份及方向性输入；保持拓扑/节点映射一致，避免重划网格造成不可比。

## 比较性能

只在同输入、同环境、已预热且科学门通过时比较耗时/分配。
先分开 Julia 加载/编译与稳态求解；给出实际重复次数和波动，不用一次耗时宣称提速。
未测量时只说明结构收益，不声称性能提升。

## 审核输出

列出基准 SHA、候选 SHA/补丁、环境、入口/配置、日志、每项 PASS/FAIL/NOT RUN 及原因。
Python 公式复现、静态阅读和技能结构校验均不等于 Julia 生产仿真验证。
没有运行证据时写“验证待补”，不把没有报告失败解释为通过。

## 来源

- [基线索引快照](https://github.com/WenStaoCHongP/JuBat/blob/59be8fd27e2a63dbc66e1678ccf4ee59e004202f/Simplify/baseline.md)
- [结构化指标快照](https://github.com/WenStaoCHongP/JuBat/blob/59be8fd27e2a63dbc66e1678ccf4ee59e004202f/Simplify/baseline/testexample/metrics.toml)
