# SOC 0.65 / 60 s 优化基线

2026-09-08 用户要求把1800 s工况改为60 s，并作为后续优化基线。本档案独立于原 `testexample` v10，不替换或撤销原门禁。

随后用户要求恢复1800 s以分析实际塑性阶段，当前 `example/testexample_soc065_1800s.jl` 已恢复1800 s。本60 s历史基线保留于 `script_snapshot.jl`（SHA-256与metrics中脚本哈希相同）；通过下面的include_string命令保留原脚本的相对路径语义。SOC=0.65、5 A、表面冷却、mix、GL/TL、J2、无预应力、fix_inner=false、原分层端点约束保持不变。

## 冻结依据

- 单次完整60 s执行，退出码0；21个时间记录点（含初始点），20/20次CZM更新收敛。
- `metrics.toml` 保存按打印精度冻结的科学指标；`run.log` 是本次完整控制台输出。
- `source_manifest.tsv` 记录当前dirty工作区源码及脚本的原始字节SHA-256，不是旧基线的LF归一化哈希。
- 当前60 s未屈服且未损伤；优化塑性/软化路径仍必须运行主动屈服、卸载、trial/commit、失败回滚的定向测试。

## 运行与判定

```powershell
$env:GKSwstype = '100'
$env:JULIA_NUM_THREADS = '1'
& 'D:\Julia-1.11.2\bin\julia.exe' --startup-file=no --threads=1 --project=. -e 'using LinearAlgebra; BLAS.set_num_threads(14); include_string(Main, read("Simplify/baseline/testexample_soc065_60s/script_snapshot.jl", String), abspath("example/testexample_soc065_1800s.jl"))'
```

比较退出码、配置、网格、记录点/更新数和全部科学结果：电压4位小数，容量4位小数，温度2位小数，分离/应力/塑性按脚本科学计数法4位小数；必须一致，不能以“接近”放行。源码哈希用于溯源，不要求优化后哈希不变。

本轮仅在内存中给 `JuBat.Solve(case)` 添加Profile宏，未修改src。76.204 s包含采样和首次编译，**仅是时间参考，不是无采样稳态性能门**。实施优化前应按相同60 s工况建立无采样性能对照，分别报告冷启动和预热后耗时，固定BLAS线程配置。

原始采样位于 `output/testexample_soc065_1800s/profile_60s_solve_20260908/`；主要结论另存 `profile_groups.txt`。嵌套采样占比不可相加。中止的1800 s与未进入Solve的采样启动均不属于本基线。
