# 发现与证据

## R1-d：load_substep geo 求解均衡化（2026-09-11，全门通过）

**改动**：`newton_raphson_czm` 求解行 `Δu = K_bc \ R_bc` → `geo_nl ? solve_equilibrated(K_bc, R_bc) : (K_bc \ R_bc)`（src/CzmSolve.jl 单点；失败语义不变：奇异即抛）。动机：均衡化（b68873e）只进了 basic，load_substep+geo 仍踩 69× 填充慢路径——四求解器重复的真实代价实例。

**改动前后探针**（nθ=24 玩具网格，geo+J2，dT=0.8/Δsoc=0.032）：
- 墙钟 **32.7→15.0 s（−54%）**——慢路径坐实且被消灭
- 迭代轨迹完全相同（46 迭代、R=6.6e-8）——均衡化只改舍入不改路径
- 位移校验和漂移 4.4e-10 相对——恰为均衡化已知偏差量级

**门**：回归 12/12（含新增 test_czm_load_substep.jl——load_substep 此前零测试覆盖，含与 basic 的生产量级载荷交叉校验 rtol 1e-8）；testexample v10 与 soc065_60s **逐位零变化**（basic 路径未触及）。

**测试设计教训**：交叉校验载荷必须取生产步增量量级（dT≈0.05/Δsoc≈0.003）——更大会使 basic 单步全量跳跃的中间态触发 J2 局部 Newton 发散（ Toys 网格上限约 dT=0.1）。

**过程乌龙两则**：①回归循环忘先 restore test/（11 连"败"实为文件不存在）；②restore 的 cp(force) 整目录替换会删掉工作区新增测试——新增测试须在 restore 之后重建或先 --sync 再 restore。
