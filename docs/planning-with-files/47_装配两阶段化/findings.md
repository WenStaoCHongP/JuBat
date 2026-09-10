# 发现与证据

## 三子批实施与门结果（2026-09-09 深夜，全部位级）

- **47-A 参考几何缓存**：`CohesiveMesh.bulk_gp_geom` 惰性字段 + `bulk_gauss_geometry` 访问器（IntQ4 原路径一次性填充 w/detJ/dNdx/dNdy）；`gl_element_residual_tangent` 增 `gp_geom` 缓存分支（`run_gp` 闭包双路驱动）。教训：struct/访问器插入位置在 gl 旧文档串与其函数之间会形成相邻双 docstring → 解析错误（挪到文档串之前解决）。
- **47-B J2 零分配**：弹性支免 `copy(C)`、`_I4` 常量提升、qbar 入参 x 视图——三处同算术位级安全。
- **47-C 两阶段化**：`BulkAssemblyCache`（ke/fe + 与 sparse() 同 pattern 的预分配 CSC + pmap 定位表）挂 `CohesiveMesh.bulk_asm`；阶段一 `Threads.@threads` + `bulk_element_kernel!`（每单元算术与旧循环逐字相同、trial 塑性槽位元素间不相交）；阶段二串行按 (e,a,b) 序累加 = Base sparse! 零初始化+顺序累加 → 逐位一致；split_KG 保留旧三元组路径。**别名陷阱修复**：返回的 CSC 共享 nzval 缓冲会被下一次装配覆写（geometric_stiffness 的 FD 测试当场抓获）——改为返回 `copy(nzval)`（~5 MB/次，相对装配成本可忽略）。

## 终局 A/B（1800 s SOC0.65，195 步与全部科学指标位级不变）

| 配置 | 墙钟 | CZM | 说明 |
|---|---:|---:|---|
| T1-a 后基线 | 293.9 s | 234.8 s | |
| +T1-c `-t1` | **229.7 s** | **171.9 s**（0.881 s/步） | 对最初 1266.6 s 累计 **5.51×** |
| +T1-c `-t8`（OPENBLAS=8） | **173.2 s** | **136.9 s**（0.702 s/步） | 累计 **7.31×**；并行阶段为纯计算，多线程下仍位级 |

60 s 工况：testexample 32.5→24.8 s（-t1）。门链：回归 11/11、testexample v10 与 soc065_60s 逐位一致。

## 结论

三批优化累计（均衡化 → T1-a → T1-c）：单线程 1266.6→229.7 s（5.51×），`-t8` 173.2 s（7.31×），全部位级无损、两份基线无需重冻结。剩余杠杆：修正 Newton（需授权+重冻结）、lu 与散射的进一步并行化、SPMe/热占比（现在 12–19%）。
