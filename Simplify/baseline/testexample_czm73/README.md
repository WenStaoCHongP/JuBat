# 73真实GP历史修复后的60秒文字基线

冻结：2026-10-22；用户明确授权审核正确的本构/历史定义修正后重冻结。当前入口为 `example/testexample.jl`；父HEAD `8ac0ae512e85131ecb5d86517148e0bf62c8405c` 加73三源码候选（CouplingState、CzmSolve、CzmPostProcess）。主仓尚未提交73。本次只复用成功运行，不运行新仿真。

## 适用工况与命令

Julia1.11.2，`-t1 --startup-file=no --project=.`，`OPENBLAS_NUM_THREADS=1`、`GKSwstype=100`；Jellyroll、5A、nθ20、60s、多SPMe/二维热，mix/basic，load_steps=10、tol=1e-3、geo_nonlinear=true、j2_plasticity=true，dt=[0.5,10]、Crank–Nicolson。完整入口定义由[保存入口原件](source_snapshot/example/testexample.jl)固定，其他默认选项由73源码清单固定。实际命令：

```powershell
& 'C:/Users/user/AppData/Local/Programs/Julia-1.11.2/bin/julia.exe' -t1 --startup-file=no --project=. example/testexample.jl
```

[run.log](run.log)与[退出码](exit_code.txt)保持已完成候选运行原字节。冻结[metrics.toml](metrics.toml)中的17科学结果、3网格和退出码共21项，按打印精度字符串精确比较，包括环向/切向应力范围、von Mises和等效塑性应变；耗时不计科学指标。

## 已审核差异，不能写成原对照21/21

[parent_comparison.json](parent_comparison.json)保留原始20/21比较（其中refrozen=false属于比较当时状态），父运行[parent_run.log](parent_run.log)保留。唯一不同是最大有效分离1.7543e-12→1.7886e-12m。原始对照仍为20/21，不因本次冻结改写。

root与独立审查确认：旧basic使用单元平均开度历史；73使用真实GP几何历史峰值再取最大值，max(norm(mean gap))与max(norm(GP gap))一般不同。聚合取GP delta_max_eff最大值；Solve在输出准备后用候选历史覆盖当前平均分离输出，非力学步从接受史读取，PostProcessing恢复单位，入口输出全时域历史峰值。生产反例覆盖非均匀GP及符号相消。18条CZM调试序列相同；其余20项按原打印精度相同。未据此声称相同接受位移或全部未打印量全精度不变。

## 源码和工件身份

[source_manifest.tsv](source_manifest.tsv)共74项：73项复用半循环最终运行前的source_snapshot原件，逐项与source_before/source_after及当前字节核对；原件保存路径逐行明确。60s入口不在原73项manifest中，故本次另存testexample.jl并证明与8ac0ae5不可变Git blob逐字节相同；这是运行后基于未变HEAD的身份核验，不冒称运行前捕获。三源码候选包含在[candidate.diff](candidate.diff)，该原始差异还含73测试变化；不将测试变化称为科学源码变化。

[source_before.json](source_before.json)、[source_after.json](source_after.json)保留来源身份原件；[artifact_manifest.csv](artifact_manifest.csv)锚定本目录工件与所引用73份源码原件。历史v11/v12/testexample及独立半循环基线不覆盖、不更改冻结数值。原实施报告中“不重冻结”等文字是用户本次追加授权前的历史状态。

该基线只用于本入口配置的后续文字回归，不等于J2完整工况、74 Newton、受损长期循环、生产精度或全面物理验收；60s D=0。半循环受损回归另行保持其独立冻结值。
