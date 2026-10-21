# v12保留为历史参照

2026-10-21，用户决定将P1独立迁移到3df40fe、移除旧混合候选，并保留v12历史而不重冻结。

此目录README.md、metrics.toml、preflight.log、run.log、source_manifest.tsv五件原字节保持。v12=testexample-20261020T161116+0800，来源为3df40fe加当时未提交的完整混合候选；仅凭HEAD不能复现。README中的“当前基线/后续必须一致”属于冻结当时的历史记录，已由本说明和总入口替代，不适用于独立迁移后的源码。

迁移前工作区（除用户指定原地保留的100GiB输出）已逐SHA备份到KB/raw/repos/jubat/diagnostics/p1_workspace_backup_20261020_181739_3df40fe/，包括51文件冻结来源及完整Git/工作文书。原输出目录D:/Graduate student/实验数据/JuBat_20261020_aa69a5b/output未复制到KB、未移动删除，清单及路径已登记；该备份不是自包含全部运行大数据包。

本次迁移对照以3df40fe+nθ20/60s为共同起点，17结果+3网格+退出码精确21/21，不要求等于v12，未修改/重冻结本目录科学指标。新正式基线尚未建立。
