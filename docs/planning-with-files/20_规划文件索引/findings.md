# planning-with-files 索引调研记录

## Requirements
- 在 `docs/planning-with-files/` 下添加索引 Markdown 文件。
- 标注每个任务/文件的大概执行内容、建立时间、修改次数等元数据。
- 命名和组织方式与现有目录一致。

## Metadata Rules
- 建立时间：已跟踪文件取 Git 最早提交日期；未跟踪文件取文件系统 `CreationTime`。
- 最近修改：已跟踪文件取 Git 最近提交日期；未跟踪文件取文件系统 `LastWriteTime`。
- 修改次数：统计包含该路径的 Git 提交数；未跟踪文件记为 `0（未跟踪）`。
- 任务说明：优先读取 `task_plan.md` 的 Goal，其次读取 Markdown 一级标题，探针脚本按文件名和首部注释概括。
- 状态：优先读取 `task_plan.md` 的 Current Phase/Status；缺失时标为“未声明”。

## Findings
- 现有目录既包含标准三文件任务，也包含探针脚本和独立设计文档，因此索引必须逐文件覆盖，不能只列 `task_plan/findings/progress`。
- 工具调用返回可在同一执行单元内解析，因此可由只读 PowerShell 生成 Markdown，再通过 `apply_patch` 写入索引，避免使用 shell 重定向写文件。
- 当前共盘点 17 个任务目录、56 个既有文件；加入总索引后为 57 个文件，其中 48 个已跟踪、9 个未跟踪或新增。
- 索引已生成目录汇总、根目录文件和 17 个逐目录明细章节，共 236 行。
- 部分旧任务没有 `Current Phase` 或标准 `Goal`，索引按“未声明”或首个文档标题回退，避免编造完成状态。
