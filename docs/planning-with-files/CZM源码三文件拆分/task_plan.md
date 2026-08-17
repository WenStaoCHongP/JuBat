# CZM 源码三文件拆分计划

## Goal

结合 `md/源码函数索引/czm.md`，将 `src/czm.jl` 按网格、核心本构和边界条件职责拆为 `CzmMesh.jl`、`Czm.jl`、`CzmBC.jl`，保持公开 API 与数值行为不变。

## Current Phase

Complete

## Phases

### Phase 1：定义、依赖与调用点审计
- [x] 核对源码索引与实际定义
- [x] 分析函数间依赖和加载顺序
- [x] 检查外部调用、include 与文档引用
- **Status:** complete

### Phase 2：实施三文件拆分
- [x] 创建 `CzmMesh.jl`、`Czm.jl`、`CzmBC.jl`
- [x] 更新模块入口加载顺序
- [x] 删除原 `src/czm.jl`
- **Status:** complete

### Phase 3：索引与架构约束更新
- [x] 拆分 `md/源码函数索引/czm.md` 的内容
- [x] 更新 `JuBat.md` 和 planning-with-files 总索引
- [x] 清理当前文档中的旧路径
- **Status:** complete

### Phase 4：验证
- [x] 运行静态定义/依赖检查
- [x] 运行 CZM 聚焦测试
- [x] 运行冻结的 `example/testexample.jl`
- **Status:** complete

## Key Questions

1. 每个 struct/function 应归属网格、核心本构还是边界条件？
2. 新文件之间需要怎样的 include 顺序才能避免前向依赖？
3. 哪些当前文档与测试需要同步新路径？

## Decisions Made

| Decision | Rationale |
|---|---|
| 保持现有函数名和顶层 export 不变 | 本任务只调整文件职责，不改变公共 API |
| 规划文件保存到本任务独立目录 | 遵守项目 planning-with-files 存放和命名约定 |

## Errors Encountered

| Error | Attempt | Resolution |
|---|---:|---|
| 首次将技能读取与无命中的记忆搜索并行执行时，整体脚本因 `rg` 返回 1 而失败 | 1 | 单独完整读取技能及模板；不重复该组合调用 |
| PowerShell 中带转义引号的组合 `rg` 正则被截断为未闭合分组 | 1 | 改用固定字符串搜索和不含引号的独立定义搜索 |
| `apply_patch` 不接受无内容变更的纯 Move hunk | 1 | 使用临时文件名并在移动时加入可识别临时标记，拆分生成时剔除该标记 |
| Windows 的 `Test-Path src/czm.jl` 大小写不敏感，会命中新 `Czm.jl`；组合检查中的 include 正则也返回 1 | 1 | 用目录项 `Name -ceq 'czm.jl'` 精确判断，并用 `Select-String -SimpleMatch` 检查 include |
| PowerShell 调用 Julia `-e` 时剥离了内层路径引号，导致 `src` 被解释为变量 | 1 | 不重复内联命令，改运行仓库现有 CZM 测试文件验证模块加载 |
| 多条 `rg` 顺序执行时，后续无命中返回 1 使组合检查整体标记失败 | 1 | 后续每个可能无命中的搜索独立执行，或在 PowerShell 中捕获输出而不依赖退出码 |
| 总索引补丁包含一个空的 Update hunk，`apply_patch` 校验拒绝 | 1 | 删除空 hunk 后重新应用完整有效补丁 |

## Constraints

- 不改变 CZM 数值公式、参数语义、公开函数名和调用方式。
- 不与当前工作树中的其他用户修改混合或回退。
- 修改后 `testexample` 的冻结数值与 PNG 哈希必须保持一致。
