# 可借鉴的方法与 JuBat 适配

采用独立撰写的操作说明，只借鉴方法；没有复制外部案例的实现或模型设置。
来源读取日期：2026-10-14。后续借鉴新增案例时核对原作者、版本和许可证；不要只依据聚合站点的摘要。

## 外部案例

| 原始来源 | 可借鉴的方法 | 对 JuBat 的适配 |
|---|---|---|
| Anthropic code-simplifier agent | 聚焦近期改动；保持行为；清理重复与过度嵌套；选择明确表达 | 用 AGENTS/科学基线约束 Julia 简化；保留必要错误门、回滚与量纲解释；不照搬其 TypeScript/React 约定 |
| Julia 官方 Performance Tips | 函数边界、类型稳定、预分配、views；先看分配和实测 | 局部 kernel 保持可理解；不全库替换 Dict；views 必须重新检查状态与历史别名；缓存要有所有者 |
| JuBat Simplify 历史案例 | 合并原位热边界入口、去掉无职责兼容转发、清理 silent catch | 以当前调用链重新证实等价性；历史 PASS 不证明现在还可删除 |

## 按收益选择模式

### 合并真正重复的实现

确认数学公式、符号、单位、异常、返回形状和写入顺序均相同。
抽出共享的局部计算或调用同一实现；保留不同 BC、本构、验证算子的独立职责。
不要把 production 和 oracle 收敛到一个实现后，又拿二者“相同”当作独立正确性证据。

### 删除空转发层

检查公开文档、include 链、方法派发与迁出 KB 的测试。
仅删除被证明无额外契约且无兼容责任的内部 wrapper；同时改完真实调用者。
保留选择型号、安装依赖等独立公开入口，即使仓库内调用不多。

### 缩短异常与控制流

把明确非法输入在入口失败；用清晰的早返回处理互斥路径。
统一已证实相同的错误构造；保留错误内容、日志和失败诊断要求。
不要用补零、旧值、截断、首元素、catch 后继续或 `get(...,默认值)` 掩盖必需数据缺失。
已有 guard 的删除仍须符合 AGENTS §9.7，不借“依赖原生报错”跳过其验证条件。

### 局部命名与数据提取

采用能区分当前/候选/已提交、GP/单元、物理/无量纲的名字。
在多个消费者之间明确只读输入与可写输出；必要 copy 是语义的一部分。
只把重复切片提取到局部变量；替换为 view 前验证下游方法接受 SubArray 且不会修改原状态。

### 局部性能改进

先定位热点。提取纯 kernel、减少重复分配或复用线程 workspace 时，保留装配与散射顺序。
避免为了“统一”新增抽象层、泛化注册器、内容哈希缓存或大范围结构体迁移。
不采用 `@fastmath` 或改写归约顺序来达到行数指标。

## 典型任务的正确处理

- “合并热边界的两个重复入口”：先核对副作用与调用者，再做窄批次 A/B 和热边界专项。
- “删除所有校验，把代码变短”：指出有状态与科学契约的校验不能批量删除；给出有证据可移除的局部候选。
- “把 copy 全换成 views”：识别 workspace/提交历史别名；只改生命周期可证明安全的切片。
- “把 k_xy 改正顺便简化”：把公式修正单列为科学变更；修复后以物理验证和获认可的基线开展简化。
- “没 Julia，先整理源码”：做静态审查或可复核补丁，列出待运行项；不要声称完成科学验收。
- “Ring 验证通过，所以生产 FEM 可以删掉另一套入口”：检查两条路径不同的材料与张量逻辑，保留独立验证意义。

## 原始来源链接

- [Anthropic code-simplifier](https://github.com/anthropics/claude-plugins-official/blob/8d392c8169665866be1d95f2bfeda5aaae7149dd/plugins/code-simplifier/agents/code-simplifier.md)，文件 blob `05e361b4ef1b688203251989707f8a924a9ed266`。
- [Julia 1.11 Performance Tips](https://docs.julialang.org/en/v1.11/manual/performance-tips/)，重点读 function barriers、preallocation、views 与 copying。
- [JuBat 历史简化与保留决策](https://github.com/WenStaoCHongP/JuBat/blob/59be8fd27e2a63dbc66e1678ccf4ee59e004202f/Simplify/baseline.md)。
