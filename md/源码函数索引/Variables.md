# Variables.jl

- **源文件**: `src/Variables.jl`
- **行数**: 278 行
- **函数/struct 计数**: 3 个函数（无独立 struct）
- **职责**: 状态变量字典的构造与更新——`StandardVariables` 全模型变量预分配、`create_element_workspace` 精简单元工作区（节省 ~60% 数组分配）、`Variable_update!` 时间步历史写回（含动态扩容）
- **相关技术文档**: `md/10_参数传递与模块架构.md`、`md/08_逐单元算法.md`

## 数据结构

本文件无独立 struct 定义。所有变量以 `Dict{String, Union{Array{Float64}, Float64}}` 形式管理。

**已知技术债**（参考 `CLAUDE.md`、`MEMORY.md`、`variables-refactor-plan.md`）：
- 字符串键 `"variables[...]"` 在 17 个文件、478 处使用
- `StandardVariables`、`create_element_workspace`、`copy_element_results`（外部）三处硬编码键表
- 当前文件包含其中两处硬编码键表（L16-L50 基础键 + L52-L86 SPMe/P2D 扩展 + L88-L126 热扩展 + L127-L146 CZM 扩展；create_element_workspace 重新枚举一遍）

## 函数清单

### `StandardVariables(case::Case, num::Int64) -> Dict` — L1-L148

按模型类型与选项构造全模型变量字典，预分配 `num` 个时间步的存储。

- **SPM/SPMe 共通基础键**（L16-L50）：33 个键，涵盖粒子浓度、表面浓度、温度、应力、位移、电压、电流、时间等
  - SPM/SPMe 模式下 `Nn = Np = 1`（L4-L8），P2D 模式下取电极节点数
  - 高斯点浓度键（L41-L44）尺寸为 `Nn*opt.Nrn*gsorder` / `Np*opt.Nrp*gsorder`
- **SPMe/P2D 扩展键**（L52-L67）：7 个电解液浓度相关键
- **P2D 专属扩展键**（L69-L86）：16 个电位/高斯点电流/高斯点 OCV 等键
- **热模型键**（L88-L126）：
  - 通用：`temperature`（L88）
  - lumped：`thermal lumped internal heat`（L91）
  - distributed2D：~30 个键（L94-L126），包括分层热源 `q_rxn_ne/pe`、`q_rev_ne/pe`、`q_ohm_s/e_ne/pe`、`q_sp`、`q_pcc/ncc`、单元电压/OCV/SOC/η/dUdT、active_mask、截止诊断等
- **CZM 扩展键**（L127-L146）：基础（`D_max`、`D_mean`、`δ_max_n`、`δ_mean_n`、`n_fractured`，L130-L134）+ 完整 cohesive 场（damage、displacement x/y、traction n/t、separation n/t，L138-L144，需要 `case.czm_mesh !== nothing`）

### `create_element_workspace(case::Case) -> Dict` — L150-L234

创建精简型单元工作区 Dict，仅包含 `SPMe_element` 调用链实际需要的键。排除 distributed2D 专用的 ~30 个 thermal2D 键（由 `CallModel_MultiSPMe` 在单元循环外独立管理）。比 `StandardVariables(case, 1)` 减少约 60% 数组分配。

- **状态提取键**（L167-L170）：4 个粒子浓度键
- **SPMe 专用键**（L173-L185）：6 个电解液浓度键（节点 + 高斯点）
- **标量键**（L187-L192）：`temperature`、`cell voltage`、`time`、`cell current` 以 `Float64` 形式（非数组）
- **SPMe_variables! 结果键**（L193-L200）：8 个交换电流密度/界面电流/过电位/OCV 键
- **高斯点浓度键**（L204-L205）：用 `opt.Nrn/Nrp`（单元数）而非 `mesh.nlen`（节点数）计算高斯点数
- **Mechanicaloutput 结果键**（L208-L219）：10 个应力/位移/耦合系数键，受 `mechanicalmodel == "full"` 门控
- **CZM 键**（L222-L231）：7 个键，受 `czm_enabled` 门控；标量 `D_max/D_mean/n_fractured` 以 `Float64` 形式

### `Variable_update!(variables_hist, variables, v::Int64) -> variables_hist` — L236-L278

将当前时间步的 `variables` 写回历史字典 `variables_hist` 的第 `v` 列。

- **动态扩容**（L238-L250）：当 `v > current_size` 时扩展 `max(1000, current_size ÷ 2)` 列，并 `@warn` 通知
- **类型分发**（L252-L276）：
  - 历史值为 `Array{Float64}` 2D：按行数匹配写入列 `[:, v]`；1D 输入自动转列；标量输入兼容单行
  - 历史值为 `Float64`：直接赋值或从数组首元素取值
  - `hist_keys` 集合过滤：仅更新历史中已存在的键（L252-L254），避免新增键导致 KeyError

## 省略项

无。所有 function 均有独立条目。

### [DEBUG]

| 行号 | 内容 | 用途推测 |
|------|------|----------|
| L244 | `@warn "时间步 $(v) 超过预分配 $(current_size)，扩展 $(expansion_size) 步（变量: $(k)）"` | 动态扩容时的用户告警；结构化但携带动态变量信息，告知用户预分配不足 |

### [PLACEHOLDER]

| 行号 | 内容 | 风险 |
|------|------|------|
| L16 | 硬编码字符串键表（StandardVariables 基础+SPMe+P2D+热+CZM 共 ~80 个键，分散在 L16-L126 的 5 个 if/elseif 分支） | **已知技术债**：键名重复在 `create_element_workspace`（L167-L231）与外部 `copy_element_results` 中再次硬编码，三处不同步风险高（参见 MEMORY.md "variables-refactor-plan.md"）。重构方向：集中到 `VARIABLE_KEYS` 常量或注册表 |
| L160 | `Nn = 1; Np = 1  # SPMe 模式`（create_element_workspace 中硬编码 SPMe 模式） | 注释明示仅支持 SPMe 模式；若用于 P2D 则 Nn/Np 错误。当前调用链仅在 per_element_spme=true 时使用，安全 |
| L243 | `expansion_size = max(1000, current_size ÷ 2)`（动态扩容幅度 magic number） | 1000 与 ÷2 是经验值；过小导致频繁扩容，过大浪费内存。注释无说明 |
| L273 | `isempty(col) \|\| (variables_hist[k] = col[1])`（Float64 历史键的数组输入兜底） | 兜底：空数组时保持原值，非空取首元素。注释缺失，语义略隐晦 |

### [COMPLEX-CHECK]

| 行号 | 内容 | 简化建议 |
|------|------|----------|
| L4 | `if case.opt.model == "SPM" \|\| case.opt.model == "SPMe"; Nn = 1; Np = 1; ...; elseif case.opt.model == "P2D"; Nn = case.mesh["negative electrode"].nlen; ...; end`（嵌套 1 层 + 多字段赋值，跨 L4-L14） | 抽出 `compute_model_dims(case) -> (Nn, Np, Ne_ngs, Ne_pgs)` helper，返回 NamedTuple；同模式在 `create_element_workspace` L158-L162 重复 |
| L52 | `if case.opt.model == "SPMe" \|\| case.opt.model == "P2D"; ...7 个键...; end; if case.opt.model == "P2D"; ...16 个键...; end`（连续两层 model 判定 + 共 ~90 行键赋值，跨 L52-L86） | 抽出 `add_spme_keys!(variables, case, num)` 与 `add_p2d_keys!(variables, case, num)` 独立函数；当前单函数 ~150 行过长 |
| L238 | `for k in keys(variables_hist); if isa(...) && ndims(...) == 2; current_size = ...; if v > current_size; expansion_size = ...; @warn ...; n_rows = ...; new_cols = ...; variables_hist[k] = hcat(...); end; end; end`（嵌套 3 层 + 多类型判定，跨 L238-L250） | 抽出 `ensure_capacity!(hist_arr, v, key_name)` helper 处理单数组扩容；当前逻辑与类型分发写回混合，可读性差 |
| L253 | `for (k, val) in pairs(variables); k in hist_keys \|\| continue; hist_val = variables_hist[k]; if isa(hist_val, Array{Float64}); nrows = ...; if isa(val, Array{Float64}); col = ...; if length(col) == nrows; ...; elseif nrows == 1 && !isempty(col); ...; end; elseif isa(val, Float64) && nrows == 1; ...; end; elseif isa(hist_val, Float64); if isa(val, Float64); ...; elseif isa(val, Array{Float64}); col = ...; isempty(col) \|\| ...; end; end; end`（嵌套 4 层 + 多 `isa` 判定 + 多 `&&`，跨 L253-L276，单条件块 ~100 字符） | 抽出 `write_back!(hist_val, val)` 多方法分派（按 hist_val 类型），消除嵌套；当前 4 层嵌套极难维护 |
