# 温度变量统一化实现计划

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 重构温度变量命名，统一内部存储（无量纲）和输出格式（有量纲），消除冗余变量和双重转换bug。

**Architecture:** 内部 `variables` 使用简短键名存储无量纲值，输出 `result` 使用 `"变量名 [单位]"` 格式存储有量纲值。删除冗余的 `"thermal2D temperature"` 和 `"T_prev"` 变量。

**Tech Stack:** Julia, Dict 数据结构

---

## 文件结构

| 文件 | 操作 | 说明 |
|------|------|------|
| `src/Variables.jl` | 修改 | 删除 L98 (`T_prev`), L99 和 L127 (`thermal2D temperature`) |
| `src/Solve.jl` | 修改 | 删除 L215, L297, L695 的写入；重构 L450-456 输出 |
| `src/PostProcessing.jl` | 修改 | 修复 L63 双重转换，统一键名 |

## 依赖分析

**`element_nodal_mean` 函数**：定义于 `src/Tools.jl:140`，在模块内部可直接调用。

**`case.param.scale` vs `case.param_dim.scale`**：
- `case.param.scale` - 归一化参数中的缩放因子（无量纲化后的）
- `case.param_dim.scale` - 维度参数中的缩放因子（物理单位）
- 两者应相同，优先使用 `case.param_dim.scale.T_ref`

**CycleData.jl 影响**：使用 `"T_nodes"` 作为内部键，不受此重构影响（仅删除 `"thermal2D temperature"` 和 `"T_prev"`）。

---

## Chunk 1: 删除冗余变量定义

### Task 1: 清理 Variables.jl 中的冗余定义

**Files:**
- Modify: `src/Variables.jl:98-99, 127`

- [ ] **Step 1: 删除 `"T_prev"` 变量定义（L98）**

删除这一行：
```julia
        variables["T_prev"] = zeros(Float64, nT, num)
```

- [ ] **Step 2: 删除第一个 `"thermal2D temperature"` 变量定义（L99）**

删除这一行：
```julia
        variables["thermal2D temperature"] = zeros(Float64, nT, num)
```

- [ ] **Step 3: 删除第二个 `"thermal2D temperature"` 变量定义（L127）**

删除这一行：
```julia
        variables["thermal2D temperature"] = zeros(Float64, nT, num)
```

- [ ] **Step 4: 验证修改**

运行 Julia 检查语法正确：
```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e "include(\"src/JuBat.jl\"); using .JuBat; println(\"✓ Variables.jl 语法正确\")"
```
Expected: `✓ Variables.jl 语法正确`

- [ ] **Step 5: 提交**

```bash
git add src/Variables.jl
git commit -m "refactor(variables): 删除冗余温度变量 T_prev 和 thermal2D temperature

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 2: 移除对冗余变量的写入

### Task 2: 清理 Solve.jl 中的冗余写入

**Files:**
- Modify: `src/Solve.jl:215, 297, 695`

- [ ] **Step 1: 删除 L215 的写入**

找到并删除：
```julia
        variables["thermal2D temperature"] = T_nodes .* case.param_dim.scale.T_ref
```

- [ ] **Step 2: 删除 L297 的写入**

找到并删除：
```julia
                variables["thermal2D temperature"] = T_nodes .* Tref
```

- [ ] **Step 3: 删除 L695 的写入**

找到并删除：
```julia
    variables["thermal2D temperature"] = T_nodes .* case.param_dim.scale.T_ref
```

- [ ] **Step 4: 验证修改**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e "include(\"src/JuBat.jl\"); using .JuBat; println(\"✓ Solve.jl 语法正确\")"
```
Expected: `✓ Solve.jl 语法正确`

- [ ] **Step 5: 提交**

```bash
git add src/Solve.jl
git commit -m "refactor(solve): 移除对 thermal2D temperature 变量的写入

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 3: 重构输出键名

### Task 3: 重构 Solve.jl 输出部分

**Files:**
- Modify: `src/Solve.jl:448-458`

- [ ] **Step 1: 重构 L450-456 的温度输出**

将：
```julia
            result["thermal2D temperature [K]"] = variables_hist["thermal2D temperature"][:, 1:v]
            result["thermal2D T_nodes history [K]"] = variables_hist["thermal2D temperature"][:, 1:v]
        end
        if case.opt.thermal_enabled
            if isa(T_nodes_carry, Array{Float64}) && length(T_nodes_carry) == case.mesh["thermal2D"].nlen
                Tref = case.param_dim.scale.T_ref
                result["thermal2D T_nodes [K]"] = T_nodes_carry .* Tref
```

替换为：
```julia
            # 节点温度时间序列（无量纲 → 有量纲）
            Tref = case.param_dim.scale.T_ref
            T_nodes_hist = variables_hist["T_nodes"][:, 1:v] .* Tref
            result["thermal2D temperature at nodes [K]"] = T_nodes_hist

            # 单元温度时间序列（节点平均 → 单元）
            mesh_th = case.mesh["thermal2D"]
            ne = size(mesh_th.element, 1)
            n_t = size(T_nodes_hist, 2)
            T_elem_hist = zeros(Float64, ne, n_t)
            for ti in 1:n_t
                T_nodes_t = variables_hist["T_nodes"][:, ti]
                T_elem_hist[:, ti] = element_nodal_mean(mesh_th, T_nodes_t)
            end
            result["thermal2D temperature [K]"] = T_elem_hist .* Tref
        end
        if case.opt.thermal_enabled
            if isa(T_nodes_carry, Array{Float64}) && length(T_nodes_carry) == case.mesh["thermal2D"].nlen
                Tref = case.param_dim.scale.T_ref
                result["thermal2D final temperature at nodes [K]"] = T_nodes_carry .* Tref
```

- [ ] **Step 2: 验证修改**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e "include(\"src/JuBat.jl\"); using .JuBat; println(\"✓ Solve.jl 输出重构完成\")"
```
Expected: `✓ Solve.jl 输出重构完成`

- [ ] **Step 3: 提交**

```bash
git add src/Solve.jl
git commit -m "refactor(solve): 重构温度输出键名，统一命名规范

- thermal2D T_nodes [K] → thermal2D final temperature at nodes [K]
- thermal2D T_nodes history [K] → 删除（冗余）
- 添加 thermal2D temperature at nodes [K]（节点温度历史）
- thermal2D temperature [K] 改为单元温度（节点平均）

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 4: 修复 PostProcessing.jl

### Task 4: 修复 PostProcessing.jl 中的双重转换和键名

**Files:**
- Modify: `src/PostProcessing.jl:63`

- [ ] **Step 1: 修复 L63 的双重转换问题**

将：
```julia
        result["thermal2D temperature [K]"] = variables["thermal2D temperature"][:, 1:v] * case.param.scale.T_ref
```

替换为：
```julia
        # 节点温度时间序列（从无量纲 T_nodes 转换）
        result["thermal2D temperature at nodes [K]"] = variables["T_nodes"][:, 1:v] * case.param_dim.scale.T_ref
```

- [ ] **Step 2: 验证修改**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia -e "include(\"src/JuBat.jl\"); using .JuBat; println(\"✓ PostProcessing.jl 修复完成\")"
```
Expected: `✓ PostProcessing.jl 修复完成`

- [ ] **Step 3: 提交**

```bash
git add src/PostProcessing.jl
git commit -m "fix(postprocessing): 修复温度双重转换bug，统一键名

- 使用 T_nodes（无量纲）替代 thermal2D temperature（有量纲）
- 修复双重 T_ref 乘法问题
- 统一键名为 thermal2D temperature at nodes [K]

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Chunk 5: 验证和最终提交

### Task 5: 验证重构完整性

**Files:**
- None (验证步骤)

- [ ] **Step 1: 检查是否还有对旧变量的引用**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && grep -rn 'variables\["thermal2D temperature"\]' src/ && grep -rn 'variables_hist\["thermal2D temperature"\]' src/ && grep -rn 'variables\["T_prev"\]' src/
```
Expected: 无输出（表示已全部清理）

- [ ] **Step 2: 验证新键名已正确添加**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && grep -rn 'thermal2D temperature at nodes' src/ && grep -rn 'thermal2D final temperature at nodes' src/
```
Expected: 输出显示新键名存在于 Solve.jl 和 PostProcessing.jl

- [ ] **Step 3: 运行示例测试**

```bash
cd "D:/OneDrive/Desktop/Jubat For Cursor/JuBat" && julia example/minimal_example.jl
```
Expected: 成功运行无报错

- [ ] **Step 4: 最终提交**

```bash
git add -A
git commit -m "refactor(temperature): 完成温度变量统一化重构

## 变更摘要
- 删除冗余变量: T_prev, thermal2D temperature
- 统一内部存储: temperature (集总), T_nodes (分布式节点)
- 统一输出格式: 变量名 [单位]
- 修复双重单位转换bug

## 键名迁移
- thermal2D T_nodes [K] → thermal2D final temperature at nodes [K]
- thermal2D temperature [K] → 单元温度（节点平均）
- 新增: thermal2D temperature at nodes [K]

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## 输出键名对照表

| 旧键名 | 新键名 | 说明 |
|--------|--------|------|
| `"temperature [K]"` | 保持不变 | 集总温度 |
| `"thermal2D temperature [K]"` | 保持键名，改为单元温度 | 节点平均到单元 |
| `"thermal2D T_nodes [K]"` | `"thermal2D final temperature at nodes [K]"` | 最终节点温度快照 |
| `"thermal2D T_nodes history [K]"` | 删除 | 冗余 |
| - | `"thermal2D temperature at nodes [K]"` | 新增：节点温度历史 |
