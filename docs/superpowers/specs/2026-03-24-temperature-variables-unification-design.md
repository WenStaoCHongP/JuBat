# 温度变量统一化设计

**日期**: 2026-03-24
**状态**: 修订中
**范围**: 重构温度变量命名，统一内部存储和输出格式

## 1. 问题分析

### 1.1 当前状态（混乱）

**内部 `variables` 字典**:
| 变量名 | 单位 | 用途 |
|--------|------|------|
| `"temperature"` | 无量纲 | 集总温度/状态向量DOF |
| `"T_nodes"` | 无量纲 | 分布式节点温度 |
| `"thermal2D temperature"` | K（混合） | 分布式节点温度 - **冗余且存在双重转换bug** |
| `"T_prev"` | 无量纲 | **未使用** |

**输出 `result` 字典**（不统一）:
| 键名 | 来源 | 问题 |
|------|------|------|
| `"temperature [K]"` | PostProcessing.jl | ✓ 格式正确，已是时间序列 |
| `"thermal2D temperature [K]"` | Solve.jl | 含义不清（节点还是单元？） |
| `"thermal2D T_nodes [K]"` | Solve.jl | ✗ 命名不一致 |
| `"thermal2D T_nodes history [K]"` | Solve.jl | ✗ 与 `"thermal2D temperature [K]"` 重复 |

### 1.2 已知 Bug

**双重单位转换问题**：
- `Solve.jl` 中 `variables["thermal2D temperature"]` 写入时已乘以 `T_ref`
- `PostProcessing.jl` 中又乘以 `T_ref`，导致值被错误放大
- **解决方案**：删除 `"thermal2D temperature"` 变量，统一使用 `"T_nodes"` 存储无量纲值

### 1.3 设计目标

1. **内部存储**：简短名称，存储无量纲值
2. **输出格式**：`"变量名 [单位]"` 格式，存储有量纲值
3. **语义清晰**：区分集总温度、分布式单元温度、分布式节点温度
4. **消除冗余**：删除不必要的重复变量

## 2. 设计方案

### 2.1 内部 `variables` 字典

保留：
| 变量名 | 含义 | 单位 | 用途 |
|--------|------|------|------|
| `"temperature"` | 集总温度 | 无量纲 (T/T_ref) | 集总模型、状态向量DOF |
| `"T_nodes"` | 节点温度场 | 无量纲 (T/T_ref) | 分布式热模型 |

删除：
- ~~`"thermal2D temperature"`~~ - 与 `"T_nodes"` 冗余，且存在双重转换bug
- ~~`"T_prev"`~~ - **确认未使用，删除**

### 2.2 输出 `result` 字典

| 键名 | 数据来源 | 形状 | 说明 |
|------|----------|------|------|
| `"temperature [K]"` | `variables["temperature"] * T_ref` | (1, n_t) | 集总温度时间序列 |
| `"thermal2D temperature [K]"` | 节点平均 → 单元 | (ne, n_t) | **分布式单元温度**时间序列 |
| `"thermal2D temperature at nodes [K]"` | `variables["T_nodes"] * T_ref` | (nT, n_t) | **分布式节点温度**时间序列 |
| `"thermal2D final temperature at nodes [K]"` | 最终状态 `T_nodes_carry * T_ref` | (nT,) | 最终节点温度快照 |

**说明**：
- `"temperature [K]"` 已包含完整时间序列，无需单独的 history 键
- `"thermal2D temperature [K]"` 和 `"thermal2D temperature at nodes [K]"` 同理

### 2.3 键名迁移路径

| 旧键名 | 新键名 | 操作 |
|--------|--------|------|
| `"thermal2D T_nodes [K]"` | `"thermal2D final temperature at nodes [K]"` | 重命名 |
| `"thermal2D T_nodes history [K]"` | 删除 | 与 `"thermal2D temperature [K]"` 重复 |
| `"thermal2D temperature [K]"` | 保持 | 但确保数据来源正确 |

### 2.4 单位转换

所有输出温度通过 `case.param_dim.scale.T_ref` 乘以无量纲值得到开尔文温度。

## 3. 需要修改的文件

| 文件 | 修改内容 |
|------|----------|
| `src/Variables.jl` | 删除 `"thermal2D temperature"` 定义（L99, L127），删除 `"T_prev"`（L98） |
| `src/Solve.jl` | 更新输出键名，移除对 `"thermal2D temperature"` 的写入，确保 `"T_nodes"` 存储无量纲值 |
| `src/PostProcessing.jl` | 统一输出键名格式，修复双重转换问题 |
| `src/CycleData.jl` | 更新 `final_state` 中的温度键名，确保使用 `"T_nodes"` |
| `src/Mechanical.jl` | 确认使用 `"T_nodes"`（当前已正确） |
| `src/Initialisation.jl` | 确认 `"T_nodes"` 返回值一致 |
| `src/CycleSolver.jl` | 更新 `"T_nodes"` 引用，确保兼容新命名 |

## 4. 验证方案

1. 运行现有测试确保功能不变
2. 检查输出 `result` 中的键名符合新规范
3. 验证温度值正确转换（无量纲 → 开尔文）
4. 特别验证不再存在双重转换问题

## 5. 向后兼容

由于这是内部代码重构，暂不添加兼容层。外部调用代码需更新键名引用：

```julia
# 旧代码
T_nodes = result["thermal2D T_nodes [K]"]

# 新代码
T_nodes = result["thermal2D final temperature at nodes [K]"]
```

如有需要，可在后续版本添加兼容性警告。
