# 温度变量统一化设计

**日期**: 2026-03-24
**状态**: 已批准
**范围**: 重构温度变量命名，统一内部存储和输出格式

## 1. 问题分析

### 1.1 当前状态（混乱）

**内部 `variables` 字典**:
| 变量名 | 单位 | 用途 |
|--------|------|------|
| `"temperature"` | 无量纲 | 集总温度/状态向量DOF |
| `"T_nodes"` | 无量纲 | 分布式节点温度 |
| `"thermal2D temperature"` | K | 分布式节点温度（有量纲）- 冗余 |
| `"T_prev"` | 无量纲 | 前一步温度 - 可能不需要 |

**输出 `result` 字典**（不统一）:
| 键名 | 来源 | 问题 |
|------|------|------|
| `"temperature [K]"` | PostProcessing.jl | ✓ 格式正确 |
| `"thermal2D temperature [K]"` | Solve.jl | ✓ 但含义不清（节点还是单元？） |
| `"thermal2D T_nodes [K]"` | Solve.jl | ✗ 命名不一致 |
| `"thermal2D T_nodes history [K]"` | Solve.jl | ✗ 命名不一致 |

### 1.2 设计目标

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
- ~~`"thermal2D temperature"`~~ - 与 `"T_nodes"` 冗余
- ~~`"T_prev"`~~ - 如确认不需要则删除

### 2.2 输出 `result` 字典

| 键名 | 数据来源 | 单位 | 说明 |
|------|----------|------|------|
| `"temperature [K]"` | `variables["temperature"]` | K | 集总温度 |
| `"temperature history [K]"` | 历史记录 | K | 集总温度时间序列 |
| `"thermal2D temperature [K]"` | 节点平均 → 单元 | K | 分布式单元温度（标量或向量） |
| `"thermal2D temperature history [K]"` | 历史记录 | K | 分布式单元温度时间序列 |
| `"thermal2D temperature at nodes [K]"` | `variables["T_nodes"]` | K | 分布式节点温度（向量） |
| `"thermal2D temperature at nodes history [K]"` | 历史记录 | K | 分布式节点温度时间序列 |

### 2.3 单位转换

所有输出温度通过 `case.param_dim.scale.T_ref` 乘以无量纲值得到开尔文温度。

## 3. 需要修改的文件

| 文件 | 修改内容 |
|------|----------|
| `src/Variables.jl` | 删除 `"thermal2D temperature"` 定义（L99, L127），删除 `"T_prev"` |
| `src/Solve.jl` | 更新输出键名，移除对 `"thermal2D temperature"` 的写入 |
| `src/PostProcessing.jl` | 统一输出键名格式 |
| `src/CycleData.jl` | 更新 `final_state` 中的温度键名 |
| `src/Mechanical.jl` | 确保使用 `"T_nodes"` |
| `src/ThermalDistributed.jl` | 如有引用则更新 |

## 4. 验证方案

1. 运行现有测试确保功能不变
2. 检查输出 `result` 中的键名符合新规范
3. 验证温度值正确转换（无量纲 → 开尔文）

## 5. 向后兼容

如果外部代码依赖旧的键名（如 `"thermal2D T_nodes [K]"`），需要通知更新。
