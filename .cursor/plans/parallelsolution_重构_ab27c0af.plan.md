---
name: Parallelsolution 重构
overview: 将 src/Parallelsolution.jl 从单文件重构为模块化结构，消除代码冗余，统一命名规范。
todos:
  - id: rename_functions
    content: 重命名函数： 调试函数加 debug 前缀， 其他函数去掉下划线前缀
    status: completed
  - id: eliminate_redundancy
    content: "消除冗余逻辑: V_branches 重复计算, 归一化逻辑分散, inactive_reason 计算, 结果写入分散"
    status: completed
  - id: verify_compatibility
    content: 验证 API 兼容性和数值结果一致性
    status: completed
  - id: update_exports
    content: 确认 JuBat.jl 中 export 语句是否需要修改
    status: completed
isProject: false
---

## 1. 重命名函数

- 调试函数: `_debug_xxx` → `debug_xxx`
- 其他函数: `_xxx` → `xxx` (去掉下划线前缀)
- 例: `_branch_voltage` → `branch_voltage`
- 例: `_newton_iteration!` → `newton_iteration!`
- 例: `_line_search` → `line_search`

1. 消除冗余逻辑
  - V_branches 重复计算 (L521-532, L535) → 提取为 `compute_initial_voltage()` 复用
  - 归一化逻辑分散 (`_initialize_currents` L179-191, L551-574) → 统一使用 `normalize_currents!()`
  - inactive_reason 计算 (L594-602) → 提取为 `compute_inactive_reason()`
2. 保持 API 兼容
  - `solve_branch_currents_newton` 签名不变
  - `JuBat.jl` 中的 export 语句不变
3. 保持单文件结构
  - 不创建 Parallelsolution/ 文件夹
  - 便于维护和 易于追踪

## 2. 函数命名规范


| 原函数名                                  | 新函数名                           | 状态    |
| ------------------------------------- | ------------------------------ | ----- |
| `_debug_check_prefactors`             | `debug_check_prefactors`       | 调试函数  |
| `_debug_check_coefficients`           | `debug_check_coefficients`     | 调试函数  |
| `_debug_check_initial_voltage`        | `debug_check_voltage`          | 调试函数  |
| `_scalarize`                          | `scalarize`                    | 辅助函数  |
| `_compute_electrochemical_prefactors` | `compute_prefactors`           | 核心计算  |
| `_compute_element_coefficients`       | `compute_element_coefficients` | 核心计算  |
| `_compute_all_coefficients`           | `compute_all_coefficients`     | 核心计算  |
| `_branch_voltage`                     | `branch_voltage`               | 模型函数  |
| `_branch_dVdI`                        | `branch_dVdI`                  | 模型函数  |
| `_newton_iteration!`                  | `newton_iteration!`            | 求解器   |
| `_line_search`                        | `line_search`                  | 求解器辅助 |
| `_initialize_currents`                | `initialize_currents`          | 初始化   |
| `_check_voltage_bounds`               | `check_voltage_bounds`         | 边界检查  |
| `_detect_cutoff_elements`             | `detect_cutoff_elements`       | 截止检测  |


## 3. 消除冗余逻辑


| 冗余位置                             | 儿理方式                            |
| -------------------------------- | ------------------------------- |
| V_branches 重复计算 (L521-532, L535) | 提取为 `compute_initial_voltage()` |
| 归一化逻辑分散 (L179-191, L551-574)     | 统一使用 `normalize_currents!()`    |
| inactive_reason 计算 (L594-602)    | 提取为 `compute_inactive_reason()` |


| active_mask 合并分散 (L494-502) | 在主函数中内联处理 |
| 初始电压计算 (L521-532) | 在 `initialize_currents` 中完成 |
| 结果写入分散 (L579-617) | 提取为 `write_results!()` |

## 4. 测试验证

1. 运行现有示例脚本确认功能正常
2. 对比重构前后的数值结果
3. 检查调试模式输出

