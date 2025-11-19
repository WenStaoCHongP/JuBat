# 调试信息说明

## 问题描述
在运行 `testexample.jl` 时出现电压变成 NaN 的错误：
```
thermal2D common voltage out of bounds: V(nd)=NaN, V(V)=NaN
```

## 已添加的调试信息

在 `src/SPMe.jl` 的 `solve_branch_currents_newton` 函数中添加了以下调试打印：

### 1. 预计算值检查（第 335-349 行）
检查关键的预计算参数是否包含 NaN/Inf：
- `prefactor_n`, `prefactor_p` - 交换电流密度的预因子
- `csn_av`, `csp_av` - 平均浓度
- `u_n_ref_val`, `u_p_ref_val` - 参考电位
- `c_sigma` - 电阻系数

如果这些值有 NaN，会打印：
```
⚠ [DEBUG] Prefactor values contain NaN/Inf:
```
并显示所有相关值和输入浓度数据。

### 2. 系数计算检查（第 362-375 行）
对每个热单元，检查计算的电化学系数是否有 NaN：
- `C1` - 开路电压相关项
- `C2` - 温度系数
- `alpha_p`, `alpha_n` - 过电位系数
- `C5` - 总电阻

如果检测到 NaN，会打印：
```
⚠ [DEBUG] Element e coeffs contains NaN/Inf:
```
并显示该单元的所有中间计算值。

### 3. 初始电压计算检查（第 400-414 行）
检查初始电压猜测值：
```
⚠ [DEBUG] Initial voltage V is NaN/Inf:
```
会显示哪个单元的分支电压是 NaN，以及对应的系数值。

### 4. 牛顿迭代过程检查（第 492-497 行）
在迭代过程中，如果 V_trial 变成 NaN：
```
⚠ [DEBUG] V_trial is NaN/Inf in iteration iter:
```
显示当前电压、增量、残差等信息。

### 5. 收敛失败后的回退检查（第 506-518 行）
如果牛顿迭代失败，使用回退方案时检查：
```
⚠ [DEBUG] V is NaN/Inf after convergence failure:
```
显示哪些单元的 C1 系数是 NaN。

### 6. 最终边界检查详细信息（第 529-538 行）
在抛出边界错误前，打印完整的诊断信息：
```
⚠ [DEBUG] Voltage out of bounds error details:
```
包括：
- 无量纲和物理电压值
- 电压上下限
- 收敛状态和迭代次数
- 电流守恒检查
- 前5个单元的系数值

## 如何使用

重新运行 `testexample.jl`，查看输出中的 `⚠ [DEBUG]` 标记。这些信息会告诉您：

1. **NaN 的源头**：最早出现 NaN 的地方
2. **问题单元**：哪个热单元的参数有问题
3. **输入数据**：导致 NaN 的输入浓度或温度值

## 可能的原因

根据调试信息，NaN 可能来自：

1. **浓度超出范围**：
   - `cn_surf` 或 `cp_surf` < 0 或 > 1
   - `ce_n_gs` 或 `ce_p_gs` <= 0

2. **温度异常**：
   - `Te_prev[e]` 是 NaN 或负值
   - 导致 `Arrhenius` 函数返回 NaN

3. **交换电流密度过小**：
   - `prefactor_n` 或 `prefactor_p` 为 0 或 NaN
   - 导致 `alpha_p`, `alpha_n` 变成 Inf

4. **电导率计算问题**：
   - `kappa` 函数返回 NaN 或 0

## 下一步

1. 运行 `testexample.jl` 并查看调试输出
2. 根据打印的信息定位 NaN 的源头
3. 检查相应的输入参数或状态变量
4. 修复根本原因（可能需要调整初始条件、参数范围或添加保护措施）
