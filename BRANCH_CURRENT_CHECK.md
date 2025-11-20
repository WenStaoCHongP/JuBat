# 分电流求解公式验证

## 📋 文档 vs 代码对比

### 文档中的公式（求分电流办法.md）

**简化后的电压方程：**
```
V_cell = C1 + C2 * [sinh^{-1}(C3*x) + sinh^{-1}(C4*x)] - C5*x
```

**系数定义：**
```
C1 = U_p - U_n + (T-T0)*dUdT_p - (T-T0)*dUdT_n + 2T(1-t+)(cs_p_av - cs_n_av)/ce0
C2 = 2T
C3 = -1/(as_p * thickness_p)
C4 = 1/(as_n * thickness_n)
C5 = R_EL + (thickness_n/sig_n + thickness_p/sig_p)/3
```

### 代码中的实现（SPMe.jl）

**电压方程（第430-433行）：**
```julia
function branch_voltage(coeff, I::Float64)
    apI = coeff.alpha_p * I
    anI = coeff.alpha_n * I
    return coeff.C1 + coeff.C2 * (asinh(apI) - asinh(anI)) - coeff.C5 * I
end
```

**系数定义（第392-396行）：**
```julia
C1 = (u_p_val(T_e) - u_n_val(T_e)) + 2.0*T_e*(1.0-param.EL.tplus)*(csp_av-csn_av)/param.EL.ce0
C2 = 2.0 * T_e
alpha_p = -1.0 / (2.0 * j0_p * param.PE.as * param.PE.thickness)
alpha_n =  1.0 / (2.0 * j0_n * param.NE.as * param.NE.thickness)
C5 = R_EL + c_sigma
```

## ❌ 发现的差异

### 差异 1：符号（加号 vs 减号）

**文档第20行：**
```
sinh^{-1}(C3*x) + sinh^{-1}(C4*x)  // 加号
```

**但文档第8行（物理公式）：**
```
η_r = 2T · [sinh^{-1}(j_p) - sinh^{-1}(j_n)]  // 减号！
```

**代码：**
```julia
asinh(apI) - asinh(anI)  // 减号
```

**结论：** 代码符合第8行的物理公式，文档第20行可能有误。代码是**正确的**。✓

### 差异 2：alpha 的定义（关键问题！）⚠️

**文档：**
```
C3 = -1/(as_p * thickness_p)
C4 = 1/(as_n * thickness_n)
```

**代码：**
```julia
alpha_p = -1.0 / (2.0 * j0_p * as_p * thickness_p)  // 多了 2*j0_p
alpha_n =  1.0 / (2.0 * j0_n * as_n * thickness_n)  // 多了 2*j0_n
```

**问题：** 代码中多了 `2*j0` 因子！

## 🔍 理论推导验证

### Butler-Volmer 方程

完整的 BV 方程：
```
j = j0 * [exp(αa*F*η/RT) - exp(-αc*F*η/RT)]
```

对称情况（αa = αc = 0.5）：
```
j = j0 * [exp(0.5*F*η/RT) - exp(-0.5*F*η/RT)]
   = 2*j0 * sinh(0.5*F*η/RT)
```

反解：
```
η = (2RT/F) * sinh^{-1}(j/(2*j0))
```

无量纲化（设置 RT/F 为参考电势）：
```
η* = 2 * sinh^{-1}(j/(2*j0))
```

### 对应到代码

```
j_p = -I/(as_p * thickness_p)  // 电流密度
η_p = 2T * sinh^{-1}(j_p / (2*j0_p))
    = 2T * sinh^{-1}(-I / (2*j0_p * as_p * thickness_p))
    = 2T * asinh(alpha_p * I)  // 其中 alpha_p = -1/(2*j0_p*as_p*thickness_p)
```

**结论：** 代码中的 `2*j0` 因子是**正确的**！文档简化公式时可能遗漏了这个因子。✓

## ✅ 公式验证结果

| 项目 | 文档 | 代码 | 状态 |
|------|------|------|------|
| C1 定义 | U_p - U_n + ... | ✓ 一致 | ✓ 正确 |
| C2 定义 | 2T | ✓ 一致 | ✓ 正确 |
| C5 定义 | R_EL + ... | ✓ 一致 | ✓ 正确 |
| asinh 符号 | + (第20行) | - (减号) | ⚠️ 文档第20行有误，应该是减号 |
| alpha 定义 | -1/(as*thick) | -1/(2*j0*as*thick) | ✓ 代码正确（BV方程） |

## 🎯 结论

**分电流求解函数的公式是正确的！**

1. ✅ 系数 C1, C2, C5 定义正确
2. ✅ alpha 包含 2*j0 因子是正确的（来自 Butler-Volmer 方程）
3. ✅ asinh 的符号（减号）是正确的（符合物理公式第8行）
4. ⚠️ 文档第20行的加号可能是笔误

**分电流求解本身没有问题！**

问题可能在其他地方（比如约束条件的实现、收敛判据等）。

## 🔍 约束条件验证

### 文档中的约束
```
sum(x_i * A_i) = R
```

### 代码中的实现（第513行）
```julia
res_I = sum(w .* I_e) - I_total
```
其中 `w = areas / A_global`

### 等价性证明
```
sum(w .* I_e) = sum((A_i/A_total) * I_e[i])
              = sum(A_i * I_e[i]) / A_total
```

要求：`sum(w .* I_e) = I_total`

等价于：`sum(A_i * I_e[i]) / A_total = I_total`

即：`sum(A_i * I_e[i]) = I_total * A_total`

**结论：** 约束条件**正确** ✓（考虑到 I_total 是归一化的电流密度）

## ⚠️ 发现的潜在问题

### 问题：分电流求解器内部的电压边界检查（第573-582行）

```julia
if !(V_MIN <= V <= V_MAX)
    throw(ErrorException("thermal2D common voltage out of bounds..."))
end
```

**影响：**
- 如果求解的公共电压 V 超出范围，**立即抛出异常**
- 这会导致整个 Solve 过程终止
- 可能就是"只运行1步"的原因

**建议：**
1. 检查初始电压是否接近边界（4.1804 vs 4.2）
2. 如果是充电，第一步可能就超限
3. 考虑放宽电压范围或改为放电

## 💡 下一步建议

重新运行测试并查看新增的调试输出：
```
[Solve] 初始化完成
  初始电压: 4.1804 V
  电压范围: [2.5, 4.2] V  ← 只差0.02V

[Solve] 迭代 1: t=0.01s, V=4.21V, dt=0.01

[INFO] 电压超出范围，终止求解  ← 可能在这里
  终止原因: 高于上限
```

或者：
```
ERROR: thermal2D common voltage out of bounds...  ← 分电流求解器抛出
```

根据输出可以确定是主循环检查还是分电流求解器检查触发的终止。
