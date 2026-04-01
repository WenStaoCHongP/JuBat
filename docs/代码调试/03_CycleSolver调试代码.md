# CycleSolver.jl 调试代码记录

> 文件: `src/CycleSolver.jl`
> 状态: A (新增文件)
> 移除行数: ~6 行

---

## 1. 状态向量长度调试打印

### 1.1 放电相位结束后的调试输出 (约原 line 307-308)

**功能**: 打印放电结束后的状态向量长度，用于调试状态传递。

```julia
            y_out = get(current_state, "y", nothing)
            V_out = get(current_state, "V", NaN)
            @printf("    → V_out=%.3fV, y_len=%d\n", V_out, y_out === nothing ? 0 : length(y_out))
```

**恢复**: 粘贴回放电相位结果输出之后。

### 1.2 静置1相位结束后的调试输出 (约原 line 335-338)

**功能**: 打印静置1结束后的状态向量长度。

```julia
            V_out = get(current_state, "V", NaN)
            @printf("    → V_out=%.3fV, y_len=%d", V_out, y_out === nothing ? 0 : length(y_out))
```

以及紧接的 `println()`。

**恢复**: 粘贴回静置1结果输出之后。

### 1.3 充电相位前的调试输出 (约原 line 367-368)

**功能**: 打印充电前的状态向量长度。

```julia
            V_in = get(current_state, "V", NaN)
            @printf("(V_in=%.3fV, y_len=%d) ", V_in, y_in === nothing ? 0 : length(y_in))
```

**恢复**: 粘贴回充电相位开始之前。

---

## 2. 说明

CycleSolver.jl 中的大部分 println/@printf 是面向用户的仿真进度输出（如 `循环 1/50`、`开始充放电循环仿真` 等），属于功能性信息，不属于调试代码。

以上 3 处是唯一需要移除的调试代码，特点是:
- 显示内部状态向量长度 (`y_len=%d`)
- 仅在开发阶段用于验证状态传递正确性
- 对用户无意义
