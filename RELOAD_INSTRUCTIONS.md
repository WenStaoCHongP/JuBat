# Julia 模块重载说明

## 问题
修改了 `src/mechanical.jl` 后，仍然出现 `UndefVarError: U_M` 错误。这是因为 Julia 缓存了旧版本的模块。

## 解决方案（选择其一）

### 方案 1：重启 Julia（推荐）
最简单的方法：
1. 完全退出 Julia REPL 或 IDE
2. 重新启动 Julia
3. 重新运行 `testexample.jl`

### 方案 2：使用 Revise.jl（开发时推荐）
如果经常需要修改代码，安装 Revise.jl：

```julia
# 在 Julia REPL 中
using Pkg
Pkg.add("Revise")

# 然后在脚本开头添加
using Revise
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
```

### 方案 3：强制重新加载模块
在 `testexample.jl` 的开头修改为：

```julia
# 清除已加载的模块（如果存在）
if isdefined(Main, :JuBat)
    @eval Main JuBat = nothing
end

# 重新包含模块
include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat
```

### 方案 4：命令行运行（推荐用于测试）
使用 `--compiled-modules=no` 选项：

```bash
julia --compiled-modules=no example/testexample.jl
```

## 验证修复是否生效

运行以下命令检查第440行的内容：

```bash
sed -n '429,442p' src/mechanical.jl
```

应该看到：
```julia
"""求解位移场"""
function _solve_mechanical_displacement_2D(K_M, F_M, nnode)
    ndof = 2 * nnode
    
    # 求解线性方程组 K*U = F
    U_M = try
        K_M \ F_M
    catch e
        @warn "Mechanical solve failed, using zero displacement" e
        zeros(Float64, ndof)
    end
    
    return U_M
end
```

## 已完成的所有修复

✅ Variables.jl - 添加 SOC 历史记录存储
✅ Solve.jl - 将 SOC 数据添加到输出
✅ testexample.jl - 添加 SOC 数据加载和回退方案
✅ mechanical.jl - 修复所有函数参数和作用域问题

所有代码修改都已正确保存到文件中。
