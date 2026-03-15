# 热验证脚本归一化修正计划

## 目标

审查并修正三个热验证脚本，确保它们与新的统一归一化方案完全一致：
- `example/热模块验证/thermal_verify.jl`
- `example/热模块验证/thermal_error_source_analysis.jl`
- `example/热模块验证/thermal_equivalent_lumped_compare.jl`

## 背景

根据 `md/01_参数定义与归一化.md`，JuBat 现在采用统一归一化方案：

### 新归一化方案（2026-03-12 更新）
1. **统一能量尺度**：P_ref = φ_ref × I_typ
2. **统一时间尺度**：t_0 = 3600 s（电化学与热模型统一）
3. **密度与比热容分离归一化**：
   - ρ* = ρ / ρ_ref
   - c* = c × ρ_ref × L³ × T_ref / (t_0 × P_ref)
4. **导热率归一化**：k* = k × L × T_ref / P_ref
5. **热源归一化**：q* = q × L³ / P_ref（即 q* = q / scale.q）

### 旧归一化方案（可能存在的问题）
- 使用传统热尺度 ρ_th, c_th, k_th
- 热模型与电化学模型使用不同时间尺度
- 体积热容未分离归一化

## 初步审查结果

### thermal_verify.jl
- ✓ Line 334, 398, 502: 使用 `q_ref = scale.q`（正确）
- ✓ Line 328-329: 使用物理时间单位（正确，Solve.jl 会自动归一化）
- ⚠️ 需要检查：update_fn 中的时间参数是否正确处理

### thermal_error_source_analysis.jl
- ✓ Line 155: 使用 `q_ref = case.param_dim.scale.q`（正确）
- ✓ Line 104-105: 使用物理时间单位（正确）
- ⚠️ 需要检查：边界条件计算是否使用正确的归一化

### thermal_equivalent_lumped_compare.jl
- ✓ Line 102: 使用 `q_ref = case.param_dim.scale.q`（正确）
- ✓ Line 63-64: 使用物理时间单位（正确）
- ⚠️ 需要检查：集总模型与分布式模型的一致性

## 任务分解

### Phase 1: 深入代码审查 [pending]
**目标**：逐行审查三个文件，识别所有与归一化相关的代码

**检查项**：
1. 热源计算与归一化
   - 是否使用 `scale.q` 进行归一化？
   - 是否有遗留的旧归一化公式？
2. 时间尺度处理
   - `opt.time` 和 `opt.dt` 是否使用物理单位？
   - update_fn 接收的时间参数是物理时间还是无量纲时间？
   - 结果输出的时间是否需要转换？
3. 边界条件
   - 传热系数 h 的使用是否正确？
   - 是否需要归一化处理？
4. 热容与导热率
   - 是否直接使用 param_dim 中的物理参数？
   - 是否有不必要的手动归一化？
5. 体积与面积计算
   - 网格积分权重是否正确？
   - 是否需要乘以长度尺度 L？

**输出**：findings.md 中记录所有发现的问题

### Phase 2: 识别问题模式 [pending]
**目标**：分类问题类型，确定修复策略

**分类**：
1. 热源归一化问题
2. 时间尺度混用问题
3. 边界条件归一化问题
4. 几何尺度问题
5. 其他维度不一致问题

**输出**：findings.md 中的问题分类表

### Phase 3: 制定修复方案 [pending]
**目标**：为每个问题制定具体修复方案

**修复策略**：
1. 统一使用 `scale.q` 进行热源归一化
2. 确保 update_fn 正确处理时间参数
3. 检查边界条件计算是否使用物理参数
4. 验证几何积分是否正确
5. 添加注释说明归一化约定

**输出**：每个文件的修复清单

### Phase 4: 实施修复 [pending]
**目标**：按照修复方案修改代码

**步骤**：
1. 备份原始文件
2. 逐文件实施修复
3. 添加注释说明归一化约定
4. 更新文档字符串

**输出**：修改后的三个文件

### Phase 5: 验证修复 [pending]
**目标**：运行修复后的脚本，验证结果正确性

**验证项**：
1. 脚本能否正常运行？
2. 输出结果的量纲是否正确？
3. 与精确解的误差是否在合理范围内？
4. 不同方法之间的一致性如何？

**输出**：验证报告

### Phase 6: 文档更新 [pending]
**目标**：更新相关文档，记录修复内容

**更新内容**：
1. 在文件头部添加归一化说明注释
2. 更新 CLAUDE.md 中的验证方案说明（如需要）
3. 记录修复历史

**输出**：更新后的文档

## 关键参考

### 归一化公式（md/01_参数定义与归一化.md）
- 功率参考：P_ref = φ_ref × I_typ
- 热源尺度：scale.q = P_ref / L³
- 时间尺度：scale.t0 = 3600 s
- 温度尺度：scale.T_ref = 298 K
- 长度尺度：scale.L = t_PE + t_SP + t_NE

### 热模型时间处理（md/05_热模型_二维分布式.md, lines 126-142）
```julia
if case.opt.model == "thermal"
    t0 = case.opt.time[1]/case.param.scale.t0      # 无量纲初始时间
    t_end = case.opt.time[end]/case.param.scale.t0   # 无量纲终止时间
    dt = case.opt.dt[1]/case.param.scale.t0        # 无量纲时间步长
end
```

### 验证脚本使用说明（md/01_参数定义与归一化.md, lines 571-587）
```julia
# 输入时间需要转换为无量纲
t_vec_nd = collect(t_phys) ./ param.scale.t0

# 热源更新函数接收无量纲时间
function update_heat_source(t_nd, vars)
    # t_nd 是无量纲时间
    vars["heat_source_fields"] = interpolate_heat_source(t_nd)
end

# 求解结果的时间需转换回物理时间
ring_sol = JuBat.Solve(case)
t_phys = ring_sol.time .* param.scale.t0
```

## 错误记录

（暂无）

## 决策记录

1. **审查策略**：先全面审查，再分类问题，最后统一修复
2. **修复原则**：最小化修改，保持代码可读性，添加注释说明
3. **验证标准**：脚本能运行 + 量纲正确 + 误差合理

## 状态

- 当前阶段：Phase 1 (深入代码审查)
- 开始时间：2026-03-15
- 预计完成：待定
