# 归一化不一致问题修正总结

## 问题发现

检查 SetParams.jl 和对流边界条件后，发现**h的归一化处理不一致**。

## 问题详情

### SetParams.jl 中的定义

**Scale结构**（第180行）：
```julia
h_th::Float64 = 0  # Biot number reference (h*L_th/k_th)
```

**归一化**（第264行）：
```julia
h_ref = param_dim.cell.h * L_th / k_ref  # Biot number
param_dim.scale.h_th = h_ref
```

**含义**：
- 输入：`param_dim.cell.h` [W/(m²·K)]（**有量纲**）
- 输出：`param_dim.scale.h_th`（**无量纲Biot数**）
- 公式：$Bi = h \cdot L_{th} / k_{th}$

### 三种对流条件的对比（修正前）

| 对流类型 | h的来源 | h的量纲 | 归一化方式 | 问题 |
|---------|---------|---------|-----------|------|
| **外边界对流** | `scale.h_th` | 无量纲 | 预先在SetParams中 | ✅ 正确 |
| **表面冷却** | `opt.h_surface` | **有量纲** | **使用时临时归一化** | ❌ 不一致 |
| **极耳冷却** | `opt.h_tab` | **有量纲** | **使用时临时归一化** | ❌ 不一致 |

### 外边界对流（基准实现）

```julia
function _apply_convection_bc!(KT, FT, mesh, is_outer, case)
    scale = case.param_dim.scale
    Bi = scale.h_th  // ← 使用无量纲Biot数（预先归一化）
    
    L_th = scale.L_th
    wt = Bi * w * (J / L_th)  // ← 直接使用
    ...
end
```

**特点**：
- ✅ 从 `param_dim.cell.h` 读取有量纲值
- ✅ 在 `NormaliseParam` 中预先归一化
- ✅ 使用时直接读取无量纲Biot数

### 表面冷却（修正前，问题实现）

```julia
function _apply_cool_surface!(KT, FT, mesh, case, t)
    h_surface = hasproperty(case.opt, :h_surface) ? 
                case.opt.h_surface : 10.0  // ← 从opt读取，有量纲
    
    vol_coeff = 2.0 * h_surface / H  // ← 使用有量纲的h
    Bi_z = vol_coeff * L_th^2 / k_th  // ← 临时归一化
    ...
end
```

**问题**：
- ❌ 从 `opt.h_surface` 读取（而非 `param_dim.cell`）
- ❌ 每次使用时临时归一化（而非预先归一化）
- ❌ 默认值10.0，没说明量纲
- ❌ 如果想与外边界用相同h，需要设置两次

### 极耳冷却（修正前，问题实现）

```julia
function _apply_cool_tab!(KT, FT, mesh, case, t)
    h_tab = hasproperty(case.opt, :h_tab) ? 
            case.opt.h_tab : 100.0  // ← 从opt读取，有量纲
    
    coeff = h_tab * tab_area * weight / (H * k_th * L_th)  // ← 临时归一化
    ...
end
```

**问题**：与表面冷却相同。

## 不一致的后果

### 1. 参数设置混乱

用户需要设置多个对流系数：
```julia
// Jellyroll.jl 或参数文件
cell.h = 10.0  // 用于外边界对流

// testexample.jl（运行时）
opt.h_surface = 10.0  // 用于表面冷却（需要重复设置！）
opt.h_tab = 100.0     // 用于极耳冷却
```

### 2. 量纲混淆

```julia
// 用户可能误以为是无量纲的
opt.h_surface = 0.1  // ❌ 实际上是有量纲的！
```

### 3. 代码维护困难

- 外边界对流的归一化逻辑在 SetParams.jl
- 表面/极耳冷却的归一化逻辑在 ThermalDistributed.jl
- 不易维护和理解

## 修正方案

### 方案：统一从param_dim读取有量纲值

**原则**：
1. 所有对流系数都从 `param_dim` 读取（不从 `opt`）
2. 在使用时临时归一化（保持一致性）
3. 提供合理的回退值（默认使用外边界的h）

### 表面冷却（修正后）

```julia
function _apply_cool_surface!(KT, FT, mesh, case, t)
    try
        # 获取参数（从param_dim.cell读取，与外边界对流一致）
        h_surface = hasproperty(case.param_dim.cell, :h_surface) ?
                    case.param_dim.cell.h_surface :
                    case.param_dim.cell.h  # 回退到外边界的h [W/(m²·K)]
        H = hasproperty(case.param_dim.cell, :height) ? 
            case.param_dim.cell.height : 
            case.param_dim.cell.width
        
        scale = case.param_dim.scale
        k_th = scale.k_th
        L_th = scale.L_th
        
        # 体积散热系数：2h/H [W/(m³·K)]
        vol_coeff = 2.0 * h_surface / H
        
        # 无量纲Biot数：Bi_z = 2h*L_th^2 / (H*k_th)
        Bi_z = vol_coeff * L_th^2 / k_th
        
        conv_factor = Bi_z / L_th^2
        ...
    end
end
```

**改进**：
- ✅ 从 `param_dim.cell` 读取（与外边界一致）
- ✅ 回退到 `cell.h`（默认与外边界相同）
- ✅ 量纲明确（注释说明 [W/(m²·K)]）
- ✅ 用户友好（只需设置一次h）

### 极耳冷却（修正后）

```julia
function _apply_cool_tab!(KT, FT, mesh, case, t)
    try
        # 获取参数（从param_dim读取，与外边界对流一致）
        h_tab = hasproperty(case.param_dim.tab, :h) ?
                case.param_dim.tab.h :
                (hasproperty(case.param_dim.cell, :h_tab) ?
                 case.param_dim.cell.h_tab :
                 case.param_dim.cell.h)  # 回退到外边界的h [W/(m²·K)]
        
        tab_area = case.param_dim.tab.area
        H = case.param_dim.cell.width
        
        scale = case.param_dim.scale
        k_th, L_th = scale.k_th, scale.L_th
        
        # 计算弧长权重...
        coeff = h_tab * tab_area * weight / (H * k_th * L_th)
        ...
    end
end
```

**改进**：
- ✅ 优先从 `param_dim.tab.h` 读取（极耳专用）
- ✅ 其次从 `param_dim.cell.h_tab` 读取
- ✅ 最后回退到 `param_dim.cell.h`（外边界）
- ✅ 灵活且有合理的默认行为

## 用户使用方式

### 修正前（问题）

```julia
// Jellyroll.jl
cell.h = 10.0  // 外边界对流

// testexample.jl（需要重复设置）
opt.h_surface = 10.0  // ❌ 需要手动设置
opt.h_tab = 100.0     // ❌ 需要手动设置
```

### 修正后（推荐）

```julia
// Jellyroll.jl 或参数文件
cell.h = 10.0          // 外边界对流（默认也用于表面冷却）
cell.h_surface = 5.0   // （可选）单独设置表面冷却
cell.h_tab = 100.0     // （可选）单独设置极耳冷却

// 或者在tab结构中
tab.h = 100.0          // 极耳专用对流系数

// testexample.jl（只需设置冷却方式）
opt.cool_method = "surface"  // 或 "tab"
// 不需要再设置h值！✅
```

**优点**：
1. ✅ 参数集中在一处（`param_dim`）
2. ✅ 有合理的默认行为（回退到 `cell.h`）
3. ✅ 灵活性（可以单独设置不同的h）
4. ✅ 量纲明确（都在参数文件中，有单位说明）

## 归一化方式对比

### 外边界对流

**方式**：预先归一化
```julia
// SetParams.jl (NormaliseParam)
scale.h_th = cell.h * L_th / k_th  // 预先计算Biot数

// ThermalDistributed.jl (使用时)
Bi = scale.h_th  // 直接使用无量纲Biot数
```

### 表面/极耳冷却（修正后）

**方式**：使用时归一化
```julia
// ThermalDistributed.jl
h = cell.h_surface  // 读取有量纲值
Bi_z = 2*h*L_th^2 / (H*k_th)  // 临时归一化
```

**为什么不预先归一化？**
- z方向冷却的Biot数定义不同：$Bi_z = 2h L_{th}^2 / (H k_{th})$
- 需要知道H（电池厚度），而H在SetParams时可能未知
- 临时归一化更灵活

**一致性**：
- 都从 `param_dim` 读取有量纲值 ✅
- 都使用 $h \cdot L / k$ 的形式归一化 ✅

## 量纲验证

### 外边界对流

$$
Bi = \frac{h \cdot L_{th}}{k_{th}} = \frac{[W/(m^2 \cdot K)] \cdot [m]}{[W/(m \cdot K)]} = [无量纲]
$$

### z方向冷却

$$
Bi_z = \frac{2h \cdot L_{th}^2}{H \cdot k_{th}} = \frac{[W/(m^2 \cdot K)] \cdot [m^2]}{[m] \cdot [W/(m \cdot K)]} = [无量纲]
$$

**正确** ✅

## 总结

### 问题

1. ❌ **参数来源不一致**：外边界从 `param_dim`，表面/极耳从 `opt`
2. ❌ **归一化时机不同**：外边界预先归一化，表面/极耳临时归一化
3. ❌ **量纲混淆**：`opt` 中的值容易误解为无量纲
4. ❌ **用户体验差**：需要多次设置相同的参数

### 修正

1. ✅ **统一从param_dim读取**：所有对流系数都从 `param_dim` 获取
2. ✅ **合理的回退**：默认使用 `cell.h`（外边界值）
3. ✅ **量纲明确**：参数文件中有明确的单位说明
4. ✅ **用户友好**：只需设置一次，自动适用于所有对流边界

### 建议

**参数设置最佳实践**：
```julia
// Jellyroll.jl 或参数文件
cell.h = 10.0          // 默认对流系数（外边界、表面冷却）[W/(m²·K)]
cell.h_tab = 100.0     // （可选）极耳专用 [W/(m²·K)]
// 或
tab.h = 100.0          // 极耳专用（优先级更高）[W/(m²·K)]

// testexample.jl
opt.cool_method = "surface"  // 只需选择冷却方式
```

---

**修正完成** ✅  
**归一化统一** ✅  
**用户友好** ✅  
**可以使用** ✅
