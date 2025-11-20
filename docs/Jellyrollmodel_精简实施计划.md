# Jellyrollmodel.jl 精简实施计划

## 总览

**目标**: 将 694 行代码精简到 ~520 行（减少 25%）  
**主要方法**: 消除重复、提取辅助函数、向量化计算  
**预估时间**: 2-4 小时  

---

## 实施路线图

### 📋 阶段1：消除重复逻辑（高优先级，1-2小时）

#### 任务1.1：精简 `material_at` 函数 ⭐⭐⭐
**时间**: 15分钟  
**减少**: 25行

**步骤**:
1. 在 `material_at` 函数前添加辅助函数：
```julia
function _find_layer_in_period(offset::Float64, order)
    acc = 0.0
    for (name, w) in order
        if offset <= acc + w
            return name, offset - acc
        end
        acc += w
    end
    last_layer, last_width = order[end]
    total_width = sum(w for (_, w) in order)
    return last_layer, offset - (total_width - last_width)
end
```

2. 替换 `material_at` 函数主体（Lines 113-139）为：
```julia
function material_at(r::Real, θ::Real, p; logic::Symbol=:spiral)
    r <= p.Rin && return :inner, 0.0
    r >= p.Rout && return :outer, 0.0
    
    offset = if logic === :spiral
        mod(r - (p.a + p.b*θ), p.t_repeat)
    elseif logic === :rings
        mod(r - p.Rin, p.t_repeat)
    else
        error("Unknown logic=$(logic). Use :spiral or :rings")
    end
    
    return _find_layer_in_period(offset, p.order)
end
```

3. 测试验证：
```julia
# 测试用例
p = jellyroll_spiral_params(param_dim)
@assert material_at(0.02, 0.0, p; logic=:spiral)[1] === material_at(0.02, 0.0, p; logic=:rings)[1]
```

**验证点**: 
- ✓ 两种逻辑模式结果一致
- ✓ 边界情况处理正确
- ✓ 代码从 ~40 行减少到 ~15 行

---

#### 任务1.2：精简 `edge_boundary` 函数 ⭐⭐⭐
**时间**: 20分钟  
**减少**: 20行

**步骤**:
1. 替换 Lines 298-352 为精简版本：
```julia
function edge_boundary(mesh, nidx::Int, param_dim; 
                       which::Symbol=:inner, 
                       theta_range::Union{Tuple{Float64,Float64},Nothing}=nothing, 
                       tol::Float64=1e-4)
    # 获取参数和偏移
    p = jellyroll_spiral_params(param_dim)
    offset = if which === :inner
        0.0
    elseif which === :outer
        p.t_repeat
    else
        error("which must be :inner or :outer")
    end
    
    # 获取节点坐标
    x, y = mesh.node[nidx, 1], mesh.node[nidx, 2]
    r = hypot(x, y)
    bval = max(p.b, 1e-12)
    
    # 获取 θ 范围
    if theta_range === nothing
        θ_min, θ_max = which === :inner ? 
            (0.0, 2.0*π) : 
            (2.0*π * max(0, Int(p.n_wind)-1), 2.0*π * Int(p.n_wind))
    else
        θ_min, θ_max = theta_range
    end
    
    # 计算并验证
    θ_cum = (r - p.a - offset) / bval
    θ_cum < θ_min || θ_cum > θ_max && return false
    
    r_theo = p.a + p.b * θ_cum + offset
    dist = hypot(x - r_theo*cos(θ_cum), y - r_theo*sin(θ_cum))
    return dist <= tol
end
```

2. 测试验证：
```julia
# 测试：边界节点识别应该一致
@assert edge_boundary(mesh, 1, param_dim; which=:inner) == 
        edge_boundary_old(mesh, 1, param_dim; which=:inner)
```

**验证点**:
- ✓ 内外边界识别正确
- ✓ 容差判断准确
- ✓ 代码从 ~55 行减少到 ~35 行

---

#### 任务1.3：精简 `jellyroll_tab_node_indices` 函数 ⭐⭐⭐
**时间**: 30-45分钟  
**减少**: 60行

**步骤**:
1. 在函数前添加两个辅助函数：
```julia
# 辅助函数1：通用节点查找
function _find_tab_nodes(mesh, tab_angles, θ_cum_nodes, θ_cum_range, 
                         delta_theta_fn, tw, Rin, Rout; reverse_range=false)
    # ... 见精简示例 ...
end

# 辅助函数2：角度增量计算
function _delta_theta_from_width(a, b, θ0, width)
    # ... 见精简示例 ...
end
```

2. 替换 `jellyroll_tab_node_indices` 主体（Lines 552-692）

3. 测试验证：
```julia
# 测试：节点索引应该一致
pos_old, neg_old = jellyroll_tab_node_indices_old(mesh, param_dim)
pos_new, neg_new = jellyroll_tab_node_indices(mesh, param_dim)
@assert Set(pos_old) == Set(pos_new)
@assert Set(neg_old) == Set(neg_new)
```

**验证点**:
- ✓ 正负极耳节点识别一致
- ✓ 边界情况处理正确
- ✓ 代码从 ~140 行减少到 ~80 行

---

### 📋 阶段2：优化计算（中优先级，30分钟）

#### 任务2.1：向量化 `jellyroll_spiral_params` ⭐⭐
**时间**: 15分钟  
**减少**: 15行

**步骤**:
1. 替换 Lines 49-98 为向量化版本（见精简示例）
2. 测试所有返回值一致性

**验证点**:
- ✓ 所有参数数值完全相同
- ✓ 代码更简洁

---

#### 任务2.2：提取 Q4 形函数为常量 ⭐⭐
**时间**: 5分钟  
**减少**: 不减少行数，但避免重复定义

**步骤**:
1. 在文件顶部（Line 40 后）添加：
```julia
# Q4 等参元形函数（常量）
const N_Q4 = (xi, eta) -> 0.25 .* [
    (1 - xi)*(1 - eta),
    (1 + xi)*(1 - eta),
    (1 + xi)*(1 + eta),
    (1 - xi)*(1 + eta)
]
```

2. 在 `jellyroll_element_layer_weights` 中删除局部定义（Lines 374-379）
3. 直接使用常量 `N_Q4`

**验证点**:
- ✓ 层权重计算结果不变

---

#### 任务2.3：向量化 `jellyroll_element_centers` ⭐
**时间**: 5分钟  
**减少**: 5行

**步骤**:
1. 替换 Lines 197-207 为：
```julia
function jellyroll_element_centers(mesh)
    ne = size(mesh.element, 1)
    return [mean(mesh.node[mesh.element[e, :], d]) 
            for e in 1:ne, d in 1:2]
end
```

**验证点**:
- ✓ 单元中心坐标完全相同

---

### 📋 阶段3：清理代码（低优先级，15分钟）

#### 任务3.1：删除未使用内容 ⭐
**时间**: 5分钟  
**减少**: 10行

**步骤**:
1. 删除 Line 8: `using Plots`
2. 评估是否删除 `pol2cart` 函数（Lines 26-32）

**验证点**:
- ✓ 代码仍能正常加载
- ✓ 所有测试通过

---

#### 任务3.2：代码格式和注释 ⭐
**时间**: 10分钟

**步骤**:
1. 添加精简后函数的文档字符串
2. 更新注释说明优化内容
3. 统一代码风格

---

## 测试计划

### 单元测试清单

```julia
# 测试文件：test/test_jellyrollmodel_refactor.jl

@testset "Jellyrollmodel 精简测试" begin
    # 1. material_at 测试
    @testset "material_at" begin
        p = jellyroll_spiral_params(param_dim)
        r, θ = 0.02, 0.5
        @test material_at(r, θ, p; logic=:spiral) == material_at_old(r, θ, p; logic=:spiral)
        @test material_at(r, θ, p; logic=:rings) == material_at_old(r, θ, p; logic=:rings)
    end
    
    # 2. edge_boundary 测试
    @testset "edge_boundary" begin
        for i in 1:10:mesh.nlen
            @test edge_boundary(mesh, i, param_dim; which=:inner) == 
                  edge_boundary_old(mesh, i, param_dim; which=:inner)
        end
    end
    
    # 3. jellyroll_tab_node_indices 测试
    @testset "tab_node_indices" begin
        pos, neg = jellyroll_tab_node_indices(mesh, param_dim)
        pos_old, neg_old = jellyroll_tab_node_indices_old(mesh, param_dim)
        @test Set(pos) == Set(pos_old)
        @test Set(neg) == Set(neg_old)
    end
    
    # 4. spiral_params 测试
    @testset "spiral_params" begin
        p = jellyroll_spiral_params(param_dim)
        p_old = jellyroll_spiral_params_old(param_dim)
        for field in fieldnames(typeof(p))
            @test getfield(p, field) ≈ getfield(p_old, field)
        end
    end
end
```

### 集成测试

运行现有的完整示例：
```bash
julia example/spme_thermal2d_example.jl
julia example/jellyroll_coupled_example.jl
```

**验证**:
- ✓ 所有示例正常运行
- ✓ 结果与精简前一致
- ✓ 无性能回退

---

## 实施时间表

| 阶段 | 任务 | 时间 | 累计 |
|------|------|------|------|
| **阶段1** | | | |
| | 1.1 material_at | 15min | 15min |
| | 1.2 edge_boundary | 20min | 35min |
| | 1.3 tab_node_indices | 45min | 80min |
| **阶段2** | | | |
| | 2.1 spiral_params | 15min | 95min |
| | 2.2 Q4 形函数 | 5min | 100min |
| | 2.3 element_centers | 5min | 105min |
| **阶段3** | | | |
| | 3.1 删除未使用 | 5min | 110min |
| | 3.2 文档和格式 | 10min | 120min |
| **测试** | | | |
| | 单元测试 | 20min | 140min |
| | 集成测试 | 10min | 150min |
| **总计** | | **2.5小时** | |

---

## 风险评估

### 低风险 ✅
- 提取辅助函数（不改变逻辑）
- 向量化计算（Julia 原生支持）
- 删除未使用代码

### 中风险 🟡
- `jellyroll_tab_node_indices` 重构（逻辑复杂）
  - **缓解**: 详细单元测试
  - **回退**: 保留原函数作为备份

### 零风险 ✅
- 提取常量定义
- 代码格式优化
- 添加注释

---

## 回滚计划

### 准备工作
1. 创建备份分支：
```bash
git checkout -b jellyrollmodel-refactor
git commit -m "Backup before Jellyrollmodel refactoring"
```

2. 保留原函数（临时）：
```julia
# 在文件末尾
function material_at_old(...) ... end
function edge_boundary_old(...) ... end
# ...
```

### 回滚步骤
如果测试失败：
1. 运行测试确定失败的函数
2. 恢复该函数到原版本
3. 重新测试
4. 分析失败原因

### 完全回滚
```bash
git checkout main
git branch -D jellyrollmodel-refactor
```

---

## 成功标准

### 必须满足
- ✅ 所有单元测试通过
- ✅ 所有集成测试通过
- ✅ 结果与原代码一致（数值误差 < 1e-12）

### 期望达成
- ✅ 代码行数减少 20-25%
- ✅ 重复代码消除 100%
- ✅ 函数平均长度 < 50 行
- ✅ 无性能回退（允许 ±5%）

### 额外收益
- ✅ 代码可读性提升
- ✅ 更易维护
- ✅ 更易测试

---

## 后续优化

完成基础精简后，可以考虑：

### 短期（1周内）
1. 添加更多单元测试
2. 性能基准测试
3. 文档更新

### 中期（1月内）
4. 拆分超长函数（如 `jellyroll_collector_seed_mesh`）
5. 提取更多通用函数
6. 优化算法（如边界搜索）

### 长期（按需）
7. 引入类型注解（性能优化）
8. 并行化计算（大规模网格）
9. GPU 加速（可选）

---

## 检查清单

### 开始前 ☑
- [ ] 创建备份分支
- [ ] 保存当前测试结果
- [ ] 准备测试数据

### 每个任务后 ☑
- [ ] 运行相关单元测试
- [ ] 检查代码格式
- [ ] 更新注释
- [ ] Commit 代码

### 完成后 ☑
- [ ] 运行所有测试
- [ ] 对比性能
- [ ] 更新文档
- [ ] Code Review
- [ ] Merge 到主分支

---

**编写日期**: 2025-11-19  
**预估时间**: 2.5 小时  
**预期成果**: 减少 ~170 行代码（25%）  
**风险等级**: 低-中
