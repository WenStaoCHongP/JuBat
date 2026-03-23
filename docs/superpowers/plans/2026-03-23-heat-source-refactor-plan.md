# 热源计算代码重构实现计划

> **Goal:** 将 Solve.jl 第633-744行的热源计算逻辑合并到 ThermalDistributed.jl，简化量纲转换，统一接口。

**Architecture:** 在 ThermalDistributed.jl 中创建统一的 `compute_heat_sources` 函数，支持单 SPMe 和多 SPMe 模式。计算过程保持无量纲，在后处理阶段转换为物理单位。

**Tech Stack:** Julia

---

## File Structure

| 文件 | 操作 | 说明 |
|------|------|------|
| `src/Variables.jl` | 修改 | 重命名热源变量（去掉单位后缀） |
| `src/ThermalDistributed.jl` | 修改 | 新增 `compute_heat_sources`，删除旧函数 |
| `src/Solve.jl` | 修改 | 替换热源计算代码为函数调用 |
| `src/PostProcessing.jl` | 修改 | 添加热源物理单位转换 |

---

## Task Structure

### Task 1: 修改 Variables.jl 热源变量定义

**Files:**
- Modify: `src/Variables.jl` 第93-130行

- [ ] **Step 1: 修改热源变量名（去掉 [W/m3] 后缀）**

将第98-108行的变量名从带单位改为不带单位：

```julia
# 修改前
variables["thermal2D Q_rxn_NE [W/m3]"] = zeros(Float64, ne, num)
# 修改后
variables["thermal2D q_rxn_ne"] = zeros(Float64, ne, num)
```

- [ ] **Step 2: 运行测试确认变量名变更不影响其他代码**

```bash
# 搜索使用这些变量的地方
grep -r "thermal2D Q_rxn_NE" src/
grep -r "thermal2D Q_rev_NE" src/
# ... 其他变量 ...
```

---

### Task 2: 在 ThermalDistributed.jl 中创建统一函数

**Files:**
- Modify: `src/ThermalDistributed.jl`

- [ ] **Step 1: 创建 `compute_heat_sources` 函数**

在 ThermalDistributed.jl 中添加新函数（约100行）：

```julia
function compute_heat_sources(case::Case, variables::Dict,
                              variables_elems::Union{Vector{<:Dict}, Nothing},
                              I_e::Vector{Float64}, T_e::Vector{Float64},
                              areas::Vector{Float64}; per_element_spme::Bool=false)
    mesh = case.mesh["thermal2D"]
    ne = size(mesh.element, 1)
    param = case.param

    # 获取层权重
    fks = jellyroll_element_properties(mesh, param)[2]

    # 从 variables 获取预分配的数组
    q_rxn_ne = variables["thermal2D q_rxn_ne"]
    q_rev_ne = variables["thermal2D q_rev_ne"]
    # ... 其他热源数组 ...

    for e in 1:ne
        # 获取电化学变量（根据 per_element_spme 判断）
        if per_element_spme
            vars_e = variables_elems[e]
            eta_n = vars_e["negative electrode overpotential"][1]
            # ... 其他变量 ...
        else
            eta_n = variables["negative electrode overpotential"][1]
            # ... 其他变量 ...
        end

        # 计算各层热源分量（无量纲）
        q_rxn_ne[e] = fks[e,1] * param.NE.as * abs(j_n) * abs(eta_n)
        # ... 其他计算 ...
    end

    # 计算总热源
    q_total = q_rxn_ne + q_rev_ne + ...
    variables["heat_source_fields"] = q_total
    variables["total heat source"] = [sum(q_total .* areas)]

    return variables
end
```

- [ ] **Step 2: 添加 `compute_heat_sources_with_czm` 包装函数**

```julia
function compute_heat_sources_with_czm(case::Case, variables::Dict,
                                        variables_elems::Union{Vector{<:Dict}, Nothing},
                                        I_e::Vector{Float64}, T_e::Vector{Float64},
                                        areas::Vector{Float64}, czm_mesh, mesh_data)
    variables = compute_heat_sources(case, variables, variables_elems, I_e, T_e, areas; per_element_spme=true)

    # 获取活跃单元
    active_elements = get_active_elements(czm_mesh, mesh_data)
    q_total = variables["heat_source_fields"]

    for e in 1:length(q_total)
        if !(e in active_elements)
            q_total[e] = 0.0
        end
    end

    variables["heat_source_fields"] = q_total
    variables["active_elements"] = active_elements

    return variables
end
```

- [ ] **Step 3: 删除旧函数**

删除以下函数：
- `heatQ_Source` (第281-297行)
- `compute_element_heat_sources` (第299-345行)
- `heatQ_Source_with_czm` (第347-375行)

---

### Task 3: 修改 Solve.jl 热源计算调用

**Files:**
- Modify: `src/Solve.jl` 第633-744行

- [ ] **Step 1: 替换热源计算代码为函数调用**

将第633-744行的热源计算代码替换为：

```julia
    # 6) 计算逐单元热源（调用 ThermalDistributed.jl 中的统一函数）
    if case.opt.czm_enabled && czm_mesh !== nothing
        variables = compute_heat_sources_with_czm(
            case, variables, variables_elems, I_e, Te_prev, areas, czm_mesh, mesh_th)
    else
        variables = compute_heat_sources(
            case, variables, variables_elems, I_e, Te_prev, areas; per_element_spme=true)
    end

    # 保存辅助变量（用于调试）
    for e in 1:ne
        vars_e = variables_elems[e]
        variables["thermal2D eta_n_e"][e] = vars_e["negative electrode overpotential"][1]
        variables["thermal2D eta_p_e"][e] = vars_e["positive electrode overpotential"][end]
        cn_surf_e = vars_e["negative particle surface lithium concentration"][1]
        cp_surf_e = vars_e["positive particle surface lithium concentration"][end]
        variables["thermal2D dUdT_n_e"][e] = param.NE.dUdT(cn_surf_e)[1]
        variables["thermal2D dUdT_p_e"][e] = param.PE.dUdT(cp_surf_e)[1]
        csn_data = vars_e["negative particle lithium concentration"]
        csp_data = vars_e["positive particle lithium concentration"]
        variables["thermal2D element soc_n"][e] = mean(vec(csn_data))
        variables["thermal2D element soc_p"][e] = mean(vec(csp_data))
    end
    variables["thermal2D element current"] = I_e
```

- [ ] **Step 2: 运行测试确认 Solve.jl 正常工作**

```bash
julia --project=. example/testexample.jl
```

---

### Task 4: 在 PostProcessing.jl 中添加热源转换

**Files:**
- Modify: `src/PostProcessing.jl`

- [ ] **Step 1: 添加热源物理单位转换代码**

在 PostProcessing 函数末尾添加：

```julia
    # 热模型后处理：将无量纲热源转换为物理单位
    if case.opt.thermalmodel == "distributed2D" && case.opt.model == "SPMe"
        q_scale = case.param_dim.scale.q

        # 分层热源（转换为物理单位）
        result["thermal2D Q_rxn_NE [W/m3]"] = variables["thermal2D q_rxn_ne"][:, 1:v] .* q_scale
        result["thermal2D Q_rev_NE [W/m3]"] = variables["thermal2D q_rev_ne"][:, 1:v] .* q_scale
        result["thermal2D Q_ohm_s_NE [W/m3]"] = variables["thermal2D q_ohm_s_ne"][:, 1:v] .* q_scale
        result["thermal2D Q_ohm_e_NE [W/m3]"] = variables["thermal2D q_ohm_e_ne"][:, 1:v] .* q_scale
        result["thermal2D Q_SP [W/m3]"] = variables["thermal2D q_sp"][:, 1:v] .* q_scale
        result["thermal2D Q_rxn_PE [W/m3]"] = variables["thermal2D q_rxn_pe"][:, 1:v] .* q_scale
        result["thermal2D Q_rev_PE [W/m3]"] = variables["thermal2D q_rev_pe"][:, 1:v] .* q_scale
        result["thermal2D Q_ohm_s_PE [W/m3]"] = variables["thermal2D q_ohm_s_pe"][:, 1:v] .* q_scale
        result["thermal2D Q_ohm_e_PE [W/m3]"] = variables["thermal2D q_ohm_e_pe"][:, 1:v] .* q_scale
        result["thermal2D Q_PCC [W/m3]"] = variables["thermal2D q_pcc"][:, 1:v] .* q_scale
        result["thermal2D Q_NCC [W/m3]"] = variables["thermal2D q_ncc"][:, 1:v] .* q_scale
    end
```

- [ ] **Step 2: 运行测试确认热源输出正确**

```bash
julia --project=. example/testexample.jl
# 检查 result 中的热源值
```

---

## Chunk 1: Variables.jl 变量重命名

- [ ] 修改 Variables.jl 第98-108行的变量名
- [ ] 搜索其他文件中引用这些变量的地方并更新

- [ ] 运行 `julia --project=. example/testexample.jl` 确认无报错

- [ ] Commit: `git commit -m "refactor(variables): rename heat source variables to dimensionless form"`

## Chunk 2: ThermalDistributed.jl 函数重构
- [ ] 创建 `compute_heat_sources` 函数
- [ ] 创建 `compute_heat_sources_with_czm` 包装函数
- [ ] 删除旧函数 (`heatQ_Source`, `compute_element_heat_sources`, `heatQ_Source_with_czm`)
- [ ] Commit: `git commit -m "refactor(thermal): add unified compute_heat_sources function"`

## Chunk 3: Solve.jl 调用点修改
- [ ] 替换第633-744行为函数调用
- [ ] 保留辅助变量保存逻辑
- [ ] 运行测试确认无报错
- [ ] Commit: `git commit -m "refactor(solve): replace heat source calculation with function call"`

## Chunk 4: PostProcessing.jl 緻加热源转换
- [ ] 添加热源物理单位转换代码
- [ ] 运行测试确认热源输出值正确
- [ ] 运行完整测试套件验证重构前后结果一致
- [ ] Commit: `git commit -m "refactor(postprocessing): add heat source unit conversion"`

## 回滚 Plan

```bash
git checkout -- src/Variables.jl src/ThermalDistributed.jl src/Solve.jl src/PostProcessing.jl
```

