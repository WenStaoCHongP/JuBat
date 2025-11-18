# 多SPMe并行架构修改计划

**目标**: 每个热单元对应一个独立的SPMe子模型，各自维护独立的电化学状态（浓度场），接受分电流和单元温度作为输入。

**日期**: 2025-11-17  
**基于**: 当前代码库（已有分流求解器 `solve_branch_currents_newton`）

---

## 一、核心架构变更

### 1.1 当前架构 vs 目标架构

#### 当前架构（单SPMe + 热场）
```
全局状态向量 yt (共享):
├─ cn_surf[1:Nrn]        # 负极粒子表面浓度
├─ cp_surf[1:Nrp]        # 正极粒子表面浓度
└─ ce[1:Nel]             # 电解液浓度

电化学求解:
SPMe(yt, I_app, T_mean) → M, K, F

热场求解:
ThermalDistributed2D(q_source) → MT, KT, FT

分流求解:
solve_branch_currents_newton → I_e[1:ne]

问题: 所有热单元共享同一个浓度场
```

#### 目标架构（多SPMe并行 + 热场）
```
独立状态向量 yt_e[e] (每单元):
├─ cn_surf_e[1:Nrn]      # 该单元的负极粒子浓度
├─ cp_surf_e[1:Nrp]      # 该单元的正极粒子浓度
└─ ce_e[1:Nel]           # 该单元的电解液浓度

并行电化学求解 (ne 个独立SPMe):
for e in 1:ne
    SPMe_e(yt_e[e], I_e[e], T_e[e]) → M_e, K_e, F_e
end

全局装配:
M_global = blockdiag(M_e[1], M_e[2], ..., M_e[ne], MT)
K_global = blockdiag(K_e[1], K_e[2], ..., K_e[ne], KT)
F_global = [F_e[1]; F_e[2]; ...; F_e[ne]; FT]

热场求解:
ThermalDistributed2D(q_source_e[1:ne]) → MT, KT, FT

分流求解:
solve_branch_currents_newton → I_e[1:ne]

优点: 每个单元独立演化，真正反映空间异质性
```

---

## 二、状态向量重构

### 2.1 状态向量结构

#### 原状态向量（单SPMe）
```julia
# yt: 1D array
yt = [cn_surf[1:Nrn]; cp_surf[1:Nrp]; ce[1:Nel]]
# 总长度: Nrn + Nrp + Nel
```

#### 新状态向量（多SPMe）
```julia
# 方案 A: 平铺向量（推荐）
yt_multi = [
    yt_e[1];  # 单元1的电化学状态 (Nrn + Nrp + Nel)
    yt_e[2];  # 单元2的电化学状态
    ...
    yt_e[ne]; # 单元ne的电化学状态
    T_nodes   # 热场节点温度 (nT)
]
# 总长度: ne * (Nrn + Nrp + Nel) + nT

# 方案 B: 结构化存储（灵活但复杂）
struct MultiSPMeState
    yt_chem::Vector{Vector{Float64}}  # yt_chem[e] 为单元e的状态
    T_nodes::Vector{Float64}           # 热场节点温度
end
```

**推荐**: 方案 A（平铺向量），与现有 Solve 框架兼容。

### 2.2 索引管理

```julia
# SetCase.jl 中新增
function SetMultiSPMeIndex(case::Case)
    ne = size(case.mesh["thermal2D"].element, 1)
    n_chem_per_elem = case.opt.Nrn + case.opt.Nrp + 
                      (case.opt.Nn + case.opt.Ns + case.opt.Np) * case.opt.gsorder
    
    # 为每个单元分配索引范围
    index_multi = Dict{String, Any}()
    for e in 1:ne
        offset = (e-1) * n_chem_per_elem
        index_multi["element $e cn_surf"] = (offset + 1):(offset + case.opt.Nrn)
        index_multi["element $e cp_surf"] = (offset + case.opt.Nrn + 1):(offset + case.opt.Nrn + case.opt.Nrp)
        # ... 其他变量
    end
    
    # 热场索引
    n_chem_total = ne * n_chem_per_elem
    nT = case.mesh["thermal2D"].nlen
    index_multi["T_nodes"] = (n_chem_total + 1):(n_chem_total + nT)
    
    case.index_multi = index_multi
    return case
end
```

---

## 三、核心函数修改

### 3.1 新增：单元级SPMe求解器

```julia
"""
    SPMe_element(case::Case, yt_e::Vector{Float64}, t::Float64, e::Int; 
                 I_e::Float64, T_e::Float64, jacobi::String="update")

为单个热单元求解SPMe模型。

输入:
- yt_e: 该单元的电化学状态向量 (长度: Nrn + Nrp + Nel)
- e: 单元编号
- I_e: 该单元的无量纲电流（由分流求解器提供）
- T_e: 该单元的无量纲温度（由热场提供）

输出:
- M_e: 该单元的质量矩阵
- K_e: 该单元的刚度矩阵
- F_e: 该单元的载荷向量
- variables_e: 该单元的电化学变量字典
"""
function SPMe_element(case::Case, yt_e::Vector{Float64}, t::Float64, e::Int; 
                      I_e::Float64, T_e::Float64, jacobi::String="update")
    param = case.param
    
    # 1) 调用 SPMe_variables，覆写 I_app 和 T_e
    variables_e = SPMe_variables(case, yt_e, t; I_app=I_e, T_e=T_e)
    
    # 2) 力学耦合（如需要）
    if case.opt.mechanicalmodel == "full"
        variables_e = Mechanicaloutput(case, variables_e)
        theta_Mn = variables_e["negative particle stress coupling diffusion coefficient"][1]
        theta_Mp = variables_e["positive particle stress coupling diffusion coefficient"][1]
    else
        theta_Mn = 0.0
        theta_Mp = 0.0
    end
    
    # 3) 粒子扩散
    csn_gs = variables_e["negative particle concentration at gauss point"]
    csp_gs = variables_e["positive particle concentration at gauss point"]
    
    if jacobi == "constant" && param.NE.M_d != []
        M_np = param.NE.M_d
        K_np = param.NE.K_d
        M_pp = param.PE.M_d
        K_pp = param.PE.K_d
    else
        mesh_np = case.mesh["negative particle"]
        mesh_pp = case.mesh["positive particle"]
        M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, theta_Mn)
        M_pp, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, theta_Mp)
    end
    M_np = M_np .* param.scale.ts_n / case.param_dim.scale.t0
    M_pp = M_pp .* param.scale.ts_p / case.param_dim.scale.t0
    
    # 4) 电解液扩散
    mesh_el = case.mesh["electrolyte"]
    M_el, K_el = ElectrolyteDiffusion(param, mesh_el, mesh_el.nlen, variables_e)
    M_el = M_el .* param.scale.te / case.param_dim.scale.t0
    
    # 5) 边界条件
    F_e = SPMe_BC(case, variables_e)
    
    # 6) 装配该单元的系统矩阵
    M_e = blockdiag(M_np, M_pp, M_el)
    K_e = blockdiag(K_np, K_pp, K_el)
    
    return M_e, K_e, F_e, variables_e
end
```

**关键点**:
- `yt_e` 是该单元的**局部**状态向量，不是全局向量
- 利用现有的 `SPMe_variables` 函数（已支持覆写 `I_app` 和 `T_e`）
- 返回该单元的局部 M/K/F

---

### 3.2 修改：CallModel 支持多SPMe模式

```julia
function CallModel(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    # 判断是否启用多SPMe模式
    multi_spme_enabled = (
        case.opt.model == "SPMe" &&
        hasproperty(case.opt, :per_element_spme) && case.opt.per_element_spme &&
        case.opt.thermalmodel == "distributed2D" &&
        haskey(case.mesh, "thermal2D")
    )
    
    if multi_spme_enabled
        return CallModel_MultiSPMe(case, yt, t, jacobi=jacobi)
    else
        # 现有逻辑（单SPMe）
        if case.opt.model == "SPM"
            M, K, F, variables = SPM(case, yt, t, jacobi=jacobi)
            y_phi = Float64[]
        elseif case.opt.model == "SPMe"
            M, K, F, variables = SPMe(case, yt, t, jacobi=jacobi)
            y_phi = Float64[]
        # ... 其他模型
        end
        
        # 热装配（现有逻辑）
        if case.opt.thermalmodel == "distributed2D"
            # ... 现有代码
        end
        
        return M, K, F, variables, y_phi
    end
end
```

---

### 3.3 新增：多SPMe模式的 CallModel

```julia
"""
    CallModel_MultiSPMe(case, yt, t; jacobi)

多SPMe并行模式的 CallModel。

状态向量结构:
yt = [yt_e[1]; yt_e[2]; ...; yt_e[ne]; T_nodes]

其中 yt_e[e] 包含单元e的电化学状态（粒子+电解液浓度）。
"""
function CallModel_MultiSPMe(case::Case, yt::Array{Float64}, t::Float64; jacobi::String)
    mesh_th = case.mesh["thermal2D"]
    ne = size(mesh_th.element, 1)
    nT = mesh_th.nlen
    
    # 1) 解析状态向量
    n_chem_per_elem = length(case.index["negative particle surface lithium concentration"]) +
                      length(case.index["positive particle surface lithium concentration"]) +
                      length(case.index["electrolyte lithium concentration"])
    
    n_chem_total = ne * n_chem_per_elem
    
    # 提取热场
    T_nodes = yt[(n_chem_total + 1):(n_chem_total + nT)]
    
    # 提取每个单元的电化学状态
    yt_chem = Vector{Vector{Float64}}(undef, ne)
    for e in 1:ne
        offset = (e-1) * n_chem_per_elem
        yt_chem[e] = yt[(offset + 1):(offset + n_chem_per_elem)]
    end
    
    # 2) 计算元素均温和面积
    areas = if haskey(case.mesh, "thermal2D element area cache")
        case.mesh["thermal2D element area cache"]
    else
        A = zeros(Float64, ne)
        ngs = length(mesh_th.gs.detJ)
        @inbounds for g in 1:ngs
            e = mesh_th.gs.ele[g]
            A[e] += mesh_th.gs.weight[g] * mesh_th.gs.detJ[g]
        end
        case.mesh["thermal2D element area cache"] = A
        A
    end
    
    Te_prev = zeros(Float64, ne)
    @inbounds for e in 1:ne
        nds = mesh_th.element[e, :]
        Te_prev[e] = sum(T_nodes[nds]) / length(nds)
    end
    
    # 3) 分流求解（获取 I_e）
    I_total = case.opt.Current(t * case.param.scale.t0) / case.param_dim.cell.I1C
    
    # 初始化 variables（用于传递给分流求解器）
    variables = StandardVariables(case, 1)
    variables["cell current"] = I_total
    variables["T_nodes"] = T_nodes
    variables["thermal2D element area"] = areas
    
    # 如果是首次调用，I_e_prev 不存在，传 nothing
    I_e_prev = haskey(case, :I_e_cache) ? case.I_e_cache : nothing
    
    variables, I_e, Vc = solve_branch_currents_newton(
        case, variables, yt, t, I_total, areas, Te_prev, I_e_prev
    )
    
    # 缓存 I_e 供下一步使用
    case.I_e_cache = copy(I_e)
    
    # 4) 并行求解每个单元的SPMe
    M_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    K_elems = Vector{SparseMatrixCSC{Float64,Int64}}(undef, ne)
    F_elems = Vector{Vector{Float64}}(undef, ne)
    variables_elems = Vector{Dict{String,Union{Array{Float64},Float64}}}(undef, ne)
    
    # 可并行化（如果Julia支持多线程）
    # Threads.@threads for e in 1:ne
    for e in 1:ne
        M_e, K_e, F_e, vars_e = SPMe_element(
            case, yt_chem[e], t, e;
            I_e = I_e[e],
            T_e = Te_prev[e],
            jacobi = jacobi
        )
        M_elems[e] = sparse(M_e)
        K_elems[e] = sparse(K_e)
        F_elems[e] = F_e
        variables_elems[e] = vars_e
    end
    
    # 5) 装配电化学全局矩阵
    M_chem = blockdiag(M_elems...)
    K_chem = blockdiag(K_elems...)
    F_chem = vcat(F_elems...)
    
    # 6) 计算逐单元热源
    q_elem = zeros(Float64, ne)
    eta_n_e = zeros(Float64, ne)
    eta_p_e = zeros(Float64, ne)
    dUdT_n_e = zeros(Float64, ne)
    dUdT_p_e = zeros(Float64, ne)
    
    for e in 1:ne
        vars_e = variables_elems[e]
        eta_n_e[e] = vars_e["negative electrode overpotential"][1]
        eta_p_e[e] = vars_e["positive electrode overpotential"][end]
        
        # dUdT（从单元变量中提取）
        cn_surf_e = vars_e["negative particle surface lithium concentration"][1]
        cp_surf_e = vars_e["positive particle surface lithium concentration"][end]
        dUdT_n_e[e] = case.param.NE.dUdT(cn_surf_e)[1]
        dUdT_p_e[e] = case.param.PE.dUdT(cp_surf_e)[1]
        
        # 热源计算（逐单元）
        T_e = Te_prev[e]
        I_e_local = I_e[e]
        
        # 反应热
        Q_rxn = abs(I_e_local * (eta_p_e[e] - eta_n_e[e]))
        
        # 可逆热
        Q_rev = abs(I_e_local) * T_e * (dUdT_p_e[e] - dUdT_n_e[e])
        
        # 欧姆热（沿用现有逻辑）
        param = case.param
        L_th = case.param_dim.scale.L_th
        kappa_ne = param.EL.kappa(param.EL.ce0, T_e) * param.NE.eps ^ param.NE.brugg
        kappa_pe = param.EL.kappa(param.EL.ce0, T_e) * param.PE.eps ^ param.PE.brugg
        kappa_sp = param.EL.kappa(param.EL.ce0, T_e) * param.SP.eps ^ param.SP.brugg
        sig_n_eff = param.NE.sig * param.NE.eps_s
        sig_p_eff = param.PE.sig * param.PE.eps_s
        
        t_n = param.NE.thickness / L_th
        t_p = param.PE.thickness / L_th
        t_sp = param.SP.thickness / L_th
        
        P_s_ne = I_e_local^2 * (t_n / sig_n_eff) / 3.0
        P_s_pe = I_e_local^2 * (t_p / sig_p_eff) / 3.0
        P_e_ne = I_e_local^2 * (t_n / kappa_ne) / 3.0
        P_e_sp = I_e_local^2 * (t_sp / kappa_sp)
        P_e_pe = I_e_local^2 * (t_p / kappa_pe) / 3.0
        
        Q_ohm = P_e_ne / t_n + P_s_ne / t_n + P_e_pe / t_p + P_s_pe / t_p + P_e_sp / t_sp
        
        Q_ele = Q_rxn + Q_rev + Q_ohm
        
        # 集流体欧姆热（使用 layer_weights）
        fks = haskey(variables, "thermal2D layer_weights") ? variables["thermal2D layer_weights"] : nothing
        if fks !== nothing
            σ_PCC = max(hasproperty(param, :PCC) && hasproperty(param.PCC, :sig) ? param.PCC.sig : 1e12, 1e-12)
            σ_NCC = max(hasproperty(param, :NCC) && hasproperty(param.NCC, :sig) ? param.NCC.sig : 1e12, 1e-12)
            Q_PCC = I_e_local^2 / (3.0 * σ_PCC)
            Q_NCC = I_e_local^2 / (3.0 * σ_NCC)
            q_elem[e] = (fks[e,1] + fks[e,2] + fks[e,3]) * Q_ele + fks[e,4] * Q_PCC + fks[e,5] * Q_NCC
        else
            q_elem[e] = Q_ele
        end
    end
    
    # 写入 variables（用于热装配）
    variables["thermal2D element current"] = I_e
    variables["thermal2D eta_n_e"] = eta_n_e
    variables["thermal2D eta_p_e"] = eta_p_e
    variables["thermal2D dUdT_n_e"] = dUdT_n_e
    variables["thermal2D dUdT_p_e"] = dUdT_p_e
    
    # 无量纲化热源
    if hasproperty(case.opt, :units_thermal) && case.opt.units_thermal == "SI"
        variables["heat_source_fields"] = q_elem
        variables["heat_source_units_code"] = 1.0
    else
        q_ref = case.param_dim.scale.q_th
        variables["heat_source_fields"] = q_elem ./ q_ref
        variables["heat_source_units_code"] = 0.0
    end
    
    # 7) 装配热学矩阵
    MT, KT, FT = ThermalDistributed2D(case, variables)
    t_ratio = case.param_dim.scale.t0 / case.param_dim.scale.t_th
    MT = MT .* t_ratio
    ThermalDistributed2D_BC(KT, FT, case, t)
    
    # 8) 全局拼装
    M = blockdiag(M_chem, sparse(MT))
    K = blockdiag(K_chem, sparse(KT))
    F = [F_chem; FT]
    
    # 9) 合并 variables（保留关键全局信息）
    # 电压取平均或用公共电压 Vc
    variables["cell voltage"] = Vc
    variables["time"] = t
    variables["temperature"] = mean(T_nodes)  # 平均温度
    
    y_phi = Float64[]
    
    return M, K, F, variables, y_phi
end
```

**关键点**:
- 解析平铺状态向量，提取每个单元的电化学状态
- 并行调用 `SPMe_element` 求解每个单元
- 使用每个单元的局部 η 和 dUdT 计算热源
- 装配全局系统矩阵

---

### 3.4 修改：Solve 主循环初始化

```julia
function Solve(case::Case)
    # ... 现有代码
    
    # 初始化
    multi_spme_enabled = (
        case.opt.model == "SPMe" &&
        hasproperty(case.opt, :per_element_spme) && case.opt.per_element_spme &&
        case.opt.thermalmodel == "distributed2D" &&
        haskey(case.mesh, "thermal2D")
    )
    
    if multi_spme_enabled
        # 多SPMe模式：构造扩展状态向量
        y0 = ModelInitialisation_MultiSPMe(case)
    else
        # 单SPMe模式（现有）
        y0 = ModelInitialisation(case)
    end
    
    # ... 现有时间推进逻辑（无需修改）
end
```

---

### 3.5 新增：多SPMe初始化

```julia
"""
    ModelInitialisation_MultiSPMe(case)

为多SPMe模式初始化状态向量。

返回:
y0 = [yt_e[1]; yt_e[2]; ...; yt_e[ne]; T_nodes]
"""
function ModelInitialisation_MultiSPMe(case::Case)
    # 1) 单个SPMe单元的初始化（使用现有函数）
    y0_single = ModelInitialisation(case)
    
    # 2) 复制到每个单元
    ne = size(case.mesh["thermal2D"].element, 1)
    n_chem = length(y0_single)
    
    y0_chem = repeat(y0_single, ne)  # 简单复制
    
    # 3) 添加热场初始温度
    nT = case.mesh["thermal2D"].nlen
    T0 = case.param.cell.T0
    T_nodes_init = fill(T0, nT)
    
    y0 = [y0_chem; T_nodes_init]
    
    return y0
end
```

**注意**: 简单复制可能不够，如果不同位置初始SOC不同，需要额外处理。

---

## 四、热源计算逻辑（已在 CallModel_MultiSPMe 中实现）

热源计算已集成在 `CallModel_MultiSPMe` 的第6步，使用每个单元的局部变量：

```julia
for e in 1:ne
    # 从该单元的 variables_elems[e] 提取
    eta_n_e[e] = vars_e["negative electrode overpotential"][1]
    eta_p_e[e] = vars_e["positive electrode overpotential"][end]
    dUdT_n_e[e] = ...
    dUdT_p_e[e] = ...
    
    # 逐单元热源
    Q_rxn = abs(I_e[e] * (eta_p_e[e] - eta_n_e[e]))
    Q_rev = abs(I_e[e]) * T_e[e] * (dUdT_p_e[e] - dUdT_n_e[e])
    Q_ohm = ...
end
```

**不再需要**修改 `heatQ_Source` 函数，因为热源在 `CallModel_MultiSPMe` 内部直接计算。

---

## 五、配置选项

### 5.1 Option.jl 新增字段

```julia
@with_kw mutable struct Option
    # ... 现有字段
    
    # 多SPMe模式开关
    per_element_spme::Bool = false     # true: 启用多SPMe模式
    
    # 调试选项
    debug_multi_spme::Bool = false     # 打印每个单元的详细信息
    debug_sample_elems::Int = 3        # 采样打印的单元数
end
```

**启用多SPMe模式**:
```julia
case.opt.per_element_spme = true
case.opt.thermalmodel = "distributed2D"
```

---

## 六、实施步骤

### 阶段 1: 单元级SPMe求解器（1-2天）
- [ ] 实现 `SPMe_element` 函数
- [ ] 测试：单个单元求解，验证与全局SPMe结果一致
- [ ] 单元测试：不同 I_e 和 T_e 输入

### 阶段 2: 多SPMe初始化与索引（1天）
- [ ] 实现 `ModelInitialisation_MultiSPMe`
- [ ] 实现 `SetMultiSPMeIndex`（可选，取决于是否需要按名称访问）
- [ ] 测试：初始化后状态向量结构正确

### 阶段 3: CallModel_MultiSPMe 实现（2-3天）
- [ ] 实现 `CallModel_MultiSPMe` 主函数
- [ ] 状态向量解析与重组
- [ ] 并行调用 `SPMe_element`
- [ ] 逐单元热源计算
- [ ] 全局装配
- [ ] 测试：零电流工况、恒流工况

### 阶段 4: CallModel 集成（半天）
- [ ] 修改 `CallModel` 添加多SPMe模式分支
- [ ] 修改 `Solve` 主循环初始化逻辑
- [ ] 测试：开关切换，确保单SPMe模式不受影响

### 阶段 5: 性能优化（1天，可选）
- [ ] 缓存温度无关的预因子（避免重复计算）
- [ ] 并行化 SPMe_element 调用（`Threads.@threads`）
- [ ] 稀疏矩阵装配优化

### 阶段 6: 验证与测试（2-3天）
- [ ] 能量守恒验证（每个单元独立守恒）
- [ ] 电流守恒验证（`sum(w .* I_e) ≈ I_total`）
- [ ] 温度分布合理性
- [ ] 浓度场空间分布（不同单元应有差异）
- [ ] 与单SPMe模式对比（均匀温度下应一致）

**总工作量**: 约 1-2 周（1人全职）

---

## 七、关键难点与解决方案

### 7.1 状态向量维度爆炸

**问题**: 状态向量从 `O(100)` 增长到 `O(100 * ne)`，对于 `ne=1000`，可能达到 `O(10^5)`。

**解决方案**:
1. **稀疏矩阵高效装配**: Julia 的 `blockdiag` 保持稀疏性
2. **并行化**: 利用 `Threads.@threads` 并行求解每个单元
3. **预条件器**: 如收敛慢，考虑 ILU 预条件
4. **降阶**: 如果 `Nrn=Nrp=1`（SPMe 简化模型），单元状态向量很小

**估算**（以 LG M50 参数为例）:
- 单SPMe: Nrn=10, Nrp=10, Nel≈40，总约60个自由度
- 多SPMe (ne=100): 60 * 100 = 6000 电化学自由度 + 约200 热自由度 = 6200
- 多SPMe (ne=1000): 60 * 1000 = 60000 电化学自由度 + 约2000 热自由度 = 62000

对于现代计算机，6万自由度的稀疏线性系统求解是可行的（约0.1-1秒/步）。

---

### 7.2 并行求解的数据竞争

**问题**: 多个线程同时调用 `SPMe_element` 可能访问共享的 `case.param`。

**解决方案**:
1. **只读参数安全**: `case.param` 中的材料参数都是只读的，线程安全
2. **局部变量**: `variables_e` 是每个线程的局部变量，无竞争
3. **预分配**: 提前分配 `M_elems`, `K_elems`, `F_elems`，每个线程写入不同位置

**代码示例**:
```julia
# 安全的并行化
Threads.@threads for e in 1:ne
    M_e, K_e, F_e, vars_e = SPMe_element(
        case, yt_chem[e], t, e;
        I_e = I_e[e],
        T_e = Te_prev[e],
        jacobi = "update"  # 避免共享 param.NE.M_d
    )
    M_elems[e] = sparse(M_e)
    K_elems[e] = sparse(K_e)
    F_elems[e] = F_e
    variables_elems[e] = vars_e
end
```

**注意**: 设置 `jacobi="update"` 避免多线程同时修改 `case.param.NE.M_d`。

---

### 7.3 分流求解器的适配

**问题**: 当前 `solve_branch_currents_newton` 假设单SPMe，如何适配多SPMe？

**解决方案**: **无需修改**。分流求解器的职责是：
- 输入: 每个单元的温度 `T_e`（已有）
- 输出: 每个单元的电流 `I_e`（已有）

在多SPMe模式下，分流求解器仍使用相同的逻辑，只是：
- 在 `CallModel_MultiSPMe` 中调用它
- 传递的 `variables` 是临时的（仅用于传递参数）
- 返回的 `I_e` 用于并行调用 `SPMe_element`

**关键点**: 分流求解器只需要全局浓度场的**代表值**（用于计算 prefactor），这可以通过取所有单元的平均或使用某个代表单元实现。

**修正建议**（在 `CallModel_MultiSPMe` 中）:
```julia
# 构造代表性的全局状态（用于分流求解器的 prefactor 计算）
yt_representative = mean(yt_chem)  # 或取某个中心单元

# 分流求解器使用代表状态
variables, I_e, Vc = solve_branch_currents_newton(
    case, variables, yt_representative, t, I_total, areas, Te_prev, I_e_prev
)
```

---

### 7.4 浓度场的空间耦合

**问题**: 不同单元的电解液浓度是否应耦合？（通过电解液的横向扩散）

**当前假设**: **无横向扩散**，每个单元的电解液浓度独立演化。

**物理合理性**:
- ✅ 合理：如果电池是**螺旋卷绕**或**叠片**结构，层间电解液扩散很慢（时间尺度远大于电化学反应）
- ❌ 不合理：如果电池是**单层大面积电极**，横向扩散不可忽略

**未来扩展**: 如需考虑横向扩散，需要：
1. 修改电解液扩散方程，添加横向（x-y 方向）项
2. 在相邻单元间添加浓度梯度的通量项
3. 这将大幅增加复杂度，建议作为 Phase 2

---

## 八、测试与验证

### 8.1 单元测试

#### 测试 1: 单元SPMe求解器
```julia
# test/test_spme_element.jl
@testset "SPMe_element" begin
    case = SetCase_SPMe_thermal2D()
    yt_e = ModelInitialisation(case)
    
    # 测试不同电流
    for I_e in [0.0, 0.5, 1.0, 2.0]
        M_e, K_e, F_e, vars_e = SPMe_element(case, yt_e, 0.0, 1; I_e=I_e, T_e=1.0)
        
        @test size(M_e, 1) == length(yt_e)
        @test vars_e["negative electrode interfacial current density"] ≈ I_e / ... atol=1e-6
    end
end
```

#### 测试 2: 多SPMe初始化
```julia
@testset "MultiSPMe_init" begin
    case = SetCase_SPMe_thermal2D()
    case.opt.per_element_spme = true
    
    y0 = ModelInitialisation_MultiSPMe(case)
    
    ne = size(case.mesh["thermal2D"].element, 1)
    n_chem = length(ModelInitialisation(case))
    nT = case.mesh["thermal2D"].nlen
    
    @test length(y0) == ne * n_chem + nT
end
```

---

### 8.2 集成测试

#### 测试 3: 均匀温度下与单SPMe对比
```julia
@testset "MultiSPMe vs SingleSPMe uniform T" begin
    # 单SPMe
    case1 = SetCase_SPMe_thermal2D()
    case1.opt.per_element_spme = false
    result1 = Solve(case1)
    
    # 多SPMe
    case2 = deepcopy(case1)
    case2.opt.per_element_spme = true
    result2 = Solve(case2)
    
    # 在均匀温度下，结果应近似一致
    @test result1["voltage"][end] ≈ result2["voltage"][end] atol=1e-3
end
```

#### 测试 4: 能量守恒
```julia
@testset "Energy conservation" begin
    case = SetCase_SPMe_thermal2D()
    case.opt.per_element_spme = true
    case.opt.debug_coupling = true
    
    result = Solve(case)
    
    # 检查能量残差（从日志中提取）
    # 或在 CallModel_MultiSPMe 中添加验证逻辑
end
```

---

### 8.3 性能基准

```julia
# benchmark/benchmark_multi_spme.jl
using BenchmarkTools

function benchmark_solve(ne)
    case = SetCase_SPMe_thermal2D()
    # 调整网格使得 ne ≈ 指定值
    case.opt.per_element_spme = true
    
    @btime Solve($case)
end

# 测试不同单元数
for ne in [10, 50, 100, 500, 1000]
    println("ne = $ne")
    benchmark_solve(ne)
end
```

**预期性能**（粗略估计）:
- ne=100: 约 10-60 秒（取决于时间步数）
- ne=1000: 约 100-600 秒

**优化目标**: 通过并行化将 ne=1000 的时间降至 50-100 秒。

---

## 九、与原修改计划的对比

| 项目 | 原修改计划 | 多SPMe架构 | 变化 |
|-----|----------|----------|------|
| 浓度场 | 全局共享 | 每单元独立 | ✅ 核心变化 |
| 状态向量 | `yt[1:n]` | `yt[1:ne*n+nT]` | ✅ 大幅扩展 |
| SPMe求解 | 单次调用 | ne次并行调用 | ✅ 新增 |
| 热源计算 | 全局η → 逐单元η | 逐单元η（从局部yt_e计算） | ✅ 更精确 |
| 分流求解 | 已实现 | 复用（需适配代表状态） | ⚠️ 微调 |
| CallModel | 单路径 | 双路径（单/多SPMe） | ✅ 新增分支 |
| 装配方式 | CallModel内 | CallModel内 | ✅ 保持 |

**关键优势**:
1. 每个单元的浓度场独立演化，真正反映空间异质性
2. 逐单元热源完全精确（不需要批处理函数的近似）
3. 易于并行化，性能可扩展

**代价**:
1. 状态向量维度增加 `ne` 倍
2. 计算量增加（但可通过并行化缓解）
3. 内存占用增加

---

## 十、实施优先级建议

### 立即实施（核心功能，2周）:
1. ✅ `SPMe_element` 函数（阶段1）
2. ✅ `ModelInitialisation_MultiSPMe`（阶段2）
3. ✅ `CallModel_MultiSPMe`（阶段3）
4. ✅ `CallModel` 集成（阶段4）
5. ✅ 基础测试（阶段6部分）

### 近期实施（优化与验证，1周）:
6. ⏭ 并行化优化（阶段5）
7. ⏭ 完整验证测试（阶段6）
8. ⏭ 性能基准测试

### 后续实施（扩展功能）:
9. ⏭ 考虑横向电解液扩散（物理扩展）
10. ⏭ 非均匀初始SOC分布
11. ⏭ 多SPMe模式的后处理可视化

---

## 十一、代码文件清单

| 文件 | 修改类型 | 描述 |
|-----|---------|------|
| `src/Option.jl` | 修改 | 添加 `per_element_spme` 字段 |
| `src/SPMe.jl` | 新增函数 | `SPMe_element` |
| `src/Solve.jl` | 新增函数 + 修改 | `CallModel_MultiSPMe`, `ModelInitialisation_MultiSPMe`, 修改 `CallModel` |
| `src/Initialisation.jl` | 新增函数 | `ModelInitialisation_MultiSPMe` |
| `test/test_multi_spme.jl` | 新增 | 单元测试与集成测试 |
| `benchmark/benchmark_multi_spme.jl` | 新增 | 性能基准测试 |

---

## 十二、总结

多SPMe并行架构是实现"每个热单元对应独立SPMe模型"目标的**正确方向**。相比原修改计划中的"批处理逐单元变量"方案，多SPMe架构：

### 优势：
✅ **物理精确性**: 每个单元的浓度场独立演化，真实反映局部电化学状态  
✅ **热源精确性**: 逐单元η/dUdT直接从局部状态计算，无需近似  
✅ **可扩展性**: 易于并行化，支持大规模单元数  
✅ **模块化**: 单元级求解器独立，易于测试和维护

### 挑战：
⚠️ **计算量**: 状态向量和计算量增加 ne 倍  
⚠️ **内存**: 需要存储 ne 个独立状态  
⚠️ **实施复杂度**: 需要重构 CallModel 和 Solve 逻辑

### 建议：
建议按照本修改计划分阶段实施，优先完成核心功能（阶段1-4），然后进行优化（阶段5）和验证（阶段6）。

**预计工作量**: 1-2周（1人全职）完成核心功能，额外1周进行优化与验证。

---

**附录**: 本修改计划可与《逐单元改进算法修改计划.md》对比，后者适用于"共享浓度场 + 逐单元变量批处理"的场景，而本计划适用于"独立浓度场 + 多SPMe并行"的场景。根据实际物理需求和计算资源选择合适方案。
