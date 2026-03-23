# 归一化实现检查计划

## 目标
检查 src 代码中是否完全实现：
1. 归一化在 `NormaliseParam` 中统一处理
2. 输出在 `PostProcessing` 中还原单位
3. 识别代码中存在的量纲混乱问题

## 检查阶段

### Phase 1: 阅读关键源文件 ✓
- [x] SetParams.jl - NormaliseParam 实现
- [x] PostProcessing.jl - 输出还原实现
- [x] ThermalDistributed.jl - 热模型归一化使用
- [x] Solve.jl - 主求解器归一化处理
- [x] Materialmatrix.jl - 材料矩阵实现
- [x] SPMe.jl - 电化学模型归一化

### Phase 2: 检查归一化实现 ✓

#### 2.1 NormaliseParam (SetParams.jl) - 热参数归一化

| 参数 | 文档公式 | 代码实现 | 状态 |
|------|----------|----------|------|
| **密度** ρ* | ρ / ρ_ref | `param.PE.rho = param_dim.PE.rho / scale.rho` (L336) | ✓ 正确 |
| **比热容** c* | c × ρ_ref × L³ × T_ref / (t₀ × P_ref) | `param.PE.heat_Q = ... * scale.rho * L^3 * T_ref / (t0 * phi * I_typ)` (L337) | ✓ 正确 |
| **导热率** k* | k × L × T_ref / P_ref | `param.PE.lambda = param_dim.PE.lambda / scale.lambda` (L335) | ✓ 正确 |
| **传热系数** h* | h × A × T_ref / P_ref | `param.cell.h = h * A * T_ref / (phi * I_typ)` (L394) | ✓ 正确 |
| **总热容** C* | C × m × T_ref / (t₀ × P_ref) | `param.cell.heat_Q = Q * m * T_ref / (t0 * phi * I_typ)` (L397) | ✓ 正确 |
| **热源尺度** q_ref | P_ref / L³ | `scale.q = P_ref / L^3` (L297) | ✓ 正确 |

#### 2.2 CZM 参数归一化 (SetParams.jl L415-436)

| 参数 | 归一化公式 | 状态 |
|------|------------|------|
| 界面换热系数 h_c0* | h_c0 × L / λ_ref | ✓ 正确 |
| 空气热导率 k_air* | k_air / λ_ref | ✓ 正确 |
| 平均自由程 λ_m* | λ_m / L | ✓ 正确 |

### Phase 3: 检查输出还原实现 ✓

#### PostProcessing.jl 输出还原

| 变量 | 还原公式 | 代码行 | 状态 |
|------|----------|--------|------|
| 时间 | t = t* × t₀ | L3 | ✓ 正确 |
| 电压 | V = V* × φ_ref | L4 | ✓ 正确 |
| 电流 | I = I* × I_1C | L5 | ✓ 正确 |
| 温度 | T = T* × T_ref | L6 | ✓ 正确 |
| 应力 | σ = σ* × E_ref | L7-10 | ✓ 正确 |
| 位移 | u = u* × r₀ | L11-12 | ✓ 正确 |
| 浓度 | c = c* × c_max | L14-17 | ✓ 正确 |

### Phase 4: 发现的问题

## 状态: completed

## 详细发现报告

见 `docs/check_normalization_findings.md`
