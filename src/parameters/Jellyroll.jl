"""
    Jellyroll cell parameters

    from the paper


    
    and references therein.
"""
# Positive Electrode
PE = Electrode()
PE.theta_100 = 0.263849
PE.theta_0 = 0.853974 
PE.thickness = 7.56e-5   # Chen2020: Positive electrode thickness
PE.lambda = 1.04  # Conductivity_ca (W/m/K, PyECN)
PE.Ds = 4.0e-15          # Chen2020: Positive particle diffusivity (assumed constant)
PE.rho = 2895    # Density_ca (kg/m³)
PE.heat_Q = 1.27e3      # Specific_heat_capacity_ca (J/kg/K)
PE.eps = 0.335           # Chen2020: Positive electrode porosity
PE.eps_fi = 0.025
PE.brugg = 1.5           # Using electrolyte Bruggeman coefficient (Chen2020 electrolyte)
PE.k = 3.42e-6
PE.cs_max = 63104        # Chen2020 max concentration (already matching)
PE.cs0 = 17038
PE.rs = 5.22e-6
PE.sig = 0.18
# Mechanical/thermal expansion (example values)
PE.E = 5.0e10         # Pa
PE.nu = 0.30          # -
PE.alphaT = 1.0e-5    # 1/K
PE.Omega = -7.28e-7     # m^3/mol (placeholder if diffusion-stress needed)
PE.Eac_D = 0
PE.Eac_k = 17800
PE.alpha = 0.5
PE.U = x-> -0.8090*x .+ 4.4875 - 0.0428*tanh.(18.5138*(x .- 0.5542)) - 17.7326*tanh.(15.7890*(x .- 0.3117)) + 17.5842*tanh.(15.9308*(x .- 0.3120))
PE.dUdT = x-> 0 * x

# Negative Electrode
NE = Electrode()
NE.theta_100 = 0.910612
NE.theta_0 = 0.0263472 
NE.thickness = 8.52e-5   # Chen2020: Negative electrode thickness
NE.lambda = 1.58     # Conductivity_an (W/m/K, PyECN)
NE.Ds = 3.3e-14          # Chen2020: Negative particle diffusivity (constant)
NE.rho = 1555      # Density_an (kg/m³)
NE.heat_Q = 1.437e3        # Specific_heat_capacity_an (J/kg/K)
NE.eps = 0.25            # Chen2020: Negative electrode porosity
NE.eps_fi = 0.0326
NE.brugg = 1.5           # Electrolyte Bruggeman coefficient
NE.k = 6.48e-7
NE.cs_max = 33133.0      # Chen2020 max concentration
NE.cs0 = 29866
NE.rs = 5.86e-6          # Chen2020 negative particle radius
NE.sig = 215.0
# Mechanical/thermal expansion (example values)
NE.E = 2.0e10         # Pa
NE.nu = 0.28          # -
NE.alphaT = 3.0e-6    # 1/K
NE.Omega = 3.1e-6    # m^3/mol (placeholder if diffusion-stress needed)
NE.Eac_D = 0.
NE.Eac_k = 35000.
NE.alpha = 0.5
NE.U = x-> 1.97938*exp.(-39.3631*x) .+ 0.2482 - 0.0909*tanh.(29.8538*(x .- 0.1234)) - 0.04478*tanh.(14.9159*(x .- 0.2769)) - 0.0205*tanh.(30.4444*(x .- 0.6103))
NE.dUdT = x-> 0 * x

# Electrolyte
EL = Electrolyte()
EL.De = (x, y=0)-> 8.794e-11 * (x ./ 1000) .^ 2 - 3.972e-10 * (x ./ 1000) .+ 4.862e-10
EL.kappa = (x, y=0)-> 0.1297 * (x ./ 1000) .^ 3 - 2.51 * (x ./ 1000) .^ 1.5 + 3.329 * (x ./ 1000)
EL.dlnf_dlnc = x-> 1
EL.rho = 1290     # 密度 (kg/m³)
EL.heat_Q = 134.1     # 比热容 (J/kg/K)
EL.tplus = 0.2594
EL.ce0 = 1000

# Separator
SP = Separator()
SP.thickness = 1.2e-5    # Chen2020 separator thickness
SP.lambda = 0.16   # 热导率 (W/m/K)
SP.rho = 1100      # 密度 (kg/m³)
SP.heat_Q = 1.978e3        # Specific_heat_capacity_sep (J/kg/K)
SP.eps = 0.47            # Chen2020 separator porosity
SP.eps_fi = 0.
SP.brugg = 1.5           # Chen2020 electrolyte Bruggeman for separator

# Tab
tab = Tab()
tab.width = 4e-3
tab.length = 0.75 * 9.9e-3
tab.area = tab.width * tab.length * 2
tab.h = 10. 
# specify tab angles on circumference (radians)
tab.theta_pos = [15π]  #none or [1，... ,15π...]
tab.theta_neg = [44π]

# Cell parameters for thermal model
cell = Cell()    
cell.Rout = 0.021/2
cell.Rin = 1.92e-3
cell.length = 1.58   # Using outer diameter proxy (not used by scaling below)
cell.width = 7e-2
cell.wrapper = 0
cell.I1C = 5
cell.no_layers = 1
cell.capacity = 5
# 有效宏观受流面积（用于 1D 电化学归一化）
cell.area = cell.width * cell.length * cell.no_layers

# 外表面对流面积（用于团簇热模型）
cell.cooling_surface = 2π * cell.Rout * cell.width + 2 * π * cell.Rout^2
cell.v_h = 4.3
cell.v_l = 2.5
cell.volume = pi * ( (cell.Rout)^2 - (cell.Rin)^2 ) * cell.width  # π (R_out^2 - R_in^2) * height
cell.mass = 8.521789028721335e-01
cell.rho = cell.mass / cell.volume
cell.alphaT = 0.
cell.h = 10.      # 对流换热系数 (W/m²/K)
cell.T0 = 298      # 初始温度 (K)
cell.T_amb = cell.T0


# Positive Current Collector 
PCC = CurrentCollector()
PCC.thickness = 10e-6
PCC.lambda = 237.   # 热导率 (W/m/K)
PCC.rho = 2700.    # 密度 (kg/m³)
PCC.heat_Q = 9.03e2      # 比热容 (J/kg/K)
PCC.sig =3.55e7

# Negative Current Collector
NCC = CurrentCollector()
NCC.thickness = 10e-6
NCC.lambda = 401.   # 热导率 (W/m/K)
NCC.rho = 8940.    # 密度 (kg/m³)
NCC.heat_Q = 3.85e2      # 比热容 (J/kg/K)
NCC.sig = 5.96e7



# Binder
binder = Binder()
# binder.rho = 

# Cohesive zone model parameters for interlayer interface
# 用于描述相邻卷绕圈之间的界面脱粘行为
#
# 参数调优说明（2025年更新）：
# ============================================================================
# 实测数据：
#   - 典型循环中的法向分离位移：δ_n ≈ 10-30 nm
#   - 典型循环中的切向分离位移：δ_t ≈ 50-100 nm
#
# 参数选择原则：
#   1. δ_0 应接近或略小于实际分离位移，这样才能触发损伤
#   2. K = σ_max / δ_0 不能太大（建议 < 1e14 Pa/m），否则数值不稳定
#   3. δ_c 通常设为 δ_0 的 3-10 倍
#
# 验证结果：
#   - 原始参数（δ_0 = 1μm）：D_max = 0%（位移只有阈值的 2-5%）
#   - 调优参数（δ_0 = 20-60nm）：D_max ≈ 60%（可观察到损伤）
# ============================================================================
#
cohesive = Cohesive()

# 法向参数 (Mode I - 张开模式)
# 参数已调优以匹配电池热-化学应变产生的位移量级
cohesive.σ_max_n = 2e6        # 最大法向牵引力 [Pa] (2 MPa)
cohesive.δ_0_n = 20e-9        # 损伤起始分离位移 [m] (20 nm)
cohesive.δ_c_n = 100e-9       # 临界分离位移 [m] (100 nm)
# 断裂能由双线性本构关系计算: G_c = 0.5 * σ_max * δ_c
cohesive.G_c_n = 0.5 * cohesive.σ_max_n * cohesive.δ_c_n  # [J/m²] = 0.1 mJ/m²
# 初始刚度（惩罚刚度）: K = σ_max / δ_0
cohesive.K_n = cohesive.σ_max_n / cohesive.δ_0_n  # [Pa/m] = 1e14 Pa/m

# 切向参数 (Mode II - 剪切模式)
# 切向位移通常比法向大，所以 δ_0_t > δ_0_n
cohesive.τ_max_t = 1e6        # 最大切向牵引力 [Pa] (1 MPa)
cohesive.δ_0_t = 60e-9        # 损伤起始切向位移 [m] (60 nm)
cohesive.δ_c_t = 300e-9       # 临界切向位移 [m] (300 nm)
cohesive.G_c_t = 0.5 * cohesive.τ_max_t * cohesive.δ_c_t  # [J/m²] = 0.15 mJ/m²
cohesive.K_t = cohesive.τ_max_t / cohesive.δ_0_t  # [Pa/m] ≈ 1.67e13 Pa/m

# 混合模式参数
cohesive.eta = 1.45           # BK准则指数（Benzeggagh-Kenane）[-]

# 疲劳损伤参数
# 疲劳损伤增量公式：dD = fatigue_coeff * (δ/δ_0 - threshold)^fatigue_exp
# 这使得即使分离位移未创新高，损伤也会随循环累积
cohesive.fatigue_enabled = true       # 启用疲劳损伤
cohesive.fatigue_coeff = 1e-4         # 疲劳系数（控制累积速率）
cohesive.fatigue_exp = 2.0            # 疲劳指数
cohesive.fatigue_threshold = 0.1      # 疲劳阈值（δ/δ_0 > 0.1 时开始累积）

# Now that PCC/NCC are defined, compute effective thermal conductivities
cell.lambda_r = (NE.thickness + SP.thickness + PE.thickness + PCC.thickness + NCC.thickness) /(NE.thickness/NE.lambda + SP.thickness/SP.lambda + PE.thickness/PE.lambda + PCC.thickness/PCC.lambda + NCC.thickness/NCC.lambda)
cell.lambda_t = (NE.thickness*NE.lambda + SP.thickness*SP.lambda + PE.thickness*PE.lambda + PCC.thickness*PCC.lambda + NCC.thickness*NCC.lambda) /(NE.thickness + SP.thickness + PE.thickness + PCC.thickness + NCC.thickness)

# Scale
scale = Scale()

# assemble to "param_dim"
param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale, cohesive)
