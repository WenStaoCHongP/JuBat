"""
    Jellyroll cell parameters

    from the paper


    
    and references therein.
"""
# Positive Electrode
PE = Electrode()
PE.theta_100 = 0.263849
PE.theta_0 = 0.853974 
PE.thickness = 75.6e-6   # Chen2020: Positive electrode thickness
PE.lambda = 0.892  # Conductivity_ca (W/m/K, PE.lambda = 1.58)
PE.Ds = 4.0e-15          # Chen2020: Positive particle diffusivity (assumed constant)
PE.rho = 3216    # Density_ca (kg/m³)
PE.heat_Q = 983      # Specific_heat_capacity_ca (J/kg/K)
PE.eps = 0.335           # Chen2020: Positive electrode porosity
PE.eps_fi = 0
PE.brugg = 1.5           # Using electrolyte Bruggeman coefficient (Chen2020 electrolyte)
PE.k = 3.42e-6
PE.cs_max = 63104        # Chen2020 max concentration (already matching)
PE.cs0 = 17038
PE.rs = 5.22e-6
PE.sig = 0.18
# Mechanical/thermal expansion (example values)
PE.E = 3.75e10         # Pa颗粒的弹性模量
PE.nu = 0.20          # -颗粒的泊松比
PE.alphaT = 1.5e-5    # 1/K
PE.E_coat = 1e9         # Pa 正极片弹性模量 (500 MPa)
PE.nu_coat = 0.3          # 正极片泊松比
PE.Omega = 7.88e-7     # m^3/mol 正极偏摩尔体积（嵌锂膨胀；2026-08-29 由负值更正）
PE.Eac_D = 0
PE.Eac_k = 17800
PE.alpha = 0.5
PE.U = x-> -0.8090*x .+ 4.4875 - 0.0428*tanh.(18.5138*(x .- 0.5542)) - 17.7326*tanh.(15.7890*(x .- 0.3117)) + 17.5842*tanh.(15.9308*(x .- 0.3120))
PE.dUdT = x-> 0 * x

# Negative Electrode
NE = Electrode()
NE.theta_100 = 0.910612
NE.theta_0 = 0.0263472 
NE.thickness = 85.2e-6   # Chen2020: Negative electrode thickness
NE.lambda = 4.058     # Conductivity_an (W/m/K, PyECN)
NE.Ds = 3.3e-14          # Chen2020: Negative particle diffusivity (constant)
NE.rho = 2213      # Density_an (kg/m³)
NE.heat_Q = 809        # Specific_heat_capacity_an (J/kg/K)
NE.eps = 0.25            # Chen2020: Negative electrode porosity
NE.eps_fi = 0
NE.brugg = 1.5           # Electrolyte Bruggeman coefficient
NE.k = 6.48e-7
NE.cs_max = 33133.0      # Chen2020 max concentration
NE.cs0 = 29866
NE.rs = 5.86e-6          # Chen2020 negative particle radius
NE.sig = 215.0
# Mechanical/thermal expansion (example values)
NE.E = 1.5e10         # Pa颗粒的弹性模量
NE.nu = 0.28          # -颗粒的泊松比
NE.alphaT = 8.0e-6    # 1/K
NE.E_coat = 1e9         # Pa 负极片弹性模量 (500 MPa)
NE.nu_coat = 0.25          # 负极片泊松比
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
EL.heat_Q = 229     # 比热容 (J/kg/K)
EL.tplus = 0.2594
EL.ce0 = 1000

# Separator
SP = Separator()
SP.thickness = 1.2e-5    # Chen2020 separator thickness
SP.lambda = 0.3344   # 热导率 (W/m/K)
SP.rho = 1548      # 密度 (kg/m³)
SP.heat_Q = 1128        # Specific_heat_capacity_sep (J/kg/K)
SP.eps = 0.47            # Chen2020 separator porosity
SP.eps_fi = 0.
SP.brugg = 1.5           # Chen2020 electrolyte Bruggeman for separator
SP.E = 750e6         # Pa 隔膜的弹性模量
SP.nu = 0.35          # 隔膜的泊松比
SP.alphaT = 30e-6    # 1/K 隔膜热膨胀系数（2026-08-31 用户指定，解除原显式置零）

# Positive Current Collector 
PCC = CurrentCollector()
PCC.thickness = 16e-6
PCC.lambda = 237.   # 热导率 (W/m/K)
PCC.rho = 2702.    # 密度 (kg/m³)
PCC.heat_Q = 8.76e2      # 比热容 (J/kg/K)
PCC.sig =3.55e7
PCC.E = 7e10         # Pa 正极集流体的弹性模量（软等效层，默认路径）
PCC.nu = 0.33          # 正极集流体的泊松比
PCC.alphaT = 23e-6   # 1/K Al 箔热膨胀系数（2026-08-31 用户指定，解除原显式置零）
# Batch 3：Al 箔塑性参数（本构模量直接用上面的 PCC.E/nu，用户参数修正后不再单列 foil 字段）
PCC.sigma_y = 60e6    # Pa（Shi 2026 屈服起始）
PCC.H = 0.0           # Pa 理想塑性（敏感性扫描 0–2 GPa）

# Negative Current Collector
NCC = CurrentCollector()
NCC.thickness = 12.00e-6
NCC.lambda = 401.   # 热导率 (W/m/K)
NCC.rho = 8933.    # 密度 (kg/m³)
NCC.heat_Q = 3.83e2      # 比热容 (J/kg/K)
NCC.sig = 5.96e7
NCC.E = 1.1e11         # Pa 负极集流体的弹性模量（软等效层，默认路径）
NCC.nu = 0.34         # 负极集流体的泊松比
NCC.alphaT = 17e-6   # 1/K Cu 箔热膨胀系数（2026-08-31 用户指定，解除原显式置零）
# Batch 3：Cu 箔塑性参数（本构模量直接用上面的 NCC.E/nu）
NCC.sigma_y = 200e6   # Pa（ED Cu 文献区间 108–441 MPa 中值）
NCC.H = 0.0           # Pa 理想塑性

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
cell.Rout = 0.0203/2
cell.Rin = 1.92e-3
cell.winding_T_ne = 3.0e6   # Pa 负极侧卷绕张力（文献整理 §10.4.1：1–5 MPa 中值）
cell.winding_T_pe = 1.0e6   # Pa 正极侧卷绕张力（0.5–2 MPa 中值）
cell.length = 1.58   # Using outer diameter proxy (not used by scaling below)
cell.width = 6.5e-2
cell.layer = 2 * (PE.thickness + NE.thickness + SP.thickness) + PCC.thickness + NCC.thickness
cell.wrapper = 0
cell.I1C = 5
cell.no_layers = 1
cell.capacity = 5
# 有效宏观受流面积（用于 1D 电化学归一化）
cell.area = cell.width * cell.length * cell.no_layers

# 外表面对流面积（用于团簇热模型）
cell.cooling_surface = 2π * cell.Rout * cell.width + 2 * π * (cell.Rout^2- cell.Rin^2)
cell.v_h = 4.3
cell.v_l = 2.5
cell.volume = pi * ( (cell.Rout)^2 - (cell.Rin)^2 ) * cell.width  # π (R_out^2 - R_in^2) * height
cell.rho = (2 * PE.rho * PE.thickness + 2 * NE.rho * NE.thickness + 2 * SP.rho * SP.thickness + PCC.rho * PCC.thickness + NCC.rho * NCC.thickness) / (2 * PE.thickness + 2 * NE.thickness + 2 * SP.thickness + PCC.thickness + NCC.thickness)
cell.mass = cell.rho * cell.volume # 质量 (kg)
cell.layer = 2 * (PE.thickness + NE.thickness + SP.thickness) + PCC.thickness + NCC.thickness
cell.heat_Q = (2 * PE.rho * PE.heat_Q * PE.thickness +2 * NE.rho * NE.heat_Q * NE.thickness +2 * SP.rho * SP.heat_Q * SP.thickness +PCC.rho * PCC.heat_Q * PCC.thickness +NCC.rho * NCC.heat_Q * NCC.thickness) / (2 * PE.rho * PE.thickness + 2 * NE.rho * NE.thickness + 2 * SP.rho * SP.thickness + PCC.rho * PCC.thickness + NCC.rho * NCC.thickness)
cell.alphaT = 0.
cell.h = 10.      # 对流换热系数 (W/m²/K)
cell.T0 = 298.15      # 初始温度 (K)
cell.T_amb = cell.T0

# Binder
binder = Binder()
# binder.rho = 

# === CZM interface parameters on current collectors ===
# --- PCC：PE-PCC 界面参数（Mode I + Mode II）---
# TODO 用户提供实测值（以下为占位，参考旧单一值 σ_max=82e6, G_c=25.3, K=2.4e17）
PCC.σ_max = 82e6       # [Pa] TODO 用户提供实测值
PCC.K_n   = 2.4e17     # [Pa/m] TODO 用户提供实测值
PCC.δ_0   = PCC.σ_max / PCC.K_n
PCC.G_c   = 25.3       # [J/m²] TODO 用户提供实测值
PCC.δ_c   = 2.0 * PCC.G_c / PCC.σ_max
# Mode II（若无独立测量，沿用 Mode I）
PCC.τ_max = PCC.σ_max
PCC.K_t   = PCC.K_n
PCC.δ_0_t = PCC.τ_max / PCC.K_t
PCC.G_c_t = PCC.G_c
PCC.δ_c_t = 2.0 * PCC.G_c_t / PCC.τ_max
PCC.eta = 1.45                 # BK 准则指数 [-]
PCC.h_c0 = 1e7                 # 界面热阻（沿用旧值）
PCC.k_air = 0.026
PCC.lambda_m = 70e-9
PCC.beta = 1.0
PCC.threshold = 70e-9

# --- NCC：NE-NCC 界面参数（Mode I + Mode II）---
NCC.σ_max = 92e6       # [Pa] TODO 用户提供实测值
NCC.K_n   = 1.2e17     # [Pa/m] TODO 用户提供实测值
NCC.δ_0   = NCC.σ_max / NCC.K_n
NCC.G_c   = 6.2        # [J/m²] TODO 用户提供实测值
NCC.δ_c   = 2.0 * NCC.G_c / NCC.σ_max
NCC.τ_max = NCC.σ_max
NCC.K_t   = NCC.K_n
NCC.δ_0_t = NCC.τ_max / NCC.K_t
NCC.G_c_t = NCC.G_c
NCC.δ_c_t = 2.0 * NCC.G_c_t / NCC.τ_max
NCC.eta = 1.45                 # BK 准则指数（与 PCC 同值）
NCC.h_c0 = 1e7                 # 界面热阻（与 PCC 同值）
NCC.k_air = 0.026
NCC.lambda_m = 70e-9
NCC.beta = 1.0
NCC.threshold = 70e-9
cell.lambda_r = (2 * NE.thickness + 2 * SP.thickness + 2 * PE.thickness + PCC.thickness + NCC.thickness) /(2 * NE.thickness/NE.lambda + 2 * SP.thickness/SP.lambda + 2 * PE.thickness/PE.lambda + PCC.thickness/PCC.lambda + NCC.thickness/NCC.lambda)
cell.lambda_t = (2 * NE.thickness*NE.lambda + 2 * SP.thickness*SP.lambda + 2 * PE.thickness*PE.lambda + PCC.thickness*PCC.lambda + NCC.thickness*NCC.lambda) /(2 * NE.thickness + 2 * SP.thickness + 2 * PE.thickness + PCC.thickness + NCC.thickness)

# Scale
scale = Scale()

# assemble to "param_dim"
param_dim = Params(PE, NE, EL, SP, cell, PCC, NCC, tab, binder, scale)
