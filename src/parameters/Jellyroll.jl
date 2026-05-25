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
PE.E = 3.75e10         # Pa
PE.nu = 0.20          # -
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
NE.E = 1.5e10         # Pa
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

# Positive Current Collector 
PCC = CurrentCollector()
PCC.thickness = 16e-6
PCC.lambda = 237.   # 热导率 (W/m/K)
PCC.rho = 2702.    # 密度 (kg/m³)
PCC.heat_Q = 8.76e2      # 比热容 (J/kg/K)
PCC.sig =3.55e7

# Negative Current Collector
NCC = CurrentCollector()
NCC.thickness = 12.00e-6
NCC.lambda = 401.   # 热导率 (W/m/K)
NCC.rho = 8933.    # 密度 (kg/m³)
NCC.heat_Q = 3.83e2      # 比热容 (J/kg/K)
NCC.sig = 5.96e7

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

# Cohesive zone model parameters for interlayer interface
# 用于描述相邻卷绕圈之间的界面脱粘行为
cohesive = Cohesive()

# 法向参数 (Mode I - 张开模式)
cohesive.σ_max_n = 92e6 # 最大法向牵引力 [Pa]
cohesive.K_n = 1.2e17 # 继续降低刚度 [Pa/m]，用于检查是否由过硬本构导致不收敛
cohesive.δ_0_n = cohesive.σ_max_n / cohesive.K_n # 损伤起始分离位移 [m]
cohesive.G_c_n = 6.2  # 断裂能 [J/m²]
cohesive.δ_c_n = 2.0 * cohesive.G_c_n / cohesive.σ_max_n  # 临界分离位移 [m]


# 切向参数 (Mode II - 剪切模式)
cohesive.τ_max_t = 92e6      # 最大切向牵引力 [Pa] (0.15 MPa)
cohesive.K_t = 1.2e17  # 继续降低刚度 [Pa/m]，用于检查是否由过硬本构导致不收敛
cohesive.δ_0_t = cohesive.τ_max_t / cohesive.K_t         # 损伤起始切向位移 [m] (35 nm)
cohesive.G_c_t = 6.2  # [J/m²] ≈ 0.0135 J/m²
cohesive.δ_c_t = 2.0 * cohesive.G_c_t / cohesive.τ_max_t        # 临界切向位移 [m] (180 nm)

# 混合模式参数
cohesive.eta = 1.45           # BK准则指数（Benzeggagh-Kenane）[-]

# Interface thermal resistance parameters (from 界面热阻曲线.py)
cohesive.h_c0 = 1e7
cohesive.k_air = 0.026
cohesive.lambda_m = 70e-9
cohesive.beta = 1.0
cohesive.threshold = 70e-9
cell.lambda_r = (2 * NE.thickness + 2 * SP.thickness + 2 * PE.thickness + PCC.thickness + NCC.thickness) /(2 * NE.thickness/NE.lambda + 2 * SP.thickness/SP.lambda + 2 * PE.thickness/PE.lambda + PCC.thickness/PCC.lambda + NCC.thickness/NCC.lambda)
cell.lambda_t = (2 * NE.thickness*NE.lambda + 2 * SP.thickness*SP.lambda + 2 * PE.thickness*PE.lambda + PCC.thickness*PCC.lambda + NCC.thickness*NCC.lambda) /(2 * NE.thickness + 2 * SP.thickness + 2 * PE.thickness + PCC.thickness + NCC.thickness)

# Scale
scale = Scale()

# assemble to "param_dim"
param_dim = Params(PE, NE, EL, SP, cell, PCC, NCC, tab, binder, scale, cohesive)
