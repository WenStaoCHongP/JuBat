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
tab.width = 40e-3
tab.length = 0.75 * 99.06e-3
tab.area = tab.width * tab.length * 2
# specify tab angles on circumference (radians)
tab.theta_pos = []
tab.theta_neg = []

# Cell parameters for thermal model
cell = Cell()    
cell.Rout = 0.021/2
cell.Rin = 1.92e-3
cell.length = 2.10   # Using outer diameter proxy (not used by scaling below)
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
cell.h = 150.      # 对流换热系数 (W/m²/K)
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

# Now that PCC/NCC are defined, compute effective thermal conductivities
cell.lambda_r = (NE.thickness + SP.thickness + PE.thickness + PCC.thickness + NCC.thickness) /
                    (NE.thickness/NE.lambda + SP.thickness/SP.lambda + PE.thickness/PE.lambda + PCC.thickness/PCC.lambda + NCC.thickness/NCC.lambda)
cell.lambda_t = (NE.thickness*NE.lambda + SP.thickness*SP.lambda + PE.thickness*PE.lambda + PCC.thickness*PCC.lambda + NCC.thickness*NCC.lambda) /
                    (NE.thickness + SP.thickness + PE.thickness + PCC.thickness + NCC.thickness)

# Scale
scale = Scale()

# assemble to "param_dim"
param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)
