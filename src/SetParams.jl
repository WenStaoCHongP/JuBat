"""
        Set up a structure for parameters of a lithium-ion cell
        The names of variables use the beginning uppercase words to specify the domain:    
            NE -- negative electrode
            PE -- positive electrode
            SP -- separator
            EL -- electrolyte
            NCC -- negative current collector
            PCC -- positive current collector
            Cell -- cell 
            Tab -- tab
            BD -- binder

        The variables are for properties of :
            theta -- stoichiometry
            thickness -- thickness [m]
            lambda -- thermal conductivities [W/(m K)]
            Ds -- solid phase diffusion coefficients [m^2/s]
            rho -- density [kg/m^3]
            heat_Q -- specific heat capacities [J/(kg K)]
            sig -- conductivity [S/m]
            eps -- porosity
            brugg -- Bruggeman coefficients
            k -- reaction rate constant [ A m^2.5/mol^1.5]
            tplus -- transference number 
            h -- heat exchange coefficient [W/(m^2 K)]
            Rs -- electrode particle radius [m]
            cs_max -- maximum concentration in solid phase [mol/m^3]
            ce0 -- electrolyte lithium-ions initial concentration [mol/m^3]
            U -- open-circuit potential (OCP) [V]
            dUdT -- entropic change of the OCP [V/K]
            theta_0 --  Theta @ 0% Lithium Concentration 
            theta_100 --  Theta @ 100% Lithium Concentration 
            alpha -- Alpha Factor
            dlnf_dlnc -- activity of electrolyte, i.e. 1 + dln(f)/dln(ce)
"""
@with_kw mutable struct Electrode
    theta_100::Float64 = 0
    theta_0::Float64 = 0
    thickness::Float64 = 0
    lambda::Float64 = 0
    Ds::Float64 = 0
    rho::Float64 = 0
    heat_Q::Float64 = 0
    eps::Float64 = 0
    eps_fi::Float64 = 0
    eps_s::Float64 = 0
    brugg::Float64 = 0
    k::Float64 = 0
    cs_max::Float64 = 0
    cs0::Float64 = 0
    rs::Float64 = 0
    as::Float64 = 0
    sig::Float64 = 0
    # Mechanical/thermal properties for stress
    E::Float64 = 0            # Young's modulus [Pa]
    nu::Float64 = 0           # Poisson's ratio [-]
    alphaT::Float64 = 0       # Thermal expansion coefficient [1/K]
    Omega::Float64 = 0        # Partial molar volume [m^3/mol]
    E_coat::Float64 = 0       # 极片（涂层）宏观弹性模量 [Pa]
    nu_coat::Float64 = 0      # 极片（涂层）宏观泊松比 [-]
    Eac_D::Float64 = 0
    Eac_k::Float64 = 0
    alpha::Float64 = 0
    U::Function = x-> 0.0
    dUdT::Function = x-> 0.0
    M_d::SparseArrays.SparseMatrixCSC{Float64, Int64} = spzeros(0,0)
    K_d::SparseArrays.SparseMatrixCSC{Float64, Int64} = spzeros(0,0)
    M_p::SparseArrays.SparseMatrixCSC{Float64, Int64} = spzeros(0,0)
    K_p::SparseArrays.SparseMatrixCSC{Float64, Int64} = spzeros(0,0)
end

@with_kw mutable struct Separator
    thickness::Float64 = 0
    lambda::Float64 = 0
    rho::Float64 = 0
    heat_Q::Float64 = 0
    eps::Float64 = 0
    eps_fi::Float64 = 0
    brugg::Float64 = 0
    E::Float64 = 0            # 弹性模量 [Pa]
    nu::Float64 = 0           # 泊松比 [-]
    alphaT::Float64 = 0       # 热膨胀系数 [1/K]
end

@with_kw mutable struct CurrentCollector
    thickness::Float64 = 0
    lambda::Float64 = 0
    rho::Float64 = 0 
    heat_Q::Float64 = 0
    sig::Float64 = 0
    E::Float64 = 0            # 弹性模量 [Pa]
    nu::Float64 = 0           # 泊松比 [-]
    alphaT::Float64 = 0       # 热膨胀系数 [1/K]
end

@with_kw mutable struct Electrolyte
    De::Function = (x,y=0)-> 0
    kappa::Function = (x,y=0)-> 0
    dlnf_dlnc::Function = (x,y=0)-> 0
    rho::Float64 = 0 
    heat_Q::Float64 = 0
    tplus::Float64 = 0
    ce0::Float64 = 0
    Eac_D::Float64 = 0
    Eac_k::Float64 = 0
end

@with_kw mutable struct Cell
    length::Float64 = 0
    width::Float64 = 0
    wrapper::Float64 = 0
    I1C::Float64 = 0
    no_layers::Int32 = 0 
    capacity::Float64 = 0 
    cooling_surface::Float64 = 0 
    area :: Float64 = 0
    v_h::Float64 = 0 
    v_l::Float64 = 0
    volume::Float64 = 0 
    rho::Float64 = 0
    mass::Float64 = 0
    alphaT::Float64 = 0
    heat_Q::Float64 = 0 
    h::Float64 = 0
    T0::Float64 = 298.
    T_amb::Float64 = 0
    # Jellyroll geometry (optional)
    Rin::Float64 = 0.0
    Rout::Float64 = 0.0
    # Effective thermal conductivities for jellyroll (optional)
    layer::Float64 = 0.0
    lambda_r::Float64 = 0.0
    lambda_t::Float64 = 0.0
    # Thermal mesh divisions for jellyroll top-view (optional)
    Nr_th::Int = 0
    Nθ_th::Int = 0
    n_windings::Int = 0
end

@with_kw mutable struct Tab
    length::Float64 = 0
    width::Float64 = 0
    area::Float64 = 0
    h::Float64 = 0
    theta_pos::Vector{Float64} = Float64[]
    theta_neg::Vector{Float64} = Float64[]
end
# param_dim.Tab.width = 0.75 * param_dim.Tab.length  
# param_dim.Tab.area = param_dim.Tab.length * param_dim.Tab.width

@with_kw mutable struct Binder
    rho::Float64 = 0
end

@with_kw mutable struct Cohesive
    # 法向 (Mode I)
    σ_max_n::Float64 = 0.0    # 最大法向牵引力 [Pa]
    δ_0_n::Float64 = 0.0      # 损伤起始分离位移 [m]
    δ_c_n::Float64 = 0.0      # 临界（完全断裂）分离位移 [m]
    G_c_n::Float64 = 0.0      # 法向断裂能 [J/m²]
    K_n::Float64 = 0.0        # 法向初始刚度（惩罚刚度）[Pa/m]
    
    # 切向 (Mode II)
    τ_max_t::Float64 = 0.0    # 最大切向牵引力 [Pa]
    δ_0_t::Float64 = 0.0      # 损伤起始切向位移 [m]
    δ_c_t::Float64 = 0.0      # 临界切向位移 [m]
    G_c_t::Float64 = 0.0      # 切向断裂能 [J/m²]
    K_t::Float64 = 0.0        # 切向初始刚度 [Pa/m]
    
    # 混合模式参数
    eta::Float64 = 1.0        # BK准则指数（Benzeggagh-Kenane）[-]
    czm_model::String = "model1"    # 模型选择:"model1"仅法向（Mode I）,"mix"混合模式（Mode I + II）

    # Interface thermal resistance parameters (physical units)
    h_c0::Float64 = 1e7        # 完好接触换热系数 [W/(m^2·K)]
    k_air::Float64 = 0.026     # 空气热导率 [W/(m·K)]
    lambda_m::Float64 = 70e-9  # 平均自由程 [m]
    beta::Float64 = 1.0        # 参数 beta [-]
    threshold::Float64 = 70e-9 # 阈值厚度 [m]

    # Viscous regularization (normalized)
    tau_visc::Float64 = 0.0     # 归一化粘性松弛时间 τ_v* = τ_v / t0 [-]
end
@with_kw mutable struct Scale
    L::Float64 = 1e-6
    r0::Float64 = 1e-6
    a0::Float64 = 1/r0
    t0::Float64 = 3600
    T_ref::Float64 = 298.
    F::Float64 = 96485.33289
    R::Float64 = 8.314
    j::Float64 = 0
    Ds_p::Float64 = 0
    Ds_n::Float64 = 0
    ts_p::Float64 = 0
    ts_n::Float64 = 0
    te::Float64 = 0
    De::Float64 = 0
    phi::Float64 = 0
    sig::Float64 = 0
    kappa::Float64 = 0
    cp_max::Float64 = 0
    cn_max::Float64 = 0
    ce::Float64 = 0
    k_p::Float64 = 0
    k_n::Float64 = 0
    I_typ::Float64 = 0
    R_cell::Float64 = 0
    E_n::Float64 = 0
    E_p::Float64 = 0
    E_coat::Float64 = 0       # 极片模量参考尺度 [Pa]
    # --- Thermal scaling (统一能量尺度) ---
    rho::Float64 = 0          # 密度尺度 = 电池平均密度 [kg/m³]
    P_ref::Float64 = 0
    lambda::Float64 = 0          # 导热率尺度参数 P_ref/(L*T_ref) [W/(m·K)]
    q::Float64 = 0               # 热源尺度参数 P_ref/L^3 [W/m³]
    h::Float64 = 0               # Biot 数 = h_cell*L/lambda_r (边界条件)
    # --- Cohesive zone model scaling ---
    σ_czm::Float64 = 0         # reference cohesive traction [Pa] (typically σ_max_n)
    δ_czm::Float64 = 0         # reference separation displacement [m] (typically δ_c_n)
    G_czm::Float64 = 0         # reference fracture energy [J/m²] (σ_czm * δ_czm)
    K_czm::Float64 = 0         # reference cohesive stiffness [Pa/m] (σ_czm / δ_czm)
end

@with_kw mutable struct Params
    PE::Electrode
    NE::Electrode
    EL::Electrolyte
    SP::Separator
    cell::Cell
    PCC::CurrentCollector
    NCC::CurrentCollector
    tab::Tab
    binder::Binder
    scale::Scale
    cohesive::Cohesive = Cohesive()  # 内聚力模型参数（可选）
end

function ChooseCell(CellType::String="LG M50")
"""
    This is a function choose a cell
    Input - CellType::String, including options of 
        1. "LG M50" for the LG M50 cells
        2. "Northrop"  for the Northrop cells from LIONSIMBA
    Output - param_dim::Params for all param_dim parameters
"""
    params_dir = joinpath(@__DIR__, "parameters")
    if CellType == "LG M50"
        include(joinpath(params_dir, "LGM50.jl"))
    elseif CellType == "Northrop"
        include(joinpath(params_dir, "Northrop.jl"))
    elseif CellType == "Enertech"
        include(joinpath(params_dir, "Enertech.jl"))
    elseif CellType == "Jellyroll"
        include(joinpath(params_dir, "Jellyroll.jl"))
    elseif CellType == "Ring"
        include(joinpath(params_dir, "Ring.jl"))
    end
    param_dim.PE.eps_s = 1 - param_dim.PE.eps - param_dim.PE.eps_fi
    param_dim.NE.eps_s = 1 - param_dim.NE.eps - param_dim.NE.eps_fi
    param_dim.PE.as = 3 * param_dim.PE.eps_s / param_dim.PE.rs
    param_dim.NE.as = 3 * param_dim.NE.eps_s / param_dim.NE.rs
    if abs(param_dim.cell.rho) < 1e-8
        param_dim.cell.rho = (
            param_dim.PE.rho * param_dim.PE.thickness + param_dim.NE.rho * param_dim.NE.thickness +
            param_dim.SP.rho * param_dim.SP.thickness + param_dim.NCC.rho * param_dim.NCC.thickness + param_dim.PCC.rho * param_dim.PCC.thickness
        ) / (param_dim.PE.thickness + param_dim.NE.thickness + param_dim.SP.thickness + param_dim.NCC.thickness + param_dim.PCC.thickness)
    else
        param_dim.cell.rho = param_dim.cell.rho
    end
    if abs(param_dim.cell.heat_Q) < 1e-8
        param_dim.cell.heat_Q = (
            param_dim.PE.heat_Q * param_dim.PE.rho * param_dim.PE.thickness + param_dim.NE.heat_Q * param_dim.NE.rho * param_dim.NE.thickness +
            param_dim.SP.heat_Q * param_dim.SP.rho * param_dim.SP.thickness + param_dim.NCC.heat_Q * param_dim.NCC.rho * param_dim.NCC.thickness +
            param_dim.PCC.heat_Q * param_dim.PCC.rho * param_dim.PCC.thickness
        ) / (param_dim.PE.rho * param_dim.PE.thickness + param_dim.NE.rho * param_dim.NE.thickness + param_dim.SP.rho * param_dim.SP.thickness + 
        param_dim.NCC.rho * param_dim.NCC.thickness  + param_dim.PCC.rho * param_dim.PCC.thickness
        ) 
    else
        param_dim.cell.heat_Q = param_dim.cell.heat_Q
    end  
    param_dim.cell.mass = param_dim.cell.rho * param_dim.cell.volume
    param_dim.scale.I_typ = param_dim.cell.I1C
    param_dim.scale.L = param_dim.PE.thickness + param_dim.NE.thickness + param_dim.SP.thickness
    param_dim.scale.j = param_dim.scale.I_typ / param_dim.scale.a0 / param_dim.scale.L / param_dim.cell.area
    param_dim.scale.ts_p = param_dim.scale.F * param_dim.PE.cs_max * param_dim.cell.area * param_dim.scale.L / param_dim.scale.I_typ
    param_dim.scale.ts_n = param_dim.scale.F * param_dim.NE.cs_max * param_dim.cell.area * param_dim.scale.L / param_dim.scale.I_typ 
    param_dim.scale.te = param_dim.scale.F * param_dim.EL.ce0 * param_dim.cell.area * param_dim.scale.L / param_dim.scale.I_typ
    param_dim.scale.Ds_p = param_dim.scale.r0^2 / param_dim.scale.ts_p
    param_dim.scale.Ds_n = param_dim.scale.r0^2 / param_dim.scale.ts_n

    param_dim.scale.De = param_dim.scale.L^2 / param_dim.scale.te
    param_dim.scale.phi = param_dim.scale.T_ref * param_dim.scale.R / param_dim.scale.F
    param_dim.scale.sig = param_dim.scale.L * param_dim.scale.I_typ / param_dim.scale.phi / param_dim.cell.area 
    param_dim.scale.kappa = param_dim.scale.L * param_dim.scale.I_typ / param_dim.scale.phi / param_dim.cell.area 
    param_dim.scale.cp_max = param_dim.PE.cs_max
    param_dim.scale.cn_max = param_dim.NE.cs_max
    param_dim.scale.ce = param_dim.EL.ce0
    param_dim.scale.E_n = param_dim.NE.cs_max * param_dim.scale.R * param_dim.scale.T_ref
    param_dim.scale.E_p = param_dim.PE.cs_max * param_dim.scale.R * param_dim.scale.T_ref
    param_dim.scale.k_p = param_dim.scale.j / param_dim.PE.cs_max / sqrt(param_dim.EL.ce0)
    param_dim.scale.k_n = param_dim.scale.j / param_dim.NE.cs_max / sqrt(param_dim.EL.ce0)
    param_dim.scale.R_cell = param_dim.scale.phi / param_dim.scale.I_typ
    param_dim.scale.rho = param_dim.cell.rho
    param_dim.scale.P_ref = param_dim.scale.phi * param_dim.scale.I_typ
    param_dim.scale.lambda = param_dim.scale.P_ref / (param_dim.scale.L * param_dim.scale.T_ref)
    param_dim.scale.h = param_dim.cell.h * param_dim.scale.L / param_dim.cell.lambda_r  # Biot 数
    param_dim.scale.q = param_dim.scale.P_ref / param_dim.scale.L^3
    param_dim.scale.σ_czm = param_dim.cohesive.σ_max_n
    param_dim.scale.δ_czm = param_dim.scale.L
    param_dim.scale.G_czm = param_dim.scale.σ_czm * param_dim.scale.δ_czm
    param_dim.scale.K_czm = param_dim.scale.σ_czm / param_dim.scale.δ_czm
    return param_dim
end

function NormaliseParam(param_dim::Params)
    """
        This is a function to normalise the parameters
        Input - param_dim::Params (with units) for a cell
        Output - param::Params (normalised)
    """
    # normalise the parameters, while thermal parameters are not covered, need revisit
    param = deepcopy(param_dim)

    # posotove electrode
    param.PE.theta_100 = param_dim.PE.theta_100
    param.PE.theta_0 = param_dim.PE.theta_0
    param.PE.cs0 = param_dim.PE.cs0 / param_dim.PE.cs_max
    param.PE.thickness = param_dim.PE.thickness / param.scale.L 
    param.PE.Ds = param_dim.PE.Ds / param.scale.Ds_p
    param.PE.eps = param_dim.PE.eps
    param.PE.eps_fi = param_dim.PE.eps_fi
    param.PE.brugg = param_dim.PE.brugg
    param.PE.k = param_dim.PE.k / param.scale.k_p
    param.PE.rs = param_dim.PE.rs / param.scale.r0
    param.PE.sig = param_dim.PE.sig / param.scale.sig
    param.PE.E = param_dim.PE.E / param_dim.scale.E_p
    param.PE.nu = param_dim.PE.nu
    param.PE.alphaT = param_dim.PE.alphaT* param.scale.T_ref
    param.PE.Omega = param_dim.PE.Omega * param_dim.PE.cs_max
    param.PE.U = x-> Base.invokelatest(param_dim.PE.U, x) / param.scale.phi
    param.PE.dUdT = x-> Base.invokelatest(param_dim.PE.dUdT, x) / param.scale.phi * param.scale.T_ref
    param.PE.as = param_dim.PE.as / param.scale.a0
    param.PE.Eac_D = param_dim.PE.Eac_D / param.scale.R / param.scale.T_ref
    param.PE.Eac_k = param_dim.PE.Eac_k / param.scale.R / param.scale.T_ref
    param.PE.lambda = param_dim.PE.lambda / param.scale.lambda
    param.PE.rho = param_dim.PE.rho / param.scale.rho
    param.PE.heat_Q = param_dim.PE.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref / (param.scale.t0 * param.scale.phi * param.scale.I_typ)

    # negative electrode
    param.NE.theta_100 = param_dim.NE.theta_100
    param.NE.theta_0 = param_dim.NE.theta_0
    param.NE.cs0 = param_dim.NE.cs0 / param_dim.NE.cs_max
    param.NE.thickness = param_dim.NE.thickness / param.scale.L
    param.NE.Ds = param_dim.NE.Ds / param.scale.Ds_n
    param.NE.eps = param_dim.NE.eps
    param.NE.eps_fi = param_dim.NE.eps_fi
    param.NE.brugg = param_dim.NE.brugg
    param.NE.k = param_dim.NE.k / param.scale.k_n
    param.NE.rs = param_dim.NE.rs / param.scale.r0
    param.NE.sig = param_dim.NE.sig / param.scale.sig
    param.NE.E = param_dim.NE.E / param_dim.scale.E_n
    param.NE.nu = param_dim.NE.nu
    param.NE.alphaT = param_dim.NE.alphaT * param.scale.T_ref
    param.NE.Omega = param_dim.NE.Omega * param_dim.NE.cs_max
    param.NE.U = x-> Base.invokelatest(param_dim.NE.U, x) / param.scale.phi
    param.NE.dUdT = x-> Base.invokelatest(param_dim.NE.dUdT, x) / param.scale.phi * param.scale.T_ref
    param.NE.as = param_dim.NE.as / param.scale.a0
    param.NE.Eac_D = param_dim.NE.Eac_D / param.scale.R / param.scale.T_ref
    param.NE.Eac_k = param_dim.NE.Eac_k / param.scale.R / param.scale.T_ref
    param.NE.lambda = param_dim.NE.lambda / param.scale.lambda
    param.NE.rho = param_dim.NE.rho / param.scale.rho
    param.NE.heat_Q = param_dim.NE.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref / (param.scale.t0 * param.scale.phi * param.scale.I_typ)

    # separator
    param.SP.thickness = param_dim.SP.thickness / param.scale.L
    param.SP.eps = param_dim.SP.eps
    param.SP.eps_fi = param_dim.SP.eps_fi
    param.SP.brugg = param_dim.SP.brugg
    param.SP.lambda = param_dim.SP.lambda / param.scale.lambda
    param.SP.rho = param_dim.SP.rho / param.scale.rho
    param.SP.heat_Q = param_dim.SP.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref / (param.scale.t0 * param.scale.phi * param.scale.I_typ)

    # positive current colloctor
    param.PCC.thickness = param_dim.PCC.thickness / param.scale.L
    param.PCC.sig =  param_dim.PCC.sig / param.scale.sig
    param.PCC.lambda = param_dim.PCC.lambda / param.scale.lambda
    param.PCC.rho = param_dim.PCC.rho / param.scale.rho
    param.PCC.heat_Q = param_dim.PCC.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref / (param.scale.t0 * param.scale.phi * param.scale.I_typ)
    # negative current colloctor
    param.NCC.thickness = param_dim.NCC.thickness / param.scale.L
    param.NCC.sig = param_dim.NCC.sig / param.scale.sig
    param.NCC.lambda = param_dim.NCC.lambda / param.scale.lambda
    param.NCC.rho = param_dim.NCC.rho / param.scale.rho
    param.NCC.heat_Q = param_dim.NCC.heat_Q * param.scale.rho * param.scale.L^3 * param.scale.T_ref / (param.scale.t0 * param.scale.phi * param.scale.I_typ)

    # electrolyte (wrap with invokelatest to avoid world-age issues for closures from parameters)
    param.EL.De = (x, y=1)-> Base.invokelatest(param_dim.EL.De, x * param.scale.ce, y * param.scale.T_ref) / param.scale.De
    param.EL.kappa = (x, y=1)-> Base.invokelatest(param_dim.EL.kappa, x * param.scale.ce, y * param.scale.T_ref) / param.scale.kappa
    param.EL.tplus = param_dim.EL.tplus
    param.EL.ce0 = param_dim.EL.ce0 / param.scale.ce

    # cell
    param.cell.cooling_surface = param_dim.cell.cooling_surface / param_dim.cell.area
    param.cell.h = param_dim.cell.h * param_dim.cell.area * param.scale.T_ref / param.scale.phi / param.scale.I_typ
    param.cell.mass = param_dim.cell.mass / param_dim.cell.mass
    param.cell.rho = param_dim.cell.rho / param.scale.rho
    param.cell.heat_Q = param_dim.cell.heat_Q * param_dim.cell.mass * param.scale.T_ref / param.scale.t0 / param.scale.phi / param.scale.I_typ
    param.cell.T_amb = param_dim.cell.T_amb / param.scale.T_ref 
    param.cell.T0 = param_dim.cell.T0 / param.scale.T_ref 
    param.cell.area = param_dim.cell.area / param_dim.cell.area
    param.cell.volume = param_dim.cell.volume / param_dim.scale.L^3
    param.cell.layer = param_dim.cell.layer / param.scale.L
    param.cell.lambda_r = param_dim.cell.lambda_r / param.scale.lambda
    param.cell.lambda_t = param_dim.cell.lambda_t / param.scale.lambda
    param.cell.Rin = param_dim.cell.Rin / param.scale.L
    param.cell.Rout = param_dim.cell.Rout / param.scale.L
    param.cell.width = param_dim.cell.width / param.scale.L

    #tab
    param.tab.length = param_dim.tab.length / param.scale.L
    param.tab.width = param_dim.tab.width / param.scale.L
    param.tab.area = param_dim.tab.area / param.scale.L^2
    param.tab.h = param_dim.tab.h * param_dim.scale.L / param_dim.cell.lambda_r
    
    # cohesive zone model 
    # 法向参数归一化
    param.cohesive.σ_max_n = param_dim.cohesive.σ_max_n / param_dim.scale.σ_czm
    param.cohesive.δ_0_n = param_dim.cohesive.δ_0_n / param_dim.scale.δ_czm
    param.cohesive.δ_c_n = param_dim.cohesive.δ_c_n / param_dim.scale.δ_czm
    param.cohesive.G_c_n = param_dim.cohesive.G_c_n / param_dim.scale.G_czm
    param.cohesive.K_n = param_dim.cohesive.K_n / param_dim.scale.K_czm
    # 切向参数归一化
    param.cohesive.τ_max_t = param_dim.cohesive.τ_max_t / param_dim.scale.σ_czm
    param.cohesive.δ_0_t = param_dim.cohesive.δ_0_t / param_dim.scale.δ_czm
    param.cohesive.δ_c_t = param_dim.cohesive.δ_c_t / param_dim.scale.δ_czm
    param.cohesive.G_c_t = param_dim.cohesive.G_c_t / param_dim.scale.G_czm
    param.cohesive.K_t = param_dim.cohesive.K_t / param_dim.scale.K_czm
    # BK指数不需要归一化（无量纲）
    param.cohesive.eta = param_dim.cohesive.eta

    # interface thermal resistance parameters (dimensionless)
    param.cohesive.h_c0 = param_dim.cohesive.h_c0 * param_dim.scale.L / param_dim.scale.lambda
    param.cohesive.k_air = param_dim.cohesive.k_air / param_dim.scale.lambda
    param.cohesive.lambda_m = param_dim.cohesive.lambda_m / param_dim.scale.L
    param.cohesive.beta = param_dim.cohesive.beta
    param.cohesive.threshold = param_dim.cohesive.threshold / param_dim.scale.L
    
    return param
end

