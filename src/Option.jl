using Parameters
"""
    CycleOption - 充放电循环参数
"""
@with_kw mutable struct CycleOption
    n_cycles::Int64 = 50                    # 循环次数
    # 各阶段时长 (s)
    t_charge::Float64 = 3600.0              # 充电时长
    t_rest1::Float64 = 600.0                # 充电后静置
    t_discharge::Float64 = 3600.0           # 放电时长
    t_rest2::Float64 = 600.0                # 放电后静置
    # 电流设置 (A)，均为正值
    I_charge::Float64 = 5.0                 # 充电电流
    I_discharge::Float64 = 5.0              # 放电电流
    # 截止条件
    V_upper::Float64 = 4.2                  # 充电截止电压
    V_lower::Float64 = 2.5                  # 放电截止电压
    # 初始状态
    SOC_init::Float64 = 0.05                # 初始SOC
    # 温度重置
    reset_T_each_cycle::Bool = true         # 每循环重置温度
    reset_T_before_charge::Bool = true      # 充电前重置温度（循环内）
    # 时间步长
    dt_cycle::Vector{Float64} = [1.0, 10.0] # [dt_min, dt_max]
end

"""
    PhaseType - 循环阶段类型枚举
"""
@enum PhaseType begin
    PHASE_CHARGE = 1      # 充电：I < 0
    PHASE_REST = 2        # 静置：I = 0
    PHASE_DISCHARGE = 3   # 放电：I > 0
end
@with_kw mutable struct CzmOptions
    enabled::Bool = false                # 是否启用CZM损伤模型
    model::String = "model1"             # "model1" or "mix"
    update_interval::Int64 = 1           # CZM损伤更新间隔（时间步数），1=每步更新
    soh_threshold::Float64 = 0.8         # SOH终止阈值，SOH≤此值时终止循环
    inner_exit_only::Bool = true         # 断裂时仅内圈单元退出电化学反应
    fix_inner::Bool = true               # CZM边界条件：true=内外圈均固定，false=仅外圈固定
    iter_method::String = "basic"        # "basic" | "load_substep" | "arc_length"
    max_iter::Int64 = 100                # CZM 牛顿迭代最大步数
    tol::Float64 = 1e-4                  # CZM 收敛容差
    load_steps::Int64 = 2                # 载荷子步数（load_substep 模式）
    arc_length_alpha::Float64 = 1.0      # 弧长法系数（arc_length 模式）
    viscous_enabled::Bool = false        # 粘性正则化开关
    viscous_tau::Float64 = 0.0           # 物理松弛时间 [s]，推荐 10~100 s
    area_loss_enabled::Bool = false      # 启用渐进式有效面积损失（D > threshold 时缩减有效面积）
    area_loss_threshold::Float64 = 0.83  # 面积开始缩减的损伤阈值
    geo_nonlinear::Bool = false          # 完全 Green-Lagrange TL 残差 + 标准初应力 K_G（Batch 2）
    winding_prestress::Bool = false      # 卷绕预应力初始应力场 σ₀(r)，缺参即 error（Batch 2'）
    j2_plasticity::Bool = false          # PCC/NCC 平面应力一致 J2 返回映射，缺 sigma_y 即 error（Batch 3）
    continuous_feedback::Bool = false    # 连续损伤–电–热反馈 + 界面热阻（Batch 6）
    friction_mu::Float64 = 0.10          # SP Coulomb 摩擦系数（Batch 8 预留，当前无消费者）
end

@with_kw mutable struct Option
#   option for a lithium-ion battery model
    Np::Int64 = 10
    Ns::Int64 = 10
    Nn::Int64 = 10
    Nrp::Int64 = 10
    Nrn::Int64 = 10
    model::String  = "SPM"
    time::Array{Float64} = [0 3600]
    meshType::String  = "L2"
    gsorder::Int64 = 4
    dimension::Int64 = 1
    #opt.load = {"constant discharge 1C for 1h"}
    Current::Function = x-> 0
    coupleMethod:: String  = "fully coupled"
    coupleOrder::Int64 = 0
    y0::Array{Float64} = []
    dt::Array{Float64} = [1, 100]
    dtType::String  = "constant" # auto or manual
    dtThreshold::Float64 = 0.01
    solveType::String  = "Crank-Nicolson" # forward, backward or Crank-Nicolson
    outputType::String  = "auto" # auto or manual
    jacobi::String = "constant" # constant or update
    thermalmodel::String  = "none" # none, lumped, distributed2D
    mechanicalmodel::String = "none" #none or full
    cite::Vector{String} = String[]
    # Thermal module
    thermal_enabled::Bool = false      # whether thermal module is active
    thermal_dim::String = "1D"        # "1D" or "2D" for distributed models
    thermalmeshType::String = "L2"    # 1D: L2/L3; 2D: Q4 (default)
    cool_method::String = "tab"      # "none", "tab" or "surface"
    collector_seeded::Bool = false     # use collector-seeded band mesh semantics (layer_weights)
    per_element_spme::Bool = false     # allow passing per-element I_app and T to SPMe
    czm::CzmOptions = CzmOptions()     # CZM (Cohesive Zone Model) 模块选项（嵌套收敛，2026-08-30 重构）
    debug_coupling::Bool = false       # print detailed logs for electro-thermal coupling
    debug_log_path::String = "output/debug.log"  # debug log file path
end