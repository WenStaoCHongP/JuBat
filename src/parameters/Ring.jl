"""
    Ring cell parameters (thermal/mechanical only).
    与 Jellyroll 归一化一致：相同 PE/NE/SP thickness → scale.L 一致；
    相同 width/length/volume 约定 → 热源量纲一致，无需桥接。
"""
# 与 Jellyroll 相同的电极/隔膜厚度，保证 scale.L = PE+NE+SP 一致
PE = Electrode()
PE.thickness = 75.6e-6
PE.cs_max = 63104
PE.theta_100 = 0.263849
PE.theta_0 = 0.853974
PE.E = 5.0e10
PE.nu = 0.30
PE.alphaT = 1.0e-5

NE = Electrode()
NE.thickness = 85.2e-6
NE.cs_max = 33133.0
NE.theta_100 = 0.910612
NE.theta_0 = 0.0263472
NE.E = 2.0e10
NE.nu = 0.28
NE.alphaT = 3.0e-6

EL = Electrolyte()
EL.ce0 = 1000

SP = Separator()
SP.thickness = 1.2e-5

PCC = CurrentCollector()
NCC = CurrentCollector()

tab = Tab()
binder = Binder()
# 界面热阻默认值（h_c0=1e7 等）已内置于 CurrentCollector 字段默认值（2026-08-30 重构）
scale = Scale()

cell = Cell()
cell.Rin = 1.92e-3
cell.Rout = 0.0203 / 2
# 与 Jellyroll 一致：width 为轴向高度，volume = π(Rout²−Rin²)*width，保证热源量纲一致
cell.width = 6.5e-2
cell.length = 1.58
cell.volume = pi * (cell.Rout^2 - cell.Rin^2) * cell.width
cell.rho = 2813
cell.heat_Q = 860
cell.h = 10
cell.T0 = 298.0
cell.T_amb = 298.0
cell.lambda_r = 1.318
cell.lambda_t = 25.26
cell.alphaT = 0.0
cell.I1C = 5.0
cell.capacity = 5.0
cell.area = cell.width * cell.length * 1
cell.no_layers = 1

param_dim = Params(PE, NE, EL, SP, cell, PCC, NCC, tab, binder, scale)
