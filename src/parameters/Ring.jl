"""
    Ring cell parameters (thermal/mechanical only)
"""

PE = Electrode()
NE = Electrode()
EL = Electrolyte()
SP = Separator()
PCC = CurrentCollector()
NCC = CurrentCollector()

tab = Tab()
binder = Binder()
cohesive = Cohesive()

cohesive.h_c0 = 1e7
cohesive.k_air = 0.026
cohesive.lambda_m = 70e-9
cohesive.beta = 1.0
cohesive.threshold = 70e-9
scale = Scale()

cell = Cell()
cell.Rin = 1.92e-3
cell.Rout = 0.021 / 2
cell.width = 1.0
cell.length = 1.0
cell.volume = pi * (cell.Rout^2 - cell.Rin^2) * cell.width
cell.rho = 2500.0
cell.heat_Q = 1000.0
cell.h = 100.0
cell.T0 = 298.0
cell.T_amb = 298.0
cell.lambda_r = 1.71
cell.lambda_t = 34.91
cell.alphaT = 0.0

# Use ring radial conductivity as thermal reference scale
PE.lambda = cell.lambda_r

NE.E = 2.0e10
NE.nu = 0.28
NE.alphaT = 3.0e-6

PE.E = 5.0e10
PE.nu = 0.30
PE.alphaT = 1.0e-5

param_dim = Params(PE, NE, EL, SP, cell, PCC, NCC, tab, binder, scale, cohesive)
