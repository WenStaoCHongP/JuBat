"""
    Ring cell parameters (thermal/mechanical only)
"""

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
cohesive = Cohesive()

cohesive.h_c0 = 1e7
cohesive.k_air = 0.026
cohesive.lambda_m = 70e-9
cohesive.beta = 1.0
cohesive.threshold = 70e-9
scale = Scale()

cell = Cell()
cell.Rin = 1.92e-3
cell.Rout = 0.0203 / 2
cell.width = 1.0
cell.length = 1.0
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
cell.area = 1.0
cell.no_layers = 1

param_dim = Params(PE, NE, EL, SP, cell, PCC, NCC, tab, binder, scale, cohesive)
