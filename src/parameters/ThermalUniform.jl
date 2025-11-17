"""
    ThermalUniform: Simple uniform material parameter set for 2D thermal verification

Properties (SI): k, rho, cp, thicknesses assembled to form cell.volume, convection h and ambient T_amb.
This file provides a compact minimal Params struct compatible with ChooseCell/SetCase.
"""

PE = JuBat.Electrode()
PE.thickness = 75e-6
PE.lambda = 1.0    # W/(m K)
PE.rho = 2500.0
PE.heat_Q = 800.0

NE = JuBat.Electrode()
NE.thickness = 75e-6
NE.lambda = 1.0
NE.rho = 2500.0
NE.heat_Q = 800.0

SP = JuBat.Separator()
SP.thickness = 50e-6
SP.lambda = 1.0
SP.rho = 1200.0
SP.heat_Q = 1000.0

EL = JuBat.Electrolyte()
EL.De = (x, y=0)-> 1.0e-10
EL.kappa = (x, y=0)-> 1.0
EL.dlnf_dlnc = (x, y=0)-> 1.0
EL.rho = 1000.0
EL.heat_Q = 2000.0
EL.tplus = 0.3
EL.ce0 = 1000.0

PCC = JuBat.CurrentCollector()
PCC.thickness = 16e-6
PCC.lambda = 200.0
PCC.rho = 2700.0
PCC.heat_Q = 900.0

NCC = JuBat.CurrentCollector()
NCC.thickness = 12e-6
NCC.lambda = 200.0
NCC.rho = 2700.0
NCC.heat_Q = 900.0

cell = JuBat.Cell()
cell.length = 0.02    # 2 cm square
cell.width  = 0.02
cell.I1C = 1.0
cell.no_layers = 1
cell.capacity = 1.0
cell.cooling_surface = 4*cell.length
cell.v_h = 4.2; cell.v_l = 2.5
cell.volume = cell.length * cell.width * (PE.thickness + SP.thickness + NE.thickness)
cell.rho = 2500.0
cell.heat_Q = 800.0
cell.h = 10.0      # W/(m^2 K) convection
cell.T0 = 298.0
cell.T_amb = 298.0
cell.area = cell.length * cell.width
cell.Rin = 0.0; cell.Rout = 0.0; cell.height = 0.0

tab = JuBat.Tab()
binder = JuBat.Binder()

scale = JuBat.Scale()

PE.cs_max = max(1.0, PE.cs_max)
NE.cs_max = max(1.0, NE.cs_max)
EL.ce0 = max(1.0, EL.ce0)

param_dim = JuBat.Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)
