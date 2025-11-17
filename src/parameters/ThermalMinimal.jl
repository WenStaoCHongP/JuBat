"""
    ThermalMinimal: Minimal parameter set for pure 2D thermal FEM validation

Goals:
- Provide only the parameters required by ThermalDistributed2D + ThermalDistributed2D_BC
- Decouple from electrochemistry as much as possible
- Keep everything in dimensional SI units here; dimensionless thermal Scheme-B is derived in ChooseCell

Minimal required fields (SI):
- Cell: T0 [K], T_amb [K], h [W/(m^2 K)], rho [kg/m^3], heat_Q (cp) [J/(kg K)]
- Layers (for k reference and optional anisotropy mixing): PE.lambda, NE.lambda, SP.lambda [W/(m K)]
- Geometric scales for the thermal model to build L_th:
  - Either cell.Rout (>0) for cylindrical/jellyroll-like, or
  - fallback uses max(cell.length, cell.width, L) with L = t_pos + t_sep + t_neg
- Thicknesses (used for L and some ohmic placeholders if needed by other utilities):
  - PE.thickness, SP.thickness, NE.thickness [m]

Notes:
- Electrochemical-specific fields are left at benign defaults.
- k_ref chosen from PE.lambda when available.
"""

PE = Electrode()
PE.thickness = 75e-6
PE.lambda = 2.0
PE.rho = 3200.0
PE.heat_Q = 700.0

NE = Electrode()
NE.thickness = 85e-6
NE.lambda = 1.5
NE.rho = 1700.0
NE.heat_Q = 700.0

SP = Separator()
SP.thickness = 12e-6
SP.lambda = 0.2
SP.rho = 400.0
SP.heat_Q = 700.0

# Electrolyte (only needed for some utilities; keep simple constants)
EL = Electrolyte()
EL.De = (x, y=0)-> 1.0e-10
EL.kappa = (x, y=0)-> 1.0
EL.dlnf_dlnc = (x, y=0)-> 1.0
EL.rho = 1200.0
EL.heat_Q = 2000.0
EL.tplus = 0.3
EL.ce0 = 1000.0

# Current collectors (kept minimal; lambda for completeness in thermal)
PCC = CurrentCollector()
PCC.thickness = 16e-6
PCC.lambda = 200.0
PCC.rho = 2700.0
PCC.heat_Q = 900.0
PCC.sig = 3.0e7

NCC = CurrentCollector()
NCC.thickness = 12e-6
NCC.lambda = 380.0
NCC.rho = 8900.0
NCC.heat_Q = 380.0
NCC.sig = 5.0e7

# Cell thermal envelope
cell = Cell()
cell.length = 0.10          # [m] fallback planar size if Rout not provided
cell.width  = 0.10          # [m]
cell.I1C = 1.0              # dummy electrochem scale (won't be used)
cell.no_layers = 1
cell.capacity = 1.0
cell.cooling_surface = cell.length * 4.0
cell.v_h = 4.2; cell.v_l = 2.5
cell.volume = cell.length * cell.width * (PE.thickness + SP.thickness + NE.thickness)
cell.rho = 2000.0           # effective density [kg/m^3]
cell.heat_Q = 800.0         # effective cp [J/(kg K)]
cell.h = 10.0               # convection coefficient [W/(m^2 K)]
cell.T0 = 298.0
cell.T_amb = 298.0
cell.area = cell.length * cell.width
# Optional cylindrical scale; keep 0 to use planar fallback in ChooseCell
cell.Rin = 0.0
cell.Rout = 0.0
cell.height = 0.0

# Tab and binder: placeholders
tab = Tab()
binder = Binder()

# Scale placeholder; ChooseCell will fill thermal Scheme-B scales later
scale = Scale()

# Minimal additional electrochem placeholders to avoid errors
PE.cs_max = max(1.0, PE.cs_max)
NE.cs_max = max(1.0, NE.cs_max)
EL.ce0 = max(1.0, EL.ce0)

param_dim = Params(PE, NE, EL, SP, cell, NCC, PCC, tab, binder, scale)
