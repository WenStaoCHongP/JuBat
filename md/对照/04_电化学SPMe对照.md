# 卷 04 — 电化学 SPMe 对照

> 理论来源：`md/04_电化学模型_SPMe.md`（377 行）
> 代码来源（7 个文件）：
> - `src/SPMe.jl`（287 行）— 主要：SPMe / SPMe_element / SPMe_variables! / SPMe_variables / SPMe_BC
> - `src/ElectrodeDiffusion.jl`（19 行）— 颗粒扩散矩阵装配
> - `src/ElectrodePotential.jl`（20 行）— 固相电位矩阵装配（P2D 用）
> - `src/ElectrolyteDiffusion.jl`（31 行）— 电解液扩散矩阵装配
> - `src/ElectrolytePotential.jl`（24 行）— 液相电位矩阵装配（P2D 用）
> - `src/SPM.jl`（82 行）— 简化版 SPM
> - `src/P2D.jl`（283 行）— 伪二维 P2D 模型（md 04 标注 "P2D 在 SPMe 框架内可选用"，作扩展对照）
> 生成日期：2026-08-02

---

## 编写说明

本卷"实现摘要"列常含简短代码片段（如 `j0=k*Arrhenius.*abs.(cs*(1-cs)*ce).^0.5`）以保证对照的可验证性，因此字数可能超过 spec §3 建议的 25 字。此为有意识取舍——优先**可追溯**而非**摘要性**。

---

## md 04 §1 模型概述

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| — | md/04 §1 三条主要假设（电解液 1D 平均 / 代表性球形颗粒 / 包含电解液动力学） | `src/SPMe.jl:209-287 (SPMe_variables)` | 单点颗粒 + 1D ce_n/ce_sp/ce_p | ✅ | md §1 列出概念假设，非公式；代码通过变量结构体现（cn/cp 单点表面浓度 + 电解液分段） |

---

## md 04 §2 控制方程

### §2.1 颗粒内扩散

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.1) | md/04 §2.1 Fick 扩散 `∂c_s/∂t = (1/r²)·∂(D_s·r²·∂c_s/∂r)/∂r` | `src/ElectrodeDiffusion.jl:1-19 (ElectrodeDiffusion)` | 弱形式装配：M = ∫Ni·Nj·r²·w·J；K = -∫dNi·dNj·Ds·r²·w·J | ✅ | md 给强形式；代码直接装配 FE 弱形式，等价 |
| (4.2) | md/04 §2.1 球心 BC `∂c_s/∂r|_{r=0} = 0` | — | 自然 BC，无需额外装配 | ✅ | md §2.1 明示；代码颗粒网格起点 r=0 不施加 Neumann 项 |
| (4.3) | md/04 §2.1 表面 BC `-D_s·∂c_s/∂r|_{r=Rs} = j/F` | `src/SPMe.jl:184-207 (SPMe_BC)` L189-191 | `flux_np[end] = -j_n*rs^2`（含 r² 因子来自 weak form） | ✅ | md 写通量 -j/F；代码已无量纲化（j 已含 1/F），并乘 r² 因子进入载荷 |
| (4.4) | md/04 §2.1 表面浓度 `c_s^surf(t) = c_s(r=Rs, t)` | `src/SPMe.jl:232-233 (SPMe_variables)` L131-132 (`SPMe_variables!`) | `cn_surf = ws["negative particle surface lithium concentration"]` | ✅ | 直接从状态向量提取末节点 |

### §2.2 电解液物质守恒

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.5) | md/04 §2.2 电解液守恒 `∂(ε_e·c_e)/∂t = ∂(D_e^eff·∂c_e/∂x)/∂x + (1-t_+^0)·a_s·j/F` | `src/ElectrolyteDiffusion.jl:1-31 (ElectrolyteDiffusion)` | M = ∫ε_e·Ni·Nj·w·J；K = -∫dNi·dNj·D_e^eff·w·J | ✅ | md 强形式；代码 FE 装配，对应 §4.2 弱形式 |
| (4.6) | md/04 §2.2 有效扩散系数 `D_e^eff = D_e·ε_e^brugg` | `src/ElectrolyteDiffusion.jl:24-26` | `De_ne = param.EL.De(ce_n_gs,T)*param.NE.eps^param.NE.brugg` | ✅ | 三区分别装配：ne/sp/pe |
| — | md/04 §2.2 区域划分（Ω_n / Ω_s / Ω_p）隔膜 a_s=0, j=0 | `src/SPMe.jl:197-201 (SPMe_BC)` | `v_sp = ...; coeff[v_sp] .*= 0` | ✅ | md §2.2 表行；代码显式将隔膜区源项置零 |
| — | md/04 §2.2 源项 `(1-t_+^0)·a_s·j/F` | `src/SPMe.jl:199-201 (SPMe_BC)` | `coeff[v_ne] .*= (1-tplus)*NE.as*j_n`；`coeff[v_pe] .*= (1-tplus)*PE.as*j_p` | ✅ | md 公式无量纲化后 j 已含 1/F；as 即 a_s |

### §2.3 电荷守恒

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.7) | md/04 §2.3 固相电荷守恒 `-∂(σ^eff·∂φ_s/∂x)/∂x = a_s·j` | `src/ElectrodePotential.jl:1-20 (ElectrodePotential)` | K = ∫dNi·dNj·σ^eff·w·J；M=0（稳态） | 🟡 | md §2.3 给普适方程；**SPMe 不直接求解 φ_s**（用代数式算 dphi_S）；此函数仅 P2D 调用（src/P2D.jl:31-32）；指向 `docs/thermal_verify/findings.md` |
| (4.8) | md/04 §2.3 液相电荷守恒（完整） `-∂(κ^eff·∂φ_e/∂x + 2(1-t_+^0)·(RT/F)·κ_D^eff·∂ln c_e/∂x)/∂x = -a_s·j` | `src/ElectrolytePotential.jl:1-24 (ElectrolytePotential)` + `src/P2D.jl:109-115 (P2D_charge_BC)` L114 | `kappa_D_eff .*= 2*T*(1-tplus).*dlnf_dlnc(ce_gs).*tau_el./ce_gs.*dcedx_gs` | 🟡 | 完整形式仅 P2D 调用；SPMe 简化为 Ohmic 项 dphi_e（见 4.13）；指向 `docs/thermal_verify/findings.md` |
| (4.9) | md/04 §2.3 液相电荷守恒（简化，忽略扩散电位） | `src/SPMe.jl:148,158 (SPMe_variables!)` | `dphi_S = I_app/3*(L_n/σ_n + L_p/σ_p)`；`dphi_e = ... - I_app*R_EL - dphi_S` | ✅ | md §2.3 列简化形式；SPMe 直接用代数 Ohmic 公式，不求解 PDE |

### §2.4 界面动力学 (Butler-Volmer)

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.10) | md/04 §2.4 过电位 `η = φ_s - φ_e - U(c_s^surf, T)` | `src/SPMe.jl:145-146 (SPMe_variables!)` | `eta_n = 2.0*T*asinh.(j_n/2.0/j0_n_av)`（反解形式） | ⚠️ | md §2.4 给定义式；SPMe 不显式计算 φ_s/φ_e，而是从已知 j 通过 asinh 反推 η；P2D (`src/P2D.jl:175-176, 257-258`) 才用定义式；指向 `docs/thermal_verify/findings.md` |
| (4.11) | md/04 §2.4 交换电流密度 `j_0 = k·c_e^α·(c_s^surf)^α·(c_s,max-c_s^surf)^α` | `src/SPMe.jl:141-142 (SPMe_variables!)` L246-247 (`SPMe_variables`) | `j0_n_gs = NE.k*Arrhenius(Eac_k,T).*abs.(cn_surf.*(1.0 .- cn_surf).*ce_n_gs).^0.5` | ✅ | 对称 α=0.5；c_s,max 已归一化为 1（cs -> cs/cs_max） |
| (4.12) | md/04 §2.4 反应通量 `j = 2·j_0·sinh(F·η/(2·R·T))` | `src/P2D.jl:249-250, 259-260 (P2D_variables)` | `j_n = j0_n.*sinh.(0.5.*eta_n./T)*2.0` | ✅ | md 公式无量纲化后 F/R=1，T*=T；代码直接用 sinh(η/(2T)) |
| — | md/04 §2.4 SPMe 反解 `η = 2·T·asinh(j/(2·j_0))` | `src/SPMe.jl:145-146 (SPMe_variables!)` L250-251 (`SPMe_variables`) | `eta_n = 2.0*T*asinh.(j_n/2.0/j0_n_av)` | ✅ | (4.12) 的解析反函数；j 已知（由 I_app 平均分配），反解 η |

---

## md 04 §3 边界条件

### §3.1 固相电位边界

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.13) | md/04 §3.1 负极集流体 `-σ^eff·∂φ_s/∂x = I_app/A_cell` | `src/P2D.jl:92,160 (P2D_charge_BC/P2D_potentials)` | `flux_ne[1] = 0.0`（直接强制：K_pot[1,1]=-1） | 🟡 | md 给 Neumann；P2D 用直接强制（参考电位）；SPMe 不直接求解；指向 `docs/thermal_verify/findings.md` |
| (4.14) | md/04 §3.1 正极集流体 `σ^eff·∂φ_s/∂x = I_app/A_cell` | `src/P2D.jl:99,161 (P2D_charge_BC/P2D_potentials)` | `flux_pe[end] = Vp0`（迭代校正） | 🟡 | md 给 Neumann；P2D 用 Vp0 迭代参考校正；指向 `docs/thermal_verify/findings.md` |

### §3.2 液相边界

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.15) | md/04 §3.2 端部绝热 `∂c_e/∂x = 0`（端部 Neumann=0） | `src/ElectrolyteDiffusion.jl:1-31 (ElectrolyteDiffusion)` | 自然 BC，装配中不施加端部 Neumann 项 | ✅ | md §3.2 列；代码隐式满足 |
| (4.16) | md/04 §3.2 参考电位 `φ_e|_{x=0} = 0` | `src/P2D.jl:130-135 (P2D_potentials)` | `K_pot[1,:]=0; K_pot[1,1]=-1` | ✅ | md 列"通常取"作参考；代码直接强制（仅 P2D 调用） |

---

## md 04 §4 离散格式

### §4.1 颗粒扩散（有限差分）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.17) | md/04 §4.1 半离散 `M_s·ċ_s = K_s·c_s + f_s` | `src/ElectrodeDiffusion.jl:14,17 (ElectrodeDiffusion)` | `M=Assemble(...Ni,Ni,r²·w·J)`；`K=Assemble(...dNi,dNj,-Ds_eff·r²·w·J)` | ✅ | md 标"有限差分"但代码用 FE 装配；数学等价（Ni·Nj 即集中质量分布） |
| (4.18) | md/04 §4.1 广义 θ 时间离散（默认 θ=0.5 即 Crank-Nicolson） | — | 不在 SPMe 模块内，由外层 ODE solver 控制 | — | 卡点：md §4.1 注"默认 opt.solveType='Crank-Nicolson' θ=0.5"；时间步进由 `src/Solve.jl`/`CrankNicolson.jl` 实现，不在本卷范围 |

### §4.2 电解液扩散（有限元）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.19) | md/04 §4.2 弱形式 `∫N^T·ε_e·ċ_e dΩ - ∫∂N^T·D_e^eff·∂c_e/∂x dΩ = ∫N^T·(1-t_+^0)·a_s·j/F dΩ` | `src/ElectrolyteDiffusion.jl:22,29 (ElectrolyteDiffusion)` | M=∫ε_e·Ni·Nj·w·J；K=-∫D_e^eff·dNi·dNj·w·J | ✅ | md §4.2 公式与代码逐项对应 |
| (4.20) | md/04 §4.2 半离散 `M_ce·ċ_e + K_ce·c_e = f_ce` | `src/ElectrolyteDiffusion.jl:22,29 (ElectrolyteDiffusion)` | M, K 由 Assemble 函数装配 | ✅ | 与 (4.19) 同装配 |
| — | md/04 §4.2 源项 `f_ce = ∫N^T·(1-t_+^0)·a_s·j/F dΩ` | `src/SPMe.jl:193-204 (SPMe_BC)` L199-201 | `coeff[v_ne] .*= (1-tplus)*NE.as*j_n`；`flux_el = Assemble1D(Vi, Ni, coeff, nlen)` | ✅ | md §4.2 公式；代码 SPMe_BC 装配 |

### §4.3 电荷守恒（有限元）

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.21) | md/04 §4.3 固相弱形式 `∫∂N^T·σ^eff·∂φ_s/∂x dΩ = ∫N^T·a_s·j dΩ + 边界通量项` | `src/ElectrodePotential.jl:18 (ElectrodePotential)` + `src/P2D.jl:88-99 (P2D_charge_BC)` | K=∫σ^eff·dNi·dNj·w·J；flux_ne=∫a_s·j·Ni·w·J | 🟡 | md 给完整形式；SPMe 不求解此方程（用代数 dphi_S）；仅 P2D 调用；指向 `docs/thermal_verify/findings.md` |
| (4.22) | md/04 §4.3 半离散（稳态） `K_s·φ_s = f_s` | `src/ElectrodePotential.jl:11 (ElectrodePotential)` | M=spzeros（零矩阵）；仅 K 非零 | ✅ | md §4.3 列稳态；代码 M=0 即稳态 |

---

## md 04 §5 端电压计算

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.23) | md/04 §5 端电压定义 `V_cell = (φ_s-φ_e)|_{x=L} - (φ_s-φ_e)|_{x=0}` | `src/P2D.jl:280 (P2D_variables)` | `variables["cell voltage"] = phis_p[end] - phis_n[1]` | ✅ | md §5 公式；P2D 直接计算；SPMe 用 (4.24) |
| (4.24) | md/04 §5 各分量组成 `V_cell = U_p - U_n + η_p - η_n + Δφ_ohm + Δφ_conc` | `src/SPMe.jl:162 (SPMe_variables!)` L268 (`SPMe_variables`) | `V_cell = u_p - u_n + eta_p - eta_n + dphi_e` | ✅ | md §5 公式；dphi_e 包含 Ohmic + 浓差（见 4.25-4.26） |
| (4.25) | md/04 §5 浓差电位 `Δφ_conc = 2·T·(1-t_+^0)·(c_sp - c_sn)/c_e0` | `src/SPMe.jl:158 (SPMe_variables!)` L264 (`SPMe_variables`) | `dphi_e = 2.0*T*(1-tplus)*(csp_av - csn_av)/ce0 - I_app*R_EL - dphi_S` | ✅ | md §5 公式前半部分；csp_av/csn_av 为电极区域 ce 平均 |
| (4.26) | md/04 §5 欧姆电压降 `Δφ_ohm = I_app·(t_n/(3·σ_n^eff) + t_p/(3·σ_p^eff) + t_n/(3·κ_n^eff) + t_sp/(κ_sp^eff) + t_p/(3·κ_p^eff))` | `src/SPMe.jl:148,155,158 (SPMe_variables!)` L254,261,264 | `dphi_S = I_app/3*(L_n/σ_n + L_p/σ_p)`；`R_EL = L_n/(3·κ_n_av) + L_sp/κ_sp_av + L_p/(3·κ_p_av)`；`dphi_e = ... - I_app*R_EL - dphi_S` | ✅ | md §5 公式拆为 dphi_S（固相）+ R_EL（液相）两部分，代码完全对应 |
| — | md/04 §5 dUdT 修正 `U(T) = U_ref + (T-T_ref)·dUdT` | `src/SPMe.jl:160-161 (SPMe_variables!)` L266-267 (`SPMe_variables`) | `u_n = NE.U(cn_surf) + (T-T0)*NE.dUdT(cn_surf)` | ✅ | md §5 未显式列；代码含熵热修正（用于热源），见卷 05 |

---

## md 04 §6 单元级 SPMe 函数

### §6.1 SPMe_element 函数

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| — | md/04 §6.1 函数签名（5 参数 + 关键字） | `src/SPMe.jl:37-93 (SPMe_element)` | `function SPMe_element(case, yt_e, t, e; I_e, T_e, jacobi="update", workspace=nothing)` | ✅ | md §6.1 签名表与代码逐字一致 |
| — | md/04 §6.1 workspace 可选预分配 | `src/SPMe.jl:39-43 (SPMe_element)` | `if workspace !== nothing: variables_e = SPMe_variables!(workspace, ...)` | ✅ | md §6.1 注释；代码三分支 if |
| — | md/04 §6.1 状态向量结构 `yt_e = [c_s_n; c_s_p; c_e]` | `src/SPMe.jl:131-135 (SPMe_variables!)` 通过 `case.index` 索引 | 索引由 `case.index[...]` 提供（详见[卷 01](./01_参数与归一化对照.md)） | ✅ | md §6.1 列 yt_e 三段结构；代码通过 case.index 分段访问 |
| — | md/04 §6.1 输出 4 元组 `(M_e, K_e, F_e, variables_e)` | `src/SPMe.jl:92 (SPMe_element)` | `return M_e, K_e, F_e, variables_e` | ✅ | md §6.1 表与代码一致 |
| — | md/04 §6.1 jacobi="constant" 时跳过 M/K 重算 | `src/SPMe.jl:62-72 (SPMe_element)` | `if jacobi == "constant" && !isempty(NE.M_d) ...` | ✅ | md §6.1 注释；代码从 param.NE.M_d 复用 |

### §6.2 SPMe_variables 函数

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| — | md/04 §6.2 函数签名（支持 I_app/T_e 覆写） | `src/SPMe.jl:209 (SPMe_variables)` | `function SPMe_variables(case, yt, t; I_app=nothing, T_e=nothing)` | ✅ | md §6.2 签名与代码一致 |
| — | md/04 §6.2 I_app 覆写（默认从 opt.Current 计算） | `src/SPMe.jl:213-217 (SPMe_variables)` L104-108 (`SPMe_variables!`) | `if isnothing(I_app): I_app = opt.Current(t*scale.t0)/scale.I_typ` | ✅ | md §6.2 注释；代码两分支 |
| — | md/04 §6.2 T_e 覆写（默认从状态向量读取） | `src/SPMe.jl:231-235 (SPMe_variables)` L128 (`SPMe_variables!`) | `T = T_e === nothing ? only(yt[index["temperature"]]) : T_e` | ✅ | md §6.2 注释；代码三元运算 |
| — | md/04 §6.2 输出 5+ 类变量（cell voltage, 过电位, etc.） | `src/SPMe.jl:269-289 (SPMe_variables)` L164-179 (`SPMe_variables!`) | 写入 14 个键：cell voltage / 过电位 / 交换电流 / OCP / etc. | ✅ | md §6.2 表 5 行（精简）；代码写入 14 键（含高斯点值） |
| — | md/04 §6.2 原位变体 `SPMe_variables!(ws, ...)` | `src/SPMe.jl:101-182 (SPMe_variables!)` | 函数签名 + 逻辑与 SPMe_variables 逐行对应 | ✅ | md §6.2 注释；代码 L101-182 整函数体 |
| — | md/04 §6.2 精简型工作区 `create_element_workspace(case)` | — | 不在 SPMe.jl | — | 卡点：md §6.2 注"~30 个键，减少 60% 分配"；此函数在 `src/Variables.jl`，不在本卷代码范围 |

### §6.3 温度依赖性

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.27) | md/04 §6.3 Arrhenius 反应速率 `k(T) = k_ref·exp(-E_a/R·(1/T - 1/T_ref))` | `src/SPMe.jl:141-142 (SPMe_variables!)` L246-247 | `j0_n_gs = NE.k*Arrhenius(Eac_k, T).*...` | ✅ | md §6.3 公式；代码调用 `Arrhenius(Eac_k, T)` 工具函数（无量纲化后 R 已吸收进 Eac_k） |
| (4.28) | md/04 §6.3 Arrhenius 扩散系数 `D_s(T) = D_s,ref·exp(-E_a/R·(1/T - 1/T_ref))` | `src/ElectrodeDiffusion.jl:15 (ElectrodeDiffusion)` | `Ds_eff = electrode.Ds*Arrhenius(electrode.Eac_D, T).*(1+(theta.*c))` | ✅ | md §6.3 公式；代码额外乘力学耦合因子 (1+theta·c)（见 §7） |
| — | md/04 §6.3 代码示例（变量名 Eac_D / Eac_k / arrhenius_D / arrhenius_k） | `src/SPMe.jl:141-142 (SPMe_variables!)` + `src/ElectrolyteDiffusion.jl:27` | 实际使用 `Arrhenius(Eac_D, T)` 工具函数（非内联 exp） | ⚠️ | md §6.3 给内联 `arrhenius_D = exp(...)` 示例；代码调用 `Arrhenius()` 函数（封装在 src/Tools.jl）；指向 `docs/thermal_verify/findings.md` |

---

## md 04 §7 力学耦合

### §7.1 应力耦合扩散系数

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.29) | md/04 §7.1 应力耦合系数 `θ_M = 2·E·Ω²/(9·(1-ν)·T)` | — | 不在 SPMe.jl 内 | — | 卡点：md §7.1 公式在 `src/Mechanical.jl::Mechanicaloutput` 中计算；本卷代码范围（SPMe.jl 等 7 文件）不含其实现，详见卷 CZM |
| — | md/04 §7.1 取值路径 `variables["negative particle stress coupling diffusion coefficient"][1]` | `src/SPMe.jl:5-6 (SPMe)` L48-49 (`SPMe_element`) | `theta_Mn = variables_e["negative particle stress coupling diffusion coefficient"][1]` | ✅ | md §7.1 代码示例与代码逐字一致 |
| — | md/04 §7.1 传递给扩散矩阵 `ElectrodeDiffusion(param.NE, mesh_np, ..., theta_Mn)` | `src/SPMe.jl:22 (SPMe)` L70 (`SPMe_element`) | `M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, theta_Mn)` | ✅ | md §7.1 代码示例与代码逐字一致；注意 csn_gs 实参顺序 |

### §7.2 对扩散方程的影响

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| (4.30) | md/04 §7.2 耦合扩散系数 `D_eff = D_s·(1 + θ_M·σ)` | `src/ElectrodeDiffusion.jl:15 (ElectrodeDiffusion)` | `Ds_eff = electrode.Ds*Arrhenius(Eac_D,T).*(1+(theta.*c))` | ⚠️ | md §7.2 公式写 `1+θ_M·σ`；代码用 `1+theta·c`（c 为浓度而非应力）；md 公式或代码符号选择不一致——但代码注释与 md §7.1 路径名"theta_M"一致；指向 `docs/thermal_verify/findings.md` |

---

## md 04 §8 代码位置

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| — | md/04 §8 表行 1 SPMe 主求解 | `src/SPMe.jl:1-35 (SPMe)` | 35 行函数 | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 2 单元级 SPMe | `src/SPMe.jl:37-93 (SPMe_element)` | 57 行函数 | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 3 变量提取 | `src/SPMe.jl:209-287 (SPMe_variables)` | 83 行函数 | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 4 变量提取（原位） | `src/SPMe.jl:101-182 (SPMe_variables!)` | 82 行函数 | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 5 边界条件 | `src/SPMe.jl:184-207 (SPMe_BC)` | 24 行函数 | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 6 颗粒扩散 | `src/ElectrodeDiffusion.jl:1-19 (ElectrodeDiffusion)` | 19 行函数 | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 7 电解液扩散 | `src/ElectrolyteDiffusion.jl:1-31 (ElectrolyteDiffusion)` | 31 行函数 | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 8 电位求解 | `src/ElectrodePotential.jl:1-20 (ElectrodePotential)` | 20 行函数（仅 P2D 调用） | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 9 电解液电位 | `src/ElectrolytePotential.jl:1-24 (ElectrolytePotential)` | 24 行函数（仅 P2D 调用） | ✅ | md §8 表与代码一致 |
| — | md/04 §8 表行 10 力学输出 | `src/Mechanical.jl::Mechanicaloutput` | 不在本卷代码范围 | — | 卡点：md §8 表行 10；详见卷 CZM |
| — | md/04 §8 表行 11 精简型工作区 | `src/Variables.jl::create_element_workspace` | 不在本卷代码范围 | — | 卡点：md §8 表行 11；详见卷 Variables |
| — | md/04 §8 表行 12 变量初始化 | `src/Variables.jl::StandardVariables` | 不在本卷代码范围 | — | 卡点：md §8 表行 12；详见卷 Variables |

---

## md 04 §代码独有 / md 缺失

| 公式编号 | 理论位置 | 代码位置 | 实现摘要 | 一致性 | 备注 |
|---|---|---|---|---|---|
| -- | 【缺失】 | `src/SPM.jl:1-31 (SPM)` | 简化版 SPM（无电解液动力学） | ✅ | md 04 主题为 SPMe；SPM 为退化版；md §1 注 "SPM 是 SPMe 的子集" |
| -- | 【缺失】 | `src/SPM.jl:34-47 (SPM_BC)` | SPM 边界条件（仅颗粒通量，无电解液源项） | ✅ | md 04 未单独描述 SPM BC；与 SPMe_BC 颗粒部分相同 |
| -- | 【缺失】 | `src/SPM.jl:49-82 (SPM_variables)` | SPM 变量提取（无电解液相关项） | ✅ | md 04 未单独描述 SPM variables；与 SPMe_variables 共用核心逻辑（j_n/j_p/eta/V_cell）但不含 dphi_e |
| -- | 【缺失】 | `src/P2D.jl:1-47 (P2D)` | P2D 主求解（求解 4 个 PDE：颗粒扩散 ×2 + 电解液扩散 + 电荷守恒） | ✅ | md 04 §1 注"SPMe 是 P2D 的简化"；P2D 完整形式见 (4.7)-(4.8) |
| -- | 【缺失】 | `src/P2D.jl:49-74 (P2D_mass_BC)` | P2D 质量 BC（含电解液源项 a_s·j_gs） | ✅ | md 04 §2.2 公式 (4.5) 源项；P2D 用高斯点 j_gs 而非平均 j |
| -- | 【缺失】 | `src/P2D.jl:76-118 (P2D_charge_BC)` | P2D 电荷 BC（含扩散电位项 kappa_D） | ✅ | md 04 §2.3 完整形式 (4.8)；代码 L109-115 装配 |
| -- | 【缺失】 | `src/P2D.jl:120-207 (P2D_potentials)` | P2D 电位迭代求解（100 次迭代上限，rel_tol=1e-9） | ⚠️ | md 04 §2.4 给显式 BV；P2D 用朗道迭代（Vp0 校正）求解非线性方程组；指向 `docs/thermal_verify/findings.md` |
| -- | 【缺失】 | `src/P2D.jl:209-283 (P2D_variables)` | P2D 变量提取（含高斯点 j0/eta/j、OCP） | ✅ | md 04 §2.4 BV 公式；P2D 在 gs 上同时计算，供 P2D_charge_BC 使用 |
| -- | 【缺失】 | `src/SPMe.jl:73-75 (SPMe_element) 时间尺度归一化` | `M_np = M_np .* scale.ts_n / param_dim.scale.t0` | ✅ | md §0 归一化方案 t_0=3600；代码在 SPMe 出口将 M 矩阵从 ts（颗粒）尺度转为 t0；详见[卷 01](./01_参数与归一化对照.md) |
| -- | 【缺失】 | `src/SPMe.jl:179 (SPMe_variables!) cell current` | `ws["cell current"] = opt.Current(t*scale.t0)/param_dim.cell.I1C` | ✅ | md 04 §6.2 表行 5 列 "cell current"；代码归一化到 I1C |
| -- | 【缺失】 | `src/SPMe.jl:90 (SPMe_element) element index` | `variables_e["element index"] = Float64(e)` | ✅ | 调试用：md 未列；代码每个 SPMe_element 输出附单元编号 |

---

## 编号总表

| 自编号 | md 出处 | 主题 |
|---|---|---|
| (4.1) | md/04 §2.1 | Fick 球形扩散方程 |
| (4.2) | md/04 §2.1 | 球心对称 BC |
| (4.3) | md/04 §2.1 | 表面流量 BC |
| (4.4) | md/04 §2.1 | 表面浓度定义 |
| (4.5) | md/04 §2.2 | 电解液守恒方程 |
| (4.6) | md/04 §2.2 | 有效扩散系数 D_e^eff |
| (4.7) | md/04 §2.3 | 固相电荷守恒 |
| (4.8) | md/04 §2.3 | 液相电荷守恒（完整） |
| (4.9) | md/04 §2.3 | 液相电荷守恒（简化） |
| (4.10) | md/04 §2.4 | 过电位定义 |
| (4.11) | md/04 §2.4 | 交换电流密度 j_0 |
| (4.12) | md/04 §2.4 | Butler-Volmer sinh |
| (4.13) | md/04 §3.1 | 负极集流体 Neumann |
| (4.14) | md/04 §3.1 | 正极集流体 Neumann |
| (4.15) | md/04 §3.2 | 电解液端部绝热 |
| (4.16) | md/04 §3.2 | 参考电位 |
| (4.17) | md/04 §4.1 | 颗粒半离散 |
| (4.18) | md/04 §4.1 | 广义 θ 时间离散 |
| (4.19) | md/04 §4.2 | 电解液弱形式 |
| (4.20) | md/04 §4.2 | 电解液半离散 |
| (4.21) | md/04 §4.3 | 固相弱形式 |
| (4.22) | md/04 §4.3 | 固相半离散（稳态） |
| (4.23) | md/04 §5 | 端电压定义（P2D） |
| (4.24) | md/04 §5 | 端电压分量组成（SPMe） |
| (4.25) | md/04 §5 | 浓差电位 Δφ_conc |
| (4.26) | md/04 §5 | 欧姆电压降 Δφ_ohm |
| (4.27) | md/04 §6.3 | Arrhenius k(T) |
| (4.28) | md/04 §6.3 | Arrhenius D_s(T) |
| (4.29) | md/04 §7.1 | 应力耦合系数 θ_M |
| (4.30) | md/04 §7.2 | 耦合扩散系数 D_eff |

---

## 草稿

### md 公式清单

**md/04_电化学模型_SPMe.md**（377 行，按 md 章节顺序）：

- §1 (L1-12): 模型概述 + 3 条假设 → 1 条对照（无公式）
- §2.1 (L17-30): Fick 扩散 + BC（球心/表面）+ 表面浓度 → (4.1)-(4.4)
- §2.2 (L32-51): 电解液守恒 + D_e^eff + 区域划分 → (4.5)-(4.6) + 2 条对照
- §2.3 (L53-75): 固相 + 液相电荷守恒（完整/简化） → (4.7)-(4.9)
- §2.4 (L77-95): 过电位 + 交换电流 + 反应通量 → (4.10)-(4.12) + 1 条对照（SPMe 反解）
- §3.1 (L103-113): 固相电位 Neumann BC → (4.13)-(4.14)
- §3.2 (L115-118): 液相绝热 + 参考电位 → (4.15)-(4.16)
- §4.1 (L124-144): 颗粒半离散 + θ 时间离散 → (4.17)-(4.18)
- §4.2 (L146-158): 电解液弱形式 + 半离散 → (4.19)-(4.20) + 1 条源项对照
- §4.3 (L160-172): 固相弱形式 + 稳态 → (4.21)-(4.22)
- §5 (L176-201): 端电压定义 + 分量组成 + 浓差 + Ohm → (4.23)-(4.26) + 1 条 dUdT 对照
- §6.1 (L207-253): SPMe_element 签名 + 状态向量 + 输出 → 5 条对照（无公式）
- §6.2 (L255-288): SPMe_variables + 原位变体 + 精简工作区 → 6 条对照（无公式）
- §6.3 (L290-311): Arrhenius k(T)/D_s(T) → (4.27)-(4.28) + 1 条对照（工具函数封装）
- §7.1 (L317-348): θ_M 公式 + 取值路径 + 传递 → (4.29) + 2 条对照
- §7.2 (L350-358): D_eff 耦合公式 → (4.30)
- §8 (L362-378): 代码位置表 12 行 → 12 条对照（3 行代码不在本卷范围）

合计：30 个自编号公式 + 约 30 条无公式对照条目

### src 函数清单

| 文件 | 行号 | 函数 | 说明 |
|---|---|---|---|
| `src/SPMe.jl` | 1-35 | `SPMe` | 主求解（电解液 + 力学耦合） |
| `src/SPMe.jl` | 37-93 | `SPMe_element` | 单元级（多 SPMe 并行） |
| `src/SPMe.jl` | 101-182 | `SPMe_variables!` | 原位变体 |
| `src/SPMe.jl` | 184-207 | `SPMe_BC` | 边界条件（颗粒通量 + 电解液源项） |
| `src/SPMe.jl` | 209-291 | `SPMe_variables` | 变量提取（创建新 Dict） |
| `src/ElectrodeDiffusion.jl` | 1-19 | `ElectrodeDiffusion` | 颗粒扩散 M/K 装配（含力学耦合因子） |
| `src/ElectrodePotential.jl` | 1-20 | `ElectrodePotential` | 固相电位 K 装配（仅 P2D 用） |
| `src/ElectrolyteDiffusion.jl` | 1-31 | `ElectrolyteDiffusion` | 电解液扩散 M/K 装配 |
| `src/ElectrolytePotential.jl` | 1-24 | `ElectrolytePotential` | 液相电位 K 装配（仅 P2D 用） |
| `src/SPM.jl` | 1-31 | `SPM` | 简化版主求解（无电解液动力学） |
| `src/SPM.jl` | 34-47 | `SPM_BC` | SPM 边界条件 |
| `src/SPM.jl` | 49-86 | `SPM_variables` | SPM 变量提取 |
| `src/P2D.jl` | 1-47 | `P2D` | P2D 主求解 |
| `src/P2D.jl` | 49-74 | `P2D_mass_BC` | 质量 BC（含电解液源项） |
| `src/P2D.jl` | 76-118 | `P2D_charge_BC` | 电荷 BC（含扩散电位） |
| `src/P2D.jl` | 120-207 | `P2D_potentials` | 电位迭代求解（朗道迭代） |
| `src/P2D.jl` | 209-287 | `P2D_variables` | P2D 变量提取 |
