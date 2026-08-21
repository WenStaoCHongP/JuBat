# tools/theory_geometry_recompute.jl
#
# 从代码实参重算 Jellyroll 螺旋几何，供 Theory/01,02,04 几何表引用（D14）。
# 只读诊断脚本：不改任何求解路径，不写文件。
#
# 运行方式: julia --startup-file=no --project=. tools/theory_geometry_recompute.jl
#
# 每一行输出都标注对应代码字段名，Theory 表格必须逐行引用该字段名，
# 使几何数值只有一个来源，避免再次漂移。

using Printf

include(joinpath(@__DIR__, "../src/JuBat.jl"))
using .JuBat

p = JuBat.ChooseCell("Jellyroll")

layer = p.cell.layer
Rin = p.cell.Rin
Rout = p.cell.Rout
b = layer / (2pi)

# theta 上界与 Jellyrollmodel.jl:40 的 theta1 完全一致（s_in=0, s_out=cell.layer）
theta_end = min((Rout - Rin - layer) / b, (Rout - Rin) / b)
N_turns = theta_end / (2pi)
# 阿基米德螺旋弧长，b << r 下取 ∫ r dθ
L_spiral = Rin * theta_end + b * theta_end^2 / 2
r_avg = L_spiral / theta_end
L_turn = 2pi * r_avg
r_out_mid = Rin + b * theta_end

println("=== 层厚（src/parameters/Jellyroll.jl） ===")
@printf("PE.thickness   = %6.2f um\n", p.PE.thickness * 1e6)
@printf("NE.thickness   = %6.2f um\n", p.NE.thickness * 1e6)
@printf("SP.thickness   = %6.2f um\n", p.SP.thickness * 1e6)
@printf("PCC.thickness  = %6.2f um\n", p.PCC.thickness * 1e6)
@printf("NCC.thickness  = %6.2f um\n", p.NCC.thickness * 1e6)

println("\n=== 螺旋几何 ===")
@printf("cell.layer (t_repeat)      = %.6e m = %.2f um\n", layer, layer * 1e6)
@printf("cell.Rin   (a)             = %.6e m = %.3f mm\n", Rin, Rin * 1e3)
@printf("cell.Rout                  = %.6e m = %.3f mm\n", Rout, Rout * 1e3)
@printf("b = cell.layer/(2pi)       = %.6e m/rad\n", b)
@printf("theta_end (Jellyrollmodel.jl:40) = %.4f rad\n", theta_end)
@printf("N_turns = theta_end/(2pi)  = %.3f\n", N_turns)
@printf("L_spiral = int r dtheta    = %.5f m\n", L_spiral)
@printf("r_avg = L_spiral/theta_end = %.6e m = %.3f mm\n", r_avg, r_avg * 1e3)
@printf("L_turn = 2pi*r_avg         = %.6e m = %.2f mm\n", L_turn, L_turn * 1e3)
@printf("L_turn/L_spiral            = %.3f %% (= 1/N_turns = %.3f %%)\n",
    100 * L_turn / L_spiral, 100 / N_turns)

println("\n=== 螺距角 gamma = t_repeat/(2*pi*r) ===")
for (nm, r) in (("r = Rin  ", Rin), ("r = r_avg", r_avg), ("r = r_out", r_out_mid))
    g = layer / (2pi * r)
    @printf("%s : gamma = %.5e rad = %.4f deg, gamma^2 = %.3e\n", nm, g, g * 180 / pi, g^2)
end

println("\n=== 周向单元长与长宽比（每层厚度方向仅 1 个 Q4）===")
for nth in (40, 80, 360)
    Lc = L_turn / nth
    @printf("n_theta=%3d: elem_len=%7.2f um | SP/NCC %5.1f:1 | PCC %5.1f:1 | PE %4.1f:1 | NE %4.1f:1\n",
        nth, Lc * 1e6, Lc / p.SP.thickness, Lc / p.PCC.thickness,
        Lc / p.PE.thickness, Lc / p.NE.thickness)
end
