# Mechanical.jl

- **源文件**：`src/Mechanical.jl`
- **行数**：339 行
- **函数/struct 计数**：8 个函数，无 struct
- **职责**：颗粒扩散应力；层分辨宏观热-扩散应力的共享恢复核、耦合在线收割和 CZM-off 固体求解工具。
- **相关技术文档**：`md/06_内聚力模型_CZM.md`、`md/15_颗粒与极片模量区分.md`

## 函数清单

### `Mechanicaloutput(case, variables)` — L1-L110

SPM/SPMe/P2D 颗粒尺度扩散应力与电化学-力学反馈入口。颗粒应力只读取 `Electrode.E/nu/Omega`，不可与涂层宏观模量混用。

### `Calstressdisp(electrode, mesh, cs, T)` — L112-L139

球形颗粒扩散应力解析计算，返回中心径向应力、表面切向应力、表面位移、应力-扩散耦合系数和 Gauss 点浓度。

### `recover_bulk_stress(node, element, material_type, u, ε0, param)` — L151-L173

共享的层分辨 Q4 平面应力恢复核：逐层调用 `moduli_of`，计算 `σ = D(ε-ε₀)` 及 von Mises 应力，返回 `scale.σ_czm` 归一空间结果。

### `macro_eigenstrain(case, variables, T_nodes)` — L182-L189

调用 `compute_czm_strain_inputs` 和 `eigenstrain_of` 生成逐力学 bulk 单元本征应变。PE/NE 分别接收本层 SOC 增量，SP/PCC/NCC 仅接收本层热应变。

### `compute_macro_stress(case, variables, T_nodes)` — L197-L204

以当前已收敛 `case.mech.u_prev` 和同一时间层载荷恢复 `xx/yy/xy/vonMises` 四分量。

### `write_macro_stress!(variables_hist, v, stress)` — L206-L212

把四个归一化应力分量写入第 `v` 个历史列。

### `export_macro_stress(case, variables, variables_hist, v, T_nodes)` — L221-L227

仅当 `case.opt.czm.enabled` 且已分配应力历史键时在线恢复并写入。本步不更新 CZM 时，`Solve.jl` 复用最近一次有效恢复状态，不写伪零。

### `thermal_diffusion_stress_2D(case, variables)` — L247-L339

CZM-off 的按需固体工具：在 `case.czm_mesh.czm_submesh.mesh_bonded` 上按逐层刚度与本征应变求线性平衡，外圈固定、内圈由 `case.opt.czm.fix_inner` 控制。直接线性求解；失败不回退为零位移。输出有量纲应力 `[Pa]` 与位移 `[m]`。

## 参数、状态与缓存契约

- 宏观涂层读取 `PE/NE.E_coat/nu_coat`，SP/PCC/NCC 读取连续层 `E/nu`；不使用已删除的 `E_eff`。
- 界面参数不经过本文件；`moduli_of(param, mt)` 直接读 `param`。
- 耦合恢复读取 `case.mech.u_prev`；不读取 `czm_layout` 或参数缓存。
- `SetCase` 后 `param` 冻结；本文件不维护装配缓存。
