# SetParams.jl

- **源文件**：`src/SetParams.jl`
- **行数**：523 行
- **函数/struct 计数**：9 个 struct + 2 个函数
- **职责**：定义有量纲参数结构，按电芯型号加载参数，建立参考尺度并生成无量纲深拷贝。
- **相关技术文档**：`md/01_参数定义与归一化.md`、`md/10_参数传递与模块架构.md`、`md/15_颗粒与极片模量区分.md`

## 数据结构

### `Electrode` — L37-L71

PE/NE 共用。区分颗粒尺度 `E/nu/Omega` 与极片涂层尺度 `E_coat/nu_coat/alphaT`，并包含电化学、传热和运行期颗粒矩阵字段。

### `Separator` — L73-L84

隔膜厚度、传热、孔隙率与连续层力学参数。

### `CurrentCollector` — L86-L115

PCC/NCC 共用，承担三个物理层级：

- 集流体连续层：`thickness/lambda/rho/heat_Q/sig/E/nu/alphaT`
- 集流体塑性：`sigma_y/H`
- 与相邻涂层的 CZM 界面：
  - Mode I：`σ_max/K_n/δ_0/G_c/δ_c`
  - Mode II：`τ_max/K_t/δ_0_t/G_c_t/δ_c_t`
  - 混合模式：`eta`
  - 界面热阻：`h_c0/k_air/lambda_m/beta/threshold`

映射固定为 `PCC ↔ :PE_PCC`、`NCC ↔ :NE_NCC`。共用的 eta/热阻值也由两个实例各自显式保存。

### `Electrolyte` — L117-L127

电解液扩散、电导率、活度、浓度和热参数。

### `Cell` — L129-L160

电芯几何、容量、热学、Jellyroll 半径/导热和卷绕张力参数。

### `Tab` — L163-L170

极耳几何、换热与正负极耳化学计量比。

### `Binder` — L174-L176

粘结剂密度。

### `Scale` — L178-L217

统一电化学、时间、热、颗粒力学、宏观极片和 CZM 参考尺度。CZM 四尺度由 PE-PCC 界面锚定：

```text
σ_czm = PCC.σ_max
δ_czm = 2 * PCC.G_c / PCC.σ_max
G_czm = σ_czm * δ_czm
K_czm = σ_czm / δ_czm
```

### `Params` — L219-L230

聚合 `PE/NE/EL/SP/cell/PCC/NCC/tab/binder/scale`。不存在独立 `cohesive` 字段。

## 函数清单

### `ChooseCell(CellType="LG M50") -> Params` — L232-L335

加载参数文件，补全电极体积分数/比表面积并建立参考尺度。缺少涂层模量或 PCC 锚点时发出明确 `@warn`；PCC 锚点存在时按上式设置 CZM 尺度。

### `NormaliseParam(param_dim) -> param` — L337-L522

`deepcopy(param_dim)` 后逐字段归一化：

- PE/NE 颗粒与涂层模量使用各自参考尺度；SP/PCC/NCC 连续层使用 `scale.E_coat`。
- PCC/NCC 的 `sigma_y/H` 除以 `scale.σ_czm`。
- 每个界面的 `σ_max/τ_max` 除以 `σ_czm`，`K_n/K_t` 除以 `K_czm`，分离量除以 `δ_czm`，断裂能除以 `G_czm`。
- `eta/beta` 保持无因次；`h_c0/k_air/lambda_m/threshold` 转到热模型 `scale.L` 空间。
- 归一化后的字段仍位于 `param.PCC/param.NCC`，求解与装配直接读取。

## 参数冻结契约

`SetCase` 完成归一化后不得修改 `param` 任一字段。`CohesiveMesh.K_bulk` 的新鲜度依赖该契约；如需新参数，应重建 Case/mesh，而不是原位改参或引入内容哈希缓存。

## 已删除接口

独立 `Cohesive` 结构、`Params.cohesive`、`czm_model/tau_visc` 材料字段及参数缓存中转层均已删除。模型与粘性配置由 `opt.czm` 单一持有。
