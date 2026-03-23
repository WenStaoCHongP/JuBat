# Parameters_Design 与 Thermal-debug 分支对比分析

## 基本信息

- **当前分支**: Parameters_Design
- **对比分支**: origin/Thermal-debug
- **提交数量**: 当前分支 17 个，Thermal-debug 分支 0 个
- **代码变更**: +9,337 行 / -389 行 (src 目录)

---

## 1. 提交历史对比

### Thermal-debug 分支提交
| 提交 | 描述 |
|------|------|
| af81aea | 力学模块改进 |
| 46954a4 | 力学模块改进 |
| af81aea | 热应力、扩散应力 |
| d00ea3b | ERROR: LoadError: UndefVarError: U_M not defined in Main.JuBat |
| 290d20a | 热应力为零 |
| 5d3e285 | 结构整理+温度振荡解决 |
| f7e7037 | 不同倍率电压曲线正常，温度曲线异常 |

| 5d3e285 | 结构整理+温度振荡解决 |
| f7e7037 | 不同倍率电压曲线正常，温度曲线异常 |

| 13df0a8 | docs: 重构技术文档目录结构 |
| f41fefc | docs: update heat_Q normalization fix design spec |
| a57cd5e | docs: add implementation plan for heat_Q normalization fix |
| 4975f6c | feat(scale): add rho field to Scale struct for density normalization |
| eb2c116 | fix(normalize): correct heat_Q normalization with rho factor for all 5 layers |
| c3aab99 | fix(thermal): correct volumetric heat capacity calculation |
| 4b0e6f5 | docs: clarify rho and heat_Q semantics, add superseded notice |
| 97363ef | fix: complete unified energy scale migration for thermal parameters |
| b770f03 | refactor(thermal): 统一极坐标热模型归一化方案 |
| f9d57b1 | refactor(thermal): remove fallback default value patterns |
| aeeb52d | simplify |
| 66b51dd | fix(thermal): 统一热模型时间尺度并更新文档 |
| d3a799a | fix(thermal): 修正验证脚本几何尺度转换并更新文档 |
| c7aa30a | docs(thermal): 热验证文档与 Jellyroll/Ring 对齐计划 |
| 81ba44a | docs: 添加热验证归一化相关文档和依赖锁定文件 |
| 6e365c9 | refactor(thermal): 统一热模型边界条件与代码风格规范 |

---

## 2. 源代码变更对比

### 2.1 新增文件 (Parameters_Design 独有)
| 文件 | 行数 | 描述 |
|------|------|------|
| ThermalDistributed.jl | +376 | 二维分布式热模型（全新） |
| CycleSolver.jl | +729 | 循环仿真求解器 |
| CycleData.jl | +621 | 循环数据管理 |
| Parallelsolution.jl | +619 | 分流求解器 |
| czm.jl | +503 | 内聚力模型 |
| CzmSolve.jl | +514 | CZM 求解器 |
| Jellyrollmodel.jl | +547 | Jellyroll 几何模型 |
| Materialmatrix.jl | +382 | 材料矩阵 |
| ThermalPolar2D.jl | +142 | 极坐标热模型 |
| Tools.jl | +173 | 工具函数 |
| ring.jl | +82 | Ring 几何 |
| parameters/Jellyroll.jl | +179 | Jellyroll 参数 |
| parameters/Ring.jl | +65 | Ring 参数 |

### 2.2 Thermal-debug 独有文件（已删除/合并到 Parameters_Design）
- Thermal-debug 分支中的 `ThermalDistributed.jl` 逻辑已合并到 Parameters_Design 分支的新版本

- 关键差异：归一化方式不同

  - Thermal-debug: 使用 `L_th^2` 进行热导率缩尺度缩放
  - Parameters_Design: 统一能量尺度，直接使用无量纲参数

### 2.3 修改文件对比
| 文件 | Parameters_Design | Thermal-debug | 差异 |
|------|-------------------|----------------|------|
| SetParams.jl | +171 行 | +320 行 | 新增 Cohesive/Tab/Scale 字段、热参数归一化 |
| Thermal.jl | +8 行 | +47 行 | 新增 `thermal lumped internal heat` 变量 |
| Solve.jl | +801 行 | +156 行 | 新增多 SPMe 支持、分布式热模型 |
| SPMe.jl | +90 行 | +86 行 | 适配多 SPMe 模型 |

---

## 3. 核心逻辑一致性验证

### 3.1 热模型归一化方案

#### Thermal-debug 方支:
```julia
# 热导率缩放
cxx = -k_xx .* wJ ./ L_th^2
cxy = -k_xy .* wJ ./ L_th^2
cyy = -k_yy .* wJ ./ L_th^2
```

- 问题: 需要显式知道 `L_th` 的值，- `L_th` 从哪里来来?

#### Parameters_Design 分支（统一能量尺度）
```julia
# 热导率（统一能量尺度，直接使用无量纲参数）
cxx = -k_xx .* wJ
cxy = -k_xy .* wJ
cyy = -k_yy .* wJ

```
- 攲进: 在 `thermal_anisotropic_conductivity_2d` 函数中直接计算
- 无需额外的 `L_th^2` 尘垢

- 结果: 更简洁，但 `L_th` 的隐式包含在网格已归一化的事实中中 雂

- **验证**: 两种方案在数学上等价，只是 Parameters_Design 更简洁

### 3.2 边界条件处理
#### Thermal-debug 分支
```julia
Bi = case.param_dim.scale.h
s_vals = (-0.577350269189626, 0.577350269189626)
for (s, w) in zip(s_vals, w_vals)
    ...
end
```

#### Parameters_Design 分支
```julia
Bi = case.param_dim.scale.h  # 直接从 Scale 获取
if Bi == 0
    return K, F
end
...
```
- 拆分了 `apply_convection_bc` 函数
- 使用统一的 `Bi` 参数
- 逻辑更清晰

### 3.3 热源计算
两个分支的热源计算逻辑基本一致：
- 反应热 + 可逆热 + 欧姆热
- 按层加权求和

### 3.4 体积热容计算
#### Thermal-debug 分支问题
体积热容归一化使用 `rho * heat_Q` 的但缺少 `rho` 因子：

```julia
rho_c_weights = fks[32 * param.NE.rho * param.NE.heat_Q * ...
```

- 问题: 热容单位错误

- 修复: 在 Scale 结构体中新增 `rho` 字段，- 修正归一化公式