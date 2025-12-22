# 力学参数导出和使用指南

## 快速开始

### 1. 安装 PyBaMM（如果还没有安装）

```bash
pip install pybamm
```

### 2. 导出 LG M50 的力学参数

```bash
python export_lgm50_mechanical_params.py
```

### 3. 查看导出结果

导出文件位于 `JuBat/output/` 目录：

```bash
# 查看 CSV 格式
cat JuBat/output/lgm50_mechanical_params.csv

# 查看 JSON 格式（详细）
cat JuBat/output/lgm50_mechanical_params.json

# 查看 Julia 代码片段
cat JuBat/output/lgm50_mechanical_params.jl
```

### 4. 更新 JuBat 参数文件

将 `lgm50_mechanical_params.jl` 中的参数复制到：
```
src/parameters/LGM50.jl
```

或者直接使用脚本中的参数（已经在当前 LGM50.jl 中更新）。

## 关键问题解答

### Q1: PyBaMM 中有化学膨胀系数吗？

**A**: 有，但不是直接以 `beta_c` 的形式。PyBaMM 使用**偏摩尔体积 (Ω)**。

- 参数名：`"Negative/Positive electrode partial molar volume [m3.mol-1]"`
- 单位：m³/mol
- 物理意义：每插入 1 mol 锂引起的体积变化

### Q2: 如何从 Ω 计算 beta_c？

**A**: 使用公式：

```
β_c = |Ω × c_s_max|
```

其中：
- Ω: 偏摩尔体积 [m³/mol]
- c_s_max: 最大锂浓度 [mol/m³]
- β_c: 化学膨胀系数（无量纲）

### Q3: LG M50 的典型值是多少？

**A**: 从 PyBaMM Chen2020 参数集计算：

| 参数 | 负极（石墨） | 正极（NMC811） |
|------|-------------|---------------|
| Ω | 3.1 × 10⁻⁶ m³/mol | -7.28 × 10⁻⁷ m³/mol |
| c_s_max | 33,133 mol/m³ | 63,104 mol/m³ |
| **β_c** | **0.103 (10.3%)** | **0.046 (4.6%)** |

这些值与文献报道一致：
- 石墨：10-13% 体积变化
- NMC：2-5% 体积变化

## 文件说明

### 新增文件

1. **export_lgm50_mechanical_params.py** ⭐
   - 从 PyBaMM 导出力学参数的主脚本
   - 计算化学膨胀系数
   - 生成多种格式的输出

2. **docs/PyBaMM_Chemical_Expansion_Parameters.md**
   - 详细解释 PyBaMM 中的化学膨胀参数
   - 理论推导和公式
   - 与文献的对比

3. **MECHANICAL_PARAMETERS_GUIDE.md** (本文件)
   - 快速使用指南

### 输出文件（运行脚本后生成）

1. **JuBat/output/lgm50_mechanical_params.csv**
   - 表格格式，易于查看
   
2. **JuBat/output/lgm50_mechanical_params.json**
   - 结构化数据，包含完整信息和说明
   
3. **JuBat/output/lgm50_mechanical_params.jl**
   - Julia 代码片段，可直接使用

### 已更新文件

1. **src/parameters/LGM50.jl**
   - 已添加 `alphaT` 和 `beta_c` 参数
   - 可以运行导出脚本后进一步验证/更新

## 使用示例

### Python 端（PyBaMM）

```python
import pybamm

# 加载参数
param = pybamm.ParameterValues("Chen2020")

# 获取偏摩尔体积
Omega_n = param["Negative electrode partial molar volume [m3.mol-1]"]
Omega_p = param["Positive electrode partial molar volume [m3.mol-1]"]

print(f"负极 Ω = {Omega_n:.6e} m³/mol")
print(f"正极 Ω = {Omega_p:.6e} m³/mol")

# 计算化学膨胀系数
cs_max_n = param["Maximum concentration in negative electrode [mol.m-3]"]
cs_max_p = param["Maximum concentration in positive electrode [mol.m-3]"]

beta_c_n = abs(Omega_n * cs_max_n)
beta_c_p = abs(Omega_p * cs_max_p)

print(f"负极 β_c = {beta_c_n:.3f} ({beta_c_n*100:.1f}% 体积变化)")
print(f"正极 β_c = {beta_c_p:.3f} ({beta_c_p*100:.1f}% 体积变化)")
```

### Julia 端（JuBat）

```julia
using JuBat

# 加载参数（已包含 beta_c）
param = JuBat.ChooseCell("LG M50")

# 查看化学膨胀系数
println("负极 β_c = ", param.NE.beta_c)
println("正极 β_c = ", param.PE.beta_c)

# 使用宏观扩散应力计算
case = JuBat.SetCase(param, opt)
variables = JuBat.diffusion_stress_2D(case, variables)
```

## 高级用法

### 列出所有力学相关参数

```bash
python export_lgm50_mechanical_params.py --list
```

这会列出 PyBaMM 中所有包含以下关键词的参数：
- molar, volume, expansion
- modulus, poisson
- thermal, mechanical
- stress, strain, elastic

### 查看 PyBaMM 所有参数

```bash
python -c "import pybamm; p=pybamm.ParameterValues('Chen2020'); print('\n'.join(sorted(p.keys())))"
```

### 在 Julia 中验证参数

```julia
# 检查体积变化是否合理
function check_volume_change(electrode)
    beta_c = electrode.beta_c
    Omega = electrode.Omega
    cs_max = electrode.cs_max
    
    # 从 Omega 重新计算
    beta_c_calc = abs(Omega * cs_max)
    
    println("参数文件中的 β_c: ", beta_c)
    println("从 Ω 计算的 β_c: ", beta_c_calc)
    println("体积变化: ", beta_c * 100, "%")
end

check_volume_change(param.NE)
check_volume_change(param.PE)
```

## 参数说明

### 已实现的参数

| 参数 | 符号 | 单位 | 来源 | 状态 |
|------|------|------|------|------|
| 偏摩尔体积 | Ω | m³/mol | PyBaMM | ✅ 已导出 |
| 最大浓度 | c_s_max | mol/m³ | PyBaMM | ✅ 已有 |
| 化学膨胀系数 | β_c | - | 计算 | ✅ 已添加 |
| 热膨胀系数 | α | 1/K | 文献 | ✅ 已添加 |
| 弹性模量 | E | Pa | 文献 | ✅ 已有 |
| 泊松比 | ν | - | 文献 | ✅ 已有 |

### 参数来源

**从 PyBaMM 直接获取**：
- Ω (偏摩尔体积)
- c_s_max (最大锂浓度)

**通过计算获得**：
- β_c = |Ω × c_s_max|

**从文献获取**：
- E (弹性模量)
- ν (泊松比)
- α (热膨胀系数)

## 验证方法

### 1. 与文献对比

| 材料 | 实验测量 | PyBaMM计算 | 匹配 |
|------|---------|-----------|------|
| 石墨 | 10-13% | ~10% | ✅ |
| NMC | 2-5% | ~5% | ✅ |

### 2. 与 PyBaMM 应力模拟对比

```python
import pybamm

# 运行带应力的模拟
model = pybamm.lithium_ion.SPM(
    options={"particle mechanics": "swelling only"}
)

sim = pybamm.Simulation(model, parameter_values="Chen2020")
sim.solve([0, 3600])

# 提取应力
stress = sim.solution["X-averaged negative particle surface tangential stress"]
```

### 3. 单位换算检查

```python
# 验证单位
Omega = 3.1e-6  # m³/mol
cs_max = 33133  # mol/m³

# 体积应变（无量纲）
epsilon = Omega * cs_max
print(f"体积应变 = {epsilon:.4f} = {epsilon*100:.2f}%")

# 应该得到约 0.1 或 10%
```

## 常见问题

### Q: 为什么负极的 β_c 这么大（10%）？

**A**: 石墨的层间距在锂化时会显著增加，这是石墨的固有特性。实验测量值为 10-13%，我们的计算值 10.3% 在合理范围内。

### Q: 正极的 Ω 为什么是负值？

**A**: 负值表示在锂化过程中晶格收缩（虽然整体还是膨胀，但相对于锂金属基准是收缩）。我们在计算 β_c 时取绝对值，因为关心的是体积变化的大小。

### Q: 为什么不直接从 PyBaMM 获取 E、ν、α？

**A**: PyBaMM 主要关注电化学，Chen2020 参数集不包含这些纯力学参数。需要从力学文献补充。推荐值：
- 石墨 E = 10-15 GPa
- NMC E = 100-200 GPa

### Q: 这些参数对应的是单晶还是复合电极？

**A**: 这些是**活性材料**的本征参数。实际复合电极（含粘结剂、导电剂、孔隙）的有效模量会更低，可以用：
```
E_eff = E_active * (1 - ε)^n
```
其中 ε 是孔隙率，n ≈ 2-3。

## 相关文档

1. **理论推导**
   - `docs/Diffusion_Stress_Macroscale_Theory.md`
   - 应力守恒方程和本构关系

2. **PyBaMM 参数详解**
   - `docs/PyBaMM_Chemical_Expansion_Parameters.md`
   - Ω 到 β_c 的转换

3. **使用指南**
   - `docs/Diffusion_Stress_Usage.md`
   - 如何在 JuBat 中使用

4. **实现总结**
   - `docs/Diffusion_Stress_Implementation_Summary.md`
   - 完整实现概览

## 参考文献

1. Chen et al. (2020), "Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models", J. Electrochem. Soc.

2. Christensen & Newman (2006), "Stress generation and fracture in lithium insertion materials", J. Solid State Electrochem.

3. Zhang et al. (2007), "Numerical simulation of intercalation-induced stress in Li-ion battery electrode particles", J. Electrochem. Soc.

---

**更新日期**: 2025-12-22  
**版本**: 1.0  
**状态**: 完成 ✅
