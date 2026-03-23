# 边界条件传热系数修正计划

**创建日期**: 2026-03-23
**状态**: in_progress

---

## 1. 目标

修正 `apply_convection_bc` 中边界条件传热系数的错误用法，从第一性原理推导并实现正确的无量纲参数。

---

## 2. 问题根因

**错误代码** (`ThermalDistributed.jl:53`)：
```julia
Bi = case.param.cell.h  # ❌ 这是集总模型参数，不是边界条件参数
```

**正确公式**（从第一性原理推导）：
$$h_{conv}^* = Bi \times k_r^* = \text{scale.h} \times \text{param.cell.lambda\_r}$$

---

## 3. 修改计划

### Phase 1: 代码修改 [in_progress]

**文件**: `src/ThermalDistributed.jl`

**修改位置 1**：第53行 `apply_convection_bc` 函数

```julia
# 修正前
Bi = case.param.cell.h  # Biot 数（统一能量尺度）

# 修正后
# 边界条件传热系数：h_conv* = Bi × k_r* = scale.h × param.cell.lambda_r
# 推导：h(T-Tamb) × L²/P_ref = h×T_ref×L²/P_ref × (T*-Tamb*)
h_conv_nd = case.param_dim.scale.h * case.param.cell.lambda_r
```

**修改位置 2**：第88行边界积分权重

```julia
# 修正前
wt = Bi * w * J

# 修正后
wt = h_conv_nd * w * J
```

### Phase 2: 删除错误引用 [pending]

**需要检查和修正的文件**：

| 文件 | 位置 | 错误内容 | 修正 |
|------|------|----------|------|
| `src/ThermalDistributed.jl` | L53 注释 | `# Biot 数（统一能量尺度）` | 删除或修正注释 |
| `md/05_热模型_二维分布式.md` | 相关章节 | 可能存在的错误引用 | 检查并修正 |
| `md/01_参数定义与归一化.md` | 边界条件章节 | 确保公式正确 | 检查一致性 |

**需要保留的正确内容**：
- `scale.h` 是真正的 Biot 数（h·L/λ_r）
- `param.cell.h` 是集总模型传热系数（h·A·T_ref/P_ref）
- 边界条件使用 $h_{conv}^* = Bi \times k_r^*$

### Phase 3: 验证修正 [pending]

**验证脚本**: `example/热模块验证/thermal_equivalent_lumped_compare.jl`

**验收标准**：

| 指标 | 修正前 | 修正后预期 |
|------|--------|------------|
| 边界温度 T_edge | ≈ 298.8 K (接近环境温度) | 接近体积平均温度 |
| 温度梯度 T_vol - T_edge | ≈ 20 K | << 1 K (Bi << 1 时) |
| 散热效率 Pout/Pin | ≈ 15-20% | ≈ 80-90% (接近集总模型) |
| 净热源 (稳态) | ≈ 0.3-0.5 W | ≈ 0 (能量守恒) |

### Phase 4: 文档更新 [pending]

**需要更新的文档**：

1. `md/05_热模型_二维分布式.md`：
   - 添加边界条件无量纲化推导
   - 明确 $h_{conv}^*$ 的定义和物理意义

2. `md/01_参数定义与归一化.md`：
   - 确保边界条件参数描述正确
   - 添加 $h_{conv}^*$ 的说明

3. `docs/thermal_verify/2D模型侧面散热低估问题分析.md`：
   - 更新问题分析和解决方案

---

## 4. 详细推导记录

### 4.1 有量纲边界条件

$$-k \frac{\partial T}{\partial n} = h(T - T_{amb})$$

### 4.2 有限元弱形式

$$\int_{\Gamma} h(T - T_{amb}) N_i \, ds = \int_{\Gamma} h \, T \, N_i \, ds - \int_{\Gamma} h \, T_{amb} \, N_i \, ds$$

### 4.3 无量纲化

$$ds = L \cdot ds^*$$
$$T = T_{ref} \cdot T^*$$

边界热流密度无量纲化（基于能量守恒，使用 $P_{ref}$）：

$$q_s^* = q_s \cdot \frac{L^2}{P_{ref}} = h(T - T_{amb}) \cdot \frac{L^2}{P_{ref}} = \frac{h \cdot T_{ref} \cdot L^2}{P_{ref}}(T^* - T_{amb}^*)$$

### 4.4 边界条件传热系数

$$h_{conv}^* = \frac{h \cdot T_{ref} \cdot L^2}{P_{ref}} = Bi \times k_r^* = \text{scale.h} \times \text{param.cell.lambda\_r}$$

### 4.5 数值验证

- $h_{conv}^* = 0.19$（正确值）
- `param.cell.h` = 3.78×10⁴（错误使用值）
- 比值差异：2×10⁵

---

## 5. 错误记录

| 错误 | 发现时间 | 状态 |
|------|----------|------|
| `apply_convection_bc` 中使用 `param.cell.h` | 2026-03-23 | 待修正 |
| 文档中将 `param.cell.h` 称为 "Biot 数" | 2026-03-23 | 待删除 |

---

## 6. 相关文件

| 文件 | 说明 |
|------|------|
| `src/ThermalDistributed.jl` | 边界条件实现（需修改） |
| `src/SetParams.jl` | 参数归一化（无需修改） |
| `md/01_参数定义与归一化.md` | 归一化文档（检查一致性） |
| `md/05_热模型_二维分布式.md` | 热模型文档（需更新） |
| `docs/thermal_verify/findings.md` | 详细发现记录 |

---

## 7. 执行进度

- [x] Phase 1.1: 从第一性原理推导正确的边界条件参数
- [x] Phase 1.2: 记录详细推导到 findings.md
- [ ] Phase 1.3: 修改 ThermalDistributed.jl 代码
- [ ] Phase 2: 删除文档中的错误引用
- [ ] Phase 3: 运行验证脚本
- [ ] Phase 4: 更新技术文档
