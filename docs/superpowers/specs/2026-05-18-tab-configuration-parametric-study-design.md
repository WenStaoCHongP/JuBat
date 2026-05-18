# 极耳工况参数研究设计

> 日期: 2026-05-18
> 分支: czm-refactor
> 状态: Draft

---

## 1. 研究目标

研究正极耳数量和周向位置对 Jellyroll 电池热-CZM 耦合响应的影响，采用控制变量法设计正交实验矩阵。

**核心问题**：
- 极耳数量（1/2/3个）如何影响温度场分布和损伤演化？
- 极耳周向位置如何影响局部热应力集中和脱粘起始？
- 双极耳配置下，极耳间距如何影响温度均匀性和损伤模式？

## 2. 控制变量法约束

### 2.1 固定参数

| 参数 | 值 | 说明 |
|------|----|------|
| `tab.width` | 4e-3 m | 极耳宽度 |
| `tab.length` | 0.75 * 9.9e-3 m | 极耳长度 |
| `tab.area` | width * length * 2 | 双面面积 |
| `tab.h` | 10 W/(m²·K) | 极耳换热系数 |
| `theta_neg` | `[44π]` | 负极耳位置（所有工况固定） |
| `cool_method` | `"tab"` | 冷却方式 |
| `theta_neg` 位置选择依据 | 负极耳在外螺旋区域，远离正极耳以最小化热耦合 | 控制变量约束说明 |
| 电化学参数 | Jellyroll 默认 | SPMe 模型参数 |
| 热参数 | Jellyroll 默认 | 各向异性导热 |
| CZM 参数 | 当前分支默认 | 双线性牵引-分离律 |

### 2.2 运行时参数（所有工况固定）

| 参数 | 值 | 说明 |
|------|----|------|
| `opt.model` | `"SPMe"` | 电化学模型 |
| `opt.thermal_enabled` | `true` | 启用热耦合 |
| `opt.thermalmodel` | `"distributed2D"` | 二维分布式热模型 |
| `opt.per_element_spme` | `true` | 逐单元 SPMe |
| `opt.czm_enabled` | `true` | 启用 CZM |
| `opt.mechanicalmodel` | `"full"` | 完整力学模型 |
| `opt.Current` | `x -> 5.0` | 1C 恒流放电 |
| `opt.time` | `[0, 3600]` | 单次放电 |
| `opt.cool_method` | `"tab"` | 极耳冷却 |
| `nθ` | 80 | 热网格周向分辨率 |
| `gsorder` | 2 | 高斯积分阶数 |
| 循环参数 | 按需 | 多循环工况另行定义 |

> 注：`n_fractured` 在单次放电中可能为零（D < 1.0），需多循环仿真才能有效对比。

### 2.3 变量

- **正极耳数量** `n_pos`: 1, 2, 3
- **正极耳周向位置** `theta_pos`: 向量，元素为角度 [rad]

## 3. 工况矩阵

正交分组法，3组实验，9个标注工况（其中7个独立）。

### 3.1 第一组：极耳数量效应

固定位置策略为**等间距分布**，仅改变极耳数量。

| Case | n_pos | theta_pos | 说明 |
|------|-------|-----------|------|
| 1 | 1 | `[θ_mid]` | 基准：正极区域中间角度 |
| 2 | 2 | `[θ_mid - Δθ, θ_mid + Δθ]` | 等间距双极耳 |
| 3 | 3 | `[θ_1, θ_mid, θ_3]` | 等间距三极耳 |

其中：
- `θ_mid` = 正极区域（内螺旋）中间角度
- `Δθ` = 正极区域角度范围的 1/4
- `θ_1`, `θ_3` = 正极区域两端各 1/6 处

**对比分析**：Case 1 vs 2 vs 3 → 极耳数量的主效应

### 3.2 第二组：极耳位置效应

数量固定为 **1个正极耳**，改变周向位置。

| Case | n_pos | theta_pos | 说明 |
|------|-------|-----------|------|
| 4 | 1 | `[θ_start]` | 正极区域起始位置（内层） |
| 5 | 1 | `[θ_mid]` | 中间位置（= Case 1） |
| 6 | 1 | `[θ_end]` | 正极区域末端位置（外层） |

其中：
- `θ_start` ≈ 内螺旋起点附近角度
- `θ_end` ≈ 内螺旋终点附近角度

**对比分析**：Case 4 vs 5 vs 6 → 极耳位置的主效应

### 3.3 第三组：双极耳间距效应

数量固定为 **2个正极耳**，改变间距。

| Case | n_pos | theta_pos | 说明 |
|------|-------|-----------|------|
| 7 | 2 | `[θ_mid - Δθ/2, θ_mid + Δθ/2]` | 紧凑间距 |
| 8 | 2 | `[θ_mid - Δθ, θ_mid + Δθ]` | 中等间距（= Case 2） |
| 9 | 2 | `[θ_start, θ_end]` | 最大间距（覆盖正极区域） |

**对比分析**：Case 7 vs 8 vs 9 → 极耳间距的主效应

### 3.4 工况交叉关系

```
独立工况: 1, 2, 3, 4, 6, 7, 9 (共7个)
重复工况: Case 1 = Case 5, Case 2 = Case 8
```

## 4. 具体角度参数

> **归一化说明**：`theta_pos` 值是无量纲的累积弧度角，归一化前后保持不变（弧度是几何不变量）。
> 角度范围 `θ_min_in` / `θ_max_in` 必须使用**归一化参数**（`param.cell.Rin`、`param.cell.layer`）从网格计算。
> 参见 `src/Jellyrollmodel.jl` 中 `jellyroll_tab_node_indices` 的 `theta_cum_in` 计算。

### 4.1 角度计算代码片段

```julia
# 从网格数据提取内螺旋角度范围
# mesh_data = jellyroll_collector_seed_mesh(param_dim; nθ=80, gsorder=2)
# param = NormaliseParam(param_dim)  # 获取归一化参数
mesh = mesh_data.mesh  # 网格对象
nn = size(mesh.node, 1)

a = param.cell.Rin                          # 归一化内半径
b = param.cell.layer / (2 * pi)             # 归一化螺旋增长率

# 计算内螺旋累积角度
theta_cum_in = [(hypot(mesh.node[i,1], mesh.node[i,2]) - a) / b for i in 1:nn]
theta_min_in, theta_max_in = extrema(theta_cum_in)

# 定义工况角度
range_in = theta_max_in - theta_min_in
θ_mid   = (theta_min_in + theta_max_in) / 2
Δθ      = range_in / 4
θ_start = theta_min_in + range_in * 0.1
θ_end   = theta_min_in + range_in * 0.9
θ_1     = theta_min_in + range_in / 6
θ_3     = theta_max_in - range_in / 6
```

### 4.2 角度标记汇总

| 角度标记 | 计算公式 | 参考值（基于 15π 附近） |
|----------|----------|------------------------|
| `θ_start` | `θ_min_in + (θ_max_in - θ_min_in) * 0.1` | 起始 10% 处 |
| `θ_mid` | `(θ_min_in + θ_max_in) / 2` | 中间位置 |
| `θ_end` | `θ_min_in + (θ_max_in - θ_min_in) * 0.9` | 末端 90% 处 |
| `Δθ` | `(θ_max_in - θ_min_in) / 4` | 范围的 1/4 |

## 5. 输出指标

### 5.1 CZM 损伤指标（主指标）

| 指标 | Key | 说明 |
|------|-----|------|
| 最大损伤 | `D_max` | 最大损伤值（0-1） |
| 平均损伤 | `D_mean` | 所有 cohesive 单元平均损伤 |
| 断裂数量 | `n_fractured` | D = 1 的单元数量 |
| 健康状态 | `soh` | 电池结构健康度 |
| 损伤空间分布 | 各节点 D 场 | 极坐标损伤云图 |

### 5.2 热场指标

| 指标 | Key | 说明 |
|------|-----|------|
| 最高温度 | `T_max` | 节点温度最大值 [K] |
| 平均温度 | `T_mean` | 体积加权平均温度 [K] |
| 最大温差 | `ΔT = T_max - T_min` | 温度不均匀度 [K] |
| 温度场 | 节点 T 场 | 极坐标温度云图 |

### 5.3 电化学性能指标

| 指标 | Key | 说明 |
|------|-----|------|
| 电压曲线 | `cell voltage [V]` | 端电压 vs 时间 |
| 容量损失 | `Q_loss` | 相对初始容量的损失 |
| 电流分布 | `element current` | 各单元电流密度分布 |

## 6. 分析方法

### 6.1 组内对比

每组内仅改变一个变量：
- 第一组（Case 1, 2, 3）：控制变量为极耳数量
- 第二组（Case 4, 5, 6）：控制变量为极耳位置
- 第三组（Case 7, 8, 9）：控制变量为极耳间距

### 6.2 跨组效应分析

通过正交分解评估各变量的主效应大小：
- 数量效应幅度 = `|max(第一组指标) - min(第一组指标)|`
- 位置效应幅度 = `|max(第二组指标) - min(第二组指标)|`
- 间距效应幅度 = `|max(第三组指标) - min(第三组指标)|`

### 6.3 可视化规范

| 图表编号 | 类型 | 内容 |
|----------|------|------|
| fig_tab1 | 多子图 | 第一组：不同极耳数量下的损伤/温度场对比 |
| fig_tab2 | 多子图 | 第二组：不同极耳位置下的损伤/温度场对比 |
| fig_tab3 | 多子图 | 第三组：不同极耳间距下的损伤/温度场对比 |
| fig_tab4 | 汇总图 | 三组工况的 D_max, T_max, ΔT, SOH 柱状对比 |

所有场图采用极坐标投影，标注极耳位置。

## 7. Julia 配置模板

```julia
# 极耳工况参数研究模板
# 在 param/Jellyroll.jl 基础上修改 theta_pos

# --- 固定参数（所有工况不变）---
tab.width = 4e-3
tab.length = 0.75 * 9.9e-3
tab.area = tab.width * tab.length * 2
tab.h = 10.0
tab.theta_neg = [44π]

# --- 第一组：极耳数量效应 ---
# Case 1: 单极耳
tab.theta_pos = [θ_mid]

# Case 2: 双极耳等间距
tab.theta_pos = [θ_mid - Δθ, θ_mid + Δθ]

# Case 3: 三极耳等间距
tab.theta_pos = [θ_1, θ_mid, θ_3]

# --- 第二组：极耳位置效应 ---
# Case 4: 起始位置
tab.theta_pos = [θ_start]

# Case 5: 中间位置 (同 Case 1, theta_pos = [θ_mid])
tab.theta_pos = [θ_mid]

# Case 6: 末端位置
tab.theta_pos = [θ_end]

# --- 第三组：双极耳间距效应 ---
# Case 7: 紧凑间距
tab.theta_pos = [θ_mid - Δθ/2, θ_mid + Δθ/2]

# Case 9: 最大间距
tab.theta_pos = [θ_start, θ_end]
```

## 8. 技术文档 `md/15_参数研究_极耳工况.md`

设计完成后将创建 `md/15_参数研究_极耳工况.md` 作为正式技术文档。该文档属于第五层（参数研究），编号 15 遵循现有 `md/` 目录编号序列（当前最大为 14_粘性正则化.md）。

---

## 附录：与已有研究的关系

| 已有研究 | 内容 | 与本研究的关系 |
|----------|------|----------------|
| `param/fig8_single_tab_debonding.py` | 单极耳不同位置损伤演化（后处理可视化） | 本研究 Case 4/5/6 的完整仿真版本 |
| `param/fig9_triple_tab_debonding.py` | 三极耳损伤演化（后处理可视化） | 本研究 Case 3 的完整仿真版本 |
| 网格敏感性分析 | 热网格和CZM网格收敛性 | 为本研究提供网格参数选择依据 |
