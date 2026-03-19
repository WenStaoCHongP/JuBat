# Parameters_Design 分支修改发现

## 1. 核心架构变更

### 1.1 归一化方案重构
**发现时间**: 2026-03-19

**关键变更**:
- Scale 结构体新增 `rho` 字段用于密度归一化
- 统一能量尺度: `P_ref = φ_ref × I_typ`
- 统一时间尺度: `t_0 = 3600` s (电化学与热模型统一)
- 热源归一化: `q* = q × L³ / P_ref`
- 体积热容归一化: `(ρc)* = ρc × L³ × T_ref / (t_0 × P_ref)`

**影响文件**:
- `src/SetParams.jl`
- `src/ThermalDistributed.jl`
- `src/ThermalPolar2D.jl`

### 1.2 Jellyroll 几何模型
**发现时间**: 2026-03-19

**新增模块**: `src/Jellyrollmodel.jl` (+547 行)

**核心功能**:
- 阿基米德螺旋线描述: `r(θ) = a + bθ`
- collector-seeded 网格生成
- 层序管理: PE → PCC → PE → SP → NE → NCC → NE → SP

### 1.3 多 SPMe 并行架构
**发现时间**: 2026-03-19

**新增模块**: `src/Parallelsolution.jl` (+619 行)

**核心功能**:
- 每个热单元独立的 SPMe 模型
- Newton-Raphson 分流求解
- 状态向量: `yt_global = [yt_chem[1]; ...; yt_chem[ne]; T_nodes]`

### 1.4 CZM 内聚力模型
**发现时间**: 2026-03-19

**新增模块**:
- `src/czm.jl` (+503 行)
- `src/CzmSolve.jl` (+514 行)

**核心功能**:
- 双线性牵引-分离律
- 热-化学载荷耦合
- 损伤状态演化

---

## 2. 热模型修正

### 2.1 热源归一化修正
**提交**: eb2c116, c3aab99

**问题**: heat_Q 归一化缺少 rho 因子

**修复**: 对所有 5 层 (NE, SP, PE, PCC, NCC) 修正热源归一化

### 2.2 时间尺度统一
**提交**: 66b51dd

**变更**: 热模型时间尺度统一为 3600s，与电化学模型一致

### 2.3 几何尺度转换
**提交**: d3a799a

**修复**: 验证脚本中的几何尺度转换，网格坐标归一化后面积需乘以 `scale.L^2`

---

## 3. 验证体系

### 3.1 新增验证脚本
| 脚本 | 用途 |
|------|------|
| thermal_verify.jl | 圆环精确解验证 |
| thermal_error_source_analysis.jl | 热误差来源分析 |
| thermal_equivalent_lumped_compare.jl | 2D 等效集总量对比 |
| jellyroll_vs_ring_thermal_compare.jl | Jellyroll vs Ring 对比 |

### 3.2 PyBaMM 对比数据
新增 `src/data/pybamm_SPMe_LGM50_*.csv` 文件用于 SPMe 验证

---

## 4. 文档体系

### 4.1 技术文档结构
采用四层结构:
- 第一层 (01-03): 参数与基础
- 第二层 (04-07): 模型实现
- 第三层 (08-10): 算法与求解
- 第四层 (11-13): 验证方案

### 4.2 设计文档
`docs/plans/` 目录包含多个实现计划和进度文档

---

## 5. 工具与调试

### 5.1 新增工具脚本 (tools/)
- `check_boundary_nodes.jl`
- `check_branch_currents.jl`
- `check_collector_mesh.jl`
- `check_thermal_kernels.jl`
- `verify_czm_*.jl` 系列

### 5.2 参数导出脚本 (param/)
- `export_lgm50t_spme.py`
- `export_pybamm_params.py`
