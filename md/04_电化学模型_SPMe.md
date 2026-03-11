# 电化学模型 (SPMe)

## 1. 模型概述

SPMe (Single Particle Model with electrolyte) 是一种简化的电池电化学模型，它在保持核心物理过程的同时降低了计算复杂度。

**主要假设**：

- 电解液浓度在空间上均匀分布（1D 平均）
- 每个电极区域使用代表性球形颗粒
- 包含电解液动力学

---

## 2. 控制方程

### 2.1 颗粒内扩散

在球形颗粒坐标系 (\(r \in [0, R_s]\)) 中，锂浓度 \(c_s(r,t)\) 满足 Fick 扩散方程：

$$
\frac{\partial c_s}{\partial t} = \frac{1}{r^2} \frac{\partial}{\partial r} \left( D_s \cdot r^2 \frac{\partial c_s}{\partial r} \right)
$$

**边界条件**：

- 球心 (\(r=0\))：\(\frac{\partial c_s}{\partial r}\big|_{r=0} = 0\) （对称性）
- 表面 (\(r=R_s\))：\(-D_s \frac{\partial c_s}{\partial r}\big|_{r=R_s} = \frac{j}{F}\)

**表面浓度**：\(c_s^{surf}(t) = c_s(r=R_s, t)\)

### 2.2 电解液物质守恒

在厚度方向 (\(x \in [0, L]\)) 中，电解液浓度 \(c_e(x,t)\) 满足：

$$
\frac{\partial (\varepsilon_e c_e)}{\partial t} = \frac{\partial}{\partial x} \left( D_e^{eff} \frac{\partial c_e}{\partial x} \right) + \frac{(1-t_+^0) a_s j}{F}
$$

其中：

- \(\varepsilon_e\)：孔隙率
- \(D_e^{eff} = D_e \cdot \varepsilon_e^{brugg}\)：有效扩散系数
- \(a_s\)：比表面积（隔膜内为 0）
- \(j\)：界面反应通量

**区域划分**：

- 负极区：\(\Omega_n = [0, L_n]\)
- 隔膜区：\(\Omega_s = [L_n, L_n+L_s]\)，\(a_s = 0, j = 0\)
- 正极区：\(\Omega_p = [L_n+L_s, L]\)

### 2.3 电荷守恒

#### 固相电荷守恒

$$
-\frac{\partial}{\partial x} \left( \sigma^{eff} \frac{\partial \phi_s}{\partial x} \right) = a_s j
$$

其中 \(\sigma^{eff} = \sigma \cdot \varepsilon^{brugg}\) 为有效固相电导率。

#### 液相电荷守恒

完整形式（包含扩散电位）：

$$
-\frac{\partial}{\partial x} \left( \kappa^{eff} \frac{\partial \phi_e}{\partial x} + 2(1-t_+^0) \frac{RT}{F} \kappa_D^{eff} \frac{\partial \ln c_e}{\partial x} \right) = -a_s j
$$

常用简化形式（忽略扩散电位项）：

$$
-\frac{\partial}{\partial x} \left( \kappa^{eff} \frac{\partial \phi_e}{\partial x} \right) = -a_s j
$$

### 2.4 界面动力学 (Butler-Volmer)

**过电位定义**：

$$
\eta = \phi_s - \phi_e - U(c_s^{surf}, T)
$$

**交换电流密度**：

$$
j_0 = k \cdot c_e^{\alpha} \cdot (c_s^{surf})^{\alpha} \cdot (c_{s,max} - c_s^{surf})^{\alpha}
$$

**反应通量**（对称情况 \(\alpha_a = \alpha_c = 0.5\)）：

$$
j = 2 j_0 \sinh\left( \frac{F\eta}{2RT} \right)
$$

---

## 3. 边界条件

### 3.1 固相电位边界

- **负极集流体 (\(x=0\))**：

$$
-\sigma^{eff} \frac{\partial \phi_s}{\partial x} = \frac{I_{app}}{A_{cell}}
$$

- **正极集流体 (\(x=L\))**：

$$
\sigma^{eff} \frac{\partial \phi_s}{\partial x} = \frac{I_{app}}{A_{cell}}
$$

### 3.2 液相边界

- **端部绝热**：\(\frac{\partial c_e}{\partial x} = 0\)
- **参考电位**：通常取 \(\phi_e|_{x=0} = 0\) 或设为平均值

---

## 4. 离散格式

### 4.1 颗粒扩散（有限差分）

**半离散形式**：

$$
M_s \dot{c}_s = K_s c_s + f_s
$$

其中：

- \(M_s\)：质量矩阵（含 \(1/r^2\) 因子）
- \(K_s\)：扩散刚度矩阵
- \(f_s\)：边界通量项

**时间离散**（后退欧拉，\(\theta=1\)）：

$$
\left( \frac{M_s}{\Delta t} - K_s \right) c_s^{n+1} = \frac{M_s}{\Delta t} c_s^n + f_s^{n+1}
$$

### 4.2 电解液扩散（有限元）

**弱形式**：

$$
\int N^T \varepsilon_e \dot{c}_e \, d\Omega - \int \frac{\partial N}{\partial x}^T D_e^{eff} \frac{\partial c_e}{\partial x} \, d\Omega = \int N^T \frac{(1-t_+^0) a_s j}{F} \, d\Omega
$$

**半离散形式**：

$$
M_{ce} \dot{c}_e + K_{ce} c_e = f_{ce}
$$

### 4.3 电荷守恒（有限元）

**固相弱形式**：

$$
\int \frac{\partial N}{\partial x}^T \sigma^{eff} \frac{\partial \phi_s}{\partial x} \, d\Omega = \int N^T a_s j \, d\Omega + \text{边界通量项}
$$

**半离散形式**（稳态）：

$$
K_s \phi_s = f_s
$$

---

## 5. 端电压计算

**端电压定义**：

$$
V_{cell} = (\phi_s - \phi_e)|_{x=L} - (\phi_s - \phi_e)|_{x=0}
$$

**各分量组成**：

$$
V_{cell} = U_p - U_n + \eta_p - \eta_n + \Delta\phi_{ohm}
$$

其中：

- \(U_p - U_n\)：开路电压差
- \(\eta_p - \eta_n\)：过电位差（反应过电位）
- \(\Delta\phi_{ohm}\)：欧姆电压降（固相+液相）

**欧姆电压降**：

$$
\Delta\phi_{ohm} = I_{app} \left( \frac{t_n}{3\sigma_n^{eff}} + \frac{t_p}{3\sigma_p^{eff}} + \frac{t_n}{3\kappa_n^{eff}} + \frac{t_{sp}}{\kappa_{sp}^{eff}} + \frac{t_p}{3\kappa_p^{eff}} \right)
$$

---

## 6. 单元级 SPMe 函数

### 6.1 SPMe_element 函数

用于逐单元电化学计算，支持多 SPMe 并行架构。

**函数签名**：
```julia
function SPMe_element(case::Case, yt_e::Array{Float64}, t::Float64, e::Int; 
    I_e::Float64, T_e::Float64, jacobi::String="update")
```

**输入参数**：

| 参数 | 类型 | 说明 |
|------|------|------|
| `case` | Case | 电池案例对象 |
| `yt_e` | Vector | 单元电化学状态向量 |
| `t` | Float64 | 当前时间 |
| `e` | Int | 单元索引 |
| `I_e` | Float64 | 单元电流（无量纲） |
| `T_e` | Float64 | 单元温度（无量纲） |
| `jacobi` | String | 雅可比更新策略 ("constant" 或 "update") |

**状态向量结构**（SPMe 模型）：

```julia
# yt_e 包含该单元的所有电化学自由度
# 结构：
yt_e = [
    c_s_n[1:Nrn],      # 负极颗粒浓度 (Nrn 个点)
    c_s_p[1:Nrp],      # 正极颗粒浓度 (Nrp 个点)
    c_e[1:Nel]         # 电解液浓度 (Nel 个点)
]
```

**注意**：SPMe 模型中，电位 \(\phi_s\) 和 \(\phi_e\) 通过代数关系求解，不包含在状态向量中。

**输出**：

| 返回值 | 说明 |
|--------|------|
| `M_e` | 单元质量矩阵 |
| `K_e` | 单元刚度矩阵 |
| `F_e` | 单元载荷向量 |
| `variables_e` | 单元变量字典 |

### 6.2 SPMe_variables 函数

**功能**：从状态向量提取电化学变量

**函数签名**：
```julia
SPMe_variables(case, yt, t; I_app, T_e)
```

**支持覆写参数**：
- `I_app`：覆写电流（默认从 opt.Current 计算）
- `T_e`：覆写温度（默认为参考温度）

**输出变量**：

| 键 | 说明 |
|----|------|
| `negative particle surface lithium concentration` | 负极表面浓度 |
| `positive particle surface lithium concentration` | 正极表面浓度 |
| `negative electrode overpotential` | 负极过电位 |
| `positive electrode overpotential` | 正极过电位 |
| `cell voltage` | 电池电压 |

### 6.3 温度依赖性

通过 Arrhenius 方程考虑温度对动力学参数的影响：

$$
k(T) = k_{ref} \exp\left(-\frac{E_a}{R}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right)
$$

$$
D_s(T) = D_{s,ref} \exp\left(-\frac{E_a}{R}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right)
$$

**代码实现**：
```julia
# 在 SPMe_element 中
Eac_D = param.NE.Eac_D  # 无量纲活化能
Eac_k = param.NE.Eac_k

# 温度修正因子
arrhenius_D = exp(-Eac_D * (1.0/T_e - 1.0))
arrhenius_k = exp(-Eac_k * (1.0/T_e - 1.0))
```

---

## 7. 力学耦合

### 7.1 应力耦合扩散系数

当启用力学模型 (`opt.mechanicalmodel == "full"`) 时，颗粒扩散方程中引入应力耦合项：

**耦合扩散系数公式**：

$$
\theta_M = \frac{2 E \Omega^2}{9 (1 - \nu) T}
$$

其中：
- \(E\)：杨氏模量
- \(\Omega\)：偏摩尔体积
- \(\nu\)：泊松比
- \(T\)：温度

**代码实现**：
```julia
# 在 SPMe 函数中
if case.opt.mechanicalmodel == "full"
    variables = Mechanicaloutput(case, variables)
    theta_Mn = variables["negative particle stress coupling diffusion coefficient"][1]
    theta_Mp = variables["positive particle stress coupling diffusion coefficient"][1]
else
    theta_Mn = 0.0
    theta_Mp = 0.0
end

# 传递给扩散矩阵计算
M_np, K_np = ElectrodeDiffusion(param.NE, mesh_np, mesh_np.nlen, csn_gs, theta_Mn)
M_pp, K_pp = ElectrodeDiffusion(param.PE, mesh_pp, mesh_pp.nlen, csp_gs, theta_Mp)
```

### 7.2 对扩散方程的影响

应力耦合修改有效扩散系数：

$$
D_{eff} = D_s \cdot (1 + \theta_M \cdot \sigma)
$$

其中 \(\sigma\) 为应力状态。

---

## 8. 代码位置

| 功能 | 文件 | 函数 |
|------|------|------|
| SPMe 主求解 | src/SPMe.jl | `SPMe` |
| 单元级 SPMe | src/SPMe.jl | `SPMe_element` |
| 变量提取 | src/SPMe.jl | `SPMe_variables` |
| 边界条件 | src/SPMe.jl | `SPMe_BC` |
| 颗粒扩散 | src/ElectrodeDiffusion.jl | `ElectrodeDiffusion` |
| 电解液扩散 | src/ElectrolyteDiffusion.jl | `ElectrolyteDiffusion` |
| 电位求解 | src/ElectrodePotential.jl | `ElectrodePotential` |
| 电解液电位 | src/ElectrolytePotential.jl | `ElectrolytePotential` |
| 力学输出 | src/Mechanical.jl | `Mechanicaloutput` |
