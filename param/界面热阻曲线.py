import numpy as np
import matplotlib.pyplot as plt

# 物理常数（单位：m, W/mK, W/m^2K）
h_c0 = 1e7
k_air = 0.026
lambda_m = 70e-9          # 平均自由程
beta = 1.0
two_beta_lambda = 2 * beta * lambda_m   # 140 nm
threshold_nm = 70e-9

# 定义 D(delta) 函数
def D(delta, delta0, delta_c):
    if delta <= delta0:
        return 0.0
    elif delta < delta_c:
        return (delta - delta0) / (delta_c - delta0)
    else:
        return 1.0

# 热阻计算函数（根据情况）
def R_case1(delta, delta0, delta_c):
    d_val = D(delta, delta0, delta_c)
    denom = h_c0 * (1 - d_val) + k_air / (delta + two_beta_lambda)
    return 1.0 / denom

def R_case2(delta, delta0, delta_c):
    if delta < delta0:
        denom = h_c0 + k_air / (delta + two_beta_lambda)
    elif delta < threshold_nm:
        d_val = D(delta, delta0, delta_c)
        denom = h_c0 * (1 - d_val) + k_air / (delta + two_beta_lambda)
    else:  # delta >= threshold_nm and delta < delta_c
        d_val = D(delta, delta0, delta_c)
        denom = h_c0 * (1 - d_val) + k_air / (delta + delta0)
    return 1.0 / denom

def R_case3(delta, delta0, delta_c):
    if delta < threshold_nm:
        denom = h_c0 + k_air / (delta + two_beta_lambda)
    elif delta < delta0:
        denom = h_c0 + k_air / (delta + delta0)
    else:  # delta >= delta0 and delta < delta_c
        d_val = D(delta, delta0, delta_c)
        denom = h_c0 * (1 - d_val) + k_air / (delta + delta0)
    return 1.0 / denom

# 三种情况的参数（单位：m）
params = [
    {'delta0': 20e-9, 'delta_c': 50e-9},   # 情况1
    {'delta0': 20e-9, 'delta_c': 100e-9},  # 情况2
    {'delta0': 80e-9, 'delta_c': 120e-9}   # 情况3
]

# 生成横坐标（nm -> m）
delta_nm = np.linspace(0, 200, 1000)
delta_m = delta_nm * 1e-9

# 计算三条曲线（只计算到各自的 delta_c）
R1 = []
R2 = []
R3 = []
for d in delta_m:
    # 情况1：只到 delta_c=50nm
    if d <= params[0]['delta_c']:
        R1.append(R_case1(d, params[0]['delta0'], params[0]['delta_c']))
    else:
        R1.append(np.nan)
    # 情况2：到 delta_c=100nm
    if d <= params[1]['delta_c']:
        R2.append(R_case2(d, params[1]['delta0'], params[1]['delta_c']))
    else:
        R2.append(np.nan)
    # 情况3：到 delta_c=120nm
    if d <= params[2]['delta_c']:
        R3.append(R_case3(d, params[2]['delta0'], params[2]['delta_c']))
    else:
        R3.append(np.nan)

# 绘图
plt.figure(figsize=(10,6))
plt.plot(delta_nm, R1, label=r'Case 1: $\delta_0=20nm,\delta_c=50nm$', lw=2)
plt.plot(delta_nm, R2, label=r'Case 2: $\delta_0=20nm,\delta_c=100nm$', lw=2)
plt.plot(delta_nm, R3, label=r'Case 3: $\delta_0=80nm,\delta_c=120nm$', lw=2)

plt.xlabel(r'Gap thickness $\delta$ (nm)')
plt.ylabel(r'Thermal resistance $R_{th}$ (m$^2$K/W)')
plt.title('Interface thermal resistance vs gap thickness')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.xlim(0, 200)
plt.tight_layout()
plt.savefig('output/界面热阻曲线.png', dpi=300)
plt.show()