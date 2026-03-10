import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# 参数
r_i = 1.0          # 内径, m
r_o = 2.0          # 外径, m
q = 1000.0         # 内热源, W/m^3
k_r = 1.71         # 径向导热系数, W/(m·K)
h = 10.0           # 对流换热系数, W/(m^2·K)
T_f = 20.0         # 环境温度, °C

# 温度函数（轴对称）
def T(r):
    return (T_f + q/(4*k_r)*(r_o**2 - r**2) 
            + q/(2*k_r)*r_i**2 * np.log(r/r_o) 
            + q/(2*h)*(r_o - r_i**2/r_o))

# 创建网格（极坐标）
theta = np.linspace(0, 2*np.pi, 200)
r = np.linspace(r_i, r_o, 200)
R, Theta = np.meshgrid(r, theta)
X = R * np.cos(Theta)
Y = R * np.sin(Theta)

# 计算温度
T_vals = T(R)

# 绘图
fig, ax = plt.subplots(figsize=(8, 6))
x_flat = X.ravel()
y_flat = Y.ravel()
t_flat = T_vals.ravel()
tri = mtri.Triangulation(x_flat, y_flat)
pc = ax.tripcolor(tri, t_flat, shading='gouraud', cmap='hot')
fig.colorbar(pc, ax=ax, label='Temperature (°C)')
ax.set_aspect('equal')
ax.set_title('Temperature distribution in anisotropic hollow cylinder\n'
             '(radial $k_r=1.71$, circumferential $k_\\theta=34.91$, axisymmetric)')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
# 标记内圆边界
circle_inner = plt.Circle((0,0), r_i, color='white', fill=False, linewidth=2)
circle_outer = plt.Circle((0,0), r_o, color='white', fill=False, linewidth=2)
ax.add_patch(circle_inner)
ax.add_patch(circle_outer)
ax.grid(False)
plt.show()