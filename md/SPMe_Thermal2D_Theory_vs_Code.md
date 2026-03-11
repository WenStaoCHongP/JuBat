# SPMe + 2D 热分布（Jellyroll）理论与实现（按有限元流程）

目标：Jellyroll 场景下，采用“每个热单元 = 一个子 SPMe 小电池”的耦合思路；其它模型与一维热不受影响，涉及代码Jellyroll.jl。

---

## 1) 几何与网格（Geometry & Mesh）

- Jellyroll 几何：内半径 Rin、外半径 Rout、高度 H；层序为 PCC → PE → SP → NE → NCC。螺旋参数 `jellyroll_spiral_params`。
- 热网格（2D，Q4）：`jellyroll_Q4_mesh(param_dim; nx, ny, gsorder, crop_to_annulus, crop_mode)` 
   `:collector_seeded`：基于两条集流体边界的螺旋曲线在角向等距采样，生成“以直代曲”的条带 Q4 单元（见下）。
- 电化学网格（1D，L2）：NE/SP/PE 三段独立于热网格；通过“单元平均+几何映射”与热网格耦合（见第6节）。

热网格划分：沿集流体边线布种 + 以直代曲条带网格（每单元包含完整层序）

- 设计目标：以集流体（NCC、PCC）边线为“导轨”布点，单元横向跨越一次完整层序（NCC→NE→SP→PE→PCC），与材料分层逻辑一致，便于并联网络与热源聚合。
- 核心思想（以直代曲）：直接在两条集流体边界的螺旋曲线 r=a+bθ 与 r=a+bθ+T_rep（T_rep=t_repeat）上等距取点，连接成两条带状“轨道”；在相邻轨道之间用直线连边形成变形的条带 Q4 单元，使每个单元横向穿过完整层序。
- 生成步骤（`jellyroll_collector_seed_mesh(...)`）：
  1) 从 `jellyroll_spiral_params` 读取 (a,b,Rin,Rout,t_repeat, boundaries)。设 s_in=0、s_out=t_repeat。
  2) 在 θ∈[θ0,θ1] 上等距采样两条曲线 r_in(θ)=a+bθ+s_in，r_out(θ)=a+bθ+s_out，其中 θ0=max(0,(Rin−a−s_in)/b)，θ1=min((Rout−a−s_out)/b,(Rout−a)/b)。
  3) 以相邻采样点配对生成 Q4 单元：[(in_i, out_i, out_{i+1}, in_{i+1})]。这是“以直代曲”的近似：边用直线连接曲线上的采样点。
  4) 为每个条带定义局部切/径向基：θ̄ = angle((x_in_i+x_out_i)/2,(y_in_i+y_out_i)/2)，e_r=[cosθ̄,sinθ̄]，e_θ=[−sinθ̄,cosθ̄]，供各向异性 K 旋转。
  5) 条带天然覆盖完整层序，赋元素层权重 f={f_NE,f_SP,f_PE,f_PCC,f_NCC}。默认 f_k = thickness_k / t_repeat；也可按局部几何做微调（例如两轨道间真实曲边长度修正）。

检查项目：果冻卷的Q4单元由SetMesh.jl 中Q4等参元映射得到，以减少代码量。

数学化定义（条带参数化）：

- 轨道采样：给定 n_θ 段，θ_j = θ0 + j(θ1−θ0)/n_θ,  j=0..n_θ。
- 节点：X_in(j) = (a+bθ_j+s_in)[cosθ_j, sinθ_j]，X_out(j) = (a+bθ_j+s_out)[cosθ_j, sinθ_j]。
- 元素：E_i = (X_in(i), X_out(i), X_out(i+1), X_in(i+1))，i=0..n_θ−1。
- 面积与法切向：可直接由 Q4 几何计算；θ̄_i ≈ (θ_i + θ_{i+1})/2。
- 等效性质聚合（条带/层权重一致）：
  - K_eff：法向（穿层、串联）调和平均，切向（沿层、并联）算术平均；再按 θ̄ 旋转至全局坐标：
    - λ_r,eff = 1 / (∑_k f_k / λ_{r,k})；λ_θ,eff = ∑_k f_k λ_{θ,k}
    - K_i = λ_r,eff e_r e_r^T + λ_θ,eff e_θ e_θ^T
  - ρ、c：体积分数加权 ρ_eff=∑_k f_k ρ_k，(ρc)_eff=∑_k f_k ρ_k c_k。
  - 内热源：先在电化学侧求得各层的体积内热源密度 Q_k（PCC、PE、SP、NE、NCC）。在映射到热单元时采用单元面积积分：

    q_e = (1/A_e) * \sum_k ∫_{A_e ∩ layer_k} Q_k(x,y) dA

    当层内 Q_k 在单元内近似常数时，上式退化为 q_e = \sum_k (A_{e,k}/A_e)·Q_k，其中 A_{e,k} 为层 k 与单元 e 的投影重叠面积。如果网格预先给出 `layer_weights f_k(e)`（即面积分数 f_k(e)=A_{e,k}/A_e），则可近似为 q_e ≈ \sum_k f_k(e)·Q_k。


## 2) 材料与参数（Materials & Parameters）

- 热材料：各层等效 ρ、c、k；可选各向异性：`K(θ)=λ_r e_r e_r^T + λ_t e_θ e_θ^T`（`jellyroll_effective_K_at`）。
- 电化学参数：cs_max、ce0、σ、κ、k、rs、as、厚度等，来源于参数文件（LGM50、Northrop、Jellyroll）。
- LGM50 vs Northrop：归一化公式一致，数值不同；细见第5节与 `src/parameters/*.jl`、`src/SetParams.jl`。

---

## 3) 控制方程（Governing Equations）

- 电化学（SPMe，一维厚度向 x∈[0,L]；负极区 Ω_n=[0,L_n]，隔膜 Ω_s=[L_n,L_n+L_s]，正极区 Ω_p=[L_n+L_s,L]）
  - 变量与常量：
    - 固相锂浓度 c_s(r,x,t)，代表性球形颗粒半径 r∈[0,R_s(x)]；电解液浓度 c_e(x,t)；
      固相/液相电位 φ_s(x,t)、φ_e(x,t)；界面反应通量 j(x,t)（A/m²-界面）；比表面积 a_s(x)；
      有效输运系数 Ds, De, σ_eff, κ_eff（Bruggeman 等效）；跨膜阳离子迁移数 t_+^0；气体常数 R，法拉第常数 F。
  - 颗粒内扩散（各电极区内独立，隔膜无固相）：
    - ∂c_s/∂t = (1/r²) ∂/∂r ( Ds · r² · ∂c_s/∂r ),  r∈(0,R_s)
    - 边界：∂c_s/∂r|_{r=0} = 0； −Ds ∂c_s/∂r|_{r=R_s} = j/F
    - 表面浓度 c_s^surf(x,t) = c_s(r=R_s,x,t)
  - 电解液物质守恒（负极/隔膜/正极三段）：
    - ∂(ε_e c_e)/∂t = ∂/∂x( D_e^eff ∂c_e/∂x ) + (1−t_+^0) a_s j / F，  x∈Ω_n∪Ω_s∪Ω_p
    - ε_e 为孔隙率；Ω_s 内 a_s=0 且 j=0（无反应），源项仅在 Ω_n、Ω_p。
  - 电荷守恒（准稳态电流连续）：
    - 固相： −∂/∂x( σ^eff ∂φ_s/∂x ) = a_s j，x∈Ω_n∪Ω_p；Ω_s 无固相（σ^eff=0）
    - 液相： −∂/∂x( κ^eff ∂φ_e/∂x + 2(1−t_+^0) (RT/F) κ_D^eff ∂ln c_e/∂x ) = −a_s j
      （常用近似：忽略扩散电位项时 −∂/∂x( κ^eff ∂φ_e/∂x ) = −a_s j）
  - 界面动力学（Butler–Volmer，α_a≈α_c≈0.5）：
    - 过电位： η = φ_s − φ_e − U(c_s^surf, T)
    - 交换电流密度： j_0 = k · (c_e)^α · (c_s^surf)^α · (c_s,max − c_s^surf)^α
    - 反应通量： j = 2 j_0 sinh( F η / (2RT) )
  - 边界与接口条件：
    - 固相电位边界（集流体）
      - x=0（负极集流体）： −σ^eff ∂φ_s/∂x = I_app/A_cell（放电为正，指向+x）
      - x=L（正极集流体）：  σ^eff ∂φ_s/∂x = I_app/A_cell
    - 液相电位与浓度：各分区内连续；端部常取绝通量 ∂c_e/∂x=0，φ_e 相对参考定零（例如端点设 φ_e=0 或平均为零）。
  - 端电压： V_cell = (φ_s − φ_e)|_{x=L} − (φ_s − φ_e)|_{x=0}

- 热（二维横截面 Ω_T⊂R²，Q4 单元）
  - 能量方程： (ρ c) ∂T/∂t = ∇·( K ∇T ) + q_total(x,y,t)
  - 热源分解（条带层权重聚合 → 2D 元素）：
    - 反应热： q_rxn = a_s j η
    - 可逆热（熵变项）： q_rev = a_s j T (∂U/∂T)
    - 欧姆热：
      - 固相： q_ohm,s ≈ i_s² / σ^eff，其中 i_s = −σ^eff ∂φ_s/∂x（1D 近似用于层平均）
      - 液相： q_ohm,e ≈ i_e² / κ^eff，其中 i_e = −κ^eff ∂φ_e/∂x（忽略扩散电位项时）
    - 各层合成：
      - 负极： Q_NE = a_s,n j_n η_n + a_s,n j_n T (∂U_n/∂T) + q_ohm,s,ne + q_ohm,e,ne
      - 隔膜： Q_SP = q_ohm,e,sp
      - 正极： Q_PE = a_s,p j_p η_p + a_s,p j_p T (∂U_p/∂T) + q_ohm,s,pe + q_ohm,e,pe
      - 集流体：
        - 正集流体（PCC）：Q_PCC ≈ I_app² / (3 · σ_PCC,eff)
        - 负集流体（NCC）：Q_NCC ≈ I_app² / (3 · σ_NCC,eff)
  - 2D 映射（collector-seeded 首选）：热源映射采用基于单元面积的积分方法：先在电化学/子单元侧求得各层（PCC/PE/SP/NE/NCC）的体积内热源密度 Q_k(x,y)（通常为分层常数或按电化学点值计算得到的层平均），然后在每个热单元 e 上对各层在该单元投影的面域进行面积积分，得到该单元来自各层的面功率，再除以单元面积得到单元平均的体积内热源。形式上，令 A_e 为单元面积，A_{e,k} := area( A_e ∩ projection(layer k) )，则单元平均体积热源为

    q_e = (1 / A_e) * \sum_k Q_k * A_{e,k}

  其中 Q_k 为层 k 的体积热源（W/m^3，或无量纲 q*），A_{e,k} 为层 k 在该单元的重叠面积（m^2）。若 Q_k 在元素内不是常数，则用数值积分替代 Q_k * A_{e,k}，即 A_{e,k}·Q_k → ∫_{A_e∩layer_k} Q_k(x,y) dA。

  得到 q_e 后，根据热模块的表达习惯，将其转为装配所需的量：若热方程以体积热源 q (W/m^3) 装配，则直接使用 q_e；若装配需要面热源（每面积功率 W/m^2），则使用 q_areal_e = q_e * H（H 为模型选定的单元高度/厚度尺度，若用无量纲则按 scale 变换）。
  - 各向异性导热张量装配（实现已更新）：
    - 首选：读取元素层权重 f_k，先在局部径/切向上聚合得到 λ_r,eff 与 λ_θ,eff（分别为调和平均与算术平均；均按 k_ref 无量纲化），再按元素极角 θ̄ 旋转得到 Kxx/Kxy/Kyy 并装配。
    - 回退：若无 f_k，则使用全局螺旋有效 λ_r_eff/λ_t_eff（同样先除以 k_ref）并按 θ̄ 旋转。
    - 代码对应：`src/ThermalDistributed.jl/ThermalDistributed2D` 的各向异性分支。
  - 边界条件（外边界 Γ_out，对流；内孔 Γ_in，默认绝热）：
    - −n·(K ∇T) = h (T − T_amb) on Γ_out；  n·(K ∇T) = 0 on Γ_in
  - 初值： T(x,y,0) = T_0（通常取 T_amb）。

---

### 3.x 控制方程的离散格式（半离散/弱式与矩阵形式）

本节给出与上述控制方程相对应的有限元/有限差分离散格式，统一写成半离散矩阵方程，时间推进通常采用后退欧拉（或 Crank–Nicolson）。

- 固相颗粒内扩散（球坐标，Fick）：
  - 方程：∂c_s/∂t = (1/r²) ∂/∂r ( D_s r² ∂c_s/∂r ),  r∈(0,R_s)
  - 边界：∂c_s/∂r|_{r=0}=0； −D_s ∂c_s/∂r|_{r=R_s} = j/F
  - 有限差分（典型二阶、径向网格 r_i=iΔr，i=0..N）：
    - 体节点 i=1..N−1：
      (c_i^{n+1}−c_i^n)/Δt = (1/r_i²) [ (D_{i+1/2} r_{i+1/2}² (c_{i+1}^{n+θ}−c_i^{n+θ})/Δr − D_{i−1/2} r_{i−1/2}² (c_i^{n+θ}−c_{i−1}^{n+θ})/Δr ) / Δr ]
      其中 θ∈[0,1]，θ=1 为后退欧拉，θ=1/2 为 Crank–Nicolson。
    - r=0 采用偶对称消去虚点（c_{−1}=c_1）或标准中心差分的对称处理；
    - r=R_s 用通量边界：−D (c_N^{n+θ}−c_{N−1}^{n+θ})/Δr = j^{n+θ}/F。
    - 写成矩阵：M_s c_s^{n+1} = M_s c_s^n + Δt (K_s c_s^{n+θ} + f_s^{n+θ})，后退欧拉时为 [(M_s/Δt) − K_s] c_s^{n+1} = (M_s/Δt) c_s^n + f_s^{n+1}。

- 电解液物质守恒（1D，L2 有限元）：
  - 方程：∂(ε_e c_e)/∂t = ∂/∂x( D_e^eff ∂c_e/∂x ) + (1−t_+^0) a_s j / F
  - 弱式：∫ Nᵀ ε_e 
        \dot c_e dΩ − ∫ (∂N/∂x)ᵀ D_e^eff (∂c_e/∂x) dΩ = ∫ Nᵀ (1−t_+^0) a_s j / F dΩ
  - 半离散：M_ce \dot c + K_ce c = f_ce，时间离散后 [(M_ce/Δt)+K_ce] c^{n+1} = (M_ce/Δt) c^n + f_ce^{n+1}。

- 固相电荷守恒（1D）：
  - 方程（负/正极区）：−∂/∂x( σ^eff ∂φ_s/∂x ) = a_s j
  - 弱式：∫ (∂N/∂x)ᵀ σ^eff (∂φ_s/∂x) dΩ = ∫ Nᵀ a_s j dΩ + 边界通量项
  - 半离散（稳态/准稳态）：K_s φ_s = f_s，其中 K_s = ∫ Bᵀ σ^eff B dΩ，f_s = ∫ Nᵀ a_s j dΩ + f_s^{(Neumann)}；
    集流体边界的电流通量 −σ^eff ∂φ_s/∂x = ± I_app/A_cell 以自然边界线性形式加入 f_s。

- 电解液电荷守恒（1D）：
  - 全式：−∂/∂x( κ^eff ∂φ_e/∂x + 2(1−t_+^0) (RT/F) κ_D^eff ∂ln c_e/∂x ) = −a_s j
  - 常用简化（忽略迁移项）：−∂/∂x( κ^eff ∂φ_e/∂x ) = −a_s j
  - 弱式（简化式）：∫ (∂N/∂x)ᵀ κ^eff (∂φ_e/∂x) dΩ = ∫ Nᵀ (a_s j) dΩ
  - 半离散：K_e φ_e = f_e，必要时把迁移项写成对 c_e 的耦合项 C(c_e)。

- Butler–Volmer 界面动力学与耦合：
  - η = φ_s − φ_e − U(c_s^surf,T)， j = 2 j_0 sinh(Fη/(2RT))， j_0 = k·c_e^α·(c_s^surf)^α·(c_s,max−c_s^surf)^α
  - 数值实现：j 通常在高斯点或节点评价并回填到右端项（f_ce、f_s）与热源（Q_rxn、Q_rev）。非线性可用外循环/内牛顿或显式滞后处理。

- 热传导（2D，Q4）：
  - 方程： (ρc) ∂T/∂t = ∇·(K ∇T) + q
  - 弱式/装配（Scheme B 无量纲，见第 8 节）：
    MT = ∬ (ρc/ρc_ref) Nᵀ N dΩ*， KT = ∬ (k/k_ref) Bᵀ B dΩ*， FT = ∬ q* N dΩ*
  - 半离散与时间离散：[(MT/Δt_th)+KT] T^{n+1} = (MT/Δt_th) T^n + FT^{n+1} + F_boundary^{n+1}。

### 3.y 并联系统求解的方程组形式（Parallel Network Equations）

下面给出更严格的分电流求解数学表述与两种推荐的数值解法（基于仓库附带笔记“求分电流办法.md”）：

设每个并联支路的分电流为 x_i（i=1..n），公共端电压为 V_cell（对所有支路相同），总电流为  I_total。每支路满足非线性端口关系

  f_i(x_i) = C_{1,i} + C_{2,i}[ asinh(C_3 x_i) + asinh(C_4 x_i) ] - C_{5,i} x_i = V_cell

其中 C_{1,i},C_{2,i},C_{5,i} 可依 i 不同（与局部 SOC、温度及材料相关），C_3,C_4 为全局常数。约束为

  Σ_{i=1}^n x_i A_i= I_total.

这构成一个 n+1 元未知量的非线性方程组（x_1..x_n, V_cell）。因为 f_i 通常单调（物理上端口电压随支路电流单调），推荐的数值解法有两种：

方法一（基于 V 的迭代法，推荐当 f_i 关于 x 单调且容易求逆时）
  1. 对于给定的候选 V，逐支路解方程 f_i(x_i)=V（数值求根）得到 x_i = g_i(V)。
  2. 构造 S(V)=Σ_i g_i(V)，解 S(V)=R（标量方程）得到 V。可用二分法或稳健的牛顿-切换策略。初始 V 可取 V0 = (1/n) Σ_i f_i(R/n)。
  3. 当 V 求出后，x_i = g_i(V)。该方法的优点是将高维问题降为标量，稳健且易于并行（每个 g_i 可并行求解）。

方法二（牛顿法求解耦合系统）
  1. 把未知向量写作 (x_1..x_n, V)ᵀ，定义残差向量

       F = [ f_1(x_1)-V,
             f_2(x_2)-V,
             ...,
             f_n(x_n)-V,
             Σ_i x_i - R ]ᵀ.

  2. 雅可比矩阵 J 的结构为带有对角块的稀疏矩阵：

       J = \begin{pmatrix}
         ∂f_1/∂x_1 & 0 & … & 0 & -1 \\
         0 & ∂f_2/∂x_2 & … & 0 & -1 \\
         … & … & … & … & … \\
         0 & 0 & … & ∂f_n/∂x_n & -1 \\
         1 & 1 & … & 1 & 0
       \end{pmatrix}

     其中

       ∂f_i/∂x_i = C_{2,i} ( C_3/√(C_3² x_i² + 1) + C_4/√(C_4² x_i² + 1) ) − C_{5,i}.

  3. 采用牛顿迭代：解线性子问题 J Δ = −F，更新 (x,V)←(x,V)+Δ，直至 ‖F‖ 小于容差。初值可取 x_i=R/A_total, V=(1/A_total)Σ f_i(R/A_total)。该方法收敛快（若初值足够好），但需要解 n+1 维线性系统，并保证 Jacobian 可逆。

实用建议：
  - 若 f_i 单调且每支路求逆容易（可用一维求根器），优先使用方法一（V-迭代），在并联支路数目很大时更高效；若 f_i 形式已知且对导数易评估，可用方法二获得二阶收敛速度。
  - 为稳健性，方法一中 S(V) 求根可先用夹逼/二分确定区间，再用牛顿加速；方法二中可对 Δ 做阻尼或使用线性求解回退策略。若出现不收敛，可回退到导纳加权近似作为初值。

实现注意（变量语义与守恒约束）:
  - 约定：在并联分配的数学表述与代码实现中，向量 x = [x_1..x_n] 明确表示每个热单元的物理电流 I_e（单位为 A 或无量纲），因此全局守恒条件应写成
  - 接口约定（代码实现）：求解器应在成功后写回：
    - `variables["thermal2D element current"] = I_e`  # 各单元电流，保证 Σ_e A_e·I_e = I_total
    - `variables["thermal2D common voltage"] = V`         # 公共端电压相等

## 4) 无量纲化与标尺（Non-dimensionalization & Scales）

- 电化学标尺（`SetParams.ChooseCell/NormaliseParam`）：
  - I_typ = I1C；t0 = 3600 s；j_ref = I_typ/(a0 L A_cell)；
  - ts_p = F cs_max,p A_cell L / I_typ；ts_n、te 类似；Ds_p = r0²/ts_p；De = L²/te；
  - ϕ_ref = R T_ref / F；σ_ref = κ_ref = L I_typ/(ϕ_ref A_cell)；
  - k_ref：k_p = j_ref/(cs_max,p √ce0)，k_n 类似；
  - U_nd = U/ϕ_ref；(dU/dT)_nd = (dU/dT)·(T_ref/ϕ_ref)。
- 热标尺（Scheme B）：
  - L_th、k_ref、(ρc)_ref、q_ref = k_ref T_ref/L_th²、t_th = (ρc)_ref L_th²/k_ref、Bi = h L_th/k_ref。
  - T* = T/T_ref，t* = t/t_th，q* = q/q_ref。

---

## 5) 边界与初值（BCs & ICs）

- 初值：T(x,0)=T_amb 或参数指定的初温。
- 边界/表面条件（本项目约定）：
  1) 上/下表面（厚度方向两面）存在对流，采用“域内等效反应项”处理；
  2) 侧面为对流边界，包括外圈 Γ_out（r=Rout）与内孔 Γ_in（r=Rin）；
  3) 极耳区域边界（见 5.x）。

写法与装配：
- 侧面对流（Robin，边界积分项）：
  - 无量纲（Scheme B）：在 Γ_side=Γ_out∪Γ_in 上，−(k/k_ref) ∂T*/∂n* = Bi_side (T* − T_amb*)，
    装配 KT += ∫_{Γ_side} Bi_side NᵀN dL*，FT += ∫_{Γ_side} Bi_side T_amb* N dL*。
  - SI：−k ∂T/∂n = h_side (T − T_amb)，装配 KT += ∫_{Γ_side} h_side NᵀN dL，FT += ∫_{Γ_side} h_side T_amb N dL。
- 上/下表面对流（域内等效“反应项”）：
  - 3D 到 2D 推导（假设厚度方向温度近似均匀）：将顶/底面对流通量 q''_tb = h_top (T−T_amb_top) + h_bot (T−T_amb_bot) 等效为 2D 方程中的体源项并移至左端，得到域内“吸收”项。
  - SI 组装（2D 每单位厚度方程）：
    KT_tb += ∬_Ω ((h_top+h_bot)/H) NᵀN dΩ，FT_tb += ∬_Ω ((h_top T_amb_top + h_bot T_amb_bot)/H) N dΩ。
  - 无量纲（Scheme B）：令 Bi_top = h_top L_th/k_ref，Bi_bot = h_bot L_th/k_ref，
    则 KT_tb += ∬_Ω* ((Bi_top+Bi_bot)/H) NᵀN dΩ*，FT_tb += ∬_Ω* ((Bi_top T_amb_top* + Bi_bot T_amb_bot*)/H) N dΩ*。
  - 注：若采用“每单位厚度”的 2D 模型，需显式除以物理厚度 H；若转入 3D 或改变缩放，需相应调整。

图例：热边界类型与装配对应关系如下图所示（外圈对流、内圈绝热，并给出恒温与辐射的标注方式，便于与实现中的 Bi、边界项装配对应起来）：

![热边界标注示意图例](../figure_thermal_bc_legend.svg)

说明：
- 侧面对流（边界积分）：`q = h_side (T − T_amb)`；无量纲为 `Bi_side = h_side L_th / k_ref`。
- 上/下表面对流（域内项）：在 2D 模型中以“体吸收”形式加入，系数为 `(h_top+h_bot)/H`（SI）或 `(Bi_top+Bi_bot)/H`（Scheme B）。
- 绝热：`q = 0 ⇔ ∂T/∂n = 0`，无边界通量项。
- 恒温（Dirichlet）：直接在边上施加 `T = T_spec` 或通过强/弱方式实现。
- 辐射（可选）：`q = εσ (T^4 − T_amb^4)`，可线性化为等效 `h_rad` 并与对流叠加到 Bi 中。

### 5.x 极耳（Tab）热边界处理

极耳区域通常位于集流体外缘特定弧段或外接矩形窗口，热边界可采用以下几种常见物理模型：

- 选区标记：
  - 几何上把极耳边界 Γ_tab 作为边界 Γ_out 的子集，或在网格生成时标注对应边段 ID（如 `boundary_id=:tab`）。
  - 极耳位置在果冻卷螺旋的起点和螺旋的终点。



- 以“极耳粘接区域温度等于极耳表面温度”处理极耳：
  - 思路：等效为受影响的条带单元节点温度为极耳表面温度。
  - 选择器：定义受极耳影响的元素集合 S_tab，可由边界标注获得。
  - 极耳发热功率由电阻发热计算得到，先简化实现:极耳温度为随时间线性增长。
- 恒温接触（Dirichlet）：
  - 直接在 Γ_tab 施加 T = T_tab（或 T* = T_tab/T_ref）。
  - 数值上采用强施加（修改方程系数）或弱惩罚法（KT += α_p ∫ NᵀN；FT += α_p ∫ T_tab N），α_p 取较大值确保约束。



实现建议（代码对应）：
- 在 `ThermalDistributed2D_BC` 中增加对 tab 边界的识别与装配：
  - 输入：`bc.tab = (type=:robin, h=h_tab, Tamb=T_amb_tab, select=selector)`；selector 可为边界 id 集合或坐标过滤器。
  - Scheme B：直接生成 Bi_tab 和 T_amb_tab* 进入边界积分；SI 模式则使用 h_tab、T_amb_tab 与物理长度 dL 进行装配。
  - 若选择 Dirichlet，则走统一的强施加路径（组装稀疏投影矩阵或行替换）。

---

## 6) 耦合策略与并联系统（Coupling & Parallel Network）

- “每热单元=子 SPMe”：
  - 子单元面积 A_e；A_global 为有效受流面积。
  - I_typ,e = I1C_total；I_app_e = I_e / I_typ,e；T_nd,e = T_e/T_ref。
  - 时间尺度保持不变：ts_e = F cs_max A_e L / I_typ,e = ts_global。
  - 步骤：分配 I_e → 逐单元求 I_app_e 与 T_nd,e → 子 SPMe → 得 q_nd,e → 热步。

### 最终目标耦合架构（与图一致）

- 数学“契约”（输入/未知/约束）：
  - 输入：I_total(t)、T_amb(t)、上一步状态（SPMe 状态、T 场）。
  - 未知：每单元电流 I_e、每单元子 SPMe 状态 y_e、元素热源 q_e、温度场 T(x,y)。
  - 约束：
    1) 端电压一致与电流守恒：Σ_e I_e*A_e = I_total，I_e = I_e(V; T_e, θ_e, SOC_e,…)，求 V。
    2) 电化学：对每单元，SPMe 方程与 BC 满足，驱动量为 I_app 与 T_nd。
    3) 热：2D 能量方程满足，源项为各元素 q_e（含 PCC/NCC）。


### 优化目标（Coupling & Parallel）

- 在每个热时间步内， I_app与 T_nd运行“子 SPMe”，得到 q_nd后再装配 2D 热；保持电化学无量纲标尺与全局一致。
- 并联分流：“以端电压 V 为未知”的单变量非线性分流（I_e=I_e(V,T,θ,SOC)）。
- 网格-材料一致性：采用 `crop_mode=:collector-seeded` 的条带网格，保证每元素包含完整层序，并附带层权重 f_k；等效 K、ρc 与热源 q 统一按 f_k 聚合，减少二次判层误差。
- 热源全量：NE/SP/PE + PCC/NCC 的反应/可逆/欧姆热全部计入；单位与热模块缩放一致（支持无量纲/SI 双路径）。
- 性能与稳定性：提供“批处理/向量化子 SPMe”接口，避免逐元素重复构造；对迭代强耦合设置最大迭代与松弛，保证步进稳定与耗时可控。

根据初始节点温度场和总电流-并联模型-计算内热源（假设电化学热源在热单元内均匀分布，经映射后取单元平均）输入到热模型，由热模型计算得到新的节点温度场，取每个热单元的节点平均作为子SPMe模型新时间步输入的温度以及总电流-并联模型-计算内热源，如此迭代，直到收敛。


### 6.x 并联关系逻辑架构图（Parallel Network Architecture）

输入初值T[节点]、I_total[总电流]、SPMe初始条件相同，I_e[单元电流]=I_total[总电流]/A_total[总面积]
I_app[单元应用电流]=I_e[单元电流]
SPMe：输入T[单元平均]、I_app[单元应用电流]，输出q_e[单元热源]
Thermal：输入q_e[单元热源]，输出T_nd[单元节点温度]
分电流求解：输入T[单元平均]、y[SPMe状态]，输出I_app[单元应用电流]
I_e[单元电流]=I_app[单元应用电流]、T[节点]=T_nd[单元节点温度]

## 7) 内热源与层到网格映射（Heat Sources & Mapping）

 - 层平均（SPM/SPMe）体热源：
  - 反应：Q_rxn = a_s j η；可逆：Q_rev = a_s j T (dU/dT)。
  - 欧姆（一致流近似）：
    - 固相：q_ohm,s ≈ I_app²/(3 σ_eff)；电解液：q_ohm,e ≈ I_app²/(3 κ_eff)（隔膜为 I_app²/κ_sp）。
  - 合层：
    - NE：Q_NE = a_s,n j_n η_n + a_s,n j_n T dU_n/dT + q_ohm,s,ne + q_ohm,e,ne。
    - SP：Q_SP = q_ohm,e,sp。
    - PE：Q_PE = a_s,p j_p η_p + a_s,p j_p T dU_p/dT + q_ohm,s,pe + q_ohm,e,pe。
    - PCC：Q_PCC = q_ohm,PCC ≈ I_app²/(3 σ_PCC,eff)。
    - NCC：Q_NCC = q_ohm,NCC ≈ I_app²/(3 σ_NCC,eff)。
 - 映射到 2D 热网格（基于单元面积积分）：
  - 基本思路：先计算每一层的体积热源密度 Q_k（W/m^3 或无量纲 q*）；再对每个热单元 e 在平面投影上按层区域做面积积分得到该单元来自各层的面功率；最后把面功率除以单元面积得到单元平均体积热源 q_e（W/m^3）。数学上：

    A_{e,k} = area( A_e ∩ projection(layer k) )

    P_{e,k} = ∫_{A_e ∩ layer_k} Q_k(x,y) dA ≈ Q_k · A_{e,k}  (若 Q_k 在层内近似常数)

    q_e = (1 / A_e) * \sum_k P_{e,k}  = \sum_k (A_{e,k}/A_e) · Q_k

  - 数值实现建议：对 A_{e,k} 的计算可采用精确的平面多边形相交（推荐）、或基于高精度规则采样的数值积分（例如高斯采样或规则网格采样）。若已有 `layer_weights f_k(e)`（由几何采样得到的面积分数），则 q_e = \sum_k f_k(e) · Q_k，与上式等价；但显式积分在处理曲边/截断单元时更准确。
  - 若 Q_k 在元素内变化显著（例如子 SPMe 给出非均匀分布），则应对 Q_k(x,y) 用高斯点或采样点进行积分：

    P_{e,k} ≈ \sum_{i∈Gauss} w_i Q_k(x_i,y_i)

  - 从体积热源到面功率的转换：若需要面密度（W/m^2），令 q_areal_e = q_e · H，其中 H 是单元的轴向高度（或模型用的厚度尺度）；在无量纲模式下按 scale.q_th 与 scale.L_th 做相应换算。
  - 装配：把 q_e（或经尺度转换的 q*）赋值到 `variables["heat_source_fields"]` 的每个单元条目，热装配时使用 FT = ∬ q* N dΩ*（或对应的 SI 量）。

---

## 8) 离散与装配（Discretization & Assembly）

- 电化学：1D L2（NE/SP/PE 分段），形成 M、K、F，边界 `SPMe_BC`；求解在 `SPMe` 中。
- 热（2D，Q4）：
  - 质量：MT = ∬ (ρc/ρc_ref) N^T N dΩ*；刚度：KT = ∬ (k/k_ref) B^T B dΩ*；右端：FT = ∬ q* N dΩ*。
  - 各向异性 K：用 `jellyroll_effective_K_at(θ)` 构造 Kxx/Kxy/Kyy 并装配。
  - 边界对流：沿边积分 KT += ∫ Bi N^T N dL*，FT += ∫ Bi T_amb* N dL*。

---

## 9) 时间推进与算法（Time Integration & Algorithm）

- 时间步：电化学以 t0 无量纲 dt；热以 t_th 无量纲 dt_th，关系 dt_th = dt · (t0/t_th)。
- 步骤（每个时间步内进行内迭代）：
  1) 初值：T^(0) ← T_old。
  2) k=0,1,…：
     - 分流：由 T^(k) 评估导纳，解 Σ_e I_e(V,T^(k))=I_total 得 V^(k)、I_e^(k)。
     - 子 SPMe：对每元素 e，I_app,nd,e^(k)、T_nd,e^(k) 作为输入，运行 SPMe 得 q_e^(k)。
     - 2D 热：以 q_e^(k) 求 T^(k+1)。
     - 收敛：||T^(k+1)−T^(k)||_∞<ε_T 且 ||I^(k)−I^(k−1)||_1/I_total<ε_I（k≥1），否则 T 松弛后继续迭代或减小时间步。
  3) 收敛后置 T_new ← T^(k+1)。


输入（Inputs）：
- `I_total(t)`、`T_amb(t)`、时间网格 `opt.time/dt`；
- 几何与参数：`param_dim`（SI 尺度/缩放）、`param`（无量纲电化学）、材料数据（σ、κ、k、ρc、dU/dT、a_s 等）；
- 网格与拓扑：`mesh_th`（Q4）、可选 `layer_weights f_k`；
- 开关：`opt.coupling_mode`、`opt.per_element_spme`、`opt.parallel_solve_V`、`opt.units_thermal`、 容差与松弛等。

输出（Outputs，典型 result/variables 键）：
- `result["thermal2D T_nodes [K]"]`、`result["thermal2D nodes xy [m]"]`；
- `result["time [s]"]`、`result["cell voltage [V]"]`；
- 跟踪：`result["thermal2D tracked element time [s]"]`、`result["thermal2D tracked element T [K]"]`（若启用跟踪与导出）；
- 诊断：`variables["thermal2D element current"]`、`variables["thermal2D common voltage"]`、`variables["thermal2D layer_weights"]`（若开启）；
- 可选：能量平衡日志、`T_mean` 时间序列、并联分流迭代状态等。

调试日志可通过 `opt.debug_coupling=true` 与 `opt.debug_sample_elems` 控制。

## 10) 校核与自检（Verification & Checks）

- 能量守恒：∑域内热源 − ∑边界散热 ≈ dE_th/dt（阈值 ~1%）。
- 快速清单：
  - 面积与电流：A_global、A_e、I_e、I_typ,e、I_app,nd,e。
  - 温度：T_e → T_nd,e。
  - SPMe 输出：η/j 标量化；dUdT 标量化；q_nd 为标量。
  - 热源映射：S_q、q_SI_e、q_areal_e（或 q*）。
  - 装配一致性：MT/KT/FT/BC 与单位/缩放一致。

---

## 11) 代码映射（Where in the Code）

- 网格与几何：`src/Jellyrollmodel.jl`（`jellyroll_Q4_mesh`、`jellyroll_spiral_params`、`material_at`、`jellyroll_element_centers`、`jellyroll_effective_K_at`、`jellyroll_element_layer_weights`）。
- 电化学：`src/SPMe.jl`（`SPMe`、`SPMe_variables`、`SPMe_BC`）。
- 热：`src/ThermalDistributed.jl`（`ThermalDistributed2D`、`ThermalDistributed2D_BC`、`heatQ_Source`、`DistributeCurrentParallel!`）。
- 参数与归一化：`src/SetParams.jl`（`ChooseCell`、`NormaliseParam`）。
- 主循环：`src/Solve.jl`（时间推进与热耦合；当 `opt.collector_seeded=true` 时自动计算并写入 `variables["thermal2D layer_weights"]`）。



## 4. 快速自检清单（Jellyroll+2D 热 + SPMe）

- 电流：
  - [ ] `I_typ_e = I1C_total`；`I_app = I_e/I_typ_e` 输入 SPMe。
- 温度：
  - [ ] `T_e` 用热单元 nodal 均温；`T_e = T_e/T_ref` 输入 SPMe。
- SPMe 输出：
  - [ ] `η_n, η_p, j_n, j_p` 为标量；若是向量，先做平均。
  - [ ] `dUdT_p, dUdT_n` 标量化（mean 或索引）。
  - [ ] `q_nd` = `Q_rxn + Q_rev + Q_ohm` 为标量。
- 热源映射：
  - [ ] `S_q = φ·I1C_total/(L·A_global)`；`q_SI_e = S_q*q_nd_e`。
  - [ ] `q_areal_e = q_SI_e*H` 注入 `variables["heat_source_fields"]`。
- 热装配：
  - [ ] 2D 热方程与 BC 的量纲一致；若使用 Scheme B 无量纲，需在 FT 把 `q_areal_e` 转为无量纲；若使用 SI，则 MT/KT/FT/BC 一致乘 H。


---

## 附录 A. 常见问题与实现细节补全

### A.1 电化学网格 vs 热网格的维度/单元不匹配如何处理？

- 现状：
  - 电化学仍沿用原1D径向/厚度方向离散（SPMe/SPM 的“L2”一维网格：负极/隔膜/正极三段），变量保存在 `variables[...]` 中，典型量是层平均或高斯点场（见 `SPMe_variables`）。
  - 热是二维 Q4 网格（`case.mesh["thermal2D"]`）。
- 处理策略：
  - 我们不做逐点场的1D→2D插值，而是先在电化学侧得到“每层的体积平均热源”（或元素均匀热源），再按面积积分得到各层对应面热源，最后除以单元面积得到单元平均体积内热源， 分配给对应 Q4 元素。
  - 代码位置：`ThermalDistributed.jl/heatQ_Source`。
  - 先计算层平均的 Q_NE、Q_SP、Q_PE、Q_PCC、Q_NCC（单位为体热源，见下文A.5公式）。

  - 这样“一维电化学 → 二维热”的维度不匹配通过“层平均+几何映射”的方式闭合，避免了不必要的插值误差，同时兼顾效率。


### A.3 Jellyroll 的归一化与一维 LGM50 等模型有何区别？

- 电化学时间尺度：
  - 子单元 1C 标尺 I_typ,e = I_1C,total 。同时子单元容量 Q_e ≈ Q_total。因此 ts_e = Q_e / I_typ,e 与全局 ts 相同。故 SPMe 的 ts_p/ts_n/te 无需随单元改变。
- 热的归一化：
  - 当前热模块按“Scheme B 无量纲”实现：热场变量、MT/KT/FT 使用 `case.param_dim.scale` 中的 (ρc_ref, k_ref, L_th, q_ref, h_th, T_ref, t_th) 进行统一缩放；边界中的 Bi = h L_th / k_ref 存为 `scale.h_th`。
  - 一维 lumped 或 LGM50 等非 Jellyroll 模式通常直接使用 lumped 热或1D热，不需要二维的几何映射与并联系统，且不会用到 `jellyroll_*` 的几何工具与分层映射。
  - 结论：电化学无量纲在各模型一致；热在 Jellyroll 使用 2D 尺度与 Bi/Fo 显式缩放，其它 1D 模式要么 lumped 要么 1D，几何与缩放更简单。

### A.4 热模型的边界条件如何设定？

- 外边界（外半径 Rout）：缺省对流边界
  - 维度化形式：−k ∂T/∂n = h (T − T_amb)。
  - 无量纲形式：−(k/k_ref) ∂T*/∂n* = Bi (T* − T_amb*)，其中 Bi = h L_th / k_ref，T* = T/T_ref。
  - 代码：`ThermalDistributed2D_BC` 自动识别 Q4 网格外边界边并组装：
    - 刚度贡献：KT += ∫ h_coeff NᵀN dL*；载荷：FT += ∫ h_coeff T_amb* N dL*；其中 h_coeff=Bi。
  - 内边界（内半径 Rin）：默认绝热；如需内侧对流，可扩展同一套路（函数已有分类 inner/outer 的脚手架）。
- 其他散热（辐射/端面/轴向）
  - 目前未启用辐射项；可选线性化写成等效 h_rad 并叠加到 Bi；端面/轴向散热在2D平面模型中等效为边界条件乘高度 H 的方案，若切 SI，要确保 MT/KT/FT 与 BC 一致乘 H（见第2节建议）。

### A.5 内热源如何计算？给出公式与代码对照

- 层平均（SPM/SPMe）公式（体热源，单位与当前热缩放一致）：
  - 反应热：Q_rxn = a_s · j · η
  - 可逆热：Q_rev = a_s · j · T · dU/dT
  - 欧姆热（简化一致流近似）：
    - 固相：q_ohm,s ≈ (I_app^2) · (t/σ_eff) / 3 ÷ t = I_app^2 / (3 σ_eff)
    - 电解液：q_ohm,e ≈ (I_app^2) · (t/κ_eff) / 3 ÷ t = I_app^2 / (3 κ_eff)（隔膜为 I_app^2 · (t_sp/κ_sp) ÷ t_sp = I_app^2 / κ_sp）
  - 分层合成：
    - 负极：Q_NE = a_s,n j_n η_n + a_s,n j_n T dU_n/dT + q_ohm,s,ne + q_ohm,e,ne
    - 隔膜：Q_SP = q_ohm,e,sp
  - 正极：Q_PE = a_s,p j_p η_p + a_s,p j_p T dU_p/dT + q_ohm,s,pe + q_ohm,e,pe
  - 正集流体：Q_PCC = q_ohm,PCC ≈ I_app^2 / (3 σ_PCC,eff)
  - 负集流体：Q_NCC = q_ohm,NCC ≈ I_app^2 / (3 σ_NCC,eff)
  - 代码对应：`ThermalDistributed.jl/heatQ_Source` 在 SPM/SPMe 分支：
    - 反应/可逆：
      - `Q_rxn_NE = as_n * jn_val * ηn_val`
      - `Q_rev_NE = as_n * jn_val * T_val * dUdT_n`（正极同理）
    - 欧姆热：依据 `I_app_val`、σ_eff、κ_eff、各层厚度 t_n,t_sp,t_p，先得每面积功率 P，再除厚度转体热源，最后求和到 Q_NE/Q_SP/Q_PE。
    - 若是 P2D，则读取高斯点 `jn, ηn, ...` 并通过 `IntV(..., mesh)` 进行加权平均，公式与上面一致但在高斯点层面累积。

- 单元映射与单位：
  - 计算得到的 Q_NE/Q_SP/Q_PE 为体热源（W/m³ 对应 SI；若使用无量纲缩放则为 q*）。
  - 在 `heatQ_Source` 末段：
    - 使用 `material_at` 将层常值赋给各 Q4 元素。
    - 若 `variables` 包含 `thermal2D element current`，则在层内按 I_e·A_e 再分配（保持层总量不变）。
    - 然后依据 `case.param_dim.scale.q_th` 转成 `variables["heat_source_fields"]` 的无量纲值，并标 `heat_source_units = :dimensionless`。
  - 在装配里（`ThermalDistributed2D`），右端 `FT = ∬ q* N dΩ*` 用的就是该无量纲热源；若切换到 SI，需要把 q_areal 或 q_vol 的单位在此处统一（见第2节）。

---

小结：以上补全确保了 1D 电化学与 2D 热的“层平均+几何映射”闭环，明确了并联电气关系的实现近似，给出了 Jellyroll 与其它 1D 模式在热缩放与几何处理上的差异，并提供了可直接对照代码的位置，便于验证与扩展。

### A.6 热网格如何划分？电化学网格如何基于热网格关联？

- 热网格（Q4）生成：
  - 函数：`jellyroll_Q4_mesh(param_dim; nx, ny, gsorder, crop_to_annulus=true, crop_mode=:inscribed)`（见 `src/Jellyrollmodel.jl`）。
  - 步骤：
    1) 在正方域 `[-Rout,Rout]×[-Rout,Rout]` 生成规则 Q4 网格（`SetMesh`）。
    2) 若 `crop_to_annulus=true`，按环域 `[Rin,Rout]` 裁剪：
       - `:inscribed`：仅保留四个节点都在环域内的单元（保证单元“完整”，数值更稳健，适合把一个热单元看作一个子电池）。
       - `:center`：按单元中心点位于环域内裁剪（可得到更满的覆盖，但可能出现跨边界的单元）。
    3) 返回带有高斯点权重与雅可比的 Q4 网格（用于装配 MT/KT/FT 与计算单元面积）。

- 电化学网格与热网格的关系：
  - 电化学仍采用一维“L2”网格（负极/隔膜/正极三段），其离散与 `Q4` 热网格是独立的；二者通过“层判定+层平均映射”耦合，而非从热网格细分出电化学网格。
  - 具体耦合：
    - 用 `jellyroll_element_centers(mesh_th)` 得到每个热单元中心 (x,y)，
    - 用 `p = jellyroll_spiral_params(param_dim)` 与 `material_at(r,θ,p; logic=:rings 或 :spiral)` 判定该中心属于 :NE/:SP/:PE/:PCC/:NCC，
    - 将电化学侧计算出的层平均热源（Q_NE/Q_SP/Q_PE …）赋给相应的热单元；若启用并联分配，再按 I_e·A_e 对层内做守恒重分配（`heatQ_Source`）。
  - 这样设计的好处：
    - 电化学保持低维、稳健与快速；热可在 2D 上捕捉空间非均匀。
    - 如果将来需要更细化的电化学空间分布，可在 `material_at` 的判定基础上，为选定区域引入“局部子模型”，无需改变热网格生成。

---

### A.7 LGM50 与 Northrop 的电化学归一化对照（一致的标尺，不同的数值）

- 共同的无量纲标尺（`src/SetParams.jl -> ChooseCell / NormaliseParam`）：
  - 电流与时间
    - I_typ = cell.I1C（整电池 1C 电流，SI 由参数文件给出）。
    - t0 = 3600 s（默认）；用于把 `opt.Current` 的时间输入换算到物理秒：t_phys = t_nd · t0。
    - 反应通量标尺 j = I_typ / (a0 · L · A_cell)，其中 a0 = 1/r0，r0 默认 1e−6 m；L = t_n + t_sp + t_p；A_cell = cell.area。
  - 扩散/电解时间尺度
    - ts_p = F · cs_max,p · A_cell · L / I_typ
    - ts_n = F · cs_max,n · A_cell · L / I_typ
    - te   = F · ce0      · A_cell · L / I_typ
    - Ds_p = r0²/ts_p，Ds_n = r0²/ts_n，De = L²/te
  - 电位与电导
    - ϕ_ref = R · T_ref / F（把 V 无量纲化）
    - σ_ref = L · I_typ / (ϕ_ref · A_cell)
    - κ_ref = L · I_typ / (ϕ_ref · A_cell)
    - 反应速率标尺：k_p = j / (cs_max,p · √ce0)，k_n 同理
  - 开路电位与熵变项
    - U_nd(θ) = U_dim(θ) / ϕ_ref
    - (dU/dT)_nd = (dU/dT)_dim · (T_ref/ϕ_ref)

- 差异点：
  - LGM50 与 Northrop 的不同之处体现在参数数值（cs_max、ce0、各层厚度、σ、κ、k 等），因此 ts_p/ts_n/te、σ_ref/κ_ref、k_p/k_n 等标尺的数值不同；但归一化公式完全一致。
  - 归一化后在 `NormaliseParam` 中：
    - 厚度、扩散、导电、反应常数等都按上述标尺转为无量纲；
    - OCP 与 dUdT 用闭包重新包装成无量纲函数；
    - `param.cell.area` 被归一化为 1（仅用于电化学尺度），而 `param_dim.cell.area` 保持 SI，供尺度计算。

代码对照：
- `src/SetParams.jl` 中 `ChooseCell` 设置 `param_dim.scale.*`；`NormaliseParam` 生成 `param`（无量纲）。
- `src/SPMe.jl/SPMe_variables` 中 `I_app = opt.Current(t*t0) / scale.I_typ` 得到无量纲电流；温度 T 亦为无量纲（T/T_ref）。

### A.8 Jellyroll + 2D 分布热 + “每热单元=子SPMe”的归一化区别

- 目标方案（物理一致的“每热单元=子SPMe”）：
  - 定义：单元面积 A_e；全局有效受流面积 A_global（通常取 `param_dim.cell.area` 或基于 Jellyroll 有效区域）。
  - 子单元 1C 标尺：I_typ,e = I_1C,total 。
  - 子单元归一化电流：I_app = I_e / I_typ,e。
  - 关键结论：时间尺度不变
    - 子单元容量 Q_e ∝ A_e，且 I_typ,e ∝ A_e ⇒ ts_e = F·cs_max·A_e·L / I_typ,e = ts_global。
    - 因此 SPMe 的 ts_p、ts_n、te 与全局相同，无需为每个单元更换扩散/电解时间标尺。
  - 温度：T_e = T_e / T_ref（T_e 为该热单元节点均温）。
  - 算法（单步强耦合示意）：
    1) 分电流非线性求解得到 I_e。
    2) 计算 I_app 并构造每单元的输入（含 T_e,e）。
    3) 用同一组无量纲尺度运行子 SPMe，得到每单元的 η、j 等，并计算 q_e,e。
    4) 用 q_ref 或 S_q 恢复到 SI 或热模块无量纲，再进入 2D 热一步。

热边界应为外边界对流，内边界绝热，极耳边界是在该条件的基础上对其进行覆盖
添加极耳的边界条件，通过识别极耳所在位置，赋予极耳所涉及的单元改变其温度与极耳温度相同，极耳温度随时间线性增长
极耳的实现形式为如图所示，tab.width = 40e-3，正极耳位于内螺旋线，负极耳位于外螺旋线上，以极耳宽度作为螺旋线长度计算角度，改变粘接到的节点的温度，默认为正极耳一个在起点，负极耳一个在终点。
在网格生成阶段为边界边添加标签/标记（如果网格工具支持），然后在 BC 中直接使用这些边界标签。

Jellyroll.jl标记正负极极耳的起点
Jellyrollmodel.jl新增极耳边界标记函数：根据参数提供的正负极极耳起点和极耳宽度，计算螺旋线的弧长，通过弧长识别角度，标记极耳边界，极耳抽象为没有厚度的弧线，正极耳只影响内螺旋线r_in(θ) = a + bθ + s_in（取 s_in = 0, 对应一周期的起始边，含 PCC 内侧）上的节点，负极耳只影响外螺旋线r_out(θ) = a + bθ + s_out（取 s_out = t_repeat）上的节点，标记对应极耳宽度影响的节点。
ThermalDistributed.jl调用极耳边界标记函数，在BC中识别极耳边界，并赋予极耳边界的节点温度与极耳温度相同，极耳温度随时间线性增长。
