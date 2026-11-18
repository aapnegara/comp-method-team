# 修改原因、边界实现与 Assignment 对应说明

[English version](README_Changes_and_Assignment_EN.md)

本文说明最终交付的物理设置、边界实现、数值配置及选择理由。

## 1. 任务范围与依据

当下任务是 Setup & boundaries：改编 mod_usr.t 的初始条件和 jet inflow，并解释 boundary conditions 与 ghost cells。
CMfAA_HPCassignment.pdf 第2问要求说明几何、初始条件、归一化；第3问要求说明边界如何填 ghost cells、每侧几层、为什么；第4问问离散方法与 stencil。
本说明重点对应初始条件、边界与虚单元，以及它们和数值重构方法之间的关系。
论文依据是 Massaglia et al. 2016, A&A 596 A12，第2.2节公式(5)-(8)与表1 Case B。

原始示例位于 /home/peder/codes/amrvac/tests/hd/jet_cloud。
它是喷流撞偏轴云团的问题，入口沿 +x；源码虽然有 IFTHREED 片段，坐标设置仍写 Cartesian_2D。
因此“改为三维”必须同时修改几何、第三分量、喷口半径、边界数组和参数文件，不能只改 setup.pl -d。
正式文件在包根目录，原始示例放在 original/ 供对照。
原因：建立可比较、可复查的独立基准，而不是覆盖已有组员工作。

## 2. 方程、变量与物理参数（老师第1问）

采用无磁场、无重力、无辐射冷却的理想气体 Euler 方程，加一个被动示踪量：

$$
\begin{aligned}
\frac{\partial \rho}{\partial t}+\nabla\cdot(\rho\mathbf{v})&=0,\\
\frac{\partial(\rho\mathbf{v})}{\partial t}
+\nabla\cdot\left(\rho\mathbf{v}\otimes\mathbf{v}+p\mathbf{I}\right)&=0,\\
\frac{\partial E}{\partial t}+\nabla\cdot\left[(E+p)\mathbf{v}\right]&=0,\\
\frac{\partial C}{\partial t}+\nabla\cdot(C\mathbf{v})&=0.
\end{aligned}
$$

其中 $\mathbf{I}$ 为单位张量，$\otimes$ 表示向量外积。示踪密度与总能量密度为：

$$
C=\rho f,\qquad
E=\frac{p}{\gamma-1}+\frac{\rho}{2}\left(v_x^2+v_y^2+v_z^2\right),
\qquad \gamma=\frac{5}{3}.
$$

守恒数组为：

$$
\mathbf{U}=(\rho,m_x,m_y,m_z,E,C)^{\mathsf T},\qquad
\mathbf{m}=\rho\mathbf{v}.
$$

物理原始量是 $\rho,v_x,v_y,v_z,p,f$；但本安装 EOS 转换不处理 tracer 槽，因此该槽必须自行保持 $C$。
为什么加示踪量：密度变化既可能来自压缩也可能来自物质混合，$f$ 能区分喷流和环境来源。

本基准参数：

$$
\rho_j=1,\qquad \eta=\frac{\rho_j}{\rho_a(0)}=0.01,\qquad
M=4,\qquad r_c=40,\qquad \alpha=2.
$$

$$
p_j=p_a=\frac{1}{\gamma}=0.6,\qquad
c_{s,j}=\sqrt{\frac{\gamma p_j}{\rho_j}}=1,\qquad
v_j=M c_{s,j}=4.
$$

为什么用喷流声速定义 Mach：避免把原 jet_cloud 的环境声速归一化直接搬过来，造成同一个 $M$ 对应不同物理量。
原 $\beta$ 云团参数、原云半径与原环境声速参数被移除，因为这里没有云团。
$r_c$ 的含义现在是环境核心半径，不是原示例云团半径；已在源码注明。

## 3. 几何、单位与初始条件（老师第2问）

三维笛卡尔域为 $x,z\in[-4,4]$、$y\in[0,8]$；喷流沿 $+y$。
长度单位为 $r_j$，密度单位为 $\rho_j$，速度单位为 $c_{s,j}$。其他归一化单位为：

$$
t_0=\frac{r_j}{c_{s,j}},\qquad
p_0=E_0=\rho_j c_{s,j}^{\,2},\qquad
M_0=\rho_j r_j^3.
$$

这里 $E_0$ 是能量密度单位，$M_0$ 是质量单位，不是 Mach 数。
代码按无量纲计算，未把真实物理单位写入 AMRVAC 的 `unit_*`；不能把日志中的单位默认值当作 pc、年或 K。
若采用论文 $r_j=100\,\mathrm{pc}$，则 $r_c=40r_j=4\,\mathrm{kpc}$；把时间换算成年还需指定物理 $c_{s,j}$。

环境球半径及初始内部状态为：

$$
\begin{aligned}
r&=\sqrt{x^2+y^2+z^2},\\
\rho_a(r)&=\frac{100}{1+(r/40)^2},\\
p&=0.6,\qquad v_x=v_y=v_z=0,\qquad C=0.
\end{aligned}
$$

注意：球半径 $r$ 用于环境，横向半径 $R=\sqrt{x^2+z^2}$ 用于喷口，不是同一个变量。
为什么删除原先预填喷流：论文初始域是环境，喷流从边界随时间进入。
为什么删除偏轴云团：论文模型的目标是分层环境内传播，不是喷流撞云。
为什么密度分层仍不加重力：该模型压力均匀，静止 Euler 动量方程没有压力梯度，解析上已经平衡；温度随 $1/\rho$ 变化。
代码在 `initonegrid_usr` 先填 $\rho$、速度和压力，再调用 `eos%to_conserved` 一次。
初始 $C=0$；不能依赖 EOS 自动把任意非零 $f$ 乘成 $\rho f$。
当前只实现 $\alpha=2$ 的 Case B；不要把 $\alpha$ 改成 $0$ 就当作论文均匀 Case A，因为原公式会多出分母 $2$，需要单独分支。

## 4. 六个物理外边界（老师第3问核心）

编号规则：对方向 $d$，最小面的编号为 $i_B=2d-1$，最大面的编号为 $i_B=2d$；代码变量名是 `iB`。

| 计算域面 | 边界编号 | 类型 |
|---|---|---|
| $x_{\min}$ | $i_B=1$ | `cont` |
| $x_{\max}$ | $i_B=2$ | `cont` |
| $y_{\min}$ | $i_B=3$ | `special`，由 `specialbound_usr` 填写 |
| $y_{\max}$ | $i_B=4$ | `cont` |
| $z_{\min}$ | $i_B=5$ | `cont` |
| $z_{\max}$ | $i_B=6$ | `cont` |

为什么把 special 从原例 case(1) 改成 case(3)：入口从 x_min 改到 y_min，必须与 .par 同时更改。
三维守恒流体变量需5个边界条目：$\rho,m_x,m_y,m_z,E$。
tracer 在本版本注册 need_bc=.false.，内置边界规则继承密度规则；special 中仍显式填 C。
因此参数写 5*'cont' 或 5*'special'；不是二维的4项，也不是因为有示踪量就随意写6项。

### 4.1 cont 是什么？
它把最近一层内部单元的值复制到每一层外部 ghost cell。
例如最大 $x$ 面，内部末单元为 $N$：

$$
\mathbf{U}_{N+1}=\mathbf{U}_{N},\qquad
\mathbf{U}_{N+2}=\mathbf{U}_{N}.
$$

最小 $x$ 面内部从 $1$ 开始：

$$
\mathbf{U}_{0}=\mathbf{U}_{1},\qquad
\mathbf{U}_{-1}=\mathbf{U}_{1}.
$$

这给出离散的法向零梯度，不是逐层镜像。
为什么用它：近似开放外边界，与论文对其他面的零梯度设定一致。
它不是严格的无反射特征边界，也不禁止回流；当边界附近法向速度指向域内时仍可能流入。
原 jet_cloud 的横向 noinflow 会裁去某些向内法向动量，故不能不加解释地用它代替论文的 cont。
真正的长期计算需扩大域或做边界敏感性测试，检查外边界是否影响喷流演化。

### 4.2 入口平面喷口外如何反射？

入口平面为 $y=0$。先对整个 `y_min` ghost 区构造镜像状态 $Q_r$。
若 ghost 的 $y$ 索引为 $1,2$，内部从 $3$ 开始，则远层 ghost $1$ 对应内部 $4$，近层 ghost $2$ 对应内部 $3$：

```text
... [ghost 1] [ghost 2] | y=0 | [inside 3] [inside 4] ...
```

一般镜像索引公式为：

$$
j_{\mathrm{mirror}}=2j_{g,\max}+1-j,
\qquad j_{g,\max}=\texttt{ixOmax2}.
$$

对应源码表达式为 `jmirror=2*ixOmax2+1-j`。
为什么不能两层都复制内部3：那是常数外推，不是反射，改变重构模板中的梯度。
对 $\rho,m_x,m_z,E,C$ 作偶对称复制；只对法向 $m_y$ 作奇对称：

$$
(m_y)_g=-(m_y)_{\mathrm{mirror}}.
$$

为什么切向动量不反号：这是无粘滑移壁的反射条件，不是 no-slip 固壁。
为什么 $E$ 不反号：法向速度反号不改变其平方，动能和内能均保持。
为什么 $C$ 直接镜像：它是随密度携带的标量密度，反射不能改变物质来源比例。
当前守恒形式处理无需对镜像状态再次 to_conserved；重复转换会错误地再次乘密度。

### 4.3 圆喷口和平滑连接
固定轴上喷流状态为：

$$
\rho=1,\qquad v_x=v_z=0,\qquad m_y=\rho_jv_j=4,
\qquad p=0.6,\qquad f=1.
$$

代码 `R=hypot(x,z)` 是 $R=\sqrt{x^2+z^2}$ 的数值稳健计算。
按论文式（8）：

$$
Q=Q_r+(Q_j-Q_r)S(R),\qquad
S(R)=\frac{1}{\cosh\!\left[\left(\frac{R}{R_s}\right)^n\right]}.
$$

各变量的平滑参数分别为：

$$
(R_s,n)=
\begin{cases}
(1.4,8), & Q=\rho,\\
(1,6), & Q\in\{v_x,m_y,v_z,p,f\}.
\end{cases}
$$

为什么两种权重：遵循论文的入口动量剖面选择，不能统一替换成相同的锐利 top-hat。
轴上 $R=0$、$S=1$，完全得到规定喷流；远处 $S\to0$，回到反射状态。
因此“喷口内固定、外面反射”是物理概述；实际过渡区既不是全固定，也不是严格纯反射。
变量 $Q_r$ 每次从当时内部场取得，而不是永远用初始环境，因为反射必须响应随时间演化的流动。

具体重建顺序：

1. 从镜像 $\mathbf{U}$ 求 $\rho_r,v_{x,r},v_{z,r},m_{y,r},p_r$ 和 $f_r=C_r/\rho_r$。
2. 按上述平滑规则混合 $\rho,v_x,v_z,m_y,p,f$。
3. 重建动量、轴向速度、总能量和示踪密度：

$$
\begin{aligned}
m_x&=\rho v_x,\qquad m_z=\rho v_z,\qquad v_y=\frac{m_y}{\rho},\\
E&=\frac{p}{\gamma-1}
+\frac{1}{2}\left[\rho\left(v_x^2+v_z^2\right)+\frac{m_y^2}{\rho}\right],\\
C&=\rho f.
\end{aligned}
$$

为什么混合 $m_y$ 而不是 $v_y$：这是论文明确指定的轴向动量平滑。
为什么不直接混合 $E$：论文平滑压力和不同形式的速度/动量，直接混合总能量一般不等价。
为什么必须包含 $v_z$：漏掉第三分量的动能会在有横向运动时算错压力。
为什么先限制幂运算底数再算 cosh：防止远场指数和 cosh 溢出。
最大 $\cosh$ 参数设为 $80$，对应权重约为 $1/\cosh(80)\approx3.61\times10^{-35}$；这是极小数值尾部截断，并非给喷口新增物理边缘。
为什么不人为加入入口扰动：论文这里没有显式入口扰动；真实三维非轴对称增长可由数值噪声启动。
本短跑不证明这些模态已经充分发展。

## 5. 虚单元层数及选择理由

虚单元层数由空间重构模板决定，与空间维度、Runge–Kutta 阶段数没有直接对应关系。当前参数配置为：

| 参数文件 | 重构方法 | 每侧虚单元层数 | 选择理由 |
|---|---|---|---|
| `smoke.par`、`demo.par`、`minmod2.par` | minmod | 2 | 两层能够覆盖分段线性重构的边界通量模板。 |
| `mp5_3.par` | MP5 | 3 | MP5 使用更宽的空间模板，本版本 AMRVAC 要求至少三层。 |

层数在 `.par` 的 `boundlist` 中由 `nghostcells` 指定。用户模块支持上述两种配置，入口物理条件与逐层镜像规则相同。

**minmod 为什么使用两层？** 在均匀网格上，受限变化量 $\delta q_i$ 与斜率 $s_i$ 为：

$$
\delta q_i=\operatorname{minmod}\left(q_i-q_{i-1},\ q_{i+1}-q_i\right),
\qquad s_i=\frac{\delta q_i}{\Delta x}.
$$

单个单元的斜率只使用相邻单元，但界面通量需要左右两个单元的重构状态。计算边界外第一层虚单元的斜率时，还需要第二层虚单元，因此采用两层。更新内部单元 $i$ 的完整空间模板可涉及 $i-2$ 到 $i+2$。

**MP5 为什么使用三层？** MP5 用更宽的模板构造界面状态，因此还需要第三层域外数据。它可在光滑区域保留更高阶的空间信息；选用三层是为满足该重构方法的模板需求。时间推进的阶段数不改变这一空间要求。

**虚单元如何填充？** 计算域外边界按第4节的 `cont` 或 `special` 规则填充；域内块边界通过相邻块数据复制或 MPI 通信填充。逐层镜像公式根据 ghost 区的索引范围工作，因此两层和三层配置都对应各自的内部镜像层。

每块采用 $8^3$ 个有效单元，便于小规模三维算例进行块分解。包含虚单元的局部矩形数组大小分别为：

$$
N_{\mathrm{local}}=
\begin{cases}
(8+2\times2)^3=12^3=1728, & \text{minmod, 2 ghost layers},\\
(8+2\times3)^3=14^3=2744, & \text{MP5, 3 ghost layers}.
\end{cases}
$$

这属于每块局部存储开销，不增加物理有效单元数。固定 `refine_max_level=1` 使算例采用单层网格，边界处理不涉及 AMR 粗细层接口。

## 6. 数值方法与参数选择（老师第4问，以及第5问的CFL部分）

基准采用固定网格有限体积、TVDLF 通量、minmod 分段线性重构与 `threestep`；`mp5_3.par` 保持通量和时间推进相同，改用 MP5 重构及三层虚单元。
本安装 threestep 默认是 ssprk3，ssprk_order 默认3；默认空间更新为 unsplit。
SSPRK(3,3) 可写：

$$
\begin{aligned}
\mathbf{U}^{(1)}&=\mathbf{U}^{n}+\Delta t\,L\!\left(\mathbf{U}^{n}\right),\\
\mathbf{U}^{(2)}&=\frac{3}{4}\mathbf{U}^{n}
+\frac{1}{4}\left[\mathbf{U}^{(1)}+\Delta t\,L\!\left(\mathbf{U}^{(1)}\right)\right],\\
\mathbf{U}^{n+1}&=\frac{1}{3}\mathbf{U}^{n}
+\frac{2}{3}\left[\mathbf{U}^{(2)}+\Delta t\,L\!\left(\mathbf{U}^{(2)}\right)\right].
\end{aligned}
$$

这里 $L$ 是空间离散算子，上标 $n$ 表示当前时间层，$(1)$、$(2)$ 表示中间阶段。
为什么选择它：结构清晰、便于先检验入口与边界；相比原例高阶 mp5/fivestep 更容易定位问题。
TVDLF 可用中心平均通量加局部最大波速控制的耗散来理解；它比接触分辨更强的HLLC更扩散。
因此输出图中的扩散混合不能全部解释为真实湍流混合。
时间形式三阶不代表整体三阶：minmod 空间平滑区通常二阶，极值/间断附近退化，整体受空间与非光滑解限制。
无额外物理源项，所以本问题无重力/冷却的分裂或非分裂源项积分，也无这些源项引入的额外时间步限制。

取 $C_{\mathrm{CFL}}=0.3$，代码设置为 `typecourant='maxsum'`。
对均匀笛卡尔网格，核心约束近似：

$$
\Delta t\leq
\frac{C_{\mathrm{CFL}}}{
\displaystyle\max_{\mathrm{cells}}\left[
\frac{|v_x|+c_s}{\Delta x}
+\frac{|v_y|+c_s}{\Delta y}
+\frac{|v_z|+c_s}{\Delta z}
\right]}.
$$

声速及 Euler 法向特征速度为：

$$
c_s=\sqrt{\frac{\gamma p}{\rho}},\qquad
\lambda\in\{v_n-c_s,\ v_n,\ v_n+c_s\}.
$$

其中 $v_n$ 是法向速度，接触、横向和示踪扰动随 $v_n$ 输运。
输出时刻和终止时刻还可能截短实际dt，因此不能只拿单一轴上状态算出的值当全部步长。
选择0.3而非原例0.8是为了保守的接口测试，不表示0.8必然错误。
small_values_method='error' 是避免发现问题后用平均修补掩盖它；输出的正值和有限性由使用者自行检查。
dtmin=1e-12 用于异常小步长终止，不是人为规定每步必须取这么小。
refine_max_level=1 是主动隔离AMR影响；日志即使列出默认细化估计器也不表示实际用了多层AMR。

网格与范围为何这样定：

| 配置 | 有效网格 | 网格间距 | 每喷流半径的单元数 | 终止时间 | 输出间隔 |
|---|---|---|---|---|---|
| smoke / minmod2 / mp5_3 | $32^3$ | $\Delta x=\Delta y=\Delta z=0.25$ | 4 | $t_{\mathrm{end}}=0.2$ | $0.1$ |
| demo | $64^3$ | $\Delta x=\Delta y=\Delta z=0.125$ | 8 | $t_{\mathrm{end}}=1$ | $0.5$ |

这些配置采用小域与较短演化时间，降低本机计算与输出成本；demo 增加入口分辨率用于观察早期演化。
输出采用 DAT + VTU，前者可重启，后者方便ParaView；交付包不预装数据，VTU由使用者运行后生成。
不同终止时刻的两次运行不是严格的网格收敛比较，不据此计算收敛阶。

## 7. 函数逐项导读

usr_init：保留原例五个回调注册与hd_activate；改Cartesian_3D，让网格与速度有三个方向。
initglobaldata_usr：集中Case B参数并检查维度、gamma、能量方程、tracer、入口位置、aux输出数与ghost层数。
  原因：配置不匹配时直接说明错误，而不是生成貌似正常但物理错误的结果。
initonegrid_usr：全域分层静止环境；调用一次EOS转换；理由见第3节。
specialbound_usr：case(3)逐层镜像、三维圆喷口平滑、完整能量重建；理由见第4节。
extra_var_output：保留 $p/\rho$ 温度代理；删除未使用的冷却模块导入，避免误以为启用了冷却。
extra_var_names_output：仍命名T，与nwauxio=1匹配。
代码中 $i,j,k$ 对应 $x,y,z$；`^L`、`^S` 是 LASY 索引宏，编译时自动展开为 `mod_usr.f`。

## 来源与版本

原始示例：AMRVAC `tests/hd/jet_cloud`，包内 `original/` 保存两份原文件。
物理依据：Massaglia et al. (2016), A&A 596 A12，第2.2节、公式(5)—(8)、表1 Case B。
作业依据：课程 `CMfAA_HPCassignment.pdf` 的第1—11问；本交付重点是第2、3问及与第4问的衔接。
原始实现见 `original/mod_usr.t`，最终实现见包根目录的 `mod_usr.t`。
实现接口以 AMRVAC commit `4ac04b2672f11ba5b8909bed14427aefd1d60425` 为依据。
其中每侧两层默认值及MP5至少三层的调整来自该版本 `src/io/mod_input_output.t`；内置外边界规则在 `src/mod_boundary_conditions.t`；tracer注册与转换约定见 `src/hd/mod_hd_phys.t` 和 `src/hd/mod_hd_eos.t`。
包内不附带求解器源码副本，编译使用接收者自己的 AMRVAC 安装。
