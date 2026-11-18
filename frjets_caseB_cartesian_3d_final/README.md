# 三维笛卡尔 FR jet：Setup & boundaries 交付

[English version](README_EN.md)

本包提供从 AMRVAC `jet_cloud` 改编的三维 Case B 设置：初始分层环境、沿 +y 的圆喷口、平滑入口与六个面的边界条件。
配套小域参数、独立截面绘图程序和实现说明。

## 包内文件

```text
frjets_caseB_cartesian_3d_final/
├── mod_usr.t
├── smoke.par
├── demo.par
├── minmod2.par
├── mp5_3.par
├── original/
│   ├── mod_usr.t
│   └── amrvac_3d.par
├── plot_demo_slice.py
├── README.md
├── README_EN.md
├── README_修改原因与Assignment.md
└── README_Changes_and_Assignment_EN.md
```

| 文件 | 用途 |
|---|---|
| `mod_usr.t` | 正式三维用户模块：参数、初始化、特殊入口、温度代理输出。 |
| `smoke.par` | 32³ 单元的小规模配置，终止 t=0.2，每隔0.1输出；方便首次运行。 |
| `demo.par` | 64³ 单元的演示配置，终止 t=1，每隔0.5输出；用于生成 demo 截面。 |
| `minmod2.par` | 32³、minmod、每侧两层虚单元；与 smoke 的配置相同，使用独立输出前缀 `output/minmod2`。 |
| `mp5_3.par` | 32³、MP5、每侧三层虚单元；终止 t=0.2，输出前缀 `output/mp5_3`。 |
| `original/mod_usr.t` | 原始 jet_cloud 用户模块，供对照阅读。 |
| `original/amrvac_3d.par` | 原始 jet_cloud 的三维参数，供对照阅读；不要与根目录的正式源码混用。 |
| `plot_demo_slice.py` | 读取一份运行生成的 VTU，画密度、示踪比例、轴向速度截面 |
| `README.md` | 文件说明、编译运行和绘图方法。 |
| `README_EN.md` | 上述文件与运行说明的英文版。 |
| `README_修改原因与Assignment.md` | 物理与数值设置的选择理由、初始条件、边界和虚单元解释。 |
| `README_Changes_and_Assignment_EN.md` | 上述修改原因与实现说明的英文版，保留公式。 |

## 物理与数值配置

三维笛卡尔域：x,z∈[-4,4]，y∈[0,8]。长度以喷流半径为单位。
喷流沿 +y，rho_j=1、M=4、eta=0.01、rc=40、alpha=2、gamma=5/3，p_j=p_a=0.6。
基准采用 TVDLF + minmod + threestep，CFL=0.3，固定单层网格，每侧两层 ghost cells。`mp5_3.par` 使用 MP5 与每侧三层 ghost cells，其他数值选项相同。
这些是本交付的运行基准；正式数值方案需要与 Yitse/Agung 统一。

代码依据的本机 AMRVAC commit 为 `4ac04b2672f11ba5b8909bed14427aefd1d60425`。
需要已经安装的 MPI-AMRVAC 及其构建依赖（Perl、make、Fortran/MPI 编译器和 MPI 启动器）。
其他 AMRVAC 版本需核对回调、EOS 与 tracer 约定。

## 编译运行

下面代码在 Linux / WSL Ubuntu 的 Bash 中执行。把路径替换为自己的解压目录和 AMRVAC 安装位置。

```bash
cd /mnt/d/2026-2027/computatianal/teamwork/frjets_caseB_cartesian_3d_final
export AMRVAC_DIR=/home/peder/codes/amrvac
export OMP_NUM_THREADS=1
set -euo pipefail
mkdir -p output

"$AMRVAC_DIR/setup.pl" -d=3 -v=3
touch amrvac.h
make -j2

# 按需要运行一项或两项；只在编译成功后执行。
mpirun -np 1 ./amrvac -i smoke.par
mpirun -np 2 ./amrvac -i demo.par

# 如需使用另外两份参数，单独选择执行：
mpirun -np 2 ./amrvac -i minmod2.par
mpirun -np 2 ./amrvac -i mp5_3.par
```

`-d=3 -v=3` 指三个空间方向和三个向量分量。`.t` 中的 LASY 宏通过构建流程展开，不能直接 gfortran 编译。
`touch amrvac.h` 提供本机构建流程所需的用户头文件；`make -j2` 控制编译并行度，`mpirun -np 2` 控制运行进程数。
修改 `mod_usr.t` 后重新 make；修改 `.par` 通常无需重新编译。
同名前缀输出可能被覆盖，重新运行前自行保存旧结果或修改 `base_filename`。

运行生成的 `amrvac`、`mod_usr.f`、`*.mod`、makefile 等是构建产物；`output/*.dat` 是快照，`output/*.vtu` 用于可视化，`output/*.log` 是求解器日志。这些文件不预装在交付包内。

## 生成 demo_slice.png

绘图需要 Python 3、NumPy 和 Matplotlib。Windows 与 WSL 的 Python 环境相互独立；Windows 已有库不代表 WSL 的 `/usr/bin/python3` 也有。绘图可以在任一侧运行，AMRVAC 的编译运行仍按上文在 WSL/Linux 中执行。

### Windows PowerShell 绘图

本机 Windows 已有包含 NumPy、Matplotlib 的 Python，可以退出 WSL 后直接运行。路径按解压位置修改：

```powershell
# 在 Windows PowerShell 中执行。
cd D:\2026-2027\computatianal\teamwork\frjets_caseB_cartesian_3d_final
python plot_demo_slice.py
```

若 PowerShell 的 `python` 指向其他环境，本机可明确指定已配置的解释器：

```powershell
& "C:\Users\Lenovo\miniconda3\python.exe" plot_demo_slice.py
```

其他机器在自己的 Windows Python 环境中安装依赖时可使用 `python -m pip install numpy matplotlib`。
WSL 的 `/mnt/d/...` 就是 Windows 的 `D:\...`，同一份 VTU 无需重新复制或运行模拟。

### WSL / Linux 绘图

推荐在独立虚拟环境中安装，避免修改系统 Python：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install numpy matplotlib
.venv/bin/python plot_demo_slice.py
```

如果创建环境时提示缺少 venv/ensurepip，Ubuntu 可先运行 `sudo apt install python3-venv`（需要时先执行 `sudo apt update`），再重试。不要在没有 pip 的系统 Python 上反复执行相同的 pip 命令。
已安装 uv 的环境也可使用 `uv venv .venv` 和 `uv pip install --python .venv/bin/python numpy matplotlib`。

本机另外已配置了包目录之外的绘图环境，可以直接使用：

```bash
/mnt/d/2026-2027/computatianal/work/frjets_plot_env/bin/python plot_demo_slice.py
```

该环境只存在于本机，不随压缩包交付。下面的可选命令使用 `python3` 举例；应替换为你已配置的解释器（Windows 为 `python`，虚拟环境为 `.venv/bin/python`）。

默认读取 `output/demo0002.vtu`，生成 `demo_slice.png`。只需要这一份 VTU，不需要 smoke、MPI 对照或其他测试文件。
也可以选择输入、图片位置和目标截面：

```bash
python3 plot_demo_slice.py output/demo0001.vtu --output demo_t05.png
python3 plot_demo_slice.py output/demo0002.vtu --z 0.5 --output demo_z05.png
```

默认 `--z 0`，选取最接近该坐标的已有单元中心平面；若距离相同，取正侧。64³配置中通常对应 z=0.0625。
图片标题读取文件实际时间，并显示实际截面坐标，不把单元中心切片冒称为插值后的 z=0 平面。
绘图程序针对本配置产生的未压缩、LittleEndian、UInt32 长度头、Float32 坐标与字段、Int32 八顶点单元的 `vtuBCCmpi` VTU，支持多个 Piece。
它要求固定均匀笛卡尔网格，改成 AMR 或其他输出格式后建议使用 ParaView，或另行修改读取方式。
程序检查的是文件结构是否适合绘图，不承担物理解的质量判定。

## 结果变量与 ParaView

在 ParaView 打开 `output/demo*.vtu`，Apply。用 Slice 设置法向 (0,0,1)、原点 (0,0,0)，可看纵向切片；隐藏原始不透明 Surface 避免挡住内部。

- `rho`：密度。
- `v1/v2/v3`：x/y/z 速度，轴向速度是 `v2`。
- `p`：压力。
- `trp1`：本机版本中仍为 C=rho*f；Calculator 用 `trp1/rho` 得到喷流比例 f。
- `T`：p/rho，无量纲温度代理，不直接对应 Kelvin。

初始内部区域都是静止环境，示踪量为零；喷流随后从底部中央进入。t=1 只展示早期注入，小域短时结果不能用于直接判定长期 FR I 形态。
脚本输出的三张图依次是密度、喷流比例和轴向速度。比较不同时间时，注意不同图片可能自动采用不同颜色范围。
