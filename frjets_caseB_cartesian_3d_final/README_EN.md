# 3D Cartesian FR jet: setup and boundaries

[中文版](README.md) | [Implementation rationale and assignment notes](README_Changes_and_Assignment_EN.md)

This package adapts the AMRVAC `jet_cloud` example to a three-dimensional setup with the physical parameters of Case B: an initially stratified atmosphere, a circular jet nozzle directed along $+y$, a smooth inlet profile, and boundary conditions on all six domain faces.

It contains small-domain parameter files, a standalone slice plotting script, and implementation notes.

## Package contents

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

| File | Purpose |
|---|---|
| `mod_usr.t` | The adapted 3D user module: parameters, initialization, special inlet boundary, and temperature-proxy output. |
| `smoke.par` | A $32^3$ configuration ending at $t=0.2$, with output every $0.1$ time units; intended for a first run. |
| `demo.par` | A $64^3$ configuration ending at $t=1$, with output every $0.5$ time units; used to generate the demo slice. |
| `minmod2.par` | $32^3$, minmod, and two ghost layers per side. It has the same configuration as `smoke.par`, with the separate output prefix `output/minmod2`. |
| `mp5_3.par` | $32^3$, MP5, and three ghost layers per side; ends at $t=0.2$ with output prefix `output/mp5_3`. |
| `original/mod_usr.t` | The original `jet_cloud` user module, retained for comparison. |
| `original/amrvac_3d.par` | The original example's 3D parameters, retained for comparison. Do not combine this file with the adapted root-level module. |
| `plot_demo_slice.py` | Reads one VTU snapshot and plots slices of density, jet tracer fraction, and axial velocity. |
| `README.md` / `README_EN.md` | Chinese / English file descriptions, build instructions, and plotting instructions. |
| `README_修改原因与Assignment.md` / `README_Changes_and_Assignment_EN.md` | Chinese / English explanations of the physical and numerical choices, initialization, boundaries, and ghost cells. |

## Physical and numerical configuration

The domain is Cartesian, with $x,z\in[-4,4]$ and $y\in[0,8]$. Lengths are measured in units of the jet radius. The jet propagates along $+y$, with

$$
\rho_j=1,\quad M=4,\quad \eta=0.01,\quad r_c=40,\quad
\alpha=2,\quad \gamma=\frac{5}{3},\quad p_j=p_a=0.6.
$$

The baseline uses TVDLF fluxes, minmod reconstruction, and `threestep` time integration, with CFL number $0.3$, one fixed grid level, and two ghost layers per side. `mp5_3.par` instead uses MP5 and three ghost layers; the remaining numerical options are the same.

These are the supplied baseline configurations. Coordinate the production numerical scheme with Yitse/Agung.

The implementation targets AMRVAC commit `4ac04b2672f11ba5b8909bed14427aefd1d60425`. An existing MPI-AMRVAC installation is required, including Perl, make, a Fortran/MPI compiler, and an MPI launcher. If using another AMRVAC version, check callback interfaces, the equation-of-state (EOS) conversions, and tracer conventions.

## Build and run

Run the following in Bash on Linux or WSL Ubuntu. Replace the package and AMRVAC paths with those on your machine. The paths shown below are examples from the development machine.

```bash
cd /mnt/d/2026-2027/computatianal/teamwork/frjets_caseB_cartesian_3d_final
export AMRVAC_DIR=/home/peder/codes/amrvac
export OMP_NUM_THREADS=1
set -euo pipefail
mkdir -p output

"$AMRVAC_DIR/setup.pl" -d=3 -v=3
touch amrvac.h
make -j2
```

After a successful build, choose the run or runs you need:

```bash
mpirun -np 1 ./amrvac -i smoke.par
mpirun -np 2 ./amrvac -i demo.par

# Alternative parameter files; run these separately if needed:
mpirun -np 2 ./amrvac -i minmod2.par
mpirun -np 2 ./amrvac -i mp5_3.par
```

`-d=3 -v=3` selects three spatial dimensions and three vector components. The `.t` file contains LASY macros expanded during the build; do not compile it directly with gfortran.

`touch amrvac.h` provides the user header required by this build workflow. `make -j2` controls compilation parallelism, while `mpirun -np 2` selects two MPI processes for the simulation. `OMP_NUM_THREADS=1` keeps the OpenMP thread count at one.

Rebuild with `make` after changing `mod_usr.t`. Editing a `.par` file usually does not require rebuilding. Reusing an output prefix may overwrite existing results, so archive them or change `base_filename` before rerunning.

Files such as `amrvac`, `mod_usr.f`, `*.mod`, and the makefile are generated build artifacts. `output/*.dat` contains AMRVAC snapshots, `output/*.vtu` is for visualization, and `output/*.log` contains solver diagnostics. These generated files are not included in the package.

## Generate demo_slice.png

Plotting requires Python 3, NumPy, and Matplotlib. Windows and WSL have separate Python environments: a library installed on Windows is not automatically available to `/usr/bin/python3` in WSL. The plotting script can run on either platform; build and run AMRVAC on Linux/WSL as described above.

### Windows PowerShell

Use a Windows Python environment with NumPy and Matplotlib installed. From PowerShell, adjust the directory to your extracted package:

```powershell
cd D:\2026-2027\computatianal\teamwork\frjets_caseB_cartesian_3d_final
python plot_demo_slice.py
```

If necessary, install the plotting dependencies into that same Python environment:

```powershell
python -m pip install numpy matplotlib
```

You can also specify an interpreter explicitly. For example, the development machine has:

```powershell
& "C:\Users\Lenovo\miniconda3\python.exe" plot_demo_slice.py
```

This interpreter path is machine-specific; replace it with your own. WSL's `/mnt/d/...` and Windows' `D:\...` refer to the same files, so the VTU data does not need to be copied or simulated again.

### WSL / Linux

Use a virtual environment to keep plotting dependencies separate from system Python:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install numpy matplotlib
.venv/bin/python plot_demo_slice.py
```

If creating the environment reports missing `venv`/`ensurepip`, Ubuntu users can install `python3-venv` with `sudo apt install python3-venv` (run `sudo apt update` first if necessary), then retry. Repeating a pip command with a Python interpreter that has no pip will not resolve the problem.

If uv is already installed, an alternative is:

```bash
uv venv .venv
uv pip install --python .venv/bin/python numpy matplotlib
.venv/bin/python plot_demo_slice.py
```

The development machine also has an external plotting environment:

```bash
/mnt/d/2026-2027/computatianal/work/frjets_plot_env/bin/python plot_demo_slice.py
```

That environment is local to the development machine and is not included in the archive. Use your own interpreter when running the commands below: `python` on Windows or `.venv/bin/python` for the virtual environment above.

### Input, output, and slice selection

By default, the script reads `output/demo0002.vtu` and writes `demo_slice.png`. It only needs that one snapshot, with no dependency on smoke runs, MPI comparisons, or validation files.

To choose another snapshot, output filename, or target plane:

```bash
python3 plot_demo_slice.py output/demo0001.vtu --output demo_t05.png
python3 plot_demo_slice.py output/demo0002.vtu --z 0.5 --output demo_z05.png
```

The default is `--z 0`. The script selects the nearest existing cell-center plane; ties are resolved toward the positive side. For the supplied $64^3$ grid, this is normally $z=0.0625$. The title uses the actual snapshot time and selected coordinate. This is a cell-center slice, not interpolation to exactly $z=0$.

The reader supports the uncompressed `vtuBCCmpi` output generated by these configurations: LittleEndian data, UInt32 length headers, Float32 coordinates and fields, and Int32 connectivity for eight-vertex cells, including multiple Piece elements.

It requires a fixed, uniform Cartesian grid. For AMR or a different output format, use ParaView or adapt the reader. The script checks whether the data layout can be plotted; it does not assess the physical quality of the solution.

## Output variables and ParaView

Open `output/demo*.vtu` in ParaView and click Apply. To view a longitudinal section, add a Slice with normal `(0,0,1)` and origin `(0,0,0)`. Hide the original opaque Surface so it does not obscure the interior.

| Field | Meaning |
|---|---|
| `rho` | Density. |
| `v1`, `v2`, `v3` | Velocities along $x$, $y$, and $z$. The axial velocity is `v2`. |
| `p` | Pressure. |
| `trp1` | In the target AMRVAC version, this remains the tracer density $C=\rho f$. Use `trp1/rho` in Calculator to obtain the jet fraction $f$. |
| `T` | The dimensionless temperature proxy $p/\rho$, not a temperature directly expressed in kelvin. |

Initially, the interior contains a stationary atmosphere with zero tracer. Jet material subsequently enters from the center of the bottom boundary. The $t=1$ output represents early injection; this small-domain, short-time setup does not establish a developed FR I morphology.

The plotting script produces density, jet fraction, and axial-velocity panels. When comparing snapshots, note that separate images may use different automatically scaled color ranges.
