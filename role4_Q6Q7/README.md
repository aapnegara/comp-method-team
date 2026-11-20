# Role 4: influence of limiters (Q6) and flux solvers (Q7)

Controlled numerical comparison on the 3D Cartesian Case-B jet
(`frjets_caseB_cartesian_3d_final/`, Zoey's setup). **Only the numerical scheme
changes between runs**; physics, grid, boundaries and CFL are identical.

## The runs

| Run | time_stepper | flux_scheme | limiter | Experiment |
|---|---|---|---|---|
| `baseline96` | threestep | tvdlf | minmod | A + B (shared reference) |
| `woodward96` | threestep | tvdlf | woodward | A: limiters (Q6) |
| `vanleer96` | threestep | tvdlf | vanleer | A |
| `koren96` | threestep | tvdlf | koren | A |
| `cada3_96` | threestep | tvdlf | cada3 | A |
| `hll96` | threestep | hll | minmod | B: solvers (Q7) |
| `hllc96` | threestep | hllc | minmod | B |
| `roe_tvd96` *(optional)* | **onestep** | **tvd** | minmod | Q7 "afterwards" |

Shared settings (copied verbatim from the 96³ run that already worked on wICE,
kept in `reference/smoke_wice_baseline96.par`): 96³ cells, blocks of 8³,
x,z ∈ [-4,4], y ∈ [0,8], `refine_max_level=1` (no AMR), CFL 0.3 (`maxsum`),
`nghostcells=2`, `small_values_method='error'`, γ = 5/3, one tracer.
Changed for all runs: `time_max=5`, log every 0.05, snapshot (.dat + .vtu) every 0.25.

## Compatibility checks (against the AMRVAC source, `src/io/mod_input_output.t`, `src/limiter/mod_limiter.t`, `src/hd/`)

- `woodward`, `vanleer`, `koren`, `cada3` are valid limiter names. AMRVAC raises
  `nghostcells` only for mp5 / ppm / weno-type limiters, so **2 ghost cells are
  sufficient for every run here** (and `mod_usr.t` accepts 2 or 3).
- `hll` and `hllc` are available for HD (`mod_hd_hllc.t`); neither changes the stencil.
- `tvd` (Roe-type TVD): AMRVAC **stops** unless `time_stepper='onestep'`, and it
  **forces `dimsplit=.true.`**. So `roe_tvd96` changes three things at once, not
  just the solver. Also, the HD Roe eigen-system (`mod_hd_roe.t`) has no explicit
  tracer wave, so tracer behaviour in that run must be checked before trusting it.
  This is why it is optional and kept outside the 7 controlled runs.
- `mod_usr.t` is an unchanged copy of `frjets_caseB_cartesian_3d_final/mod_usr.t`.

## Files

| File | Purpose |
|---|---|
| `reference/smoke_wice_baseline96.par` | The `.par` that already ran successfully on wICE (untouched). |
| `gen_pars.sh` | Regenerates `pars/*.par` from the reference and prints the validation diff. |
| `pars/*.par` | The 7 (+1) parameter files. |
| `mod_usr.t` | Identical copy of the team's user module. |
| `build.sh` | Compiles one `amrvac` executable for all runs into `build/`. |
| `job.slurm` | Generic SLURM job (wICE, `lp_edu_amrvac`, `foss/2025a`, 4 tasks). |
| `submit_all.sh` | Submits one job per run, each in `runs/<name>/`; `--test`, `--dry-run`, `--force`. |
| `analyze.py` | Diagnostics, CSV tables and comparison figures into `analysis/`. |

## How to run on wICE

```bash
# 0. Get this folder onto the cluster (source on DATA, runs on SCRATCH)
cd $VSC_DATA/massaglia-jets/role4 && git pull       # or copy role4_Q6Q7/ there
cp -r role4_Q6Q7 $VSC_SCRATCH/massaglia-jets/role4/
cd $VSC_SCRATCH/massaglia-jets/role4/role4_Q6Q7

# 1. Build once (light compile, fine on the login node)
./build.sh
#    ...or reuse the binary that already worked (mod_usr.t is identical):
#    export AMRVAC_EXE=$VSC_SCRATCH/massaglia-jets/role4/test/frjets_caseB_cartesian_3d_final/amrvac

# 2. Validate the parameter files: baseline96 may differ from the reference
#    ONLY in time_max and dtsave; every other run ONLY in its scheme line.
./gen_pars.sh
diff pars/baseline96.par pars/hllc96.par

# 3. Dry run: shows the sbatch commands without submitting
./submit_all.sh --dry-run

# 4. Short test: baseline to t = 0.5 (30 min limit)
./submit_all.sh --test
squeue -M wice -u $USER
```

**Check the test before going further:**

1. `runs/test_baseline96/timing.txt` shows `exit_status=0`.
2. `tail runs/test_baseline96/output/test_baseline96.log` reaches t = 0.5.
3. `python analyze.py --ymax 3` → `analysis/slices_A_limiters.png`: inflow at
   y = 0 around x = 0; v_y close to 4 in the nozzle; f close to 1 in jet material;
   no NaNs and ρ > 0 in `summary.csv`.
4. Use `wall_seconds` to estimate the cost of the full runs: t = 5 costs roughly
   10× the t = 0.5 test. Adjust `WALLTIME` if needed (default 4 h).

Then submit the main runs:

```bash
./submit_all.sh                       # the 7 controlled runs
./submit_all.sh roe_tvd96             # optional, after the 7 are done
NTASKS=8 WALLTIME=08:00:00 ./submit_all.sh koren96   # example override
PARTITION=batch_sapphirerapids ./submit_all.sh       # default is batch_icelake
```

**If a run reaches its walltime:** `job.slurm` passes an extra `walltime.par`
with `wall_time_max` = SLURM walltime − 5 min. AMRVAC then leaves its time loop
by itself and writes the current state as a final snapshot, so nothing is lost.
`timing.txt` shows `t_reached` < 5. Continue from the latest `.dat` with:

```bash
WALLTIME=01:00:00 ./submit_all.sh --resume hllc96
```

`analyze.py` adds up `wall_seconds` over all parts of a resumed run.

Each run lives in `runs/<name>/` (`<name>.par`, `output/`, `<name>.stdout`,
`timing.txt`, `slurm-<id>.out`). A non-empty `output/` is never overwritten
unless you pass `--force`. Your earlier smoke/baseline outputs in `test/` are
not touched. **Keep the same `NTASKS` for all 7 runs**, so wall times can be compared.

## Analysis

```bash
module load matplotlib   # or a venv with numpy + matplotlib
python analyze.py --ymax 4
```

| Output | Content |
|---|---|
| `summary.csv` | Per run: scheme, exit status, failure message, t reached, number of steps, wall time, cell updates/s/core, final y_head, y_bow, ρ and p min/max, V_jet, V_mix. |
| `timeseries.csv` | The same diagnostics for every snapshot. |
| `yhead_<exp>.png` | y_head(t) and y_bow(t). |
| `axis_<exp>.png` | On-axis ρ, v_y, f and p at a common time. |
| `slices_<exp>.png` | z = 0 slices of log ρ, f, v_y and log p, with the same colour scale per column. |

Definitions are identical for every run (see the docstring in `analyze.py`):
f = trp1/ρ; y_head = largest y with f > 0.5; y_bow = largest y with
|p/0.6 − 1| > 0.05; V_mix = volume with 0.1 < f < 0.9.

V_mix measures the **numerical** smearing of the jet/ambient interface at this
resolution. It is not physical turbulent mixing.

Rough expectation: the 1D ram-pressure balance gives
v_head ≈ v_j √η / (1 + √η) ≈ 0.36, so y_head(5) ≈ 2 r_j. The jet stays well
inside the domain, and the y_max boundary should not matter by t = 5.

## What to look for (discussion for Q6/Q7)

- **Limiters (A):** minmod is the most diffusive. vanleer, woodward and koren/cada3
  should give a thinner transition layer (smaller V_mix), sharper shocks and
  more small-scale structure in the cocoon. They may also show oscillations or
  smaller ρ_min / p_min, and a higher risk of `small value` errors.
- **Solvers (B):** tvdlf → hll → hllc means progressively less dissipation.
  HLLC resolves the contact wave, so the jet/ambient density jump should become
  visibly sharper. Compare the cost per step against the gain in sharpness.
- A crash is a result: report the time and location from `<name>.stdout`, and do
  not switch off `small_values_method` silently.
- Runs at a single resolution show differences between schemes, not convergence.
  Link this to Yitske's resolution study (Q5) for the order-of-accuracy discussion.
