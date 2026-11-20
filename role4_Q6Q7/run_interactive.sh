#!/usr/bin/env bash
# Run Role-4 runs ONE AFTER ANOTHER inside an interactive allocation, e.g.
#   srun -M wice -A lp_edu_amrvac -p batch_icelake --nodes=1 --ntasks=4 -t 01:45:00 --pty bash -l
#   bash run_interactive.sh                  (the 7 main runs)
#   bash run_interactive.sh hllc96 koren96   (selected runs)
# Uses the same job.slurm as the batch jobs, so outputs/timing.txt are identical.
# Env: PER_RUN_MIN=25 (AMRVAC stops cleanly and saves after PER_RUN_MIN-5 min)
set -uo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
[ -n "${SLURM_JOB_ID:-}" ] || { echo "ERROR: not inside a SLURM allocation (never run this on a login node)"; exit 1; }
AMRVAC_EXE=${AMRVAC_EXE:-$ROOT/build/amrvac}
[ -x "$AMRVAC_EXE" ] || { echo "ERROR: no executable $AMRVAC_EXE (run ./build.sh)"; exit 1; }
NT=${SLURM_NTASKS:-4}
runs=("$@"); [ ${#runs[@]} -eq 0 ] && runs=(baseline96 woodward96 vanleer96 koren96 cada3_96 hll96 hllc96)

for name in "${runs[@]}"; do
  dir=$ROOT/runs/$name
  if [ -n "$(ls -A "$dir/output" 2>/dev/null)" ]; then
    echo "SKIP $name: $dir/output is not empty"; continue
  fi
  mkdir -p "$dir/output" && cp "$ROOT/pars/$name.par" "$dir/$name.par"
  echo ">>> $(date +%H:%M:%S)  $name  (ntasks=$NT)"
  ( cd "$dir" && RUN=$name AMRVAC_EXE=$AMRVAC_EXE SLURM_NTASKS=$NT \
      WALLTIME_MIN=${PER_RUN_MIN:-25} RESUME=0 bash "$ROOT/job.slurm" ) \
    | grep -E "exit_status|wall_seconds|t_reached" | sed 's/^/    /'
done
echo "All done. Check:  grep -H -E 'exit_status|wall_seconds|t_reached' runs/*/timing.txt"
