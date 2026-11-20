#!/usr/bin/env bash
#SBATCH --clusters=wice
#SBATCH --account=lp_edu_amrvac
#SBATCH --ntasks=4
#SBATCH --time=04:00:00
#SBATCH --job-name=role4
# Do not sbatch this directly: use submit_all.sh, which sets RUN, AMRVAC_EXE,
# WALLTIME_MIN, RESUME, --chdir=runs/<RUN> and a per-run --output file.
#
# Safe stop: AMRVAC gets wall_time_max = (SLURM walltime - 5 min), so it leaves
# the time loop by itself and WRITES the current state as a final snapshot
# before SLURM would kill it. Continue such a run later with RESUME=1
# (submit_all.sh --resume <name>): AMRVAC '-resume' restarts from the latest .dat.
set -uo pipefail
: "${RUN:?RUN not set}" "${AMRVAC_EXE:?AMRVAC_EXE not set}"

module purge
module load foss/2025a
export OMP_NUM_THREADS=1

mkdir -p output

# Extra par file read AFTER <RUN>.par: only sets the wall-clock stop (no physics/numerics).
WALLTIME_MIN=${WALLTIME_MIN:-240}
STOP_H=$(awk -v m="$WALLTIME_MIN" 'BEGIN{h=(m-5)/60; if(h<0.05)h=0.05; printf "%.3f", h}')
printf '&stoplist\n wall_time_max=%sd0\n/\n' "$STOP_H" > walltime.par

EXTRA=()
if [ "${RESUME:-0}" = 1 ]; then
  EXTRA=(-resume)
  [ -f "output/${RUN}.log" ] && cp "output/${RUN}.log" "output/${RUN}.log.before_resume_$(date +%s)"
  echo "----- resume $(date -Is)" >> timing.txt
else
  : > timing.txt
fi
{
  echo "run=$RUN  jobid=${SLURM_JOB_ID:-}  ntasks=${SLURM_NTASKS:-}  nodes=${SLURM_JOB_NODELIST:-}"
  echo "exe=$AMRVAC_EXE  amrvac_commit=$(git -C "${AMRVAC_DIR:-/nonexistent}" log -1 --format=%h 2>/dev/null || echo unknown)"
  echo "start=$(date -Is)"
  echo "wall_time_max_hours=$STOP_H  resume=${RESUME:-0}"
} >> timing.txt

t0=$(date +%s)
mpirun -np "${SLURM_NTASKS}" "$AMRVAC_EXE" -i "${RUN}.par" walltime.par "${EXTRA[@]}" >> "${RUN}.stdout" 2>&1
status=$?
t1=$(date +%s)

{
  echo "end=$(date -Is)"
  echo "exit_status=$status"
  echo "wall_seconds=$((t1-t0))"
  echo "t_reached=$(awk 'NF>2 && $1 ~ /^[0-9]+$/ {t=$2} END{print t}' "output/${RUN}.log" 2>/dev/null)"
} >> timing.txt
cat timing.txt
exit $status
