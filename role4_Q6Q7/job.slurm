#!/usr/bin/env bash
#SBATCH --clusters=wice
#SBATCH --account=lp_edu_amrvac
#SBATCH --ntasks=4
#SBATCH --time=04:00:00
#SBATCH --job-name=role4
# Do not sbatch this directly: use submit_all.sh, which sets RUN, AMRVAC_EXE,
# --chdir=runs/<RUN> and a per-run --output file.
set -uo pipefail
: "${RUN:?RUN not set}" "${AMRVAC_EXE:?AMRVAC_EXE not set}"

module purge
module load foss/2025a
export OMP_NUM_THREADS=1

mkdir -p output
{
  echo "run=$RUN  jobid=${SLURM_JOB_ID:-}  ntasks=${SLURM_NTASKS:-}  nodes=${SLURM_JOB_NODELIST:-}"
  echo "exe=$AMRVAC_EXE  amrvac_commit=$(git -C "${AMRVAC_DIR:-/nonexistent}" log -1 --format=%h 2>/dev/null || echo unknown)"
  echo "start=$(date -Is)"
} > timing.txt

t0=$(date +%s)
mpirun -np "${SLURM_NTASKS}" "$AMRVAC_EXE" -i "${RUN}.par" > "${RUN}.stdout" 2>&1
status=$?
t1=$(date +%s)

{
  echo "end=$(date -Is)"
  echo "exit_status=$status"
  echo "wall_seconds=$((t1-t0))"
} >> timing.txt
cat timing.txt
exit $status
