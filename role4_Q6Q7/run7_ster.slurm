#!/bin/bash
#SBATCH -p ivs-short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J role4
# Role 4 (Q6/Q7) on the IvS cluster: the 7 runs one after another, 4 MPI tasks each.
# Lives in (and runs from) the folder that holds amrvac + the 7 .par files.
# If you needed 'module load ...' to compile amrvac here, add the same lines below.

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p output
echo "host=$(hostname)  start=$(date)"
for r in baseline96 woodward96 vanleer96 koren96 cada3_96 hll96 hllc96; do
  s=$(date +%s)
  mpirun -np 4 ./amrvac -i $r.par > $r.stdout 2>&1
  echo "$r exit=$? wall=$(( $(date +%s) - s ))s"
done
echo "end=$(date)"
