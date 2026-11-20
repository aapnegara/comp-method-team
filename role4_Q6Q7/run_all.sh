#!/usr/bin/env bash
# Run all 7 Role-4 schemes one after another with the SAME core count.
# Usage:  NP=4 ./run_all.sh            (all runs)
#         NP=4 ./run_all.sh Q6_koren   (selected runs)
# A failing run does not stop the others; see output/<name>.stdout.
set -u
NP=${NP:-4}
mkdir -p output
runs=("$@"); [ ${#runs[@]} -eq 0 ] && runs=(Q6_woodward Q6_vanleer Q6_koren Q6_cada3 Q7_hll Q7_hllc Q7_roe_tvd)
echo "run,np,status,wall_s" > output/timing.csv
for r in "${runs[@]}"; do
  echo ">>> $r  (np=$NP)"
  t0=$(date +%s)
  mpirun -np "$NP" ./amrvac -i common.par "schemes/$r.par" > "output/$r.stdout" 2>&1
  st=$?; t1=$(date +%s)
  echo "$r,$NP,$st,$((t1-t0))" >> output/timing.csv
  [ $st -eq 0 ] && echo "    done in $((t1-t0)) s" || echo "    FAILED (exit $st) after $((t1-t0)) s -> output/$r.stdout"
done
