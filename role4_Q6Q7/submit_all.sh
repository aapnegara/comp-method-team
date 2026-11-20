#!/usr/bin/env bash
# Submit Role-4 runs to wICE, one job per run, each in its own directory:
#   runs/<name>/{<name>.par, output/, <name>.stdout, timing.txt, slurm-<jobid>.out}
#
# Usage:
#   ./submit_all.sh --test                 short baseline check (t=0.5, 30 min)
#   ./submit_all.sh                        the 7 main runs
#   ./submit_all.sh hllc96 koren96         selected runs
#   ./submit_all.sh roe_tvd96              optional Roe run
#   ./submit_all.sh --resume hllc96        continue a run that stopped at its walltime
# Options:  --dry-run (print, do not submit)   --force (start again from t=0 in a used run dir)
#           --resume (continue from the latest snapshot; keeps existing output)
# Env:      PARTITION=batch_icelake  NTASKS=4  WALLTIME=04:00:00  AMRVAC_EXE=<path>  (default: ./build/amrvac)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
PARTITION=${PARTITION:-batch_icelake}   # keep ALL runs on one partition (fair cost comparison)
NTASKS=${NTASKS:-4}
WALLTIME=${WALLTIME:-04:00:00}
AMRVAC_EXE=${AMRVAC_EXE:-$ROOT/build/amrvac}
DRY=0; FORCE=0; TEST=0; RESUME=0; runs=()
for a in "$@"; do
  case $a in
    --dry-run) DRY=1 ;; --force) FORCE=1 ;; --test) TEST=1 ;; --resume) RESUME=1 ;;
    *) runs+=("$a") ;;
  esac
done
[ ${#runs[@]} -eq 0 ] && runs=(baseline96 woodward96 vanleer96 koren96 cada3_96 hll96 hllc96)
[ $TEST -eq 1 ] && runs=(test_baseline96)
[ -x "$AMRVAC_EXE" ] || { echo "ERROR: executable not found: $AMRVAC_EXE (run ./build.sh or set AMRVAC_EXE)"; exit 1; }

for name in "${runs[@]}"; do
  dir=$ROOT/runs/$name
  if [ $RESUME -eq 1 ]; then
    ls "$dir"/output/"$name"[0-9][0-9][0-9][0-9].dat >/dev/null 2>&1 || { echo "SKIP $name: nothing to resume"; continue; }
  elif [ -d "$dir/output" ] && [ -n "$(ls -A "$dir/output" 2>/dev/null)" ] && [ $FORCE -eq 0 ]; then
    echo "SKIP $name: $dir/output is not empty (use --force, or move it away first)"; continue
  fi
  mkdir -p "$dir/output"
  wt=$WALLTIME
  if [ $RESUME -eq 1 ]; then
    :   # keep the existing <name>.par
  elif [ "$name" = test_baseline96 ]; then
    sed -e "s|^ base_filename=.*| base_filename='output/test_baseline96'|" \
        -e "s|^ time_max=.*| time_max=0.5d0|" \
        -e "s|^ dtsave(2)=.*| dtsave(2)=0.1d0|" "$ROOT/pars/baseline96.par" > "$dir/$name.par"
    wt=00:30:00
  else
    [ -f "$ROOT/pars/$name.par" ] || { echo "ERROR: no pars/$name.par"; exit 1; }
    cp "$ROOT/pars/$name.par" "$dir/$name.par"
  fi
  wt_min=$(awk -F: '{ if (NF==3) print $1*60+$2; else if (NF==2) print $1; else print $1 }' <<< "${wt#*-}")
  cmd=(sbatch --partition="$PARTITION" --ntasks="$NTASKS" --time="$wt" --job-name="r4_$name" --chdir="$dir"
       --output="$dir/slurm-%j.out" --export=ALL,RUN="$name",AMRVAC_EXE="$AMRVAC_EXE",WALLTIME_MIN="$wt_min",RESUME="$RESUME"
       "$ROOT/job.slurm")
  if [ $DRY -eq 1 ]; then echo "[dry-run] ${cmd[*]}"; else echo "submit $name"; "${cmd[@]}"; fi
done
echo "Monitor with:  squeue -M wice -u \$USER"
