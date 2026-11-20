#!/usr/bin/env bash
# Run this on YOUR LAPTOP, inside your clone of comp-method-team.
# Copies the 7 Role-4 .par files + run7_ster.slurm to the IvS cluster.
#
#   bash role4_Q6Q7/push_to_ster.sh r1081254@<ster-host>
#   bash role4_Q6Q7/push_to_ster.sh r1081254@<ster-host> --submit   (also sbatch it)
#
# It takes the files straight from the GitHub branch, so your own checkout
# and working files are not changed.
set -euo pipefail
REMOTE=${1:?usage: push_to_ster.sh user@host [--submit]}
SUBMIT=${2:-}
BRANCH=claude/vibrant-carson-bwcfes
REMOTE_DIR=/STER/agungp/amrvac_runs/frjets_caseB_cartesian_3d_final
RUNS="baseline96 woodward96 vanleer96 koren96 cada3_96 hll96 hllc96"

cd "$(git rev-parse --show-toplevel)"
git fetch -q origin "$BRANCH"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
for r in $RUNS; do git show "FETCH_HEAD:role4_Q6Q7/pars/$r.par" > "$TMP/$r.par"; done
git show FETCH_HEAD:role4_Q6Q7/run7_ster.slurm > "$TMP/run7.slurm"
echo "Copying to $REMOTE:$REMOTE_DIR :"; ls "$TMP"

# Protect old outputs that the new baseline96 would overwrite, then copy.
ssh "$REMOTE" "cd $REMOTE_DIR && test -x amrvac || echo 'WARNING: no ./amrvac in $REMOTE_DIR (compile it first)';
  if ls output/baseline96* >/dev/null 2>&1; then
    d=output/old_baseline96_\$(date +%Y%m%d_%H%M%S); mkdir -p \$d && mv output/baseline96* \$d/ && echo \"moved old output/baseline96* to \$d\"; fi"
scp -q "$TMP"/*.par "$TMP/run7.slurm" "$REMOTE:$REMOTE_DIR/"
echo "Done."

if [ "$SUBMIT" = "--submit" ]; then
  ssh "$REMOTE" "cd $REMOTE_DIR && sbatch run7.slurm && squeue -u \$USER"
else
  echo "Next, on the cluster:  cd $REMOTE_DIR && sbatch run7.slurm"
fi
