#!/usr/bin/env bash
# Compile ONE amrvac executable for all Role-4 runs (mod_usr.t is identical for all).
# Light compile; OK on a login node. Usage: ./build.sh
set -euo pipefail
cd "$(dirname "$0")"
: "${AMRVAC_DIR:?AMRVAC_DIR is not set}"
module load foss/2025a 2>/dev/null || true
mkdir -p build && cp mod_usr.t build/
cd build
"$AMRVAC_DIR/setup.pl" -d=3 -v=3
touch amrvac.h
make -j4
echo "Built: $PWD/amrvac   (AMRVAC commit: $(git -C "$AMRVAC_DIR" log -1 --format=%h 2>/dev/null || echo unknown))"
