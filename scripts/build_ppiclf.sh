#!/bin/bash
# Run from inside a case directory containing SIZE and ppiclf/.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
CASEDIR="$(pwd)"
OUT="$CASEDIR/ppiclf_build"
for f in SIZE ppiclf/PPICLF_USER.h ppiclf/ppiclf_user.f; do
  [ -f "$CASEDIR/$f" ] || { echo "Missing $f: run from a case directory"; exit 1; }
done
[ -f "$EMBERSIM_ROOT/external/ppiclF/Makefile" ] || { echo "Run scripts/setup.sh first"; exit 1; }

rm -rf "$OUT" && mkdir -p "$OUT"
cp -r "$EMBERSIM_ROOT/external/ppiclF/." "$OUT/"
rm -rf "$OUT/.git"
patch -p1 -d "$OUT" < "$EMBERSIM_ROOT/patches/ppiclf-post-rk3.patch"
cp ppiclf/PPICLF_USER.h ppiclf/ppiclf_user.f SIZE "$OUT/source/"

cd "$OUT"
make FFLAGS="-cpp -std=legacy -mcmodel=medium -g -O" CPFLAGS="-E -mcmodel=medium -g -O"
find "$OUT" -name 'libppiclF.a'
