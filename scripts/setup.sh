#!/bin/bash
# One-time setup: fetch pinned dependencies and build Nek5000 tools.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
git submodule update --init external/Nek5000 external/KTH_Toolbox external/ppiclF
source scripts/env.sh
( cd external/Nek5000/tools && ./maketools genmap gmsh2nek )
