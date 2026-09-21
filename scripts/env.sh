# Usage: source scripts/env.sh
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export EMBERSIM_ROOT="$(dirname "$_here")"

# Toolchain used for the thesis runs (NCI Gadi). Edit for other systems.
if command -v module >/dev/null 2>&1; then
  module load intel-mpi/2021.17.2
fi

export NEK_SOURCE_ROOT="$EMBERSIM_ROOT/external/Nek5000"
export TOOLBOX_SRC="$EMBERSIM_ROOT/external/KTH_Toolbox"
