#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUTPUT_DIR="${1:-/projectnb/dietzelab/guYANG/Soil_Param/tsoil_meeting_figures_VER502}"
VERTICAL_POSITION="${2:-502}"

if [[ $# -ge 2 ]]; then
  shift 2
elif [[ $# -eq 1 ]]; then
  shift 1
fi

Rscript "${SCRIPT_DIR}/generate_tsoil_meeting_figures.R" \
  "--depth=${VERTICAL_POSITION}" \
  "--output-dir=${OUTPUT_DIR}" \
  "$@"

echo ""
echo "All figures and tables were written to: ${OUTPUT_DIR}"

