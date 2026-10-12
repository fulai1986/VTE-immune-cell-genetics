#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?Set PROJECT_DIR to the analysis workspace}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/02_make_celltype_pseudobulk_expression.R"
LOG_DIR="$PROJECT_DIR/scratch/logs"
mkdir -p "$LOG_DIR"

CELLTYPES=(
  "CD4 NC"
  "CD4 ET"
  "CD4 SOX4"
  "CD8 ET"
  "CD8 NC"
  "CD8 S100B"
  "NK"
  "NK R"
  "B IN"
  "B Mem"
  "Plasma"
  "Mono C"
  "Mono NC"
  "DC"
)

for celltype in "${CELLTYPES[@]}"; do
  safe_name="$(echo "$celltype" | sed 's/[^A-Za-z0-9]/_/g')"
  summary="$PROJECT_DIR/results/qc/${safe_name}.all_pools.pseudobulk_summary.tsv"
  if [[ -s "$summary" ]]; then
    echo "[skip] $celltype already has $summary"
    continue
  fi
  echo "[run] $celltype"
  Rscript "$SCRIPT" --cell-type "$celltype" --pools all --min-cells-per-donor 10 \
    2>&1 | tee "$LOG_DIR/pseudobulk_${safe_name}.all_pools.log"
done

echo "[done] all requested pseudobulk jobs"
