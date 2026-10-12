#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?Set PROJECT_DIR to the analysis workspace}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
RUNNER="${SCRIPT_DIR}/08_run_chr_full_eqtl.py"

CELL_TYPES=(
  "CD4 NC"
  "CD8 ET"
  "NK"
  "CD8 NC"
  "B IN"
  "CD4 ET"
  "B Mem"
  "CD8 S100B"
  "Mono C"
  "Mono NC"
  "NK R"
  "DC"
  "CD4 SOX4"
  "Plasma"
)

safe_name() {
  printf "%s" "$1" | sed -E 's/[^A-Za-z0-9]+/_/g; s/^_+//; s/_+$//'
}

cd "$PROJECT_DIR"

for cell_type in "${CELL_TYPES[@]}"; do
  cell_safe="$(safe_name "$cell_type")"
  for chrom in $(seq 1 22); do
    out="results/eqtl/${cell_safe}/chr${chrom}.full.cis_eqtl.tsv.gz"
    top="results/eqtl/${cell_safe}/chr${chrom}.full.top_by_gene.tsv"
    if [[ -s "$out" && -s "$top" ]]; then
      echo "[skip] ${cell_type} chr${chrom}"
      continue
    fi
    echo "[run] ${cell_type} chr${chrom}"
    "$PYTHON" -u "$RUNNER" --cell-type "$cell_type" --chrom "$chrom" --progress-every 20000
  done
done
