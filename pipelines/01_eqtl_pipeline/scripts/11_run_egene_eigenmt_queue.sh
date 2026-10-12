#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?Set PROJECT_DIR to the analysis workspace}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/10_run_egene_eigenmt_cell_chr_plink.R"
SUMMARY_SCRIPT="${SCRIPT_DIR}/12_summarize_egene_eigenmt_cell.R"
OUT_DIR="${PROJECT_DIR}/manuscript_results/tables/egene_eigenmt"
LOG_DIR="${PROJECT_DIR}/manuscript_results/logs/egene_eigenmt"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"

CELL_TYPES="${CELL_TYPES:-CD4 NC;CD4 ET;CD4 SOX4;CD8 NC;CD8 ET;CD8 S100B;NK;NK R;B IN;B Mem;Mono C;Mono NC;DC;Plasma}"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
GENE_BATCH_SIZE="${GENE_BATCH_SIZE:-10000}"
EIGENMT_BLOCK_SIZE="${EIGENMT_BLOCK_SIZE:-200}"
EIGENMT_VAR_THRESH="${EIGENMT_VAR_THRESH:-0.99}"

safe_name() {
  echo "$1" | sed -E 's/[^A-Za-z0-9]+/_/g'
}

IFS=';' read -r -a CELLS <<< "${CELL_TYPES}"
for CELL in "${CELLS[@]}"; do
  CELL_SAFE="$(safe_name "${CELL}")"
  echo "[$(date)] Starting cell: ${CELL}"
  for CHR in ${CHROMS}; do
    OUT_FILE="${OUT_DIR}/${CELL_SAFE}.chr${CHR}.eigenmt_like.by_gene.tsv"
    LOG_FILE="${LOG_DIR}/${CELL_SAFE}.chr${CHR}.eigenmt_like.log"
    if [[ -s "${OUT_FILE}" ]]; then
      echo "[$(date)] Skip existing ${CELL} chr${CHR}: ${OUT_FILE}"
      continue
    fi
    echo "[$(date)] Run ${CELL} chr${CHR}"
    Rscript "${SCRIPT}" \
      --cell-type "${CELL}" \
      --chrom "${CHR}" \
      --gene-batch-size "${GENE_BATCH_SIZE}" \
      --eigenmt-block-size "${EIGENMT_BLOCK_SIZE}" \
      --eigenmt-var-thresh "${EIGENMT_VAR_THRESH}" \
      > "${LOG_FILE}" 2>&1
  done
  Rscript "${SUMMARY_SCRIPT}" --cell-type "${CELL}" \
    > "${LOG_DIR}/${CELL_SAFE}.cellwise_summary.log" 2>&1
  echo "[$(date)] Finished cell: ${CELL}"
done
