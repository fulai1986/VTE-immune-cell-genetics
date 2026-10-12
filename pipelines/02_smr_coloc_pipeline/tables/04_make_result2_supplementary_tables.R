#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

get_project_dir <- function() {
  env_dir <- Sys.getenv("PROJECT_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = FALSE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
    return(normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

project_dir <- get_project_dir()
supp_dir <- file.path(project_dir, "manuscript_results/tables/supplementary")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

cell_summary_file <- file.path(
  project_dir,
  "results/smr_analysis/qc/finngen_r13_discovery_onek1k_params.cellwise_bh_fdr_heidi001_summary.tsv"
)
tier_file <- file.path(
  project_dir,
  "results/smr_analysis/qc/finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf_beta_or.tsv"
)
all_signal_file <- file.path(
  project_dir,
  "manuscript_results/tables/supplementary/Table_S4_finngen_smr_coloc_all_fdr05_signals.tsv"
)

cell_order <- c(
  "CD4_NC", "CD8_ET", "NK", "CD8_NC", "B_IN", "CD4_ET", "B_Mem",
  "Mono_C", "CD8_S100B", "Mono_NC", "NK_R", "DC", "CD4_SOX4", "Plasma"
)
tier_order <- c("Tier 1", "Tier 2", "Tier 3", "Not prioritized")

cell_summary <- fread(cell_summary_file)
tier_dt <- fread(tier_file)

tier_dt[, `:=`(
  cell_label = gsub("_", " ", cell),
  tier = factor(tier, levels = tier_order)
)]

tier_counts <- dcast(
  tier_dt[, .N, by = .(cell, tier)],
  cell ~ tier,
  value.var = "N",
  fill = 0
)
for (col in tier_order) {
  if (!col %in% names(tier_counts)) {
    tier_counts[, (col) := 0L]
  }
}

setnames(
  cell_summary,
  old = c(
    "cell",
    "n_smr_rows",
    "n_cell_BH_FDR_lt_0.05",
    "n_cell_BH_FDR_lt_0.05_and_HEIDI_ge_0.01",
    "n_cell_BH_FDR_lt_0.05_and_HEIDI_ge_0.05",
    "min_p_SMR",
    "min_q_SMR_cell_BH",
    "top_gene",
    "top_snp",
    "top_p_SMR",
    "top_q_SMR_cell_BH",
    "top_p_HEIDI"
  ),
  new = c(
    "cell",
    "SMR-tested probes",
    "SMR FDR < 0.05 signals",
    "SMR FDR < 0.05 and HEIDI P >= 0.01 signals",
    "SMR FDR < 0.05 and HEIDI P >= 0.05 signals",
    "Minimum SMR P",
    "Minimum cell-wise SMR FDR",
    "Top gene",
    "Top SNP coordinate (hg19)",
    "Top SMR P",
    "Top cell-wise SMR FDR",
    "Top HEIDI P"
  )
)
cell_summary[, `Cell type` := gsub("_", " ", cell)]
cell_summary <- merge(cell_summary, tier_counts, by = "cell", all.x = TRUE)
for (col in tier_order) {
  cell_summary[is.na(get(col)), (col) := 0L]
}
cell_summary[, cell_order_idx := match(cell, cell_order)]

table_s4 <- cell_summary[
  order(
    -`SMR FDR < 0.05 signals`,
    -`Tier 1`,
    -`SMR FDR < 0.05 and HEIDI P >= 0.01 signals`,
    cell_order_idx
  ),
  .(
  `Cell type`,
  `SMR-tested probes`,
  `SMR FDR < 0.05 signals`,
  `SMR FDR < 0.05 and HEIDI P >= 0.01 signals`,
  `Tier 1 signals` = `Tier 1`,
  `Tier 2 signals` = `Tier 2`,
  `Tier 3 signals` = `Tier 3`,
  `Not prioritized signals` = `Not prioritized`,
  `Minimum SMR P`,
  `Minimum cell-wise SMR FDR`,
  `Top gene`,
  `Top SNP coordinate (hg19)`,
  `Top SMR P`,
  `Top cell-wise SMR FDR`,
  `Top HEIDI P`
)]

table_s4_file <- file.path(supp_dir, "Table_S4_cell_wise_smr_discovery_summary.tsv")
fwrite(table_s4, table_s4_file, sep = "\t")

table_s5 <- fread(all_signal_file)
table_s5[, `:=`(
  Tier = factor(Tier, levels = tier_order),
  tier_order_idx = as.integer(factor(Tier, levels = tier_order)),
  cell_order_idx = match(gsub(" ", "_", `Cell type`), cell_order)
)]
setorder(table_s5, tier_order_idx, Gene, cell_order_idx)
table_s5[, c("tier_order_idx", "cell_order_idx") := NULL]

table_s5_file <- file.path(supp_dir, "Table_S5_all_finngen_smr_fdr_significant_signals.tsv")
fwrite(table_s5, table_s5_file, sep = "\t")

readme <- file.path(supp_dir, "README_supplementary_tables.md")
readme_lines <- c(
  "# Supplementary Tables",
  "",
  "## Result 1: cell-type-specific eQTL",
  "",
  "- `Table_S1_cell_type_eqtl_analysis_summary.tsv`: Cell numbers, donor coverage, tested genes, and significant eGenes for each PBMC cell type.",
  "- `Table_S2_significant_egenes_by_cell_type.tsv.gz`: Significant eGenes defined by eigenMT-like gene-level correction followed by cell-type-wise BH-FDR < 0.05.",
  "- `Table_S3_onek1k_eqtl_concordance_benchmark.tsv`: Concordance between project pseudobulk eQTLs and OneK1K ROUND=1 primary eSNP-eGene benchmark sets.",
  "",
  "## Result 2: FinnGen SMR and colocalization",
  "",
  "- `Table_S4_cell_wise_smr_discovery_summary.tsv`: Cell-wise SMR discovery summary across 14 PBMC cell types, including numbers of FDR-significant, HEIDI-supported, and tiered signals.",
  "- `Table_S5_all_finngen_smr_fdr_significant_signals.tsv`: Full list of all 19 FinnGen SMR cell-wise FDR < 0.05 cell-gene pairs, including SMR, HEIDI, coloc.abf, rsID/hg19 coordinate, and meta-GWAS validation fields.",
  ""
)
writeLines(readme_lines, readme)

cat("Wrote ", table_s4_file, "\n", sep = "")
cat("Wrote ", table_s5_file, "\n", sep = "")
cat("Updated ", readme, "\n", sep = "")
