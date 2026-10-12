#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- Sys.getenv("PROJECT_DIR")
if (!nzchar(project_dir)) {
  project_dir <- Sys.getenv("PROJECT_DIR")
}
project_dir <- normalizePath(project_dir, mustWork = TRUE)

result2_table <- file.path(
  project_dir,
  "manuscript_results/tables/Table_figure3_smr_coloc_prioritized_tier1_3.tsv"
)
meta_smr_table <- file.path(
  project_dir,
  "results/smr_analysis/qc/finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_validation.cellwise_fdr.tsv"
)
meta_coloc_table <- file.path(
  project_dir,
  "results/smr_analysis/qc/finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_coloc_abf_summary.tsv"
)
rsid_table <- file.path(
  project_dir,
  "manuscript_results/tables/top_snp_rsid_mapping_hg19.tsv"
)

out_dir <- file.path(
  project_dir,
  "manuscript_results/completed_results/result4_meta_gwas_validation/supplementary"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fg <- fread(result2_table)
meta <- fread(meta_smr_table)
coloc <- fread(meta_coloc_table)
rsid <- fread(rsid_table)

fg[, pair_key := gsub(" ", "_", `Cell-gene pair`)]
fg[, pair_key := sub("_-_", " - ", pair_key)]
fg[, pair_key := gsub(" - ", "_-_", pair_key)]
fg[, pair_key := paste0(gsub(" ", "_", `Cell type`), "_-_", Gene)]

meta[, pair_key := paste0(cell, "_-_", Gene)]
coloc[, pair_key := paste0(cell, "_-_", Gene)]

setnames(rsid, c("topSNP", "rsid"), c("snp_hg19", "rsid"))

meta_cols <- c(
  "pair_key",
  "topSNP",
  "meta_ukb_approx_topSNP",
  "meta_ukb_approx_b_SMR",
  "meta_ukb_approx_se_SMR",
  "meta_ukb_approx_OR",
  "meta_ukb_approx_OR_95CI",
  "meta_ukb_approx_p_SMR",
  "meta_ukb_approx_q_SMR_cell_BH",
  "meta_ukb_approx_p_HEIDI",
  "meta_ukb_approx_same_direction",
  "meta_ukb_approx_validation_pass_cellwise_fdr"
)

coloc_cols <- c(
  "pair_key",
  "coloc_nsnps",
  "PP.H3.abf",
  "PP.H4.abf",
  "PP4_over_PP3_PP4",
  "coloc_abf_positive",
  "validation_and_coloc_positive"
)

tab <- merge(
  fg,
  meta[, ..meta_cols],
  by = "pair_key",
  all.x = TRUE,
  sort = FALSE
)
tab <- merge(
  tab,
  coloc[, ..coloc_cols],
  by = "pair_key",
  all.x = TRUE,
  sort = FALSE
)

tab <- merge(
  tab,
  rsid[, .(snp_hg19, finngen_top_snp_rsid = rsid)],
  by.x = "topSNP",
  by.y = "snp_hg19",
  all.x = TRUE,
  sort = FALSE
)
tab <- merge(
  tab,
  rsid[, .(snp_hg19, meta_top_snp_rsid = rsid)],
  by.x = "meta_ukb_approx_topSNP",
  by.y = "snp_hg19",
  all.x = TRUE,
  sort = FALSE
)

tier_order <- c("Tier 1", "Tier 2", "Tier 3")
pair_order <- c(
  "CD8 ET - FADS1",
  "NK - FADS1",
  "NK - MAFF",
  "CD8 ET - SLC44A2",
  "CD8 ET - TMEM258",
  "CD8 NC - TMEM258",
  "CD4 NC - TSPAN15",
  "CD4 ET - FADS1",
  "NK - SLC44A2"
)
tab[, `Cell-gene pair` := factor(`Cell-gene pair`, levels = pair_order)]
setorder(tab, `Cell-gene pair`)
tab[, `Cell-gene pair` := as.character(`Cell-gene pair`)]

fmt_num <- function(x, digits = 6) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "fg", digits = digits))
}
fmt_p <- function(x) {
  fifelse(
    is.na(x),
    NA_character_,
    fifelse(x < 0.001, formatC(x, format = "e", digits = 3), sprintf("%.4f", x))
  )
}
pass_fail <- function(x) {
  fifelse(is.na(x), NA_character_, fifelse(x, "Pass", "Fail"))
}

out <- tab[, .(
  Tier,
  `Cell-gene pair`,
  Gene,
  `Cell type`,
  `FinnGen top SNP (rsID)` = finngen_top_snp_rsid,
  `FinnGen top SNP (hg19)` = topSNP,
  `FinnGen SMR beta (log OR)` = fmt_num(`SMR beta (log OR)`),
  `FinnGen OR (95% CI)` = `SMR OR (95% CI)`,
  `FinnGen SMR P` = fmt_p(`SMR P`),
  `FinnGen cell-wise FDR` = fmt_p(`Cell-wise SMR FDR`),
  `FinnGen HEIDI P` = fmt_p(`HEIDI P`),
  `FinnGen HEIDI status` = fifelse(`HEIDI P` >= 0.01, "Pass", "Fail"),
  `FinnGen coloc PP4` = fmt_num(`coloc PP4`),
  `Meta-GWAS top SNP (rsID)` = meta_top_snp_rsid,
  `Meta-GWAS top SNP (hg19)` = meta_ukb_approx_topSNP,
  `Meta-GWAS SMR beta (log OR)` = fmt_num(meta_ukb_approx_b_SMR),
  `Meta-GWAS OR (95% CI)` = meta_ukb_approx_OR_95CI,
  `Meta-GWAS SMR P` = fmt_p(meta_ukb_approx_p_SMR),
  `Meta-GWAS cell-wise FDR` = fmt_p(meta_ukb_approx_q_SMR_cell_BH),
  `Meta-GWAS HEIDI P` = fmt_p(meta_ukb_approx_p_HEIDI),
  `Meta-GWAS HEIDI status` = fifelse(meta_ukb_approx_p_HEIDI >= 0.01, "Pass", "Fail"),
  `Meta-GWAS coloc SNPs` = coloc_nsnps,
  `Meta-GWAS coloc PP3` = fmt_num(PP.H3.abf),
  `Meta-GWAS coloc PP4` = fmt_num(PP.H4.abf),
  `Meta-GWAS PP4/(PP3+PP4)` = fmt_num(PP4_over_PP3_PP4),
  `Direction consistent` = pass_fail(meta_ukb_approx_same_direction),
  `Meta-GWAS SMR-HEIDI validation` = pass_fail(meta_ukb_approx_validation_pass_cellwise_fdr),
  `Meta-GWAS coloc positive` = pass_fail(coloc_abf_positive),
  `Overall replication status` = pass_fail(validation_and_coloc_positive)
)]

out_path <- file.path(out_dir, "Table_S9_independent_meta_gwas_validation_of_finngen_prioritized_smr_coloc_signals.tsv")
fwrite(out, out_path, sep = "\t")

message("Wrote ", out_path)
