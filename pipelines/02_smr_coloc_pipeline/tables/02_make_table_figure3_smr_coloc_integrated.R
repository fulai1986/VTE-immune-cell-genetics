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
tier_file <- file.path(
  project_dir,
  "results/smr_analysis/qc/finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf_beta_or.tsv"
)
hit_file <- file.path(
  project_dir,
  "results/smr_analysis/qc/finngen_r13_discovery_onek1k_params.cellwise_bh_fdr_lt_0.05_hits.tsv"
)
rsid_file <- file.path(project_dir, "manuscript_results/tables/top_snp_rsid_mapping_hg19.tsv")
out_dir <- file.path(project_dir, "manuscript_results/tables")
supp_dir <- file.path(out_dir, "supplementary")

dt <- fread(tier_file)
hits <- fread(hit_file)[, .(cell, Gene, topSNP, se_SMR)]
dt <- merge(dt, hits, by = c("cell", "Gene", "topSNP"), all.x = TRUE)
rsid_map <- fread(rsid_file)[, .(topSNP, rsid, CHR_hg19, BP_hg19)]
dt <- merge(dt, rsid_map, by = "topSNP", all.x = TRUE)
if (anyNA(dt$se_SMR)) {
  stop("Missing se_SMR values after merging SMR hit file")
}
if (anyNA(dt$rsid)) {
  stop("Missing rsID values after merging rsID map")
}

tier_order <- c("Tier 1", "Tier 2", "Tier 3", "Not prioritized")
dt[, tier := factor(tier, levels = tier_order)]
dt[, cell_label := gsub("_", " ", cell)]
dt[, cell_gene_pair := paste0(cell_label, " - ", Gene)]
dt[, direction := fifelse(b_SMR_finngen > 0, "Risk", "Protective")]
dt[, `:=`(
  or_low = exp(b_SMR_finngen - 1.96 * se_SMR),
  or_high = exp(b_SMR_finngen + 1.96 * se_SMR),
  SMR_FDR = smr_q_SMR_cell_BH,
  HEIDI_P = smr_p_HEIDI,
  coloc_PP4 = PP.H4.abf,
  coloc_ratio = PP4_over_PP3_PP4
)]
dt[, OR_95CI := sprintf("%.2f (%.2f-%.2f)", OR_SMR_finngen, or_low, or_high)]

main <- dt[tier %in% c("Tier 1", "Tier 2", "Tier 3")]
main <- main[order(as.integer(tier), Gene, cell_label)]

main_out <- main[, .(
  Tier = as.character(tier),
  `Cell-gene pair` = cell_gene_pair,
  Gene,
  `Cell type` = cell_label,
  `Top SNP (rsID)` = rsid,
  `SMR beta (log OR)` = b_SMR_finngen,
  `SMR OR (95% CI)` = OR_95CI,
  Direction = direction,
  `SMR P` = smr_p_SMR,
  `Cell-wise SMR FDR` = SMR_FDR,
  `HEIDI P` = HEIDI_P,
  `coloc PP4` = coloc_PP4,
  `PP4/(PP3+PP4)` = coloc_ratio
)]

supp <- dt[order(as.integer(tier), Gene, cell_label)]
supp_out <- supp[, .(
  Tier = as.character(tier),
  Reason = reason,
  `Cell-gene pair` = cell_gene_pair,
  Gene,
  `Cell type` = cell_label,
  Chr = chr,
  `Top SNP (rsID)` = rsid,
  `Top SNP coordinate (hg19)` = topSNP,
  `Top SNP chr (hg19)` = CHR_hg19,
  `Top SNP position (hg19)` = BP_hg19,
  `FinnGen SMR beta (log OR)` = b_SMR_finngen,
  `FinnGen SMR SE` = se_SMR,
  `FinnGen SMR OR` = OR_SMR_finngen,
  `FinnGen SMR OR lower 95% CI` = or_low,
  `FinnGen SMR OR upper 95% CI` = or_high,
  Direction = direction,
  `GWAS P at top eQTL` = smr_p_GWAS,
  `eQTL P at top eQTL` = smr_p_eQTL,
  `SMR P` = smr_p_SMR,
  `Cell-wise SMR FDR` = SMR_FDR,
  `HEIDI P` = HEIDI_P,
  `coloc number of SNPs` = coloc_nsnps,
  `coloc PP3` = PP.H3.abf,
  `coloc PP4` = coloc_PP4,
  `PP4/(PP3+PP4)` = coloc_ratio,
  `Meta same direction` = meta_same_direction,
  `Meta SMR P < 0.05` = meta_nominal_p_SMR_lt_0.05,
  `Meta HEIDI P >= 0.01` = meta_heidi001_pass,
  `Meta HEIDI P >= 0.05` = meta_heidi005_pass,
  `Meta SMR beta` = b_SMR_meta,
  `Meta SMR OR` = OR_SMR_meta,
  `Meta SMR P` = meta_p_SMR,
  `Meta HEIDI P` = meta_p_HEIDI
)]

main_file <- file.path(out_dir, "Table_figure3_smr_coloc_prioritized_tier1_3.tsv")
supp_file <- file.path(supp_dir, "Table_S4_finngen_smr_coloc_all_fdr05_signals.tsv")
fwrite(main_out, main_file, sep = "\t")
fwrite(supp_out, supp_file, sep = "\t")

cat("Wrote ", main_file, "\n", sep = "")
cat("Wrote ", supp_file, "\n", sep = "")
