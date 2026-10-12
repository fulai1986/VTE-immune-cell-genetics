#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
tab_dir <- file.path(project_dir, "manuscript_results/tables")
egene_dir <- file.path(tab_dir, "egene_eigenmt")
out_dir <- file.path(tab_dir, "supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_order <- c(
  "CD4 NC", "CD8 ET", "NK", "CD8 NC", "B IN", "CD4 ET", "B Mem",
  "Mono C", "CD8 S100B", "Mono NC", "NK R", "DC", "CD4 SOX4", "Plasma"
)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

summary_file <- file.path(tab_dir, "figure2_eqtl_egene_summary_for_plot.tsv")
summary_dt <- fread(summary_file)
summary_dt[, cell_type := factor(cell_type, levels = cell_order)]
setorder(summary_dt, cell_type)

supp_s1 <- summary_dt[, .(
  cell_type = as.character(cell_type),
  cells_retained,
  donors_contributing_to_pseudobulk = donors_retained,
  genes_with_expression = genes_expression,
  genes_considered_for_cis_eqtl = genes_considered,
  cis_pairs_tested,
  chromosomes_completed,
  genes_tested_for_egene = genes_tested_eigenmt,
  significant_egenes = egene_eigenmt_fdr05,
  egene_proportion = egene_fraction,
  median_meff_fraction,
  lead_p_lt_1e4,
  lead_p_lt_1e6,
  lead_p_lt_1e10,
  median_neglog10_lead_p = median_neglog10p,
  q90_neglog10_lead_p = q90_neglog10p,
  max_neglog10_lead_p = max_neglog10p,
  egene_definition = "P_gene = min(P_lead x M_eff, 1); cell-type-wise BH-FDR < 0.05"
)]
fwrite(
  supp_s1,
  file.path(out_dir, "Table_S1_cell_type_eqtl_analysis_summary.tsv"),
  sep = "\t"
)

egene_files <- file.path(egene_dir, paste0(safe_name(cell_order), ".eigenmt_like.cellwise.by_gene.tsv.gz"))
if (!all(file.exists(egene_files))) {
  stop("Missing eGene files: ", paste(egene_files[!file.exists(egene_files)], collapse = ", "))
}

supp_s2 <- rbindlist(lapply(egene_files, function(f) {
  fread(f)[egene_eigenmt_fdr05 == TRUE, .(
    cell_type,
    chrom,
    gene_id,
    gene_symbol,
    gene_start,
    gene_end,
    lead_snp = top_snp,
    lead_snp_pos = top_snp_pos,
    lead_beta = top_beta,
    lead_se = top_se,
    lead_p = p_lead,
    n_cis_snps = n_cis_snps_qfile,
    n_ld_snps,
    meff = meff_eigenmt_like,
    meff_fraction = meff_eigenmt_fraction,
    p_gene_bonferroni = p_gene_bonferroni_eigenmt,
    p_gene_sidak = p_gene_sidak_eigenmt,
    cell_wise_fdr_bonferroni = q_eigenmt_bonferroni_cell,
    cell_wise_fdr_sidak = q_eigenmt_sidak_cell
  )]
}), use.names = TRUE)
supp_s2[, cell_type := factor(cell_type, levels = cell_order)]
setorder(supp_s2, cell_type, cell_wise_fdr_bonferroni, lead_p, gene_symbol)
supp_s2[, cell_type := as.character(cell_type)]
fwrite(
  supp_s2,
  file.path(out_dir, "Table_S2_significant_egenes_by_cell_type.tsv.gz"),
  sep = "\t"
)

concordance_file <- file.path(tab_dir, "figure2c_official_top_random500_concordance_summary.tsv")
concordance <- fread(concordance_file)
concordance[, cell_type := factor(cell_type, levels = cell_order)]
concordance[, benchmark_set := fifelse(
  sample_set == "official_top500_by_p",
  "strongest_500_ROUND1_primary_eSNP_eGene_pairs",
  "random_500_ROUND1_primary_eSNP_eGene_pairs"
)]
setorder(concordance, cell_type, benchmark_set)
supp_s3 <- concordance[, .(
  cell_type = as.character(cell_type),
  benchmark_set,
  pairs_sampled = official_pairs_sampled,
  pairs_matched_in_pseudobulk_eqtl = pairs_found_in_ours,
  matched_pair_fraction = pair_found_rate,
  direction_comparable_pairs = direction_comparable,
  direction_concordant_pairs = direction_concordant,
  direction_concordance = direction_concordance_rate,
  nominal_support_p_lt_0_05 = our_nominal_p05_among_found,
  median_pseudobulk_p_among_matched_pairs = median_our_p_among_found,
  spearman_rho_neglog10p = spearman_neglog10p_corr_found,
  pearson_r_neglog10p = pearson_neglog10p_corr_found
)]
fwrite(
  supp_s3,
  file.path(out_dir, "Table_S3_onek1k_eqtl_concordance_benchmark.tsv"),
  sep = "\t"
)

readme <- c(
  "# Supplementary tables for Figure 2 eQTL results",
  "",
  "Table S1. Cell-type eQTL analysis summary.",
  "One row per PBMC cell type. Significant eGenes are defined as P_gene = min(P_lead x M_eff, 1), followed by cell-type-wise BH-FDR < 0.05.",
  "",
  "Table S2. Significant eGenes by cell type.",
  "One row per significant cell-type-gene pair. Lead SNP and beta correspond to the pseudobulk cis-eQTL lead association for that gene in that cell type.",
  "",
  "Table S3. OneK1K eQTL concordance benchmark.",
  "Concordance between our pseudobulk eQTL results and OneK1K ROUND=1 primary eSNP-eGene pairs, including strongest 500 and random 500 benchmarks per cell type."
)
writeLines(readme, file.path(out_dir, "README_supplementary_tables.md"))

cat("Wrote supplementary tables to", out_dir, "\n")
cat("S1 rows:", nrow(supp_s1), "\n")
cat("S2 rows:", nrow(supp_s2), "\n")
cat("S3 rows:", nrow(supp_s3), "\n")
