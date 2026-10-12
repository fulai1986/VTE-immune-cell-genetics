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
qc_dir <- file.path(project_dir, "results/smr_analysis/qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

smr_file <- file.path(
  qc_dir,
  "finngen_r13_discovery_onek1k_params.cellwise_bh_fdr_lt_0.05_hits.tsv"
)
coloc_file <- file.path(
  qc_dir,
  "finngen_r13_discovery_fdr05_hits.coloc_abf_summary.tsv"
)
meta_file <- file.path(
  qc_dir,
  "finngen_r13_discovery_onek1k_params.fdr05_hits.meta_validation.tsv"
)

out_file <- file.path(
  qc_dir,
  "finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf.tsv"
)
out_beta_or_file <- file.path(
  qc_dir,
  "finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf_beta_or.tsv"
)

smr <- fread(smr_file)
coloc <- fread(coloc_file)

keep_coloc <- c(
  "cell", "Gene", "chr", "topSNP", "coloc_status", "coloc_error", "coloc_nsnps",
  "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
  "PP4_over_PP3_PP4", "coloc_abf_positive_pp4_0.8"
)
keep_coloc <- intersect(keep_coloc, names(coloc))
dt <- merge(
  smr,
  coloc[, ..keep_coloc],
  by = c("cell", "Gene", "chr", "topSNP"),
  all.x = TRUE
)

dt[, coloc_positive := fifelse(
  coloc_status == "OK" &
    is.finite(PP.H4.abf) & PP.H4.abf >= 0.8 &
    is.finite(PP4_over_PP3_PP4) & PP4_over_PP3_PP4 >= 0.8,
  TRUE,
  FALSE
)]
dt[, heidi001_pass := is.finite(p_HEIDI) & p_HEIDI >= 0.01]

dt[, tier := fifelse(
  heidi001_pass & coloc_positive,
  "Tier 1",
  fifelse(
    !heidi001_pass & coloc_positive,
    "Tier 2",
    fifelse(heidi001_pass & !coloc_positive, "Tier 3", "Not prioritized")
  )
)]
dt[, reason := fifelse(
  tier == "Tier 1",
  "SMR FDR<0.05, HEIDI>=0.01, coloc.abf positive",
  fifelse(
    tier == "Tier 2",
    "SMR FDR<0.05, HEIDI<0.01, coloc.abf positive",
    fifelse(
      tier == "Tier 3",
      "SMR FDR<0.05, HEIDI>=0.01, coloc.abf not positive",
      "SMR FDR<0.05 without HEIDI/coloc support"
    )
  )
)]

dt[, `:=`(
  b_SMR_finngen = b_SMR,
  OR_SMR_finngen = exp(b_SMR),
  smr_p_GWAS = p_GWAS,
  smr_p_eQTL = p_eQTL,
  smr_p_SMR = p_SMR,
  smr_q_SMR_cell_BH = q_SMR_cell_BH,
  smr_p_HEIDI = p_HEIDI
)]

if (file.exists(meta_file)) {
  meta <- fread(meta_file)
  meta_keep <- c(
    "cell", "Gene",
    "same_direction", "meta_nominal_p_SMR_lt_0.05",
    "meta_heidi001_pass", "meta_heidi005_pass",
    "b_SMR_meta", "p_SMR_meta", "p_HEIDI_meta"
  )
  meta_keep <- intersect(meta_keep, names(meta))
  dt <- merge(dt, meta[, ..meta_keep], by = c("cell", "Gene"), all.x = TRUE)
  dt[, `:=`(
    OR_SMR_meta = exp(b_SMR_meta),
    meta_p_SMR = p_SMR_meta,
    meta_p_HEIDI = p_HEIDI_meta
  )]
} else {
  dt[, `:=`(
    same_direction = NA,
    meta_nominal_p_SMR_lt_0.05 = NA,
    meta_heidi001_pass = NA,
    meta_heidi005_pass = NA,
    b_SMR_meta = NA_real_,
    OR_SMR_meta = NA_real_,
    meta_p_SMR = NA_real_,
    meta_p_HEIDI = NA_real_
  )]
}

tier_levels <- c("Tier 1", "Tier 2", "Tier 3", "Not prioritized")
dt[, tier := factor(tier, levels = tier_levels)]
setorder(dt, tier, cell, Gene)

out_cols <- c(
  "tier", "reason", "cell", "Gene", "chr", "topSNP",
  "b_SMR_finngen", "OR_SMR_finngen",
  "smr_p_GWAS", "smr_p_eQTL", "smr_p_SMR", "smr_q_SMR_cell_BH", "smr_p_HEIDI",
  "coloc_nsnps", "PP.H3.abf", "PP.H4.abf", "PP4_over_PP3_PP4",
  "same_direction", "meta_nominal_p_SMR_lt_0.05",
  "meta_heidi001_pass", "meta_heidi005_pass",
  "b_SMR_meta", "OR_SMR_meta", "meta_p_SMR", "meta_p_HEIDI"
)
out_cols <- intersect(out_cols, names(dt))

fwrite(dt[, ..out_cols], out_file, sep = "\t")
fwrite(dt[, ..out_cols], out_beta_or_file, sep = "\t")

message("Wrote ", out_file)
message("Wrote ", out_beta_or_file)
