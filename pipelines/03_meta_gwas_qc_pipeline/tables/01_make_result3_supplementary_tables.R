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
    return(normalizePath(file.path(script_dir, "..", "..", "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

project_dir <- get_project_dir()
out_dir <- file.path(
  project_dir,
  "manuscript_results/completed_results/result3_meta_gwas_qc/supplementary"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(file.path(project_dir, "gwas/config/gwas_dataset_manifest.tsv"))
std <- fread(file.path(project_dir, "gwas/qc/discovery_gwas_standardization_summary.tsv"))
metal <- fread(file.path(project_dir, "gwas/qc/vte_discovery_finngen_independent_ivw_ukb_approx_logor.metal_input_prepare_summary.tsv"))

meta_ids <- c("decode_chb_cvdc_dbds", "decode_iceland", "decode_intermountain", "ukb_i9_vte")
cohort_labels <- c(
  decode_chb_cvdc_dbds = "deCODE DK",
  decode_iceland = "deCODE Iceland",
  decode_intermountain = "deCODE Intermountain",
  ukb_i9_vte = "UK Biobank"
)
source_labels <- c(
  decode_chb_cvdc_dbds = "deCODE public summary statistics",
  decode_iceland = "deCODE public summary statistics",
  decode_intermountain = "deCODE public summary statistics",
  ukb_i9_vte = "OpenGWAS / UK Biobank"
)
metal_ids <- c(
  decode_chb_cvdc_dbds = "decode_chb_cvdc_dbds",
  decode_iceland = "decode_iceland",
  decode_intermountain = "decode_intermountain",
  ukb_i9_vte = "ukb_i9_vte_approx_logor"
)

tab_s6 <- manifest[dataset_id %in% meta_ids]
tab_s6[, Dataset := cohort_labels[dataset_id]]
tab_s6[, Source := source_labels[dataset_id]]
tab_s6[, Role := "FinnGen-independent meta-GWAS input"]
tab_s6[, Trait := trait]
tab_s6[, Ancestry := ancestry]
tab_s6[, Cases := cases]
tab_s6[, Controls := controls]
tab_s6[, `Total N` := cases + controls]
tab_s6[, `Genome build` := "GRCh37/hg19"]
tab_s6 <- merge(
  tab_s6,
  std[, .(dataset_id, `Raw variants` = raw_rows, `Standardized variants` = written_rows)],
  by = "dataset_id",
  all.x = TRUE
)
metal_join <- copy(metal)
metal_join[, dataset_id := names(metal_ids)[match(dataset_id, metal_ids)]]
tab_s6 <- merge(
  tab_s6,
  metal_join[, .(dataset_id, `METAL input variants` = metal_rows, `Duplicate markers removed` = duplicate_marker_rows_removed)],
  by = "dataset_id",
  all.x = TRUE
)
tab_s6[, Notes := fifelse(
  dataset_id == "ukb_i9_vte",
  "UK Biobank VTE summary statistics included in FinnGen-independent meta-GWAS; FinnGen was not included.",
  notes
)]
tab_s6 <- tab_s6[match(meta_ids, dataset_id)]
tab_s6 <- tab_s6[, .(
  Dataset, Source, Role, Trait, Ancestry, Cases, Controls, `Total N`,
  `Genome build`, `Raw variants`, `Standardized variants`,
  `METAL input variants`, `Duplicate markers removed`, Notes
)]
fwrite(tab_s6, file.path(out_dir, "Table_S6_meta_gwas_cohort_information.tsv"), sep = "\t", quote = FALSE, na = "")

tab_s7 <- data.table(
  Metric = c(
    "Primary meta-GWAS variants",
    "OneK1K LD reference variants",
    "OneK1K LD-matched downstream variants",
    "OneK1K LD-matched retained fraction",
    "Genome-wide significant variants",
    "Genome-wide significant variants with HetP > 0.01",
    "Genome-wide significant retained fraction after heterogeneity filter"
  ),
  Value = c("6,891,026", "7,248,782", "5,783,976", "83.9%", "3,428", "2,626", "76.6%"),
  Notes = c(
    "Variants available in at least three contributing cohorts (n cohorts >= 3).",
    "OneK1K PLINK LD reference; duplicate positions = 2,806.",
    "Primary meta-GWAS variants harmonized to the OneK1K LD reference for downstream SMR/coloc validation.",
    "5,783,976 / 6,891,026 variants retained.",
    "Variants with P < 5e-8 in the primary meta-GWAS; includes two P-value underflow records plotted at the y-axis cap.",
    "Genome-wide significant variants passing the predefined heterogeneity screen.",
    "2,626 / 3,428 genome-wide significant variants retained after HetP > 0.01."
  )
)
fwrite(tab_s7, file.path(out_dir, "Table_S7_meta_gwas_variant_coverage_and_qc_summary.tsv"), sep = "\t", quote = FALSE, na = "")

message("Wrote Tables S6-S7 to ", out_dir)
