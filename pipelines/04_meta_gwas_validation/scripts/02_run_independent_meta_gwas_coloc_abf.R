#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
})

project_dir <- Sys.getenv("PROJECT_DIR")

candidate_path <- file.path(
  project_dir,
  "results/smr_analysis/qc/finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_validation.cellwise_fdr.tsv"
)
gwas_path <- file.path(
  project_dir,
  "gwas/harmonized/oneK1K_ld/vte_discovery_finngen_independent_ivw_ukb_approx_logor.primary_ncohort3plus.onek1k_ld_harmonized.smr.ma"
)
eqtl_qc_root <- file.path(project_dir, "results/qc/eqtl")
qfile_root <- file.path(project_dir, "results/smr")
out_root <- file.path(project_dir, "results/coloc/meta_ncohort3plus_ukb_approx_logor_tier1_3_abf")
qc_dir <- file.path(project_dir, "results/smr_analysis/qc")
report_dir <- file.path(project_dir, "reports/smr_analysis")
display_dir <- file.path(project_dir, "manuscript_results/completed_results/result2_smr_coloc/source_tables")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(display_dir, recursive = TRUE, showWarnings = FALSE)

discovery_cases <- 18569 + 11700 + 6981 + 4620
discovery_controls <- 213503 + 311745 + 51252 + 356574
discovery_s <- discovery_cases / (discovery_cases + discovery_controls)

safe_id <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
to_num <- function(x) suppressWarnings(as.numeric(x))
maf <- function(freq) pmin(freq, 1 - freq)
truthy <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes")

read_eqtl_n <- function(cell, chrom) {
  summary_path <- file.path(eqtl_qc_root, cell, sprintf("chr%d.full.summary.tsv", chrom))
  summary <- fread(summary_path)
  as.numeric(summary[metric == "donors", value][1])
}

read_gene_qfile <- function(cell, chrom, gene) {
  qfile <- file.path(
    qfile_root,
    cell,
    sprintf("chr%d", chrom),
    sprintf("%s.chr%d.r2_0.8.maf_0.01.allele_id.qfile.txt", cell, chrom)
  )
  if (!file.exists(qfile)) {
    stop("Missing qfile: ", qfile)
  }
  q <- fread(qfile)
  q[Gene == gene]
}

read_gwas_chr <- function(chrom) {
  g <- fread(gwas_path, showProgress = FALSE)
  g[startsWith(SNP, paste0(chrom, ":"))][
    ,
    BP := as.integer(tstrsplit(SNP, ":", fixed = TRUE)[[2]])
  ]
}

candidates <- fread(candidate_path)[tier %in% c("Tier 1", "Tier 2", "Tier 3")]
candidates[, candidate_id := sprintf("%02d_%s_%s_chr%s", .I, cell, Gene, chr)]

gwas_cache <- new.env(parent = emptyenv())
summary_rows <- vector("list", nrow(candidates))

for (i in seq_len(nrow(candidates))) {
  candidate <- candidates[i]
  cell <- candidate$cell
  chrom <- as.integer(candidate$chr)
  gene <- candidate$Gene
  candidate_id <- safe_id(candidate$candidate_id)
  message(sprintf("[%d/%d] %s %s chr%d", i, nrow(candidates), cell, gene, chrom))

  chrom_key <- as.character(chrom)
  if (!exists(chrom_key, envir = gwas_cache, inherits = FALSE)) {
    assign(chrom_key, read_gwas_chr(chrom), envir = gwas_cache)
  }

  gwas_chr <- get(chrom_key, envir = gwas_cache, inherits = FALSE)
  eqtl <- read_gene_qfile(cell, chrom, gene)
  eqtl_n <- read_eqtl_n(cell, chrom)

  status <- "OK"
  error <- ""
  coloc_summary <- c(
    nsnps = NA_real_,
    PP.H0.abf = NA_real_,
    PP.H1.abf = NA_real_,
    PP.H2.abf = NA_real_,
    PP.H3.abf = NA_real_,
    PP.H4.abf = NA_real_
  )
  merged <- data.table()

  tryCatch({
    if (!nrow(eqtl)) {
      stop("No eQTL rows for gene")
    }

    eqtl_keep <- eqtl[
      ,
      .(
        SNP,
        BP_eqtl = as.integer(BP),
        A1_eqtl = A1,
        A2_eqtl = A2,
        freq_eqtl = to_num(Freq),
        beta_eqtl = to_num(b),
        se_eqtl = to_num(se),
        p_eqtl = to_num(p)
      )
    ]
    gwas_keep <- gwas_chr[
      ,
      .(
        SNP,
        BP_gwas = as.integer(BP),
        A1_gwas = A1,
        A2_gwas = A2,
        freq_gwas = to_num(freq),
        beta_gwas = to_num(b),
        se_gwas = to_num(se),
        p_gwas = to_num(p),
        n_gwas = to_num(n)
      )
    ]

    merged <- merge(eqtl_keep, gwas_keep, by = "SNP", all = FALSE)
    merged <- merged[
      is.finite(beta_eqtl) &
        is.finite(se_eqtl) & se_eqtl > 0 &
        is.finite(beta_gwas) &
        is.finite(se_gwas) & se_gwas > 0 &
        is.finite(freq_eqtl) & freq_eqtl > 0 & freq_eqtl < 1 &
        is.finite(freq_gwas) & freq_gwas > 0 & freq_gwas < 1
    ]
    merged <- unique(merged, by = "SNP")
    if (nrow(merged) < 50) {
      stop(sprintf("Too few overlapping SNPs: %d", nrow(merged)))
    }

    dataset_eqtl <- list(
      beta = merged$beta_eqtl,
      varbeta = merged$se_eqtl^2,
      snp = merged$SNP,
      position = merged$BP_eqtl,
      MAF = maf(merged$freq_eqtl),
      N = eqtl_n,
      type = "quant"
    )
    dataset_gwas <- list(
      beta = merged$beta_gwas,
      varbeta = merged$se_gwas^2,
      snp = merged$SNP,
      position = merged$BP_gwas,
      MAF = maf(merged$freq_gwas),
      N = round(stats::median(merged$n_gwas, na.rm = TRUE)),
      s = discovery_s,
      type = "cc"
    )
    coloc_result <- coloc.abf(dataset1 = dataset_gwas, dataset2 = dataset_eqtl)
    coloc_summary <- coloc_result$summary
  }, error = function(e) {
    status <<- "FAIL"
    error <<- conditionMessage(e)
  })

  candidate_dir <- file.path(out_root, candidate_id)
  dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
  if (nrow(merged)) {
    fwrite(merged, file.path(candidate_dir, "coloc_input_merged.tsv"), sep = "\t")
  }

  pp3 <- unname(coloc_summary["PP.H3.abf"])
  pp4 <- unname(coloc_summary["PP.H4.abf"])
  ratio <- ifelse(is.finite(pp4) && is.finite(pp3) && (pp3 + pp4) > 0, pp4 / (pp3 + pp4), NA_real_)

  summary_rows[[i]] <- data.table(
    candidate_id = candidate_id,
    tier = candidate$tier,
    cell = cell,
    Gene = gene,
    chr = chrom,
    topSNP_finngen = candidate$topSNP,
    b_SMR_finngen = to_num(candidate$b_SMR_finngen),
    OR_SMR_finngen = to_num(candidate$OR_SMR_finngen),
    finngen_q_SMR_cell_BH = to_num(candidate$smr_q_SMR_cell_BH),
    finngen_p_HEIDI = to_num(candidate$smr_p_HEIDI),
    finngen_coloc_PP4 = to_num(candidate$coloc_abf_PP4),
    meta_topSNP = candidate$meta_ukb_approx_topSNP,
    meta_b_SMR = to_num(candidate$meta_ukb_approx_b_SMR),
    meta_OR = to_num(candidate$meta_ukb_approx_OR),
    meta_OR_95CI = candidate$meta_ukb_approx_OR_95CI,
    meta_p_SMR = to_num(candidate$meta_ukb_approx_p_SMR),
    meta_q_SMR_cell_BH = to_num(candidate$meta_ukb_approx_q_SMR_cell_BH),
    meta_p_HEIDI = to_num(candidate$meta_ukb_approx_p_HEIDI),
    meta_same_direction = truthy(candidate$meta_ukb_approx_same_direction),
    meta_validation_pass_cellwise_fdr = truthy(candidate$meta_ukb_approx_validation_pass_cellwise_fdr),
    eqtl_n = eqtl_n,
    gwas_n_median_in_coloc = if (nrow(merged)) round(stats::median(merged$n_gwas, na.rm = TRUE)) else NA_real_,
    coloc_status = status,
    coloc_error = error,
    coloc_nsnps = unname(coloc_summary["nsnps"]),
    PP.H0.abf = unname(coloc_summary["PP.H0.abf"]),
    PP.H1.abf = unname(coloc_summary["PP.H1.abf"]),
    PP.H2.abf = unname(coloc_summary["PP.H2.abf"]),
    PP.H3.abf = pp3,
    PP.H4.abf = pp4,
    PP4_over_PP3_PP4 = ratio
  )
}

summary <- rbindlist(summary_rows, fill = TRUE)
summary[, coloc_abf_positive := coloc_status == "OK" & PP.H4.abf >= 0.8 & PP4_over_PP3_PP4 >= 0.8]
summary[, validation_and_coloc_positive := meta_validation_pass_cellwise_fdr == TRUE & coloc_abf_positive == TRUE]
setorder(summary, tier, -validation_and_coloc_positive, meta_p_SMR)

summary_path <- file.path(qc_dir, "finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_coloc_abf_summary.tsv")
positive_path <- file.path(qc_dir, "finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_coloc_abf_positive.tsv")
display_path <- file.path(display_dir, "tier1_3_meta_ncohort3plus_ukb_approx_logor_coloc_abf_display.tsv")
report_path <- file.path(report_dir, "FINNGEN_R13_DISCOVERY_TIER1_3_META_NCOHORT3PLUS_UKB_APPROX_LOGOR_COLOC_ABF_SUMMARY.md")

fwrite(summary, summary_path, sep = "\t", na = "NA")
fwrite(summary[coloc_abf_positive == TRUE], positive_path, sep = "\t", na = "NA")

display <- summary[
  ,
  .(
    Tier = tier,
    `Cell-gene pair` = paste(gsub("_", " ", cell), Gene, sep = " - "),
    `Meta OR (95% CI)` = meta_OR_95CI,
    `Meta p_SMR` = meta_p_SMR,
    `Meta cell-wise FDR` = meta_q_SMR_cell_BH,
    `Meta p_HEIDI` = meta_p_HEIDI,
    `Meta SMR validation` = ifelse(meta_validation_pass_cellwise_fdr, "Pass", "Not pass"),
    `Meta coloc nsnps` = coloc_nsnps,
    `Meta coloc PP4` = PP.H4.abf,
    `PP4/(PP3+PP4)` = PP4_over_PP3_PP4,
    `Meta coloc positive` = ifelse(coloc_abf_positive, "Yes", "No"),
    `Validation + coloc` = ifelse(validation_and_coloc_positive, "Yes", "No")
  )
]
fwrite(display, display_path, sep = "\t", na = "NA")

counts_tier <- summary[
  ,
  .(
    n = .N,
    meta_smr_validation = sum(meta_validation_pass_cellwise_fdr, na.rm = TRUE),
    coloc_abf_positive = sum(coloc_abf_positive, na.rm = TRUE),
    validation_and_coloc_positive = sum(validation_and_coloc_positive, na.rm = TRUE)
  ),
  by = tier
][order(tier)]

report_lines <- c(
  "# Meta-GWAS coloc.abf for FinnGen Tier 1-3 candidates: UKB approximate log(OR)",
  "",
  "This run re-tested colocalization in the FinnGen-independent meta-GWAS validation dataset generated after converting UKB `ukb-d-I9_VTE` from signed Z scores to approximate log(OR).",
  "",
  "Inputs:",
  "",
  "- GWAS: `gwas/harmonized/oneK1K_ld/vte_discovery_finngen_independent_ivw_ukb_approx_logor.primary_ncohort3plus.onek1k_ld_harmonized.smr.ma`.",
  "- eQTL: matching OneK1K pseudobulk qfile rows for each cell-gene candidate.",
  "- Candidate set: FinnGen discovery Tier 1-3 candidates.",
  "",
  "Settings:",
  "",
  "- Method: `coloc.abf` from coloc.",
  sprintf("- GWAS case fraction: %d / (%d + %d) = %.5f.", discovery_cases, discovery_cases, discovery_controls, discovery_s),
  "- GWAS N: median SNP-level N within each coloc region, from the harmonized meta-GWAS file.",
  "- Positive rule: `PP.H4.abf >= 0.8` and `PP4/(PP3+PP4) >= 0.8`.",
  "",
  "## Counts",
  "",
  sprintf("- Candidates tested: %d", nrow(summary)),
  sprintf("- Successful coloc.abf runs: %d", sum(summary$coloc_status == "OK")),
  sprintf("- Meta coloc.abf positive: %d", sum(summary$coloc_abf_positive, na.rm = TRUE)),
  sprintf("- Meta SMR validation and coloc.abf positive: %d", sum(summary$validation_and_coloc_positive, na.rm = TRUE)),
  "",
  "| Tier | N | Meta SMR validation | Meta coloc.abf positive | Both |",
  "|---|---:|---:|---:|---:|",
  sprintf(
    "| %s | %d | %d | %d | %d |",
    counts_tier$tier,
    counts_tier$n,
    counts_tier$meta_smr_validation,
    counts_tier$coloc_abf_positive,
    counts_tier$validation_and_coloc_positive
  ),
  "",
  "## Outputs",
  "",
  sprintf("- Summary: `%s`", sub(paste0(project_dir, "/"), "", summary_path)),
  sprintf("- Positives: `%s`", sub(paste0(project_dir, "/"), "", positive_path)),
  sprintf("- Display table: `%s`", sub(paste0(project_dir, "/"), "", display_path)),
  sprintf("- Per-candidate merged inputs: `%s`", sub(paste0(project_dir, "/"), "", out_root))
)
writeLines(report_lines, report_path)

message("Wrote ", summary_path)
message("Wrote ", positive_path)
message("Wrote ", display_path)
message("Wrote ", report_path)
