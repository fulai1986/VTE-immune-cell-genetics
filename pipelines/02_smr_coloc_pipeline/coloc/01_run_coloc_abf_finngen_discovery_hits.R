#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(coloc)
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

candidate_path <- file.path(project_dir, "results/smr_analysis/qc/finngen_r13_discovery_onek1k_params.cellwise_bh_fdr_lt_0.05_hits.tsv")
gwas_path <- file.path(project_dir, "gwas/replication/finngen_r13/harmonized_oneK1K_ld/finngen_R13_I9_VTE.hg19.onek1k_ld_harmonized.smr.ma")
eqtl_qc_root <- file.path(project_dir, "results/qc/eqtl")
qfile_root <- file.path(project_dir, "results/smr")
out_root <- file.path(project_dir, "results/coloc/finngen_r13_discovery_abf")
qc_dir <- file.path(project_dir, "results/smr_analysis/qc")
report_dir <- file.path(project_dir, "reports/smr_analysis")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, ...) {
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

safe_id <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

read_eqtl_n <- function(cell, chrom) {
  summary_path <- file.path(eqtl_qc_root, cell, sprintf("chr%d.full.summary.tsv", chrom))
  summary <- read_tsv(summary_path)
  as.numeric(summary$value[summary$metric == "donors"][1])
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
  con <- file(qfile, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
  gene_i <- match("Gene", header)
  rows <- list()
  n <- 0L
  repeat {
    lines <- readLines(con, n = 200000)
    if (!length(lines)) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    keep <- vapply(parts, function(x) length(x) >= gene_i && identical(x[[gene_i]], gene), logical(1))
    if (any(keep)) {
      for (x in parts[keep]) {
        n <- n + 1L
        rows[[n]] <- x
      }
    }
  }
  if (!length(rows)) {
    return(data.frame())
  }
  matrix_rows <- do.call(rbind, rows)
  df <- as.data.frame(matrix_rows, stringsAsFactors = FALSE)
  names(df) <- header
  df
}

read_gwas_chr <- function(chrom) {
  con <- file(gwas_path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(trimws(readLines(con, n = 1)), "\\s+")[[1]]
  snp_i <- match("SNP", header)
  rows <- list()
  n <- 0L
  repeat {
    lines <- readLines(con, n = 300000)
    if (!length(lines)) break
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    parts <- strsplit(lines, "\\s+")
    keep <- vapply(
      parts,
      function(x) length(x) >= snp_i && startsWith(x[[snp_i]], paste0(chrom, ":")),
      logical(1)
    )
    if (any(keep)) {
      for (x in parts[keep]) {
        n <- n + 1L
        rows[[n]] <- x
      }
    }
  }
  if (!length(rows)) {
    return(data.frame())
  }
  matrix_rows <- do.call(rbind, rows)
  df <- as.data.frame(matrix_rows, stringsAsFactors = FALSE)
  names(df) <- header
  df$BP <- vapply(strsplit(df$SNP, ":", fixed = TRUE), function(x) x[[2]], character(1))
  df
}

to_num <- function(x) suppressWarnings(as.numeric(x))
maf <- function(freq) pmin(freq, 1 - freq)

candidates <- read_tsv(candidate_path)
candidates$candidate_id <- sprintf("%02d_%s_%s_chr%s", seq_len(nrow(candidates)), candidates$cell, candidates$Gene, candidates$chr)

gwas_cache <- new.env(parent = emptyenv())
summary_rows <- list()

for (i in seq_len(nrow(candidates))) {
  candidate <- candidates[i, ]
  cell <- candidate$cell
  chrom <- as.integer(candidate$chr)
  gene <- candidate$Gene
  candidate_id <- safe_id(candidate$candidate_id)
  message(sprintf("[%d/%d] %s %s chr%d", i, nrow(candidates), cell, gene, chrom))

  if (!exists(as.character(chrom), envir = gwas_cache, inherits = FALSE)) {
    assign(as.character(chrom), read_gwas_chr(chrom), envir = gwas_cache)
  }
  gwas_chr <- get(as.character(chrom), envir = gwas_cache, inherits = FALSE)
  eqtl <- read_gene_qfile(cell, chrom, gene)
  eqtl_n <- read_eqtl_n(cell, chrom)

  status <- "OK"
  error <- ""
  coloc_summary <- rep(NA_real_, 6)
  names(coloc_summary) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
  merged <- data.frame()

  tryCatch({
    if (!nrow(eqtl)) {
      stop("No eQTL rows for gene")
    }
    eqtl_keep <- eqtl[, c("SNP", "BP", "A1", "A2", "Freq", "b", "se", "p")]
    names(eqtl_keep) <- c("SNP", "BP_eqtl", "A1_eqtl", "A2_eqtl", "freq_eqtl", "beta_eqtl", "se_eqtl", "p_eqtl")
    gwas_keep <- gwas_chr[, c("SNP", "BP", "A1", "A2", "freq", "b", "se", "p", "n")]
    names(gwas_keep) <- c("SNP", "BP_gwas", "A1_gwas", "A2_gwas", "freq_gwas", "beta_gwas", "se_gwas", "p_gwas", "n_gwas")
    merged <- merge(eqtl_keep, gwas_keep, by = "SNP", all = FALSE)
    if (!nrow(merged)) {
      stop("No overlapping SNPs")
    }
    num_cols <- c("BP_eqtl", "freq_eqtl", "beta_eqtl", "se_eqtl", "p_eqtl", "BP_gwas", "freq_gwas", "beta_gwas", "se_gwas", "p_gwas", "n_gwas")
    for (col in num_cols) {
      merged[[col]] <- to_num(merged[[col]])
    }
    merged <- merged[
      is.finite(merged$beta_eqtl) &
        is.finite(merged$se_eqtl) & merged$se_eqtl > 0 &
        is.finite(merged$beta_gwas) &
        is.finite(merged$se_gwas) & merged$se_gwas > 0 &
        is.finite(merged$freq_eqtl) & merged$freq_eqtl > 0 & merged$freq_eqtl < 1 &
        is.finite(merged$freq_gwas) & merged$freq_gwas > 0 & merged$freq_gwas < 1,
    ]
    merged <- merged[!duplicated(merged$SNP), ]
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
      N = 500186,
      s = 27928 / 500186,
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
    write.table(
      merged,
      file.path(candidate_dir, "coloc_input_merged.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  row <- data.frame(
    candidate_id = candidate_id,
    cell = cell,
    Gene = gene,
    chr = chrom,
    probeID = candidate$probeID,
    Probe_bp = candidate$Probe_bp,
    topSNP = candidate$topSNP,
    smr_p_GWAS = candidate$p_GWAS,
    smr_p_eQTL = candidate$p_eQTL,
    smr_p_SMR = candidate$p_SMR,
    smr_q_SMR_cell_BH = candidate$q_SMR_cell_BH,
    smr_p_HEIDI = candidate$p_HEIDI,
    eqtl_n = eqtl_n,
    coloc_status = status,
    coloc_error = error,
    coloc_nsnps = unname(coloc_summary["nsnps"]),
    PP.H0.abf = unname(coloc_summary["PP.H0.abf"]),
    PP.H1.abf = unname(coloc_summary["PP.H1.abf"]),
    PP.H2.abf = unname(coloc_summary["PP.H2.abf"]),
    PP.H3.abf = unname(coloc_summary["PP.H3.abf"]),
    PP.H4.abf = unname(coloc_summary["PP.H4.abf"]),
    PP4_over_PP3_PP4 = ifelse(
      is.finite(unname(coloc_summary["PP.H4.abf"])) &&
        is.finite(unname(coloc_summary["PP.H3.abf"])) &&
        (unname(coloc_summary["PP.H3.abf"]) + unname(coloc_summary["PP.H4.abf"])) > 0,
      unname(coloc_summary["PP.H4.abf"]) / (unname(coloc_summary["PP.H3.abf"]) + unname(coloc_summary["PP.H4.abf"])),
      NA_real_
    ),
    coloc_abf_positive_pp4_0.8 = is.finite(unname(coloc_summary["PP.H4.abf"])) && unname(coloc_summary["PP.H4.abf"]) >= 0.8,
    stringsAsFactors = FALSE
  )
  summary_rows[[i]] <- row
}

summary <- do.call(rbind, summary_rows)
summary_path <- file.path(qc_dir, "finngen_r13_discovery_fdr05_hits.coloc_abf_summary.tsv")
write.table(summary, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)

positive <- summary[summary$coloc_status == "OK" & summary$PP.H4.abf >= 0.8 & summary$PP4_over_PP3_PP4 >= 0.8, ]
positive_path <- file.path(qc_dir, "finngen_r13_discovery_fdr05_hits.coloc_abf_positive.tsv")
write.table(positive, positive_path, sep = "\t", row.names = FALSE, quote = FALSE)

report_path <- file.path(report_dir, "FINNGEN_R13_DISCOVERY_FDR05_HITS_COLOC_ABF_SUMMARY.md")
lines <- c(
  "# FinnGen R13 Discovery FDR<0.05 Hits: coloc.abf Summary",
  "",
  "Inputs:",
  "",
  "- Discovery GWAS: FinnGen R13 `I9_VTE`, GRCh37-lifted and OneK1K-LD harmonized.",
  "- eQTL: matching OneK1K pseudobulk qfile rows for each cell-gene candidate.",
  "- Candidate set: all FinnGen discovery SMR hits with cell-wise BH-FDR < 0.05.",
  "",
  "coloc settings:",
  "",
  "- `coloc.abf` from coloc 5.2.3.",
  "- GWAS dataset: case-control, `N=500186`, `s=27928/500186`.",
  "- eQTL dataset: quantitative, cell-specific donor N from `results/qc/eqtl/<cell>/chr*.full.summary.tsv`.",
  "- Positive rule used for this summary: `PP.H4.abf >= 0.8` and `PP4/(PP3+PP4) >= 0.8`.",
  "",
  "## Counts",
  "",
  sprintf("- Candidates tested: %d", nrow(summary)),
  sprintf("- Successful coloc runs: %d", sum(summary$coloc_status == "OK")),
  sprintf("- coloc.abf positive: %d", nrow(positive)),
  "",
  "## Outputs",
  "",
  sprintf("- Summary: `%s`", sub(paste0("^", project_dir, "/"), "", summary_path)),
  sprintf("- Positive hits: `%s`", sub(paste0("^", project_dir, "/"), "", positive_path)),
  sprintf("- Per-candidate merged inputs: `%s`", sub(paste0("^", project_dir, "/"), "", out_root))
)
writeLines(lines, report_path)

message("Wrote ", summary_path)
message("Wrote ", positive_path)
message("Wrote ", report_path)
