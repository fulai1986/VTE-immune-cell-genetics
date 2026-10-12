#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
in_dir <- file.path(project_dir, "manuscript_results/tables/egene_eigenmt")
out_dir <- in_dir

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) return(default)
  args[[hit + 1L]]
}

cell_type <- get_arg("--cell-type", "CD8 ET")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
cell_safe <- safe_name(cell_type)

files <- Sys.glob(file.path(in_dir, paste0(cell_safe, ".chr*.eigenmt_like.by_gene.tsv")))
if (!length(files)) stop("No eigenMT-like by-gene files found for ", cell_type)

x <- rbindlist(lapply(files, fread), fill = TRUE)
x <- unique(x, by = c("cell_type", "chrom", "gene_id"))
x[, q_eigenmt_bonferroni_cell := p.adjust(p_gene_bonferroni_eigenmt, method = "BH")]
x[, q_eigenmt_sidak_cell := p.adjust(p_gene_sidak_eigenmt, method = "BH")]
x[, egene_eigenmt_fdr05 := q_eigenmt_bonferroni_cell < 0.05]
setorder(x, q_eigenmt_bonferroni_cell, chrom, gene_start, gene_id)

summary <- x[, .(
  cell_type = cell_type[1],
  chromosomes_completed = uniqueN(chrom),
  genes_tested = .N,
  median_n_cis_snps = median(n_cis_snps_qfile, na.rm = TRUE),
  median_meff_eigenmt_like = median(meff_eigenmt_like, na.rm = TRUE),
  median_meff_eigenmt_fraction = median(meff_eigenmt_fraction, na.rm = TRUE),
  nominal_gene_p05 = sum(p_gene_bonferroni_eigenmt < 0.05, na.rm = TRUE),
  egene_eigenmt_fdr05 = sum(egene_eigenmt_fdr05, na.rm = TRUE),
  min_p_gene = min(p_gene_bonferroni_eigenmt, na.rm = TRUE),
  min_q_cell = min(q_eigenmt_bonferroni_cell, na.rm = TRUE)
)]

prefix <- file.path(out_dir, paste0(cell_safe, ".eigenmt_like.cellwise"))
fwrite(x, paste0(prefix, ".by_gene.tsv.gz"), sep = "\t")
fwrite(summary, paste0(prefix, ".summary.tsv"), sep = "\t")
message("Wrote ", paste0(prefix, ".by_gene.tsv.gz"))
message("Wrote ", paste0(prefix, ".summary.tsv"))
