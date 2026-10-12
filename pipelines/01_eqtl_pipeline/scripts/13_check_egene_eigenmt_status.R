#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
in_dir <- file.path(project_dir, "manuscript_results/tables/egene_eigenmt")
out_file <- file.path(in_dir, "egene_eigenmt_status.tsv")

cells <- c(
  "CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 NC", "CD8 ET", "CD8 S100B",
  "NK", "NK R", "B IN", "B Mem", "Mono C", "Mono NC", "DC", "Plasma"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

rows <- lapply(cells, function(cell) {
  cell_safe <- safe_name(cell)
  chr_files <- Sys.glob(file.path(in_dir, paste0(cell_safe, ".chr*.eigenmt_like.by_gene.tsv")))
  chroms <- as.integer(sub(".*\\.chr([0-9]+)\\.eigenmt_like\\.by_gene\\.tsv$", "\\1", chr_files))
  summary_file <- file.path(in_dir, paste0(cell_safe, ".eigenmt_like.cellwise.summary.tsv"))
  if (file.exists(summary_file)) {
    s <- fread(summary_file)
  } else {
    s <- data.table(
      genes_tested = NA_integer_,
      egene_eigenmt_fdr05 = NA_integer_,
      median_meff_eigenmt_fraction = NA_real_
    )
  }
  data.table(
    cell_type = cell,
    chromosomes_completed = length(unique(chroms)),
    completed_chroms = paste(sort(unique(chroms)), collapse = ","),
    genes_tested_current = s$genes_tested[1],
    egene_eigenmt_fdr05_current = s$egene_eigenmt_fdr05[1],
    median_meff_fraction_current = s$median_meff_eigenmt_fraction[1]
  )
})

status <- rbindlist(rows, fill = TRUE)
fwrite(status, out_file, sep = "\t")
print(status)
message("Wrote ", out_file)
