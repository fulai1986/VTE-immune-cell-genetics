if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")
library(data.table)
library(Matrix)

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)

read_project_paths <- function(path) {
  x <- fread(path)
  setNames(x$path, x$key)
}

paths <- read_project_paths(file.path(project_dir, "config", "project_paths.tsv"))
tar_file <- paths[["raw_counts_tar"]]
raw_dir <- file.path(project_dir, "scratch", "GSE196735_RAW")
out_dir <- file.path(project_dir, "results", "pseudobulk")
qc_dir <- file.path(project_dir, "results", "qc")
report_dir <- file.path(project_dir, "reports")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

available_pools <- c(1:39, 41:65, 67:77)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) return(default)
  args[[hit + 1L]]
}

parse_pools <- function(x) {
  if (is.null(x) || x == "") return(1:3)
  if (tolower(x) == "all") return(available_pools)
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  pools <- integer()
  for (part in parts) {
    part <- trimws(part)
    if (grepl("-", part, fixed = TRUE)) {
      ends <- as.integer(strsplit(part, "-", fixed = TRUE)[[1]])
      pools <- c(pools, seq.int(ends[[1]], ends[[2]]))
    } else {
      pools <- c(pools, as.integer(part))
    }
  }
  intersect(unique(pools), available_pools)
}

target_cell_type <- get_arg("--cell-type", "NK")
pools_to_use <- parse_pools(get_arg("--pools", "1-3"))
min_cells_per_donor <- as.integer(get_arg("--min-cells-per-donor", "10"))
min_donors_expressed <- as.integer(get_arg("--min-donors-expressed", "20"))

cell_meta_file <- file.path(project_dir, "metadata", "cell_metadata_science14_author_style.tsv.gz")

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
pool_label <- function(pools) {
  if (setequal(pools, available_pools)) return("all_pools")
  paste0("pools_", min(pools), "_", max(pools))
}
rank_int <- function(x) {
  ok <- is.finite(x)
  y <- rep(NA_real_, length(x))
  if (sum(ok) >= 3) {
    r <- rank(x[ok], ties.method = "average")
    y[ok] <- qnorm((r - 0.5) / length(r))
  }
  y
}

extract_pool_files <- function(pools) {
  all_files <- untar(tar_file, list = TRUE)
  patterns <- paste0("(Pool|Sample)(", paste(pools, collapse = "|"), ")_")
  keep <- all_files[grepl(patterns, all_files) & grepl("RawCounts|Individual_Barcodes|GenotypeSamples", all_files)]
  missing <- keep[!file.exists(file.path(raw_dir, keep))]
  if (length(missing) > 0) {
    message("Extracting ", length(missing), " files from GEO tar...")
    untar(tar_file, files = missing, exdir = raw_dir)
  }
}

pool_file_table <- function() {
  files <- list.files(raw_dir, full.names = TRUE)
  data.table(path = files, file = basename(files))[
    grepl("RawCounts|Individual_Barcodes|GenotypeSamples", file)
  ][
    , pool := as.integer(sub(".*(?:Pool|Sample)([0-9]+).*", "\\1", file))
  ][
    , type := fifelse(grepl("RawCounts", file), "counts",
               fifelse(grepl("Individual_Barcodes", file), "barcode", "sample_list"))
  ]
}

sum_counts_by_donor <- function(count_mat, cell_to_donor) {
  donor <- as.character(cell_to_donor[colnames(count_mat)])
  donor_levels <- unique(donor)
  design <- sparseMatrix(
    i = seq_along(donor),
    j = match(donor, donor_levels),
    x = 1,
    dims = c(length(donor), length(donor_levels)),
    dimnames = list(NULL, donor_levels)
  )
  count_mat %*% design
}

read_counts_for_cells <- function(path, wanted_barcodes) {
  header <- names(fread(path, nrows = 0, check.names = FALSE))
  gene_col <- header[1]
  cols <- intersect(header[-1], wanted_barcodes)
  if (length(cols) == 0) return(NULL)
  x <- fread(path, select = c(gene_col, cols), check.names = FALSE, showProgress = TRUE)
  genes <- x[[1]]
  mat <- as.matrix(x[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(as.character(genes))
  mat
}

message("Target cell type: ", target_cell_type)
message("Pools: ", paste(pools_to_use, collapse = ","))
extract_pool_files(pools_to_use)
files_dt <- pool_file_table()

cell_meta <- fread(cell_meta_file, select = c(
  "barcode",
  "donor_id",
  "pool_number",
  "cell_type_cellxgene",
  "predicted_celltype_l2",
  "science14_cell_type",
  "classification_note",
  "marker_score_for_assigned_label"
))
setnames(cell_meta,
         old = c("donor_id", "pool_number"),
         new = c("individual_id", "pool"))
cell_meta[, pool := as.integer(pool)]
cell_meta[, barcode_h5ad := barcode]
cell_meta <- cell_meta[science14_cell_type == target_cell_type & pool %in% pools_to_use]
if (nrow(cell_meta) == 0) stop("No cells found for target_cell_type in selected pools.")

donor_counts <- cell_meta[, .N, by = individual_id]
keep_donors <- donor_counts[N >= min_cells_per_donor, individual_id]
cell_meta <- cell_meta[individual_id %in% keep_donors]

message("Cells retained: ", nrow(cell_meta))
message("Donors retained: ", uniqueN(cell_meta$individual_id))

pb_list <- list()
for (p in pools_to_use) {
  message("Processing pool ", p)
  counts_file <- files_dt[pool == p & type == "counts", path][1]
  if (is.na(counts_file)) stop("Missing RawCounts for pool ", p)

  cm <- cell_meta[pool == p]
  if (nrow(cm) == 0) {
    message("  no retained cells for this pool")
    next
  }

  count_barcodes <- if ("barcode_h5ad" %in% names(cm)) cm$barcode_h5ad else cm$barcode
  mat <- read_counts_for_cells(counts_file, count_barcodes)
  if (is.null(mat)) {
    message("  no count columns matched")
    next
  }

  cell_to_donor <- setNames(cm$individual_id, count_barcodes)
  pb <- sum_counts_by_donor(mat, cell_to_donor)
  pb_list[[as.character(p)]] <- pb
  rm(mat, pb)
  gc()
}

if (length(pb_list) == 0) stop("No pseudo-bulk matrix created.")

all_genes <- Reduce(union, lapply(pb_list, rownames))
all_donors <- Reduce(union, lapply(pb_list, colnames))
pb_all <- Matrix(0, nrow = length(all_genes), ncol = length(all_donors), sparse = TRUE)
rownames(pb_all) <- all_genes
colnames(pb_all) <- all_donors

for (pb in pb_list) {
  pb_all[rownames(pb), colnames(pb)] <- pb_all[rownames(pb), colnames(pb)] + pb
}

out_prefix <- paste0(safe_name(target_cell_type), ".", pool_label(pools_to_use))

raw_out <- file.path(out_dir, paste0(out_prefix, ".raw_counts.tsv.gz"))
raw_dt <- as.data.table(as.matrix(pb_all), keep.rownames = "Probe")
setnames(raw_dt, c("Probe", colnames(pb_all)))
fwrite(raw_dt, raw_out, sep = "\t")

lib_size <- Matrix::colSums(pb_all)
cpm <- t(t(pb_all) / lib_size * 1e6)
expr <- log2(as.matrix(cpm) + 1)
keep_gene <- rowSums(expr > 0) >= min_donors_expressed
expr <- expr[keep_gene, , drop = FALSE]
expr <- t(apply(expr, 1, rank_int))
rownames(expr) <- rownames(pb_all)[keep_gene]
colnames(expr) <- colnames(pb_all)

expr_out <- file.path(out_dir, paste0(out_prefix, ".expression_rankint.tsv.gz"))
expr_dt <- as.data.table(expr, keep.rownames = "Probe")
setnames(expr_dt, c("Probe", colnames(expr)))
fwrite(expr_dt, expr_out, sep = "\t")

summary <- data.table(
  metric = c("target_cell_type", "pools", "cells_retained", "donors_retained", "genes_raw", "genes_expression"),
  value = c(
    target_cell_type,
    pool_label(pools_to_use),
    nrow(cell_meta),
    uniqueN(cell_meta$individual_id),
    nrow(pb_all),
    nrow(expr)
  )
)
summary_out <- file.path(qc_dir, paste0(out_prefix, ".pseudobulk_summary.tsv"))
fwrite(summary, summary_out, sep = "\t")

donor_qc <- data.table(
  individual_id = colnames(pb_all),
  library_size = as.numeric(lib_size[colnames(pb_all)]),
  cells_in_celltype = as.integer(donor_counts[match(colnames(pb_all), individual_id), N])
)
donor_qc_out <- file.path(qc_dir, paste0(out_prefix, ".pseudobulk_donor_qc.tsv"))
fwrite(donor_qc, donor_qc_out, sep = "\t")

report_out <- file.path(report_dir, paste0("PSEUDOBULK_", out_prefix, ".md"))
writeLines(c(
  paste0("# Pseudobulk ", target_cell_type),
  "",
  paste0("- Pools: ", pool_label(pools_to_use)),
  paste0("- Cells retained: ", nrow(cell_meta)),
  paste0("- Donors retained: ", uniqueN(cell_meta$individual_id)),
  paste0("- Raw genes: ", nrow(pb_all)),
  paste0("- Expression genes after filtering: ", nrow(expr)),
  paste0("- Raw counts: `", sub(paste0(project_dir, "/"), "", raw_out, fixed = TRUE), "`"),
  paste0("- RankINT expression: `", sub(paste0(project_dir, "/"), "", expr_out, fixed = TRUE), "`"),
  paste0("- Summary: `", sub(paste0(project_dir, "/"), "", summary_out, fixed = TRUE), "`"),
  paste0("- Donor QC: `", sub(paste0(project_dir, "/"), "", donor_qc_out, fixed = TRUE), "`")
), report_out)

message("Wrote raw pseudo-bulk counts: ", raw_out)
message("Wrote normalized expression: ", expr_out)
message("Wrote summary: ", summary_out)
message("Wrote donor QC: ", donor_qc_out)
message("Wrote report: ", report_out)
print(summary)
