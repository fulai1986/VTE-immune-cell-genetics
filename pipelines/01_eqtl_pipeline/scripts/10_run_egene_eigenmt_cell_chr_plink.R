#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
eqtl_dir <- file.path(project_dir, "results/eqtl")
smr_dir <- file.path(project_dir, "results/smr")
ld_dir <- file.path(project_dir, "results/ld_reference/onek1k")
plink2 <- Sys.getenv("PLINK2_BIN", unset = file.path(project_dir, "tools/plink2/plink2"))
out_dir <- file.path(project_dir, "manuscript_results/tables/egene_eigenmt")
log_dir <- file.path(project_dir, "manuscript_results/logs/egene_eigenmt/plink_blocks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) return(default)
  args[[hit + 1L]]
}

cell_type <- get_arg("--cell-type", "CD8 ET")
chrom <- as.integer(get_arg("--chrom", "22"))
eigenmt_block_size <- as.integer(get_arg("--eigenmt-block-size", "200"))
eigenmt_var_thresh <- as.numeric(get_arg("--eigenmt-var-thresh", "0.99"))
limit_genes <- as.integer(get_arg("--limit-genes", "0"))

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
cell_safe <- safe_name(cell_type)

calc_pcs_for_variance <- function(mat, var_thresh) {
  if (ncol(mat) == 1L) return(1L)
  ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  ev[ev < 0 & ev > -1e-8] <- 0
  ev <- pmax(ev, 0)
  total <- sum(ev)
  if (!is.finite(total) || total <= 0) return(NA_integer_)
  which(cumsum(sort(ev, decreasing = TRUE)) / total >= var_thresh)[1]
}

message("Loading eQTL top-by-gene and qfile for ", cell_type, " chr", chrom)
top_file <- file.path(eqtl_dir, cell_safe, paste0("chr", chrom, ".full.top_by_gene.tsv"))
qfile <- file.path(
  smr_dir, cell_safe, paste0("chr", chrom),
  paste0(cell_safe, ".chr", chrom, ".r2_0.8.maf_0.01.allele_id.qfile.txt")
)
ld_prefix <- file.path(
  ld_dir, paste0("chr", chrom),
  paste0("onek1k.chr", chrom, ".r2_0.8.maf_0.01.allele_id")
)

top <- fread(top_file)
q <- unique(fread(qfile, select = c("SNP", "Probe", "BP")))
genes_available <- intersect(top$gene_id, unique(q$Probe))
top <- top[gene_id %in% genes_available]
setorder(top, chrom, gene_start, gene_id)
if (limit_genes > 0L) top <- top[seq_len(min(limit_genes, nrow(top)))]
gene_ids <- top$gene_id
q <- q[Probe %in% gene_ids]
setorder(q, Probe, BP, SNP)
gene_snps <- split(q$SNP, q$Probe)

tmp_root <- tempfile(pattern = paste0(cell_safe, "_chr", chrom, "_eigenmt_"))
dir.create(tmp_root, recursive = TRUE)
on.exit(unlink(tmp_root, recursive = TRUE, force = TRUE), add = TRUE)

run_block <- function(snps, gene, block_idx) {
  snps <- unique(snps)
  if (!length(snps)) return(list(n_ld_snps = 0L, pcs = NA_integer_))
  if (length(snps) == 1L) return(list(n_ld_snps = 1L, pcs = 1L))

  block_dir <- file.path(tmp_root, paste0(gsub("[^A-Za-z0-9]+", "_", gene), "_b", block_idx))
  extract_file <- paste0(block_dir, ".extract.txt")
  out_prefix <- block_dir
  writeLines(snps, extract_file)
  cmd <- c(
    shQuote(plink2),
    "--bfile", shQuote(ld_prefix),
    "--extract", shQuote(extract_file),
    "--r-unphased square",
    "--threads 1",
    "--out", shQuote(out_prefix)
  )
  log_file <- file.path(log_dir, paste0(cell_safe, ".chr", chrom, ".", gsub("[^A-Za-z0-9]+", "_", gene), ".block", block_idx, ".plink.log"))
  status <- system(paste(cmd, collapse = " "), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (status != 0L) {
    return(list(n_ld_snps = 0L, pcs = NA_integer_))
  }
  mat_file <- paste0(out_prefix, ".unphased.vcor1")
  var_file <- paste0(out_prefix, ".unphased.vcor1.vars")
  if (!file.exists(mat_file) || !file.exists(var_file)) {
    return(list(n_ld_snps = 0L, pcs = NA_integer_))
  }
  vars <- readLines(var_file, warn = FALSE)
  nvars <- length(vars)
  if (nvars == 0L) return(list(n_ld_snps = 0L, pcs = NA_integer_))
  if (nvars == 1L) return(list(n_ld_snps = 1L, pcs = 1L))
  mat <- as.matrix(fread(mat_file, header = FALSE))
  storage.mode(mat) <- "double"
  if (nrow(mat) != nvars || ncol(mat) != nvars) {
    return(list(n_ld_snps = nvars, pcs = NA_integer_))
  }
  mat[!is.finite(mat)] <- 0
  diag(mat) <- 1
  list(n_ld_snps = nvars, pcs = calc_pcs_for_variance(mat, eigenmt_var_thresh))
}

message("Computing PLINK-backed eigenMT-like Meff for ", length(gene_ids), " genes")
rows <- vector("list", length(gene_ids))
block_counter <- 0L
for (ii in seq_along(gene_ids)) {
  gene <- gene_ids[ii]
  snps <- unique(gene_snps[[gene]])
  blocks <- split(snps, ceiling(seq_along(snps) / eigenmt_block_size))
  block_results <- lapply(seq_along(blocks), function(bb) {
    block_counter <<- block_counter + 1L
    run_block(blocks[[bb]], gene, bb)
  })
  n_ld <- sum(vapply(block_results, `[[`, numeric(1), "n_ld_snps"), na.rm = TRUE)
  pcs <- sum(vapply(block_results, `[[`, numeric(1), "pcs"), na.rm = TRUE)
  if (ii %% 25L == 0L) {
    message("genes_done=", ii, "/", length(gene_ids), "; blocks_done=", block_counter)
  }
  p_lead <- top$p[ii]
  rows[[ii]] <- data.table(
    cell_type = cell_type,
    chrom = chrom,
    gene_id = gene,
    gene_symbol = top$gene_symbol[ii],
    gene_start = top$gene_start[ii],
    gene_end = top$gene_end[ii],
    top_snp = top$snp_id[ii],
    top_snp_pos = top$pos[ii],
    top_beta = top$beta[ii],
    top_se = top$se[ii],
    p_lead = p_lead,
    n_cis_snps_qfile = length(snps),
    n_ld_snps = n_ld,
    meff_eigenmt_like = pcs,
    eigenmt_block_size = eigenmt_block_size,
    eigenmt_var_thresh = eigenmt_var_thresh,
    eigenmt_blocks = length(blocks),
    meff_eigenmt_fraction = pcs / length(snps),
    p_gene_bonferroni_eigenmt = pmin(p_lead * pcs, 1),
    p_gene_sidak_eigenmt = 1 - exp(pcs * log1p(-p_lead))
  )
}

out <- rbindlist(rows, fill = TRUE)
summary <- out[, .(
  cell_type = cell_type[1],
  chrom = chrom[1],
  genes = .N,
  backend = "plink_ld_reference",
  total_blocks = sum(eigenmt_blocks, na.rm = TRUE),
  median_n_cis_snps = median(n_cis_snps_qfile, na.rm = TRUE),
  median_meff_eigenmt_like = median(meff_eigenmt_like, na.rm = TRUE),
  median_meff_eigenmt_fraction = median(meff_eigenmt_fraction, na.rm = TRUE),
  nominal_egene_p_gene_0.05 = sum(p_gene_bonferroni_eigenmt < 0.05, na.rm = TRUE)
)]

prefix <- file.path(out_dir, paste0(cell_safe, ".chr", chrom, ".eigenmt_like"))
fwrite(out, paste0(prefix, ".by_gene.tsv"), sep = "\t")
fwrite(summary, paste0(prefix, ".summary.tsv"), sep = "\t")
message("Wrote ", paste0(prefix, ".by_gene.tsv"))
message("Wrote ", paste0(prefix, ".summary.tsv"))
