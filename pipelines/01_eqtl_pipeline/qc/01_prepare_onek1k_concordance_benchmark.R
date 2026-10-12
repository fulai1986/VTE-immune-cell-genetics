#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
official_dir <- file.path(project_dir, "reference/onek1k_official/esnp")
eqtl_dir <- file.path(project_dir, "results/eqtl")
out_dir <- file.path(project_dir, "manuscript_results/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_map <- data.table(
  official_cell_id = c(
    "cd4nc", "cd4et", "cd4sox4", "cd8nc", "cd8et", "cd8s100b",
    "nk", "nkr", "bin", "bmem", "monoc", "mononc", "dc", "plasma"
  ),
  cell_type = c(
    "CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 NC", "CD8 ET", "CD8 S100B",
    "NK", "NK R", "B IN", "B Mem", "Mono C", "Mono NC", "DC", "Plasma"
  ),
  cell_safe = c(
    "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET", "CD8_S100B",
    "NK", "NK_R", "B_IN", "B_Mem", "Mono_C", "Mono_NC", "DC", "Plasma"
  )
)

sample_official <- function(row, n = 500L) {
  path <- file.path(official_dir, paste0(row$official_cell_id, "_esnp_table.tsv.gz"))
  x <- fread(path, showProgress = FALSE, nThread = 1)
  x <- x[as.integer(ROUND) == 1L]
  x[, P_VALUE_NUM := as.numeric(P_VALUE)]
  setorder(x, P_VALUE_NUM)
  top <- x[seq_len(min(n, .N))]
  top[, sample_set := "official_top500_by_p"]
  set.seed(20260812L + match(row$official_cell_id, cell_map$official_cell_id))
  random <- x[sample.int(.N, min(n, .N))]
  random[, sample_set := "official_random500"]
  out <- rbindlist(list(top, random), use.names = TRUE, fill = TRUE)
  out[, `:=`(
    cell_type_project = row$cell_type,
    cell_safe = row$cell_safe,
    official_cell_id = row$official_cell_id,
    CHR_INT = as.integer(CHR),
    POS_INT = as.integer(POS),
    official_neglog10p = -log10(pmax(P_VALUE_NUM, 1e-300))
  )]
  out[, official_row_id := paste(sample_set, official_cell_id, GENE, CHR_INT, POS_INT, A1, A2, sep = "|")]
  out[, match_key := paste(cell_safe, CHR_INT, POS_INT, GENE, sep = "|")]
  out
}

official_sample <- rbindlist(lapply(seq_len(nrow(cell_map)), function(i) {
  sample_official(cell_map[i])
}), use.names = TRUE, fill = TRUE)

scan_one_cell_chr <- function(cell_safe, chrom, keys) {
  path <- file.path(eqtl_dir, cell_safe, paste0("chr", chrom, ".full.cis_eqtl.tsv.gz"))
  if (!file.exists(path) || length(keys) == 0L) return(NULL)
  message("Scanning ", cell_safe, " chr", chrom, " for ", length(keys), " official sampled pairs")
  use_cols <- c("chrom", "pos", "snp_id", "ref", "alt", "gene_symbol", "beta", "se", "p")
  x <- fread(path, select = use_cols, showProgress = FALSE, nThread = 1)
  x[, match_key := paste(cell_safe, chrom, pos, gene_symbol, sep = "|")]
  x[match_key %in% keys]
}

hits <- list()
idx <- 1L
for (cs in unique(official_sample$cell_safe)) {
  for (chrom in sort(unique(official_sample[cell_safe == cs]$CHR_INT))) {
    keys <- official_sample[cell_safe == cs & CHR_INT == chrom, unique(match_key)]
    hit <- scan_one_cell_chr(cs, chrom, keys)
  if (!is.null(hit) && nrow(hit) > 0L) {
      hits[[idx]] <- hit
      idx <- idx + 1L
    }
  }
}
project_hits <- rbindlist(hits, use.names = TRUE, fill = TRUE)
project_hits <- unique(project_hits, by = "match_key")

cmp <- merge(official_sample, project_hits, by = "match_key", all.x = TRUE)
cmp[, pair_found := is.finite(beta)]

cmp[, project_beta_on_official_a2 := fifelse(
  A2 == alt,
  as.numeric(beta),
  fifelse(A2 == ref, -as.numeric(beta), NA_real_)
)]
cmp[, direction_comparable := is.finite(project_beta_on_official_a2) & is.finite(as.numeric(SPEARMANS_RHO))]
cmp[, direction_concordant := direction_comparable & sign(as.numeric(SPEARMANS_RHO)) == sign(project_beta_on_official_a2)]
cmp[, ours_neglog10p := -log10(pmax(as.numeric(p), 1e-300))]
cmp[, our_nominal_p05 := as.numeric(p) < 0.05]

summary_dt <- cmp[, .(
  official_pairs_sampled = .N,
  pairs_found_in_ours = sum(pair_found),
  pair_found_rate = mean(pair_found),
  direction_comparable = sum(direction_comparable),
  direction_concordant = sum(direction_concordant),
  direction_concordance_rate = ifelse(sum(direction_comparable) > 1,
    mean(direction_concordant[direction_comparable]),
    NA_real_
  ),
  our_nominal_p05_among_found = ifelse(sum(pair_found) > 0,
    mean(our_nominal_p05[pair_found], na.rm = TRUE),
    NA_real_
  ),
  median_our_p_among_found = ifelse(sum(pair_found) > 0,
    median(as.numeric(p)[pair_found], na.rm = TRUE),
    NA_real_
  ),
  spearman_neglog10p_corr_found = ifelse(sum(pair_found) > 2,
    suppressWarnings(cor(official_neglog10p[pair_found], ours_neglog10p[pair_found], method = "spearman", use = "complete.obs")),
    NA_real_
  ),
  pearson_neglog10p_corr_found = ifelse(sum(pair_found) > 2,
    suppressWarnings(cor(official_neglog10p[pair_found], ours_neglog10p[pair_found], method = "pearson", use = "complete.obs")),
    NA_real_
  )
), by = .(cell_type = cell_type_project, cell_safe, sample_set)]

cell_order <- cell_map$cell_type
summary_dt[, cell_type := factor(cell_type, levels = cell_order)]
summary_dt[, sample_set := factor(sample_set, levels = c("official_top500_by_p", "official_random500"))]
setorder(summary_dt, cell_type, sample_set)

detail_path <- file.path(out_dir, "figure2c_official_top_random500_concordance_detail.tsv.gz")
summary_path <- file.path(out_dir, "figure2c_official_top_random500_concordance_summary.tsv")
fwrite(cmp, detail_path, sep = "\t")
fwrite(summary_dt, summary_path, sep = "\t")
message("Wrote ", detail_path)
message("Wrote ", summary_path)
print(summary_dt)
