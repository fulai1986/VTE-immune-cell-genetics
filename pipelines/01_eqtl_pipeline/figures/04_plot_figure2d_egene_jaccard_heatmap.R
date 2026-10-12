#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
egene_dir <- file.path(project_dir, "manuscript_results/tables/egene_eigenmt")
fig_dir <- file.path(project_dir, "manuscript_results/figures")
tab_dir <- file.path(project_dir, "manuscript_results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cell_order <- c(
  "CD4 NC", "CD4 ET", "CD4 SOX4",
  "CD8 ET", "CD8 NC", "CD8 S100B",
  "NK", "NK R",
  "Plasma", "B Mem", "B IN",
  "Mono C", "Mono NC", "DC"
)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
files <- file.path(egene_dir, paste0(safe_name(cell_order), ".eigenmt_like.cellwise.by_gene.tsv.gz"))
if (!all(file.exists(files))) {
  stop("Missing files: ", paste(files[!file.exists(files)], collapse = ", "))
}

egene <- rbindlist(lapply(files, function(f) {
  fread(
    f,
    select = c("cell_type", "gene_id", "gene_symbol", "top_snp", "top_beta", "egene_eigenmt_fdr05")
  )[egene_eigenmt_fdr05 == TRUE]
}), use.names = TRUE)
egene <- unique(egene)
egene[, cell_type := factor(cell_type, levels = cell_order)]

gene_sets <- lapply(cell_order, function(ct) unique(egene[as.character(cell_type) == ct, gene_id]))
names(gene_sets) <- cell_order

jaccard <- matrix(NA_real_, nrow = length(cell_order), ncol = length(cell_order), dimnames = list(cell_order, cell_order))
overlap_n <- matrix(NA_integer_, nrow = length(cell_order), ncol = length(cell_order), dimnames = list(cell_order, cell_order))
for (i in seq_along(cell_order)) {
  for (j in seq_along(cell_order)) {
    inter <- intersect(gene_sets[[i]], gene_sets[[j]])
    union <- union(gene_sets[[i]], gene_sets[[j]])
    jaccard[i, j] <- length(inter) / length(union)
    overlap_n[i, j] <- length(inter)
  }
}

beta_corr <- matrix(NA_real_, nrow = length(cell_order), ncol = length(cell_order), dimnames = list(cell_order, cell_order))
beta_pair_n <- matrix(NA_integer_, nrow = length(cell_order), ncol = length(cell_order), dimnames = list(cell_order, cell_order))
for (i in seq_along(cell_order)) {
  for (j in seq_along(cell_order)) {
    xi <- egene[as.character(cell_type) == cell_order[i], .(gene_id, top_snp_i = top_snp, beta_i = top_beta)]
    xj <- egene[as.character(cell_type) == cell_order[j], .(gene_id, top_snp_j = top_snp, beta_j = top_beta)]
    m <- merge(xi, xj, by = "gene_id")
    m <- m[top_snp_i == top_snp_j & is.finite(beta_i) & is.finite(beta_j)]
    beta_pair_n[i, j] <- nrow(m)
    beta_corr[i, j] <- if (nrow(m) >= 10) suppressWarnings(cor(m$beta_i, m$beta_j, method = "spearman")) else NA_real_
  }
}

fwrite(as.data.table(jaccard, keep.rownames = "cell_type"), file.path(tab_dir, "figure2e_egene_jaccard_matrix.tsv"), sep = "\t")
fwrite(as.data.table(overlap_n, keep.rownames = "cell_type"), file.path(tab_dir, "figure2e_egene_overlap_count_matrix.tsv"), sep = "\t")
fwrite(as.data.table(beta_corr, keep.rownames = "cell_type"), file.path(tab_dir, "figure2e_lead_beta_same_top_snp_spearman_matrix.tsv"), sep = "\t")
fwrite(as.data.table(beta_pair_n, keep.rownames = "cell_type"), file.path(tab_dir, "figure2e_lead_beta_same_top_snp_pair_count_matrix.tsv"), sep = "\t")

draw_colorbar <- function(x0, y0, w, h, pal, title, labels, ticks) {
  n <- length(pal)
  ys <- seq(y0, y0 + h, length.out = n + 1)
  for (i in seq_len(n)) {
    rect(x0, ys[i], x0 + w, ys[i + 1], col = pal[i], border = NA)
  }
  segments(x0, y0, x0, y0 + h, col = "#606060", lwd = 0.35)
  segments(x0 + w, y0, x0 + w, y0 + h, col = "#606060", lwd = 0.35)
  text(x0, y0 + h + 0.030, labels = title, adj = c(0, 0), cex = 0.62)
  for (i in seq_along(ticks)) {
    y <- y0 + h * ticks[i]
    segments(x0 + w, y, x0 + w + 0.010, y, col = "#404040", lwd = 0.45)
    text(x0 + w + 0.017, y, labels = labels[i], adj = c(0, 0.5), cex = 0.62)
  }
}

draw_heatmap <- function(mat, file, title = NULL, legend_title, pal, limits = c(0, 1), labels = TRUE, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(file, width = 4.7, height = 4.5, units = "in", res = 600)
  } else {
    pdf(file, width = 4.7, height = 4.5, useDingbats = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 0, 0), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  n <- nrow(mat)
  x0 <- 0.23
  x1 <- 0.79
  y0 <- 0.18
  y1 <- 0.82
  cell_w <- (x1 - x0) / n
  cell_h <- (y1 - y0) / n
  pal_n <- length(pal)
  map_col <- function(v) {
    if (!is.finite(v)) return("#F2F2F2")
    scaled <- (v - limits[1]) / (limits[2] - limits[1])
    pal[pmax(1, pmin(pal_n, floor(scaled * (pal_n - 1)) + 1))]
  }

  if (!is.null(title) && nzchar(title)) {
    text((x0 + x1) / 2, 0.92, title, cex = 0.78, font = 2)
  }
  for (i in seq_len(n)) {
    y_top <- y1 - (i - 1) * cell_h
    y_bottom <- y_top - cell_h
    text(x0 - 0.014, (y_top + y_bottom) / 2, rownames(mat)[i], adj = c(1, 0.5), cex = 0.58)
    text(x0 + (i - 0.5) * cell_w, y1 + 0.014, colnames(mat)[i], srt = 45, adj = c(0, 0.5), cex = 0.54, xpd = NA)
    for (j in seq_len(n)) {
      x_left <- x0 + (j - 1) * cell_w
      x_right <- x_left + cell_w
      val <- mat[i, j]
      rect(x_left, y_bottom, x_right, y_top, col = map_col(val), border = "white", lwd = 0.8)
      if (labels && i != j && is.finite(val) && val >= 0.20) {
        text((x_left + x_right) / 2, (y_top + y_bottom) / 2, sprintf("%.2f", val),
             cex = 0.42, col = ifelse(val >= 0.45, "white", "#111111"))
      }
    }
  }
  rect(x0, y0, x1, y1, border = "#4A4A4A", lwd = 0.5)
  draw_colorbar(0.84, 0.36, 0.018, 0.25, pal, legend_title, c(as.character(limits[1]), as.character(mean(limits)), as.character(limits[2])), c(0, 0.5, 1))
}

warm_pal <- colorRampPalette(c("#FFF7EC", "#FEE8C8", "#FDBB84", "#FC8D59", "#E34A33", "#B30000"))(101)
div_pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)

draw_heatmap(
  jaccard,
  file.path(fig_dir, "figure2e_egene_jaccard_heatmap.png"),
  "eGene sharing between cell types",
  "Jaccard\nindex",
  warm_pal,
  limits = c(0, 1),
  labels = TRUE,
  device = "png"
)
draw_heatmap(
  jaccard,
  file.path(fig_dir, "figure2e_egene_jaccard_heatmap.pdf"),
  "eGene sharing between cell types",
  "Jaccard\nindex",
  warm_pal,
  limits = c(0, 1),
  labels = TRUE,
  device = "pdf"
)
draw_heatmap(
  jaccard,
  file.path(fig_dir, "figure2e_egene_jaccard_heatmap_notitle.png"),
  NULL,
  "eGene\nJaccard\nindex",
  warm_pal,
  limits = c(0, 1),
  labels = TRUE,
  device = "png"
)
draw_heatmap(
  jaccard,
  file.path(fig_dir, "figure2e_egene_jaccard_heatmap_notitle.pdf"),
  NULL,
  "eGene\nJaccard\nindex",
  warm_pal,
  limits = c(0, 1),
  labels = TRUE,
  device = "pdf"
)
draw_heatmap(
  beta_corr,
  file.path(fig_dir, "figure2e_lead_beta_same_top_snp_spearman_heatmap.png"),
  "Lead-eQTL beta correlation",
  "Spearman\nrho",
  div_pal,
  limits = c(-1, 1),
  labels = FALSE,
  device = "png"
)
draw_heatmap(
  beta_corr,
  file.path(fig_dir, "figure2e_lead_beta_same_top_snp_spearman_heatmap.pdf"),
  "Lead-eQTL beta correlation",
  "Spearman\nrho",
  div_pal,
  limits = c(-1, 1),
  labels = FALSE,
  device = "pdf"
)

cat("Wrote Figure 2E heatmaps\n")
