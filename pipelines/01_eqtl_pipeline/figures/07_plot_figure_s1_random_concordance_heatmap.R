#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
tab_file <- file.path(project_dir, "manuscript_results/tables/figure2c_official_top_random500_concordance_summary.tsv")
out_dir <- file.path(project_dir, "manuscript_results/figures/supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_order <- c(
  "CD4 NC", "CD8 ET", "NK", "CD8 NC", "B IN", "CD4 ET", "B Mem",
  "Mono C", "CD8 S100B", "Mono NC", "NK R", "DC", "CD4 SOX4", "Plasma"
)

dt <- fread(tab_file)
dt[, cell_type := factor(cell_type, levels = cell_order)]
dt[, benchmark := fifelse(sample_set == "official_top500_by_p", "Strongest", "Random")]
dt[, benchmark := factor(benchmark, levels = c("Strongest", "Random"))]
dt <- dt[benchmark == "Random"]

plot_dt <- rbindlist(list(
  dt[, .(cell_type, benchmark, pairs = official_pairs_sampled, metric = "Direction", value = direction_concordance_rate)],
  dt[, .(cell_type, benchmark, pairs = official_pairs_sampled, metric = "P < 0.05", value = our_nominal_p05_among_found)],
  dt[, .(cell_type, benchmark, pairs = official_pairs_sampled, metric = "Spearman rho", value = spearman_neglog10p_corr_found)]
), use.names = TRUE)
plot_dt[, metric := factor(metric, levels = c("Direction", "P < 0.05", "Spearman rho"))]
setorder(plot_dt, cell_type, benchmark, metric)

heat_dt <- dcast(plot_dt, cell_type ~ benchmark + metric, value.var = "value")
mat_cols <- c(
  "Random_Direction", "Random_P < 0.05", "Random_Spearman rho"
)
mat <- as.matrix(heat_dt[, ..mat_cols])
rownames(mat) <- as.character(heat_dt$cell_type)
colnames(mat) <- c("Direction", "P < 0.05", "rho")

pair_lookup <- unique(dt[, .(cell_type, benchmark, pairs = official_pairs_sampled)])
random_pair <- pair_lookup[benchmark == "Random"][match(rownames(mat), as.character(cell_type)), pairs]

pal <- colorRampPalette(c("#FFF7EC", "#FEE8C8", "#FDBB84", "#FC8D59", "#E34A33", "#B30000"))(101)
to_col <- function(x) pal[pmax(1, pmin(101, floor(x * 100) + 1))]

draw_colorbar <- function(x0, y0, w, h) {
  ys <- seq(y0, y0 + h, length.out = length(pal) + 1)
  for (i in seq_along(pal)) {
    rect(x0, ys[i], x0 + w, ys[i + 1], col = pal[i], border = NA)
  }
  rect(x0, y0, x0 + w, y0 + h, border = "#555555", lwd = 0.4)
  text(x0, y0 + h + 0.028, "Metric\nvalue", adj = c(0, 0), cex = 0.66)
  for (v in c(0, 0.5, 1)) {
    y <- y0 + h * v
    segments(x0 + w, y, x0 + w + 0.010, y, col = "#333333", lwd = 0.4)
    text(x0 + w + 0.016, y, v, adj = c(0, 0.5), cex = 0.62)
  }
}

draw_heatmap <- function(file, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(file, width = 4.8, height = 5.2, units = "in", res = 600)
  } else {
    pdf(file, width = 4.8, height = 5.2, useDingbats = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 0, 0), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  n_row <- nrow(mat)
  n_col <- ncol(mat)
  x0 <- 0.31
  x1 <- 0.72
  y0 <- 0.16
  y1 <- 0.83
  cw <- (x1 - x0) / n_col
  rh <- (y1 - y0) / n_row

  text((x0 + x1) / 2, 0.940, "Random pairs", cex = 0.74, font = 2)
  text((x0 + x1) / 2, 0.902, "up to 500; all available if <500", cex = 0.50)

  for (j in seq_len(n_col)) {
    text(x0 + (j - 0.5) * cw, y1 + 0.027, colnames(mat)[j],
         adj = c(0.5, 0.5), cex = 0.56, xpd = NA)
  }
  for (i in seq_len(n_row)) {
    y_top <- y1 - (i - 1) * rh
    y_bot <- y_top - rh
    text(x0 - 0.018, (y_top + y_bot) / 2, rownames(mat)[i], adj = c(1, 0.5), cex = 0.66)
    text(x1 + 0.020, (y_top + y_bot) / 2,
         sprintf("n=%s", random_pair[i]), adj = c(0, 0.5), cex = 0.54, col = "#333333")
    for (j in seq_len(n_col)) {
      x_left <- x0 + (j - 1) * cw
      x_right <- x_left + cw
      val <- mat[i, j]
      rect(x_left, y_bot, x_right, y_top, col = to_col(val), border = "white", lwd = 0.8)
      lab <- if (j %in% c(1, 2)) sprintf("%.0f%%", val * 100) else sprintf("%.2f", val)
      text((x_left + x_right) / 2, (y_top + y_bot) / 2, lab,
           cex = 0.50, font = 2, col = ifelse(val >= 0.68, "white", "#111111"))
    }
  }
  rect(x0, y0, x1, y1, border = "#4A4A4A", lwd = 0.5)
  text(x1 + 0.020, y1 + 0.027, "Pairs used", adj = c(0, 0.5), cex = 0.58)
  draw_colorbar(0.89, 0.42, 0.020, 0.24)
}

png_file <- file.path(out_dir, "Figure_S1_onek1k_eqtl_concordance_benchmark.png")
pdf_file <- file.path(out_dir, "Figure_S1_onek1k_eqtl_concordance_benchmark.pdf")
draw_heatmap(png_file, "png")
draw_heatmap(pdf_file, "pdf")

summary_file <- file.path(project_dir, "manuscript_results/tables/supplementary/Figure_S1_onek1k_eqtl_concordance_benchmark_plot_data.tsv")
fwrite(plot_dt, summary_file, sep = "\t")

cat("Wrote", png_file, "\n")
cat("Wrote", pdf_file, "\n")
cat("Wrote", summary_file, "\n")
