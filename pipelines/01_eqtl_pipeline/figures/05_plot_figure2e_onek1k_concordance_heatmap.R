#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
summary_file <- file.path(project_dir, "manuscript_results/tables/figure2c_official_top_random500_concordance_summary.tsv")
egene_file <- file.path(project_dir, "manuscript_results/tables/figure2_eqtl_egene_summary_for_plot.tsv")
out_dir <- file.path(project_dir, "manuscript_results/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary <- fread(summary_file)
egene <- fread(egene_file)

dt <- summary[sample_set == "official_top500_by_p"]
dt <- merge(
  dt,
  egene[, .(cell_type, cells_retained, donors_retained, egene_eigenmt_fdr05)],
  by = "cell_type",
  all.x = TRUE
)

setorder(dt, -egene_eigenmt_fdr05, -donors_retained, -cells_retained)

cells <- dt$cell_type
metric_raw <- cbind(
  Recovery = 100 * dt$pair_found_rate,
  Direction = 100 * dt$direction_concordance_rate,
  `P < 0.05` = 100 * dt$our_nominal_p05_among_found,
  `Spearman rho` = dt$spearman_neglog10p_corr_found
)
metric_scaled <- metric_raw
metric_scaled[, 1:3] <- metric_scaled[, 1:3] / 100

metric_cols <- colorRampPalette(c("#FFF7EC", "#FEE8C8", "#FDBB84", "#FC8D59", "#E34A33", "#B30000"))(101)
to_col <- function(x, pal) pal[pmax(1, pmin(101, floor(x * 100) + 1))]

draw_colorbar <- function(x0, y0, w, h, pal, title, labels, ticks) {
  n <- length(pal)
  ys <- seq(y0, y0 + h, length.out = n + 1)
  for (i in seq_len(n)) {
    rect(x0, ys[i], x0 + w, ys[i + 1], col = pal[i], border = NA)
  }
  segments(x0, y0, x0, y0 + h, col = "#606060", lwd = 0.35)
  segments(x0 + w, y0, x0 + w, y0 + h, col = "#606060", lwd = 0.35)
  text(x0, y0 + h + 0.028, labels = title, adj = c(0, 0), cex = 0.66, col = "#202020")
  for (i in seq_along(ticks)) {
    y <- y0 + h * ticks[i]
    segments(x0 + w, y, x0 + w + 0.010, y, col = "#404040", lwd = 0.45)
    text(x0 + w + 0.017, y, labels = labels[i], adj = c(0, 0.5), cex = 0.64, col = "#202020")
  }
}

draw_plot <- function(file, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(file, width = 5.9, height = 4.35, units = "in", res = 600)
  } else {
    pdf(file, width = 5.9, height = 4.35, useDingbats = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 0, 0), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  n <- length(cells)
  y_top <- 0.88
  y_bottom <- 0.19
  row_h <- (y_top - y_bottom) / n
  label_x <- 0.245
  metric_x <- 0.285
  metric_w <- 0.115
  metric_gap <- 0.006

  for (i in seq_len(n)) {
    y1 <- y_top - i * row_h
    y2 <- y_top - (i - 1) * row_h
    ymid <- (y1 + y2) / 2

    text(label_x, ymid, cells[i], adj = c(1, 0.5), cex = 0.78, col = "#111111")

    for (j in seq_len(ncol(metric_scaled))) {
      x1 <- metric_x + (j - 1) * (metric_w + metric_gap)
      x2 <- x1 + metric_w
      val <- metric_scaled[i, j]
      rect(x1, y1, x2, y2, col = to_col(val, metric_cols), border = "white", lwd = 0.9)
      lab <- if (j <= 3) sprintf("%.0f%%", metric_raw[i, j]) else sprintf("%.2f", metric_raw[i, j])
      txt_col <- if (val >= 0.68) "white" else "#151515"
      text((x1 + x2) / 2, ymid, lab, col = txt_col, cex = 0.65, font = 2)
    }
  }

  text(
    metric_x, 0.955,
    "Strongest 500 ROUND=1 primary eSNP-eGene pairs",
    adj = c(0, 0.5), cex = 0.70, col = "#222222"
  )

  x_centers <- metric_x + (seq_len(ncol(metric_scaled)) - 1) * (metric_w + metric_gap) + metric_w / 2
  for (k in seq_along(colnames(metric_raw))) {
    text(x_centers[k], y_bottom - 0.036, colnames(metric_raw)[k], srt = 32, adj = c(1, 0.5), cex = 0.70, xpd = NA, col = "#111111")
  }

  draw_colorbar(0.820, 0.43, 0.018, 0.22, metric_cols, "Metric\nvalue", c("0", "0.5", "1"), c(0, 0.5, 1))
}

png_file <- file.path(out_dir, "figure2c_official_validation_summary_heatmap_polished_no_donor.png")
pdf_file <- file.path(out_dir, "figure2c_official_validation_summary_heatmap_polished_no_donor.pdf")
draw_plot(png_file, "png")
draw_plot(pdf_file, "pdf")
cat("Wrote", png_file, "\n")
cat("Wrote", pdf_file, "\n")
