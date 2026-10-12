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
cell_colors <- c(
  "CD4 NC" = "#8E2C7D",
  "CD4 ET" = "#B06AA7",
  "CD4 SOX4" = "#C8B3D9",
  "CD8 ET" = "#1F77B4",
  "CD8 NC" = "#4C8DCC",
  "CD8 S100B" = "#8DBCE6",
  "NK" = "#3BAF57",
  "NK R" = "#93C47D",
  "Plasma" = "#B8D99A",
  "B Mem" = "#F4E04D",
  "B IN" = "#F2B233",
  "Mono C" = "#F28E2B",
  "Mono NC" = "#E8591C",
  "DC" = "#D71920"
)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
files <- file.path(egene_dir, paste0(safe_name(cell_order), ".eigenmt_like.cellwise.by_gene.tsv.gz"))
if (!all(file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("Missing eGene files: ", paste(missing, collapse = ", "))
}

egene <- rbindlist(lapply(files, function(f) {
  x <- fread(f, select = c("cell_type", "gene_id", "gene_symbol", "egene_eigenmt_fdr05"))
  x[egene_eigenmt_fdr05 == TRUE, .(cell_type, gene_id, gene_symbol)]
}), use.names = TRUE)
egene <- unique(egene)
egene[, cell_type := factor(cell_type, levels = cell_order)]

cell_totals <- egene[, .(n_egenes = uniqueN(gene_id)), by = cell_type]
cell_totals <- merge(data.table(cell_type = factor(cell_order, levels = cell_order)), cell_totals, by = "cell_type", all.x = TRUE)
cell_totals[is.na(n_egenes), n_egenes := 0L]

gene_sets <- egene[, .(
  n_cell_types = uniqueN(cell_type),
  cell_types = paste(as.character(sort(cell_type)), collapse = "|"),
  gene_symbol = gene_symbol[1]
), by = gene_id]
intersections <- gene_sets[, .(
  n_genes = .N,
  n_cell_types = first(n_cell_types)
), by = cell_types]
setorder(intersections, -n_genes, n_cell_types, cell_types)
intersections[, rank := .I]

top_n <- 30L
top_intersections <- intersections[1:min(top_n, .N)]
fwrite(intersections, file.path(tab_dir, "figure2d_egene_all_intersections.tsv"), sep = "\t")
fwrite(top_intersections, file.path(tab_dir, "figure2d_egene_top30_intersections_science_style.tsv"), sep = "\t")

membership <- matrix(FALSE, nrow = length(cell_order), ncol = nrow(top_intersections),
                     dimnames = list(cell_order, paste0("I", seq_len(nrow(top_intersections)))))
for (j in seq_len(nrow(top_intersections))) {
  members <- strsplit(top_intersections$cell_types[j], "\\|", fixed = FALSE)[[1]]
  membership[members, j] <- TRUE
}

draw_upset <- function(file, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(file, width = 7.2, height = 4.3, units = "in", res = 600)
  } else {
    pdf(file, width = 7.2, height = 4.3, useDingbats = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 0, 0), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  n_col <- ncol(membership)
  n_row <- nrow(membership)

  matrix_x0 <- 0.22
  matrix_x1 <- 0.76
  matrix_y0 <- 0.18
  matrix_y1 <- 0.52
  top_y0 <- 0.56
  top_y1 <- 0.86
  right_x0 <- 0.78
  right_x1 <- 0.94

  xs <- seq(matrix_x0, matrix_x1, length.out = n_col)
  ys <- seq(matrix_y1, matrix_y0, length.out = n_row)
  dx <- if (n_col > 1) diff(xs)[1] else 0.015
  dy <- if (n_row > 1) abs(diff(ys)[1]) else 0.020

  row_shade <- c("#FFFFFF", "#F6F6F6")
  for (i in seq_len(n_row)) {
    rect(matrix_x0 - dx * 0.55, ys[i] - dy * 0.45, matrix_x1 + dx * 0.55, ys[i] + dy * 0.45,
         col = row_shade[(i %% 2) + 1], border = NA)
  }

  for (i in seq_len(n_row)) {
    text(matrix_x0 - 0.018, ys[i], cell_order[i], adj = c(1, 0.5), cex = 0.62, col = "#333333")
  }

  bg_col <- "#E2E2E2"
  active_col <- "#243B8F"
  for (j in seq_len(n_col)) {
    active_rows <- which(membership[, j])
    if (length(active_rows) > 1) {
      segments(xs[j], min(ys[active_rows]), xs[j], max(ys[active_rows]), col = active_col, lwd = 1.0)
    }
    points(rep(xs[j], n_row), ys, pch = 16, cex = 0.42, col = bg_col)
    points(rep(xs[j], length(active_rows)), ys[active_rows], pch = 16, cex = 0.50, col = active_col)
  }

  top_max <- max(top_intersections$n_genes) * 1.10
  for (j in seq_len(n_col)) {
    h <- top_intersections$n_genes[j] / top_max * (top_y1 - top_y0)
    col <- if (top_intersections$n_cell_types[j] == 1) cell_colors[strsplit(top_intersections$cell_types[j], "\\|")[[1]]] else "#303030"
    rect(xs[j] - dx * 0.32, top_y0, xs[j] + dx * 0.32, top_y0 + h, col = col, border = NA)
  }
  axis_vals <- pretty(c(0, top_max), n = 4)
  axis_vals <- axis_vals[axis_vals >= 0 & axis_vals <= top_max]
  for (v in axis_vals) {
    y <- top_y0 + v / top_max * (top_y1 - top_y0)
    segments(matrix_x0 - dx * 0.55, y, matrix_x1 + dx * 0.55, y, col = "#E7E7E7", lwd = 0.45)
    text(matrix_x0 - 0.020, y, format(round(v), big.mark = ","), adj = c(1, 0.5), cex = 0.54, col = "#333333")
  }
  segments(matrix_x0 - dx * 0.55, top_y0, matrix_x1 + dx * 0.55, top_y0, col = "#333333", lwd = 0.5)
  segments(matrix_x0 - dx * 0.55, top_y0, matrix_x0 - dx * 0.55, top_y1, col = "#333333", lwd = 0.5)
  text(matrix_x0 - 0.067, (top_y0 + top_y1) / 2, "Intersection eGenes", srt = 90, cex = 0.56)

  right_max <- max(cell_totals$n_egenes) * 1.08
  for (i in seq_len(n_row)) {
    w <- cell_totals$n_egenes[i] / right_max * (right_x1 - right_x0)
    rect(right_x0, ys[i] - dy * 0.32, right_x0 + w, ys[i] + dy * 0.32,
         col = cell_colors[as.character(cell_totals$cell_type[i])], border = NA)
  }
  axis_vals_r <- pretty(c(0, right_max), n = 4)
  axis_vals_r <- axis_vals_r[axis_vals_r >= 0 & axis_vals_r <= right_max]
  for (v in axis_vals_r) {
    x <- right_x0 + v / right_max * (right_x1 - right_x0)
    segments(x, matrix_y0 - dy * 0.45, x, matrix_y1 + dy * 0.45, col = "#E7E7E7", lwd = 0.45)
    text(x, matrix_y0 - 0.035, format(round(v), big.mark = ","), adj = c(0.5, 1), cex = 0.52, col = "#333333")
  }
  segments(right_x0, matrix_y0 - dy * 0.45, right_x1, matrix_y0 - dy * 0.45, col = "#333333", lwd = 0.5)
  text((right_x0 + right_x1) / 2, matrix_y0 - 0.075, "Total eGenes", cex = 0.60)

}

png_file <- file.path(fig_dir, "figure2d_egene_sharing_upset_science_style.png")
pdf_file <- file.path(fig_dir, "figure2d_egene_sharing_upset_science_style.pdf")
draw_upset(png_file, "png")
draw_upset(pdf_file, "pdf")
cat("Wrote", png_file, "\n")
cat("Wrote", pdf_file, "\n")
