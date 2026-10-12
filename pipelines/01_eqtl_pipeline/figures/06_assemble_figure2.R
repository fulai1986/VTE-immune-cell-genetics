#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(grid)
  library(png)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
fig_dir <- file.path(project_dir, "manuscript_results/figures")

panels <- list(
  A = file.path(fig_dir, "figure2a_eqtl_sample_basis_split_axes.png"),
  B = file.path(fig_dir, "figure2b_egene_counts.png"),
  C = file.path(fig_dir, "figure2d_egene_sharing_upset_science_style.png"),
  D = file.path(fig_dir, "figure2e_egene_jaccard_heatmap_notitle.png"),
  E = file.path(fig_dir, "figure2c_official_validation_summary_heatmap_polished_no_donor.png")
)

stopifnot(all(file.exists(unlist(panels))))

trim_white <- function(img, tol = 0.985, pad = 20L) {
  rgb <- img[, , 1:3, drop = FALSE]
  alpha_mask <- if (dim(img)[3] == 4) img[, , 4] > 0.01 else TRUE
  nonwhite <- apply(rgb < tol, c(1, 2), any) & alpha_mask
  rows <- which(rowSums(nonwhite) > 0)
  cols <- which(colSums(nonwhite) > 0)
  if (!length(rows) || !length(cols)) {
    return(img)
  }
  r0 <- max(1L, min(rows) - pad)
  r1 <- min(dim(img)[1], max(rows) + pad)
  c0 <- max(1L, min(cols) - pad)
  c1 <- min(dim(img)[2], max(cols) + pad)
  img[r0:r1, c0:c1, , drop = FALSE]
}

read_panel <- function(path, pad = 20L) {
  rasterGrob(trim_white(readPNG(path), pad = pad), interpolate = TRUE)
}

panel_grob <- function(path, label) {
  g <- read_panel(path)
  g$vp <- viewport(
    x = unit(0.525, "npc"),
    y = unit(0.475, "npc"),
    width = unit(0.94, "npc"),
    height = unit(0.90, "npc")
  )
  grobTree(
    g,
    textGrob(label, x = unit(0.018, "npc"), y = unit(0.965, "npc"),
             just = c("left", "top"), gp = gpar(fontsize = 15, fontface = "bold"))
  )
}

panel_grob_scaled <- function(path, label, scale = 1) {
  g <- read_panel(path)
  g$vp <- viewport(
    x = unit(0.525, "npc"),
    y = unit(0.475, "npc"),
    width = unit(0.94 * scale, "npc"),
    height = unit(0.90 * scale, "npc")
  )
  grobTree(
    g,
    textGrob(label, x = unit(0.018, "npc"), y = unit(0.965, "npc"),
             just = c("left", "top"), gp = gpar(fontsize = 15, fontface = "bold"))
  )
}

out_png <- file.path(fig_dir, "figure2_combined_preview.png")
out_pdf <- file.path(fig_dir, "figure2_combined_preview.pdf")

draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 13,
    heights = unit(c(1.00, 1.18), "null"),
    widths = unit(rep(1, 13), "null")
  )))

  # Top row: analysis basis and detected eGenes.
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:6))
  grid.draw(panel_grob(panels$A, "A"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 6:13))
  grid.draw(panel_grob(panels$B, "B"))
  popViewport()

  # Bottom row: broad sharing, pairwise sharing, and OneK1K reference concordance.
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:5))
  grid.draw(panel_grob(panels$C, "C"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 6:9))
  grid.draw(panel_grob_scaled(panels$D, "D", scale = 0.92))
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 10:13))
  grid.draw(panel_grob_scaled(panels$E, "E", scale = 0.96))
  popViewport()
}

png(out_png, width = 15.8, height = 8.8, units = "in", res = 300)
draw_combined()
dev.off()

pdf(out_pdf, width = 15.8, height = 8.8, useDingbats = FALSE)
draw_combined()
dev.off()

cat("Wrote", out_png, "\n")
cat("Wrote", out_pdf, "\n")
