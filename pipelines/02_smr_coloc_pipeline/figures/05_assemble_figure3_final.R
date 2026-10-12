#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(grid)
  library(png)
})

get_project_dir <- function() {
  env_dir <- Sys.getenv("PROJECT_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = FALSE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
    return(normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

project_dir <- get_project_dir()
fig_dir <- file.path(project_dir, "manuscript_results/figures")
completed_dir <- file.path(
  project_dir,
  "manuscript_results/completed_results/result2_smr_coloc/main"
)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(completed_dir, recursive = TRUE, showWarnings = FALSE)

panel_paths <- list(
  A = file.path(fig_dir, "figure3a_alt1_ranked_stacked_bar.png"),
  B = file.path(fig_dir, "figure3b_finngen_smr_tier1_3_forest_full.png"),
  C = file.path(fig_dir, "figure3d_cross_cell_smr_heatmap_with_official_effect_lines.png"),
  D = file.path(fig_dir, "figure4_tier1_coloc/figure4_tier1_coloc_reference_style_combined.png")
)

missing <- names(panel_paths)[!file.exists(unlist(panel_paths))]
if (length(missing)) {
  stop("Missing panel image(s): ", paste(missing, collapse = ", "))
}

read_panel <- function(path) rasterGrob(readPNG(path), interpolate = TRUE)

draw_panel <- function(path, label, x, y, w, h) {
  pushViewport(viewport(x = x, y = y, width = w, height = h, just = c("left", "bottom")))
  grid.draw(read_panel(path))
  grid.text(
    label,
    x = unit(0.01, "npc"),
    y = unit(0.99, "npc"),
    just = c("left", "top"),
    gp = gpar(fontface = "bold", fontsize = 20)
  )
  popViewport()
}

out_png <- file.path(fig_dir, "figure3_smr_coloc_final.png")
out_pdf <- file.path(fig_dir, "figure3_smr_coloc_final.pdf")

draw_figure <- function() {
  grid.newpage()
  draw_panel(panel_paths$A, "A", 0.02, 0.60, 0.43, 0.36)
  draw_panel(panel_paths$B, "B", 0.47, 0.58, 0.51, 0.39)
  draw_panel(panel_paths$C, "C", 0.02, 0.04, 0.45, 0.48)
  draw_panel(panel_paths$D, "D", 0.51, 0.04, 0.47, 0.49)
}

png(out_png, width = 3600, height = 2400, res = 300)
draw_figure()
dev.off()

pdf(out_pdf, width = 12, height = 8, useDingbats = FALSE)
draw_figure()
dev.off()

file.copy(out_png, file.path(completed_dir, "Figure3_smr_coloc.png"), overwrite = TRUE)
file.copy(out_pdf, file.path(completed_dir, "Figure3_smr_coloc.pdf"), overwrite = TRUE)

message("Wrote ", out_png)
message("Wrote ", out_pdf)
