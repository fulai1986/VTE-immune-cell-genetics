#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(png)
  library(scales)
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
heatmap_file <- file.path(project_dir, "manuscript_results/tables/figure3d_cross_cell_smr_heatmap_plot_data.tsv")
effect_dir <- file.path(project_dir, "manuscript_results/figures/smr_official_plots")
out_dir <- file.path(project_dir, "manuscript_results/figures")

cell_order <- c(
  "CD4 NC", "CD8 ET", "NK", "CD8 NC", "B IN", "CD4 ET", "B Mem",
  "Mono C", "CD8 S100B", "Mono NC", "NK R", "DC", "CD4 SOX4", "Plasma"
)
gene_order <- c("FADS1", "SLC44A2", "MAFF", "TMEM258", "TSPAN15")

hm <- fread(heatmap_file)
hm <- hm[as.character(Gene) %in% gene_order]
hm[, Gene_chr := as.character(Gene)]
hm[, cell_chr := as.character(cell_type)]
hm[, Gene := factor(Gene_chr, levels = rev(gene_order))]
hm[, cell_type := factor(cell_chr, levels = cell_order)]
hm[, fdr_sig := fifelse(is.na(fdr_sig), FALSE, fdr_sig)]
hm[, star_label := fifelse(fdr_sig, "*", "")]
hm[, beta_clip := pmax(pmin(b_SMR, 0.65), -0.65)]

base_theme <- theme_classic(base_size = 8, base_family = "sans")

heatmap_plot <- ggplot(hm, aes(x = cell_type, y = Gene)) +
  geom_tile(
    aes(fill = beta_clip),
    color = "white",
    linewidth = 0.65,
    width = 0.96,
    height = 0.9
  ) +
  geom_tile(
    data = hm[is.na(beta_clip)],
    fill = "#F0F0F0",
    color = "white",
    linewidth = 0.65,
    width = 0.96,
    height = 0.9
  ) +
  geom_text(
    data = hm[fdr_sig == TRUE],
    aes(label = star_label),
    size = 3.8,
    fontface = "bold",
    color = "#111111",
    vjust = 0.55
  ) +
  geom_tile(
    data = hm[(cell_chr == "CD8 ET" & Gene_chr == "SLC44A2") | (cell_chr == "NK" & Gene_chr == "FADS1")],
    fill = NA,
    color = "#1B1B1B",
    linewidth = 0.9,
    width = 0.98,
    height = 0.92
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.65, 0.65),
    oob = squish,
    na.value = "#F0F0F0",
    name = expression(beta[SMR]~"(log OR)")
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  guides(fill = guide_colorbar(
    title.position = "top",
    barwidth = unit(0.34, "cm"),
    barheight = unit(2.0, "cm"),
    ticks.colour = "#303030"
  )) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    axis.text.x = element_text(size = 7.2, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8.4, color = "black", face = "italic"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    plot.margin = margin(2, 8, 2, 6)
  )

read_effect <- function(path) {
  rasterGrob(readPNG(path), interpolate = TRUE)
}

effect_slc44a2 <- read_effect(file.path(effect_dir, "CD8_ET_SLC44A2_SMR_effect_line.official_peqtl1e4.png"))
effect_fads1 <- read_effect(file.path(effect_dir, "NK_FADS1_SMR_effect_line.official_peqtl1e4.png"))

out_png <- file.path(out_dir, "figure3d_cross_cell_smr_heatmap_with_official_effect_lines.png")
out_pdf <- file.path(out_dir, "figure3d_cross_cell_smr_heatmap_with_official_effect_lines.pdf")

draw_figure <- function() {
  grid.newpage()

  pushViewport(viewport(x = 0.50, y = 0.69, width = 0.90, height = 0.52))
  grid.draw(ggplotGrob(heatmap_plot))
  popViewport()

  pushViewport(viewport(x = 0.31, y = 0.22, width = 0.32, height = 0.36))
  grid.draw(effect_slc44a2)
  popViewport()

  pushViewport(viewport(x = 0.70, y = 0.22, width = 0.32, height = 0.36))
  grid.draw(effect_fads1)
  popViewport()

  # Anchors point to the highlighted heatmap tiles:
  # CD8 ET-SLC44A2 is column 2 and SLC44A2 row; NK-FADS1 is column 3 and FADS1 row.
  grid.lines(
    x = unit(c(0.238, 0.260, 0.315), "npc"),
    y = unit(c(0.765, 0.510, 0.395), "npc"),
    gp = gpar(col = "#4A4A4A", lwd = 1.0),
    arrow = arrow(type = "closed", length = unit(0.055, "inches"))
  )
  grid.lines(
    x = unit(c(0.284, 0.610, 0.695), "npc"),
    y = unit(c(0.835, 0.505, 0.395), "npc"),
    gp = gpar(col = "#4A4A4A", lwd = 1.0),
    arrow = arrow(type = "closed", length = unit(0.055, "inches"))
  )
}

png(out_png, width = 3600, height = 3000, res = 500, bg = "white")
draw_figure()
dev.off()

pdf(out_pdf, width = 7.2, height = 6.0, useDingbats = FALSE)
draw_figure()
dev.off()

cat(out_png, "\n", out_pdf, "\n", sep = "")
