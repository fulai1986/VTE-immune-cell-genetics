#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(gridExtra)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
input_file <- file.path(project_dir, "manuscript_results/tables/figure2_eqtl_egene_summary_for_plot.tsv")
out_dir <- file.path(project_dir, "manuscript_results/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_file)
dt <- dt[order(-cells_retained)]
dt[, cell_type := factor(cell_type, levels = cell_type)]
dt[, log10_cells_retained := log10(cells_retained)]

cell_palette <- c(
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

base_theme <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.text = element_text(color = "black", size = 7),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    axis.title = element_text(size = 8.5),
    plot.margin = margin(4, 6, 4, 6)
  )

p_cells <- ggplot(dt, aes(x = cell_type, y = log10_cells_retained)) +
  geom_col(aes(fill = cell_type), width = 0.62, alpha = 0.92) +
  scale_fill_manual(values = cell_palette, guide = "none") +
  scale_y_continuous(
    breaks = c(3, 4, 5, 6),
    limits = c(0, 6),
    expand = expansion(mult = c(0.01, 0.08))
  ) +
  labs(y = expression(log[10]~"(cells retained)"), x = NULL) +
  base_theme

p_donors <- ggplot(dt, aes(x = cell_type, y = donors_retained)) +
  geom_segment(aes(y = 0, yend = donors_retained, xend = cell_type),
               linewidth = 0.35, color = "#B8B8B8") +
  geom_point(aes(fill = cell_type), shape = 21, size = 2.3, stroke = 0.45,
             color = "#333333") +
  geom_text(aes(label = donors_retained),
            nudge_y = 28, size = 2.15, vjust = 0, color = "#333333") +
  scale_fill_manual(values = cell_palette, guide = "none") +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    limits = c(0, 1050),
    breaks = c(0, 250, 500, 750, 1000),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(y = "Donors contributing to pseudobulk", x = NULL) +
  base_theme

fig_split <- arrangeGrob(
  p_cells,
  p_donors,
  ncol = 2,
  widths = c(1, 1)
)

split_pdf <- file.path(out_dir, "figure2a_eqtl_sample_basis_split_axes.pdf")
split_png <- file.path(out_dir, "figure2a_eqtl_sample_basis_split_axes.png")
ggsave(split_pdf, fig_split, width = 6.8, height = 3.4, device = "pdf")
ggsave(split_png, fig_split, width = 6.8, height = 3.4, dpi = 600)

dt[, cells_log10 := log10(cells_retained)]

p_landscape <- ggplot(dt, aes(x = donors_retained, y = cells_retained)) +
  geom_point(aes(size = egene_eigenmt_fdr05), shape = 21, stroke = 0.5,
             fill = "#4C78A8", color = "#1B3553", alpha = 0.85) +
  geom_text(aes(label = cell_type), size = 2.4, nudge_y = 0.055 * max(dt$cells_retained),
            check_overlap = TRUE) +
  scale_y_log10(labels = label_number(big.mark = ",")) +
  scale_x_continuous(limits = c(0, 1050), breaks = c(0, 250, 500, 750, 1000)) +
  scale_size_area(max_size = 8, name = "eGenes") +
  labs(
    x = "Donors contributing to pseudobulk",
    y = "Cells retained for pseudobulk eQTL"
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.text = element_text(color = "black", size = 7),
    axis.title = element_text(size = 8.5),
    legend.position = "right",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7)
  )

landscape_pdf <- file.path(out_dir, "figure2a_eqtl_sample_basis_power_landscape.pdf")
landscape_png <- file.path(out_dir, "figure2a_eqtl_sample_basis_power_landscape.png")
ggsave(landscape_pdf, p_landscape, width = 4.9, height = 3.6, device = "pdf")
ggsave(landscape_png, p_landscape, width = 4.9, height = 3.6, dpi = 600)

cat("Wrote", split_pdf, "\n")
cat("Wrote", split_png, "\n")
cat("Wrote", landscape_pdf, "\n")
cat("Wrote", landscape_png, "\n")
