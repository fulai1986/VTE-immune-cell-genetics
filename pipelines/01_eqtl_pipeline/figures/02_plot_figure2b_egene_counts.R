#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) dirname(normalizePath(sub("^--file=", "", script_file[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
input_file <- file.path(project_dir, "manuscript_results/tables/figure2_eqtl_egene_summary_for_plot.tsv")
out_dir <- file.path(project_dir, "manuscript_results/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_file)
dt <- dt[order(-egene_eigenmt_fdr05)]
dt[, cell_type := factor(cell_type, levels = cell_type)]
dt[, egene_label := format(egene_eigenmt_fdr05, big.mark = ",")]

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

fig <- ggplot(dt, aes(x = cell_type, y = egene_eigenmt_fdr05)) +
  geom_col(aes(fill = cell_type), width = 0.68, alpha = 0.94) +
  geom_text(
    aes(label = egene_label),
    vjust = -0.25,
    size = 2.45,
    color = "#333333"
  ) +
  scale_fill_manual(values = cell_palette, guide = "none") +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Number of eGenes"
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 8.5, margin = margin(r = 6)),
    plot.margin = margin(6, 8, 6, 6)
  )

pdf_file <- file.path(out_dir, "figure2b_egene_counts.pdf")
png_file <- file.path(out_dir, "figure2b_egene_counts.png")
ggsave(pdf_file, fig, width = 6.4, height = 3.9, device = "pdf")
ggsave(png_file, fig, width = 6.4, height = 3.9, dpi = 600)

cat("Wrote", pdf_file, "\n")
cat("Wrote", png_file, "\n")
