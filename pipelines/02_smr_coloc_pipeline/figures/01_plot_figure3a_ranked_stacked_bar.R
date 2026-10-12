#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
input_file <- file.path(project_dir, "results/smr_analysis/qc/finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf_beta_or.tsv")
out_dir <- file.path(project_dir, "manuscript_results/figures")
table_dir <- file.path(project_dir, "manuscript_results/tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

cell_order <- c(
  "CD4 NC", "CD8 ET", "NK", "CD8 NC", "B IN", "CD4 ET", "B Mem",
  "Mono C", "CD8 S100B", "Mono NC", "NK R", "DC", "CD4 SOX4", "Plasma"
)
cell_safe_to_plot <- setNames(cell_order, gsub(" ", "_", cell_order))
tier_order <- c("Tier 1", "Tier 2", "Tier 3", "Not prioritized")

tier_palette <- c(
  "Tier 1" = "#8E2C7D",
  "Tier 2" = "#D95F02",
  "Tier 3" = "#E6AB02",
  "Not prioritized" = "#BDBDBD"
)

dt <- fread(input_file)
dt[, cell_type := fifelse(cell %in% names(cell_safe_to_plot), cell_safe_to_plot[cell], gsub("_", " ", cell))]
dt[, tier := factor(tier, levels = tier_order)]

summary_dt <- dt[, .N, by = .(cell_type, tier)]
complete_grid <- CJ(cell_type = cell_order, tier = tier_order)
complete_grid[, `:=`(
  cell_type = factor(cell_type, levels = cell_order),
  tier = factor(tier, levels = tier_order)
)]
summary_dt <- merge(complete_grid, summary_dt, by = c("cell_type", "tier"), all.x = TRUE)
summary_dt[is.na(N), N := 0L]
summary_dt[, cell_total := sum(N), by = cell_type]
summary_dt[, prioritized_total := sum(N[tier %in% c("Tier 1", "Tier 2", "Tier 3")]), by = cell_type]

fwrite(summary_dt, file.path(table_dir, "figure3a_finngen_smr_signal_distribution_alternatives.tsv"), sep = "\t")

theme_smr <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8.5, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.32, "cm"),
    plot.margin = margin(8, 8, 6, 6)
  )

rank_order <- unique(summary_dt[cell_total > 0][order(-cell_total, cell_type), as.character(cell_type)])
rank_dt <- copy(summary_dt[cell_total > 0])
rank_dt[, cell_ranked := factor(as.character(cell_type), levels = rev(rank_order))]

fig_ranked <- ggplot(rank_dt, aes(x = N, y = cell_ranked, fill = tier)) +
  geom_col(width = 0.66, color = "white", linewidth = 0.25) +
  geom_text(
    data = unique(rank_dt[, .(cell_ranked, cell_total)]),
    aes(x = cell_total + 0.12, y = cell_ranked, label = cell_total),
    inherit.aes = FALSE,
    size = 2.7,
    hjust = 0,
    color = "#2B2B2B"
  ) +
  scale_fill_manual(values = tier_palette, drop = FALSE) +
  scale_x_continuous(
    breaks = 0:6,
    limits = c(0, 6.7),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Number of cell-gene pairs with cell-wise BH-FDR < 0.05",
    y = NULL
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_smr +
  theme(
    legend.position = "top",
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

dot_dt <- copy(summary_dt)
dot_dt[, cell_type := factor(as.character(cell_type), levels = cell_order)]
dot_dt[, tier := factor(as.character(tier), levels = rev(tier_order))]
dot_dt[, shown_count := fifelse(N > 0, as.character(N), "")]

fig_matrix <- ggplot(dot_dt, aes(x = cell_type, y = tier)) +
  geom_tile(fill = "#F7F7F7", color = "white", linewidth = 0.5) +
  geom_point(aes(size = N, fill = tier), shape = 21, color = "#303030", stroke = 0.25, alpha = ifelse(dot_dt$N > 0, 0.95, 0)) +
  geom_text(aes(label = shown_count), size = 2.35, color = "white", fontface = "bold") +
  scale_fill_manual(values = tier_palette, drop = FALSE) +
  scale_size_area(max_size = 9, limits = c(0, max(dot_dt$N)), breaks = c(1, 2, 3), name = "Pairs") +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_fixed(ratio = 0.95) +
  guides(
    fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4, alpha = 1)),
    size = guide_legend(order = 2)
  ) +
  theme_smr +
  theme(
    axis.text.x = element_text(size = 7, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    panel.background = element_rect(fill = "white", color = NA)
  )

rank_pdf <- file.path(out_dir, "figure3a_alt1_ranked_stacked_bar.pdf")
rank_png <- file.path(out_dir, "figure3a_alt1_ranked_stacked_bar.png")
matrix_pdf <- file.path(out_dir, "figure3a_alt2_tier_dot_matrix.pdf")
matrix_png <- file.path(out_dir, "figure3a_alt2_tier_dot_matrix.png")

ggsave(rank_pdf, fig_ranked, width = 4.7, height = 2.7, device = "pdf", useDingbats = FALSE)
ggsave(rank_png, fig_ranked, width = 4.7, height = 2.7, dpi = 600)
ggsave(matrix_pdf, fig_matrix, width = 6.6, height = 3.0, device = "pdf", useDingbats = FALSE)
ggsave(matrix_png, fig_matrix, width = 6.6, height = 3.0, dpi = 600)

cat("Wrote", rank_pdf, "\n")
cat("Wrote", rank_png, "\n")
cat("Wrote", matrix_pdf, "\n")
cat("Wrote", matrix_png, "\n")
