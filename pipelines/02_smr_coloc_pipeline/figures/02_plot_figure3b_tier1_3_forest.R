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
tier_file <- file.path(project_dir, "results/smr_analysis/qc/finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf_beta_or.tsv")
hit_file <- file.path(project_dir, "results/smr_analysis/qc/finngen_r13_discovery_onek1k_params.cellwise_bh_fdr_lt_0.05_hits.tsv")
out_dir <- file.path(project_dir, "manuscript_results/figures")
table_dir <- file.path(project_dir, "manuscript_results/tables")

tier_order <- c("Tier 1", "Tier 2", "Tier 3")
tier_bg <- c(
  "Tier 1" = "#F5EDF3",
  "Tier 2" = "#FBF0DF",
  "Tier 3" = "#F7F2DE"
)
tier_col <- c(
  "Tier 1" = "#8E3B46",
  "Tier 2" = "#C76E00",
  "Tier 3" = "#7A8425"
)

dt <- fread(tier_file)[tier %in% tier_order]
hits <- fread(hit_file)[, .(cell, Gene, topSNP, se_SMR)]
dt <- merge(dt, hits, by = c("cell", "Gene", "topSNP"), all.x = TRUE)
if (any(is.na(dt$se_SMR))) stop("Missing se_SMR after merge")

dt[, tier := factor(tier, levels = tier_order)]
dt[, Cell := gsub("_", " ", cell)]
dt[, `:=`(
  ci_low = exp(b_SMR_finngen - 1.96 * se_SMR),
  ci_high = exp(b_SMR_finngen + 1.96 * se_SMR),
  OR = OR_SMR_finngen,
  or_label = sprintf("%.2f (%.2f-%.2f)", OR_SMR_finngen,
                     exp(b_SMR_finngen - 1.96 * se_SMR),
                     exp(b_SMR_finngen + 1.96 * se_SMR)),
  fdr_label = ifelse(smr_q_SMR_cell_BH < 0.001,
                     formatC(smr_q_SMR_cell_BH, format = "e", digits = 1),
                     sprintf("%.3f", smr_q_SMR_cell_BH)),
  heidi_label = ifelse(smr_p_HEIDI >= 0.01, "Pass", "Fail"),
  pp4_label = sprintf("%.3f", PP.H4.abf)
)]

dt <- dt[order(as.integer(tier), Gene, Cell)]
dt[, y := .N:1]
dt[, pair_label := paste0(Cell, " - ", Gene)]

bg_dt <- dt[, .(ymin = min(y) - 0.48, ymax = max(y) + 0.48), by = tier]
tier_labels <- dt[, .(y = mean(y)), by = tier]

fwrite(
  dt[, .(
    tier, Cell, Gene, topSNP, b_SMR_finngen, se_SMR, OR, ci_low, ci_high,
    or_label, smr_p_SMR, smr_q_SMR_cell_BH, smr_p_HEIDI, PP.H4.abf
  )],
  file.path(table_dir, "figure3b_finngen_smr_tier1_3_forest_full_plot_data.tsv"),
  sep = "\t"
)

forest_low <- 0.40
forest_high <- 1.50
forest_left <- 0.345
forest_right <- 0.585
map_or <- function(x) {
  forest_left + (log(x) - log(forest_low)) /
    (log(forest_high) - log(forest_low)) * (forest_right - forest_left)
}

dt[, `:=`(
  x_or_plot = map_or(OR),
  x_ci_low_plot = map_or(pmax(ci_low, forest_low)),
  x_ci_high_plot = map_or(pmin(ci_high, forest_high))
)]

axis_ticks <- c(0.5, 0.75, 1.0, 1.25, 1.5)
axis_dt <- data.table(
  tick = axis_ticks,
  x = map_or(axis_ticks),
  label = c("0.50", "0.75", "1.00", "1.25", "1.50")
)

x_tier <- 0.030
x_pair <- 0.140
x_or <- 0.690
x_fdr <- 0.805
x_heidi <- 0.875
x_pp4 <- 0.940
x_min <- 0
x_max <- 1

fig <- ggplot(dt, aes(y = y)) +
  geom_rect(
    data = bg_dt,
    aes(xmin = x_min, xmax = x_max, ymin = ymin, ymax = ymax, fill = tier),
    inherit.aes = FALSE,
    alpha = 0.82,
    color = NA
  ) +
  geom_hline(
    data = bg_dt,
    aes(yintercept = ymin),
    inherit.aes = FALSE,
    linewidth = 0.28,
    color = "white"
  ) +
  geom_segment(
    data = axis_dt,
    aes(x = x, xend = x, y = 0.40, yend = 0.55),
    inherit.aes = FALSE,
    linewidth = 0.32,
    color = "#222222"
  ) +
  geom_segment(
    aes(x = x_ci_low_plot, xend = x_ci_high_plot, y = y, yend = y, color = tier),
    linewidth = 0.78
  ) +
  geom_segment(
    x = map_or(1),
    xend = map_or(1),
    y = 0.55,
    yend = max(dt$y) + 0.47,
    linewidth = 0.38,
    linetype = "dashed",
    color = "#4A4A4A",
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(x = x_or_plot, fill = tier),
    shape = 21,
    size = 3.15,
    stroke = 0.58,
    color = "#171717"
  ) +
  geom_text(
    data = tier_labels,
    aes(x = x_tier, y = y, label = as.character(tier), color = tier),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.70,
    fontface = "bold"
  ) +
  geom_text(aes(x = x_pair, label = pair_label), hjust = 0, size = 2.72, color = "#111111") +
  geom_text(aes(x = x_or, label = or_label), hjust = 0, size = 2.38, color = "#111111") +
  geom_text(aes(x = x_fdr, label = fdr_label), hjust = 0, size = 2.38, color = "#111111") +
  geom_text(aes(x = x_heidi, label = heidi_label), hjust = 0, size = 2.38, color = "#2E6F40") +
  geom_text(aes(x = x_pp4, label = pp4_label), hjust = 0, size = 2.38, color = "#111111") +
  annotate("text", x = x_pair, y = max(dt$y) + 0.78, label = "Cell-gene pair", hjust = 0, size = 2.65, fontface = "bold") +
  annotate("text", x = x_or, y = max(dt$y) + 0.78, label = "OR (95% CI)", hjust = 0, size = 2.65, fontface = "bold") +
  annotate("text", x = x_fdr, y = max(dt$y) + 0.78, label = "FDR", hjust = 0, size = 2.65, fontface = "bold") +
  annotate("text", x = x_heidi, y = max(dt$y) + 0.78, label = "HEIDI", hjust = 0, size = 2.65, fontface = "bold") +
  annotate("text", x = x_pp4, y = max(dt$y) + 0.78, label = "PP4", hjust = 0, size = 2.65, fontface = "bold") +
  geom_text(
    data = axis_dt,
    aes(x = x, y = 0.18, label = label),
    inherit.aes = FALSE,
    size = 2.35,
    color = "#111111"
  ) +
  annotate("text", x = (forest_left + forest_right) / 2, y = -0.35,
           label = "OR for genetically predicted expression",
           size = 2.7, hjust = 0.5) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(expand = expansion(add = c(0.70, 1.05))) +
  scale_fill_manual(values = tier_bg, guide = "none") +
  scale_color_manual(values = tier_col, guide = "none") +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(8, 8, 12, 5)
  )

pdf_file <- file.path(out_dir, "figure3b_finngen_smr_tier1_3_forest_full.pdf")
png_file <- file.path(out_dir, "figure3b_finngen_smr_tier1_3_forest_full.png")
ggsave(pdf_file, fig, width = 6.8, height = 3.65, device = "pdf", useDingbats = FALSE)
ggsave(png_file, fig, width = 6.8, height = 3.65, dpi = 600, bg = "white")

cat("Wrote", png_file, "\n")
cat("Wrote", pdf_file, "\n")
