#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

project_dir <- Sys.getenv("PROJECT_DIR")
if (!nzchar(project_dir)) {
  project_dir <- Sys.getenv("PROJECT_DIR")
}
project_dir <- normalizePath(project_dir)

input_path <- file.path(
  project_dir,
  "results",
  "smr_analysis",
  "qc",
  "finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_validation.cellwise_fdr.tsv"
)
coloc_path <- file.path(
  project_dir,
  "results",
  "smr_analysis",
  "qc",
  "finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_coloc_abf_summary.tsv"
)
finngen_table1_path <- file.path(
  project_dir,
  "manuscript_results",
  "tables",
  "Table_figure3_smr_coloc_prioritized_tier1_3.tsv"
)
out_dir <- file.path(
  project_dir,
  "manuscript_results",
  "completed_results",
  "result4_meta_gwas_validation"
)
main_dir <- file.path(out_dir, "main")
source_dir <- file.path(out_dir, "source_tables")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_path)
coloc <- fread(coloc_path, select = c("cell", "Gene", "PP.H4.abf", "PP4_over_PP3_PP4", "coloc_abf_positive", "validation_and_coloc_positive"))
setnames(coloc, c("PP.H4.abf", "PP4_over_PP3_PP4"), c("meta_coloc_PP4", "meta_PP4_over_PP3_PP4"))
dt <- merge(dt, coloc, by = c("cell", "Gene"), all.x = TRUE, sort = FALSE)
finngen_table1 <- fread(finngen_table1_path, select = c("Cell-gene pair", "SMR OR (95% CI)"))
setnames(finngen_table1, c("Cell-gene pair", "SMR OR (95% CI)"), c("pair_display", "finngen_or_ci"))

parse_or_ci <- function(x) {
  m <- regexec("^([0-9.]+) \\(([0-9.]+)-([0-9.]+)\\)$", x)
  parts <- regmatches(x, m)
  out <- t(vapply(parts, function(p) {
    if (length(p) != 4) {
      return(c(or = NA_real_, lcl = NA_real_, ucl = NA_real_))
    }
    as.numeric(p[2:4])
  }, numeric(3)))
  out <- as.data.table(out)
  setnames(out, c("or", "lcl", "ucl"))
  out
}
finngen_ci <- parse_or_ci(finngen_table1$finngen_or_ci)
finngen_table1 <- cbind(finngen_table1, finngen_ci)

tier_order <- c("Tier 1", "Tier 2", "Tier 3")
pair_order <- c(
  "CD8_ET - FADS1",
  "NK - FADS1",
  "NK - MAFF",
  "CD8_ET - SLC44A2",
  "CD8_ET - TMEM258",
  "CD8_NC - TMEM258",
  "CD4_NC - TSPAN15",
  "CD4_ET - FADS1",
  "NK - SLC44A2"
)
dt[, pair := paste(gsub("_", " ", cell), "-", Gene)]
dt[, pair_key := paste(cell, "-", Gene)]
dt <- merge(dt, finngen_table1[, .(pair_display, finngen_or_ci, finngen_or_table = or, finngen_lcl_table = lcl, finngen_ucl_table = ucl)],
            by.x = "pair", by.y = "pair_display", all.x = TRUE, sort = FALSE)
dt[, tier := factor(tier, levels = tier_order)]
dt[, pair_key := factor(pair_key, levels = pair_order)]
setorder(dt, pair_key)
dt[, y := rev(seq_len(.N))]

dt[, `:=`(
  finngen_or_plot = fifelse(is.finite(finngen_or_table), finngen_or_table, OR_SMR_finngen),
  finngen_lcl = fifelse(is.finite(finngen_lcl_table), finngen_lcl_table, NA_real_),
  finngen_ucl = fifelse(is.finite(finngen_ucl_table), finngen_ucl_table, NA_real_),
  meta_lcl = exp(meta_ukb_approx_b_SMR - 1.96 * meta_ukb_approx_se_SMR),
  meta_ucl = exp(meta_ukb_approx_b_SMR + 1.96 * meta_ukb_approx_se_SMR)
)]

format_p <- function(x) {
  ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x))
}

dt[, `:=`(
  finngen_or_label = finngen_or_ci,
  meta_or_label = meta_ukb_approx_OR_95CI,
  meta_fdr_label = format_p(meta_ukb_approx_q_SMR_cell_BH),
  heidi_label = ifelse(meta_ukb_approx_p_HEIDI >= 0.01, "Pass", "Fail"),
  pp4_label = sprintf("%.3f", meta_coloc_PP4),
  overall_pass = validation_and_coloc_positive == TRUE
)]

write_cols <- c(
  "tier", "cell", "Gene", "pair", "OR_SMR_finngen", "finngen_or_label", "finngen_or_plot", "finngen_lcl", "finngen_ucl",
  "meta_ukb_approx_OR", "meta_or_label", "meta_lcl", "meta_ucl", "meta_ukb_approx_q_SMR_cell_BH",
  "meta_ukb_approx_p_HEIDI", "meta_coloc_PP4", "meta_ukb_approx_validation_pass_cellwise_fdr",
  "coloc_abf_positive", "validation_and_coloc_positive"
)
fwrite(dt[, ..write_cols], file.path(source_dir, "Table_result4_tier_stratified_replication_forest_values.tsv"), sep = "\t")

or_to_x <- function(or, center, half_width = 0.72, log_min = log(0.45), log_max = log(1.55)) {
  pmin(pmax(center + ((log(or) - log_min) / (log_max - log_min) - 0.5) * 2 * half_width, center - half_width), center + half_width)
}

left_center <- 3.10
right_center <- 5.12
half_width <- 0.72
log_min <- log(0.45)
log_max <- log(1.55)

plot_dt <- rbindlist(list(
  dt[, .(
    y,
    source = "FinnGen",
    or = finngen_or_plot,
    lcl = finngen_lcl,
    ucl = finngen_ucl,
    x = or_to_x(finngen_or_plot, left_center, half_width, log_min, log_max),
    xmin = or_to_x(finngen_lcl, left_center, half_width, log_min, log_max),
    xmax = or_to_x(finngen_ucl, left_center, half_width, log_min, log_max)
  )],
  dt[, .(
    y,
    source = "Independent meta-GWAS",
    or = meta_ukb_approx_OR,
    lcl = meta_lcl,
    ucl = meta_ucl,
    x = or_to_x(meta_ukb_approx_OR, right_center, half_width, log_min, log_max),
    xmin = or_to_x(meta_lcl, right_center, half_width, log_min, log_max),
    xmax = or_to_x(meta_ucl, right_center, half_width, log_min, log_max)
  )]
))
plot_dt[, source := factor(source, levels = c("FinnGen", "Independent meta-GWAS"))]

bg <- dt[, .(ymin = min(y) - 0.43, ymax = max(y) + 0.43), by = tier]
bg[, fill := c("#F5EDF3", "#FBF0DF", "#F7F2DE")[as.integer(tier)]]
tier_labels <- dt[, .(y = mean(y)), by = tier]

tick_or <- c(0.5, 1.0, 1.5)
ticks <- rbindlist(list(
  data.table(source = "FinnGen", or = tick_or, x = or_to_x(tick_or, left_center, half_width, log_min, log_max), center = left_center),
  data.table(source = "Independent meta-GWAS", or = tick_or, x = or_to_x(tick_or, right_center, half_width, log_min, log_max), center = right_center)
))
null_lines <- data.table(
  x = c(or_to_x(1, left_center, half_width, log_min, log_max), or_to_x(1, right_center, half_width, log_min, log_max)),
  xend = c(or_to_x(1, left_center, half_width, log_min, log_max), or_to_x(1, right_center, half_width, log_min, log_max)),
  y = 0.35,
  yend = 9.55
)

source_cols <- c("FinnGen" = "#8A1538", "Independent meta-GWAS" = "#2F5D7C")
p <- ggplot() +
  geom_rect(
    data = bg,
    aes(xmin = 0.05, xmax = 9.35, ymin = ymin, ymax = ymax),
    fill = bg$fill,
    colour = NA,
    alpha = 0.82
  ) +
  geom_segment(
    data = null_lines,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.35,
    linetype = "dashed",
    colour = "#6B7280"
  ) +
  geom_segment(
    data = plot_dt,
    aes(x = xmin, xend = xmax, y = y, yend = y, colour = source),
    linewidth = 0.82,
    lineend = "round"
  ) +
  geom_point(
    data = plot_dt[source == "FinnGen"],
    aes(x = x, y = y, fill = source),
    shape = 16,
    size = 3.45,
    colour = "#2B2B2B",
    stroke = 0.35
  ) +
  geom_point(
    data = plot_dt[source == "Independent meta-GWAS"],
    aes(x = x, y = y, fill = source),
    shape = 15,
    size = 3.45,
    colour = "#2B2B2B",
    stroke = 0.35
  ) +
  scale_colour_manual(values = source_cols, guide = "none") +
  scale_fill_manual(values = source_cols, guide = "none") +
  geom_text(data = tier_labels, aes(x = 0.18, y = y, label = as.character(tier)), hjust = 0, size = 3.35, fontface = "bold", colour = "#6B1238") +
  geom_text(data = dt, aes(x = 1.05, y = y, label = pair), hjust = 0, size = 3.45, colour = "#1F2933") +
  geom_text(data = dt, aes(x = 6.58, y = y, label = meta_fdr_label), hjust = 0.5, size = 3.1, colour = "#1F2933") +
  geom_text(data = dt, aes(x = 7.26, y = y, label = heidi_label), hjust = 0.5, size = 3.15, colour = ifelse(dt$meta_ukb_approx_p_HEIDI >= 0.01, "#2E6F40", "#9A3412")) +
  geom_text(data = dt, aes(x = 7.88, y = y, label = pp4_label), hjust = 0.5, size = 3.1, colour = ifelse(dt$coloc_abf_positive, "#1F2933", "#9A3412")) +
  geom_point(data = dt[overall_pass == TRUE], aes(x = 8.52, y = y), shape = 16, size = 3.05, colour = "#1F2933", stroke = 0.65) +
  geom_point(data = dt[overall_pass == FALSE], aes(x = 8.52, y = y), shape = 1, size = 3.05, colour = "#1F2933", stroke = 0.65) +
  geom_text(data = data.table(x = c(left_center, right_center, mean(c(6.58, 8.52))), y = 10.36, label = c("FinnGen discovery", "Independent meta-GWAS", "Replication evidence")),
            aes(x = x, y = y, label = label), hjust = 0.5, fontface = "bold", size = 3.45, colour = "#111827") +
  geom_text(data = data.table(x = 1.05, y = 9.90, label = "Cell-gene pair"),
            aes(x = x, y = y, label = label), hjust = 0, size = 3.2, colour = "#4B5563") +
  geom_text(data = data.table(x = c(left_center, right_center, 6.58, 7.26, 7.88, 8.52), y = 9.90, label = c("OR (95% CI)", "OR (95% CI)", "FDR", "HEIDI", "PP4", "Overall")),
            aes(x = x, y = y, label = label), hjust = 0.5, size = 3.2, colour = "#4B5563") +
  geom_segment(data = ticks, aes(x = x, xend = x, y = 0.12, yend = 0.24), linewidth = 0.35, colour = "#111827") +
  geom_segment(data = data.table(center = c(left_center, right_center)), aes(x = or_to_x(0.5, center, half_width, log_min, log_max), xend = or_to_x(1.5, center, half_width, log_min, log_max), y = 0.12, yend = 0.12), linewidth = 0.35, colour = "#111827") +
  geom_text(data = ticks, aes(x = x, y = -0.10, label = sprintf("%.1f", or)), size = 2.9, colour = "#111827") +
  annotate("text", x = mean(c(left_center, right_center)), y = -0.58, label = "OR for genetically predicted expression", size = 3.45, colour = "#111827") +
  coord_cartesian(xlim = c(0.05, 9.40), ylim = c(-0.78, 10.62), clip = "off") +
  theme_void(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(14, 14, 18, 14)
  )

png_path <- file.path(main_dir, "Figure5_tier_stratified_replication_forest.png")
pdf_path <- file.path(main_dir, "Figure5_tier_stratified_replication_forest.pdf")
ggsave(png_path, p, width = 10.2, height = 5.7, units = "in", dpi = 600, bg = "white")
ggsave(pdf_path, p, width = 10.2, height = 5.7, units = "in", device = "pdf", bg = "white")

message("Wrote ", png_path)
message("Wrote ", pdf_path)
