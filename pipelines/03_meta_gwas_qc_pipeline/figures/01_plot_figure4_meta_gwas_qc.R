#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(png)
  library(scales)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
get_project_dir <- function() {
  env_project <- Sys.getenv("PROJECT_DIR")
  if (nzchar(env_project)) {
    return(normalizePath(env_project))
  }
  normalizePath(file.path(dirname(script_path), "..", "..", ".."))
}

project_dir <- get_project_dir()

prefix <- "vte_discovery_finngen_independent_ivw_ukb_approx_logor"
meta_path <- file.path(project_dir, "gwas", "meta", paste0(prefix, ".primary_ncohort3plus.tsv.gz"))
run_dir <- Sys.getenv("STANDARD_LDSC_RUN_DIR", unset = file.path(project_dir, "gwas/ldsc_standard"))
ldsc_matched_path <- file.path(run_dir, "prepared/vte_meta_hm3_from_raw.ldscore_matched.tsv.gz")
ldsc_qc_path <- file.path(project_dir, "manuscript_results/completed_results/result3_meta_gwas_qc/supplementary/Table_S8_ldsc_calibration_summary.tsv")
out_dir <- file.path(
  project_dir,
  "manuscript_results",
  "completed_results",
  "result3_meta_gwas_qc"
)
main_dir <- file.path(out_dir, "main")
source_dir <- file.path(out_dir, "source_tables")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

unicode_pdf_device <- function(filename, width, height, bg, ...) {
  if (Sys.info()[["sysname"]] == "Darwin") {
    grDevices::quartz(type = "pdf", file = filename, width = width, height = height, bg = bg)
  } else {
    grDevices::cairo_pdf(filename = filename, width = width, height = height, bg = bg, ...)
  }
}

panel_col <- "#8A1538"
accent_col <- "#B55A30"
muted_col <- "#6B7280"
light_bg <- "#F6F6F6"
line_grey <- "#7C8490"
dark_grey <- "#30343B"
light_grey <- "#D9DDE3"

theme_fig4 <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.3),
      plot.title = element_blank(),
      plot.margin = margin(7, 7, 7, 7),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1)
    )
}

panel_label <- function(p, label) {
  p + annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = label,
    hjust = -0.35,
    vjust = 1.25,
    fontface = "bold",
    size = 7.1
  )
}

make_manhattan_dt <- function(dt, cap_y = 50, max_background = 300000L) {
  set.seed(20260815)
  chr_lengths <- dt[, .(chr_max = as.numeric(max(BP, na.rm = TRUE))), by = CHR][order(CHR)]
  chr_lengths[, offset := shift(cumsum(chr_max), fill = 0)]
  centers <- chr_lengths[, .(CHR, center = offset + chr_max / 2)]
  x <- merge(dt, chr_lengths[, .(CHR, offset)], by = "CHR", sort = FALSE)
  x[, pos_cum := as.numeric(BP) + offset]
  x[, logp := -log10(pmax(p, .Machine$double.xmin))]
  x[, logp_capped := pmin(logp, cap_y)]
  keep <- x[p <= 1e-5]
  bg <- x[p > 1e-5]
  if (nrow(bg) > max_background) {
    bg <- bg[sample(.N, max_background)]
  }
  plot_dt <- rbindlist(list(keep, bg), use.names = TRUE)
  plot_dt[, chr_group := ifelse(CHR %% 2 == 1, "odd", "even")]
  list(plot_dt = plot_dt, centers = centers)
}

make_qq_dt <- function(p, max_points = 250000L) {
  p <- sort(p[is.finite(p) & p > 0 & p <= 1])
  n <- length(p)
  observed <- -log10(p)
  expected <- -log10((seq_len(n) - 0.5) / n)
  idx <- unique(c(
    seq_len(min(n, 20000L)),
    round(seq(20001L, n, length.out = min(max_points, n)))
  ))
  ci_low_p <- qbeta(0.975, idx, n - idx + 1)
  ci_high_p <- qbeta(0.025, idx, n - idx + 1)
  data.table(
    expected = expected[idx],
    observed = observed[idx],
    lower = -log10(ci_low_p),
    upper = -log10(ci_high_p)
  )
}

message("Reading meta-GWAS...")
meta <- fread(meta_path, select = c("CHR", "BP", "p", "het_p"))
meta[, `:=`(
  CHR = as.integer(CHR),
  BP = as.integer(BP),
  p = as.numeric(p),
  het_p = as.numeric(het_p)
)]
meta <- meta[CHR %between% c(1L, 22L) & is.finite(BP) & is.finite(p) & p >= 0 & p <= 1]

message("Reading LDSC matched HapMap3 SNPs...")
ldsc <- fread(ldsc_matched_path, select = c("SNP", "P"))
munged <- fread(file.path(run_dir, "munged/vte_meta_standard.sumstats.gz"), select = c("SNP", "Z"))
munged <- munged[is.finite(Z)]
ldsc_primary <- merge(ldsc, munged[Z^2 < 80], by = "SNP")
stopifnot(nrow(ldsc_primary) == nrow(munged[Z^2 < 80]))
s8 <- fread(ldsc_qc_path)
standard <- s8[`Analysis set` == "Standard LDSC regression set"]
stopifnot(nrow(standard) == 1L, nrow(munged) == standard$SNPs)
primary_qc <- data.table(n_snps=standard$SNPs, mean_chi2=standard$`Mean chi-square`,
 lambda_gc=standard$`Lambda GC`, ldsc_intercept=standard$`LDSC intercept`,
 ldsc_intercept_se=standard$`LDSC intercept SE`, attenuation_ratio=standard$`Attenuation ratio`,
 attenuation_ratio_se=standard$`Attenuation ratio SE`)

meta_n <- 6891026L
harmonized_n <- 5783976L
gws_n <- 3428L
gws_het_pass_n <- 2626L
hm3_n <- as.integer(primary_qc$n_snps)
harmonized_pct <- harmonized_n / meta_n
het_pass_pct <- gws_het_pass_n / gws_n

fwrite(
  data.table(
    metric = c(
      "meta_gwas_ncohorts_ge3_variants",
      "oneK1K_ld_harmonized_variants",
      "harmonized_fraction",
      "genome_wide_significant_variants",
      "genome_wide_significant_hetp_gt_0.01",
      "gws_hetp_gt_0.01_fraction",
      "hapmap3_eur_ldsc_qc_snps",
      "hapmap3_allele_mismatches",
      "mean_chi2",
      "lambda_gc",
      "ldsc_intercept",
      "ldsc_intercept_se",
      "attenuation_ratio"
    ),
    value = c(
      meta_n,
      harmonized_n,
      sprintf("%.3f", harmonized_pct),
      gws_n,
      gws_het_pass_n,
      sprintf("%.3f", het_pass_pct),
      hm3_n,
      0,
      sprintf("%.3f", primary_qc$mean_chi2),
      sprintf("%.3f", primary_qc$lambda_gc),
      sprintf("%.3f", primary_qc$ldsc_intercept),
      sprintf("%.3f", primary_qc$ldsc_intercept_se),
      sprintf("%.3f", primary_qc$attenuation_ratio)
    )
  ),
  file.path(source_dir, "Table_result3_figure4_plot_values.tsv"),
  sep = "\t"
)

cohorts <- data.table(
  cohort = c("deCODE DK", "deCODE Iceland", "deCODE Intermountain", "UK Biobank"),
  sample = c(
    "18,569 cases / 213,503 controls",
    "11,700 cases / 311,745 controls",
    "6,981 cases / 51,252 controls",
    "4,620 cases / 356,574 controls"
  ),
  x = 0.055,
  y = c(0.78, 0.61, 0.44, 0.27)
)
p_a <- ggplot() +
  geom_rect(aes(xmin = 0.02, xmax = 0.47, ymin = 0.09, ymax = 0.90), fill = light_bg, colour = NA) +
  geom_rect(
    data = cohorts,
    aes(xmin = x, xmax = x + 0.37, ymin = y - 0.062, ymax = y + 0.062),
    fill = "white",
    colour = "#9CA3AF",
    linewidth = 0.34
  ) +
  geom_text(
    data = cohorts,
    aes(x = x + 0.018, y = y + 0.022, label = cohort),
    hjust = 0,
    fontface = "bold",
    size = 3.75
  ) +
  geom_text(
    data = cohorts,
    aes(x = x + 0.018, y = y - 0.025, label = sample),
    hjust = 0,
    size = 3.10,
    colour = dark_grey
  ) +
  geom_segment(aes(x = 0.49, xend = 0.55, y = 0.84, yend = 0.84), linewidth = 0.42, colour = line_grey) +
  geom_segment(aes(x = 0.49, xend = 0.55, y = 0.16, yend = 0.16), linewidth = 0.42, colour = line_grey) +
  geom_segment(aes(x = 0.55, xend = 0.55, y = 0.84, yend = 0.16), linewidth = 0.42, colour = line_grey) +
  geom_segment(
    aes(x = 0.55, xend = 0.63, y = 0.50, yend = 0.50),
    arrow = arrow(length = unit(0.10, "inches"), type = "closed"),
    linewidth = 0.42,
    colour = line_grey
  ) +
  geom_rect(aes(xmin = 0.64, xmax = 0.98, ymin = 0.40, ymax = 0.70), fill = panel_col, colour = NA) +
  geom_text(aes(x = 0.81, y = 0.62, label = "VTE meta-GWAS"), colour = "white", fontface = "bold", size = 4.25) +
  geom_text(aes(x = 0.81, y = 0.535, label = "41,870 cases / 933,074 controls"), colour = "white", size = 3.05) +
  geom_text(aes(x = 0.81, y = 0.455, label = "FinnGen-independent"), colour = "white", fontface = "bold", size = 3.05) +
  geom_text(aes(x = 0.81, y = 0.29, label = "FinnGen excluded"), colour = panel_col, fontface = "bold", size = 3.25) +
  geom_text(aes(x = 0.81, y = 0.22, label = "Independent of discovery GWAS"), colour = dark_grey, size = 2.75) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_size = 10.5) +
  theme(plot.margin = margin(6, 8, 3, 12)) +
  annotate("text", x = -0.035, y = 1.00, label = "A", hjust = 0, vjust = 1, fontface = "bold", size = 7.1)

p_b <- ggplot() +
  geom_rect(aes(xmin = 0.20, xmax = 0.80, ymin = 0.70, ymax = 0.91), fill = panel_col, colour = NA) +
  geom_text(aes(x = 0.50, y = 0.845, label = "Primary meta-GWAS"), colour = "white", fontface = "bold", size = 3.75) +
  geom_text(aes(x = 0.50, y = 0.785, label = "n cohorts ≥ 3"), colour = "white", size = 2.85) +
  geom_text(aes(x = 0.50, y = 0.735, label = "6,891,026 variants"), colour = "white", size = 3.00) +
  geom_segment(aes(x = 0.50, xend = 0.50, y = 0.68, yend = 0.62), linewidth = 0.38, colour = line_grey) +
  geom_segment(aes(x = 0.25, xend = 0.75, y = 0.62, yend = 0.62), linewidth = 0.38, colour = line_grey) +
  geom_segment(
    aes(x = 0.25, xend = 0.25, y = 0.62, yend = 0.51),
    arrow = arrow(length = unit(0.085, "inches"), type = "closed"),
    linewidth = 0.38,
    colour = line_grey
  ) +
  geom_segment(
    aes(x = 0.75, xend = 0.75, y = 0.62, yend = 0.51),
    arrow = arrow(length = unit(0.085, "inches"), type = "closed"),
    linewidth = 0.38,
    colour = line_grey
  ) +
  geom_rect(aes(xmin = 0.04, xmax = 0.46, ymin = 0.25, ymax = 0.50), fill = "#ECEFF3", colour = "#BFC4CC", linewidth = 0.34) +
  geom_text(aes(x = 0.25, y = 0.425, label = "Matched to OneK1K\nLD reference"), colour = dark_grey, fontface = "bold", size = 3.00, lineheight = 0.95) +
  geom_text(aes(x = 0.25, y = 0.335, label = "5,783,976 variants"), colour = dark_grey, size = 2.85) +
  geom_text(aes(x = 0.25, y = 0.275, label = "83.9% retained"), colour = panel_col, fontface = "bold", size = 2.95) +
  geom_rect(aes(xmin = 0.56, xmax = 0.95, ymin = 0.25, ymax = 0.50), fill = light_bg, colour = "#BFC4CC", linewidth = 0.34) +
  geom_text(aes(x = 0.755, y = 0.435, label = "Genome-wide\nsignificant"), fontface = "bold", size = 2.95, lineheight = 0.95) +
  geom_text(aes(x = 0.755, y = 0.350, label = "3,428 variants"), size = 2.85) +
  geom_segment(
    aes(x = 0.755, xend = 0.755, y = 0.315, yend = 0.245),
    arrow = arrow(length = unit(0.075, "inches"), type = "closed"),
    linewidth = 0.34,
    colour = line_grey
  ) +
  geom_rect(aes(xmin = 0.56, xmax = 0.95, ymin = 0.06, ymax = 0.22), fill = "white", colour = "#BFC4CC", linewidth = 0.34) +
  geom_text(aes(x = 0.755, y = 0.170, label = "HetP > 0.01"), size = 2.75, fontface = "bold") +
  geom_text(aes(x = 0.755, y = 0.120, label = "2,626 variants retained"), size = 2.60) +
  geom_text(aes(x = 0.755, y = 0.075, label = "76.6% retained"), colour = panel_col, fontface = "bold", size = 2.80) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_size = 10.5) +
  theme(plot.margin = margin(6, 8, 3, 12)) +
  annotate("text", x = -0.035, y = 1.00, label = "B", hjust = 0, vjust = 1, fontface = "bold", size = 7.1)

man <- make_manhattan_dt(meta, cap_y = 50)
plot_dt <- man$plot_dt
centers <- man$centers
plot_dt[, signal := ifelse(p <= 5e-8, "GWS", ifelse(CHR %% 2 == 1, "odd", "even"))]
threshold_label_x <- min(plot_dt$pos_cum, na.rm = TRUE) + 0.018 * diff(range(plot_dt$pos_cum, na.rm = TRUE))
p_c_no_label <- ggplot(plot_dt, aes(x = pos_cum, y = logp_capped, colour = chr_group)) +
  geom_point(data = plot_dt[signal != "GWS"], aes(colour = signal), size = 0.30, alpha = 0.70, stroke = 0) +
  geom_point(data = plot_dt[signal == "GWS"], aes(colour = signal), size = 0.42, alpha = 0.88, stroke = 0) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", linewidth = 0.35, colour = "#303030") +
  scale_colour_manual(values = c(odd = "#5F6670", even = "#C8CDD4", GWS = panel_col), guide = "none") +
  scale_x_continuous(
    breaks = centers[, center],
    labels = centers[, CHR],
    expand = expansion(mult = c(0.004, 0.004))
  ) +
  scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.025))) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  annotate("text", x = Inf, y = 48.1, label = "3,428 GWS variants", hjust = 1.02, vjust = 1, size = 3.05) +
  annotate("text", x = threshold_label_x, y = -log10(5e-8) + 1.35, label = "P == 5 %*% 10^{-8}", parse = TRUE, hjust = 0, vjust = 0, size = 2.45, colour = "#4B5563") +
  theme_fig4(10.5) +
  theme(axis.text.x = element_text(size = 6.3, colour = "black"))
p_c <- p_c_no_label +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.42, vjust = 1.30, fontface = "bold", size = 7.1)

qq_dt <- make_qq_dt(ldsc_primary$P)
max_xy <- ceiling(max(c(qq_dt$expected, qq_dt$observed), na.rm = TRUE))
p_d_no_label <- ggplot(qq_dt, aes(x = expected, y = observed)) +
  geom_point(size = 0.42, alpha = 0.55, colour = dark_grey, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, colour = panel_col) +
  coord_cartesian(xlim = c(0, max(7, max(qq_dt$expected))), ylim = c(0, max_xy)) +
  annotate(
    "text",
    x = 0.35,
    y = max_xy - 0.30,
    hjust = 0,
    vjust = 1,
    label = sprintf(
      "HapMap3/EUR SNPs used for LDSC\nQQ display: Chi-square < 80 (n = %s)\nStandard LDSC: all n = %s\nλGC = %.3f\nLDSC intercept = %.3f ± %.3f\nRatio = %.3f ± %.3f",
      format(nrow(ldsc_primary), big.mark = ",", scientific = FALSE),
      format(hm3_n, big.mark = ",", scientific = FALSE),
      primary_qc$lambda_gc,
      primary_qc$ldsc_intercept,
      primary_qc$ldsc_intercept_se,
      primary_qc$attenuation_ratio,
      primary_qc$attenuation_ratio_se
    ),
    colour = dark_grey,
    size = 2.65
  ) +
  labs(x = expression(Expected~~-log[10](P)), y = expression(Observed~~-log[10](P))) +
  theme_fig4(10.5)
p_d <- p_d_no_label +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.52, vjust = 1.30, fontface = "bold", size = 7.1)

top_row <- arrangeGrob(p_a, p_b, ncol = 2, widths = c(1.04, 0.96))
bottom_row <- arrangeGrob(p_c, p_d, ncol = 2, widths = c(1.45, 0.95))
figure4 <- arrangeGrob(top_row, bottom_row, ncol = 1, heights = c(0.82, 1.68))

png_path <- file.path(main_dir, "Figure4_meta_gwas_construction_qc.png")
pdf_path <- file.path(main_dir, "Figure4_meta_gwas_construction_qc.pdf")
ab_vector_pdf_path <- file.path(main_dir, "Figure4AB_meta_gwas_construction_filtering_vector.pdf")
ggsave(png_path, figure4, width = 10.8, height = 8.0, units = "in", dpi = 350, bg = "white")
ggsave(pdf_path, figure4, width = 10.8, height = 8.0, units = "in", device = unicode_pdf_device, bg = "white")
ggsave(ab_vector_pdf_path, top_row, width = 10.8, height = 2.65, units = "in", device = unicode_pdf_device, bg = "white")

single_panel_dir <- file.path(main_dir, "single_panels_no_labels")
dir.create(single_panel_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(single_panel_dir, "Figure4C_manhattan_no_panel_label.png"), p_c_no_label, width = 6.85, height = 4.05, units = "in", dpi = 600, bg = "white")
ggsave(file.path(single_panel_dir, "Figure4C_manhattan_no_panel_label.pdf"), p_c_no_label, width = 6.85, height = 4.05, units = "in", device = unicode_pdf_device, bg = "white")
ggsave(file.path(single_panel_dir, "Figure4D_qq_ldsc_no_panel_label.png"), p_d_no_label, width = 4.35, height = 4.05, units = "in", dpi = 600, bg = "white")
ggsave(file.path(single_panel_dir, "Figure4D_qq_ldsc_no_panel_label.pdf"), p_d_no_label, width = 4.35, height = 4.05, units = "in", device = unicode_pdf_device, bg = "white")

message("Wrote ", png_path)
message("Wrote ", pdf_path)
message("Wrote ", ab_vector_pdf_path)

hybrid_dir <- file.path(main_dir, "hybrid_vector_panels")
dir.create(hybrid_dir, recursive = TRUE, showWarnings = FALSE)
c_raster_path <- file.path(hybrid_dir, "Figure4C_manhattan_raster.png")
d_raster_path <- file.path(hybrid_dir, "Figure4D_qq_ldsc_raster.png")

ggsave(c_raster_path, p_c, width = 6.85, height = 4.05, units = "in", dpi = 600, bg = "white")
ggsave(d_raster_path, p_d, width = 4.35, height = 4.05, units = "in", dpi = 600, bg = "white")

panel_c_raster <- rasterGrob(readPNG(c_raster_path), interpolate = TRUE)
panel_d_raster <- rasterGrob(readPNG(d_raster_path), interpolate = TRUE)
bottom_row_hybrid <- arrangeGrob(
  panel_c_raster,
  panel_d_raster,
  ncol = 2,
  widths = c(1.45, 0.95)
)
figure4_hybrid <- arrangeGrob(
  top_row,
  bottom_row_hybrid,
  ncol = 1,
  heights = c(0.82, 1.68)
)

hybrid_pdf_path <- file.path(main_dir, "Figure4_meta_gwas_construction_qc_hybrid_vector_AB_raster_CD.pdf")
hybrid_svg_path <- file.path(main_dir, "Figure4_meta_gwas_construction_qc_hybrid_vector_AB_raster_CD.svg")
ggsave(hybrid_pdf_path, figure4_hybrid, width = 10.8, height = 8.0, units = "in", device = unicode_pdf_device, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(hybrid_svg_path, figure4_hybrid, width = 10.8, height = 8.0, units = "in", device = "svg", bg = "white")
  message("Wrote ", hybrid_svg_path)
} else {
  message("Skipped SVG output because the svglite package is not installed.")
}

message("Wrote ", hybrid_pdf_path)
