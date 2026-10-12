#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
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
summary_path <- file.path(project_dir, "results/smr_analysis/qc/finngen_r13_discovery_tier1.coloc_abf_summary.tsv")
abf_root <- file.path(project_dir, "results/coloc/finngen_r13_discovery_abf")
ld_root <- file.path(project_dir, "results/coloc/finngen_r13_discovery_susie")
out_dir <- file.path(project_dir, "manuscript_results/figures/figure4_tier1_coloc")
table_dir <- file.path(project_dir, "manuscript_results/tables")
rsid_path <- file.path(project_dir, "manuscript_results/tables/top_snp_rsid_mapping_hg19.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

neglog10 <- function(p) -log10(pmax(as.numeric(p), .Machine$double.xmin))

locate_coloc_dir <- function(root, row) {
  candidates <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  pattern <- sprintf("_%s_%s_chr%d$", row$cell, row$Gene, as.integer(row$chr))
  hit <- candidates[grepl(pattern, basename(candidates))]
  if (length(hit) != 1) {
    stop("Cannot locate directory for ", row$cell, " ", row$Gene, " chr", row$chr, " under ", root)
  }
  hit
}

read_ld_r2 <- function(ld_dir, top_snp) {
  vars_path <- file.path(ld_dir, "onek1k_ld.unphased.vcor1.vars")
  mat_path <- file.path(ld_dir, "onek1k_ld.unphased.vcor1")
  if (!file.exists(vars_path) || !file.exists(mat_path)) {
    return(NULL)
  }
  vars <- fread(vars_path, header = FALSE)$V1
  top_idx <- match(top_snp, vars)
  if (is.na(top_idx)) {
    return(data.table(SNP = vars, ld_r2 = NA_real_))
  }
  con <- file(mat_path, open = "r")
  on.exit(close(con), add = TRUE)
  line <- NULL
  for (i in seq_len(top_idx)) {
    line <- readLines(con, n = 1)
  }
  r <- as.numeric(strsplit(line, "\\s+")[[1]])
  data.table(SNP = vars, ld_r2 = pmin(pmax(r^2, 0), 1))
}

ld_bin_palette <- c(
  "0-0.2" = "#253494",
  "0.2-0.4" = "#74CFE3",
  "0.4-0.6" = "#238B45",
  "0.6-0.8" = "#FDBF2D",
  "0.8-1" = "#E31A1C"
)

bin_ld <- function(x) {
  cut(
    x,
    breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
    labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
    right = FALSE
  )
}

thin_points <- function(x, max_n = 4200L) {
  if (nrow(x) <= max_n) return(x)
  keep <- unique(c(
    x[order(p_gwas)][seq_len(min(300L, .N)), SNP],
    x[order(p_eqtl)][seq_len(min(300L, .N)), SNP],
    x[order(-ld_r2)][seq_len(min(300L, .N)), SNP],
    x[is_top == TRUE, SNP]
  ))
  top <- x[SNP %in% keep]
  rest <- x[!SNP %in% keep]
  if (nrow(rest) > max_n - nrow(top)) {
    set.seed(20260813)
    rest <- rest[sample(.N, max_n - nrow(top))]
  }
  rbind(top, rest, fill = TRUE)
}

make_panel <- function(row) {
  abf_dir <- locate_coloc_dir(abf_root, row)
  ld_dir <- locate_coloc_dir(ld_root, row)
  x <- fread(file.path(abf_dir, "coloc_input_merged.tsv"))
  ld <- read_ld_r2(ld_dir, row$topSNP)
  if (!is.null(ld)) {
    x <- merge(x, ld, by = "SNP", all.x = TRUE)
  } else {
    x[, ld_r2 := NA_real_]
  }
  x[, `:=`(
    bp = as.numeric(BP_gwas),
    pos_mb = as.numeric(BP_gwas) / 1e6,
    gwas_logp = neglog10(p_gwas),
    eqtl_logp = neglog10(p_eqtl),
    is_top = SNP == row$topSNP,
    ld_bin = bin_ld(ld_r2)
  )]
  x[is.na(ld_bin), ld_bin := "0-0.2"]
  x <- x[is.finite(pos_mb) & is.finite(gwas_logp) & is.finite(eqtl_logp)]
  xp <- thin_points(x)
  top <- x[SNP == row$topSNP]
  top_bp_numeric <- if (nrow(top)) top$bp[1] else as.numeric(strsplit(row$topSNP, ":", fixed = TRUE)[[1]][2])
  top_pos <- top_bp_numeric / 1e6
  top_label <- if (!is.na(row$rsid) && nzchar(row$rsid)) row$rsid else paste0("chr", as.integer(row$chr), ":", as.integer(round(top_bp_numeric)))
  xlim <- range(x$pos_mb, na.rm = TRUE)

  base_theme <- theme_classic(base_size = 10.5, base_family = "sans") +
    theme(
      axis.text = element_text(size = 9.2, color = "black"),
      axis.title = element_text(size = 10.5, color = "black"),
      plot.margin = margin(2, 3, 2, 3),
      legend.position = "none"
    )

  shared_scales <- list(
    scale_color_manual(values = ld_bin_palette, drop = FALSE),
    coord_cartesian(xlim = xlim)
  )

  p_left <- ggplot(xp, aes(gwas_logp, eqtl_logp)) +
    geom_point(aes(color = ld_bin), size = 1.25, alpha = 0.84) +
    geom_point(data = x[is_top == TRUE], shape = 23, fill = "#C51B8A", color = "#4A004A", size = 3.25, stroke = 0.65) +
    geom_text(
      data = x[is_top == TRUE],
      aes(label = top_label),
      size = 3.2,
      hjust = 0.5,
      vjust = -0.55
    ) +
    scale_color_manual(values = ld_bin_palette, drop = FALSE) +
    labs(x = expression("GWAS "*-log[10](P)), y = expression("eQTL "*-log[10](P))) +
    base_theme +
    theme(
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4))
    )

  p_gwas <- ggplot(xp, aes(pos_mb, gwas_logp)) +
    geom_point(aes(color = ld_bin), size = 1.25, alpha = 0.84) +
    geom_point(data = x[is_top == TRUE], shape = 23, fill = "#C51B8A", color = "#4A004A", size = 3.15, stroke = 0.65) +
    geom_vline(xintercept = top_pos, color = "#4A004A", linetype = 2, linewidth = 0.34) +
    shared_scales +
    labs(x = NULL, y = expression("GWAS "*-log[10](P))) +
    base_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_eqtl <- ggplot(xp, aes(pos_mb, eqtl_logp)) +
    geom_point(aes(color = ld_bin), size = 1.25, alpha = 0.84) +
    geom_point(data = x[is_top == TRUE], shape = 23, fill = "#C51B8A", color = "#4A004A", size = 3.15, stroke = 0.65) +
    geom_vline(xintercept = top_pos, color = "#4A004A", linetype = 2, linewidth = 0.34) +
    shared_scales +
    labs(x = sprintf("chr%d position (Mb)", as.integer(row$chr)), y = expression("eQTL "*-log[10](P))) +
    base_theme +
    theme(axis.title.x = element_text(margin = margin(t = 4)))

  title <- sprintf("%s - %s   PP4 = %.3f", gsub("_", " ", row$cell), row$Gene, row$finngen_PP.H4.abf)
  right <- arrangeGrob(p_gwas, p_eqtl, ncol = 1, heights = c(1, 1))
  body <- arrangeGrob(p_left, right, ncol = 2, widths = c(1.15, 1))
  grob <- arrangeGrob(
    textGrob(title, x = unit(0.02, "npc"), just = "left", gp = gpar(fontsize = 11.5, fontface = "bold")),
    body,
    ncol = 1,
    heights = c(0.11, 1)
  )

  plot_data <- x[, .(
    cell = row$cell, Gene = row$Gene, chr = row$chr, topSNP = row$topSNP, topSNP_rsid = row$rsid,
    SNP, bp, pos_mb, gwas_logp, eqtl_logp, p_gwas, p_eqtl, ld_r2, ld_bin, is_top
  )]
  list(grob = grob, data = plot_data)
}

summary <- fread(summary_path)
rsid_map <- fread(rsid_path)[, .(topSNP, rsid, CHR_hg19, BP_hg19)]
summary <- merge(summary, rsid_map, by = "topSNP", all.x = TRUE)
if (anyNA(summary$rsid)) {
  stop("Missing rsID mapping for Figure 4 top SNPs")
}
summary <- summary[order(match(paste(cell, Gene), c("CD8_ET SLC44A2", "CD8_ET FADS1", "NK FADS1", "NK MAFF")))]
panels <- lapply(seq_len(nrow(summary)), function(i) make_panel(summary[i]))

plot_data <- rbindlist(lapply(panels, `[[`, "data"), fill = TRUE)
fwrite(plot_data, file.path(table_dir, "figure4_tier1_coloc_reference_style_plot_data.tsv"), sep = "\t")

manifest <- data.table()
for (i in seq_len(nrow(summary))) {
  row <- summary[i]
  prefix <- sprintf("%s_%s_chr%d.reference_style_coloc", row$cell, row$Gene, as.integer(row$chr))
  png_file <- file.path(out_dir, paste0(prefix, ".png"))
  pdf_file <- file.path(out_dir, paste0(prefix, ".pdf"))
  ggsave(png_file, panels[[i]]$grob, width = 5.0, height = 4.6, dpi = 600, bg = "white")
  ggsave(pdf_file, panels[[i]]$grob, width = 5.0, height = 4.6, device = "pdf", useDingbats = FALSE, bg = "white")
  manifest <- rbind(
    manifest,
    data.table(cell = row$cell, Gene = row$Gene, chr = row$chr, topSNP = row$topSNP,
               topSNP_rsid = row$rsid, PP.H4.abf = row$finngen_PP.H4.abf, png = png_file, pdf = pdf_file)
  )
}

legend_grob <- {
  ld_df <- data.table(ld_bin = factor(names(ld_bin_palette), levels = names(ld_bin_palette)), x = seq_along(ld_bin_palette), y = 1)
  lg <- ggplot(ld_df, aes(x, y, fill = ld_bin)) +
    geom_tile(width = 0.9, height = 0.5, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = ld_bin_palette, name = expression(LD~r^2), drop = FALSE) +
    guides(fill = guide_legend(nrow = 1, title.position = "left")) +
    theme_void(base_size = 10.5) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10.5),
      legend.text = element_text(size = 9.2),
      legend.key.width = unit(0.55, "cm"),
      legend.key.height = unit(0.30, "cm")
    )
  g <- ggplotGrob(lg)
  g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}

combined <- arrangeGrob(
  arrangeGrob(panels[[1]]$grob, panels[[2]]$grob, panels[[3]]$grob, panels[[4]]$grob, ncol = 2),
  legend_grob,
  ncol = 1,
  heights = c(1, 0.05)
)

combined_png <- file.path(out_dir, "figure4_tier1_coloc_reference_style_combined.png")
combined_pdf <- file.path(out_dir, "figure4_tier1_coloc_reference_style_combined.pdf")
ggsave(combined_png, combined, width = 10.2, height = 9.6, dpi = 600, bg = "white")
ggsave(combined_pdf, combined, width = 10.2, height = 9.6, device = "pdf", useDingbats = FALSE, bg = "white")

manifest_path <- file.path(out_dir, "figure4_tier1_coloc_reference_style_manifest.tsv")
fwrite(manifest, manifest_path, sep = "\t")

cat("Wrote", combined_png, "\n")
cat("Wrote", combined_pdf, "\n")
cat("Wrote", manifest_path, "\n")
cat("Wrote", file.path(table_dir, "figure4_tier1_coloc_reference_style_plot_data.tsv"), "\n")
