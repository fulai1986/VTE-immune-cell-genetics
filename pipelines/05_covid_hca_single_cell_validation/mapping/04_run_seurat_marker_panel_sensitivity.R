#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("PROJECT_DIR"))) stop("Set PROJECT_DIR to the analysis workspace before running this script.")

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
})

base_dir <- file.path(Sys.getenv("PROJECT_DIR"), "covid_validation/hca_onek1k_mapping")
input_dir <- file.path(base_dir, "seurat_label_transfer")
out_dir <- file.path(base_dir, "seurat_label_transfer", "results_marker_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260827)

read_expr <- function(mtx_path, features_path, metadata_path) {
  mat <- readMM(mtx_path)
  features <- fread(features_path)
  meta <- fread(metadata_path)
  rownames(mat) <- features$gene
  colnames(mat) <- meta$cell_id
  list(mat = as(mat, "dgCMatrix"), meta = as.data.frame(meta))
}

ref_in <- read_expr(
  file.path(input_dir, "onek1k_reference_logexpr.mtx"),
  file.path(input_dir, "shared_genes.tsv"),
  file.path(input_dir, "onek1k_reference_metadata.tsv")
)
qry_in <- read_expr(
  file.path(input_dir, "hca_healthy_query_logexpr.mtx"),
  file.path(input_dir, "shared_genes.tsv"),
  file.path(input_dir, "hca_healthy_query_metadata.tsv")
)

marker_features <- fread(file.path(base_dir, "results", "hca_onek1k_mapping_feature_panel.tsv"))$feature_gene
features <- intersect(marker_features, rownames(ref_in$mat))
features <- intersect(features, rownames(qry_in$mat))
message("Marker features used: ", length(features))

rownames(ref_in$meta) <- ref_in$meta$cell_id
rownames(qry_in$meta) <- qry_in$meta$cell_id

ref <- CreateSeuratObject(counts = ref_in$mat[features, ], meta.data = ref_in$meta, assay = "RNA")
qry <- CreateSeuratObject(counts = qry_in$mat[features, ], meta.data = qry_in$meta, assay = "RNA")

ref <- SetAssayData(ref, assay = "RNA", layer = "data", new.data = ref_in$mat[features, ])
qry <- SetAssayData(qry, assay = "RNA", layer = "data", new.data = qry_in$mat[features, ])
VariableFeatures(ref) <- features
VariableFeatures(qry) <- features

ref <- ScaleData(ref, features = features, verbose = FALSE)
qry <- ScaleData(qry, features = features, verbose = FALSE)
ref <- RunPCA(ref, features = features, npcs = 20, verbose = FALSE)
qry <- RunPCA(qry, features = features, npcs = 20, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = ref,
  query = qry,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  dims = 1:20,
  features = features,
  verbose = TRUE
)

pred <- TransferData(
  anchorset = anchors,
  refdata = ref$science14_cell_type,
  dims = 1:20,
  verbose = TRUE
)

pred_dt <- as.data.table(pred, keep.rownames = "cell_id")
qry_meta <- as.data.table(qry@meta.data)
if (!"cell_id" %in% names(qry_meta)) {
  qry_meta[, cell_id := rownames(qry@meta.data)]
}
cell_out <- merge(qry_meta, pred_dt, by = "cell_id", all.x = TRUE, sort = FALSE)

consensus <- cell_out[, {
  tab <- sort(table(predicted.id), decreasing = TRUE)
  top_label <- names(tab)[1]
  second_label <- if (length(tab) > 1) names(tab)[2] else NA_character_
  top_n <- as.integer(tab[1])
  second_n <- if (length(tab) > 1) as.integer(tab[2]) else 0L
  top_fraction <- top_n / .N
  second_fraction <- second_n / .N
  .(
    sampled_healthy_cells = .N,
    healthy_samples = uniqueN(sample_id),
    seurat_marker_top_OneK1K_label = top_label,
    seurat_marker_top_label_cells = top_n,
    seurat_marker_top_label_fraction = top_fraction,
    seurat_marker_second_OneK1K_label = second_label,
    seurat_marker_second_label_cells = second_n,
    seurat_marker_second_label_fraction = second_fraction,
    seurat_marker_top_minus_second_fraction = top_fraction - second_fraction,
    seurat_marker_median_prediction_score_max = median(prediction.score.max, na.rm = TRUE)
  )
}, by = full_clustering]

setorder(consensus, -sampled_healthy_cells)
fwrite(cell_out, file.path(out_dir, "hca_healthy_seurat_marker_label_transfer_cell_level.tsv.gz"), sep = "\t")
fwrite(consensus, file.path(out_dir, "hca_cluster_seurat_marker_label_transfer_consensus.tsv"), sep = "\t")
fwrite(data.table(metric = c("reference_cells", "query_cells", "features", "anchors"),
                  value = c(ncol(ref), ncol(qry), length(features), length(anchors@anchors))),
       file.path(out_dir, "seurat_marker_label_transfer_qc.tsv"), sep = "\t")

print(consensus)
