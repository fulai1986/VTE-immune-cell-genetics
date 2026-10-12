#!/usr/bin/env python3
"""Centroid-based OneK1K-compatible mapping for Stephenson/HCA healthy PBMC cells.

This intentionally avoids requiring scanpy/anndata. It reads the HCA h5ad file
with h5py, builds downsampled OneK1K science14 reference centroids from raw
counts, assigns HCA Healthy cells to the nearest OneK1K centroid, and summarizes
the result at the HCA full_clustering level.
"""

from __future__ import annotations

import argparse
import gzip
import re
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
DEFAULT_RAW_DIR = PROJECT / "scratch" / "GSE196735_RAW"
DEFAULT_METADATA = PROJECT / "metadata" / "cell_metadata_science14_author_style.tsv.gz"
DEFAULT_HCA_H5AD = PROJECT / "covid_validation" / "stephenson_hca" / "covid_portal.h5ad"
DEFAULT_OUTDIR = PROJECT / "covid_validation" / "hca_onek1k_mapping" / "results"


SCORE_SETS = {
    "CD4 NC": ["CD3D", "CD3E", "TRAC", "CD4", "CCR7", "SELL", "LRRN3", "IL7R"],
    "CD4 ET": ["CD3D", "CD3E", "TRAC", "CD4", "KLRB1", "GZMK", "TNFSF13B", "IL7R", "CCL5"],
    "CD4 SOX4": ["CD3D", "CD3E", "TRAC", "CD4", "SOX4", "ID2", "SELL"],
    "CD8 NC": ["CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "LTB", "CCR7", "SELL", "PASK"],
    "CD8 ET": ["CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "GNLY", "NKG7", "GZMB", "GZMH", "PRF1", "KLRB1"],
    "CD8 S100B": ["CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "S100B", "KLRB1", "LTB", "IL7R", "GZMK"],
    "NK": ["NCAM1", "KLRD1", "KLRF1", "NKG7", "GNLY", "GZMA", "GZMB", "PRF1"],
    "NK R": ["NCAM1", "KLRD1", "GZMK", "XCL1", "XCL2", "IL7R"],
    "B IN": ["MS4A1", "CD79A", "CD79B", "TCL1A", "FCER2", "IL4R", "BANK1"],
    "B Mem": ["MS4A1", "CD79A", "CD79B", "CD27", "TNFRSF13B", "BANK1"],
    "Plasma": ["TNFRSF17", "JCHAIN", "MZB1", "XBP1", "IGHG1", "IGHA1"],
    "Mono C": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "VCAN", "LST1"],
    "Mono NC": ["FCGR3A", "MS4A7", "LST1", "CTSS", "LYZ", "IFITM3"],
    "DC": ["FCER1A", "CST3", "CLEC10A", "LILRA4", "IRF7", "SERPINF1", "ITGAX"],
}

EXTRA_FEATURES = [
    "SLC44A2",
    "FADS1",
    "MAFF",
    "TMEM258",
    "TSPAN15",
    "PPBP",
    "PF4",
    "GP9",
    "ITGA2B",
    "ITGB3",
    "GP1BA",
    "GP1BB",
    "SELP",
    "TREML1",
    "CD2",
    "CD7",
    "CD247",
    "CD14",
    "FCGR3A",
    "HLA-DRA",
    "HLA-DPB1",
]

PLATELET_CLUSTERS = {"Platelets"}
HCA_SPECIFIC_CLUSTERS = {"Platelets", "RBC", "HSC"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--hca-h5ad", type=Path, default=DEFAULT_HCA_H5AD)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--max-cells-per-type", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260820)
    parser.add_argument("--reuse-reference", action="store_true")
    parser.add_argument("--min-cluster-cells", type=int, default=50)
    parser.add_argument("--high-top-fraction", type=float, default=0.80)
    parser.add_argument("--intermediate-top-fraction", type=float, default=0.60)
    parser.add_argument("--high-score", type=float, default=0.25)
    parser.add_argument("--intermediate-score", type=float, default=0.10)
    parser.add_argument("--high-margin", type=float, default=0.20)
    parser.add_argument("--intermediate-margin", type=float, default=0.10)
    return parser.parse_args()


def feature_panel() -> list[str]:
    genes = []
    seen = set()
    for label in SCORE_SETS:
        for gene in SCORE_SETS[label]:
            if gene not in seen:
                genes.append(gene)
                seen.add(gene)
    for gene in EXTRA_FEATURES:
        if gene not in seen:
            genes.append(gene)
            seen.add(gene)
    return genes


def decode_h5ad_column(obs_group: h5py.Group, column: str) -> np.ndarray:
    values = obs_group[column][:]
    categories_path = f"__categories/{column}"
    if categories_path in obs_group:
        categories = [
            x.decode("utf-8") if isinstance(x, bytes) else str(x)
            for x in obs_group[categories_path][:]
        ]
        return np.array([categories[int(x)] if int(x) >= 0 else None for x in values], dtype=object)
    return np.array([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values], dtype=object)


def pool_number_from_raw_file(path: Path) -> int | None:
    match = re.search(r"_Pool(\d+)_RawCounts\.csv\.gz$", path.name)
    return int(match.group(1)) if match else None


def available_barcodes_by_pool(raw_dir: Path) -> dict[int, set[str]]:
    available: dict[int, set[str]] = {}
    for raw_file in sorted(raw_dir.glob("*_Pool*_RawCounts.csv.gz")):
        pool = pool_number_from_raw_file(raw_file)
        if pool is None:
            continue
        with gzip.open(raw_file, "rt") as handle:
            header = handle.readline().rstrip("\n").split(",")
        available[pool] = {x.strip('"') for x in header[1:]}
    return available


def load_and_sample_onek1k_metadata(
    path: Path, raw_dir: Path, max_cells: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(path, sep="\t", usecols=["barcode", "pool_number", "science14_cell_type"])
    meta = meta[meta["science14_cell_type"].notna()].copy()
    meta = meta[~meta["science14_cell_type"].eq("unclassified")].copy()
    meta["pool_number"] = meta["pool_number"].astype(int)
    available = available_barcodes_by_pool(raw_dir)
    meta["present_in_raw_counts"] = [
        str(barcode) in available.get(int(pool), set())
        for barcode, pool in zip(meta["barcode"], meta["pool_number"])
    ]
    availability_summary = (
        meta.groupby("science14_cell_type")
        .agg(
            metadata_cells=("barcode", "size"),
            raw_count_available_cells=("present_in_raw_counts", "sum"),
        )
        .reset_index()
    )
    meta = meta[meta["present_in_raw_counts"]].copy()
    sampled = (
        meta.groupby("science14_cell_type", group_keys=False)
        .sample(n=max_cells, replace=False, random_state=seed)
        if (meta.groupby("science14_cell_type").size() >= max_cells).all()
        else meta.groupby("science14_cell_type", group_keys=False).apply(
            lambda x: x.sample(n=min(len(x), max_cells), replace=False, random_state=seed)
        )
    )
    sampled = sampled.reset_index(drop=True)
    sampled["reference_cell_id"] = (
        sampled["pool_number"].astype(str) + ":" + sampled["barcode"].astype(str)
    )
    return sampled, availability_summary


def read_onek1k_reference_matrix(
    raw_dir: Path, selected: pd.DataFrame, features: list[str]
) -> tuple[pd.DataFrame, np.ndarray, set[str]]:
    selected = selected.reset_index(drop=True).copy()
    selected["selected_index"] = np.arange(len(selected), dtype=np.int64)
    feature_to_idx = {gene: i for i, gene in enumerate(features)}
    counts = np.zeros((len(selected), len(features)), dtype=np.float32)
    totals = np.zeros(len(selected), dtype=np.float64)
    seen_features: set[str] = set()

    by_pool: dict[int, dict[str, int]] = {}
    for row in selected.itertuples(index=False):
        by_pool.setdefault(int(row.pool_number), {})[str(row.barcode)] = int(row.selected_index)

    raw_files = sorted(raw_dir.glob("*_Pool*_RawCounts.csv.gz"))
    for raw_file in raw_files:
        pool = pool_number_from_raw_file(raw_file)
        if pool not in by_pool:
            continue
        barcode_to_idx = by_pool[pool]
        with gzip.open(raw_file, "rt") as handle:
            header = handle.readline().rstrip("\n").split(",")
        header_barcodes = [x.strip('"') for x in header[1:]]
        use_positions = [0]
        selected_columns = []
        selected_indices = []
        for pos, barcode in enumerate(header_barcodes, start=1):
            if barcode in barcode_to_idx:
                use_positions.append(pos)
                selected_columns.append(barcode)
                selected_indices.append(barcode_to_idx[barcode])
        if not selected_indices:
            continue

        for chunk in pd.read_csv(
            raw_file,
            compression="gzip",
            header=0,
            index_col=0,
            usecols=use_positions,
            chunksize=1000,
        ):
            values = chunk.to_numpy(dtype=np.float32, copy=False)
            totals[np.array(selected_indices)] += values.sum(axis=0)
            feature_rows = chunk.index.intersection(features)
            if len(feature_rows):
                row_pos = chunk.index.get_indexer(feature_rows)
                for gene, local_row in zip(feature_rows, row_pos):
                    counts[np.array(selected_indices), feature_to_idx[gene]] = values[local_row, :]
                    seen_features.add(str(gene))

    valid = totals > 0
    if not valid.all():
        missing = selected.loc[~valid, ["reference_cell_id", "science14_cell_type"]]
        raise RuntimeError(f"OneK1K selected cells with zero total counts: {missing.head().to_dict('records')}")

    norm = np.log1p((counts / totals[:, None]) * 10000.0)
    centroids = []
    for label, group in selected.groupby("science14_cell_type", sort=True):
        idx = group.index.to_numpy()
        centroids.append(pd.Series(norm[idx, :].mean(axis=0), index=features, name=label))
    ref = pd.DataFrame(centroids)
    ref.index.name = "OneK1K_label"
    return ref, norm, seen_features


def read_hca_healthy_feature_matrix(
    h5ad: Path, features: list[str]
) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    with h5py.File(h5ad, "r") as handle:
        obs = handle["obs"]
        var_names = [
            x.decode("utf-8") if isinstance(x, bytes) else str(x)
            for x in handle["var"]["names"][:]
        ]
        gene_to_col = {gene: i for i, gene in enumerate(var_names)}
        common = [gene for gene in features if gene in gene_to_col]
        common_cols = np.array([gene_to_col[gene] for gene in common], dtype=np.int64)
        col_to_feature = np.full(len(var_names), -1, dtype=np.int32)
        col_to_feature[common_cols] = np.arange(len(common), dtype=np.int32)

        status = decode_h5ad_column(obs, "Status")
        healthy_rows = np.where(status == "Healthy")[0]
        selected_obs = pd.DataFrame(
            {
                "hca_cell_id": decode_h5ad_column(obs, "_index")[healthy_rows],
                "sample_id": decode_h5ad_column(obs, "sample_id")[healthy_rows],
                "Status": status[healthy_rows],
                "Collection_Day": decode_h5ad_column(obs, "Collection_Day")[healthy_rows],
                "full_clustering": decode_h5ad_column(obs, "full_clustering")[healthy_rows],
                "Worst_Clinical_Status": decode_h5ad_column(obs, "Worst_Clinical_Status")[healthy_rows],
                "Sex": decode_h5ad_column(obs, "Sex")[healthy_rows],
                "Age": decode_h5ad_column(obs, "Age")[healthy_rows],
                "Site": decode_h5ad_column(obs, "Site")[healthy_rows],
            }
        )

        mat = np.zeros((len(healthy_rows), len(common)), dtype=np.float32)
        x = handle["X"]
        data = x["data"]
        indices = x["indices"]
        indptr = x["indptr"]
        for out_i, row_i in enumerate(healthy_rows):
            start = int(indptr[row_i])
            end = int(indptr[row_i + 1])
            cols = indices[start:end]
            feature_pos = col_to_feature[cols]
            keep = feature_pos >= 0
            if keep.any():
                mat[out_i, feature_pos[keep]] = data[start:end][keep]
    return selected_obs, mat, common


def standardize_by_reference(
    ref: pd.DataFrame, hca: np.ndarray, common_features: list[str]
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    ref_common = ref.loc[:, common_features].to_numpy(dtype=np.float64)
    mu = ref_common.mean(axis=0)
    sd = ref_common.std(axis=0, ddof=1)
    keep = np.isfinite(sd) & (sd > 1e-8)
    kept_features = [gene for gene, ok in zip(common_features, keep) if ok]
    ref_z = (ref_common[:, keep] - mu[keep]) / sd[keep]
    hca_z = (hca[:, keep].astype(np.float64) - mu[keep]) / sd[keep]
    # Cross-study expression scales are not identical; a few near-invariant
    # reference genes can otherwise dominate cosine scores.
    ref_z = np.clip(np.nan_to_num(ref_z, nan=0.0, posinf=10.0, neginf=-10.0), -10.0, 10.0)
    hca_z = np.clip(np.nan_to_num(hca_z, nan=0.0, posinf=10.0, neginf=-10.0), -10.0, 10.0)
    return ref_z, hca_z, kept_features


def cosine_scores(query: np.ndarray, ref: np.ndarray) -> np.ndarray:
    query = np.clip(np.nan_to_num(np.asarray(query, dtype=np.float64), nan=0.0), -10.0, 10.0)
    ref = np.clip(np.nan_to_num(np.asarray(ref, dtype=np.float64), nan=0.0), -10.0, 10.0)
    query_denom = np.linalg.norm(query, axis=1, keepdims=True)
    ref_denom = np.linalg.norm(ref, axis=1, keepdims=True)
    query_norm = np.divide(query, query_denom, out=np.zeros_like(query), where=query_denom > 1e-12)
    ref_norm = np.divide(ref, ref_denom, out=np.zeros_like(ref), where=ref_denom > 1e-12)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        scores = query_norm @ ref_norm.T
    return np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)


def call_confidence(
    cluster: str,
    healthy_cells: int,
    top_fraction: float,
    median_score: float,
    fraction_margin: float,
    args: argparse.Namespace,
) -> str:
    if cluster in PLATELET_CLUSTERS:
        return "Not applicable"
    if cluster in HCA_SPECIFIC_CLUSTERS:
        return "HCA-specific"
    if healthy_cells < args.min_cluster_cells:
        return "Low-cell"
    if (
        top_fraction >= args.high_top_fraction
        and median_score >= args.high_score
        and fraction_margin >= args.high_margin
    ):
        return "High"
    if (
        top_fraction >= args.intermediate_top_fraction
        and median_score >= args.intermediate_score
        and fraction_margin >= args.intermediate_margin
    ):
        return "Intermediate"
    return "Unresolved"


def marker_set_scores(cluster_means: pd.DataFrame) -> pd.DataFrame:
    marker_genes = [gene for gene in feature_panel() if gene in cluster_means.columns]
    values = cluster_means.loc[:, marker_genes].to_numpy(dtype=np.float64)
    mu = np.nanmean(values, axis=0)
    sd = np.nanstd(values, axis=0, ddof=1)
    sd[sd == 0] = np.nan
    z = (values - mu) / sd
    zdf = pd.DataFrame(z, index=cluster_means.index, columns=marker_genes)
    rows = []
    for cluster in zdf.index:
        scores = {}
        for label, genes in SCORE_SETS.items():
            present = [gene for gene in genes if gene in zdf.columns]
            scores[label] = float(np.nanmean(zdf.loc[cluster, present])) if present else np.nan
        ordered = sorted(scores.items(), key=lambda item: -np.inf if np.isnan(item[1]) else item[1], reverse=True)
        rows.append(
            {
                "full_clustering": cluster,
                "marker_top_label": ordered[0][0],
                "marker_top_score": ordered[0][1],
                "marker_second_label": ordered[1][0],
                "marker_second_score": ordered[1][1],
                "marker_score_margin": ordered[0][1] - ordered[1][1],
            }
        )
    return pd.DataFrame(rows)


def make_consensus(
    healthy_obs: pd.DataFrame,
    predictions: pd.DataFrame,
    hca_matrix: np.ndarray,
    common_features: list[str],
    ref_labels: list[str],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pred = pd.concat([healthy_obs.reset_index(drop=True), predictions], axis=1)
    cluster_means = pd.DataFrame(hca_matrix, columns=common_features).groupby(
        healthy_obs["full_clustering"].to_numpy()
    ).mean()
    marker_scores = marker_set_scores(cluster_means)

    rows = []
    similarity_rows = []
    for cluster, group in pred.groupby("full_clustering", sort=False):
        counts = group["predicted_OneK1K_label"].value_counts()
        fractions = counts / counts.sum()
        top_label = counts.index[0]
        top_count = int(counts.iloc[0])
        top_fraction = float(fractions.iloc[0])
        second_label = counts.index[1] if len(counts) > 1 else ""
        second_count = int(counts.iloc[1]) if len(counts) > 1 else 0
        second_fraction = float(fractions.iloc[1]) if len(counts) > 1 else 0.0
        median_score = float(group["prediction_score"].median())
        median_margin = float(group["prediction_score_margin"].median())
        compatible = "" if cluster in HCA_SPECIFIC_CLUSTERS else str(top_label)
        confidence = call_confidence(
            cluster,
            len(group),
            top_fraction,
            median_score,
            top_fraction - second_fraction,
            args,
        )
        if cluster in PLATELET_CLUSTERS:
            note = "HCA platelet compartment; not a OneK1K eQTL exposure cell"
        elif cluster in HCA_SPECIFIC_CLUSTERS:
            note = "HCA-specific non-PBMC/non-exposure compartment"
        elif confidence == "Unresolved":
            note = "Mixed OneK1K predictions; keep HCA full_clustering for functional analyses"
        elif "prolif" in cluster.lower():
            note = "Proliferative state; use cautiously even if nearest identity is assigned"
        else:
            note = "OneK1K-compatible centroid-based assignment"
        rows.append(
            {
                "full_clustering": cluster,
                "healthy_cells": len(group),
                "healthy_samples": group["sample_id"].nunique(),
                "OneK1K_compatible": compatible,
                "mapping_confidence": confidence,
                "top_OneK1K_label": top_label,
                "top_label_cells": top_count,
                "top_label_fraction": top_fraction,
                "second_OneK1K_label": second_label,
                "second_label_cells": second_count,
                "second_label_fraction": second_fraction,
                "top_minus_second_fraction": top_fraction - second_fraction,
                "median_prediction_score": median_score,
                "median_prediction_score_margin": median_margin,
                "mapping_note": note,
            }
        )

        for label in ref_labels:
            similarity_rows.append(
                {
                    "full_clustering": cluster,
                    "OneK1K_label": label,
                    "mean_cell_prediction_score": float(group[f"score::{label}"].mean()),
                    "median_cell_prediction_score": float(group[f"score::{label}"].median()),
                }
            )

    consensus = pd.DataFrame(rows).merge(marker_scores, on="full_clustering", how="left")
    consensus = consensus.sort_values(["healthy_cells"], ascending=False)
    similarity = pd.DataFrame(similarity_rows)
    return consensus, similarity


def build_all_obs_mapping(h5ad: Path, consensus: pd.DataFrame) -> pd.DataFrame:
    cluster_map = consensus.set_index("full_clustering")
    with h5py.File(h5ad, "r") as handle:
        obs = handle["obs"]
        full = decode_h5ad_column(obs, "full_clustering")
        out = pd.DataFrame(
            {
                "hca_cell_id": decode_h5ad_column(obs, "_index"),
                "sample_id": decode_h5ad_column(obs, "sample_id"),
                "Status": decode_h5ad_column(obs, "Status"),
                "Collection_Day": decode_h5ad_column(obs, "Collection_Day"),
                "full_clustering": full,
            }
        )
    out["OneK1K_compatible"] = out["full_clustering"].map(cluster_map["OneK1K_compatible"]).fillna("")
    out["mapping_confidence"] = out["full_clustering"].map(cluster_map["mapping_confidence"]).fillna("Unmapped")
    return out


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    features = feature_panel()
    pd.DataFrame({"feature_gene": features}).to_csv(
        args.outdir / "hca_onek1k_mapping_feature_panel.tsv", sep="\t", index=False
    )

    selected, availability_summary = load_and_sample_onek1k_metadata(
        args.metadata, args.raw_dir, args.max_cells_per_type, args.seed
    )
    availability_summary.to_csv(
        args.outdir / "onek1k_metadata_raw_count_availability_summary.tsv", sep="\t", index=False
    )
    selected.groupby("science14_cell_type").size().rename("n_reference_cells").reset_index().to_csv(
        args.outdir / "onek1k_reference_sampling_summary.tsv", sep="\t", index=False
    )

    ref_path = args.outdir / "onek1k_reference_centroids_logcp10k.tsv"
    if args.reuse_reference and ref_path.exists():
        ref = pd.read_csv(ref_path, sep="\t", index_col=0)
        seen_reference_features = set(ref.columns)
    else:
        ref, _, seen_reference_features = read_onek1k_reference_matrix(args.raw_dir, selected, features)
        ref.to_csv(ref_path, sep="\t")

    hca_obs, hca_matrix, hca_features = read_hca_healthy_feature_matrix(args.hca_h5ad, features)
    common_features = [gene for gene in hca_features if gene in seen_reference_features]
    hca_common = hca_matrix[:, [hca_features.index(gene) for gene in common_features]]
    ref_z, hca_z, kept_features = standardize_by_reference(ref, hca_common, common_features)
    ref_labels = ref.index.tolist()
    scores = cosine_scores(hca_z, ref_z)
    order = np.argsort(-scores, axis=1)
    top_idx = order[:, 0]
    second_idx = order[:, 1]

    predictions = pd.DataFrame(
        {
            "predicted_OneK1K_label": [ref_labels[i] for i in top_idx],
            "prediction_score": scores[np.arange(scores.shape[0]), top_idx],
            "second_OneK1K_label": [ref_labels[i] for i in second_idx],
            "second_prediction_score": scores[np.arange(scores.shape[0]), second_idx],
        }
    )
    predictions["prediction_score_margin"] = (
        predictions["prediction_score"] - predictions["second_prediction_score"]
    )
    for i, label in enumerate(ref_labels):
        predictions[f"score::{label}"] = scores[:, i]

    healthy_cell_predictions = pd.concat([hca_obs.reset_index(drop=True), predictions], axis=1)
    healthy_cell_predictions.to_csv(
        args.outdir / "hca_healthy_onek1k_label_transfer_cell_level.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    cluster_counts = hca_obs.groupby("full_clustering").agg(
        healthy_cells=("hca_cell_id", "size"), healthy_samples=("sample_id", "nunique")
    )
    cluster_counts.sort_values("healthy_cells", ascending=False).reset_index().to_csv(
        args.outdir / "hca_healthy_cluster_counts.tsv", sep="\t", index=False
    )

    consensus, similarity = make_consensus(
        hca_obs,
        predictions,
        hca_common[:, [common_features.index(gene) for gene in kept_features]],
        kept_features,
        ref_labels,
        args,
    )
    consensus.to_csv(args.outdir / "hca_cluster_onek1k_mapping_consensus.tsv", sep="\t", index=False)
    similarity.to_csv(args.outdir / "hca_cluster_onek1k_similarity_matrix.tsv", sep="\t", index=False)

    all_obs = build_all_obs_mapping(args.hca_h5ad, consensus)
    all_obs.to_csv(
        args.outdir / "hca_obs_with_onek1k_compatible_annotation.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    qc = pd.DataFrame(
        [
            {"metric": "reference_cell_types", "value": len(ref_labels)},
            {"metric": "selected_reference_cells", "value": len(selected)},
            {"metric": "requested_features", "value": len(features)},
            {"metric": "reference_features_seen", "value": len(seen_reference_features)},
            {"metric": "hca_features_found", "value": len(hca_features)},
            {"metric": "mapping_features_used", "value": len(kept_features)},
            {"metric": "hca_healthy_cells", "value": len(hca_obs)},
            {"metric": "hca_healthy_clusters", "value": hca_obs["full_clustering"].nunique()},
        ]
    )
    qc.to_csv(args.outdir / "hca_onek1k_centroid_mapping_qc.tsv", sep="\t", index=False)

    focus = consensus[
        consensus["full_clustering"].isin(
            ["CD8.TE", "CD8.EM", "NK_16hi", "NK_56hi", "NK_prolif", "CD14_mono", "CD16_mono", "Platelets"]
        )
    ]
    print("QC")
    print(qc.to_string(index=False))
    print("\nFocused HCA cluster consensus")
    print(
        focus[
            [
                "full_clustering",
                "healthy_cells",
                "OneK1K_compatible",
                "mapping_confidence",
                "top_OneK1K_label",
                "top_label_fraction",
                "second_OneK1K_label",
                "second_label_fraction",
                "median_prediction_score",
                "marker_top_label",
                "mapping_note",
            ]
        ].to_string(index=False)
    )
    print(f"\nWrote outputs to {args.outdir}")


if __name__ == "__main__":
    main()
