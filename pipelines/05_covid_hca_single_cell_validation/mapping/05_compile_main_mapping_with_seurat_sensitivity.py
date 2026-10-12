#!/usr/bin/env python3
"""Combine the lineage-aware primary mapping with Seurat sensitivity results."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
BASE = PROJECT / "covid_validation/hca_onek1k_mapping"
CONSENSUS = BASE / "results/hca_cluster_onek1k_mapping_consensus.curated.tsv"
SEURAT = (
    BASE
    / "seurat_label_transfer/results_marker_panel/"
    "hca_cluster_seurat_marker_label_transfer_consensus.tsv"
)
OUT = BASE / "results/hca_to_onek1k_mapping_main_with_seurat_sensitivity.tsv"


def main() -> None:
    primary = pd.read_csv(CONSENSUS, sep="\t")
    seurat = pd.read_csv(SEURAT, sep="\t")

    keep = primary[
        [
            "full_clustering",
            "healthy_cells",
            "healthy_samples",
            "curated_OneK1K_compatible",
            "curated_mapping_confidence",
            "top_OneK1K_label",
            "top_label_fraction",
            "marker_top_label",
            "curated_mapping_note",
        ]
    ].copy()
    keep = keep.rename(
        columns={
            "full_clustering": "HCA cluster",
            "healthy_cells": "Healthy cells",
            "healthy_samples": "Healthy samples",
            "curated_OneK1K_compatible": "Main OneK1K-compatible label",
            "curated_mapping_confidence": "Main confidence",
            "top_OneK1K_label": "Centroid top label",
            "top_label_fraction": "Centroid top fraction",
            "marker_top_label": "Marker sanity label",
            "curated_mapping_note": "Note",
        }
    )

    sensitivity = seurat[
        [
            "full_clustering",
            "seurat_marker_top_OneK1K_label",
            "seurat_marker_top_label_fraction",
            "seurat_marker_second_OneK1K_label",
            "seurat_marker_median_prediction_score_max",
        ]
    ].rename(
        columns={
            "full_clustering": "HCA cluster",
            "seurat_marker_top_OneK1K_label": "Seurat top label",
            "seurat_marker_top_label_fraction": "Seurat top fraction",
            "seurat_marker_second_OneK1K_label": "Seurat second label",
            "seurat_marker_median_prediction_score_max": "Seurat median score",
        }
    )

    out = keep.merge(sensitivity, on="HCA cluster", how="left", validate="one_to_one")
    if out["HCA cluster"].duplicated().any():
        raise RuntimeError("Duplicate HCA clusters after mapping merge")
    if len(out) != len(primary):
        raise RuntimeError("Mapping row count changed during Seurat merge")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False, na_rep="NA")
    print(f"clusters={len(out)}")
    print(f"seurat_results={out['Seurat top label'].notna().sum()}")
    print(OUT)


if __name__ == "__main__":
    main()
