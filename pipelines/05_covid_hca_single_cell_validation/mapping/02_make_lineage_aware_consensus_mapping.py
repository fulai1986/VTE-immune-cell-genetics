#!/usr/bin/env python3
"""Add lineage-aware curated OneK1K-compatible labels for all HCA clusters.

The centroid-based centroid mapping is useful evidence, but broad CD4/CD8 states can
be pulled toward cytotoxic or naive programs when using a small feature panel.
This script keeps the original algorithmic result and adds a curated mapping
that respects the HCA full_clustering lineage names.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


BASE = Path(os.environ["PROJECT_DIR"]).expanduser().resolve() / "covid_validation/hca_onek1k_mapping"
RESULTS = BASE / "results"


CURATED = {
    # CD4 T-cell states. HCA has detailed disease/reference states; OneK1K has
    # broader CD4 NC, CD4 ET, and CD4 SOX4 compartments.
    "CD4.Naive": ("CD4 NC", "High", "CD4 naive/central-memory-like state"),
    "CD4.CM": ("CD4 NC", "High", "CD4 central-memory-like state"),
    "CD4.IL22": ("CD4 ET", "Intermediate", "CD4 effector/helper state"),
    "CD4.Tfh": ("CD4 ET", "Intermediate", "CD4 helper effector-like state"),
    "Treg": ("CD4 ET", "Intermediate", "CD4 regulatory state mapped to nearest CD4 effector-like OneK1K compartment"),
    "CD4.EM": ("CD4 ET", "Intermediate", "CD4 effector-memory state"),
    "CD4.Th1": ("CD4 ET", "Intermediate", "CD4 helper effector state"),
    "CD4.Th2": ("CD4 ET", "Low-cell", "Low-cell CD4 helper state"),
    "CD4.Th17": ("CD4 ET", "Low-cell", "Very low-cell CD4 helper state"),
    "CD4.Prolif": ("CD4 ET", "Low-cell", "Proliferative CD4 state; identity mapped but state retained"),
    # CD8 and unconventional T-cell states.
    "CD8.Naive": ("CD8 NC", "High", "CD8 naive/central-memory-like state"),
    "CD8.TE": ("CD8 ET", "High", "CD8 terminal effector state; centroid-based mapping also high confidence"),
    "CD8.EM": ("CD8 ET", "High", "CD8 effector-memory state; centroid-based mapping also high confidence"),
    "CD8.Prolif": ("CD8 ET", "Intermediate", "Proliferative CD8 state; identity mapped but state retained"),
    "MAIT": ("CD8 S100B", "Intermediate", "MAIT-like T state closest to OneK1K CD8 S100B-like compartment"),
    "gdT": ("CD8 ET", "Intermediate", "Cytotoxic unconventional T state mapped to CD8 effector-like compartment"),
    "NKT": ("CD8 ET", "Intermediate", "Mixed T/NK state; retain HCA full_clustering for functional analysis"),
    # NK/ILC states.
    "NK_16hi": ("NK", "High", "CD16-high cytotoxic NK state; centroid-based mapping high confidence"),
    "NK_56hi": ("NK R", "Intermediate", "CD56-bright-like NK state; directionally maps to OneK1K NK R but mixed with NK"),
    "NK_prolif": ("NK", "Intermediate", "Proliferative NK state; identity mapped but state retained"),
    "ILC1_3": ("NK R", "Intermediate", "ILC/NK-related state mapped to nearest NK_R-like compartment"),
    "ILC2": ("NK R", "Low-cell", "Low-cell ILC state mapped to nearest NK_R-like compartment"),
    # B-cell and plasma states.
    "B_naive": ("B IN", "High", "Naive/intermediate B-cell state"),
    "B_immature": ("B IN", "Intermediate", "Immature B-cell state mapped to B IN"),
    "B_switched_memory": ("B Mem", "High", "Memory B-cell state"),
    "B_non-switched_memory": ("B Mem", "High", "Memory B-cell state"),
    "B_exhausted": ("B Mem", "High", "Memory/exhausted B-cell state"),
    "Plasma_cell_IgA": ("Plasma", "High", "Plasma cell state; marker sanity supports Plasma"),
    "Plasma_cell_IgG": ("Plasma", "High", "Plasma cell state"),
    "Plasma_cell_IgM": ("Plasma", "Low-cell", "Low-cell plasma cell state"),
    "Plasmablast": ("Plasma", "High", "Plasmablast/plasma state"),
    # Monocytes.
    "CD83_CD14_mono": ("Mono C", "High", "Classical/inflammatory CD14 monocyte state"),
    "CD14_mono": ("Mono C", "High", "Classical monocyte state"),
    "CD16_mono": ("Mono NC", "High", "Non-classical monocyte state"),
    "C1_CD16_mono": ("Mono NC", "High", "C1-positive non-classical monocyte state"),
    "Mono_prolif": ("Mono C", "Low-cell", "Low-cell proliferative monocyte state"),
    # DC states.
    "DC2": ("DC", "High", "Dendritic cell state"),
    "DC3": ("DC", "Intermediate", "Dendritic/monocyte-adjacent state mapped to DC"),
    "pDC": ("DC", "High", "Plasmacytoid dendritic cell state mapped to OneK1K DC"),
    "DC1": ("DC", "High", "Dendritic cell state"),
    "ASDC": ("DC", "Low-cell", "Low-cell dendritic cell state"),
    "DC_prolif": ("DC", "Low-cell", "Low-cell proliferative dendritic state"),
    # Not OneK1K eQTL exposure compartments.
    "Platelets": ("", "Not applicable", "Platelet compartment; downstream thromboinflammatory effector, not OneK1K exposure"),
    "RBC": ("", "HCA-specific", "RBC-like compartment; no OneK1K 14-cell exposure mapping"),
    "HSC": ("", "HCA-specific", "HSC-like compartment; no OneK1K 14-cell exposure mapping"),
}


def main() -> None:
    consensus_path = RESULTS / "hca_cluster_onek1k_mapping_consensus.tsv"
    obs_path = RESULTS / "hca_obs_with_onek1k_compatible_annotation.tsv.gz"
    consensus = pd.read_csv(consensus_path, sep="\t")
    missing = sorted(set(consensus["full_clustering"]) - set(CURATED))
    if missing:
        raise RuntimeError(f"Missing curated mapping for HCA clusters: {missing}")

    curated_rows = []
    for cluster in consensus["full_clustering"]:
        label, confidence, note = CURATED[cluster]
        curated_rows.append(
            {
                "full_clustering": cluster,
                "curated_OneK1K_compatible": label,
                "curated_mapping_confidence": confidence,
                "curated_mapping_note": note,
            }
        )
    curated = pd.DataFrame(curated_rows)
    out = consensus.merge(curated, on="full_clustering", how="left")
    out.to_csv(RESULTS / "hca_cluster_onek1k_mapping_consensus.curated.tsv", sep="\t", index=False)

    obs = pd.read_csv(obs_path, sep="\t")
    obs = obs.drop(columns=["OneK1K_compatible", "mapping_confidence"], errors="ignore")
    obs = obs.merge(curated, on="full_clustering", how="left")
    obs.to_csv(
        RESULTS / "hca_obs_with_onek1k_compatible_annotation.curated.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    print("Curated mapping counts")
    print(out["curated_mapping_confidence"].value_counts().to_string())
    print("\nCD4/CD8/NK focused curated mapping")
    focus = out[out["full_clustering"].str.contains("CD4|CD8|NK|MAIT|gdT|Treg|ILC", regex=True)]
    print(
        focus[
            [
                "full_clustering",
                "healthy_cells",
                "curated_OneK1K_compatible",
                "curated_mapping_confidence",
                "top_OneK1K_label",
                "top_label_fraction",
                "marker_top_label",
                "curated_mapping_note",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
