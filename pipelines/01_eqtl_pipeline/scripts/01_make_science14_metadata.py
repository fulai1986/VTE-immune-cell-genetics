#!/usr/bin/env python3
"""Create Science14 author-style cell metadata for OneK1K PBMC eQTL."""

from __future__ import annotations

import csv
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
PATHS_FILE = PROJECT_DIR / "config" / "project_paths.tsv"
METADATA_DIR = PROJECT_DIR / "metadata"
QC_DIR = PROJECT_DIR / "results" / "qc"
REPORT_DIR = PROJECT_DIR / "reports"


MARKER_GENES = [
    "CD3D",
    "CD4",
    "CCR7",
    "SELL",
    "LRRN3",
    "KLRB1",
    "GZMK",
    "TNFSF13B",
    "IL7R",
    "SOX4",
    "ID2",
    "CD8A",
    "LTB",
    "PASK",
    "GNLY",
    "NKG7",
    "S100B",
    "SLC4A10",
    "NCAM1",
    "GZMA",
    "GZMB",
    "XCL1",
    "XCL2",
    "MS4A1",
    "TCL1A",
    "FCER2",
    "IL4R",
    "CD27",
    "TNFRSF13B",
    "TNFRSF17",
    "JCHAIN",
    "CD14",
    "LYZ",
    "FCGR3A",
    "FCER1A",
    "CST3",
    "SERPINF1",
]


SCORE_SETS = {
    "CD4 NC": ["CD3D", "CD4", "CCR7", "SELL", "LRRN3"],
    "CD4 ET": ["CD3D", "CD4", "KLRB1", "GZMK", "TNFSF13B", "IL7R"],
    "CD4 SOX4": ["CD3D", "CD4", "SOX4", "ID2", "SELL"],
    "CD8 NC": ["CD3D", "CD8A", "LTB", "CCR7", "PASK"],
    "CD8 ET": ["CD3D", "CD8A", "GNLY", "NKG7", "KLRB1"],
    "CD8 S100B": ["CD3D", "CD8A", "S100B", "KLRB1", "LTB", "IL7R", "GZMK"],
    "NK": ["NCAM1", "GZMA", "GZMB"],
    "NK R": ["NCAM1", "GZMK", "XCL1", "XCL2"],
    "B IN": ["MS4A1", "TCL1A", "FCER2", "IL4R"],
    "B Mem": ["MS4A1", "CD27", "TNFRSF13B"],
    "Plasma": ["TNFRSF17", "JCHAIN"],
    "Mono C": ["CD14", "LYZ"],
    "Mono NC": ["FCGR3A"],
    "DC": ["FCER1A", "CST3", "SERPINF1"],
}


def decode_array(values) -> np.ndarray:
    arr = np.asarray(values)
    if arr.dtype.kind in {"S", "O"}:
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])
    return arr.astype(str)


def read_paths() -> dict[str, Path]:
    with PATHS_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: Path(row["path"]) for row in reader}


def read_cat_obs(h5: h5py.File, col: str) -> np.ndarray:
    base = f"obs/{col}"
    if f"{base}/codes" in h5:
        codes = h5[f"{base}/codes"][:].astype(int)
        cats = decode_array(h5[f"{base}/categories"][:])
        return cats[codes]
    return decode_array(h5[base][:])


def read_var_gene_symbols(h5: h5py.File) -> np.ndarray:
    base = "var/feature_name"
    if f"{base}/codes" in h5:
        codes = h5[f"{base}/codes"][:].astype(int)
        cats = decode_array(h5[f"{base}/categories"][:])
        return cats[codes]
    return decode_array(h5[base][:])


def scaled_targets(targets: pd.DataFrame, labels: list[str], observed_total: int) -> dict[str, int]:
    original = targets.set_index("science14_cell_type").loc[labels, "original_cells"].astype(float)
    raw = original / original.sum() * observed_total
    n = np.floor(raw).astype(int)
    remainder = int(observed_total - n.sum())
    if remainder > 0:
        order = np.argsort(-(raw.to_numpy() - n.to_numpy()))
        for idx in order[:remainder]:
            n.iloc[idx] += 1
    return {label: int(n.loc[label]) for label in labels}


def assign_top(
    labels: np.ndarray,
    notes: np.ndarray,
    available: np.ndarray,
    scores: dict[str, np.ndarray],
    label: str,
    n: int,
    note: str,
) -> np.ndarray:
    if n <= 0:
        return available
    candidate = np.where(available & np.isfinite(scores[label]))[0]
    if candidate.size == 0:
        return available
    ordered = candidate[np.argsort(-scores[label][candidate])]
    keep = ordered[: min(n, ordered.size)]
    labels[keep] = label
    notes[keep] = note
    available[keep] = False
    return available


def score_set(marker_expr: np.ndarray, marker_names: list[str], geneset: list[str]) -> np.ndarray:
    idx = [marker_names.index(g) for g in geneset if g in marker_names]
    if not idx:
        return np.full(marker_expr.shape[0], np.nan)
    return marker_expr[:, idx].mean(axis=1)


def main() -> int:
    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    paths = read_paths()
    h5ad_file = paths["cellxgene_h5ad"]
    table_s7_xlsx = paths["science_tables_xlsx"]

    print("Reading Science Table S7 targets...")
    s7 = pd.read_excel(table_s7_xlsx, sheet_name="Table.S7", skiprows=2)
    s7.columns = [str(c).strip() for c in s7.columns]
    s7 = s7[s7["Cell Type"].notna()].copy()
    s7["Number of Cells"] = pd.to_numeric(s7["Number of Cells"], errors="coerce").fillna(0).astype(int)
    targets = (
        s7.groupby("Cell Type", as_index=False)["Number of Cells"]
        .sum()
        .rename(columns={"Cell Type": "science14_cell_type", "Number of Cells": "original_cells"})
    )

    print("Reading h5ad metadata...")
    with h5py.File(h5ad_file, "r") as h5:
        barcode = read_cat_obs(h5, "barcode")
        donor_id = read_cat_obs(h5, "donor_id")
        pool_number = read_cat_obs(h5, "pool_number")
        cell_type_cxg = read_cat_obs(h5, "cell_type")
        pred_l2 = read_cat_obs(h5, "predicted.celltype.l2")
        genes = np.char.upper(read_var_gene_symbols(h5).astype(str))

        gene_to_idx = {gene: i for i, gene in enumerate(genes)}
        marker_names = [g for g in MARKER_GENES if g in gene_to_idx]
        missing = [g for g in MARKER_GENES if g not in gene_to_idx]
        if missing:
            print("Missing marker genes:", ", ".join(missing))

        print("Extracting marker expression from h5ad CSR matrix...")
        n_cells = len(barcode)
        marker_expr = np.zeros((n_cells, len(marker_names)), dtype=np.float32)
        target_lookup = {gene_to_idx[g]: j for j, g in enumerate(marker_names)}
        indptr = h5["X/indptr"][:]
        chunk_cells = 50_000
        for start in range(0, n_cells, chunk_cells):
            end = min(n_cells, start + chunk_cells)
            start_nnz = int(indptr[start])
            end_nnz = int(indptr[end])
            if end_nnz > start_nnz:
                idx = h5["X/indices"][start_nnz:end_nnz]
                val = h5["X/data"][start_nnz:end_nnz]
                nnz_per_cell = np.diff(indptr[start : end + 1]).astype(int)
                cell_ids = np.repeat(np.arange(start, end), nnz_per_cell)
                hit_mask = np.isin(idx, list(target_lookup.keys()))
                if hit_mask.any():
                    cols = np.array([target_lookup[int(i)] for i in idx[hit_mask]], dtype=int)
                    marker_expr[cell_ids[hit_mask], cols] = val[hit_mask]
            if start == 0 or end == n_cells or start % 200_000 == 0:
                print(f"Processed cells {start + 1}-{end}")

    scores = {label: score_set(marker_expr, marker_names, geneset) for label, geneset in SCORE_SETS.items()}

    science_labels = np.full(len(barcode), "unclassified", dtype=object)
    notes = np.full(len(barcode), "not_in_author_style_gate", dtype=object)

    pred = pred_l2.astype(str)
    t_gate = np.isin(
        pred,
        [
            "CD4 TCM",
            "CD4 Naive",
            "CD4 TEM",
            "Treg",
            "CD4 CTL",
            "CD4 Proliferating",
            "CD8 TEM",
            "CD8 Naive",
            "CD8 TCM",
            "CD8 Proliferating",
            "MAIT",
            "gdT",
            "dnT",
        ],
    )
    b_gate = np.isin(pred, ["B naive", "B intermediate", "B memory", "Plasmablast"])
    nk_gate = np.isin(pred, ["NK", "NK_CD56bright", "NK Proliferating", "ILC"])
    mono_c_gate = pred == "CD14 Mono"
    mono_nc_gate = pred == "CD16 Mono"
    dc_gate = np.isin(pred, ["cDC2", "cDC1", "pDC", "ASDC"])

    print("Assigning T-cell branch...")
    t_labels = ["CD4 SOX4", "CD8 S100B", "CD4 ET", "CD8 NC", "CD8 ET", "CD4 NC"]
    t_targets = scaled_targets(targets, t_labels, int(t_gate.sum()))
    available = t_gate.copy()
    for lab in t_labels[:-1]:
        available = assign_top(
            science_labels, notes, available, scores, lab, t_targets[lab], "T_branch_TableS7_scaled_marker_rank"
        )
    science_labels[available] = "CD4 NC"
    notes[available] = "T_branch_remaining_CD4_NC"

    print("Assigning B-cell branch...")
    b_labels = ["Plasma", "B Mem", "B IN"]
    b_targets = scaled_targets(targets, b_labels, int(b_gate.sum()))
    available = b_gate.copy()
    for lab in b_labels[:-1]:
        available = assign_top(
            science_labels, notes, available, scores, lab, b_targets[lab], "B_branch_TableS7_scaled_marker_rank"
        )
    science_labels[available] = "B IN"
    notes[available] = "B_branch_remaining_B_IN"

    print("Assigning NK branch...")
    nk_targets = scaled_targets(targets, ["NK R", "NK"], int(nk_gate.sum()))
    available = nk_gate.copy()
    available = assign_top(
        science_labels, notes, available, scores, "NK R", nk_targets["NK R"], "NK_branch_TableS7_scaled_marker_rank"
    )
    science_labels[available] = "NK"
    notes[available] = "NK_branch_remaining_NK"

    science_labels[mono_c_gate] = "Mono C"
    notes[mono_c_gate] = "direct_l2_CD14_mono"
    science_labels[mono_nc_gate] = "Mono NC"
    notes[mono_nc_gate] = "direct_l2_CD16_mono"
    science_labels[dc_gate] = "DC"
    notes[dc_gate] = "direct_l2_DC"

    marker_score = np.full(len(barcode), np.nan)
    for lab, values in scores.items():
        mask = science_labels == lab
        marker_score[mask] = values[mask]

    dt = pd.DataFrame(
        {
            "barcode": barcode,
            "donor_id": donor_id,
            "pool_number": pool_number,
            "cell_type_cellxgene": cell_type_cxg,
            "predicted_celltype_l2": pred_l2,
            "science14_cell_type": science_labels,
            "classification_note": notes,
            "marker_score_for_assigned_label": marker_score,
        }
    )

    out_file = METADATA_DIR / "cell_metadata_science14_author_style.tsv.gz"
    score_file = QC_DIR / "science14_marker_scores.tsv.gz"
    counts_file = QC_DIR / "science14_celltype_counts.tsv"
    donor_counts_file = QC_DIR / "science14_donor_cell_counts.tsv.gz"
    compare_file = QC_DIR / "science_table_s7_vs_author_style_hierarchical_counts.tsv"
    report_file = REPORT_DIR / "STEP1_SCIENCE14_METADATA.md"

    print("Writing outputs...")
    dt.to_csv(out_file, sep="\t", index=False)
    pd.DataFrame({"barcode": barcode, **scores}).to_csv(score_file, sep="\t", index=False)

    counts = dt.groupby("science14_cell_type", as_index=False).size().rename(columns={"size": "author_style_cells"})
    counts = counts.sort_values("author_style_cells", ascending=False)
    counts.to_csv(counts_file, sep="\t", index=False)

    donor_counts = (
        dt[dt["science14_cell_type"] != "unclassified"]
        .groupby(["donor_id", "science14_cell_type"], as_index=False)
        .size()
        .rename(columns={"size": "N"})
        .sort_values(["science14_cell_type", "donor_id"])
    )
    donor_counts.to_csv(donor_counts_file, sep="\t", index=False)

    compare = targets.merge(counts, on="science14_cell_type", how="outer").fillna(0)
    compare["original_cells"] = compare["original_cells"].astype(int)
    compare["author_style_cells"] = compare["author_style_cells"].astype(int)
    compare["diff"] = compare["author_style_cells"] - compare["original_cells"]
    compare["pct_of_original"] = np.where(
        compare["original_cells"] > 0, 100 * compare["author_style_cells"] / compare["original_cells"], np.nan
    )
    compare["original_pct"] = 100 * compare["original_cells"] / compare["original_cells"].sum()
    compare["author_style_pct"] = 100 * compare["author_style_cells"] / compare["author_style_cells"].sum()
    compare["pct_point_diff"] = compare["author_style_pct"] - compare["original_pct"]
    order_labels = [
        "CD4 NC",
        "CD4 ET",
        "CD4 SOX4",
        "CD8 ET",
        "CD8 NC",
        "CD8 S100B",
        "NK",
        "NK R",
        "B IN",
        "B Mem",
        "Plasma",
        "Mono C",
        "Mono NC",
        "DC",
        "unclassified",
    ]
    compare["_order"] = compare["science14_cell_type"].map({v: i for i, v in enumerate(order_labels)})
    compare = compare.sort_values(["_order", "science14_cell_type"]).drop(columns=["_order"])
    compare.to_csv(compare_file, sep="\t", index=False)

    report_file.write_text(
        "\n".join(
            [
                "# Step 1 Science14 Metadata",
                "",
                f"- Cells classified: {len(dt)}",
                f"- Donors represented: {dt['donor_id'].nunique()}",
                f"- Science14 classes plus unclassified: {dt['science14_cell_type'].nunique()}",
                f"- Main metadata: `{out_file.relative_to(PROJECT_DIR)}`",
                f"- Cell-type counts: `{counts_file.relative_to(PROJECT_DIR)}`",
                f"- Donor cell counts: `{donor_counts_file.relative_to(PROJECT_DIR)}`",
                f"- Table S7 comparison: `{compare_file.relative_to(PROJECT_DIR)}`",
                "",
                "## Notes",
                "",
                "- Classification logic follows the prior author-style hierarchical Science14 script.",
                "- Outputs are contained inside the new eQTL project directory.",
            ]
        )
        + "\n"
    )

    print(compare)
    print(f"Wrote {out_file}")
    print(f"Wrote {counts_file}")
    print(f"Wrote {donor_counts_file}")
    print(f"Wrote {compare_file}")
    print(f"Wrote {report_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
