#!/usr/bin/env python3
"""Prepare sampled OneK1K and HCA Healthy matrices for Seurat label transfer.

The HCA covid_portal.h5ad stores normalized expression values in X, whereas the
OneK1K GEO matrices are raw counts. This script outputs a compact, shared-gene
matrix set for a Seurat sensitivity analysis; the R step uses log-normalized
expression layers rather than treating HCA X as raw UMI counts.
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
from scipy import sparse
from scipy.io import mmwrite


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
DEFAULT_RAW_DIR = PROJECT / "scratch" / "GSE196735_RAW"
DEFAULT_METADATA = PROJECT / "metadata" / "cell_metadata_science14_author_style.tsv.gz"
DEFAULT_HCA_H5AD = PROJECT / "covid_validation" / "stephenson_hca" / "covid_portal.h5ad"
DEFAULT_OUTDIR = PROJECT / "covid_validation" / "hca_onek1k_mapping" / "seurat_label_transfer"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--hca-h5ad", type=Path, default=DEFAULT_HCA_H5AD)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--max-ref-cells-per-type", type=int, default=1000)
    parser.add_argument("--max-query-cells-per-cluster", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260827)
    return parser.parse_args()


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


def read_onek1k_gene_list(raw_dir: Path) -> list[str]:
    first = sorted(raw_dir.glob("*_Pool*_RawCounts.csv.gz"))[0]
    genes = []
    for chunk in pd.read_csv(first, compression="gzip", header=0, index_col=0, usecols=[0], chunksize=10000):
        genes.extend([str(x) for x in chunk.index])
    return genes


def read_hca_gene_list(h5ad: Path) -> list[str]:
    with h5py.File(h5ad, "r") as handle:
        return [
            x.decode("utf-8") if isinstance(x, bytes) else str(x)
            for x in handle["var"]["names"][:]
        ]


def select_onek1k_cells(metadata: Path, raw_dir: Path, max_cells: int, seed: int) -> pd.DataFrame:
    meta = pd.read_csv(metadata, sep="\t", usecols=["barcode", "pool_number", "science14_cell_type"])
    meta = meta[meta["science14_cell_type"].notna()].copy()
    meta = meta[~meta["science14_cell_type"].eq("unclassified")].copy()
    meta["pool_number"] = meta["pool_number"].astype(int)
    available = available_barcodes_by_pool(raw_dir)
    meta = meta[
        [
            str(barcode) in available.get(int(pool), set())
            for barcode, pool in zip(meta["barcode"], meta["pool_number"])
        ]
    ].copy()
    sampled = meta.groupby("science14_cell_type", group_keys=False).apply(
        lambda x: x.sample(n=min(len(x), max_cells), replace=False, random_state=seed)
    )
    sampled = sampled.reset_index(drop=True)
    sampled["cell_id"] = (
        "OneK1K_pool" + sampled["pool_number"].astype(str) + "_" + sampled["barcode"].astype(str)
    )
    return sampled


def select_hca_healthy_cells(h5ad: Path, max_cells: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    with h5py.File(h5ad, "r") as handle:
        obs = handle["obs"]
        status = decode_h5ad_column(obs, "Status")
        full = decode_h5ad_column(obs, "full_clustering")
        sample_id = decode_h5ad_column(obs, "sample_id")
        cell_id = decode_h5ad_column(obs, "_index")
    healthy_idx = np.where(status == "Healthy")[0]
    rows = []
    for cluster in pd.Series(full[healthy_idx]).dropna().unique():
        idx = healthy_idx[full[healthy_idx] == cluster]
        if len(idx) > max_cells:
            idx = rng.choice(idx, size=max_cells, replace=False)
        for i in idx:
            rows.append(
                {
                    "hca_row_index": int(i),
                    "cell_id": f"HCA_{cell_id[i]}",
                    "full_clustering": full[i],
                    "sample_id": sample_id[i],
                    "Status": status[i],
                }
            )
    return pd.DataFrame(rows)


def write_onek1k_matrix(
    raw_dir: Path,
    selected: pd.DataFrame,
    common_genes: list[str],
    outdir: Path,
) -> None:
    gene_to_row = {gene: i for i, gene in enumerate(common_genes)}
    selected = selected.reset_index(drop=True).copy()
    selected["col_index"] = np.arange(selected.shape[0], dtype=np.int64)
    by_pool: dict[int, dict[str, int]] = {}
    for row in selected.itertuples(index=False):
        by_pool.setdefault(int(row.pool_number), {})[str(row.barcode)] = int(row.col_index)

    data_parts = []
    row_parts = []
    col_parts = []
    totals = np.zeros(selected.shape[0], dtype=np.float64)

    for raw_file in sorted(raw_dir.glob("*_Pool*_RawCounts.csv.gz")):
        pool = pool_number_from_raw_file(raw_file)
        if pool not in by_pool:
            continue
        barcode_to_col = by_pool[pool]
        with gzip.open(raw_file, "rt") as handle:
            header = handle.readline().rstrip("\n").split(",")
        header_barcodes = [x.strip('"') for x in header[1:]]
        use_positions = [0]
        selected_cols = []
        for pos, barcode in enumerate(header_barcodes, start=1):
            if barcode in barcode_to_col:
                use_positions.append(pos)
                selected_cols.append(barcode_to_col[barcode])
        if not selected_cols:
            continue
        selected_cols_array = np.array(selected_cols, dtype=np.int64)
        for chunk in pd.read_csv(
            raw_file,
            compression="gzip",
            header=0,
            index_col=0,
            usecols=use_positions,
            chunksize=1000,
        ):
            totals[selected_cols_array] += chunk.to_numpy(dtype=np.float32, copy=False).sum(axis=0)
            use_rows = [gene for gene in chunk.index if gene in gene_to_row]
            if not use_rows:
                continue
            sub = chunk.loc[use_rows].to_numpy(dtype=np.float32, copy=False)
            rr, cc = np.nonzero(sub)
            if rr.size:
                row_parts.append(np.array([gene_to_row[str(use_rows[i])] for i in rr], dtype=np.int32))
                col_parts.append(selected_cols_array[cc].astype(np.int32))
                data_parts.append(sub[rr, cc].astype(np.float32))

    counts = sparse.coo_matrix(
        (
            np.concatenate(data_parts),
            (np.concatenate(row_parts), np.concatenate(col_parts)),
        ),
        shape=(len(common_genes), selected.shape[0]),
    ).tocsr()
    scale = np.divide(10000.0, totals, out=np.zeros_like(totals), where=totals > 0)
    log_norm = counts @ sparse.diags(scale)
    log_norm.data = np.log1p(log_norm.data)
    mmwrite(str(outdir / "onek1k_reference_logexpr.mtx"), log_norm.astype(np.float32))
    selected[["cell_id", "science14_cell_type", "pool_number", "barcode"]].to_csv(
        outdir / "onek1k_reference_metadata.tsv", sep="\t", index=False
    )


def write_hca_query_matrix(
    h5ad: Path,
    selected: pd.DataFrame,
    common_genes: list[str],
    hca_genes: list[str],
    outdir: Path,
) -> None:
    hca_gene_to_common = {gene: i for i, gene in enumerate(common_genes)}
    hca_col_to_row = np.full(len(hca_genes), -1, dtype=np.int32)
    for col, gene in enumerate(hca_genes):
        if gene in hca_gene_to_common:
            hca_col_to_row[col] = hca_gene_to_common[gene]

    data_parts = []
    row_parts = []
    col_parts = []
    selected_rows = selected["hca_row_index"].to_numpy(dtype=np.int64)
    with h5py.File(h5ad, "r") as handle:
        x = handle["X"]
        data = x["data"]
        indices = x["indices"]
        indptr = x["indptr"]
        for out_col, hca_row in enumerate(selected_rows):
            start = int(indptr[hca_row])
            end = int(indptr[hca_row + 1])
            cols = indices[start:end]
            mapped = hca_col_to_row[cols]
            keep = mapped >= 0
            if keep.any():
                row_parts.append(mapped[keep].astype(np.int32))
                col_parts.append(np.full(int(keep.sum()), out_col, dtype=np.int32))
                data_parts.append(data[start:end][keep].astype(np.float32))

    mat = sparse.coo_matrix(
        (
            np.concatenate(data_parts),
            (np.concatenate(row_parts), np.concatenate(col_parts)),
        ),
        shape=(len(common_genes), selected.shape[0]),
    ).tocsr()
    mmwrite(str(outdir / "hca_healthy_query_logexpr.mtx"), mat.astype(np.float32))
    selected[["cell_id", "full_clustering", "sample_id", "Status", "hca_row_index"]].to_csv(
        outdir / "hca_healthy_query_metadata.tsv", sep="\t", index=False
    )


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    onek1k_genes = read_onek1k_gene_list(args.raw_dir)
    hca_genes = read_hca_gene_list(args.hca_h5ad)
    common_genes = [gene for gene in onek1k_genes if gene in set(hca_genes)]
    pd.DataFrame({"gene": common_genes}).to_csv(args.outdir / "shared_genes.tsv", sep="\t", index=False)

    ref_meta = select_onek1k_cells(args.metadata, args.raw_dir, args.max_ref_cells_per_type, args.seed)
    query_meta = select_hca_healthy_cells(args.hca_h5ad, args.max_query_cells_per_cluster, args.seed)

    write_onek1k_matrix(args.raw_dir, ref_meta, common_genes, args.outdir)
    write_hca_query_matrix(args.hca_h5ad, query_meta, common_genes, hca_genes, args.outdir)

    qc = pd.DataFrame(
        [
            {"metric": "shared_genes", "value": len(common_genes)},
            {"metric": "reference_cells", "value": ref_meta.shape[0]},
            {"metric": "query_cells", "value": query_meta.shape[0]},
            {"metric": "query_clusters", "value": query_meta["full_clustering"].nunique()},
            {"metric": "max_ref_cells_per_type", "value": args.max_ref_cells_per_type},
            {"metric": "max_query_cells_per_cluster", "value": args.max_query_cells_per_cluster},
        ]
    )
    qc.to_csv(args.outdir / "seurat_input_qc.tsv", sep="\t", index=False)
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()
