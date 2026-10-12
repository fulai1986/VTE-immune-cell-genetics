#!/usr/bin/env python3
"""Platelet MSigDB Hallmark programs for CD8 ET-SLC44A2 contextualization."""

from __future__ import annotations

from collections import defaultdict
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
BASE = PROJECT / "covid_validation/stephenson_hca"
RESULTS = BASE / "result5_expression_context/results"
H5AD = BASE / "covid_portal.h5ad"
MAP = PROJECT / "covid_validation/hca_onek1k_mapping/results/hca_all_cells_onek1k_compatible_annotation.curated.tsv.gz"
GMT = PROJECT / "resources/msigdb/h.all.v2025.1.Hs.symbols.gmt"

PAIR_SCORES = RESULTS / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv"
LIANA_FOCUS = RESULTS / "liana_all_eligible_cd8et_monocyte_platelet_focus_edges.tsv.gz"
LIANA_AXIS = RESULTS / "liana_all_eligible_cd8et_monocyte_platelet_axis_scores.tsv"

OUT_EXPR = RESULTS / "covid_d0_platelet_selected_hallmark_gene_expression.tsv.gz"
OUT_SCORES = RESULTS / "covid_d0_platelet_selected_hallmark_scores.tsv"
OUT_SLC44A2_ASSOC = RESULTS / "cd8et_slc44a2_positive_fraction_platelet_hallmark_associations.tsv"
OUT_COMM_ASSOC = RESULTS / "platelet_hallmark_vs_liana_communication_associations.tsv"

MIN_PLATELET_CELLS = 20

SELECTED_HALLMARKS = [
    "HALLMARK_COAGULATION",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_HYPOXIA",
    "HALLMARK_APOPTOSIS",
]


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


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce").astype(float)
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / sd


def bh_fdr(pvalues: pd.Series) -> np.ndarray:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * ok.sum() / np.arange(1, ok.sum() + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked, 1.0)
    return q


def load_hallmarks() -> dict[str, list[str]]:
    hallmarks: dict[str, list[str]] = {}
    with GMT.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] in SELECTED_HALLMARKS:
                hallmarks[fields[0]] = fields[2:]
    missing = sorted(set(SELECTED_HALLMARKS) - set(hallmarks))
    if missing:
        raise RuntimeError(f"Missing Hallmark gene sets: {missing}")
    return hallmarks


def load_platelet_meta() -> pd.DataFrame:
    mapped = pd.read_csv(MAP, sep="\t")
    mapped["cell_index"] = np.arange(len(mapped), dtype=np.int64)
    with h5py.File(H5AD, "r") as handle:
        obs = handle["obs"]
        clinical = pd.DataFrame(
            {
                "hca_cell_id": decode_h5ad_column(obs, "_index"),
                "Sex": decode_h5ad_column(obs, "Sex"),
                "Age": decode_h5ad_column(obs, "Age"),
                "Site": decode_h5ad_column(obs, "Site"),
            }
        )
    meta = mapped.merge(clinical, on="hca_cell_id", how="left")
    meta = meta[
        meta["Status"].eq("Covid")
        & meta["Collection_Day"].eq("D0")
        & meta["full_clustering"].eq("Platelets")
    ].copy()
    return meta.sort_values("cell_index")


def extract_platelet_expression(hallmarks: dict[str, list[str]]) -> pd.DataFrame:
    genes = sorted({gene for gene_list in hallmarks.values() for gene in gene_list})
    meta = load_platelet_meta()
    with h5py.File(H5AD, "r") as handle:
        var_names = [
            x.decode("utf-8") if isinstance(x, bytes) else str(x)
            for x in handle["var"]["names"][:]
        ]
        gene_idx = {gene: var_names.index(gene) for gene in genes if gene in var_names}
        target_cols = set(gene_idx.values())
        idx_to_gene = {idx: gene for gene, idx in gene_idx.items()}

        row_idx = meta["cell_index"].to_numpy(dtype=np.int64)
        row_sample = meta["sample_id"].astype(str).to_numpy()
        n_cells = defaultdict(int)
        sums = defaultdict(float)
        positives = defaultdict(int)
        x = handle["X"]
        data = x["data"]
        indices = x["indices"]
        indptr = x["indptr"]

        target_col_array = np.array(sorted(target_cols), dtype=indices.dtype)
        for i, row in enumerate(row_idx):
            sample = row_sample[i]
            n_cells[sample] += 1
            start = int(indptr[row])
            end = int(indptr[row + 1])
            cols = indices[start:end]
            if len(cols) == 0:
                continue
            keep = np.isin(cols, target_col_array)
            if not keep.any():
                continue
            vals = data[start:end][keep]
            for col, val in zip(cols[keep], vals):
                gene = idx_to_gene[int(col)]
                key = (sample, gene)
                sums[key] += float(val)
                if float(val) > 0:
                    positives[key] += 1

    rows = []
    for sample, ncells in n_cells.items():
        for gene in gene_idx:
            key = (sample, gene)
            rows.append(
                {
                    "sample_id": sample,
                    "gene": gene,
                    "platelet_n_cells": ncells,
                    "mean_expr": sums[key] / ncells,
                    "pct_expr": positives[key] / ncells,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_EXPR, sep="\t", index=False, compression="gzip")
    return out


def compute_scores(expr: pd.DataFrame, hallmarks: dict[str, list[str]]) -> pd.DataFrame:
    expr = expr[expr["platelet_n_cells"] >= MIN_PLATELET_CELLS].copy()
    expr["gene_z"] = expr.groupby("gene")["mean_expr"].transform(zscore)
    rows = []
    for hallmark, genes in hallmarks.items():
        present = sorted(set(genes).intersection(expr["gene"]))
        sub = expr[expr["gene"].isin(present)].copy()
        score = (
            sub.groupby("sample_id", as_index=False)
            .agg(
                platelet_hallmark_score=("gene_z", "mean"),
                platelet_n_cells=("platelet_n_cells", "first"),
                n_genes_present=("gene", "nunique"),
                genes_present=("gene", lambda x: ",".join(sorted(set(map(str, x))))),
            )
            .assign(hallmark=hallmark)
        )
        rows.append(score)
    scores = pd.concat(rows, ignore_index=True)
    clinical = load_platelet_meta()[["sample_id", "Sex", "Age", "Site"]].drop_duplicates()
    scores = scores.merge(clinical, on="sample_id", how="left")
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    return scores


def slc44a2_fraction() -> pd.DataFrame:
    pair = pd.read_csv(PAIR_SCORES, sep="\t")
    metrics = pair[pair["cell"].eq("CD8_ET") & pair["Gene"].eq("SLC44A2")][
        ["sample_id", "candidate_pct_expr"]
    ].drop_duplicates()
    metrics = metrics.rename(columns={"candidate_pct_expr": "cd8et_slc44a2_positive_fraction"})
    q1, q2 = metrics["cd8et_slc44a2_positive_fraction"].quantile([1 / 3, 2 / 3])
    metrics["slc44a2_fraction_group"] = "Middle"
    metrics.loc[metrics["cd8et_slc44a2_positive_fraction"] <= q1, "slc44a2_fraction_group"] = "Low"
    metrics.loc[metrics["cd8et_slc44a2_positive_fraction"] >= q2, "slc44a2_fraction_group"] = "High"
    return metrics


def summarize_slc44a2_associations(scores: pd.DataFrame) -> pd.DataFrame:
    dat = scores.merge(slc44a2_fraction(), on="sample_id", how="inner")
    rows = []
    for hallmark, sub in dat.groupby("hallmark", sort=True):
        rho, p = stats.spearmanr(
            sub["cd8et_slc44a2_positive_fraction"], sub["platelet_hallmark_score"]
        )
        low = sub[sub["slc44a2_fraction_group"].eq("Low")]["platelet_hallmark_score"]
        high = sub[sub["slc44a2_fraction_group"].eq("High")]["platelet_hallmark_score"]
        _, mw_p = stats.mannwhitneyu(high, low, alternative="two-sided")
        rows.append(
            {
                "hallmark": hallmark,
                "n": len(sub),
                "n_low": len(low),
                "n_high": len(high),
                "n_genes_present": int(sub["n_genes_present"].iloc[0]),
                "rho_vs_slc44a2_fraction": float(rho),
                "p_vs_slc44a2_fraction": float(p),
                "low_median": float(low.median()),
                "high_median": float(high.median()),
                "delta_high_minus_low": float(high.median() - low.median()),
                "high_vs_low_p": float(mw_p),
            }
        )
    out = pd.DataFrame(rows)
    out["FDR_vs_slc44a2_fraction"] = bh_fdr(out["p_vs_slc44a2_fraction"])
    out["FDR_high_vs_low"] = bh_fdr(out["high_vs_low_p"])
    out = out.sort_values(["FDR_vs_slc44a2_fraction", "p_vs_slc44a2_fraction"])
    out.to_csv(OUT_SLC44A2_ASSOC, sep="\t", index=False)
    return out


def liana_axis_with_summed_rank() -> pd.DataFrame:
    axis = pd.read_csv(LIANA_AXIS, sep="\t")
    focus = pd.read_csv(LIANA_FOCUS, sep="\t")
    focus["cellphone_pvals"] = pd.to_numeric(focus["cellphone_pvals"], errors="coerce")
    focus["magnitude_rank"] = pd.to_numeric(focus["magnitude_rank"], errors="coerce").clip(lower=1e-300)
    focus["rank_strength"] = -np.log10(focus["magnitude_rank"])
    summed = (
        focus[focus["cellphone_pvals"] < 0.05]
        .groupby(["sample_id", "axis"], as_index=False)["rank_strength"]
        .sum()
        .rename(columns={"rank_strength": "sum_rank_strength_sig"})
    )
    out = axis.merge(summed, on=["sample_id", "axis"], how="left")
    out["sum_rank_strength_sig"] = out["sum_rank_strength_sig"].fillna(0.0)
    return out


def summarize_communication_associations(scores: pd.DataFrame) -> pd.DataFrame:
    axis = liana_axis_with_summed_rank()
    dat = axis.merge(scores, on="sample_id", how="inner")
    rows = []
    for (hallmark, axis_name), sub in dat.groupby(["hallmark", "axis"], sort=True):
        for metric, label in [
            ("n_cellphone_p_lt_0.05", "Significant LR pairs"),
            ("sum_rank_strength_sig", "Summed rank evidence"),
        ]:
            rho, p = stats.spearmanr(sub["platelet_hallmark_score"], sub[metric])
            rows.append(
                {
                    "hallmark": hallmark,
                    "axis": axis_name,
                    "metric": metric,
                    "metric_label": label,
                    "n": len(sub),
                    "rho_platelet_hallmark_vs_communication": float(rho),
                    "p_platelet_hallmark_vs_communication": float(p),
                }
            )
    out = pd.DataFrame(rows)
    out["FDR"] = bh_fdr(out["p_platelet_hallmark_vs_communication"])
    out = out.sort_values(["FDR", "p_platelet_hallmark_vs_communication"])
    out.to_csv(OUT_COMM_ASSOC, sep="\t", index=False)
    return out


def main() -> None:
    hallmarks = load_hallmarks()
    expr = extract_platelet_expression(hallmarks)
    scores = compute_scores(expr, hallmarks)
    slc = summarize_slc44a2_associations(scores)
    comm = summarize_communication_associations(scores)

    print("Platelet Hallmark association with CD8 ET SLC44A2-positive fraction:")
    print(slc.to_string(index=False))
    print("\nPlatelet Hallmark association with LIANA communication axes:")
    print(comm.head(30).to_string(index=False))
    print("\nOutputs:")
    for path in [OUT_EXPR, OUT_SCORES, OUT_SLC44A2_ASSOC, OUT_COMM_ASSOC]:
        print(path)


if __name__ == "__main__":
    main()
