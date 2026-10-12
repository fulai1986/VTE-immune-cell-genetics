#!/usr/bin/env python3
"""COVID-only same-cell program associations for Tier 1-3 VTE cell-gene pairs.

Unit of analysis is donor/sample-level pseudobulk within each
OneK1K-compatible cell type. COVID D0 cells are assigned by the fixed
Healthy-derived lineage-aware mapping; no single-cell-level correlation is used.
"""

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
H5AD = BASE / "covid_portal.h5ad"
MAP = PROJECT / "covid_validation/hca_onek1k_mapping/results/hca_all_cells_onek1k_compatible_annotation.curated.tsv.gz"
TIER = PROJECT / "results/smr_analysis/qc/finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf_beta_or.tsv"
OUT = BASE / "result5_expression_context/results"

MIN_CELLS_PER_SAMPLE_CELL = 20

PROGRAMS = {
    "cytotoxicity": ["NKG7", "GNLY", "GZMB", "GZMH", "GZMA", "PRF1", "CST7", "CTSW"],
    "T_cell_activation": ["CD69", "IL2RA", "TNFRSF9", "ICOS", "CD40LG", "FOS", "JUN", "JUNB", "NFKBIA"],
    "IFN_response": ["IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "MX2", "OAS1", "OAS2", "IFI6", "IRF7", "STAT1", "GBP1"],
    "inflammatory_response": ["TNF", "IL1B", "NFKBIA", "CXCL8", "CCL3", "CCL4", "CCL5", "JUNB", "FOS"],
    "migration_chemokine": ["CCL3", "CCL4", "CCL5", "CXCR3", "CCR7", "SELL", "S1PR1", "ITGAL", "ITGB2"],
    "exhaustion": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1"],
    "proliferation": ["MKI67", "TOP2A", "STMN1", "TUBA1B", "TYMS", "PCNA", "MCM5"],
    "metabolism": ["GAPDH", "LDHA", "PGK1", "ENO1", "PKM", "SLC2A3", "HIF1A", "MT-CO1", "MT-CO2", "MT-ATP6"],
}


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


def load_pairs() -> pd.DataFrame:
    pairs = pd.read_csv(TIER, sep="\t")
    pairs = pairs[pairs["tier"].isin(["Tier 1", "Tier 2", "Tier 3"])].copy()
    pairs["cell_label"] = pairs["cell"].str.replace("_", " ", regex=False)
    return pairs[
        [
            "tier",
            "cell",
            "cell_label",
            "Gene",
            "b_SMR_finngen",
            "OR_SMR_finngen",
            "smr_q_SMR_cell_BH",
            "smr_p_HEIDI",
            "PP.H4.abf",
        ]
    ].drop_duplicates().sort_values(["tier", "cell", "Gene"])


def load_meta(source_labels: set[str]) -> pd.DataFrame:
    mapped = pd.read_csv(MAP, sep="\t")
    mapped = mapped.rename(
        columns={
            "curated_OneK1K_compatible": "OneK1K_compatible",
            "curated_mapping_confidence": "mapping_confidence",
        }
    )
    mapped["cell_index"] = np.arange(len(mapped), dtype=np.int64)
    with h5py.File(H5AD, "r") as handle:
        obs = handle["obs"]
        clinical = pd.DataFrame(
            {
                "hca_cell_id": decode_h5ad_column(obs, "_index"),
                "Sex": decode_h5ad_column(obs, "Sex"),
                "Age": decode_h5ad_column(obs, "Age"),
                "Site": decode_h5ad_column(obs, "Site"),
                "Status_on_day_collection_summary": decode_h5ad_column(obs, "Status_on_day_collection_summary"),
                "Worst_Clinical_Status": decode_h5ad_column(obs, "Worst_Clinical_Status"),
            }
        )
    meta = mapped.merge(clinical, on="hca_cell_id", how="left")
    use = meta[
        meta["Status"].eq("Covid")
        & meta["Collection_Day"].eq("D0")
        & meta["OneK1K_compatible"].isin(source_labels)
        & ~meta["mapping_confidence"].isin(["Low-cell", "HCA-specific", "Not applicable"])
    ].copy()
    return use.sort_values("cell_index")


def extract_expression(meta: pd.DataFrame, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    with h5py.File(H5AD, "r") as handle:
        var_names = [
            x.decode("utf-8") if isinstance(x, bytes) else str(x)
            for x in handle["var"]["names"][:]
        ]
        gene_idx = {gene: var_names.index(gene) for gene in genes if gene in var_names}
        missing = sorted(set(genes) - set(gene_idx))
        if missing:
            print("Missing genes skipped:", ",".join(missing))
        target_cols = set(gene_idx.values())
        idx_to_gene = {idx: gene for gene, idx in gene_idx.items()}

        n_cells = defaultdict(int)
        sums = defaultdict(float)
        positives = defaultdict(int)
        row_idx = meta["cell_index"].to_numpy(dtype=np.int64)
        row_sample = meta["sample_id"].astype(str).to_numpy()
        row_cell = meta["OneK1K_compatible"].astype(str).to_numpy()

        x = handle["X"]
        data = x["data"]
        indices = x["indices"]
        indptr = x["indptr"]
        for i, row in enumerate(row_idx):
            sample = row_sample[i]
            cell = row_cell[i]
            n_cells[(sample, cell)] += 1
            start = int(indptr[row])
            end = int(indptr[row + 1])
            cols = indices[start:end]
            if len(cols) == 0:
                continue
            keep = np.isin(cols, list(target_cols))
            if not keep.any():
                continue
            for col, val in zip(cols[keep], data[start:end][keep]):
                gene = idx_to_gene[int(col)]
                key = (sample, cell, gene)
                sums[key] += float(val)
                if float(val) > 0:
                    positives[key] += 1

    rows = []
    for (sample, cell), nc in n_cells.items():
        if nc < MIN_CELLS_PER_SAMPLE_CELL:
            continue
        for gene in gene_idx:
            key = (sample, cell, gene)
            rows.append(
                {
                    "sample_id": sample,
                    "cell_label": cell,
                    "gene": gene,
                    "n_target_cells": nc,
                    "mean_expr": sums[key] / nc,
                    "pct_expr": positives[key] / nc,
                }
            )
    expr = pd.DataFrame(rows)
    sample_total = (
        meta.drop_duplicates("cell_index")
        .groupby("sample_id")
        .size()
        .rename("total_mapped_source_cells")
        .reset_index()
    )
    sample_meta = (
        meta.drop_duplicates("sample_id")[
            ["sample_id", "Sex", "Age", "Site", "Status_on_day_collection_summary", "Worst_Clinical_Status"]
        ]
        .merge(sample_total, on="sample_id", how="left")
    )
    return expr, sample_meta


def make_candidate_activity(expr: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pair in pairs.itertuples(index=False):
        sub = expr[(expr["cell_label"].eq(pair.cell_label)) & (expr["gene"].eq(pair.Gene))].copy()
        if sub.empty:
            continue
        sub["candidate_mean_z"] = zscore(sub["mean_expr"])
        sub["candidate_pct_z"] = zscore(sub["pct_expr"])
        sub["candidate_activity"] = (sub["candidate_mean_z"] + sub["candidate_pct_z"]) / 2
        for row in sub.itertuples(index=False):
            rows.append(
                {
                    "tier": pair.tier,
                    "cell": pair.cell,
                    "cell_label": pair.cell_label,
                    "Gene": pair.Gene,
                    "sample_id": row.sample_id,
                    "n_target_cells": row.n_target_cells,
                    "candidate_mean_expr": row.mean_expr,
                    "candidate_pct_expr": row.pct_expr,
                    "candidate_activity": row.candidate_activity,
                    "b_SMR_finngen": pair.b_SMR_finngen,
                    "OR_SMR_finngen": pair.OR_SMR_finngen,
                    "smr_q_SMR_cell_BH": pair.smr_q_SMR_cell_BH,
                    "smr_p_HEIDI": pair.smr_p_HEIDI,
                    "PP.H4.abf": getattr(pair, "_8"),
                }
            )
    return pd.DataFrame(rows)


def make_program_scores(expr: pd.DataFrame, pair_activity: pd.DataFrame) -> pd.DataFrame:
    out = pair_activity.copy()
    program_frames = []
    for cell_label, cell_expr in expr.groupby("cell_label"):
        mat = cell_expr.pivot_table(index="sample_id", columns="gene", values="mean_expr", aggfunc="mean")
        for program, genes in PROGRAMS.items():
            present = [gene for gene in genes if gene in mat.columns]
            if len(present) < 3:
                continue
            zmat = mat[present].copy()
            for gene in zmat.columns:
                zmat[gene] = zscore(zmat[gene])
            tmp = zmat.mean(axis=1, skipna=True).rename("program_score").reset_index()
            tmp["cell_label"] = cell_label
            tmp["program"] = program
            tmp["n_program_genes_present"] = len(present)
            tmp["program_genes_present"] = ",".join(present)
            program_frames.append(tmp)
    programs = pd.concat(program_frames, ignore_index=True)
    out = out.merge(programs, on=["sample_id", "cell_label"], how="left")
    return out


def ols_adjusted(sub: pd.DataFrame) -> tuple[float, float, float, int, float]:
    dat = sub[
        [
            "program_score",
            "candidate_activity",
            "Sex",
            "Age",
            "Site",
            "total_mapped_source_cells",
            "n_target_cells",
        ]
    ].replace([np.inf, -np.inf], np.nan).dropna()
    if len(dat) < 20 or dat["candidate_activity"].std(ddof=1) <= 0:
        return np.nan, np.nan, np.nan, len(dat), np.nan
    design = pd.concat(
        [
            pd.Series(1.0, index=dat.index, name="Intercept"),
            dat[["candidate_activity", "total_mapped_source_cells", "n_target_cells"]],
            pd.get_dummies(dat[["Sex", "Age", "Site"]].astype(str), drop_first=True, dtype=float),
        ],
        axis=1,
    )
    design = design.loc[:, design.std(axis=0).fillna(0).gt(0) | (design.columns == "Intercept")]
    x = design.to_numpy(dtype=float)
    y = dat["program_score"].to_numpy(dtype=float)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ beta
    rank = np.linalg.matrix_rank(x)
    df_resid = len(y) - rank
    if df_resid <= 0:
        return np.nan, np.nan, np.nan, len(dat), np.nan
    sigma2 = float((resid @ resid) / df_resid)
    cov_beta = sigma2 * np.linalg.pinv(x.T @ x)
    term_idx = list(design.columns).index("candidate_activity")
    b = float(beta[term_idx])
    se = float(np.sqrt(cov_beta[term_idx, term_idx]))
    tval = b / se if se > 0 else np.nan
    p = float(2 * stats.t.sf(abs(tval), df_resid)) if np.isfinite(tval) else np.nan
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = float(1 - (resid @ resid) / ss_tot) if ss_tot > 0 else np.nan
    return b, se, p, len(dat), r2


def summarize_associations(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, sub in scores.groupby(["tier", "cell", "cell_label", "Gene", "program"], sort=False):
        tier, cell, cell_label, gene, program = key
        dat = sub[["candidate_activity", "program_score"]].dropna()
        if len(dat) >= 20 and dat["candidate_activity"].nunique() > 1 and dat["program_score"].nunique() > 1:
            rho, sp = stats.spearmanr(dat["candidate_activity"], dat["program_score"])
        else:
            rho, sp = np.nan, np.nan
        b, se, p, n_model, r2 = ols_adjusted(sub)
        rows.append(
            {
                "tier": tier,
                "cell": cell,
                "cell_label": cell_label,
                "Gene": gene,
                "program": program,
                "n_spearman": len(dat),
                "spearman_rho": rho,
                "spearman_p": sp,
                "n_adjusted_model": n_model,
                "adjusted_beta_candidate_activity": b,
                "adjusted_se": se,
                "adjusted_p": p,
                "adjusted_r_squared": r2,
                "n_program_genes_present": sub["n_program_genes_present"].dropna().iloc[0]
                if sub["n_program_genes_present"].notna().any()
                else np.nan,
                "program_genes_present": sub["program_genes_present"].dropna().iloc[0]
                if sub["program_genes_present"].notna().any()
                else "",
            }
        )
    out = pd.DataFrame(rows)
    out["spearman_q_BH_all"] = bh_fdr(out["spearman_p"])
    out["adjusted_q_BH_all"] = bh_fdr(out["adjusted_p"])
    for program in sorted(out["program"].dropna().unique()):
        idx = out.index[out["program"].eq(program)]
        out.loc[idx, "spearman_q_BH_within_program"] = bh_fdr(out.loc[idx, "spearman_p"])
        out.loc[idx, "adjusted_q_BH_within_program"] = bh_fdr(out.loc[idx, "adjusted_p"])
    return out.sort_values(["adjusted_p", "spearman_p"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pairs = load_pairs()
    source_labels = set(pairs["cell_label"])
    genes = sorted(set(pairs["Gene"]) | {gene for genes in PROGRAMS.values() for gene in genes})
    meta = load_meta(source_labels)
    expr, sample_meta = extract_expression(meta, genes)
    expr.to_csv(OUT / "covid_d0_tier1_3_same_cell_gene_expression.curated_mapping.tsv.gz", sep="\t", index=False, compression="gzip")
    sample_meta.to_csv(OUT / "covid_d0_same_cell_program_sample_metadata.tsv", sep="\t", index=False)

    pair_activity = make_candidate_activity(expr, pairs)
    scores = make_program_scores(expr, pair_activity).merge(sample_meta, on="sample_id", how="left")
    scores.to_csv(OUT / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv", sep="\t", index=False)

    assoc = summarize_associations(scores)
    assoc.to_csv(OUT / "tier1_3_pairs_covid_d0_same_cell_program_associations.curated_mapping.tsv", sep="\t", index=False)

    print("Same-cell score rows:", scores.shape)
    print("Same-cell associations:", assoc.shape)
    print("\nTop adjusted same-cell associations")
    print(assoc.head(30).to_string(index=False))
    print("\nTop Spearman same-cell associations")
    print(assoc.sort_values(["spearman_p", "adjusted_p"]).head(30).to_string(index=False))


if __name__ == "__main__":
    main()
