#!/usr/bin/env python3
"""CD8 ET-SLC44A2 COVID-only Hallmark pathway associations.

The unit of analysis is donor/sample-level pseudobulk. For each COVID D0
sample, the script computes CD8 ET Hallmark scores from normalized h5ad X
values, then tests each Hallmark score against CD8 ET SLC44A2 mean expression
and SLC44A2-positive fraction separately.
"""

from __future__ import annotations

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
GMT = PROJECT / "resources/msigdb/h.all.v2025.1.Hs.symbols.gmt"
RESULT_DIR = BASE / "result5_expression_context/results"

PAIR_SCORES = RESULT_DIR / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv"
OUT_ASSOC = RESULT_DIR / "cd8_et_slc44a2_covid_d0_hallmark_expression_vs_fraction.tsv"
OUT_SCORE = RESULT_DIR / "cd8_et_covid_d0_hallmark_scores.tsv"
OUT_CONCORDANCE = RESULT_DIR / "cd8_et_slc44a2_covid_d0_hallmark_expression_fraction_concordance.tsv"

MIN_CELLS_PER_SAMPLE = 20
MIN_HALLMARK_GENES = 10


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


def read_hallmark_gmt(path: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].replace("HALLMARK_", "")
            gene_sets[name] = sorted(set(parts[2:]))
    return gene_sets


def zscore_columns(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame.astype(float)
    means = values.mean(axis=0)
    sds = values.std(axis=0, ddof=1)
    keep = sds.replace(0, np.nan).notna()
    return (values.loc[:, keep] - means.loc[keep]) / sds.loc[keep]


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


def load_cd8_et_meta() -> pd.DataFrame:
    mapped = pd.read_csv(MAP, sep="\t")
    mapped = mapped.rename(
        columns={
            "curated_OneK1K_compatible": "OneK1K_compatible",
            "curated_mapping_confidence": "mapping_confidence",
        }
    )
    mapped["cell_index"] = np.arange(len(mapped), dtype=np.int64)
    use = mapped[
        mapped["Status"].eq("Covid")
        & mapped["Collection_Day"].eq("D0")
        & mapped["OneK1K_compatible"].eq("CD8 ET")
        & ~mapped["mapping_confidence"].isin(["Low-cell", "HCA-specific", "Not applicable"])
    ].copy()
    return use.sort_values("cell_index")


def extract_sample_gene_means(meta: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    sample_ids = sorted(meta["sample_id"].astype(str).unique())
    sample_to_i = {sample: i for i, sample in enumerate(sample_ids)}

    with h5py.File(H5AD, "r") as handle:
        var_names = [
            x.decode("utf-8") if isinstance(x, bytes) else str(x)
            for x in handle["var"]["names"][:]
        ]
        gene_to_var = {gene: i for i, gene in enumerate(var_names)}
        genes_present = [gene for gene in genes if gene in gene_to_var]
        var_to_gene_i = np.full(len(var_names), -1, dtype=np.int32)
        for j, gene in enumerate(genes_present):
            var_to_gene_i[gene_to_var[gene]] = j

        sums = np.zeros((len(sample_ids), len(genes_present)), dtype=np.float64)
        n_cells = np.zeros(len(sample_ids), dtype=np.int64)

        x = handle["X"]
        data = x["data"]
        indices = x["indices"]
        indptr = x["indptr"]

        row_idx = meta["cell_index"].to_numpy(dtype=np.int64)
        row_sample_i = meta["sample_id"].astype(str).map(sample_to_i).to_numpy(dtype=np.int32)

        for row, si in zip(row_idx, row_sample_i):
            n_cells[si] += 1
            start = int(indptr[row])
            end = int(indptr[row + 1])
            cols = indices[start:end]
            mapped_cols = var_to_gene_i[cols]
            keep = mapped_cols >= 0
            if keep.any():
                np.add.at(sums[si], mapped_cols[keep], data[start:end][keep].astype(float))

    keep_samples = n_cells >= MIN_CELLS_PER_SAMPLE
    means = pd.DataFrame(
        sums[keep_samples] / n_cells[keep_samples, None],
        index=np.array(sample_ids)[keep_samples],
        columns=genes_present,
    )
    means.index.name = "sample_id"
    return means


def make_hallmark_scores(gene_means: pd.DataFrame, hallmark: dict[str, list[str]]) -> pd.DataFrame:
    gene_z = zscore_columns(gene_means)
    rows = []
    for name, genes in hallmark.items():
        present = [gene for gene in genes if gene in gene_z.columns]
        if len(present) < MIN_HALLMARK_GENES:
            continue
        score = gene_z[present].mean(axis=1)
        rows.append(
            pd.DataFrame(
                {
                    "sample_id": score.index,
                    "hallmark": name,
                    "hallmark_score": score.to_numpy(),
                    "n_genes_present": len(present),
                    "genes_present": ",".join(present),
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def fit_association(dat: pd.DataFrame, exposure: str) -> dict[str, float]:
    cols = [
        "hallmark_score",
        exposure,
        "Sex",
        "Age",
        "Site",
        "total_mapped_source_cells",
        "n_target_cells",
    ]
    use = dat.dropna(subset=cols).copy()
    use[f"{exposure}_z"] = stats.zscore(use[exposure].astype(float), nan_policy="omit")
    use["log_total_mapped_source_cells"] = np.log10(use["total_mapped_source_cells"].astype(float) + 1)
    use["log_n_target_cells"] = np.log10(use["n_target_cells"].astype(float) + 1)

    y = use["hallmark_score"].astype(float)
    design = pd.concat(
        [
            use[[f"{exposure}_z", "log_total_mapped_source_cells", "log_n_target_cells"]],
            pd.get_dummies(use[["Sex", "Age", "Site"]].astype(str), drop_first=True),
        ],
        axis=1,
    )
    design.insert(0, "const", 1.0)
    design = design.astype(float)
    xmat = design.to_numpy(dtype=float)
    yvec = y.to_numpy(dtype=float)
    beta, _, rank, _ = np.linalg.lstsq(xmat, yvec, rcond=None)
    resid = yvec - xmat @ beta
    df_resid = max(len(yvec) - rank, 1)
    sigma2 = float((resid @ resid) / df_resid)
    xtx_inv = np.linalg.pinv(xmat.T @ xmat)
    se = np.sqrt(np.diag(xtx_inv) * sigma2)
    tvals = beta / se
    pvals = 2 * stats.t.sf(np.abs(tvals), df=df_resid)
    coef = pd.Series(beta, index=design.columns)
    coef_se = pd.Series(se, index=design.columns)
    coef_t = pd.Series(tvals, index=design.columns)
    coef_p = pd.Series(pvals, index=design.columns)
    rho, rho_p = stats.spearmanr(use[exposure].astype(float), y)
    term = f"{exposure}_z"
    return {
        "n": float(len(use)),
        "spearman_rho": float(rho),
        "spearman_p": float(rho_p),
        "beta_z": float(coef[term]),
        "se": float(coef_se[term]),
        "t": float(coef_t[term]),
        "p": float(coef_p[term]),
    }


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    hallmark = read_hallmark_gmt(GMT)
    hallmark_genes = sorted(set(gene for genes in hallmark.values() for gene in genes))

    meta = load_cd8_et_meta()
    gene_means = extract_sample_gene_means(meta, hallmark_genes)
    scores = make_hallmark_scores(gene_means, hallmark)
    scores.to_csv(OUT_SCORE, sep="\t", index=False)

    pair_scores = pd.read_csv(PAIR_SCORES, sep="\t")
    candidate = pair_scores[
        pair_scores["cell"].eq("CD8_ET") & pair_scores["Gene"].eq("SLC44A2")
    ][
        [
            "sample_id",
            "candidate_mean_expr",
            "candidate_pct_expr",
            "n_target_cells",
            "Sex",
            "Age",
            "Site",
            "total_mapped_source_cells",
        ]
    ].drop_duplicates()

    merged = scores.merge(candidate, on="sample_id", how="inner")
    rows = []
    for hallmark_name, sub in merged.groupby("hallmark", sort=True):
        for exposure in ["candidate_mean_expr", "candidate_pct_expr"]:
            result = fit_association(sub, exposure)
            result.update(
                {
                    "cell": "CD8_ET",
                    "cell_label": "CD8 ET",
                    "Gene": "SLC44A2",
                    "hallmark": hallmark_name,
                    "exposure": "mean_expression"
                    if exposure == "candidate_mean_expr"
                    else "positive_fraction",
                    "n_genes_present": int(sub["n_genes_present"].iloc[0]),
                }
            )
            rows.append(result)

    assoc = pd.DataFrame(rows)
    assoc["FDR_within_exposure"] = assoc.groupby("exposure")["p"].transform(lambda x: bh_fdr(x))
    assoc["Spearman_FDR_within_exposure"] = assoc.groupby("exposure")["spearman_p"].transform(
        lambda x: bh_fdr(x)
    )
    assoc = assoc.sort_values(["exposure", "FDR_within_exposure", "p", "hallmark"])
    assoc.to_csv(OUT_ASSOC, sep="\t", index=False)

    wide = assoc.pivot(index="hallmark", columns="exposure")
    concordance = pd.DataFrame(
        {
            "hallmark": wide.index,
            "beta_z.mean_expression": wide[("beta_z", "mean_expression")],
            "p.mean_expression": wide[("p", "mean_expression")],
            "FDR.mean_expression": wide[("FDR_within_exposure", "mean_expression")],
            "beta_z.positive_fraction": wide[("beta_z", "positive_fraction")],
            "p.positive_fraction": wide[("p", "positive_fraction")],
            "FDR.positive_fraction": wide[("FDR_within_exposure", "positive_fraction")],
        }
    ).reset_index(drop=True)
    concordance["same_direction_beta"] = np.sign(concordance["beta_z.mean_expression"]) == np.sign(
        concordance["beta_z.positive_fraction"]
    )
    concordance["both_nominal_p_lt_0.05"] = (concordance["p.mean_expression"] < 0.05) & (
        concordance["p.positive_fraction"] < 0.05
    )
    concordance["both_fdr_lt_0.05"] = (concordance["FDR.mean_expression"] < 0.05) & (
        concordance["FDR.positive_fraction"] < 0.05
    )
    concordance["min_fdr"] = concordance[["FDR.mean_expression", "FDR.positive_fraction"]].max(axis=1)
    concordance = concordance.sort_values(
        ["both_fdr_lt_0.05", "both_nominal_p_lt_0.05", "min_fdr"],
        ascending=[False, False, True],
    )
    concordance.to_csv(OUT_CONCORDANCE, sep="\t", index=False)

    print(f"Hallmark sets: {len(hallmark)}")
    print(f"Samples scored: {gene_means.shape[0]}")
    print(f"Hallmark genes present: {gene_means.shape[1]}")
    print(OUT_ASSOC)
    print(OUT_CONCORDANCE)
    print(concordance.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
