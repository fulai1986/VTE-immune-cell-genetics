#!/usr/bin/env python3
"""LIANA all-eligible-donor test for CD8 ET-SLC44A2 communication axes.

This includes COVID D0 donors with sufficient CD8 ET, monocyte, and platelet
cells. LIANA is run separately per donor; donor-level
communication scores are then associated with CD8 ET SLC44A2 mean expression and
positive-cell fraction.
"""

from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import h5py
import liana as li
import numpy as np
import pandas as pd
from scipy import sparse, stats


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
BASE = PROJECT / "covid_validation/stephenson_hca"
H5AD = BASE / "covid_portal.h5ad"
MAP = PROJECT / "covid_validation/hca_onek1k_mapping/results/hca_all_cells_onek1k_compatible_annotation.curated.tsv.gz"
OUT = BASE / "result5_expression_context/results"
PAIR_SCORES = OUT / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv"

MIN_CD8_ET = 100
MIN_MONOCYTE = 100
MIN_PLATELET = 20
GROUPS = ["CD8 ET", "Mono C", "Mono NC", "Platelets"]


def read_var_names(handle: h5py.File) -> list[str]:
    return [
        x.decode("utf-8") if isinstance(x, bytes) else str(x)
        for x in handle["var"]["names"][:]
    ]


def subset_csr(handle: h5py.File, row_idx: np.ndarray) -> sparse.csr_matrix:
    x = handle["X"]
    data = x["data"]
    indices = x["indices"]
    indptr = x["indptr"]
    shape = tuple(x.attrs["shape"])
    n_vars = shape[1]

    chunks_data = []
    chunks_indices = []
    new_indptr = [0]
    for row in row_idx:
        start = int(indptr[row])
        end = int(indptr[row + 1])
        row_data = data[start:end]
        row_indices = indices[start:end]
        chunks_data.append(row_data)
        chunks_indices.append(row_indices)
        new_indptr.append(new_indptr[-1] + len(row_data))

    out_data = np.concatenate(chunks_data).astype(np.float32) if chunks_data else np.array([], dtype=np.float32)
    out_indices = (
        np.concatenate(chunks_indices).astype(np.int32) if chunks_indices else np.array([], dtype=np.int32)
    )
    return sparse.csr_matrix(
        (out_data, out_indices, np.asarray(new_indptr, dtype=np.int64)),
        shape=(len(row_idx), n_vars),
    )


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


def build_meta() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(MAP, sep="\t")
    meta = meta.rename(
        columns={
            "curated_OneK1K_compatible": "OneK1K_compatible",
            "curated_mapping_confidence": "mapping_confidence",
        }
    )
    meta["cell_index"] = np.arange(len(meta), dtype=np.int64)
    meta["comm_group"] = meta["OneK1K_compatible"]
    meta.loc[meta["full_clustering"].eq("Platelets"), "comm_group"] = "Platelets"
    covid = meta[
        meta["Status"].eq("Covid")
        & meta["Collection_Day"].eq("D0")
        & meta["comm_group"].isin(GROUPS)
        & (
            meta["full_clustering"].eq("Platelets")
            | ~meta["mapping_confidence"].isin(["Low-cell", "HCA-specific", "Not applicable"])
        )
    ].copy()
    covid.loc[covid["full_clustering"].eq("Platelets"), "mapping_confidence"] = "HCA-specific"
    counts = covid.groupby(["sample_id", "comm_group"]).size().unstack(fill_value=0)
    counts["Mono_total"] = counts.get("Mono C", 0) + counts.get("Mono NC", 0)
    eligible = counts[
        (counts.get("CD8 ET", 0) >= MIN_CD8_ET)
        & (counts["Mono_total"] >= MIN_MONOCYTE)
        & (counts.get("Platelets", 0) >= MIN_PLATELET)
    ].copy()
    eligible["total_focus"] = eligible.get("CD8 ET", 0) + eligible["Mono_total"] + eligible.get("Platelets", 0)
    covid = covid[covid["sample_id"].isin(eligible.index)].copy()
    return covid.sort_values(["sample_id", "cell_index"]), eligible.sort_values("total_focus", ascending=False)


def run_liana_donor(donor: str, meta: pd.DataFrame) -> pd.DataFrame:
    sub = meta[meta["sample_id"].eq(donor)].copy()
    with h5py.File(H5AD, "r") as handle:
        x = subset_csr(handle, sub["cell_index"].to_numpy(dtype=np.int64))
        var_names = read_var_names(handle)
    if not np.isfinite(x.data).all():
        x.data[~np.isfinite(x.data)] = 0

    adata = ad.AnnData(
        X=x,
        obs=sub[["hca_cell_id", "sample_id", "comm_group", "full_clustering"]].set_index("hca_cell_id"),
        var=pd.DataFrame(index=var_names),
    )
    li.mt.rank_aggregate(
        adata,
        groupby="comm_group",
        resource_name="consensus",
        expr_prop=0.05,
        min_cells=20,
        use_raw=False,
        n_perms=100,
        n_jobs=1,
        inplace=True,
        verbose=False,
    )
    res = adata.uns["liana_res"].copy()
    res.insert(0, "sample_id", donor)
    return res


def add_axis(df: pd.DataFrame) -> pd.DataFrame:
    focus = df.copy()
    focus["axis"] = np.nan
    focus.loc[
        focus["source"].eq("CD8 ET") & focus["target"].isin(["Mono C", "Mono NC"]),
        "axis",
    ] = "CD8_ET_to_monocyte"
    focus.loc[
        focus["source"].isin(["Mono C", "Mono NC"]) & focus["target"].eq("Platelets"),
        "axis",
    ] = "monocyte_to_platelet"
    focus.loc[
        focus["source"].eq("CD8 ET") & focus["target"].eq("Platelets"),
        "axis",
    ] = "CD8_ET_to_platelet"
    return focus.dropna(subset=["axis"]).copy()


def make_axis_scores(focus: pd.DataFrame) -> pd.DataFrame:
    use = focus.copy()
    use["magnitude_rank"] = pd.to_numeric(use["magnitude_rank"], errors="coerce").clip(lower=1e-300)
    use["rank_strength"] = -np.log10(use["magnitude_rank"])
    use["lrscore"] = pd.to_numeric(use["lrscore"], errors="coerce")
    use["lr_means"] = pd.to_numeric(use["lr_means"], errors="coerce")
    use["cellphone_pvals"] = pd.to_numeric(use["cellphone_pvals"], errors="coerce")

    rows = []
    for (sample, axis), sub in use.groupby(["sample_id", "axis"]):
        sig = sub[sub["cellphone_pvals"] < 0.05]
        if sig.empty:
            sig = sub.nsmallest(min(50, len(sub)), "magnitude_rank")
        top = sub.nsmallest(min(50, len(sub)), "magnitude_rank")
        rows.append(
            {
                "sample_id": sample,
                "axis": axis,
                "n_lr_pairs": len(sub),
                "n_cellphone_p_lt_0.05": int((sub["cellphone_pvals"] < 0.05).sum()),
                "mean_lrscore_sig": sig["lrscore"].mean(),
                "mean_lr_means_sig": sig["lr_means"].mean(),
                "mean_rank_strength_sig": sig["rank_strength"].mean(),
                "mean_lrscore_top50": top["lrscore"].mean(),
                "mean_lr_means_top50": top["lr_means"].mean(),
                "mean_rank_strength_top50": top["rank_strength"].mean(),
            }
        )
    return pd.DataFrame(rows)


def load_slc44a2_metrics() -> pd.DataFrame:
    pair = pd.read_csv(PAIR_SCORES, sep="\t")
    return pair[pair["cell"].eq("CD8_ET") & pair["Gene"].eq("SLC44A2")][
        [
            "sample_id",
            "candidate_mean_expr",
            "candidate_pct_expr",
            "candidate_activity",
            "n_target_cells",
            "Sex",
            "Age",
            "Site",
            "total_mapped_source_cells",
        ]
    ].drop_duplicates()


def residualize(y: pd.Series, covars: pd.DataFrame) -> pd.Series:
    x = pd.concat(
        [
            np.log10(covars[["total_mapped_source_cells", "n_target_cells"]].astype(float) + 1),
            pd.get_dummies(covars[["Sex", "Age", "Site"]].astype(str), drop_first=True),
        ],
        axis=1,
    )
    x.insert(0, "const", 1.0)
    x = x.astype(float).to_numpy()
    yy = y.astype(float).to_numpy()
    keep = np.isfinite(yy) & np.isfinite(x).all(axis=1)
    beta, _, _, _ = np.linalg.lstsq(x[keep], yy[keep], rcond=None)
    out = np.full(len(yy), np.nan)
    out[keep] = yy[keep] - x[keep] @ beta
    return pd.Series(out, index=y.index)


def association_tests(axis_scores: pd.DataFrame, metrics: pd.DataFrame) -> pd.DataFrame:
    merged = axis_scores.merge(metrics, on="sample_id", how="inner")
    score_cols = [
        "mean_lrscore_sig",
        "mean_lr_means_sig",
        "mean_rank_strength_sig",
        "mean_lrscore_top50",
        "mean_lr_means_top50",
        "mean_rank_strength_top50",
    ]
    exposures = ["candidate_mean_expr", "candidate_pct_expr"]
    rows = []
    for axis, sub in merged.groupby("axis"):
        for score in score_cols:
            for exposure in exposures:
                dat = sub.dropna(subset=[score, exposure]).copy()
                if len(dat) < 10:
                    continue
                rho, p = stats.spearmanr(dat[exposure], dat[score])
                y_resid = residualize(dat[score], dat)
                x_resid = residualize(dat[exposure], dat)
                ok = y_resid.notna() & x_resid.notna()
                adj_rho, adj_p = stats.spearmanr(x_resid[ok], y_resid[ok])
                rows.append(
                    {
                        "axis": axis,
                        "communication_score": score,
                        "exposure": "mean_expression" if exposure == "candidate_mean_expr" else "positive_fraction",
                        "n": len(dat),
                        "spearman_rho": rho,
                        "spearman_p": p,
                        "adjusted_residual_spearman_rho": adj_rho,
                        "adjusted_residual_spearman_p": adj_p,
                    }
                )
    out = pd.DataFrame(rows)
    out["FDR_spearman"] = bh_fdr(out["spearman_p"])
    out["FDR_adjusted_residual_spearman"] = bh_fdr(out["adjusted_residual_spearman_p"])
    return out.sort_values(["FDR_adjusted_residual_spearman", "FDR_spearman", "axis"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    meta, eligible = build_meta()
    eligible.to_csv(OUT / "liana_all_eligible_donor_cell_counts.tsv", sep="\t")
    print(f"Eligible donors: {len(eligible)}", flush=True)
    print(eligible.head(20).to_string(), flush=True)

    donor_files = []
    raw_dir = OUT / "liana_per_donor"
    raw_dir.mkdir(parents=True, exist_ok=True)
    for i, donor in enumerate(eligible.index, start=1):
        out_file = raw_dir / f"{donor}.liana_rank_aggregate.tsv.gz"
        donor_files.append(out_file)
        if out_file.exists():
            print(f"[{i}/{len(eligible)}] skip existing {donor}", flush=True)
            continue
        print(f"[{i}/{len(eligible)}] running {donor}", flush=True)
        res = run_liana_donor(donor, meta)
        res.to_csv(out_file, sep="\t", index=False)

    raw = pd.concat((pd.read_csv(path, sep="\t") for path in donor_files), ignore_index=True)
    focus = add_axis(raw)
    focus.to_csv(OUT / "liana_all_eligible_cd8et_monocyte_platelet_focus_edges.tsv.gz", sep="\t", index=False)
    axis_scores = make_axis_scores(focus)
    axis_scores.to_csv(OUT / "liana_all_eligible_cd8et_monocyte_platelet_axis_scores.tsv", sep="\t", index=False)
    assoc = association_tests(axis_scores, load_slc44a2_metrics())
    assoc.to_csv(OUT / "liana_cd8et_slc44a2_expression_fraction_axis_associations.tsv", sep="\t", index=False)

    print("\nAssociation results:")
    print(assoc.to_string(index=False))


if __name__ == "__main__":
    main()
