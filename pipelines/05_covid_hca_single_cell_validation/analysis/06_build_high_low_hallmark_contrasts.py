#!/usr/bin/env python3
"""Build CD8 ET and platelet Hallmark contrasts for SLC44A2 tertile groups."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
RESULTS = PROJECT / "covid_validation/stephenson_hca/result5_expression_context/results"
EXPOSURE = RESULTS / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv"
CD8_SCORES = RESULTS / "cd8_et_covid_d0_hallmark_scores.tsv"
PLATELET_SCORES = RESULTS / "covid_d0_platelet_selected_hallmark_scores.tsv"
CD8_OUT = RESULTS / "cd8et_slc44a2_pseudobulk_mean_hallmark_high_low.tsv"
PLATELET_OUT = RESULTS / "platelet_hallmark_high_low_cd8et_slc44a2_28v28.tsv"


def bh_fdr(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return out


def exposure_groups() -> pd.DataFrame:
    frame = pd.read_csv(EXPOSURE, sep="\t")
    frame = frame.loc[
        frame["cell"].eq("CD8_ET") & frame["Gene"].eq("SLC44A2"),
        ["sample_id", "candidate_mean_expr"],
    ].drop_duplicates("sample_id")
    lower, upper = frame["candidate_mean_expr"].quantile([1 / 3, 2 / 3])
    frame["group"] = "Middle"
    frame.loc[frame["candidate_mean_expr"] <= lower, "group"] = "Low"
    frame.loc[frame["candidate_mean_expr"] >= upper, "group"] = "High"
    return frame


def summarize(
    score_file: Path,
    outcome: str,
    value: str,
    output: Path,
) -> pd.DataFrame:
    frame = pd.read_csv(score_file, sep="\t").merge(
        exposure_groups()[["sample_id", "group"]], on="sample_id", how="inner"
    )
    frame = frame[frame["group"].isin(["Low", "High"])].copy()
    rows = []
    for name, sub in frame.groupby(outcome, sort=False):
        low = sub.loc[sub["group"].eq("Low"), value].dropna().to_numpy(float)
        high = sub.loc[sub["group"].eq("High"), value].dropna().to_numpy(float)
        p = stats.mannwhitneyu(low, high, alternative="two-sided").pvalue
        genes = int(sub["n_genes_present"].max()) if "n_genes_present" in sub else np.nan
        rows.append(
            {
                outcome: name,
                "n_low": len(low),
                "n_high": len(high),
                "median_low": np.median(low),
                "median_high": np.median(high),
                "delta_high_minus_low": np.median(high) - np.median(low),
                "p": p,
                "n_genes_present": genes,
            }
        )
    out = pd.DataFrame(rows)
    out["FDR"] = bh_fdr(out["p"])
    out.to_csv(output, sep="\t", index=False)
    return out


def main() -> None:
    cd8 = summarize(CD8_SCORES, "hallmark", "hallmark_score", CD8_OUT)
    platelet = summarize(
        PLATELET_SCORES,
        "hallmark",
        "platelet_hallmark_score",
        PLATELET_OUT,
    )
    if len(cd8) != 50 or len(platelet) != 9:
        raise RuntimeError(f"Unexpected Hallmark counts: CD8={len(cd8)}, platelet={len(platelet)}")
    print(CD8_OUT)
    print(PLATELET_OUT)


if __name__ == "__main__":
    main()
