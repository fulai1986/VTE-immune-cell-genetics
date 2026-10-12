#!/usr/bin/env python3
"""Summarize per-donor LIANA metrics and High-minus-Low contrasts."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
RESULTS = PROJECT / "covid_validation/stephenson_hca/result5_expression_context/results"
EXPOSURE = RESULTS / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv"
LIANA_EDGES = RESULTS / "liana_all_eligible_cd8et_monocyte_platelet_focus_edges.tsv.gz"
DONOR_OUT = RESULTS / "figure6_liana_donor_metrics.tsv"
CONTRAST_OUT = RESULTS / "figure5c_cd8et_monocyte_platelet_bootstrap_contrasts.tsv"

AXES = ["CD8_ET_to_monocyte", "CD8_ET_to_platelet", "monocyte_to_platelet"]
AXIS_LABELS = {
    "CD8_ET_to_monocyte": "CD8 ET -> monocyte",
    "CD8_ET_to_platelet": "CD8 ET -> platelet",
    "monocyte_to_platelet": "monocyte -> platelet",
}
METRICS = {
    "n_cellphone_p_lt_0.05": "Significant LR pairs",
    "sum_rank_strength_sig": "Summed rank evidence",
}


def bh_fdr(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    return out


def groups() -> pd.DataFrame:
    frame = pd.read_csv(EXPOSURE, sep="\t")
    frame = frame.loc[
        frame["cell"].eq("CD8_ET") & frame["Gene"].eq("SLC44A2"),
        ["sample_id", "candidate_mean_expr"],
    ].drop_duplicates("sample_id")
    lower, upper = frame["candidate_mean_expr"].quantile([1 / 3, 2 / 3])
    frame["slc44a2_group"] = "Middle"
    frame.loc[frame["candidate_mean_expr"] <= lower, "slc44a2_group"] = "Low"
    frame.loc[frame["candidate_mean_expr"] >= upper, "slc44a2_group"] = "High"
    return frame


def donor_metrics() -> pd.DataFrame:
    edges = pd.read_csv(LIANA_EDGES, sep="\t")
    edges = edges[edges["axis"].isin(AXES)].copy()
    edges["cellphone_pvals"] = pd.to_numeric(edges["cellphone_pvals"], errors="coerce")
    edges["magnitude_rank"] = pd.to_numeric(edges["magnitude_rank"], errors="coerce").clip(1e-300)
    edges["rank_strength"] = -np.log10(edges["magnitude_rank"])
    edges["is_significant"] = edges["cellphone_pvals"] < 0.05

    sample_ids = sorted(set(groups()["sample_id"]).intersection(edges["sample_id"]))
    grid = pd.MultiIndex.from_product(
        [sample_ids, AXES], names=["sample_id", "axis"]
    ).to_frame(index=False)
    count = (
        edges.groupby(["sample_id", "axis"], as_index=False)["is_significant"]
        .sum()
        .rename(columns={"is_significant": "n_cellphone_p_lt_0.05"})
    )
    strength = (
        edges[edges["is_significant"]]
        .groupby(["sample_id", "axis"], as_index=False)["rank_strength"]
        .sum()
        .rename(columns={"rank_strength": "sum_rank_strength_sig"})
    )
    out = grid.merge(count, on=["sample_id", "axis"], how="left").merge(
        strength, on=["sample_id", "axis"], how="left"
    )
    out[list(METRICS)] = out[list(METRICS)].fillna(0)
    out = out.merge(groups(), on="sample_id", how="left")
    out = out[out["slc44a2_group"].isin(["Low", "High"])].copy()
    out.to_csv(DONOR_OUT, sep="\t", index=False)
    return out


def bootstrap_median_difference(low: np.ndarray, high: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    draws = 20_000
    low_index = rng.integers(0, len(low), size=(draws, len(low)))
    high_index = rng.integers(0, len(high), size=(draws, len(high)))
    differences = np.median(high[high_index], axis=1) - np.median(low[low_index], axis=1)
    return tuple(np.quantile(differences, [0.025, 0.975]))


def main() -> None:
    donor = donor_metrics()
    rows = []
    seed = 20260909
    for metric, metric_label in METRICS.items():
        for axis in AXES:
            sub = donor[donor["axis"].eq(axis)]
            low = sub.loc[sub["slc44a2_group"].eq("Low"), metric].to_numpy(float)
            high = sub.loc[sub["slc44a2_group"].eq("High"), metric].to_numpy(float)
            ci_low, ci_high = bootstrap_median_difference(low, high, seed)
            seed += 1
            rows.append(
                {
                    "axis": axis,
                    "axis_label": AXIS_LABELS[axis],
                    "metric": metric,
                    "metric_label": metric_label,
                    "n_low": len(low),
                    "n_high": len(high),
                    "median_low": np.median(low),
                    "median_high": np.median(high),
                    "delta_high_minus_low": np.median(high) - np.median(low),
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "p": stats.mannwhitneyu(low, high, alternative="two-sided").pvalue,
                }
            )
    out = pd.DataFrame(rows)
    out["FDR"] = bh_fdr(out["p"])
    out.to_csv(CONTRAST_OUT, sep="\t", index=False)
    print(DONOR_OUT)
    print(CONTRAST_OUT)


if __name__ == "__main__":
    main()
