#!/usr/bin/env python3
"""Continuous donor expression and collection-day categorical severity sensitivity.

Reads existing donor-level scores/LIANA outputs; creates derived files only.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
RESULTS = PROJECT / "covid_validation/stephenson_hca/result5_expression_context/results"
OUT = PROJECT / "manuscript_results/current_sensitivity"
AXIS_FILE = RESULTS / "liana_all_eligible_cd8et_monocyte_platelet_axis_scores.tsv"
EDGES_FILE = RESULTS / "liana_all_eligible_cd8et_monocyte_platelet_focus_edges.tsv.gz"
EXPRESSION_FILE = RESULTS / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv"
METADATA_FILE = RESULTS / "covid_d0_same_cell_program_sample_metadata.tsv"
AXES = ["CD8_ET_to_monocyte", "CD8_ET_to_platelet", "monocyte_to_platelet"]
METRICS = ["n_cellphone_p_lt_0.05", "sum_rank_strength_sig"]
SEVERITY = ["Asymptomatic", "Mild", "Moderate", "Severe", "Critical"]
METRIC_LABEL = {METRICS[0]: "Number of significant LR pairs", METRICS[1]: "Summed LIANA rank evidence"}

def partial_spearman(frame: pd.DataFrame, metric: str, coding: str) -> tuple[float, float, float]:
    """Residual correlation of midranks; t P accounts for severity model dimension."""
    rx = stats.rankdata(frame.candidate_mean_expr.to_numpy(float))
    ry = stats.rankdata(frame[metric].to_numpy(float))
    if coding == "categorical":
        codes = pd.Categorical(frame.severity, categories=SEVERITY).codes
        controls = np.column_stack([np.ones(len(frame))] + [(codes == j).astype(float) for j in range(1, 5)])
    else:
        controls = np.column_stack([np.ones(len(frame)), frame.severity_rank.to_numpy(float)])
    ex = rx - controls @ np.linalg.lstsq(controls, rx, rcond=None)[0]
    ey = ry - controls @ np.linalg.lstsq(controls, ry, rcond=None)[0]
    rho = float(np.corrcoef(ex, ey)[0, 1])
    df = len(frame) - np.linalg.matrix_rank(controls) - 1
    t = rho * np.sqrt(df / (1 - rho * rho))
    p = float(2 * stats.t.sf(abs(t), df))
    return rho, p, float(df)


def bootstrap_ci(frame: pd.DataFrame, metric: str, seed: int, n_boot: int = 5000) -> tuple[float, float]:
    """Bootstrap donors within D0 severity strata, retaining observed stratum sizes."""
    rng = np.random.default_rng(seed)
    strata = [np.flatnonzero(frame.severity.to_numpy() == x) for x in SEVERITY]
    values = []
    for _ in range(n_boot):
        picks = np.concatenate([rng.choice(ix, size=len(ix), replace=True) for ix in strata])
        b = frame.iloc[picks]
        rho, _, _ = partial_spearman(b, metric, "categorical")
        if np.isfinite(rho):
            values.append(rho)
    return tuple(np.quantile(values, [0.025, 0.975]))


def make_input() -> pd.DataFrame:
    axes = pd.read_csv(AXIS_FILE, sep="\t")
    edges = pd.read_csv(EDGES_FILE, sep="\t")
    edges["cellphone_pvals"] = pd.to_numeric(edges.cellphone_pvals, errors="coerce")
    edges["magnitude_rank"] = pd.to_numeric(edges.magnitude_rank, errors="coerce").clip(lower=1e-300)
    summed = (edges[edges.cellphone_pvals.lt(0.05)]
              .assign(rank_strength=lambda x: -np.log10(x.magnitude_rank))
              .groupby(["sample_id", "axis"], as_index=False).rank_strength.sum()
              .rename(columns={"rank_strength": "sum_rank_strength_sig"}))
    axes = axes.merge(summed, on=["sample_id", "axis"], how="left", validate="one_to_one")
    axes["sum_rank_strength_sig"] = axes.sum_rank_strength_sig.fillna(0)
    expression = pd.read_csv(EXPRESSION_FILE, sep="\t")
    expression = expression.loc[expression.cell.eq("CD8_ET") & expression.Gene.eq("SLC44A2"), ["sample_id", "candidate_mean_expr"]].drop_duplicates("sample_id")
    metadata = pd.read_csv(METADATA_FILE, sep="\t")
    assert not metadata.sample_id.duplicated().any()
    frame = axes.merge(expression, on="sample_id", how="inner", validate="many_to_one")
    frame = frame.merge(metadata[["sample_id", "Status_on_day_collection_summary"]], on="sample_id", how="left", validate="many_to_one")
    frame["severity"] = frame.Status_on_day_collection_summary.astype("string").str.strip()
    assert frame.severity.isin(SEVERITY).all()
    frame["severity_rank"] = frame.severity.map({x: i for i, x in enumerate(SEVERITY)}).astype(int)
    assert frame.sample_id.nunique() == 54
    assert frame.groupby("axis").sample_id.nunique().eq(54).all()
    assert not frame.duplicated(["sample_id", "axis"]).any()
    return frame



def main():
    OUT.mkdir(parents=True, exist_ok=True)
    liana = make_input()
    expression = pd.read_csv(EXPRESSION_FILE, sep="\t")
    expression = expression.loc[expression.cell.eq("CD8_ET") & expression.Gene.eq("SLC44A2"), ["sample_id", "candidate_mean_expr"]].drop_duplicates()
    assert not expression.sample_id.duplicated().any()
    metadata = pd.read_csv(METADATA_FILE, sep="\t")
    assert not metadata.sample_id.duplicated().any()
    metadata["severity"] = metadata.Status_on_day_collection_summary.astype("string").str.strip()
    exposure = expression.merge(metadata[["sample_id", "severity"]], on="sample_id", validate="one_to_one")
    assert len(exposure) == 84 and exposure.severity.isin(SEVERITY).all()
    targets = []
    for family, filename, column in [
        ("CD8 ET Hallmark", "cd8_et_covid_d0_hallmark_scores.tsv", "hallmark_score"),
        ("Platelet Hallmark", "covid_d0_platelet_selected_hallmark_scores.tsv", "platelet_hallmark_score")]:
        data = pd.read_csv(RESULTS / filename, sep="\t").merge(exposure, on="sample_id", validate="many_to_one")
        for name, sub in data.groupby("hallmark", sort=True):
            targets.append((family, name, sub.sort_values("sample_id"), column, None))
    for mi, metric in enumerate(METRICS):
        for ai, axis in enumerate(AXES):
            # Preserve the original donor order for exact seeded bootstrap replay.
            sub = liana.loc[liana.axis.eq(axis)].copy()
            targets.append(("LIANA communication", axis + " | " + METRIC_LABEL[metric], sub, metric, 20260930 + 10 * mi + ai))
    rows = []
    for family, readout, sub, col, seed in targets:
        assert not sub.sample_id.duplicated().any()
        assert sub.severity.isin(SEVERITY).all()
        assert np.isfinite(sub[["candidate_mean_expr", col]].to_numpy(float)).all()
        rho, p = stats.spearmanr(sub.candidate_mean_expr, sub[col])
        adjusted, adjusted_p, df = partial_spearman(sub, col, "categorical")
        assert df == len(sub) - 6
        z = np.arctanh(np.clip(rho, -.999999, .999999))
        width = stats.norm.ppf(.975) / np.sqrt(len(sub) - 3)
        low, high = bootstrap_ci(sub, col, seed) if seed is not None else (np.nan, np.nan)
        rows.append(dict(family=family, readout=readout, n_donors=len(sub), unadjusted_rho=rho,
            unadjusted_p=p, unadjusted_ci_low=np.tanh(z-width), unadjusted_ci_high=np.tanh(z+width),
            severity_adjusted_rho=adjusted, severity_adjusted_p=adjusted_p,
            adjusted_ci_low=low, adjusted_ci_high=high, residual_df=int(df)))
    out = pd.DataFrame(rows)
    expected = {"CD8 ET Hallmark": (50, 84), "LIANA communication": (6, 54), "Platelet Hallmark": (9, 76)}
    for family, sub in out.groupby("family"):
        tests, donors = expected[family]
        assert len(sub) == tests and sub.n_donors.eq(donors).all()
        out.loc[sub.index, "unadjusted_fdr"] = stats.false_discovery_control(sub.unadjusted_p.to_numpy())
        out.loc[sub.index, "severity_adjusted_fdr"] = stats.false_discovery_control(sub.severity_adjusted_p.to_numpy())
        out.loc[sub.index, "fdr_tests"] = tests
    out["fdr_tests"] = out.fdr_tests.astype(int)
    out.to_csv(OUT / "continuous_severity_results.tsv", sep="\t", index=False)
    li = out.loc[out.family.eq("LIANA communication")].copy()
    li[["axis", "metric"]] = li.readout.str.split(" | ", regex=False, expand=True)
    li["metric"] = li.metric.replace({"Number of significant LR pairs": "Significant LR pairs"})
    li = li.rename(columns={"severity_adjusted_rho":"adjusted_rho", "severity_adjusted_fdr":"fdr_six", "adjusted_ci_low":"ci_low", "adjusted_ci_high":"ci_high"})
    li.to_csv(OUT / "liana_adjusted_display.tsv", sep="\t", index=False)
    print(out.groupby("family").size().to_string())

if __name__ == "__main__":
    main()
