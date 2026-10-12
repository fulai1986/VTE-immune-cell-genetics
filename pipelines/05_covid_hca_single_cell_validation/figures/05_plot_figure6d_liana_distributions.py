#!/usr/bin/env python3
"""Standalone donor-level LIANA distribution panel for Figure 6D."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
BASE = PROJECT / "covid_validation/stephenson_hca/result5_expression_context"
RESULTS = BASE / "results"
FIGURES = BASE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

DONOR_METRICS = RESULTS / "figure6_liana_donor_metrics.tsv"
CONTRASTS = RESULTS / "figure5c_cd8et_monocyte_platelet_bootstrap_contrasts.tsv"
OUT = FIGURES / "figure6d_liana_donor_level_communication_scores"

AXIS_ORDER = ["CD8_ET_to_monocyte", "CD8_ET_to_platelet", "monocyte_to_platelet"]
AXIS_LABEL = {
    "CD8_ET_to_monocyte": "CD8 ET \u2192 monocyte",
    "CD8_ET_to_platelet": "CD8 ET \u2192 platelet",
    "monocyte_to_platelet": "Monocyte \u2192 platelet",
}
METRIC_ORDER = ["n_cellphone_p_lt_0.05", "sum_rank_strength_sig"]
METRIC_LABEL = {
    "n_cellphone_p_lt_0.05": "Number of LR pairs",
    "sum_rank_strength_sig": "Summed LIANA rank evidence",
}

LOW = "#7A838D"
HIGH = "#C2183F"
LOW_FILL = "#EEF1F4"
HIGH_FILL = "#F9CDD6"
TEXT = "#111827"
MUTED = "#65707C"
GRID = "#E5EAF0"
CARD_BG = "#F7F8FA"
HEADER = "#E9EDF2"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 8.3,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 7.3,
            "axes.linewidth": 0.85,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def ensure_donor_metrics() -> pd.DataFrame:
    if not DONOR_METRICS.exists():
        raise FileNotFoundError(
            f"Missing {DONOR_METRICS}. Run analysis/07_summarize_liana_high_low_contrasts.py first."
        )
    return pd.read_csv(DONOR_METRICS, sep="\t")


def format_fdr(x: float) -> str:
    return f"FDR = {x:.2e}"


def draw_violin(ax: plt.Axes, donor: pd.DataFrame, axis_name: str, metric: str) -> None:
    vals = []
    for group in ["Low", "High"]:
        vals.append(
            donor.loc[
                donor["axis"].eq(axis_name) & donor["slc44a2_group"].eq(group),
                metric,
            ]
            .astype(float)
            .to_numpy()
        )

    parts = ax.violinplot(
        vals,
        positions=[0, 1],
        widths=0.74,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body, color in zip(parts["bodies"], [LOW, HIGH]):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.17)

    for i, (arr, color) in enumerate(zip(vals, [LOW, HIGH])):
        rng = np.random.default_rng(2026 + i)
        jitter = rng.normal(0, 0.045, len(arr))
        ax.scatter(np.full(len(arr), i) + jitter, arr, s=11, color=color, alpha=0.55, linewidth=0, zorder=2)
        q1, med, q3 = np.percentile(arr, [25, 50, 75])
        ax.plot([i, i], [q1, q3], color="#222222", lw=1.0, zorder=3)
        ax.scatter([i], [med], marker="s", s=24, facecolor="white", edgecolor="#222222", linewidth=0.85, zorder=4)

    ymax = max(np.nanmax(vals[0]), np.nanmax(vals[1]))
    ymin = min(np.nanmin(vals[0]), np.nanmin(vals[1]))
    pad = max((ymax - ymin) * 0.24, ymax * 0.08, 1.0)
    ax.set_ylim(max(0, ymin - pad * 0.25), ymax + pad)

    y_bracket = ymax + pad * 0.38
    ax.plot([0, 0, 1, 1], [y_bracket * 0.985, y_bracket, y_bracket, y_bracket * 0.985], color="#4B5563", lw=0.75)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Low", "High"], fontsize=8)
    ax.grid(axis="y", color=GRID, lw=0.65)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=2.5, width=0.75)


def draw_card(fig: plt.Figure, slot, donor: pd.DataFrame, contrasts: pd.DataFrame, axis_name: str) -> None:
    card = fig.add_subplot(slot)
    card.axis("off")
    card.add_patch(Rectangle((0, 0), 1, 1, transform=card.transAxes, facecolor=CARD_BG, edgecolor="#E1E5EA", lw=0.85))
    card.add_patch(Rectangle((0, 0.875), 1, 0.125, transform=card.transAxes, facecolor=HEADER, edgecolor="none"))
    card.text(0.5, 0.938, AXIS_LABEL[axis_name], ha="center", va="center", fontsize=10.2, fontweight="bold", color=TEXT)

    ax_left = card.inset_axes([0.105, 0.125, 0.375, 0.585])
    ax_right = card.inset_axes([0.575, 0.125, 0.375, 0.585])
    for ax, metric in [(ax_left, METRIC_ORDER[0]), (ax_right, METRIC_ORDER[1])]:
        row = contrasts.loc[contrasts["axis"].eq(axis_name) & contrasts["metric"].eq(metric)].iloc[0]
        draw_violin(ax, donor, axis_name, metric)

    for x, metric in [(0.287, METRIC_ORDER[0]), (0.767, METRIC_ORDER[1])]:
        row = contrasts.loc[contrasts["axis"].eq(axis_name) & contrasts["metric"].eq(metric)].iloc[0]
        sig = float(row["FDR"]) < 0.05
        card.text(
            x,
            0.805,
            METRIC_LABEL[metric],
            ha="center",
            va="center",
            fontsize=8.0,
            color=TEXT,
        )
        card.text(
            x,
            0.752,
            format_fdr(float(row["FDR"])),
            ha="center",
            va="center",
            fontsize=7.8,
            color=HIGH if sig else MUTED,
            fontweight="bold" if sig else "normal",
        )

    ax_left.set_ylabel("")
    ax_right.set_ylabel("")


def main() -> None:
    setup_style()
    donor = ensure_donor_metrics()
    contrasts = pd.read_csv(CONTRASTS, sep="\t")
    contrasts = contrasts[contrasts["metric"].isin(METRIC_ORDER)].copy()

    fig = plt.figure(figsize=(12.0, 3.05))
    gs = fig.add_gridspec(1, 3, left=0.035, right=0.985, top=0.82, bottom=0.11, wspace=0.025)
    fig.text(0.035, 0.91, "Donor-level communication scores (LIANA)", ha="left", va="center", fontsize=13.2, fontweight="bold")

    for i, axis_name in enumerate(AXIS_ORDER):
        draw_card(fig, gs[0, i], donor, contrasts, axis_name)

    fig.text(
        0.018,
        0.435,
        "Donor-level score",
        rotation=90,
        ha="center",
        va="center",
        fontsize=9.0,
        color=TEXT,
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{OUT}.{ext}", dpi=600, bbox_inches="tight")
    print(f"{OUT}.png")


if __name__ == "__main__":
    main()
