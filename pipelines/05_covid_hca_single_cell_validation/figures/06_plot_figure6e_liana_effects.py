#!/usr/bin/env python3
"""Standalone raw effect-size forest/table panel for LIANA contrasts."""

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

CONTRASTS = RESULTS / "figure5c_cd8et_monocyte_platelet_bootstrap_contrasts.tsv"
OUT_DATA = RESULTS / "figure6e_liana_raw_effect_size.tsv"
OUT = FIGURES / "figure6e_liana_effect_size"

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

TEXT = "#111827"
MUTED = "#747D87"
HIGH = "#C2183F"
GRID = "#E2E7ED"
HEADER = "#F1F3F5"
ROW_SHADE = "#F8F9FA"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 9.0,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.2,
            "axes.linewidth": 0.85,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def build_plot_data() -> pd.DataFrame:
    contrasts = pd.read_csv(CONTRASTS, sep="\t")
    rows = []
    for axis in AXIS_ORDER:
        for metric in METRIC_ORDER:
            row = contrasts.loc[contrasts["axis"].eq(axis) & contrasts["metric"].eq(metric)].iloc[0]
            rows.append(
                {
                    "axis": axis,
                    "axis_label": AXIS_LABEL[axis],
                    "metric": metric,
                    "metric_label": METRIC_LABEL[metric],
                    "n_low": int(row["n_low"]),
                    "n_high": int(row["n_high"]),
                    "delta_high_minus_low": float(row["delta_high_minus_low"]),
                    "ci_low": float(row["ci_low"]),
                    "ci_high": float(row["ci_high"]),
                    "FDR": float(row["FDR"]),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DATA, sep="\t", index=False)
    return out


def fmt_effect(row) -> str:
    return f"\u0394 = {row.delta_high_minus_low:+.1f}"


def fmt_fdr(x: float) -> str:
    return f"{x:.2e}"


def main() -> None:
    setup_style()
    dat = build_plot_data()
    n = len(dat)
    y = np.arange(n)[::-1]

    fig = plt.figure(figsize=(8.25, 3.05))
    ax_bg = fig.add_axes([0, 0, 1, 1])
    ax_bg.axis("off")
    ax_bg.add_patch(Rectangle((0.025, 0.105), 0.95, 0.765, transform=fig.transFigure, facecolor="white", edgecolor="#D9DEE6", lw=0.8))
    ax_bg.add_patch(Rectangle((0.025, 0.792), 0.95, 0.078, transform=fig.transFigure, facecolor=HEADER, edgecolor="none"))

    ax_left = fig.add_axes([0.045, 0.195, 0.335, 0.545])
    ax = fig.add_axes([0.392, 0.195, 0.305, 0.545])
    ax_right = fig.add_axes([0.720, 0.195, 0.225, 0.545])

    ax_left.set_xlim(0, 1)
    ax_left.set_ylim(-0.6, n - 0.4)
    ax_left.axis("off")
    ax_right.set_xlim(0, 1)
    ax_right.set_ylim(-0.6, n - 0.4)
    ax_right.axis("off")

    ax.set_xlim(-42, 65)
    ax.set_ylim(-0.6, n - 0.4)
    ax.axvline(0, color="#4B5563", lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax.grid(axis="x", color=GRID, lw=0.75, zorder=0)
    ax.set_xlabel("High \u2212 Low difference (bootstrap 95% CI)", fontsize=8.8, labelpad=6)
    ax.set_yticks([])
    ax.set_xticks([-40, -20, 0, 20, 40, 60])
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)

    fig.text(0.045, 0.832, "Communication axis", ha="left", va="center", fontsize=8.6, fontweight="bold", color=TEXT)
    fig.text(0.220, 0.832, "Metric", ha="left", va="center", fontsize=8.6, fontweight="bold", color=TEXT)
    fig.text(0.433, 0.832, "Higher in Low", ha="center", va="center", fontsize=7.5, color=TEXT)
    fig.text(0.634, 0.832, "Higher in High", ha="center", va="center", fontsize=7.5, color=TEXT)
    fig.text(0.725, 0.832, "\u0394 (High \u2212 Low)", ha="left", va="center", fontsize=8.6, fontweight="bold", color=TEXT)
    fig.text(0.875, 0.832, "FDR", ha="left", va="center", fontsize=8.6, fontweight="bold", color=TEXT)
    arrow_y = n - 0.23
    ax.annotate("", xy=(-35, arrow_y), xytext=(-5, arrow_y), arrowprops=dict(arrowstyle="<-", lw=0.7, color=TEXT), clip_on=False)
    ax.annotate("", xy=(35, arrow_y), xytext=(5, arrow_y), arrowprops=dict(arrowstyle="->", lw=0.7, color=TEXT), clip_on=False)

    # Light row shading by communication axis.
    for block_start in [0, 4]:
        top_y = y[block_start] + 0.46
        bottom_y = y[min(block_start + 1, n - 1)] - 0.46
        y0_fig = fig.transFigure.inverted().transform(ax.transData.transform((0, bottom_y)))[1]
        y1_fig = fig.transFigure.inverted().transform(ax.transData.transform((0, top_y)))[1]
        ax_bg.add_patch(Rectangle((0.02, y0_fig), 0.96, y1_fig - y0_fig, transform=fig.transFigure, facecolor=ROW_SHADE, edgecolor="none", zorder=-1))

    last_axis = None
    for yi, row in zip(y, dat.itertuples(index=False)):
        sig = float(row.FDR) < 0.05
        color = HIGH if sig else MUTED
        if row.axis != last_axis:
            ax_left.text(0.00, yi + 0.23, row.axis_label, ha="left", va="center", fontsize=8.6, fontweight="bold", color=TEXT)
            last_axis = row.axis
        ax_left.text(0.50, yi, row.metric_label, ha="left", va="center", fontsize=7.6, color=TEXT, clip_on=False)
        ax.hlines(yi, row.ci_low, row.ci_high, color=color, lw=1.25, zorder=2)
        ax.scatter(row.delta_high_minus_low, yi, s=36, color=color, edgecolor="white", linewidth=0.55, zorder=3)
        ax_right.text(0.02, yi, fmt_effect(row), ha="left", va="center", fontsize=8.0, color=TEXT)
        ax_right.text(
            0.66,
            yi,
            fmt_fdr(float(row.FDR)),
            ha="left",
            va="center",
            fontsize=8.0,
            color=HIGH if sig else MUTED,
            fontweight="bold" if sig else "normal",
        )

    for sep in [1.5, 3.5]:
        for a in [ax_left, ax, ax_right]:
            a.axhline(sep, color="#DDE3EA", lw=0.75, zorder=0)

    fig.text(0.025, 0.955, "Effect sizes of communication changes (High vs Low)", ha="left", va="top", fontsize=11.6, fontweight="bold", color=TEXT)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{OUT}.{ext}", dpi=600, bbox_inches="tight")
    print(f"{OUT}.png")
    print(OUT_DATA)


if __name__ == "__main__":
    main()
