#!/usr/bin/env python3
"""Figure 6E-style six-row panel for severity-adjusted continuous LIANA associations."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
HERE = PROJECT / "manuscript_results/current_sensitivity"
HERE.mkdir(parents=True, exist_ok=True)
INPUT = HERE / "liana_adjusted_display.tsv"
OUT = HERE / "Figure_S4_liana_severity_adjusted_associations"
ORDER = ["CD8_ET_to_monocyte", "CD8_ET_to_platelet", "monocyte_to_platelet"]
AXIS = {
    "CD8_ET_to_monocyte": "CD8 ET → monocyte",
    "CD8_ET_to_platelet": "CD8 ET → platelet",
    "monocyte_to_platelet": "Monocyte → platelet",
}
METRIC = ["Significant LR pairs", "Summed LIANA rank evidence"]
TEXT = "#111827"
RED = "#C2183F"
GRID = "#E2E7ED"


def main() -> None:
    dat = pd.read_csv(INPUT, sep="\t")
    assert len(dat) == 6 and dat.n_donors.eq(54).all()
    dat = dat.set_index(["axis", "metric"]).loc[[(a, m) for a in ORDER for m in METRIC]].reset_index()
    assert dat.fdr_six.lt(.05).all()

    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8.5, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(8.6, 3.9), facecolor="white")
    bg = fig.add_axes([0, 0, 1, 1]); bg.axis("off")
    bg.add_patch(Rectangle((.018, .06), .964, .81, transform=fig.transFigure,
                           facecolor="white", edgecolor="#D9DEE6", lw=.8))
    bg.add_patch(Rectangle((.018, .79), .964, .08, transform=fig.transFigure,
                           facecolor="#F1F3F5", edgecolor="none"))

    left = fig.add_axes([.04, .23, .35, .49])
    forest = fig.add_axes([.40, .23, .30, .49])
    right = fig.add_axes([.72, .23, .24, .49])
    for ax in [left, right]:
        ax.set_xlim(0, 1); ax.set_ylim(-.55, 5.55); ax.axis("off")
    forest.set_xlim(-.08, .78); forest.set_ylim(-.55, 5.55)
    forest.axvline(0, color="#7B8490", lw=.9, ls=(0, (4, 2)), zorder=1)
    forest.grid(axis="x", color=GRID, lw=.7, zorder=0)
    forest.set_yticks([])
    forest.set_xticks([0, .2, .4, .6], ["0", "0.2", "0.4", "0.6"])
    forest.set_xlabel("Severity-adjusted partial Spearman ρ\n(bootstrap 95% CI)", fontsize=8.5)
    for side in ["top", "left", "right"]:
        forest.spines[side].set_visible(False)

    fig.text(.04, .83, "Communication axis", va="center", weight="bold", color=TEXT)
    fig.text(.235, .83, "Metric", va="center", weight="bold", color=TEXT)
    fig.text(.535, .83, "Positive association", ha="center", va="center", color=TEXT, fontsize=8)
    fig.text(.72, .83, "Adjusted ρ", va="center", weight="bold", color=TEXT)
    fig.text(.87, .83, "FDR", va="center", weight="bold", color=TEXT)
    fig.text(.018, .97,
             "LIANA communication: severity-adjusted associations",
             va="top", weight="bold", fontsize=10.5, color=TEXT)
    fig.text(.018, .925, "With continuous donor-mean CD8 ET SLC44A2 expression",
             va="top", fontsize=8.9, color="#4B5563")

    for j, row in enumerate(dat.itertuples(index=False)):
        y = 5 - j
        if j in (0, 4):
            p0 = fig.transFigure.inverted().transform(left.transData.transform((0, y - .47)))[1]
            p1 = fig.transFigure.inverted().transform(left.transData.transform((0, y + .47)))[1]
            bg.add_patch(Rectangle((.023, p0), .954, p1-p0, transform=fig.transFigure,
                                   facecolor="#F8F9FA", edgecolor="none", zorder=-1))
        if j % 2 == 0:
            left.text(0, y - .5, AXIS[row.axis], va="center", ha="left", fontsize=8.3,
                      weight="bold", color=TEXT)
        left.text(.56, y, "Number of LR pairs" if row.metric == METRIC[0] else row.metric,
                  va="center", ha="left", fontsize=7.4, color=TEXT)
        forest.hlines(y, row.ci_low, row.ci_high, color=RED, lw=1.6, zorder=2)
        forest.vlines([row.ci_low, row.ci_high], y-.065, y+.065, color=RED, lw=1)
        forest.scatter(row.adjusted_rho, y, s=42, color=RED, edgecolor="white", lw=.5, zorder=3)
        right.text(.02, y, f"ρ = {row.adjusted_rho:.3f}", va="center", ha="left", fontsize=8, color=TEXT)
        right.text(.63, y, f"{row.fdr_six:.2e}", va="center", ha="left", fontsize=8,
                   color=RED, weight="bold")

    for sep in [1.5, 3.5]:
        for ax in [left, forest, right]:
            ax.axhline(sep, color="#DDE3EA", lw=.75, zorder=0)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT.with_suffix("." + ext), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(OUT.with_suffix(".png"))


if __name__ == "__main__":
    main()
