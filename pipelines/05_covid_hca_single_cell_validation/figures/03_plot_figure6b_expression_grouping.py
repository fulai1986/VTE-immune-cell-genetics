#!/usr/bin/env python3
"""Publication-style SLC44A2 CD8 ET donor grouping panel."""

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

PAIR_SCORES = RESULTS / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv"
OUT = FIGURES / "figure6b_slc44a2_expression_grouping_publication"

LOW = "#707A84"
MID = "#C3CAD1"
HIGH = "#B5163A"
GRID = "#E8ECEF"
TEXT = "#111827"
MUTED = "#5B6570"
HEADER = "#EEF1F4"
LOW_BAND = "#ECEFF2"
HIGH_BAND = "#F9DCE3"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 7.8,
            "axes.linewidth": 0.85,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_data() -> pd.DataFrame:
    pair = pd.read_csv(PAIR_SCORES, sep="\t")
    dat = (
        pair[pair["cell"].eq("CD8_ET") & pair["Gene"].eq("SLC44A2")]
        [["sample_id", "candidate_mean_expr"]]
        .drop_duplicates()
        .sort_values("candidate_mean_expr")
        .reset_index(drop=True)
    )
    q_low, q_high = dat["candidate_mean_expr"].quantile([1 / 3, 2 / 3])
    dat["group"] = "Middle"
    dat.loc[dat["candidate_mean_expr"] <= q_low, "group"] = "Low"
    dat.loc[dat["candidate_mean_expr"] >= q_high, "group"] = "High"
    dat["rank"] = np.arange(1, len(dat) + 1)
    return dat


def main() -> None:
    setup_style()
    dat = load_data()
    n_low = int(dat["group"].eq("Low").sum())
    n_high = int(dat["group"].eq("High").sum())
    low_end = dat.index[dat["group"].eq("Low")].max() + 1
    high_start = dat.index[dat["group"].eq("High")].min() + 1

    fig = plt.figure(figsize=(3.75, 2.85))
    ax_bg = fig.add_axes([0, 0, 1, 1])
    ax_bg.axis("off")
    ax_bg.add_patch(
        Rectangle(
            (0.055, 0.055),
            0.89,
            0.88,
            transform=fig.transFigure,
            facecolor="white",
            edgecolor="#D8DEE6",
            linewidth=1.0,
        )
    )
    ax_bg.add_patch(
        Rectangle(
            (0.055, 0.835),
            0.89,
            0.10,
            transform=fig.transFigure,
            facecolor=HEADER,
            edgecolor="none",
        )
    )
    fig.text(
        0.50,
        0.884,
        r"$\it{SLC44A2}$ expression in CD8 ET cells",
        ha="center",
        va="center",
        fontsize=11.2,
        fontweight="bold",
        color=TEXT,
    )

    ax = fig.add_axes([0.20, 0.39, 0.70, 0.36])
    colors = dat["group"].map({"Low": LOW, "Middle": MID, "High": HIGH})
    ax.plot(dat["rank"], dat["candidate_mean_expr"], color="#AAB2BB", lw=0.9, zorder=1)
    ax.scatter(
        dat["rank"],
        dat["candidate_mean_expr"],
        s=18,
        c=colors,
        edgecolor="white",
        linewidth=0.35,
        zorder=2,
    )
    ax.axvline(low_end + 0.5, color="#7E8792", lw=0.9, ls=(0, (4, 3)))
    ax.axvline(high_start - 0.5, color="#7E8792", lw=0.9, ls=(0, (4, 3)))
    ax.set_xlim(0, len(dat) + 2)
    ax.set_ylim(0.15, 0.77)
    ax.set_xticks([1, 20, 40, 60, 80])
    ax.set_yticks([0.2, 0.4, 0.6])
    ax.set_xlabel("")
    ax.set_ylabel("Pseudobulk expression", labelpad=7, fontsize=8.7)
    ax.grid(axis="y", color=GRID, lw=0.85)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.85, color="#222222")

    fig.canvas.draw()
    xlabel_y = 0.285
    band_y = 0.225
    label_y = 0.168
    n_y = 0.112

    def data_x_to_fig(x_value: float) -> float:
        return fig.transFigure.inverted().transform(ax.transData.transform((x_value, ax.get_ylim()[0])))[0]

    low_left = data_x_to_fig(1)
    low_right = data_x_to_fig(low_end)
    high_left = data_x_to_fig(high_start)
    high_right = data_x_to_fig(len(dat))
    low_mid = (low_left + low_right) / 2
    high_mid = (high_left + high_right) / 2
    fig.text(0.55, xlabel_y, "Donors ranked by expression", ha="center", va="center", fontsize=8.9, color=TEXT)
    ax_bg.add_patch(Rectangle((low_left, band_y), low_right - low_left, 0.027, transform=fig.transFigure, facecolor=LOW_BAND, edgecolor="none"))
    ax_bg.add_patch(Rectangle((high_left, band_y), high_right - high_left, 0.027, transform=fig.transFigure, facecolor=HIGH_BAND, edgecolor="none"))
    fig.text(low_mid, label_y, "Low", ha="center", va="center", fontsize=9.2, fontweight="bold", color=TEXT)
    fig.text(low_mid, n_y, f"n = {n_low}", ha="center", va="center", fontsize=8.8, color=TEXT)
    fig.text(high_mid, label_y, "High", ha="center", va="center", fontsize=9.2, fontweight="bold", color=TEXT)
    fig.text(high_mid, n_y, f"n = {n_high}", ha="center", va="center", fontsize=8.8, color=TEXT)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{OUT}.{ext}", dpi=600, bbox_inches="tight")
    print(f"{OUT}.png")


if __name__ == "__main__":
    main()
