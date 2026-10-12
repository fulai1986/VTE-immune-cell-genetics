#!/usr/bin/env python3
"""Figure 6A: HCA detectability of replicated Tier 1 signals."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
BASE = PROJECT / "covid_validation/stephenson_hca/result5_expression_context"
RESULTS = BASE / "results"
FIGURES = BASE / "figures"
DATA = RESULTS / "figure5a_replicated_tier1_hca_expression_qc_plot_data.tsv"
OUT = FIGURES / "figure6a_expression_detectability_qc"

PAIR_ORDER = ["CD8 ET - SLC44A2", "CD8 ET - FADS1", "NK - FADS1"]
PAIR_LABELS = {
    "CD8 ET - SLC44A2": r"CD8 ET - $\it{SLC44A2}$",
    "CD8 ET - FADS1": r"CD8 ET - $\it{FADS1}$",
    "NK - FADS1": r"NK - $\it{FADS1}$",
}
ACCENT = "#9B1B30"
GREY = "#737B84"
PALE_ACCENT = "#F7E8ED"
TEXT = "#20242A"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 8.2,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 9.5,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_data() -> pd.DataFrame:
    df = pd.read_csv(DATA, sep="\t")
    df["pair"] = pd.Categorical(df["pair"], categories=PAIR_ORDER, ordered=True)
    return df.sort_values("pair").reset_index(drop=True)


def draw_interval_axis(
    axes: tuple[plt.Axes, plt.Axes],
    df: pd.DataFrame,
    y: np.ndarray,
    median_col: str,
    q1_col: str,
    q3_col: str,
    xlims: tuple[tuple[float, float], tuple[float, float]],
    xticks: tuple[list[float], list[float]],
    xlabel: str,
    value_fmt: str,
    value_dx: float,
) -> None:
    for ax, xlim, ticks in zip(axes, xlims, xticks):
        for yi, row in zip(y, df.itertuples(index=False)):
            pair = str(row.pair)
            is_lead = pair == "CD8 ET - SLC44A2"
            color = ACCENT if is_lead else GREY
            q1 = getattr(row, q1_col)
            q3 = getattr(row, q3_col)
            med = getattr(row, median_col)
            left = max(q1, xlim[0])
            right = min(q3, xlim[1])
            if left <= right:
                ax.hlines(yi, left, right, color=color, lw=1.7, zorder=3)
            if xlim[0] <= q1 <= xlim[1]:
                ax.vlines(q1, yi - 0.065, yi + 0.065, color=color, lw=1.7, zorder=3)
            if xlim[0] <= q3 <= xlim[1]:
                ax.vlines(q3, yi - 0.065, yi + 0.065, color=color, lw=1.7, zorder=3)
            if xlim[0] <= med <= xlim[1]:
                ax.scatter(
                    med,
                    yi,
                    s=38 if is_lead else 32,
                    color=color,
                    edgecolor="white",
                    linewidth=0.55,
                    zorder=4,
                )
                text_x = min(med + value_dx, xlim[1] - (xlim[1] - xlim[0]) * 0.02)
                ax.text(
                    text_x,
                    yi + 0.22,
                    value_fmt.format(med),
                    ha="left",
                    va="bottom",
                    fontsize=7.8,
                    color=TEXT,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.35},
                    zorder=5,
                )

        ax.set_xlim(*xlim)
        ax.set_xticks(ticks)
        ax.grid(axis="x", color="#DDE2E7", lw=0.75)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)

    axes[1].spines["left"].set_visible(False)
    axes[1].tick_params(axis="y", left=False, labelleft=False)
    axes[0].set_xlabel(xlabel)
    axes[1].set_xlabel("")

    d = 0.012
    kwargs = dict(transform=axes[0].transAxes, color="#4F545A", clip_on=False, lw=0.8)
    axes[0].plot((1 - d, 1 + d), (-d, +d), **kwargs)
    axes[0].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    kwargs = dict(transform=axes[1].transAxes, color="#4F545A", clip_on=False, lw=0.8)
    axes[1].plot((-d, +d), (-d, +d), **kwargs)
    axes[1].plot((-d, +d), (1 - d, 1 + d), **kwargs)


def main() -> None:
    setup_style()
    df = load_data()
    y = np.arange(len(df))[::-1]

    fig = plt.figure(figsize=(6.15, 2.08))
    gs = fig.add_gridspec(
        1,
        5,
        width_ratios=[1.02, 0.88, 1.12, 0.75, 1.00],
        wspace=0.07,
        left=0.055,
        right=0.985,
        bottom=0.29,
        top=0.95,
    )
    ax_label = fig.add_subplot(gs[0, 0])
    ax_expr_low = fig.add_subplot(gs[0, 1], sharey=ax_label)
    ax_expr_high = fig.add_subplot(gs[0, 2], sharey=ax_label)
    ax_frac_low = fig.add_subplot(gs[0, 3], sharey=ax_label)
    ax_frac_high = fig.add_subplot(gs[0, 4], sharey=ax_label)

    for ax in [ax_label, ax_expr_low, ax_expr_high, ax_frac_low, ax_frac_high]:
        for yi, pair in zip(y, df["pair"].astype(str)):
            fill = PALE_ACCENT if pair == "CD8 ET - SLC44A2" else "#FAFAFA"
            ax.axhspan(yi - 0.42, yi + 0.42, color=fill, zorder=0)
        ax.set_ylim(-0.55, len(df) - 0.45)

    ax_label.set_xlim(0, 1)
    ax_label.set_xticks([])
    ax_label.set_yticks(y)
    ax_label.set_yticklabels([])
    for yi, row in zip(y, df.itertuples(index=False)):
        pair = str(row.pair)
        ax_label.text(0.00, yi + 0.08, PAIR_LABELS[pair], ha="left", va="center", fontsize=9.4, color=TEXT)
        ax_label.text(
            0.00,
            yi - 0.22,
            f"n={int(row.n_samples)} donors",
            ha="left",
            va="center",
            fontsize=7.3,
            color="#60666D",
        )
    for spine in ax_label.spines.values():
        spine.set_visible(False)
    ax_label.tick_params(left=False, bottom=False)

    draw_interval_axis(
        (ax_expr_low, ax_expr_high),
        df,
        y,
        median_col="median_mean_expr",
        q1_col="q1_mean_expr",
        q3_col="q3_mean_expr",
        xlims=((0.0, 0.095), (0.30, 0.52)),
        xticks=([0, 0.04, 0.08], [0.35, 0.45]),
        xlabel="Median pseudobulk\nexpression",
        value_fmt="{:.3f}",
        value_dx=0.006,
    )
    draw_interval_axis(
        (ax_frac_low, ax_frac_high),
        df,
        y,
        median_col="median_pct_expr",
        q1_col="q1_pct_expr",
        q3_col="q3_pct_expr",
        xlims=((0.0, 7.0), (18.0, 36.0)),
        xticks=([0, 3, 6], [20, 30]),
        xlabel="Median positive\ncells (%)",
        value_fmt="{:.1f}",
        value_dx=0.45,
    )
    for ax in [ax_expr_low, ax_expr_high, ax_frac_low, ax_frac_high]:
        ax.tick_params(axis="y", labelleft=False, left=False)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{OUT}.{ext}", dpi=450 if ext == "png" else 300)
    plt.close(fig)
    print(f"Wrote {OUT}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
