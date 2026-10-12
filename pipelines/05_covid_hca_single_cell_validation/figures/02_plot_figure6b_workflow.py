#!/usr/bin/env python3
"""Export standalone workflow panels for Result 5."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle, Rectangle


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
BASE = PROJECT / "covid_validation/stephenson_hca/result5_expression_context"
RESULTS = BASE / "results"
FIGURES = BASE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

PAIR_SCORES = RESULTS / "tier1_3_pairs_covid_d0_same_cell_program_scores.curated_mapping.tsv"

LOW = "#7A838D"
HIGH = "#B5163A"
MUTED = "#56616D"
GRID = "#E3E8EE"
LOW_FILL = "#ECEFF2"
HIGH_FILL = "#FFE3E5"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def assign_groups() -> pd.DataFrame:
    pair = pd.read_csv(PAIR_SCORES, sep="\t")
    metrics = pair[pair["cell"].eq("CD8_ET") & pair["Gene"].eq("SLC44A2")][
        ["sample_id", "candidate_mean_expr"]
    ].drop_duplicates()
    q_low, q_high = metrics["candidate_mean_expr"].quantile([1 / 3, 2 / 3])
    metrics["slc44a2_group"] = "Middle"
    metrics.loc[metrics["candidate_mean_expr"] <= q_low, "slc44a2_group"] = "Low"
    metrics.loc[metrics["candidate_mean_expr"] >= q_high, "slc44a2_group"] = "High"
    return metrics.sort_values("candidate_mean_expr").reset_index(drop=True)


def add_panel_box(ax, title: str, header_color: str = "#E9EDF2") -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(
        FancyBboxPatch(
            (0.02, 0.02),
            0.96,
            0.92,
            boxstyle="round,pad=0.012,rounding_size=0.015",
            facecolor="white",
            edgecolor="#D9DEE5",
            linewidth=0.8,
        )
    )
    ax.add_patch(Rectangle((0.02, 0.82), 0.96, 0.12, transform=ax.transAxes, color=header_color, ec="none"))
    ax.text(0.50, 0.88, title, ha="center", va="center", fontsize=10.5, fontweight="bold")


def draw_expression_panel(ax, groups: pd.DataFrame) -> None:
    add_panel_box(ax, "SLC44A2 expression in CD8 ET cells", "#E9EDF2")
    inner = ax.inset_axes([0.16, 0.32, 0.76, 0.46])
    x = list(range(len(groups)))
    colors = groups["slc44a2_group"].map({"Low": LOW, "Middle": "#B7BDC5", "High": HIGH}).fillna("#B7BDC5")
    inner.plot(x, groups["candidate_mean_expr"], color="#AEB6BF", lw=0.85, zorder=1)
    inner.scatter(x, groups["candidate_mean_expr"], s=16, c=colors, edgecolor="white", linewidth=0.25, zorder=2)
    q1 = groups.index[groups["slc44a2_group"].eq("Low")].max()
    q2 = groups.index[groups["slc44a2_group"].eq("High")].min()
    inner.axvline(q1 + 0.5, color="#7E8792", lw=0.8, ls="--")
    inner.axvline(q2 - 0.5, color="#7E8792", lw=0.8, ls="--")
    inner.set_xlabel("Donors ranked by expression", fontsize=8.4, labelpad=6)
    inner.set_ylabel("Pseudobulk expression", fontsize=8.4, labelpad=7)
    inner.tick_params(labelsize=7.2, length=2)
    inner.spines["top"].set_visible(False)
    inner.spines["right"].set_visible(False)
    inner.grid(axis="y", color=GRID, lw=0.7)
    inner.set_xlim(-2, len(groups) + 1)

    ax.add_patch(Rectangle((0.205, 0.160), 0.255, 0.035, color=LOW_FILL, ec="none"))
    ax.add_patch(Rectangle((0.585, 0.160), 0.255, 0.035, color=HIGH_FILL, ec="none"))
    n_low = int(groups["slc44a2_group"].eq("Low").sum())
    n_high = int(groups["slc44a2_group"].eq("High").sum())
    ax.text(0.333, 0.130, "Low", ha="center", va="top", fontsize=9.0, fontweight="bold")
    ax.text(0.333, 0.070, f"n = {n_low}", ha="center", va="top", fontsize=8.8)
    ax.text(0.713, 0.130, "High", ha="center", va="top", fontsize=9.0, fontweight="bold")
    ax.text(0.713, 0.070, f"n = {n_high}", ha="center", va="top", fontsize=8.8)


def draw_workflow_network(ax) -> None:
    add_panel_box(ax, "Cell-cell communication\n(LIANA, donor-level)", "#E9EDF2")
    pos = {"CD8 ET": (0.50, 0.68), "Monocyte": (0.25, 0.39), "Platelet": (0.75, 0.39)}
    colors = {"CD8 ET": "#8BC6B2", "Monocyte": "#9E9AD1", "Platelet": "#DF5961"}
    edges = [
        ("CD8 ET", "Monocyte", 1.7, "#8A939D", -0.12),
        ("CD8 ET", "Platelet", 3.7, HIGH, 0.10),
        ("Monocyte", "Platelet", 2.9, HIGH, 0.00),
    ]
    for s, t, lw, col, rad in edges:
        ax.add_patch(
            FancyArrowPatch(
                pos[s],
                pos[t],
                arrowstyle="-|>",
                mutation_scale=14,
                linewidth=lw,
                color=col,
                alpha=0.82,
                connectionstyle=f"arc3,rad={rad}",
                shrinkA=24,
                shrinkB=24,
            )
        )
    for name, (x, y) in pos.items():
        ax.add_patch(Circle((x, y), 0.066, color=colors[name], ec="#425466", lw=0.9, zorder=3))
    ax.text(0.50, 0.79, "CD8 ET", ha="center", va="center", fontsize=9.2, fontweight="bold")
    ax.text(0.25, 0.245, "Monocyte", ha="center", va="center", fontsize=8.9)
    ax.text(0.75, 0.245, "Platelet", ha="center", va="center", fontsize=8.9)

    ax.text(0.08, 0.135, "Predefined axes", fontsize=8.3, fontweight="bold", ha="left", va="center")
    ax.text(0.08, 0.078, "CD8 ET \u2192 monocyte", fontsize=7.35, ha="left", va="center", color="#56616D")
    ax.text(0.405, 0.078, "CD8 ET \u2192 platelet", fontsize=7.35, ha="left", va="center", color=HIGH)
    ax.text(0.72, 0.078, "monocyte \u2192 platelet", fontsize=7.35, ha="left", va="center", color=HIGH)


def save_single(draw_func, out_prefix: Path, groups=None, figsize=(3.0, 2.35)) -> None:
    fig, ax = plt.subplots(figsize=figsize)
    if groups is None:
        draw_func(ax)
    else:
        draw_func(ax, groups)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{out_prefix}.{ext}", dpi=450, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    setup_style()
    groups = assign_groups()
    save_single(draw_expression_panel, FIGURES / "figure6b_slc44a2_expression_grouping", groups, figsize=(3.65, 2.75))
    save_single(draw_workflow_network, FIGURES / "figure6b_liana_axis_schematic", figsize=(3.45, 2.35))


if __name__ == "__main__":
    main()
