#!/usr/bin/env python3
"""Platelet Hallmark raw/adjusted comparison in the established CD8 heatmap style."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
HERE = PROJECT / "manuscript_results/current_sensitivity"
HERE.mkdir(parents=True, exist_ok=True)
SOURCE = HERE / "continuous_severity_results.tsv"
ORDER = [
    ("HALLMARK_COAGULATION", "Coagulation", "Immune"),
    ("HALLMARK_COMPLEMENT", "Complement", "Immune"),
    ("HALLMARK_IL6_JAK_STAT3_SIGNALING", "IL6–JAK–STAT3 signaling", "Immune"),
    ("HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response", "Immune"),
    ("HALLMARK_INTERFERON_ALPHA_RESPONSE", "Interferon alpha response", "Immune"),
    ("HALLMARK_INTERFERON_GAMMA_RESPONSE", "Interferon gamma response", "Immune"),
    ("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF–NFκB signaling", "Immune"),
    ("HALLMARK_APOPTOSIS", "Apoptosis", "Stress"),
    ("HALLMARK_HYPOXIA", "Hypoxia", "Stress"),
]
COLORS = {"Immune": "#C44E52", "Stress": "#D9822B"}


def main() -> None:
    HERE.mkdir(exist_ok=True)
    source = pd.read_csv(SOURCE, sep="\t")
    source = source.loc[source.family.eq("Platelet Hallmark")].set_index("readout")
    assert len(source) == len(ORDER) == 9
    data = source.loc[[row[0] for row in ORDER]].reset_index()
    data.insert(1, "pathway_label", [row[1] for row in ORDER])
    data.insert(2, "category", [row[2] for row in ORDER])
    data["change_in_rho"] = data.severity_adjusted_rho - data.unadjusted_rho
    assert data.n_donors.eq(76).all()
    data[["readout", "pathway_label", "category", "n_donors", "unadjusted_rho",
          "severity_adjusted_rho", "change_in_rho", "severity_adjusted_fdr"]].to_csv(
              HERE / "platelet_comparison_heatmap_display_data.tsv", sep="\t", index=False)

    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    n = len(data)
    top, bottom = .85, .18
    fig = plt.figure(figsize=(10.4, 4.45), facecolor="white")
    gs = fig.add_gridspec(1, 4, width_ratios=[.13, 2.25, 1.9, .9],
                          left=.39, right=.96, bottom=bottom, top=top, wspace=.13)
    group = fig.add_subplot(gs[0, 0])
    heat = fig.add_subplot(gs[0, 1], sharey=group)
    delta_ax = fig.add_subplot(gs[0, 2], sharey=group)
    fdr_ax = fig.add_subplot(gs[0, 3], sharey=group)
    values = data[["unadjusted_rho", "severity_adjusted_rho"]].to_numpy(float)
    im = heat.imshow(values, aspect="auto", cmap="RdBu_r",
                     norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), interpolation="none")

    for i, row in enumerate(data.itertuples(index=False)):
        group.add_patch(plt.Rectangle((0, i-.5), 1, 1, color=COLORS[row.category], lw=0))
        for j in range(2):
            v = values[i, j]
            heat.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8.5,
                      color="white" if v >= .6 else "#17202A", weight="bold" if j == 1 else "normal")
        delta = float(row.change_in_rho)
        color = "#A71935" if delta >= 0 else "#377DA8"
        delta_ax.barh(i, delta, color=color, height=.56)
        delta_ax.text(delta + (.006 if delta >= 0 else -.006), i, f"{delta:+.3f}",
                      ha="left" if delta >= 0 else "right", va="center", fontsize=8, color="#26323B")
        fdr_ax.text(.04, i, f"{row.severity_adjusted_fdr:.2e}", ha="left", va="center",
                    fontsize=8.2, color="#A71935")
        y = top - (i+.5)*(top-bottom)/n
        fig.text(.375, y, row.pathway_label, ha="right", va="center", fontsize=8.5, color="#17202A")

    for first, last, label, category in [
        (0, 6, "Immune / inflammatory /\ncoagulation", "Immune"),
        (7, 8, "Stress / apoptosis", "Stress")]:
        y = top - ((first+last)/2+.5)*(top-bottom)/n
        fig.text(.22, y, label, ha="right", va="center", fontsize=8.2,
                 color=COLORS[category], weight="bold")
    for ax in [group, heat, delta_ax, fdr_ax]:
        ax.set_ylim(n-.5, -.5)
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.axhline(6.5, color="#B9C1C8", lw=.8)
    group.set_xlim(0, 1); group.set_xticks([])
    heat.set_xticks([0, 1], ["Unadjusted", "D0 severity\nadjusted"], fontsize=8.7)
    heat.xaxis.tick_top(); heat.tick_params(axis="x", length=0, pad=9)
    for i in range(n+1):
        heat.axhline(i-.5, color="white", lw=1.5)
    heat.axvline(.5, color="white", lw=1.5)
    delta_ax.set_xlim(-.1, .1)
    delta_ax.axvline(0, color="#9BA4AB", lw=.8)
    delta_ax.set_xticks([-.1, 0, .1], ["−0.1", "0", "+0.1"], fontsize=7.8)
    delta_ax.set_title("Change in ρ\n(adjusted − raw)", fontsize=8.7, pad=10)
    fdr_ax.set_xlim(0, 1); fdr_ax.set_xticks([])
    fdr_ax.set_title("Adjusted\nFDR", fontsize=8.7, pad=10, loc="left")
    cax = fig.add_axes([.43, .085, .26, .018])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-1, -.5, 0, .5, 1])
    cb.set_label("Spearman ρ", fontsize=8, labelpad=2)
    cb.ax.tick_params(labelsize=7.5, length=2)
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"Figure_S5_platelet_hallmark_severity_sensitivity.{ext}", dpi=300)
    plt.close(fig)
    print("Platelet comparison heatmap saved")


if __name__ == "__main__":
    main()
