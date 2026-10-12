#!/usr/bin/env python3
"""Compact Figure-6C-style heatmap for continuous CD8 ET associations."""

from __future__ import annotations

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
DATA = HERE / "continuous_severity_results.tsv"
SELECTED = [
    "HEME_METABOLISM", "PI3K_AKT_MTOR_SIGNALING", "TGF_BETA_SIGNALING",
    "P53_PATHWAY", "APOPTOSIS", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "HYPOXIA",
    "COAGULATION", "COMPLEMENT", "IL6_JAK_STAT3_SIGNALING",
    "INFLAMMATORY_RESPONSE", "TNFA_SIGNALING_VIA_NFKB",
]
LABELS = ["Heme metabolism", "PI3K–AKT–mTOR signaling", "TGF-β signaling",
          "P53 pathway", "Apoptosis", "Reactive oxygen species", "Hypoxia",
          "Coagulation", "Complement", "IL6–JAK–STAT3 signaling",
          "Inflammatory response", "TNF–NFκB signaling"]
GROUPS = ["Metabolic", "Growth/signaling", "Growth/signaling", "Stress", "Stress", "Stress", "Stress",
          "Immune/coagulation", "Immune/coagulation", "Immune/coagulation", "Immune/coagulation", "Immune/coagulation"]
GROUP_COLORS = {"Metabolic": "#4F6D7A", "Growth/signaling": "#8172B3",
                "Stress": "#D9822B", "Immune/coagulation": "#C44E52"}
GROUP_LABELS = [
    ("Metabolic\nremodeling", 0, 0, "Metabolic"),
    ("Growth-factor\nand signaling", 1, 2, "Growth/signaling"),
    ("Stress / apoptosis /\nDNA damage", 3, 6, "Stress"),
    ("Immune / inflammatory /\ncoagulation", 7, 11, "Immune/coagulation"),
]


def fdr_label(q: float) -> str:
    return f"{q:.1e}" if q < .001 else f"{q:.3f}"


def main() -> None:
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    data = pd.read_csv(DATA, sep="\t").set_index("readout").loc[SELECTED].reset_index()
    assert data.n_donors.eq(84).all()
    values = data[["unadjusted_rho", "severity_adjusted_rho"]].to_numpy(float)
    delta = data.severity_adjusted_rho.to_numpy(float) - data.unadjusted_rho.to_numpy(float)
    n = len(data)

    fig = plt.figure(figsize=(10.4, 5.25))
    gs = fig.add_gridspec(1, 4, width_ratios=[.13, 2.25, 1.9, .9],
                          left=.39, right=.96, bottom=.15, top=.90, wspace=.13)
    ax_group = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1], sharey=ax_group)
    ax_delta = fig.add_subplot(gs[0, 2], sharey=ax_group)
    ax_fdr = fig.add_subplot(gs[0, 3], sharey=ax_group)

    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im = ax_heat.imshow(values, aspect="auto", cmap="RdBu_r", norm=norm, interpolation="none")
    for i in range(n):
        ax_group.add_patch(plt.Rectangle((0, i-.5), 1, 1, color=GROUP_COLORS[GROUPS[i]], lw=0))
        for j in range(2):
            ax_heat.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center",
                         color="white" if values[i, j] >= .6 else "#17202A", fontsize=8.3,
                         fontweight="bold" if j == 1 else "normal")
        color = "#A71935" if delta[i] >= 0 else "#377DA8"
        ax_delta.barh(i, delta[i], color=color, height=.55)
        ax_delta.text(delta[i] + (.012 if delta[i] >= 0 else -.012), i, f"{delta[i]:+.3f}",
                      va="center", ha="left" if delta[i] >= 0 else "right", fontsize=8,
                      color="#26323B")
        ax_fdr.text(.04, i, fdr_label(float(data.severity_adjusted_fdr.iloc[i])),
                    va="center", ha="left", color="#A71935" if data.severity_adjusted_fdr.iloc[i] < .05 else "#677381",
                    fontsize=8.3)

    for ax in [ax_group, ax_heat, ax_delta, ax_fdr]:
        ax.set_ylim(n-.5, -.5)
        ax.set_yticks(np.arange(n))
        ax.tick_params(axis="y", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    ax_group.set_xlim(0, 1)
    ax_group.set_xticks([])
    ax_group.tick_params(axis="y", labelleft=False)
    ax_heat.tick_params(axis="y", labelleft=False)
    for i, label in enumerate(LABELS):
        ypos = .90 - (i + .5) * (.90 - .15) / n
        fig.text(.375, ypos, label, ha="right", va="center", fontsize=8.3, color="#17202A")
    for label, first, last, group in GROUP_LABELS:
        midpoint = (first + last) / 2
        ypos = .90 - (midpoint + .5) * (.90 - .15) / n
        fig.text(.22, ypos, label, ha="right", va="center", fontsize=8.1,
                 color=GROUP_COLORS[group], fontweight="bold")
    ax_heat.set_xticks([0, 1], ["Unadjusted", "D0 severity\nadjusted"], fontsize=8.4)
    ax_heat.xaxis.tick_top()
    ax_heat.tick_params(axis="x", length=0, pad=9)
    for i in range(n+1):
        ax_heat.axhline(i-.5, color="white", linewidth=1.5)
    ax_heat.axvline(.5, color="white", linewidth=1.5)
    ax_delta.set_xlim(-.17, .17)
    ax_delta.tick_params(axis="y", labelleft=False)
    ax_delta.axvline(0, color="#9BA4AB", lw=.8)
    ax_delta.set_xticks([-.1, 0, .1], ["−0.1", "0", "+0.1"], fontsize=7.5)
    ax_delta.set_title("Change in ρ\n(adjusted − raw)", fontsize=8.3, pad=10)
    ax_fdr.set_xlim(0, 1)
    ax_fdr.tick_params(axis="y", labelleft=False)
    ax_fdr.set_xticks([])
    ax_fdr.set_title("Adjusted\nFDR", fontsize=8.3, pad=10, loc="left")
    for i in range(1, n):
        if GROUPS[i] != GROUPS[i-1]:
            for ax in [ax_group, ax_heat, ax_delta, ax_fdr]:
                ax.axhline(i-.5, color="#B9C1C8", linewidth=.8)

    cax = fig.add_axes([.43, .072, .26, .016])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-1, -.5, 0, .5, 1])
    cb.set_label("Spearman ρ", fontsize=7.5, labelpad=2)
    cb.ax.tick_params(labelsize=7, length=2)
    fig.savefig(HERE / "Figure_S3_cd8et_hallmark_severity_sensitivity.png", dpi=300)
    fig.savefig(HERE / "Figure_S3_cd8et_hallmark_severity_sensitivity.pdf")
    plt.close(fig)
    print(HERE / "Figure_S3_cd8et_hallmark_severity_sensitivity.png")


if __name__ == "__main__":
    main()
