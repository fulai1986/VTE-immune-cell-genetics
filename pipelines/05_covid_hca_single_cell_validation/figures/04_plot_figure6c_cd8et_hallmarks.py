#!/usr/bin/env python3
"""Figure 6C: Hallmark programs by CD8 ET SLC44A2 high/low groups."""

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
DATA = RESULTS / "cd8et_slc44a2_pseudobulk_mean_hallmark_high_low.tsv"
OUT = FIGURES / "figure6c_cd8et_slc44a2_hallmark_high_low_heatmap"

SELECTED = [
    "TNFA_SIGNALING_VIA_NFKB",
    "INFLAMMATORY_RESPONSE",
    "IL6_JAK_STAT3_SIGNALING",
    "COMPLEMENT",
    "COAGULATION",
    "HYPOXIA",
    "REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "APOPTOSIS",
    "P53_PATHWAY",
    "TGF_BETA_SIGNALING",
    "PI3K_AKT_MTOR_SIGNALING",
    "HEME_METABOLISM",
]

CATEGORY = {
    "HEME_METABOLISM": "Metabolic\nremodeling",
    "PI3K_AKT_MTOR_SIGNALING": "Growth-factor\nand signaling",
    "TGF_BETA_SIGNALING": "Growth-factor\nand signaling",
    "P53_PATHWAY": "Stress /\napoptosis /\nDNA damage",
    "APOPTOSIS": "Stress /\napoptosis /\nDNA damage",
    "REACTIVE_OXYGEN_SPECIES_PATHWAY": "Stress /\napoptosis /\nDNA damage",
    "HYPOXIA": "Stress /\napoptosis /\nDNA damage",
    "COAGULATION": "Immune /\ninflammatory /\ncoagulation",
    "COMPLEMENT": "Immune /\ninflammatory /\ncoagulation",
    "IL6_JAK_STAT3_SIGNALING": "Immune /\ninflammatory /\ncoagulation",
    "INFLAMMATORY_RESPONSE": "Immune /\ninflammatory /\ncoagulation",
    "TNFA_SIGNALING_VIA_NFKB": "Immune /\ninflammatory /\ncoagulation",
}

CATEGORY_COLOR = {
    "Metabolic\nremodeling": "#4F6D7A",
    "Growth-factor\nand signaling": "#8172B3",
    "Stress /\napoptosis /\nDNA damage": "#D9822B",
    "Immune /\ninflammatory /\ncoagulation": "#C44E52",
}


def clean_name(name: str) -> str:
    replacements = {
        "TNFA_SIGNALING_VIA_NFKB": "TNF-NFkB signaling",
        "INFLAMMATORY_RESPONSE": "Inflammatory response",
        "IL6_JAK_STAT3_SIGNALING": "IL6-JAK-STAT3 signaling",
        "COMPLEMENT": "Complement",
        "COAGULATION": "Coagulation",
        "HYPOXIA": "Hypoxia",
        "REACTIVE_OXYGEN_SPECIES_PATHWAY": "Reactive oxygen species",
        "APOPTOSIS": "Apoptosis",
        "P53_PATHWAY": "P53 pathway",
        "TGF_BETA_SIGNALING": "TGF-beta signaling",
        "PI3K_AKT_MTOR_SIGNALING": "PI3K-AKT-mTOR signaling",
        "HEME_METABOLISM": "Heme metabolism",
    }
    return replacements.get(name, name.replace("_", " ").title())


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def main() -> None:
    setup_style()
    df = pd.read_csv(DATA, sep="\t")
    df = df[df["hallmark"].isin(SELECTED)].copy()
    df["hallmark"] = pd.Categorical(df["hallmark"], categories=SELECTED[::-1], ordered=True)
    df = df.sort_values("hallmark").reset_index(drop=True)
    df["label"] = df["hallmark"].astype(str).map(clean_name)
    df["category"] = df["hallmark"].astype(str).map(CATEGORY)
    df["fdr_label"] = df["FDR"].map(lambda x: f"{x:.1e}" if x < 0.001 else f"{x:.3f}")

    mat = df[["median_low", "median_high"]].to_numpy()
    delta = df["delta_high_minus_low"].to_numpy()
    y = np.arange(len(df))

    fig = plt.figure(figsize=(8.4, 4.55))
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[2.35, 1.78, 1.35, 0.92],
        wspace=0.12,
        left=0.055,
        right=0.975,
        top=0.84,
        bottom=0.20,
    )
    ax_left = fig.add_subplot(gs[0, 0])
    ax_hm = fig.add_subplot(gs[0, 1])
    ax_delta = fig.add_subplot(gs[0, 2])
    ax_fdr = fig.add_subplot(gs[0, 3])

    for ax in [ax_left, ax_hm, ax_delta, ax_fdr]:
        ax.set_ylim(len(df) - 0.5, -0.5)

    for i in range(1, len(df)):
        ax_left.hlines(i - 0.5, xmin=0.38, xmax=0.98, color="#CFCFCF", lw=0.6, linestyles=(0, (3, 2)), zorder=0)

    vmax = float(np.nanmax(np.abs(mat)))
    im = ax_hm.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax_hm.set_xticks([0, 1])
    n_low = int(df["n_low"].iloc[0])
    n_high = int(df["n_high"].iloc[0])
    ax_hm.set_xticklabels([f"Low (n={n_low})", f"High (n={n_high})"])
    ax_hm.set_yticks([])
    ax_hm.set_yticklabels([])
    ax_hm.tick_params(axis="both", length=0, left=False, right=False)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            text_color = "white" if abs(mat[i, j]) >= 0.23 else "#20242A"
            ax_hm.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8.3, color=text_color)
    ax_hm.set_title("Median Hallmark score", fontsize=10.5, fontweight="bold", pad=18)
    ax_hm.plot([-0.48, 1.48], [-0.90, -0.90], color="#20242A", lw=0.8, clip_on=False)
    ax_hm.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax_hm.set_yticks(np.arange(-0.5, len(df), 1), minor=True)
    ax_hm.grid(which="minor", color="white", linewidth=1.4)
    ax_hm.tick_params(which="minor", bottom=False, left=False)

    cax = fig.add_axes([0.395, 0.078, 0.21, 0.020])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("Hallmark score", fontsize=7.5, labelpad=1)
    cbar.ax.tick_params(labelsize=7, length=2, pad=1)

    ax_delta.barh(y, delta, color="#9B1B30", height=0.58)
    ax_delta.axvline(0, color="#4F545A", lw=0.8)
    ax_delta.set_xlim(0, 0.82)
    ax_delta.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax_delta.set_title("High-Low (Delta)", fontsize=10.0, fontweight="bold", pad=18)
    ax_delta.set_xlabel("Delta (High-Low)", fontsize=8)
    ax_delta.set_yticks([])
    ax_delta.tick_params(axis="y", left=False, right=False, labelleft=False, length=0)
    ax_delta.tick_params(axis="x", labelsize=7)
    ax_delta.grid(axis="x", color="#E2E5E9", lw=0.7)
    ax_delta.spines["top"].set_visible(False)
    ax_delta.spines["right"].set_visible(False)
    ax_delta.spines["left"].set_visible(False)
    for yi, dval in zip(y, delta):
        ax_delta.text(dval + 0.035, yi, f"{dval:.2f}", ha="left", va="center", fontsize=7.8, color="#20242A")

    ax_fdr.set_xlim(0, 1)
    ax_fdr.set_xticks([])
    ax_fdr.set_yticks([])
    ax_fdr.tick_params(axis="y", left=False, right=False, labelleft=False, length=0)
    ax_fdr.set_title("FDR\nBenjamini-Hochberg\nq-value", fontsize=8.8, fontweight="bold", pad=12, linespacing=1.25)
    for yi, fdr in zip(y, df["fdr_label"]):
        ax_fdr.text(0.50, yi, fdr, ha="center", va="center", fontsize=8.0, color="#20242A")
    for spine in ax_fdr.spines.values():
        spine.set_visible(False)

    ax_left.set_xlim(0, 1)
    ax_left.set_xticks([])
    ax_left.set_yticks([])
    ax_left.tick_params(axis="y", left=False, right=False, labelleft=False, length=0)
    ax_left.tick_params(axis="x", bottom=False, labelbottom=False)
    group_runs = []
    current_cat = None
    start = None
    for idx, cat in enumerate(df["category"]):
        if cat != current_cat:
            if current_cat is not None:
                group_runs.append((current_cat, start, idx - 1))
            current_cat = cat
            start = idx
    group_runs.append((current_cat, start, len(df) - 1))
    for cat, start_idx, end_idx in group_runs:
        y0 = start_idx - 0.39
        height = end_idx - start_idx + 0.78
        color = CATEGORY_COLOR[cat]
        ax_left.add_patch(plt.Rectangle((0.315, y0), 0.018, height, color=color, lw=0))
        ax_left.text(
            0.24,
            y0 + height / 2,
            cat,
            ha="right",
            va="center",
            fontsize=7.2,
            color=color,
            fontweight="bold" if height > 1.4 else "normal",
            linespacing=0.90,
            clip_on=False,
        )
    ax_left.text(0.47, -0.82, "Pathway", ha="left", va="bottom", fontsize=9.5, fontweight="bold", clip_on=False)
    for yi, label in enumerate(df["label"]):
        ax_left.text(0.38, yi, label, ha="left", va="center", fontsize=8.6, color="#20242A")
    for spine in ax_left.spines.values():
        spine.set_visible(False)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{OUT}.{ext}", dpi=450 if ext == "png" else 300)
    plt.close(fig)
    print(f"Wrote {OUT}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
