#!/usr/bin/env python3
"""Plot Supplementary Figure S2A-B for HCA to OneK1K-compatible mapping."""

import os
from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
import pandas as pd


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
MAPPING_DIR = PROJECT / "covid_validation/hca_onek1k_mapping"
RESULTS = MAPPING_DIR / "results"
OUT = PROJECT / "manuscript_results/completed_results/result5_covid_single_cell_validation/supplementary"
OUT.mkdir(parents=True, exist_ok=True)

CONSENSUS = RESULTS / "hca_to_onek1k_mapping_main_with_seurat_sensitivity.tsv"
COVID_COUNTS = RESULTS / "hca_covid_cluster_onek1k_compatible_counts.tsv"

CELL_ORDER = [
    "CD4 NC",
    "CD4 ET",
    "CD4 SOX4",
    "CD8 ET",
    "CD8 NC",
    "CD8 S100B",
    "NK",
    "NK R",
    "B IN",
    "B Mem",
    "Plasma",
    "Mono C",
    "Mono NC",
    "DC",
    "HCA-specific",
    "Not mapped",
]

ONEK1K_CELL_ORDER = [
    "CD4 NC",
    "CD4 ET",
    "CD4 SOX4",
    "CD8 ET",
    "CD8 NC",
    "CD8 S100B",
    "NK",
    "NK R",
    "B IN",
    "B Mem",
    "Plasma",
    "Mono C",
    "Mono NC",
    "DC",
]

HCA_CLUSTER_ORDER = [
    "CD4.Naive", "CD4.CM", "CD4.Tfh", "CD4.IL22", "Treg", "CD4.EM",
    "CD4.Th1", "CD4.Prolif", "CD4.Th2", "CD4.Th17", "CD8.TE", "gdT",
    "CD8.EM", "CD8.Prolif", "NKT", "CD8.Naive", "MAIT", "NK_16hi",
    "NK_prolif", "NK_56hi", "ILC1_3", "ILC2", "B_naive", "B_immature",
    "B_switched_memory", "B_non-switched_memory", "B_exhausted", "Plasmablast",
    "Plasma_cell_IgG", "Plasma_cell_IgA", "Plasma_cell_IgM", "CD14_mono",
    "CD83_CD14_mono", "Mono_prolif", "CD16_mono", "C1_CD16_mono", "DC1",
    "DC2", "DC3", "pDC", "ASDC", "DC_prolif",
]

CONF_COLORS = {
    "High": "#A0173D",
    "Intermediate": "#E7903F",
    "HCA-specific": "#5E6F82",
    "Not applicable": "#5E6F82",
    "Unresolved": "#8C8C8C",
}

LABEL_COLORS = {
    "CD4 NC": "#8A2D78",
    "CD4 ET": "#B260A5",
    "CD4 SOX4": "#C7A6D8",
    "CD8 ET": "#2F75B5",
    "CD8 NC": "#4F86C6",
    "CD8 S100B": "#8BB9E3",
    "NK": "#3BAA55",
    "NK R": "#86B96F",
    "B IN": "#F2AF32",
    "B Mem": "#E7C437",
    "Plasma": "#B89CCB",
    "Mono C": "#E77728",
    "Mono NC": "#E24A2A",
    "DC": "#C51B2D",
    "HCA-specific": "#6B7280",
    "Not mapped": "#9CA3AF",
}

LINEAGE_COLORS = {
    "T cell": "#6E58A8",
    "NK": "#2A7F62",
    "B cell": "#D9A514",
    "Myeloid": "#C95743",
    "Platelet/RBC": "#5E6F82",
    "Other": "#9E9E9E",
}


def set_style():
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
        }
    )


def rounded_box(ax, xy, w, h, text, fc, ec="#D5DAE1", color="#111827", fs=10, weight="bold"):
    patch = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.035",
        fc=fc,
        ec=ec,
        lw=1.2,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + w / 2,
        xy[1] + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fs,
        color=color,
        fontweight=weight,
        linespacing=1.25,
    )
    return patch


def arrow(ax, start, end, color="#6B7280", lw=1.5, rad=0.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            lw=lw,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
        )
    )


def panel_a():
    set_style()
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.02, 0.95, "A", fontsize=18, fontweight="bold", ha="left", va="top")
    ax.text(
        0.10,
        0.94,
        "HCA healthy-to-OneK1K-compatible cell-type mapping workflow",
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="top",
    )

    rounded_box(
        ax,
        (0.06, 0.56),
        0.21,
        0.22,
        "OneK1K PBMC\nscience14 reference\n14 eQTL cell types",
        fc="#F4F6F8",
        ec="#8FA3B8",
        fs=9.5,
    )
    rounded_box(
        ax,
        (0.06, 0.22),
        0.21,
        0.22,
        "HCA healthy cells\nStephenson et al.\ncluster annotations",
        fc="#F4F6F8",
        ec="#8FA3B8",
        fs=9.5,
    )
    rounded_box(
        ax,
        (0.39, 0.38),
        0.24,
        0.25,
        "Lineage-aware\nconsensus mapping\ncentroid + marker support",
        fc="#FFF6E6",
        ec="#D7A95C",
        fs=9.5,
    )
    rounded_box(
        ax,
        (0.74, 0.48),
        0.20,
        0.21,
        "Fixed mapping\napplied to\nHCA COVID cells",
        fc="#FBE8EE",
        ec="#C24764",
        fs=9.5,
    )
    rounded_box(
        ax,
        (0.74, 0.19),
        0.20,
        0.17,
        "Platelets retained as\nHCA-specific effector\ncompartment",
        fc="#EEF2F6",
        ec="#7B8794",
        fs=8.8,
        color="#263241",
    )

    arrow(ax, (0.28, 0.67), (0.39, 0.53))
    arrow(ax, (0.28, 0.33), (0.39, 0.48))
    arrow(ax, (0.63, 0.51), (0.74, 0.59))
    arrow(ax, (0.63, 0.43), (0.74, 0.29), color="#7B8794")

    ax.plot([0.70, 0.96], [0.79, 0.79], color="#C24764", lw=1.2)
    ax.text(
        0.83,
        0.84,
        "COVID cells were not used to define cell identities",
        ha="center",
        va="center",
        fontsize=8.5,
        color="#A0173D",
        fontweight="bold",
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"Figure_S2A_hca_onek1k_mapping_workflow.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def infer_lineage(label):
    label = str(label)
    if label.startswith("CD4") or label.startswith("CD8"):
        return "T cell"
    if label.startswith("NK"):
        return "NK"
    if label.startswith("B") or label == "Plasma":
        return "B cell"
    if label.startswith("Mono") or label == "DC":
        return "Myeloid"
    if label in ["HCA-specific", "Not mapped"]:
        return "Platelet/RBC"
    return "Other"


def clean_label(x):
    if pd.isna(x) or str(x).strip() == "":
        return "Not mapped"
    return str(x)


def size_from_count(x):
    if x <= 0:
        return 24
    return 28 + 38 * math.log10(x + 10)


def panel_b():
    set_style()
    main = pd.read_csv(CONSENSUS, sep="\t")
    covid = pd.read_csv(COVID_COUNTS, sep="\t")

    df = main.merge(
        covid[["full_clustering", "covid_cells"]],
        left_on="HCA cluster",
        right_on="full_clustering",
        how="left",
    )
    df["covid_cells"] = df["covid_cells"].fillna(0).astype(int)
    df["label"] = df["Main OneK1K-compatible label"].apply(clean_label)
    df.loc[df["Main confidence"].isin(["HCA-specific", "Not applicable"]), "label"] = "HCA-specific"
    df["confidence"] = df["Main confidence"].replace({"Not applicable": "HCA-specific"}).fillna("Unresolved")
    df["lineage"] = df["label"].apply(infer_lineage)
    df["label"] = pd.Categorical(df["label"], CELL_ORDER, ordered=True)
    df["lineage"] = pd.Categorical(df["lineage"], ["T cell", "NK", "B cell", "Myeloid", "Platelet/RBC", "Other"], ordered=True)
    df = df.sort_values(["lineage", "label", "Healthy cells", "HCA cluster"], ascending=[True, True, False, True])
    df = df.reset_index(drop=True)

    fig_h = max(8.2, 0.18 * len(df) + 2.0)
    fig, ax = plt.subplots(figsize=(8.4, fig_h))
    ax.set_facecolor("white")

    x_map = {ct: i for i, ct in enumerate(CELL_ORDER)}
    y = list(range(len(df)))

    for i, row in df.iterrows():
        if i % 2 == 0:
            ax.add_patch(Rectangle((-0.5, i - 0.5), len(CELL_ORDER), 1, fc="#FAFAFA", ec="none", zorder=0))

    for i, row in df.iterrows():
        ax.scatter(
            x_map[str(row["label"])],
            i,
            s=size_from_count(row["covid_cells"]),
            c=CONF_COLORS.get(row["confidence"], "#8C8C8C"),
            edgecolors="#2F3640",
            linewidths=0.45,
            alpha=0.95,
            zorder=3,
        )

    ax.set_xlim(-0.7, len(CELL_ORDER) - 0.3)
    ax.set_ylim(len(df) - 0.5, -0.5)
    ax.set_xticks(range(len(CELL_ORDER)))
    ax.set_xticklabels(CELL_ORDER, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(y)
    ax.set_yticklabels(df["HCA cluster"])
    ax.tick_params(axis="both", length=0)
    ax.set_xlabel("OneK1K-compatible label")
    ax.set_ylabel("HCA cluster")
    ax.grid(axis="x", color="#E5E7EB", lw=0.8)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.text(-0.09, 1.035, "B", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")
    ax.text(
        0.0,
        1.03,
        "HCA cluster correspondence to OneK1K-compatible cell types",
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="top",
    )

    x0 = len(CELL_ORDER) - 0.1
    legend_y = 1
    for conf in ["High", "Intermediate", "HCA-specific", "Unresolved"]:
        ax.scatter(x0 + 1.0, legend_y, s=90, c=CONF_COLORS[conf], edgecolors="#2F3640", linewidths=0.45, clip_on=False)
        ax.text(x0 + 1.35, legend_y, conf, va="center", fontsize=8, clip_on=False)
        legend_y += 1.15
    ax.text(x0 + 0.95, 0.0, "Mapping\nconfidence", fontsize=8.5, fontweight="bold", ha="left", va="bottom", clip_on=False)

    size_y = legend_y + 1.1
    ax.text(x0 + 0.95, size_y - 0.4, "COVID\ncells", fontsize=8.5, fontweight="bold", ha="left", va="bottom", clip_on=False)
    for cnt in [1000, 10000, 50000]:
        size_y += 1.15
        ax.scatter(x0 + 1.0, size_y, s=size_from_count(cnt), c="white", edgecolors="#2F3640", linewidths=0.6, clip_on=False)
        ax.text(x0 + 1.35, size_y, f"{cnt:,}", va="center", fontsize=8, clip_on=False)

    ax.text(
        0.0,
        -0.055,
        "Platelets are retained as an HCA-specific downstream effector compartment and are not forced into the OneK1K eQTL exposure framework.",
        transform=ax.transAxes,
        fontsize=7.5,
        color="#4B5563",
        ha="left",
        va="top",
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"Figure_S2B_hca_cluster_onek1k_correspondence.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def bezier_band(ax, x0, x1, y0_low, y0_high, y1_low, y1_high, color, alpha=0.72):
    c = (x1 - x0) * 0.45
    verts = [
        (x0, y0_low),
        (x0 + c, y0_low),
        (x1 - c, y1_low),
        (x1, y1_low),
        (x1, y1_high),
        (x1 - c, y1_high),
        (x0 + c, y0_high),
        (x0, y0_high),
        (x0, y0_low),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), fc=color, ec="none", alpha=alpha, zorder=1))


def panel_b_sankey():
    set_style()
    main = pd.read_csv(CONSENSUS, sep="\t")
    covid = pd.read_csv(COVID_COUNTS, sep="\t")
    df = main.merge(
        covid[["full_clustering", "covid_cells"]],
        left_on="HCA cluster",
        right_on="full_clustering",
        how="left",
    )
    df["covid_cells"] = df["covid_cells"].fillna(0).astype(int)
    df["label"] = df["Main OneK1K-compatible label"].apply(clean_label)
    df.loc[df["Main confidence"].isin(["HCA-specific", "Not applicable"]), "label"] = "HCA-specific"
    df = df[
        (~df["HCA cluster"].isin(["Platelets", "RBC"]))
        & (~df["label"].isin(["HCA-specific", "Not mapped"]))
    ].copy()
    df["confidence"] = df["Main confidence"].replace({"Not applicable": "HCA-specific"}).fillna("Unresolved")
    df["lineage"] = df["label"].apply(infer_lineage)
    df["lineage_order"] = df["lineage"].map(
        {"T cell": 0, "NK": 1, "B cell": 2, "Myeloid": 3, "Platelet/RBC": 4, "Other": 5}
    )
    df["label_order"] = df["label"].map({v: i for i, v in enumerate(ONEK1K_CELL_ORDER)})
    df = df.sort_values(["lineage_order", "label_order", "covid_cells", "HCA cluster"], ascending=[True, True, False, True])

    total = df["covid_cells"].sum()
    scale = 0.82 / total
    gap_left = 0.0025
    gap_right = 0.009

    y_top = 0.91
    left_pos = {}
    y = y_top
    for _, row in df.iterrows():
        h = max(row["covid_cells"] * scale, 0.0014)
        left_pos[row["HCA cluster"]] = (y - h, y)
        y -= h + gap_left

    right_pos = {}
    label_ranges = {}
    y = y_top
    for label in ONEK1K_CELL_ORDER:
        sub = df[df["label"] == label].sort_values(["lineage_order", "covid_cells"], ascending=[True, False])
        if sub.empty:
            continue
        block_top = y
        for _, row in sub.iterrows():
            h = max(row["covid_cells"] * scale, 0.0014)
            right_pos[row["HCA cluster"]] = (y - h, y)
            y -= h + gap_left
        block_low = y + gap_left
        label_ranges[label] = (block_low, block_top, int(sub["covid_cells"].sum()))
        y -= gap_right

    fig, ax = plt.subplots(figsize=(10.6, 8.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.015, 0.985, "B", fontsize=18, fontweight="bold", ha="left", va="top")
    ax.text(
        0.075,
        0.975,
        "HCA cluster correspondence to OneK1K-compatible cell types",
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="top",
    )

    x_left0, x_left1 = 0.245, 0.270
    x_right0, x_right1 = 0.755, 0.780
    band_x0, band_x1 = x_left1, x_right0

    for _, row in df.iterrows():
        l0, l1 = left_pos[row["HCA cluster"]]
        r0, r1 = right_pos[row["HCA cluster"]]
        color = LABEL_COLORS.get(row["label"], "#9CA3AF")
        bezier_band(ax, band_x0, band_x1, l0, l1, r0, r1, color=color, alpha=0.66)
        ax.add_patch(Rectangle((x_left0, l0), x_left1 - x_left0, l1 - l0, fc=color, ec="white", lw=0.25, zorder=2))
        ax.add_patch(Rectangle((x_right0, r0), x_right1 - x_right0, r1 - r0, fc=color, ec="white", lw=0.25, zorder=2))

    ax.text((x_left0 + x_left1) / 2, 0.94, "HCA cluster", ha="center", va="bottom", fontsize=9.5, fontweight="bold")
    ax.text((x_right0 + x_right1) / 2, 0.94, "OneK1K-compatible label", ha="center", va="bottom", fontsize=9.5, fontweight="bold")

    min_label_height = 0.012
    for _, row in df.iterrows():
        l0, l1 = left_pos[row["HCA cluster"]]
        mid = (l0 + l1) / 2
        if row["covid_cells"] >= 180 or (l1 - l0) >= min_label_height:
            ax.text(
                x_left0 - 0.012,
                mid,
                row["HCA cluster"],
                ha="right",
                va="center",
                fontsize=6.6,
                color="#1F2937",
            )

    for label, (low, high, count) in label_ranges.items():
        mid = (low + high) / 2
        ax.text(x_right1 + 0.012, mid, label, ha="left", va="center", fontsize=8.0, color="#111827", fontweight="bold")
        ax.text(x_right1 + 0.118, mid, f"{count:,}", ha="left", va="center", fontsize=7.0, color="#6B7280")

    ax.text(
        x_right1 + 0.118,
        0.94,
        "COVID cells",
        ha="left",
        va="bottom",
        fontsize=8,
        color="#6B7280",
        fontweight="bold",
    )

    ax.text(
        0.245,
        0.055,
        "Flows are scaled by COVID cell counts after applying the predefined healthy HCA-to-OneK1K-compatible mapping.",
        ha="left",
        va="center",
        fontsize=7.5,
        color="#4B5563",
    )
    ax.text(
        0.245,
        0.025,
        "Platelet and RBC-like compartments are excluded from this correspondence panel and documented separately as HCA-specific annotations.",
        ha="left",
        va="center",
        fontsize=7.5,
        color="#4B5563",
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / f"Figure_S2B_hca_cluster_onek1k_correspondence_sankey.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def panel_b_sankey_reference():
    """Render the correspondence panel in the compact table-alluvial layout."""
    set_style()
    main = pd.read_csv(CONSENSUS, sep="\t")
    covid = pd.read_csv(COVID_COUNTS, sep="\t")
    df = main.merge(
        covid[["full_clustering", "covid_cells"]],
        left_on="HCA cluster",
        right_on="full_clustering",
        how="left",
    )
    df["covid_cells"] = df["covid_cells"].fillna(0).astype(int)
    df["label"] = df["Main OneK1K-compatible label"].apply(clean_label)
    df.loc[df["Main confidence"].isin(["HCA-specific", "Not applicable"]), "label"] = "HCA-specific"
    df = df[
        (~df["HCA cluster"].isin(["Platelets", "RBC"]))
        & (~df["label"].isin(["HCA-specific", "Not mapped"]))
    ].copy()
    df["cluster_order"] = df["HCA cluster"].map({v: i for i, v in enumerate(HCA_CLUSTER_ORDER)})
    df["label_order"] = df["label"].map({v: i for i, v in enumerate(ONEK1K_CELL_ORDER)})
    df = df.sort_values(["label_order", "cluster_order", "HCA cluster"]).reset_index(drop=True)

    # Square-root compression keeps rare clusters visible while retaining the
    # abundance hierarchy; exact COVID-cell totals are shown in the right column.
    weights = df["covid_cells"].clip(lower=1).pow(0.5)
    weight_map = dict(zip(df["HCA cluster"], weights))

    n_rows = len(df)
    left_top, left_bottom = 0.895, 0.105
    row_step = (left_top - left_bottom) / (n_rows - 1)
    max_weight = max(weight_map.values())
    max_left_h = row_step * 0.76
    min_left_h = 0.0028
    left_pos = {}
    for i, row in df.iterrows():
        center = left_top - i * row_step
        h = min_left_h + (max_left_h - min_left_h) * weight_map[row["HCA cluster"]] / max_weight
        left_pos[row["HCA cluster"]] = (center - h / 2, center + h / 2)

    present_labels = [label for label in ONEK1K_CELL_ORDER if (df["label"] == label).any()]
    right_top, right_bottom = 0.895, 0.105
    group_gap = 0.010
    usable_height = right_top - right_bottom - group_gap * (len(present_labels) - 1)
    right_scale = usable_height / weights.sum()
    right_pos = {}
    label_ranges = {}
    y = right_top
    for label in present_labels:
        sub = df[df["label"] == label]
        block_top = y
        for _, row in sub.iterrows():
            h = weight_map[row["HCA cluster"]] * right_scale
            right_pos[row["HCA cluster"]] = (y - h, y)
            y -= h
        block_low = y
        label_ranges[label] = (block_low, block_top, int(sub["covid_cells"].sum()))
        y -= group_gap

    fig, ax = plt.subplots(figsize=(10.8, 9.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.052, 0.982,
        "HCA cluster correspondence to OneK1K-compatible cell types",
        fontsize=12, fontweight="bold", ha="left", va="top",
    )

    header_fc = "#EEF1F5"
    header_ec = "#E3E7EC"
    ax.add_patch(
        FancyBboxPatch(
            (0.052, 0.914), 0.775, 0.035,
            boxstyle="round,pad=0.003,rounding_size=0.006",
            fc=header_fc, ec=header_ec, lw=0.5,
        )
    )
    ax.add_patch(
        FancyBboxPatch(
            (0.838, 0.914), 0.145, 0.035,
            boxstyle="round,pad=0.003,rounding_size=0.006",
            fc=header_fc, ec=header_ec, lw=0.5,
        )
    )
    ax.text(0.205, 0.9315, "HCA cluster", fontsize=10, fontweight="bold", ha="center", va="center")
    ax.text(0.675, 0.9315, "OneK1K-compatible label", fontsize=10, fontweight="bold", ha="center", va="center")
    ax.text(0.9105, 0.9315, "COVID cells (n)", fontsize=10, fontweight="bold", ha="center", va="center")

    x_left0, x_left1 = 0.188, 0.211
    x_right0, x_right1 = 0.721, 0.744
    for _, row in df.iterrows():
        l0, l1 = left_pos[row["HCA cluster"]]
        r0, r1 = right_pos[row["HCA cluster"]]
        color = LABEL_COLORS.get(row["label"], "#9CA3AF")
        bezier_band(ax, x_left1, x_right0, l0, l1, r0, r1, color=color, alpha=0.67)
        ax.add_patch(
            Rectangle((x_left0, l0), x_left1 - x_left0, l1 - l0,
                      fc=color, ec="white", lw=0.45, zorder=3)
        )
        ax.add_patch(
            Rectangle((x_right0, r0), x_right1 - x_right0, r1 - r0,
                      fc=color, ec="white", lw=0.45, zorder=3)
        )
        ax.text(
            x_left0 - 0.010, (l0 + l1) / 2, row["HCA cluster"],
            ha="right", va="center", fontsize=7.1, color="#172033",
        )

    for label, (low, high, count) in label_ranges.items():
        center = (low + high) / 2
        count_box_h = max(0.026, min(0.066, high - low))
        ax.text(
            x_right1 + 0.012, center, label,
            ha="left", va="center", fontsize=8.8,
            color="#172033", fontweight="bold",
        )
        ax.add_patch(
            FancyBboxPatch(
                (0.842, center - count_box_h / 2), 0.136, count_box_h,
                boxstyle="round,pad=0.004,rounding_size=0.006",
                fc="#F3F5F8", ec="none", zorder=0,
            )
        )
        ax.text(0.910, center, f"{count:,}", ha="center", va="center", fontsize=8.5, color="#172033")

    ax.add_patch(
        FancyBboxPatch(
            (0.045, 0.020), 0.932, 0.062,
            boxstyle="round,pad=0.008,rounding_size=0.010",
            fc="#EEF1F5", ec="none",
        )
    )
    ax.text(
        0.058, 0.060,
        "Flow widths are square-root scaled by COVID cell abundance after applying the predefined healthy HCA-to-OneK1K-compatible mapping.",
        ha="left", va="center", fontsize=7.5, color="#46536A",
    )
    ax.text(
        0.058, 0.037,
        "Platelet, RBC-like, and other HCA-specific compartments are excluded from this correspondence panel and documented separately.",
        ha="left", va="center", fontsize=7.5, color="#46536A",
    )

    output_base = OUT / "Figure_S2_hca_cluster_onek1k_correspondence"
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(output_base.with_suffix(f".{ext}"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    panel_a()
    panel_b()
    panel_b_sankey_reference()
    print(f"Wrote Supplementary Figure S2A-B panels to {OUT}")


if __name__ == "__main__":
    main()
