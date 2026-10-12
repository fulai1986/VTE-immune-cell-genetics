#!/usr/bin/env python3
"""Plot Supplementary Figure S2A: healthy-bridge cell mapping workflow."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
OUT = PROJECT / "manuscript_results/completed_results/result5_covid_single_cell_validation/supplementary"
OUT.mkdir(parents=True, exist_ok=True)


def set_style():
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def draw_panel():
    set_style()
    fig, ax = plt.subplots(figsize=(10.8, 4.9))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ink = "#172033"
    muted = "#566277"
    border = "#D9DEE6"
    arrow_color = "#778292"
    burgundy = "#A51D45"

    def card(x, y, w, h, title, body, fc, ec=border, title_color=ink):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.008,rounding_size=0.014",
                fc=fc,
                ec=ec,
                lw=0.9,
            )
        )
        ax.text(
            x + w / 2,
            y + h * 0.73,
            title,
            ha="center",
            va="center",
            fontsize=10.0,
            fontweight="bold",
            color=title_color,
        )
        ax.text(
            x + w / 2,
            y + h * 0.35,
            body,
            ha="center",
            va="center",
            fontsize=8.0,
            color=muted,
            linespacing=1.35,
        )

    def section_header(x, y, w, text):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                0.060,
                boxstyle="round,pad=0.003,rounding_size=0.008",
                fc="#EEF1F5",
                ec="#E3E7EC",
                lw=0.5,
            )
        )
        ax.text(
            x + w / 2,
            y + 0.030,
            text,
            ha="center",
            va="center",
            fontsize=9.3,
            fontweight="bold",
            color=ink,
        )

    def flow_arrow(x0, y0, x1, y1, color=arrow_color):
        ax.add_patch(
            FancyArrowPatch(
                (x0, y0),
                (x1, y1),
                arrowstyle="-|>",
                mutation_scale=13,
                color=color,
                lw=1.7,
                connectionstyle="arc3,rad=0.0",
            )
        )

    ax.text(0.012, 0.985, "A", fontsize=18, fontweight="bold", ha="left", va="top")
    ax.text(
        0.073,
        0.982,
        "Healthy-reference bridge for OneK1K-compatible annotation of HCA COVID cells",
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="top",
        color=ink,
    )

    x1, w1 = 0.050, 0.225
    x2, w2 = 0.315, 0.255
    x3, w3 = 0.610, 0.170
    x4, w4 = 0.820, 0.155
    header_y = 0.845
    section_header(x1, header_y, w1, "Reference inputs")
    section_header(x2, header_y, w2, "Identity mapping in healthy cells")
    section_header(x3, header_y, w3, "Frozen correspondence")
    section_header(x4, header_y, w4, "COVID application")

    card(
        x1,
        0.585,
        w1,
        0.205,
        "OneK1K PBMC reference",
        "science14 framework: 14 cell types\n14,000 sampled reference cells\n77 shared immune-identity genes",
        "#EAF1F8",
        ec="#83A4C2",
    )
    card(
        x1,
        0.320,
        w1,
        0.205,
        "HCA healthy bridge",
        "100,593 cells from 24 samples\n45 original HCA clusters\nCOVID cells excluded",
        "#EDF5F0",
        ec="#7EAD98",
    )
    ax.text(
        x1 + w1 / 2,
        0.555,
        "+",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
        color=arrow_color,
    )

    ax.add_patch(
        FancyBboxPatch(
            (x2, 0.320),
            w2,
            0.470,
            boxstyle="round,pad=0.008,rounding_size=0.014",
            fc="#FFF6E7",
            ec="#D6A654",
            lw=1.0,
        )
    )
    ax.text(
        x2 + w2 / 2,
        0.735,
        "Lineage-aware consensus mapping",
        ha="center",
        va="center",
        fontsize=9.3,
        fontweight="bold",
        color=ink,
    )
    mapping_steps = [
        ("1", "Nearest-centroid cosine similarity"),
        ("2", "Lineage and marker concordance"),
        ("3", "Cluster-level consensus + confidence"),
    ]
    for (number, label), y in zip(mapping_steps, [0.645, 0.555, 0.465]):
        ax.add_patch(
            FancyBboxPatch(
                (x2 + 0.020, y - 0.035),
                0.035,
                0.070,
                boxstyle="round,pad=0.002,rounding_size=0.010",
                fc="#D9A443",
                ec="none",
            )
        )
        ax.text(
            x2 + 0.0375,
            y,
            number,
            ha="center",
            va="center",
            fontsize=8.2,
            fontweight="bold",
            color="white",
        )
        ax.text(x2 + 0.066, y, label, ha="left", va="center", fontsize=7.25, color=ink)

    ax.plot([x2 + 0.022, x2 + w2 - 0.022], [0.420, 0.420], color="#E4D3B1", lw=0.8)
    ax.text(
        x2 + w2 / 2,
        0.386,
        "HCA full_clustering retained",
        ha="center",
        va="center",
        fontsize=8.0,
        fontweight="bold",
        color="#76571E",
    )
    ax.add_patch(
        FancyBboxPatch(
            (x2 + 0.036, 0.337),
            w2 - 0.072,
            0.038,
            boxstyle="round,pad=0.003,rounding_size=0.008",
            fc="#F1E5CD",
            ec="none",
        )
    )
    ax.text(
        x2 + w2 / 2,
        0.356,
        "Seurat label-transfer sensitivity analysis",
        ha="center",
        va="center",
        fontsize=7.3,
        color="#76571E",
    )

    card(
        x3,
        0.470,
        w3,
        0.250,
        "Cluster-level mapping",
        "42 HCA clusters mapped\nwithin the 14-class framework\n13 labels observed in COVID",
        "#F2EDF6",
        ec="#A48AB7",
    )
    card(
        x4,
        0.470,
        w4,
        0.250,
        "Apply fixed labels",
        "COVID cells receive the\npredefined correspondence\nwithout redefining identity",
        "#FBECEF",
        ec="#C64A69",
        title_color=burgundy,
    )

    flow_arrow(x1 + w1, 0.555, x2 - 0.008, 0.555)
    flow_arrow(x2 + w2, 0.555, x3 - 0.008, 0.555)
    flow_arrow(x3 + w3, 0.555, x4 - 0.008, 0.555, color=burgundy)

    ax.add_patch(
        FancyBboxPatch(
            (x4 - 0.005, 0.755),
            w4 + 0.010,
            0.055,
            boxstyle="round,pad=0.004,rounding_size=0.010",
            fc="#F8DDE4",
            ec="none",
        )
    )
    ax.text(
        x4 + w4 / 2,
        0.782,
        "COVID did not define cell identities",
        ha="center",
        va="center",
        fontsize=7.4,
        fontweight="bold",
        color=burgundy,
    )

    section_header(0.050, 0.245, 0.925, "Downstream use")
    ax.add_patch(
        FancyBboxPatch(
            (0.050, 0.125),
            0.445,
            0.095,
            boxstyle="round,pad=0.006,rounding_size=0.012",
            fc="#F3F5F8",
            ec=border,
            lw=0.7,
        )
    )
    ax.text(0.070, 0.190, "Genetics-guided analysis", fontsize=7.9, fontweight="bold", color=ink, ha="left", va="center")
    ax.text(
        0.070, 0.148,
        "Use OneK1K-compatible labels to align\ncell-gene pairs across datasets",
        fontsize=6.9, color=muted, ha="left", va="center", linespacing=1.25,
    )

    ax.add_patch(
        FancyBboxPatch(
            (0.520, 0.125),
            0.455,
            0.095,
            boxstyle="round,pad=0.006,rounding_size=0.012",
            fc="#F3F5F8",
            ec=border,
            lw=0.7,
        )
    )
    ax.text(0.540, 0.190, "Functional and communication analysis", fontsize=7.9, fontweight="bold", color=ink, ha="left", va="center")
    ax.text(
        0.540, 0.148,
        "Retain HCA subclusters; platelets remain\nan HCA-specific effector compartment",
        fontsize=6.9, color=muted, ha="left", va="center", linespacing=1.25,
    )

    ax.add_patch(
        FancyBboxPatch(
            (0.045, 0.025),
            0.930,
            0.060,
            boxstyle="round,pad=0.007,rounding_size=0.010",
            fc="#EEF1F5",
            ec="none",
        )
    )
    ax.text(
        0.060,
        0.055,
        "Healthy HCA cells served as the annotation bridge; the correspondence was fixed before analysis of the COVID compartment.",
        fontsize=7.7,
        color="#46536A",
        ha="left",
        va="center",
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(
            OUT / f"Figure_S2A_hca_onek1k_mapping_workflow.{ext}",
            dpi=600,
            bbox_inches="tight",
        )
    plt.close(fig)


if __name__ == "__main__":
    draw_panel()
    print(f"Wrote Figure S2A to {OUT}")
