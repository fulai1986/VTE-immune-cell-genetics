#!/usr/bin/env python3
"""Plot HCA COVID expression detectability for replicated Tier 1 signals."""

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
WORK = BASE / "covid_validation/stephenson_hca/result5_expression_context"
MAPPING_DIR = BASE / "covid_validation/hca_onek1k_mapping"

EXPR_FILE = WORK / "results/covid_d0_tier1_3_same_cell_gene_expression.curated_mapping.tsv.gz"
REPL_FILE = (
    BASE
    / "manuscript_results/completed_results/result4_meta_gwas_validation/source_tables/"
    / "Table_result4_tier_stratified_replication_forest_values.tsv"
)
MAP_FILE = MAPPING_DIR / "results/hca_to_onek1k_mapping_main_with_seurat_sensitivity.tsv"

OUT_DIR = WORK / "figures"
TAB_DIR = WORK / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

OUT_PREFIX = "figure5a_replicated_tier1_hca_expression_qc"


def format_pair(cell: str, gene: str) -> str:
    return f"{cell.replace('_', ' ')} - {gene}"


def main() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 9,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    repl = pd.read_csv(REPL_FILE, sep="\t")
    tier1_repl = repl.query(
        "tier == 'Tier 1' and validation_and_coloc_positive == True"
    ).copy()
    tier1_repl["pair"] = tier1_repl.apply(
        lambda x: format_pair(x["cell"], x["Gene"]), axis=1
    )

    desired_order = ["CD8 ET - SLC44A2", "CD8 ET - FADS1", "NK - FADS1"]
    tier1_repl["pair"] = pd.Categorical(
        tier1_repl["pair"], categories=desired_order, ordered=True
    )
    tier1_repl = tier1_repl.sort_values("pair")

    expr = pd.read_csv(EXPR_FILE, sep="\t")
    expr["pair"] = expr.apply(lambda x: format_pair(x["cell_label"], x["gene"]), axis=1)
    expr = expr[expr["pair"].isin(desired_order)].copy()
    expr["pair"] = pd.Categorical(expr["pair"], categories=desired_order, ordered=True)
    expr["pct_expr_percent"] = expr["pct_expr"] * 100.0

    plot_table = (
        expr.groupby("pair", observed=True)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_cells=("n_target_cells", "median"),
            median_mean_expr=("mean_expr", "median"),
            q1_mean_expr=("mean_expr", lambda x: np.quantile(x, 0.25)),
            q3_mean_expr=("mean_expr", lambda x: np.quantile(x, 0.75)),
            median_pct_expr=("pct_expr_percent", "median"),
            q1_pct_expr=("pct_expr_percent", lambda x: np.quantile(x, 0.25)),
            q3_pct_expr=("pct_expr_percent", lambda x: np.quantile(x, 0.75)),
        )
        .reset_index()
    )
    plot_table["detectability"] = np.where(
        (plot_table["median_mean_expr"] >= 0.1)
        & (plot_table["median_pct_expr"] >= 10.0),
        "Functionally analyzable",
        "Sparse in HCA scRNA-seq",
    )

    table_out = TAB_DIR / f"{OUT_PREFIX}_plot_data.tsv"
    plot_table.to_csv(table_out, sep="\t", index=False)

    mapping = pd.read_csv(MAP_FILE, sep="\t")
    mapping = mapping[
        mapping["Main OneK1K-compatible label"].isin(["CD8 ET", "NK"])
        & mapping["HCA cluster"].isin(["CD8.TE", "CD8.EM", "CD8.Prolif", "NK_16hi", "NK_prolif"])
    ].copy()
    mapping["cluster_label"] = mapping.apply(
        lambda x: f"{x['HCA cluster']} ({x['Main confidence']})", axis=1
    )

    colors = {
        "CD8 ET - SLC44A2": "#9b1b30",
        "CD8 ET - FADS1": "#4b6f9e",
        "NK - FADS1": "#4b6f9e",
    }
    status_colors = {
        "Functionally analyzable": "#9b1b30",
        "Sparse in HCA scRNA-seq": "#6f7782",
    }

    fig = plt.figure(figsize=(7.1, 4.6), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=3,
        height_ratios=[0.95, 2.35],
        width_ratios=[1.2, 1.65, 1.65],
        hspace=0.42,
        wspace=0.56,
    )

    ax_map = fig.add_subplot(gs[0, :])
    ax_map.axis("off")
    ax_map.text(
        0.00,
        0.97,
        "HCA Healthy-to-OneK1K-compatible annotation bridge",
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
        transform=ax_map.transAxes,
    )
    ax_map.text(
        0.00,
        0.70,
        "Healthy HCA clusters were mapped to OneK1K-compatible labels, then fixed for COVID donor-level pseudobulk expression.",
        ha="left",
        va="top",
        fontsize=8.8,
        color="#30343b",
        transform=ax_map.transAxes,
    )

    box_y = {"CD8 ET": 0.30, "NK": 0.30}
    x0s = {"CD8 ET": 0.05, "NK": 0.55}
    for label in ["CD8 ET", "NK"]:
        x0 = x0s[label]
        ax_map.text(
            x0,
            0.46,
            label,
            ha="left",
            va="center",
            fontsize=10,
            fontweight="bold",
            color="#20242a",
            transform=ax_map.transAxes,
        )
        sub = mapping[mapping["Main OneK1K-compatible label"] == label]
        cluster_text = "  |  ".join(sub["cluster_label"].tolist())
        ax_map.text(
            x0 + 0.11,
            box_y[label],
            cluster_text,
            ha="left",
            va="center",
            fontsize=8.2,
            color="#30343b",
            transform=ax_map.transAxes,
            bbox=dict(
                boxstyle="round,pad=0.28,rounding_size=0.08",
                fc="#f1f3f5",
                ec="#d9dee5",
                lw=0.8,
            ),
        )

    ax_mean = fig.add_subplot(gs[1, 0:2])
    ax_pct = fig.add_subplot(gs[1, 2], sharey=ax_mean)

    pairs = desired_order
    y_positions = np.arange(len(pairs))[::-1]
    rng = np.random.default_rng(18)

    for idx, pair in enumerate(pairs):
        y = y_positions[idx]
        sub = expr[expr["pair"].astype(str) == pair]
        jitter = rng.normal(0, 0.045, size=len(sub))
        c = colors[pair]

        ax_mean.scatter(
            sub["mean_expr"],
            np.full(len(sub), y) + jitter,
            s=17,
            color=c,
            alpha=0.33 if pair != "CD8 ET - SLC44A2" else 0.42,
            linewidth=0,
            rasterized=True,
        )
        pt = plot_table[plot_table["pair"].astype(str) == pair].iloc[0]
        ax_mean.errorbar(
            pt["median_mean_expr"],
            y,
            xerr=[
                [pt["median_mean_expr"] - pt["q1_mean_expr"]],
                [pt["q3_mean_expr"] - pt["median_mean_expr"]],
            ],
            fmt="o",
            ms=7.5,
            color=c,
            mec="black",
            mew=0.7,
            capsize=3.5,
            lw=1.6,
            zorder=5,
        )

        ax_pct.scatter(
            sub["pct_expr_percent"],
            np.full(len(sub), y) + jitter,
            s=17,
            color=c,
            alpha=0.33 if pair != "CD8 ET - SLC44A2" else 0.42,
            linewidth=0,
            rasterized=True,
        )
        ax_pct.errorbar(
            pt["median_pct_expr"],
            y,
            xerr=[
                [pt["median_pct_expr"] - pt["q1_pct_expr"]],
                [pt["q3_pct_expr"] - pt["median_pct_expr"]],
            ],
            fmt="o",
            ms=7.5,
            color=c,
            mec="black",
            mew=0.7,
            capsize=3.5,
            lw=1.6,
            zorder=5,
        )
        ax_pct.text(
            41.5,
            y,
            pt["detectability"],
            ha="left",
            va="center",
            fontsize=8.7,
            color=status_colors[pt["detectability"]],
            fontweight="bold" if pt["detectability"] == "Functionally analyzable" else "normal",
            clip_on=False,
        )

    y_labels = [
        r"CD8 ET - $\it{SLC44A2}$",
        r"CD8 ET - $\it{FADS1}$",
        r"NK - $\it{FADS1}$",
    ]
    ax_mean.set_yticks(y_positions)
    ax_mean.set_yticklabels(y_labels)
    ax_mean.set_xlabel("Mean expression per donor")
    ax_mean.set_xlim(-0.01, 0.55)
    ax_mean.set_title("Expression abundance", loc="left", fontsize=10.5, fontweight="bold")

    ax_pct.set_xlabel("Positive-cell fraction (%)")
    ax_pct.set_xlim(-1, 40)
    ax_pct.set_title("Expression prevalence", loc="left", fontsize=10.5, fontweight="bold")
    plt.setp(ax_pct.get_yticklabels(), visible=False)
    ax_pct.tick_params(axis="y", length=0)
    ax_pct.axvline(5, ls=(0, (3, 3)), lw=0.9, color="#9aa1aa")
    ax_pct.text(5.5, 2.45, "5%", ha="left", va="center", fontsize=7.5, color="#6f7782")

    for ax in [ax_mean, ax_pct]:
        ax.set_ylim(-0.55, len(pairs) - 0.45)
        ax.grid(axis="x", color="#e2e5ea", lw=0.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)

    for idx, pair in enumerate(pairs):
        y = y_positions[idx]
        pt = plot_table[plot_table["pair"].astype(str) == pair].iloc[0]
        ax_mean.text(
            0.535,
            y + 0.18,
            f"n={int(pt['n_samples'])}, median cells={int(round(pt['median_cells']))}",
            ha="right",
            va="bottom",
            fontsize=7.6,
            color="#555b64",
        )

    fig.text(0.018, 0.985, "A", ha="left", va="top", fontsize=18, fontweight="bold")
    fig.suptitle(
        "HCA COVID-19 expression detectability of replicated Tier 1 VTE regulatory signals",
        x=0.05,
        y=0.995,
        ha="left",
        fontsize=11.5,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.17, right=0.88, top=0.88, bottom=0.14)

    for ext, dpi in [("png", 450), ("pdf", 450), ("svg", 450)]:
        fig.savefig(OUT_DIR / f"{OUT_PREFIX}.{ext}", dpi=dpi)
    plt.close(fig)

    print(f"Wrote: {OUT_DIR / (OUT_PREFIX + '.png')}")
    print(f"Wrote: {OUT_DIR / (OUT_PREFIX + '.pdf')}")
    print(f"Wrote: {OUT_DIR / (OUT_PREFIX + '.svg')}")
    print(f"Wrote: {table_out}")


if __name__ == "__main__":
    main()
