# Result 2 SMR and coloc pipeline

This folder contains the code used for Result 2: FinnGen VTE SMR discovery, coloc.abf support, evidence tier assignment, and Figure 3/Table 1 generation.

## Inputs

Set the project root before running:

```bash
export PROJECT_DIR=/path/to/onek1k_pbmc_eqtl
```

Expected project resources are listed in `config/project_paths.example.tsv`. The workflow assumes that cell-type eQTL files have already been exported in SMR-compatible format and that the FinnGen VTE GWAS has been harmonized to GRCh37/OneK1K LD SNPs.

## Software

Required command-line tools and packages:

- SMR executable
- Python packages: `pandas`, `python-docx`
- R packages: `data.table`, `coloc`, `ggplot2`, `grid`, `png`, `scales`

## Key analysis settings

- SMR discovery GWAS: FinnGen R13 VTE.
- Cell types: 14 PBMC classes from `config/science14_celltypes.tsv`.
- SMR settings: `--smr`, `--peqtl-smr 0.05`, `--diff-freq-prop 0.99`.
- Multiple testing: BH-FDR within each cell type.
- Discovery threshold: cell-wise SMR FDR < 0.05.
- HEIDI support threshold: P_HEIDI >= 0.01.
- coloc.abf support threshold: PP4 >= 0.8 and PP4 / (PP3 + PP4) >= 0.8.
- Evidence tiers:
  - Tier 1: SMR FDR < 0.05, P_HEIDI >= 0.01, and coloc.abf positive.
  - Tier 2: SMR FDR < 0.05, P_HEIDI < 0.01, and coloc.abf positive.
  - Tier 3: SMR FDR < 0.05, P_HEIDI >= 0.01, but coloc.abf not positive.
  - Not prioritized: remaining SMR FDR < 0.05 pairs.

## Workflow

Run SMR per cell type:

```bash
tail -n +2 "$PROJECT_DIR/config/science14_celltypes.tsv" | cut -f1 | sed '/^$/d; s/ /_/g' | while read -r cell; do
  python3 scripts/01_run_finngen_discovery_smr_for_cell.py --cell "$cell" --thread-num 4
done
```

Summarize SMR discovery hits:

```bash
python3 scripts/02_summarize_finngen_discovery_smr_all_cells.py
```

Run coloc.abf for SMR FDR-positive pairs and assign tiers:

```bash
Rscript coloc/01_run_coloc_abf_finngen_discovery_hits.R
Rscript coloc/02_assign_evidence_tiers.R
```

Generate Table 1 and supplementary tables:

```bash
python3 tables/01_make_top_snp_rsid_mapping.py
Rscript tables/02_make_table_figure3_smr_coloc_integrated.R
python3 tables/03_make_table1_smr_coloc_three_line.py
Rscript tables/04_make_result2_supplementary_tables.R
```

Generate Figure 3:

```bash
Rscript figures/01_plot_figure3a_ranked_stacked_bar.R
Rscript figures/02_plot_figure3b_tier1_3_forest.R
Rscript figures/03_plot_figure3c_cross_cell_smr_heatmap_effect_lines.R
Rscript figures/04_plot_figure3d_tier1_coloc_reference_style.R
Rscript figures/05_assemble_figure3_final.R
```

Main outputs are written to:

- `manuscript_results/figures/figure3_smr_coloc_final.png`
- `manuscript_results/figures/figure3_smr_coloc_final.pdf`
- `manuscript_results/tables/Table_figure3_smr_coloc_prioritized_tier1_3.tsv`
- `manuscript_results/tables/supplementary/Table_S4_finngen_smr_coloc_all_fdr05_signals.tsv`
- `manuscript_results/tables/supplementary/Table_S5_finngen_smr_discovery_by_cell.tsv`
