# Result 4 independent meta-GWAS validation pipeline

This folder contains the code used for Result 4: validation of FinnGen-prioritized SMR-coloc signals in the FinnGen-independent VTE meta-GWAS.

## Inputs

Set the project root before running:

```bash
export PROJECT_DIR=/path/to/onek1k_pbmc_eqtl
```

Expected resources:

- FinnGen discovery Tier 1-3 candidate table under `results/smr_analysis/qc/`
- FinnGen-independent VTE meta-GWAS summary statistics harmonized to the OneK1K LD reference:
  `gwas/harmonized/oneK1K_ld/vte_discovery_finngen_independent_ivw_ukb_approx_logor.primary_ncohort3plus.onek1k_ld_harmonized.smr.ma`
- OneK1K cell-type eQTL SMR BESD files under `results/smr/`
- OneK1K LD reference files under `results/ld_reference/onek1k/`
- Local SMR binary under `tools/smr/source/SMR/build/RelWithDebInfo/smr`
- R packages: `data.table`, `coloc`, `ggplot2`, `grid`

## Main analysis settings

- Discovery set: FinnGen VTE SMR-coloc Tier 1-3 signals.
- Validation set: FinnGen-independent VTE meta-GWAS with `n_cohorts >= 3`.
- SMR validation threshold: same direction, cell-wise BH-FDR < 0.05, and HEIDI P >= 0.01.
- Colocalization method: `coloc.abf`.
- Overall replication status: SMR-HEIDI validation plus coloc-positive support.

## Workflow

Run targeted SMR validation of FinnGen Tier 1-3 signals in the independent meta-GWAS:

```bash
python3 scripts/00_prepare_replication_candidates.py
python3 scripts/01_run_independent_meta_gwas_smr_validation.py
python3 scripts/01b_adjust_replication_smr_fdr.py
```

Run coloc.abf for the same Tier 1-3 signal set using independent meta-GWAS summary statistics:

```bash
Rscript scripts/02_run_independent_meta_gwas_coloc_abf.R
```

Generate Figure 5:

```bash
Rscript figures/01_plot_figure5_tier_stratified_replication_forest.R
```

Generate Supplementary Table S9:

```bash
Rscript tables/01_make_table_s9_independent_meta_gwas_validation.R
```

## Main outputs

- `manuscript_results/completed_results/result4_meta_gwas_validation/main/Figure5_tier_stratified_replication_forest.png`
- `manuscript_results/completed_results/result4_meta_gwas_validation/main/Figure5_tier_stratified_replication_forest.pdf`
- `manuscript_results/completed_results/result4_meta_gwas_validation/supplementary/Table_S9_independent_meta_gwas_validation_of_finngen_prioritized_smr_coloc_signals.tsv`

## Source outputs

- `manuscript_results/completed_results/result4_meta_gwas_validation/source_tables/Table_result4_meta_gwas_smr_validation_full.tsv`
- `manuscript_results/completed_results/result4_meta_gwas_validation/source_tables/Table_result4_meta_gwas_coloc_abf_full.tsv`
- `manuscript_results/completed_results/result4_meta_gwas_validation/source_tables/Table_result4_tier_stratified_replication_forest_values.tsv`

## Key Result 4 summary

- Tier 1: 4/4 signals showed independent meta-GWAS SMR-HEIDI validation; 3/4 also showed coloc-positive support.
- Tier 2: 2/3 signals showed independent meta-GWAS SMR-HEIDI validation; 2/3 also showed coloc-positive support.
- Tier 3: 2/2 signals showed independent meta-GWAS SMR-HEIDI validation; 1/2 also showed coloc-positive support.
