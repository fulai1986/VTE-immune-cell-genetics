# Part 1: OneK1K PBMC pseudobulk eQTL pipeline

This folder contains the code used for the first results section: cell-type-specific donor-level pseudobulk cis-eQTL mapping in OneK1K PBMCs, eGene definition, OneK1K concordance QC, and Figure 2/Supplementary Figure S1 generation.

The analysis is organized as a staged workflow. Raw or controlled-access data are not included.

## Directory layout

- `config/`: example path configuration.
- `scripts/`: core eQTL pipeline.
- `qc/`: OneK1K concordance benchmark and supplementary table generation.
- `figures/`: Figure 2 and Figure S1 plotting scripts.

## Required inputs

Set `PROJECT_DIR` to the analysis workspace. Create `$PROJECT_DIR/config/project_paths.tsv` from `config/project_paths.example.tsv` and update paths for the local environment. Required resources include:

- OneK1K raw count archive.
- OneK1K CELLxGENE h5ad or equivalent cell metadata source.
- Science supplementary tables / annotation references used to define the 14 PBMC cell types.
- GENCODE v19 gene annotation.
- Michigan imputation server chromosome zip archives.
- Genotype sample-to-donor map.
- Genotype PC table.
- OneK1K reference eSNP-eGene tables for the concordance benchmark.

Encrypted imputation archives should be accessed using environment variables rather than storing passwords in code:

```bash
export IMPUTE_PASS_A='...'
export IMPUTE_PASS_B='...'
```

Alternatively, set `IMPUTE_PASSWORD_SOURCE` to a local helper file containing `PASS_A` and `PASS_B`. Do not commit that helper file.

## Core eQTL workflow

Run from this folder unless otherwise noted.

```bash
python3 scripts/00_check_inputs.py
python3 scripts/01_make_science14_metadata.py
bash scripts/03_run_all_pseudobulk.sh
python3 scripts/04_summarize_pseudobulk.py
python3 scripts/05_prepare_eqtl_covariates.py
python3 scripts/06_inspect_imputed_vcf.py
python3 scripts/07_build_vcf_sample_map.py
bash scripts/09_run_full_eqtl_queue.sh
python3 ../00_reference_preparation/16_make_onek1k_ld_reference.py --allele-specific-ids
python3 ../00_reference_preparation/18_prepare_smr_qfile_besd.py --make-besd
bash scripts/11_run_egene_eigenmt_queue.sh
Rscript scripts/13_check_egene_eigenmt_status.R
```

The cis-eQTL model uses donor-level pseudobulk expression within each cell type. Expression and imputed ALT dosage are residualized on covariates, followed by linear association testing within a +/-1 Mb cis-window around gene boundaries. Variant filters are `R2 >= 0.8` and `MAF >= 0.01`.

eGenes are defined by a Meff/eigenMT-like gene-level correction:

```text
P_gene = min(P_lead x M_eff, 1)
```

followed by cell-type-wise BH-FDR < 0.05.

## QC and supplementary outputs

```bash
Rscript qc/01_prepare_onek1k_concordance_benchmark.R
Rscript qc/02_make_eqtl_supplementary_tables.R
```

## Figure generation

```bash
Rscript figures/01_plot_figure2a_sample_basis.R
Rscript figures/02_plot_figure2b_egene_counts.R
Rscript figures/03_plot_figure2c_egene_sharing_upset.R
Rscript figures/04_plot_figure2d_egene_jaccard_heatmap.R
Rscript figures/05_plot_figure2e_onek1k_concordance_heatmap.R
Rscript figures/06_assemble_figure2.R
Rscript figures/07_plot_figure_s1_random_concordance_heatmap.R
```

## Main outputs

- `results/pseudobulk/`: donor-level pseudobulk counts and rank-inverse-normalized expression.
- `metadata/covariates/`: cell-type-specific donor covariates.
- `results/eqtl/`: full cis-eQTL results and top association per gene.
- `manuscript_results/tables/egene_eigenmt/`: gene-level eGene calls.
- `manuscript_results/figures/figure2.*`: Figure 2.
- `manuscript_results/figures/supplementary/Figure_S1_onek1k_eqtl_concordance_benchmark.*`: Figure S1.
- `manuscript_results/tables/supplementary/`: Supplementary Tables S1-S3.

## Notes

This workflow uses donor-level pseudobulk eQTL mapping rather than a single-cell mixed model. eQTL beta values are ALT dosage effects.
