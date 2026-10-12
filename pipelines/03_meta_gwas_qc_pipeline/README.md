# Result 3 meta-GWAS QC pipeline

This folder contains the code used for Result 3: construction and quality assessment of the FinnGen-independent VTE meta-GWAS used as an external validation resource.

## Inputs

Set the project root before running:

```bash
export PROJECT_DIR=/path/to/onek1k_pbmc_eqtl
```

Expected resources:

- raw or preformatted VTE GWAS summary statistics under `gwas/download/` and `gwas/formatted/`
- GWAS dataset manifest following `config/gwas_dataset_manifest.example.tsv`
- METAL executable at `tools/metal/METAL/build/metal/metal`
- OneK1K LD reference under `results/ld_reference/onek1k/`
- HapMap3/EUR LDSC reference under `tools/ldsc_reference/eur_w_ld_chr/`

## Main analysis settings

- Meta-GWAS role: FinnGen-independent external validation resource.
- Included datasets: deCODE DK, deCODE Iceland, deCODE Intermountain, and UK Biobank.
- Excluded dataset: FinnGen.
- Primary analysis: variants present in at least three cohorts (`n_cohorts >= 3`).
- Meta-analysis method: inverse-variance fixed-effect meta-analysis using METAL with `SCHEME STDERR`.
- UK Biobank handling: `ukb-d-I9_VTE` was converted from signed Z statistics to approximate log(OR) effects before meta-analysis.
- Downstream harmonization: primary meta-GWAS variants were matched to the OneK1K LD reference for SMR/coloc-compatible analyses.
- LDSC QC: HapMap3/EUR SNPs were used for standard weighted two-step LDSC calibration metrics.

## Workflow

Standardize source GWAS files:

```bash
python3 scripts/01_standardize_discovery_gwas.py
```

Prepare METAL inputs for the deCODE and original UKB formatted files:

```bash
python3 scripts/02_prepare_metal_inputs.py
```

Convert the UKB VTE summary statistics to approximate log(OR) scale:

```bash
python3 scripts/03_convert_ukb_i9_vte_to_approx_logor.py
```

Run the final FinnGen-independent meta-GWAS with UKB approximate log(OR), create the `n_cohorts >= 3` and `n_cohorts = 4` outputs, and harmonize to OneK1K LD:

```bash
python3 scripts/04_run_meta_gwas_ukb_approx_logor.py
```

Prepare HapMap3/EUR LDSC-matched summary statistics and compute LDSC QC metrics:

```bash
export LDSC_DIR=/path/to/CBIIT-ldsc
export LDSC_PYTHON=/path/to/ldsc-python3.9/bin/python
python3 qc/03_run_standard_ldsc.py
```

Generate Figure 4 and Result 3 source tables:

```bash
Rscript tables/01_make_result3_supplementary_tables.R
python3 tables/02_make_table_s8_standard_ldsc.py
Rscript figures/01_plot_figure4_meta_gwas_qc.R
```

## Main outputs

- `manuscript_results/completed_results/result3_meta_gwas_qc/main/Figure4_meta_gwas_construction_qc_final.jpg`
- `manuscript_results/completed_results/result3_meta_gwas_qc/supplementary/Table_S6_meta_gwas_cohort_information.tsv`
- `manuscript_results/completed_results/result3_meta_gwas_qc/supplementary/Table_S7_meta_gwas_variant_coverage_and_qc_summary.tsv`
- `manuscript_results/completed_results/result3_meta_gwas_qc/supplementary/Table_S8_ldsc_calibration_summary.tsv`

## Key Result 3 values

- Total meta-GWAS sample size: 41,870 cases and 933,074 controls.
- Primary `n_cohorts >= 3` variants: 6,891,026.
- OneK1K LD-matched variants: 5,783,976.
- Genome-wide significant variants: 3,428.
- Genome-wide significant variants with HetP > 0.01: 2,626.
- HapMap3/EUR SNPs used for LDSC after QC: 1,117,715.
- LDSC intercept: 1.0308 +/- 0.0089.
- Attenuation ratio: 0.1127 +/- 0.0326.

## Standard LDSC provenance and reuse

The verified run used CBIIT/ldsc, branch ldsc39, commit
`6c673952cee74bd5c57aef1555a03b1c015399a0` (reported LDSC 3.0.1), Python 3.9,
NumPy 1.23.5, pandas 1.5.3, SciPy 1.10.1 and bitarray 2.6.0. Install this external
tool in a separate environment. The runner reads the primary meta-GWAS,
allele-matches HapMap3/European references, munges with signed Z and merge-alleles,
and uses the same EUR LD scores for reference and weights. The reference excludes
the MHC. Default sample-size filtering, two-step cutoff 30 and 200 jackknife
blocks are retained. No extra chi-square > 80 exclusion is applied to the fit.
The QQ plot uses the munged SNPs with chi-square < 80 (1,117,556 for this dataset).
Total case-plus-control N limits interpretation of observed-scale heritability.

The default run directory is `gwas/ldsc_standard/`; override it with
`STANDARD_LDSC_RUN_DIR`. The runner refuses nonempty run directories. To reuse
completed results, run `meta-tables` instead of `meta-qc`. The S8 exporter reads
preparation counts, the LDSC log and the munged SNP file; it does not hard-code
calibration estimates. Figure 4 reads the same standard run and S8 source TSV.
The author-arranged manuscript Figure 4 image is preserved separately.
