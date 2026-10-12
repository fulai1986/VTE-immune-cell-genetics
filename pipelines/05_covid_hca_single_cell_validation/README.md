# Result 5 COVID/HCA single-cell validation pipeline

This directory contains the finalized analysis and figure-generation code for Result 5 and Figure 6. The archive follows the same structure as the completed Result 1-4 pipelines.

## Analysis scope

- HCA Healthy cells establish the bridge from HCA `full_clustering` labels to OneK1K-compatible cell types.
- COVID cells do not define cell identity; the fixed Healthy-derived mapping is applied to COVID cells.
- The primary annotation is a lineage-aware consensus mapping. Marker-panel Seurat label transfer is retained as a sensitivity analysis.
- Platelets remain an HCA-specific downstream effector compartment and are not forced into the OneK1K 14-cell taxonomy.
- Functional follow-up is restricted to CD8 ET-SLC44A2 after strict genetic replication and HCA expression-detectability filtering.
- The main exposure is donor-level CD8 ET pseudobulk SLC44A2 expression. Low and High groups are the lower and upper tertiles; the middle tertile is omitted from group contrasts.
- LIANA is run per donor with the consensus resource, a 5% expression threshold, and at least 20 cells per sender/receiver group.

Set `PROJECT_DIR` to the analysis workspace before running any script. Code and data may reside in separate directories. Required external resources are listed in `config/project_paths.example.tsv`.

Python dependencies are listed in `config/python_requirements.txt`. The Seurat sensitivity step additionally requires R packages `Matrix`, `Seurat`, and `data.table`.

## Directory structure

- `mapping/`: Healthy-HCA to OneK1K-compatible mapping and Seurat sensitivity.
- `analysis/`: candidate detectability, Hallmark scoring, donor-level LIANA, and platelet programs.
- `figures/`: Figure 6 panel generators and Supplementary Figures S2-S5.
- `tables/`: Supplementary Tables S10-S13 builders.

## Mapping workflow

```bash
python3 mapping/01_build_onek1k_centroids_and_map_hca_healthy.py
python3 mapping/02_make_lineage_aware_consensus_mapping.py
python3 mapping/03_prepare_seurat_label_transfer_inputs.py
Rscript mapping/04_run_seurat_marker_panel_sensitivity.R
python3 mapping/05_compile_main_mapping_with_seurat_sensitivity.py
python3 mapping/06_export_fixed_mapping_and_covid_counts.py
```

The centroid-based script constructs balanced OneK1K reference profiles and maps Healthy HCA cells using shared expression features. The lineage-aware step combines centroid evidence, marker sanity checks, and lineage identity. Seurat marker-panel transfer is reported as sensitivity evidence rather than replacing the primary mapping.

## Main analysis workflow

```bash
python3 analysis/01_prepare_covid_donor_pseudobulk.py
python3 analysis/02_build_candidate_detectability_qc.py
python3 analysis/03_score_cd8et_hallmark_programs.py
python3 analysis/04_run_donor_level_liana.py
python3 analysis/05_score_platelet_hallmark_programs.py
python3 analysis/06_build_high_low_hallmark_contrasts.py
python3 analysis/07_summarize_liana_high_low_contrasts.py
```

The final grouped analyses use two-sided Mann-Whitney tests with BH-FDR within each prespecified outcome family. LIANA effect intervals are bootstrap 95% confidence intervals for the raw High-minus-Low median difference. The two retained communication metrics are the number of consensus-supported LR pairs passing the predefined significance threshold and summed LIANA rank evidence among those pairs.

## Figure workflow

```bash
python3 figures/01_plot_figure6a_candidate_detectability.py
python3 figures/02_plot_figure6b_workflow.py
python3 figures/03_plot_figure6b_expression_grouping.py
python3 figures/04_plot_figure6c_cd8et_hallmarks.py
python3 figures/05_plot_figure6d_liana_distributions.py
python3 figures/06_plot_figure6e_liana_effects.py
python3 figures/07_plot_figure6f_platelet_hallmarks.py
python3 figures/08_plot_supp_figure_s2_mapping.py
python3 figures/09_plot_supp_figure_s2a_workflow.py
python3 analysis/08_continuous_severity_sensitivity.py
python3 tables/05_build_table_s13_severity_adjusted.py
python3 figures/10_plot_figure_s3_cd8et_severity.py
python3 figures/11_plot_figure_s4_liana_severity.py
python3 figures/12_plot_figure_s5_platelet_severity.py
```

The final Figure 6 manuscript composite was assembled from statistical panel exports and author-prepared workflow artwork. This code repository generates the component panels and source values. The submitted composite accompanies the manuscript.

## Supplementary tables

```bash
python3 tables/01_build_table_s10_hca_onek1k_mapping.py
python3 tables/02_build_table_s11_candidate_detectability.py
python3 analysis/08_continuous_severity_sensitivity.py
python3 tables/05_build_table_s13_severity_adjusted.py
python3 figures/10_plot_figure_s3_cd8et_severity.py
python3 figures/11_plot_figure_s4_liana_severity.py
python3 figures/12_plot_figure_s5_platelet_severity.py
python3 tables/03_prepare_table_s12_figure6_results.py
python3 tables/04_format_table_s12_figure6_results.py
```

- Table S10: full HCA-to-OneK1K-compatible mapping and Seurat sensitivity evidence.
- Table S11: strict Tier 1 replication and HCA detectability evidence.
- Table S12: all grouped results underlying Figure 6C-F.
- Table S13: all 65 continuous associations, before and after collection-day severity adjustment.

The XLSX formatter uses the public `openpyxl` package. TSV files are the source tables.

## Final outputs

- Figure 6 component exports under `covid_validation/stephenson_hca/result5_expression_context/figures/` in the workspace.
- `supplementary/Figure_S2_hca_cluster_onek1k_correspondence.{png,pdf,svg}`
- `manuscript_results/current_sensitivity/Figure_S3_cd8et_hallmark_severity_sensitivity.{png,pdf}`
- `manuscript_results/current_sensitivity/Figure_S4_liana_severity_adjusted_associations.{png,pdf,svg}`
- `manuscript_results/current_sensitivity/Figure_S5_platelet_hallmark_severity_sensitivity.{png,pdf}`
- `supplementary/Table_S10_hca_to_onek1k_compatible_mapping.{tsv,xlsx}`
- `supplementary/Table_S11_replicated_tier1_candidate_detectability_in_hca_covid.{tsv,xlsx}`
- `manuscript_results/current_sensitivity/Table_S13_continuous_slc44a2_severity_adjusted.{tsv,xlsx}`
- `supplementary/Table_S12_complete_slc44a2_high_low_results_underlying_figure6.{tsv,xlsx}`

## Continuous severity sensitivity

The statistical unit is the donor. All eligible D0 donors are included (84 CD8
ET, 54 LIANA, 76 platelet), including the middle expression tertile. Exposure and
outcome midranks are residualized on an intercept and four indicators for the
five collection-day severity categories. Residual correlations are partial
Spearman coefficients; two-sided t P values use n-6 degrees of freedom. BH FDR
is calculated separately across 50 CD8 ET, six LIANA and nine platelet readouts.
LIANA percentile CIs use 5,000 within-severity donor bootstrap draws with the
original fixed seeds and donor order. Unadjusted CIs use the Fisher-z method.
The script checks donor uniqueness, outcome counts and valid severity categories.
No site, age, sex or detection-depth adjustment is included in this analysis.
Scores use the original standardized donor expression algorithm, not GSVA/GSEA.
