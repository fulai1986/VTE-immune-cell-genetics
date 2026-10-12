# Release validation

## Scope

The publication package was assembled from the archived Result 1-5 code.
Original analysis files and manuscript results were not overwritten. Changes
cover workspace paths, entry points, dependency documentation, missing export
steps, and public-library table formatting. Scientific thresholds and the
archived bootstrap algorithm were retained.

Source syntax and workflow references were checked for 49 Python, 29 R, and
three Bash pipeline files across 17 workflow stages.
Five automated tests cover workspace initialization and preservation, numerical
helpers, tertile assignment, fixed-mapping export, and workflow settings.
Twenty-two analysis, table, and plotting scripts were executed successfully in
a separate validation workspace using existing intermediate inputs.

## Numerical comparisons

The following regenerated outputs matched their archived counterparts within
floating-point tolerance:

- Candidate-set replication FDR for the nine prioritized associations.
- All 50 CD8 ET Hallmark High-versus-Low contrasts.
- All nine selected platelet Hallmark contrasts.
- HCA cluster-to-OneK1K COVID cell counts.
- Figure 5 replication source values.
- Archived unadjusted continuous-expression sensitivity results (formerly S12).

Table S11 retains full numerical precision in regenerated exports, whereas the
archived table rounds FDR, HEIDI P, PP4, expression, and detectability summaries.
These are presentation-precision differences, not different analysis estimates.

## Outstanding exact-reproduction discrepancy

The archived LIANA summary script, with its existing 20,000 bootstrap draws and
fixed seeds, does not exactly reproduce two stored confidence-interval lower
bounds for CD8 ET to monocyte:

| Metric | Archived lower bound | Regenerated lower bound |
| --- | ---: | ---: |
| Significant LR-pair count | -37.0 | -39.0 |
| Summed rank evidence | -19.57662430045808 | -19.939511089292186 |

The upper bounds, medians, median differences, P values, FDR values, and other
four communication contrasts match. Both affected intervals continue to include
zero. Current Table S12 inherits these same two lower-bound differences. The precise
origin of the discrepancy has not been established; it should not be attributed
to random variation without verifying the original computation. No seed or
statistical parameter was adjusted to force agreement with stored results.

The original figure and tables remain unchanged. Exact regeneration of those
two reported endpoints still requires reconciliation before claiming complete
numerical reproducibility. Do not replace the manuscript outputs with rerun
exports without reviewing this difference.

## Limits

This was not a complete raw-data rerun of genome-wide eQTL mapping, meta-GWAS,
all SMR/colocalization jobs, or single-cell communication inference. Installation
in a clean external environment and the hosted CI run have not yet been tested.
No GitHub upload has been performed. The package contains code, configuration,
and documentation, not participant-level data or third-party binaries.

## Current manuscript synchronization, 2026-10-09

The release now uses standard LDSC, grouped Table S12, continuous severity-adjusted
Table S13, and Figures S3-S5. Validation used an isolated workspace with read-only
links to existing donor summaries; source analysis files were not overwritten.

- All 65 S13 outcomes matched the current manuscript TSV, including unadjusted and
  adjusted rho, P, FDR, sample counts and all six adjusted LIANA bootstrap CIs.
  FDR differences were at most 1.36e-20 (floating-point rounding).
- All three S3-S5 PNG exports were byte-identical to the current manuscript images.
- S8 was re-exported from the completed standard LDSC log, preparation QC and
  munged SNP file; every TSV cell matched the current S8 source table.
- The renamed S12 builder reproduced all 65 rows from stored grouped contrasts.
  This does not resolve the two historical bootstrap endpoints noted above when
  recomputing grouped contrasts from donor metrics.
- Seven unit tests passed, including independent residualization and S8 checks
  distinguishing the full regression from the chi-square < 80 display subset.
- Syntax checks passed for 55 Python, 28 R and three Bash pipeline scripts; 19
  workflow stages reference existing scripts.

The standard LDSC runner was connected and syntax-checked; the expensive original
GWAS/LDSC analysis was not repeated. Current statistical figures were rebuilt from
retained inputs. No GitHub upload or commit was performed.

Figure 4 was also regenerated successfully from existing primary meta-GWAS and
standard LDSC files in the isolated workspace. The updated plot reads standard
S8 estimates and the actual munged/display SNP sets. Its PDF exports now use
Quartz on macOS or Cairo elsewhere to preserve symbols. The publication's
author-arranged composite was not overwritten.

The supplementary-table packager was run on the current 13 workbooks and its
archive retained all workbook bytes exactly.
