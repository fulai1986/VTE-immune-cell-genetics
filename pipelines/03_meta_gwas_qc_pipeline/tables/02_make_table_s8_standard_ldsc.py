#!/usr/bin/env python3
"""Build Table S8 from the completed standard LDSC run, without rerunning GWAS/LDSC.

The result directory can be overridden with STANDARD_LDSC_RUN_DIR. The output
directory can be overridden with TABLE_S8_OUT_DIR for a non-destructive check.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import re
from pathlib import Path


ROOT = Path(os.environ["PROJECT_DIR"]).resolve()
RUN = Path(os.environ.get("STANDARD_LDSC_RUN_DIR", ROOT / "gwas/ldsc_standard"))
OUT = Path(os.environ.get(
    "TABLE_S8_OUT_DIR",
    ROOT / "manuscript_results/completed_results/result3_meta_gwas_qc/supplementary",
))
LOG = RUN / "results/vte_meta_standard_ldsc.log"
PREP = RUN / "prepared/preparation_qc.tsv"
MUNGED = RUN / "munged/vte_meta_standard.sumstats.gz"


def log_match(lines: str, pattern: str) -> tuple[str, ...]:
    found = re.findall(pattern, lines, flags=re.MULTILINE)
    if len(found) != 1:
        raise ValueError(f"Expected one LDSC log match for {pattern!r}; found {len(found)}")
    value = found[0]
    return value if isinstance(value, tuple) else (value,)


def main() -> None:
    log = LOG.read_text()
    n_regression = int(log_match(log, r"^After merging with regression SNP LD, (\d+) SNPs remain\.$")[0])
    mean_chi2 = float(log_match(log, r"^Mean Chi\^2: ([0-9.]+)$")[0])
    lambda_gc = float(log_match(log, r"^Lambda GC: ([0-9.]+)$")[0])
    intercept, intercept_se = map(float, log_match(log, r"^Intercept: ([0-9.]+) \(([0-9.]+)\)$"))
    ratio, ratio_se = map(float, log_match(log, r"^Ratio: ([0-9.]+) \(([0-9.]+)\)$"))

    with PREP.open(newline="") as handle:
        prep = {row["metric"]: int(row["value"]) for row in csv.DictReader(handle, delimiter="\t")}
    n_display = 0
    n_munged = 0
    with gzip.open(MUNGED, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                z = float(row["Z"])
            except (TypeError, ValueError):
                continue
            if not math.isfinite(z):
                continue
            n_munged += 1
            if z ** 2 < 80:
                n_display += 1
    if n_munged != n_regression:
        raise ValueError(f"Munged SNP count {n_munged} differs from LDSC regression count {n_regression}")

    header = [
        "Analysis set", "SNPs", "Mean chi-square", "Lambda GC", "LDSC intercept",
        "LDSC intercept SE", "Attenuation ratio", "Attenuation ratio SE", "Notes",
    ]
    rows = [
        ["Primary meta-GWAS variants", prep["gwas_rows"], "", "", "", "", "", "",
         "Variants in the n_cohorts >= 3 meta-GWAS source; standard LDSC was not run on this unrestricted set."],
        ["HapMap3/EUR reference-matched variants", prep["matched_rows"], "", "", "", "", "", "",
         "Position- and allele-matched variants before standard summary-statistic processing."],
        ["Standard LDSC regression set", n_regression, f"{mean_chi2:.4f}", f"{lambda_gc:.4f}",
         f"{intercept:.4f}", f"{intercept_se:.4f}", f"{ratio:.4f}", f"{ratio_se:.4f}",
         "HapMap3/EUR-matched SNPs after munge_sumstats.py filtering; European LD reference excludes the MHC; ldsc.py --h2 used its default two-step estimator and a 200-block jackknife. No separate chi-square > 80 exclusion."],
        ["Figure 4D QQ display subset", n_display, "", "", "", "", "", "",
         "Visual subset of the standard LDSC SNP set with chi-square < 80; no separate LDSC fit. Figure 4D reports calibration metrics from the full standard LDSC regression set above."],
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "Table_S8_ldsc_calibration_summary.tsv"
    with target.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    print(target)


if __name__ == "__main__":
    main()
