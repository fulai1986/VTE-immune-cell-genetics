#!/usr/bin/env python3
"""Approximate log(OR) effects for OpenGWAS ukb-d-I9_VTE.

The local OpenGWAS VCF for ukb-d-I9_VTE is labelled StudyType=Continuous.
This script preserves the original allele orientation but uses the signed
Z-score and case/control counts to approximate a liability/logistic-scale
log(OR) and SE:

    Z = ES / SE
    SE_logOR = sqrt((1/N_case + 1/N_control) / (2 * EAF * (1 - EAF)))
    beta_logOR = Z * SE_logOR

This is a sensitivity conversion, not a substitute for native logistic GWAS.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
GWAS_DIR = PROJECT_DIR / "gwas"
INFILE = GWAS_DIR / "formatted" / "ukb_i9_vte.standardized.tsv.gz"
OUTFILE = GWAS_DIR / "formatted" / "ukb_i9_vte.approx_logor.standardized.tsv.gz"
QC_FILE = GWAS_DIR / "qc" / "ukb_i9_vte.approx_logor_conversion_qc.tsv"
TARGET_FILE = GWAS_DIR / "qc" / "ukb_i9_vte.approx_logor_target_snp_check.tsv"
REPORT = GWAS_DIR / "reports" / "UKB_I9_VTE_APPROX_LOGOR_CONVERSION_QC.md"

N_CASE = 4620
N_CONTROL = 356574

TARGET_POSITIONS = {
    ("10", "71218646"): "rs17490626",
    ("11", "61560081"): "rs174538",
    ("11", "61567753"): "rs174544",
    ("11", "61579760"): "rs174555",
    ("11", "61581656"): "rs174559",
    ("11", "61582708"): "rs174561",
    ("19", "10724441"): "rs78579480",
    ("19", "10727014"): "rs10420809",
    ("22", "38557326"): "rs84473",
}


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def main() -> int:
    GWAS_DIR.joinpath("qc").mkdir(parents=True, exist_ok=True)
    GWAS_DIR.joinpath("reports").mkdir(parents=True, exist_ok=True)

    counts = {
        "input_rows": 0,
        "written_rows": 0,
        "skip_bad_numeric": 0,
        "skip_bad_frequency": 0,
        "skip_zero_se": 0,
    }
    targets: list[dict[str, str]] = []

    with gzip.open(INFILE, "rt", newline="") as fin, gzip.open(OUTFILE, "wt", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames or []
        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            counts["input_rows"] += 1
            beta = parse_float(row["b"])
            se = parse_float(row["se"])
            freq = parse_float(row["freq"])
            if not all(math.isfinite(x) for x in (beta, se, freq)):
                counts["skip_bad_numeric"] += 1
                continue
            if se <= 0:
                counts["skip_zero_se"] += 1
                continue
            if freq <= 0 or freq >= 1:
                counts["skip_bad_frequency"] += 1
                continue

            z = beta / se
            se_logor = math.sqrt((1.0 / N_CASE + 1.0 / N_CONTROL) / (2.0 * freq * (1.0 - freq)))
            beta_logor = z * se_logor

            out = dict(row)
            out["b"] = f"{beta_logor:.8g}"
            out["se"] = f"{se_logor:.8g}"
            writer.writerow(out)
            counts["written_rows"] += 1

            key = (row["CHR"], row["BP"])
            if key in TARGET_POSITIONS:
                targets.append(
                    {
                        "expected_rsid": TARGET_POSITIONS[key],
                        "SNP": row["SNP"],
                        "CHR": row["CHR"],
                        "BP": row["BP"],
                        "A1": row["A1"],
                        "A2": row["A2"],
                        "freq": row["freq"],
                        "original_ES": row["b"],
                        "original_SE": row["se"],
                        "Z": f"{z:.8g}",
                        "approx_logOR": f"{beta_logor:.8g}",
                        "approx_SE_logOR": f"{se_logor:.8g}",
                        "p": row["p"],
                    }
                )

    with QC_FILE.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in counts.items():
            writer.writerow([key, value])
        writer.writerow(["n_case", N_CASE])
        writer.writerow(["n_control", N_CONTROL])

    with TARGET_FILE.open("w", newline="") as handle:
        fieldnames = [
            "expected_rsid",
            "SNP",
            "CHR",
            "BP",
            "A1",
            "A2",
            "freq",
            "original_ES",
            "original_SE",
            "Z",
            "approx_logOR",
            "approx_SE_logOR",
            "p",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(targets, key=lambda x: (int(x["CHR"]), int(x["BP"]))))

    lines = [
        "# UKB I9_VTE approximate log(OR) conversion QC",
        "",
        "## Conclusion",
        "",
        "Generated an approximate log(OR)-scale sensitivity file for `ukb-d-I9_VTE` using signed Z-scores, case/control counts, and effect allele frequency.",
        "This file should be treated as sensitivity input, not as a native logistic GWAS replacement.",
        "",
        "## Inputs and outputs",
        "",
        f"- Input standardized file: `{INFILE.relative_to(PROJECT_DIR)}`",
        f"- Output approximate log(OR) file: `{OUTFILE.relative_to(PROJECT_DIR)}`",
        f"- QC table: `{QC_FILE.relative_to(PROJECT_DIR)}`",
        f"- Target SNP check: `{TARGET_FILE.relative_to(PROJECT_DIR)}`",
        "",
        "## Formula",
        "",
        "`Z = ES / SE`",
        "",
        "`SE_logOR = sqrt((1/N_case + 1/N_control) / (2 * EAF * (1 - EAF)))`",
        "",
        "`beta_logOR = Z * SE_logOR`",
        "",
        f"N_case = {N_CASE:,}; N_control = {N_CONTROL:,}.",
    ]
    REPORT.write_text("\n".join(lines) + "\n")

    print(f"Wrote {OUTFILE}")
    print(f"Wrote {QC_FILE}")
    print(f"Wrote {TARGET_FILE}")
    print(f"Wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
