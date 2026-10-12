#!/usr/bin/env python3
"""Create METAL-ready input files from standardized discovery GWAS tables."""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
GWAS_DIR = PROJECT_DIR / "gwas"
FORMATTED_DIR = GWAS_DIR / "formatted"
METAL_INPUT_DIR = GWAS_DIR / "meta" / "input"
QC_DIR = GWAS_DIR / "qc"
REPORT_DIR = GWAS_DIR / "reports"

DATASETS = [
    "decode_chb_cvdc_dbds",
    "decode_iceland",
    "decode_intermountain",
    "ukb_i9_vte",
]

METAL_COLUMNS = ["MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value", "Freq1", "Weight"]


def marker_name(row: dict[str, str]) -> str:
    a = sorted([row["A1"], row["A2"]])
    return f"{row['CHR']}:{row['BP']}:{a[0]}:{a[1]}"


def as_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return float("inf")


def prepare_one(dataset: str) -> dict[str, int | str]:
    in_path = FORMATTED_DIR / f"{dataset}.standardized.tsv.gz"
    out_path = METAL_INPUT_DIR / f"{dataset}.metal.txt.gz"
    best: dict[str, dict[str, str]] = {}
    raw_rows = 0
    duplicate_marker_rows = 0

    with gzip.open(in_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw_rows += 1
            marker = marker_name(row)
            previous = best.get(marker)
            if previous is not None:
                duplicate_marker_rows += 1
                if as_float(row["p"]) >= as_float(previous["p"]):
                    continue
            best[marker] = row

    with gzip.open(out_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=METAL_COLUMNS)
        writer.writeheader()
        for marker in sorted(best, key=lambda x: tuple([int(x.split(":")[0]), int(x.split(":")[1]), x])):
            row = best[marker]
            writer.writerow(
                {
                    "MarkerName": marker,
                    "Allele1": row["A1"],
                    "Allele2": row["A2"],
                    "Effect": row["b"],
                    "StdErr": row["se"],
                    "P-value": row["p"],
                    "Freq1": row["freq"],
                    "Weight": row["n"],
                }
            )

    return {
        "dataset_id": dataset,
        "standardized_rows": raw_rows,
        "metal_rows": len(best),
        "duplicate_marker_rows_removed": duplicate_marker_rows,
        "metal_input": str(out_path.relative_to(PROJECT_DIR)),
    }


def write_qc(rows: list[dict[str, int | str]]) -> None:
    path = QC_DIR / "metal_input_prepare_summary.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    report = REPORT_DIR / "METAL_INPUT_PREPARE_QC.md"
    lines = [
        "# METAL Input Prepare QC",
        "",
        "| Dataset | Standardized rows | METAL rows | Duplicate marker rows removed | METAL input |",
        "|---|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['dataset_id']} | {row['standardized_rows']} | {row['metal_rows']} | "
            f"{row['duplicate_marker_rows_removed']} | `{row['metal_input']}` |"
        )
    lines += [
        "",
        "Marker construction:",
        "",
        "- `MarkerName = CHR:BP:allele_low:allele_high`.",
        "- `Allele1` remains the cohort effect allele.",
        "- `Allele2` remains the cohort non-effect allele.",
        "- Within each cohort, duplicate `MarkerName` rows are collapsed by keeping the row with the smallest P value.",
        "- This keeps marker IDs stable across cohorts while allowing METAL to flip effects when Allele1/Allele2 are reversed.",
    ]
    report.write_text("\n".join(lines) + "\n")


def main() -> int:
    METAL_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for dataset in DATASETS:
        print(f"Preparing METAL input for {dataset}", flush=True)
        rows.append(prepare_one(dataset))
    write_qc(rows)
    print(f"Wrote METAL inputs to {METAL_INPUT_DIR}")
    print(f"Wrote QC to {QC_DIR / 'metal_input_prepare_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
