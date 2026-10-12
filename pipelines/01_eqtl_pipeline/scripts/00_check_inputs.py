#!/usr/bin/env python3
"""Check local inputs for the OneK1K PBMC Science14 eQTL project."""

from __future__ import annotations

import csv
import gzip
import sys
import zipfile
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
PATHS_FILE = PROJECT_DIR / "config" / "project_paths.tsv"
QC_DIR = PROJECT_DIR / "results" / "qc"
REPORT_DIR = PROJECT_DIR / "reports"


REQUIRED_KEYS = [
    "raw_counts_tar",
    "cellxgene_h5ad",
    "science_tables_xlsx",
    "science_sm_txt",
    "imputation_dir",
    "genotype_sample_map",
    "genotype_pc_table",
    "gene_annotation",
]


def read_paths(path: Path) -> dict[str, Path]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: Path(row["path"]) for row in reader}


def read_header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return handle.readline().rstrip("\n").split("\t")


def chr_zip_status(imputation_dir: Path) -> tuple[int, list[str], list[str]]:
    present = []
    missing = []
    for chrom in range(1, 23):
        zip_path = imputation_dir / f"chr_{chrom}.zip"
        if not zip_path.exists():
            missing.append(str(chrom))
            continue
        try:
            with zipfile.ZipFile(zip_path) as zf:
                names = set(zf.namelist())
            expected = {f"chr{chrom}.dose.vcf.gz", f"chr{chrom}.info.gz"}
            if expected.issubset(names):
                present.append(str(chrom))
            else:
                missing.append(f"{chrom}:missing_expected_members")
        except zipfile.BadZipFile:
            missing.append(f"{chrom}:bad_zip")
    return len(present), present, missing


def main() -> int:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    paths = read_paths(PATHS_FILE)
    rows = []
    failures = []

    for key in REQUIRED_KEYS:
        path = paths.get(key)
        exists = bool(path and path.exists())
        rows.append(
            {
                "check": f"path:{key}",
                "status": "PASS" if exists else "FAIL",
                "detail": str(path) if path else "missing_from_project_paths.tsv",
            }
        )
        if not exists:
            failures.append(f"Missing path for {key}: {path}")

    imputation_dir = paths.get("imputation_dir")
    if imputation_dir and imputation_dir.exists():
        n_chr, present, missing = chr_zip_status(imputation_dir)
        status = "PASS" if n_chr == 22 else "FAIL"
        rows.append(
            {
                "check": "imputation_chr_zip_archives",
                "status": status,
                "detail": f"present={','.join(present)}; missing_or_bad={','.join(missing) if missing else 'none'}",
            }
        )
        if status != "PASS":
            failures.append("Not all chromosome imputation zip archives contain expected dose/info files.")

    expected_headers = {
        "genotype_sample_map": {"SampleID", "donor_id"},
        "genotype_pc_table": {"donor_id", "SampleID", "PC1", "PC10"},
        "gene_annotation": {"Probe", "chr", "start", "end", "gene_symbol"},
    }
    for key, expected in expected_headers.items():
        path = paths.get(key)
        if not path or not path.exists():
            continue
        header = set(read_header(path))
        missing = sorted(expected - header)
        status = "PASS" if not missing else "FAIL"
        rows.append(
            {
                "check": f"header:{key}",
                "status": status,
                "detail": "missing=none" if not missing else "missing=" + ",".join(missing),
            }
        )
        if missing:
            failures.append(f"{key} is missing columns: {', '.join(missing)}")

    methods_path = paths.get("science_sm_txt")
    if methods_path and methods_path.exists():
        text = methods_path.read_text(errors="replace")
        phrases = [
            "average expression of each gene per person",
            "sex, age, six genotype-based",
            "two PEER factors",
            "+/- 1Mb",
        ]
        missing = [phrase for phrase in phrases if phrase not in text]
        rows.append(
            {
                "check": "science_method_text_evidence",
                "status": "PASS" if not missing else "WARN",
                "detail": "all_key_phrases_found" if not missing else "missing_phrases=" + "; ".join(missing),
            }
        )

    out_tsv = QC_DIR / "step0_input_inventory.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    report = REPORT_DIR / "STEP0_INPUT_CHECK.md"
    n_pass = sum(row["status"] == "PASS" for row in rows)
    n_warn = sum(row["status"] == "WARN" for row in rows)
    n_fail = sum(row["status"] == "FAIL" for row in rows)
    lines = [
        "# Step 0 Input Check",
        "",
        f"- Checks passed: {n_pass}",
        f"- Warnings: {n_warn}",
        f"- Failures: {n_fail}",
        f"- Inventory table: `{out_tsv.relative_to(PROJECT_DIR)}`",
        "",
        "## Notes",
        "",
        "- The Science supplementary methods support a donor-level per-cell-type expression strategy.",
        "- Validate input formats before running the full autosomal analysis.",
    ]
    if failures:
        lines.extend(["", "## Failures", ""])
        lines.extend([f"- {failure}" for failure in failures])
    report.write_text("\n".join(lines) + "\n")

    print(f"Wrote {out_tsv}")
    print(f"Wrote {report}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
