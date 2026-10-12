#!/usr/bin/env python3
"""Build VCF sample to donor map for imputed OneK1K dosage files."""

from __future__ import annotations

import csv
import gzip
import os
import re
import subprocess
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
CONFIG_DIR = PROJECT_DIR / "config"
METADATA_DIR = PROJECT_DIR / "metadata"
QC_DIR = PROJECT_DIR / "results" / "qc"
REPORT_DIR = PROJECT_DIR / "reports"


def read_paths() -> dict[str, Path]:
    with (CONFIG_DIR / "project_paths.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: Path(row["path"]) for row in reader}


def read_password_a() -> str:
    env_a = os.environ.get("IMPUTE_PASS_A")
    if env_a:
        return env_a
    password_source = os.environ.get("IMPUTE_PASSWORD_SOURCE")
    if not password_source:
        raise RuntimeError("Set IMPUTE_PASS_A or IMPUTE_PASSWORD_SOURCE before reading encrypted VCF archives.")
    text = Path(password_source).read_text()
    match = re.search(r'^PASS_A\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    if not match:
        raise RuntimeError(f"Could not read PASS_A from {password_source}")
    return match.group(1)


def read_chr1_samples(imputation_dir: Path) -> list[str]:
    cmd = [
        "unzip",
        "-P",
        read_password_a(),
        "-p",
        str(imputation_dir / "chr_1.zip"),
        "chr1.dose.vcf.gz",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    assert proc.stdout is not None
    handle = gzip.GzipFile(fileobj=proc.stdout)
    try:
        for raw in handle:
            line = raw.decode("utf-8", errors="replace").rstrip("\n")
            if line.startswith("#CHROM"):
                return line.split("\t")[9:]
    finally:
        try:
            handle.close()
        except Exception:
            pass
        proc.terminate()
    raise RuntimeError("No #CHROM header found in chr1 dose VCF")


def normalize_sample_id(vcf_sample: str) -> str:
    value = vcf_sample.removeprefix("OK1K_")
    if re.fullmatch(r"\d+", value):
        return str(int(value))
    if re.fullmatch(r"\d+_2", value):
        left, suffix = value.split("_", 1)
        return f"{int(left)}_{suffix}"
    return value


def main() -> int:
    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    paths = read_paths()

    vcf_samples = read_chr1_samples(paths["imputation_dir"])
    vcf_map = pd.DataFrame({"vcf_sample": vcf_samples})
    vcf_map["SampleID"] = vcf_map["vcf_sample"].map(normalize_sample_id)

    sample_map = pd.read_csv(paths["genotype_sample_map"], sep="\t", dtype=str)
    out = vcf_map.merge(sample_map[["SampleID", "donor_id", "geo_gender", "gsm_accession"]], on="SampleID", how="left")
    out["has_donor_id"] = out["donor_id"].notna() & (out["donor_id"] != "")

    out_path = METADATA_DIR / "imputed_vcf_sample_to_donor_map.tsv"
    out.to_csv(out_path, sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "vcf_samples", "value": len(out)},
            {"metric": "mapped_sample_ids", "value": int(out["donor_id"].notna().sum())},
            {"metric": "mapped_nonempty_donor_ids", "value": int(out["has_donor_id"].sum())},
            {"metric": "unique_mapped_donors", "value": int(out.loc[out["has_donor_id"], "donor_id"].nunique())},
            {"metric": "samples_without_donor_id", "value": int((~out["has_donor_id"]).sum())},
        ]
    )
    summary_path = QC_DIR / "imputed_vcf_sample_map_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    report_path = REPORT_DIR / "STEP4_IMPUTED_VCF_SAMPLE_MAP.md"
    report_path.write_text(
        "\n".join(
            [
                "# Step 4 Imputed VCF Sample Map",
                "",
                f"- VCF samples: {len(out)}",
                f"- Samples with non-empty donor_id: {int(out['has_donor_id'].sum())}",
                f"- Unique mapped donors: {int(out.loc[out['has_donor_id'], 'donor_id'].nunique())}",
                f"- Samples without donor_id: {int((~out['has_donor_id']).sum())}",
                f"- Map: `{out_path.relative_to(PROJECT_DIR)}`",
                f"- Summary: `{summary_path.relative_to(PROJECT_DIR)}`",
                "",
                "Samples without donor_id are retained in the map for provenance but should be excluded from eQTL dosage matrices.",
            ]
        )
        + "\n"
    )
    print(f"Wrote {out_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
