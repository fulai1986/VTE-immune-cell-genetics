#!/usr/bin/env python3
"""Inspect encrypted Michigan imputation VCF archives without writing passwords to outputs."""

from __future__ import annotations

import csv
import gzip
import contextlib
import re
import subprocess
import os
import sys
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
CONFIG_DIR = PROJECT_DIR / "config"
QC_DIR = PROJECT_DIR / "results" / "qc"
REPORT_DIR = PROJECT_DIR / "reports"


def read_paths() -> dict[str, Path]:
    with (CONFIG_DIR / "project_paths.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: Path(row["path"]) for row in reader}


def read_existing_passwords() -> tuple[str, str]:
    env_a = os.environ.get("IMPUTE_PASS_A")
    env_b = os.environ.get("IMPUTE_PASS_B")
    if env_a and env_b:
        return env_a, env_b
    password_source = os.environ.get("IMPUTE_PASSWORD_SOURCE")
    if not password_source:
        raise RuntimeError("Set IMPUTE_PASS_A/IMPUTE_PASS_B or IMPUTE_PASSWORD_SOURCE before reading encrypted VCF archives.")
    text = Path(password_source).read_text()
    pass_a = re.search(r'^PASS_A\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    pass_b = re.search(r'^PASS_B\s*=\s*"([^"]+)"', text, flags=re.MULTILINE)
    if not pass_a or not pass_b:
        raise RuntimeError(f"Could not read PASS_A/PASS_B from {password_source}")
    return pass_a.group(1), pass_b.group(1)


def password_for_chrom(chrom: int, pass_a: str, pass_b: str) -> str:
    return pass_a if chrom in {1, 10, 11, 12} else pass_b


def stream_zip_member(zip_path: Path, member: str, password: str):
    unzip = subprocess.Popen(
        ["unzip", "-P", password, "-p", str(zip_path), member],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert unzip.stdout is not None
    gz = gzip.GzipFile(fileobj=unzip.stdout)
    return unzip, gz


def stop_process(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    with contextlib.suppress(Exception):
        proc.terminate()
        proc.wait(timeout=5)
    if proc.poll() is None:
        with contextlib.suppress(Exception):
            proc.kill()
            proc.wait(timeout=5)


def inspect_chrom(chrom: int, imputation_dir: Path) -> dict[str, object]:
    pass_a, pass_b = read_existing_passwords()
    password = password_for_chrom(chrom, pass_a, pass_b)
    zip_path = imputation_dir / f"chr_{chrom}.zip"

    unzip, handle = stream_zip_member(zip_path, f"chr{chrom}.dose.vcf.gz", password)
    n_meta = 0
    info_r2 = "not_found"
    samples = []
    try:
        for raw in handle:
            line = raw.decode("utf-8", errors="replace").rstrip("\n")
            if line.startswith("##"):
                n_meta += 1
                if line.startswith("##INFO=<ID=R2"):
                    info_r2 = line
                continue
            if line.startswith("#CHROM"):
                parts = line.split("\t")
                samples = parts[9:]
                break
    finally:
        with contextlib.suppress(Exception):
            handle.close()
        stop_process(unzip)
    if not samples:
        stderr = unzip.stderr.read() if unzip.stderr else b""
        raise RuntimeError(stderr.decode("utf-8", errors="replace"))

    unzip_info, info_handle = stream_zip_member(zip_path, f"chr{chrom}.info.gz", password)
    info_header = ""
    first_variant = ""
    try:
        info_header = info_handle.readline().decode("utf-8", errors="replace").rstrip("\n")
        first_variant = info_handle.readline().decode("utf-8", errors="replace").rstrip("\n")
    finally:
        with contextlib.suppress(Exception):
            info_handle.close()
        stop_process(unzip_info)

    return {
        "chrom": chrom,
        "zip_path": str(zip_path),
        "dose_member": f"chr{chrom}.dose.vcf.gz",
        "info_member": f"chr{chrom}.info.gz",
        "vcf_meta_lines_seen": n_meta,
        "vcf_sample_count": len(samples),
        "first_sample": samples[0] if samples else "",
        "last_sample": samples[-1] if samples else "",
        "info_r2_defined": info_r2 != "not_found",
        "info_header": info_header,
        "first_info_variant": first_variant,
    }


def main(argv: list[str]) -> int:
    chrom = int(argv[1]) if len(argv) > 1 else 1
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    paths = read_paths()
    result = inspect_chrom(chrom, paths["imputation_dir"])

    sample_map = pd.read_csv(paths["genotype_sample_map"], sep="\t", dtype={"SampleID": str, "donor_id": str})
    result["genotype_sample_map_rows"] = len(sample_map)
    result["genotype_sample_map_unique_donors"] = sample_map["donor_id"].nunique()

    out = QC_DIR / f"chr{chrom}_imputed_vcf_inventory.tsv"
    pd.DataFrame([result]).to_csv(out, sep="\t", index=False)

    report = REPORT_DIR / f"STEP4_CHR{chrom}_IMPUTED_VCF_CHECK.md"
    report.write_text(
        "\n".join(
            [
                f"# Step 4 chr{chrom} Imputed VCF Check",
                "",
                f"- VCF samples: {result['vcf_sample_count']}",
                f"- First sample: `{result['first_sample']}`",
                f"- Last sample: `{result['last_sample']}`",
                f"- INFO R2 field present: {result['info_r2_defined']}",
                f"- Genotype sample map rows: {result['genotype_sample_map_rows']}",
                f"- Unique mapped donors: {result['genotype_sample_map_unique_donors']}",
                f"- Inventory table: `{out.relative_to(PROJECT_DIR)}`",
                "",
                "Passwords were supplied through the local environment and were not written to this report.",
            ]
        )
        + "\n"
    )
    print(f"Wrote {out}")
    print(f"Wrote {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
