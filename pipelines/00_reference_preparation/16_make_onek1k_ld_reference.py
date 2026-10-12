#!/usr/bin/env python3
"""Create OneK1K PLINK LD reference files from imputed VCF zip archives.

The eQTL pipeline uses VCF DS as ALT-allele dosage. PLINK2 imports
dosage=DS and writes .bim allele columns as ALT, REF for these VCFs, which
keeps the LD reference allele order aligned with the SMR export plan.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
CONFIG_DIR = PROJECT_DIR / "config"
METADATA_DIR = PROJECT_DIR / "metadata"
OUT_DIR = PROJECT_DIR / "results" / "ld_reference" / "onek1k"
PLINK2 = Path(os.environ.get("PLINK2_BIN", PROJECT_DIR / "tools" / "plink2" / "plink2"))


def read_paths() -> dict[str, Path]:
    with (CONFIG_DIR / "project_paths.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: Path(row["path"]) for row in reader}


def read_passwords() -> tuple[str, str]:
    try:
        return os.environ["IMPUTE_PASS_A"], os.environ["IMPUTE_PASS_B"]
    except KeyError as exc:
        raise RuntimeError("Set IMPUTE_PASS_A and IMPUTE_PASS_B before reading imputation archives.") from exc


def password_for_chrom(chrom: int, pass_a: str, pass_b: str) -> str:
    return pass_a if chrom in {1, 10, 11, 12} else pass_b


def write_keep_file() -> Path:
    sample_map = METADATA_DIR / "imputed_vcf_sample_to_donor_map.tsv"
    keep_path = OUT_DIR / "keep_mapped_vcf_samples.tsv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with sample_map.open() as inp, keep_path.open("w", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for row in reader:
            if str(row["has_donor_id"]).lower() in {"true", "1"}:
                writer.writerow(["0", row["vcf_sample"]])
    return keep_path


def run_command(cmd: list[str], log_path: Path | None = None) -> None:
    print(" ".join(cmd[:1] + ["<args omitted>"] if cmd and cmd[0] == "unzip" else cmd), flush=True)
    with (log_path.open("w") if log_path else subprocess.DEVNULL) as log:
        proc = subprocess.run(cmd, stdout=log if log_path else None, stderr=subprocess.STDOUT if log_path else None)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)


def extract_vcf(chrom: int, imputation_dir: Path, tmp_dir: Path, pass_a: str, pass_b: str, overwrite: bool) -> Path:
    tmp_vcf = tmp_dir / f"chr{chrom}.dose.vcf.gz"
    if tmp_vcf.exists() and not overwrite:
        return tmp_vcf
    tmp_dir.mkdir(parents=True, exist_ok=True)
    zip_path = imputation_dir / f"chr_{chrom}.zip"
    member = f"chr{chrom}.dose.vcf.gz"
    cmd = ["unzip", "-P", password_for_chrom(chrom, pass_a, pass_b), "-p", str(zip_path), member]
    print(f"Extracting chr{chrom} dose VCF to temporary file", flush=True)
    with tmp_vcf.open("wb") as out:
        proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"Failed to extract {member} from {zip_path}: {proc.stderr.decode(errors='replace')}")
    return tmp_vcf


def make_ld_reference(chrom: int, args: argparse.Namespace, keep_path: Path, pass_a: str, pass_b: str) -> dict[str, str | int]:
    paths = read_paths()
    chrom_dir = OUT_DIR / f"chr{chrom}"
    tmp_dir = OUT_DIR / "tmp"
    qc_dir = OUT_DIR / "qc"
    chrom_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    id_suffix = ".allele_id" if args.allele_specific_ids else ""
    out_prefix = chrom_dir / f"onek1k.chr{chrom}.r2_{args.min_r2:g}.maf_{args.min_maf:g}{id_suffix}"
    bed = Path(str(out_prefix) + ".bed")
    bim = Path(str(out_prefix) + ".bim")
    fam = Path(str(out_prefix) + ".fam")
    if bed.exists() and bim.exists() and fam.exists() and not args.overwrite:
        return {"chrom": chrom, "status": "exists", "prefix": str(out_prefix)}

    tmp_vcf = extract_vcf(chrom, paths["imputation_dir"], tmp_dir, pass_a, pass_b, args.overwrite)
    log_path = qc_dir / f"chr{chrom}.plink2.log"
    expr = f"(R2 >= {args.min_r2}) && (MAF >= {args.min_maf})"
    cmd = [
        str(PLINK2),
        "--vcf",
        str(tmp_vcf),
        "dosage=DS",
        "--const-fid",
        "0",
        "--keep",
        str(keep_path),
        "--import-max-alleles",
        "2",
        "--snps-only",
        "just-acgt",
    ]
    if args.allele_specific_ids:
        cmd.extend(["--set-all-var-ids", "@:#:$r:$a"])
    cmd.extend(
        [
        "--extract-if-info",
        expr,
        "--maf",
        str(args.min_maf),
        "--threads",
        str(args.threads),
        "--make-bed",
        "--out",
        str(out_prefix),
        ]
    )
    print(f"Running PLINK2 chr{chrom}", flush=True)
    run_command(cmd, log_path)
    if args.delete_tmp:
        tmp_vcf.unlink(missing_ok=True)
    return {"chrom": chrom, "status": "written", "prefix": str(out_prefix)}


def write_manifest(rows: list[dict[str, str | int]], args: argparse.Namespace) -> Path:
    manifest = OUT_DIR / "ld_reference_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        fieldnames = ["chrom", "status", "prefix", "min_r2", "min_maf", "sample_keep_file", "plink2"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out.update(
                {
                    "min_r2": args.min_r2,
                    "min_maf": args.min_maf,
                    "sample_keep_file": str(OUT_DIR / "keep_mapped_vcf_samples.tsv"),
                    "plink2": str(PLINK2),
                }
            )
            writer.writerow(out)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom", type=int, action="append", help="Chromosome to process; may be repeated. Default: 1-22.")
    parser.add_argument("--min-r2", type=float, default=0.8, help="Minimum INFO/R2 filter.")
    parser.add_argument("--min-maf", type=float, default=0.01, help="Minimum MAF filter after sample filtering.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-tmp", action="store_true", help="Keep temporary extracted VCF.gz files.")
    parser.add_argument(
        "--pos-only-ids",
        dest="allele_specific_ids",
        action="store_false",
        help="Keep original chr:pos variant IDs. Not recommended for SMR because duplicate positions can occur.",
    )
    parser.set_defaults(allele_specific_ids=True)
    args = parser.parse_args()
    args.delete_tmp = not args.keep_tmp

    if not PLINK2.exists():
        raise FileNotFoundError(f"PLINK2 not found: {PLINK2}")

    pass_a, pass_b = read_passwords()
    keep_path = write_keep_file()
    chroms = args.chrom or list(range(1, 23))
    rows = [make_ld_reference(chrom, args, keep_path, pass_a, pass_b) for chrom in chroms]
    manifest = write_manifest(rows, args)
    print(f"Wrote manifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
