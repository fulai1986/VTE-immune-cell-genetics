#!/usr/bin/env python3
"""Standardize discovery VTE GWAS files before METAL meta-analysis.

Outputs one gzipped table per cohort with columns suitable for allele-aware
harmonization and METAL input generation.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path
from statistics import NormalDist


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
GWAS_DIR = PROJECT_DIR / "gwas"
DOWNLOAD_DIR = GWAS_DIR / "download"
FORMATTED_DIR = GWAS_DIR / "formatted"
QC_DIR = GWAS_DIR / "qc"
REPORT_DIR = GWAS_DIR / "reports"
MANIFEST = GWAS_DIR / "config" / "gwas_dataset_manifest.tsv"

OUT_COLUMNS = [
    "SNP",
    "CHR",
    "BP",
    "A1",
    "A2",
    "freq",
    "maf",
    "b",
    "se",
    "p",
    "n",
    "case_n",
    "control_n",
    "INFO",
    "source",
]

ACGT = {"A", "C", "G", "T"}
AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
NORM = NormalDist()


def read_manifest() -> dict[str, dict[str, str]]:
    with MANIFEST.open() as handle:
        return {row["dataset_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def is_autosome(chrom: str) -> bool:
    return chrom.isdigit() and 1 <= int(chrom) <= 22


def is_snv(a1: str, a2: str) -> bool:
    return len(a1) == 1 and len(a2) == 1 and a1 in ACGT and a2 in ACGT and a1 != a2


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def clean_freq(value: str) -> float:
    freq = parse_float(value)
    if not math.isfinite(freq):
        return math.nan
    if freq > 1:
        freq = freq / 100.0
    return freq


def p_to_se(beta: float, pvalue: float) -> float:
    if not math.isfinite(beta) or not math.isfinite(pvalue):
        return math.nan
    if beta == 0 or pvalue <= 0 or pvalue >= 1:
        return math.nan
    z_abs = -NORM.inv_cdf(pvalue / 2.0)
    if not math.isfinite(z_abs) or z_abs <= 0:
        return math.nan
    return abs(beta) / z_abs


def p_from_lp(lp: float) -> float:
    if not math.isfinite(lp) or lp < 0:
        return math.nan
    if lp > 323:
        return 0.0
    return 10.0 ** (-lp)


def init_counts(dataset_id: str) -> dict[str, int | str]:
    return {
        "dataset_id": dataset_id,
        "raw_rows": 0,
        "written_rows": 0,
        "skip_non_autosome": 0,
        "skip_non_snv": 0,
        "skip_ambiguous": 0,
        "skip_bad_freq": 0,
        "skip_maf_lt_0.01": 0,
        "skip_info_lt_0.8": 0,
        "skip_bad_effect": 0,
    }


def write_summary(rows: list[dict[str, int | str]]) -> None:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    path = QC_DIR / "discovery_gwas_standardization_summary.tsv"
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    report = REPORT_DIR / "DISCOVERY_GWAS_STANDARDIZATION_QC.md"
    lines = [
        "# Discovery GWAS Standardization QC",
        "",
        "| Dataset | Raw rows | Written | Non-autosome | Non-SNV | Ambiguous | Bad freq | MAF<0.01 | INFO<0.8 | Bad effect |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['dataset_id']} | {row['raw_rows']} | {row['written_rows']} | "
            f"{row['skip_non_autosome']} | {row['skip_non_snv']} | {row['skip_ambiguous']} | "
            f"{row['skip_bad_freq']} | {row['skip_maf_lt_0.01']} | {row['skip_info_lt_0.8']} | "
            f"{row['skip_bad_effect']} |"
        )
    lines += [
        "",
        "Filters used for the primary METAL input:",
        "",
        "- Autosomes only.",
        "- Biallelic A/C/G/T SNVs only.",
        "- A/T and C/G ambiguous SNPs removed.",
        "- Effect allele frequency must be in [0, 1]; values greater than 1 are divided by 100 before checking.",
        "- MAF >= 0.01.",
        "- INFO >= 0.8 when INFO is available.",
        "- Valid beta, SE, and P required.",
    ]
    report.write_text("\n".join(lines) + "\n")


def standardize_decode(dataset_id: str, manifest: dict[str, str]) -> dict[str, int | str]:
    counts = init_counts(dataset_id)
    in_path = PROJECT_DIR / manifest["local_raw_path"]
    out_path = FORMATTED_DIR / f"{dataset_id}.standardized.tsv.gz"
    case_n = manifest["cases"]
    control_n = manifest["controls"]
    total_n = str(int(case_n) + int(control_n))

    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=OUT_COLUMNS)
        writer.writeheader()
        for row in reader:
            counts["raw_rows"] += 1
            chrom = row["CHR"]
            if not is_autosome(chrom):
                counts["skip_non_autosome"] += 1
                continue
            a1 = row["A1"].upper()
            a2 = row["A0"].upper()
            if not is_snv(a1, a2):
                counts["skip_non_snv"] += 1
                continue
            if (a1, a2) in AMBIGUOUS:
                counts["skip_ambiguous"] += 1
                continue
            freq = clean_freq(row["A1_FREQ"])
            if not math.isfinite(freq) or freq < 0 or freq > 1:
                counts["skip_bad_freq"] += 1
                continue
            maf = min(freq, 1.0 - freq)
            if maf < 0.01:
                counts["skip_maf_lt_0.01"] += 1
                continue
            info = parse_float(row["INFO"])
            if math.isfinite(info) and info < 0.8:
                counts["skip_info_lt_0.8"] += 1
                continue
            beta = parse_float(row["BETA"])
            pvalue = parse_float(row["PVAL"])
            se = p_to_se(beta, pvalue)
            if not all(math.isfinite(x) for x in [beta, pvalue, se]) or pvalue <= 0 or pvalue > 1 or se <= 0:
                counts["skip_bad_effect"] += 1
                continue
            writer.writerow(
                {
                    "SNP": row["RS"] if row["RS"] else f"{chrom}:{row['BP']}:{a2}:{a1}",
                    "CHR": chrom,
                    "BP": row["BP"],
                    "A1": a1,
                    "A2": a2,
                    "freq": f"{freq:.8g}",
                    "maf": f"{maf:.8g}",
                    "b": f"{beta:.8g}",
                    "se": f"{se:.8g}",
                    "p": f"{pvalue:.8g}",
                    "n": total_n,
                    "case_n": case_n,
                    "control_n": control_n,
                    "INFO": f"{info:.8g}" if math.isfinite(info) else "NA",
                    "source": dataset_id,
                }
            )
            counts["written_rows"] += 1
    return counts


def standardize_ukb(dataset_id: str, manifest: dict[str, str]) -> dict[str, int | str]:
    counts = init_counts(dataset_id)
    in_path = PROJECT_DIR / manifest["local_raw_path"]
    out_path = FORMATTED_DIR / f"{dataset_id}.standardized.tsv.gz"
    case_n = manifest["cases"]
    control_n = manifest["controls"]

    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt", newline="") as fout:
        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=OUT_COLUMNS)
        writer.writeheader()
        for line in fin:
            if line.startswith("#"):
                continue
            counts["raw_rows"] += 1
            parts = line.rstrip("\n").split("\t")
            chrom, bp, rsid, ref, alt = parts[:5]
            if not is_autosome(chrom):
                counts["skip_non_autosome"] += 1
                continue
            a1 = alt.upper()
            a2 = ref.upper()
            if not is_snv(a1, a2):
                counts["skip_non_snv"] += 1
                continue
            if (a1, a2) in AMBIGUOUS:
                counts["skip_ambiguous"] += 1
                continue
            fmt = parts[8].split(":")
            vals = parts[9].split(":")
            fdict = dict(zip(fmt, vals))
            freq = clean_freq(fdict.get("AF", ""))
            if not math.isfinite(freq) or freq < 0 or freq > 1:
                counts["skip_bad_freq"] += 1
                continue
            maf = min(freq, 1.0 - freq)
            if maf < 0.01:
                counts["skip_maf_lt_0.01"] += 1
                continue
            beta = parse_float(fdict.get("ES", ""))
            se = parse_float(fdict.get("SE", ""))
            lp = parse_float(fdict.get("LP", ""))
            pvalue = p_from_lp(lp)
            n = fdict.get("SS", "")
            if not all(math.isfinite(x) for x in [beta, se, pvalue]) or se <= 0 or pvalue < 0 or pvalue > 1:
                counts["skip_bad_effect"] += 1
                continue
            writer.writerow(
                {
                    "SNP": rsid if rsid != "." else f"{chrom}:{bp}:{ref}:{alt}",
                    "CHR": chrom,
                    "BP": bp,
                    "A1": a1,
                    "A2": a2,
                    "freq": f"{freq:.8g}",
                    "maf": f"{maf:.8g}",
                    "b": f"{beta:.8g}",
                    "se": f"{se:.8g}",
                    "p": f"{pvalue:.8g}",
                    "n": n,
                    "case_n": case_n,
                    "control_n": control_n,
                    "INFO": "NA",
                    "source": dataset_id,
                }
            )
            counts["written_rows"] += 1
    return counts


def main() -> int:
    FORMATTED_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = read_manifest()
    rows = []
    for dataset_id in ["decode_chb_cvdc_dbds", "decode_iceland", "decode_intermountain"]:
        print(f"Standardizing {dataset_id}", flush=True)
        rows.append(standardize_decode(dataset_id, manifest[dataset_id]))
    print("Standardizing ukb_i9_vte", flush=True)
    rows.append(standardize_ukb("ukb_i9_vte", manifest["ukb_i9_vte"]))
    write_summary(rows)
    print(f"Wrote standardized files to {FORMATTED_DIR}")
    print(f"Wrote QC summary to {QC_DIR / 'discovery_gwas_standardization_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
