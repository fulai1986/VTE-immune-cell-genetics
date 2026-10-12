#!/usr/bin/env python3
"""Prepare HapMap3-matched summary statistics for LDSC QC."""

from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
GWAS_PATH = (
    PROJECT_DIR
    / "gwas"
    / "meta"
    / "vte_discovery_finngen_independent_ivw_ukb_approx_logor.primary_ncohort3plus.tsv.gz"
)
REF_DIR = PROJECT_DIR / "tools" / "ldsc_reference" / "eur_w_ld_chr"
OUT_DIR = PROJECT_DIR / "gwas" / "ldsc"
SUMSTATS = OUT_DIR / "vte_meta_ukb_approx_logor_ncohort3plus.hm3.sumstats.gz"
MATCHED = OUT_DIR / "vte_meta_ukb_approx_logor_ncohort3plus.hm3.ldscore_matched.tsv.gz"
QC = OUT_DIR / "vte_meta_ukb_approx_logor_ncohort3plus.ldsc_prepare_qc.tsv"

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def comp(allele: str) -> str:
    return allele.translate(COMPLEMENT)


def allele_match(gwas_a1: str, gwas_a2: str, ref_a1: str, ref_a2: str) -> tuple[bool, bool, str]:
    """Return matched, flip_effect, status."""
    gwas_a1 = gwas_a1.upper()
    gwas_a2 = gwas_a2.upper()
    ref_a1 = ref_a1.upper()
    ref_a2 = ref_a2.upper()
    if gwas_a1 == ref_a1 and gwas_a2 == ref_a2:
        return True, False, "exact"
    if gwas_a1 == ref_a2 and gwas_a2 == ref_a1:
        return True, True, "swapped"
    if comp(gwas_a1) == ref_a1 and comp(gwas_a2) == ref_a2:
        return True, False, "strand"
    if comp(gwas_a1) == ref_a2 and comp(gwas_a2) == ref_a1:
        return True, True, "strand_swapped"
    return False, False, "allele_mismatch"


def load_hm3_alleles() -> dict[str, tuple[str, str]]:
    hm3 = {}
    with (REF_DIR / "w_hm3.snplist").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            hm3[row["SNP"]] = (row["A1"].upper(), row["A2"].upper())
    return hm3


def load_ld_reference() -> tuple[dict[tuple[str, str], dict[str, str]], dict[str, int]]:
    hm3 = load_hm3_alleles()
    ref = {}
    counts = {
        "ldscore_rows": 0,
        "hm3_rows": len(hm3),
        "ldscore_hm3_rows": 0,
        "duplicate_chr_bp": 0,
    }
    for chrom in range(1, 23):
        path = REF_DIR / f"{chrom}.l2.ldscore.gz"
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                counts["ldscore_rows"] += 1
                snp = row["SNP"]
                if snp not in hm3:
                    continue
                key = (row["CHR"], row["BP"])
                if key in ref:
                    counts["duplicate_chr_bp"] += 1
                    continue
                a1, a2 = hm3[snp]
                ref[key] = {
                    "SNP": snp,
                    "CHR": row["CHR"],
                    "BP": row["BP"],
                    "A1_ref": a1,
                    "A2_ref": a2,
                    "MAF_ref": row["MAF"],
                    "L2": row["L2"],
                }
                counts["ldscore_hm3_rows"] += 1
    return ref, counts


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ref, counts = load_ld_reference()
    counts.update(
        {
            "gwas_rows": 0,
            "matched_rows": 0,
            "no_hm3_position": 0,
            "allele_mismatch": 0,
            "bad_numeric": 0,
            "exact": 0,
            "swapped": 0,
            "strand": 0,
            "strand_swapped": 0,
        }
    )

    sumstats_fields = ["SNP", "A1", "A2", "N", "Z", "P"]
    matched_fields = [
        "SNP",
        "CHR",
        "BP",
        "A1",
        "A2",
        "N",
        "Z",
        "P",
        "BETA",
        "SE",
        "CHI2",
        "L2",
        "MAF_ref",
        "match_status",
        "original_SNP",
        "original_A1",
        "original_A2",
    ]
    seen = set()
    with gzip.open(GWAS_PATH, "rt", newline="") as fin, gzip.open(SUMSTATS, "wt", newline="") as fout, gzip.open(
        MATCHED, "wt", newline=""
    ) as mout:
        reader = csv.DictReader(fin, delimiter="\t")
        sw = csv.DictWriter(fout, delimiter="\t", fieldnames=sumstats_fields)
        mw = csv.DictWriter(mout, delimiter="\t", fieldnames=matched_fields)
        sw.writeheader()
        mw.writeheader()
        for row in reader:
            counts["gwas_rows"] += 1
            key = (row["CHR"], row["BP"])
            ref_row = ref.get(key)
            if ref_row is None:
                counts["no_hm3_position"] += 1
                continue
            matched, flip, status = allele_match(row["A1"], row["A2"], ref_row["A1_ref"], ref_row["A2_ref"])
            if not matched:
                counts["allele_mismatch"] += 1
                continue
            beta = parse_float(row["b"])
            se = parse_float(row["se"])
            p = parse_float(row["p"])
            n = parse_float(row["n"])
            l2 = parse_float(ref_row["L2"])
            if not all(math.isfinite(x) for x in [beta, se, p, n, l2]) or se <= 0 or p <= 0 or p > 1:
                counts["bad_numeric"] += 1
                continue
            if ref_row["SNP"] in seen:
                continue
            seen.add(ref_row["SNP"])
            beta = -beta if flip else beta
            z = beta / se
            out = {
                "SNP": ref_row["SNP"],
                "A1": ref_row["A1_ref"],
                "A2": ref_row["A2_ref"],
                "N": f"{n:.0f}",
                "Z": f"{z:.10g}",
                "P": f"{p:.10g}",
            }
            sw.writerow(out)
            mw.writerow(
                {
                    **out,
                    "CHR": row["CHR"],
                    "BP": row["BP"],
                    "BETA": f"{beta:.10g}",
                    "SE": f"{se:.10g}",
                    "CHI2": f"{z * z:.10g}",
                    "L2": ref_row["L2"],
                    "MAF_ref": ref_row["MAF_ref"],
                    "match_status": status,
                    "original_SNP": row["SNP"],
                    "original_A1": row["A1"],
                    "original_A2": row["A2"],
                }
            )
            counts["matched_rows"] += 1
            counts[status] += 1

    with QC.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in counts.items():
            writer.writerow([key, value])

    print(f"Wrote {SUMSTATS}")
    print(f"Wrote {MATCHED}")
    print(f"Wrote {QC}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
