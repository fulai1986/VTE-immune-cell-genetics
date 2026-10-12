#!/usr/bin/env python3
"""Harmonize lifted FinnGen R13 I9_VTE GWAS to OneK1K LD/eQTL IDs for SMR."""

from __future__ import annotations

import csv
import gzip
import math
from collections import defaultdict
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
GWAS_DIR = PROJECT_DIR / "gwas"
FINNGEN_DIR = GWAS_DIR / "replication" / "finngen_r13"
INPUT = FINNGEN_DIR / "liftover_grch37" / "finngen_R13_I9_VTE.hg19.tsv.gz"
HARMONIZED_DIR = FINNGEN_DIR / "harmonized_oneK1K_ld"
QC_DIR = GWAS_DIR / "qc"
REPORT_DIR = GWAS_DIR / "reports"
LD_DIR = PROJECT_DIR / "results" / "ld_reference" / "onek1k"

PHENOCODE = "I9_VTE"
PREFIX = "finngen_R13_I9_VTE"
SAMPLE_N = 27_928 + 472_258

SMR_COLUMNS = ["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]
ANNOTATED_COLUMNS = [
    "SNP",
    "CHR",
    "BP",
    "A1",
    "A2",
    "freq",
    "b",
    "se",
    "p",
    "n",
    "rsids",
    "nearest_genes",
    "af_alt_cases",
    "af_alt_controls",
    "original_chrom",
    "original_pos",
    "liftover_strand",
    "source_ref",
    "source_alt",
    "source_af_alt",
    "source_beta",
    "harmonization",
]

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def comp(allele: str) -> str:
    return allele.upper().translate(COMPLEMENT)


def parse_ld_id(snp_id: str) -> tuple[str, int, str, str] | None:
    fields = snp_id.split(":")
    if len(fields) != 4:
        return None
    chrom, bp, ref, alt = fields
    try:
        return chrom, int(bp), ref.upper(), alt.upper()
    except ValueError:
        return None


def load_ld_map() -> tuple[dict[tuple[str, int], list[tuple[str, str, str]]], dict[str, int]]:
    ld_map: dict[tuple[str, int], list[tuple[str, str, str]]] = defaultdict(list)
    counts = {
        "ld_variants": 0,
        "ld_bad_ids": 0,
        "ld_duplicate_positions": 0,
    }
    for chrom in range(1, 23):
        bim = LD_DIR / f"chr{chrom}" / f"onek1k.chr{chrom}.r2_0.8.maf_0.01.allele_id.bim"
        with bim.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                parsed = parse_ld_id(parts[1])
                if parsed is None:
                    counts["ld_bad_ids"] += 1
                    continue
                ld_chrom, bp, ref, alt = parsed
                key = (ld_chrom, bp)
                if ld_map[key]:
                    counts["ld_duplicate_positions"] += 1
                ld_map[key].append((ref, alt, parts[1]))
                counts["ld_variants"] += 1
    return ld_map, counts


def parse_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return math.nan


def fmt(value: float) -> str:
    return f"{value:.10g}"


def harmonize_row(row: dict[str, str], ld_records: list[tuple[str, str, str]]) -> tuple[dict[str, str], str] | tuple[None, str]:
    source_ref = row["ref"].upper()
    source_alt = row["alt"].upper()
    beta = parse_float(row["beta"])
    se = parse_float(row["sebeta"])
    p = parse_float(row["pval"])
    af_alt = parse_float(row["af_alt"])

    if not all(math.isfinite(x) for x in [beta, se, p, af_alt]) or p <= 0 or p > 1 or se <= 0:
        return None, "bad_effect"
    if af_alt < -1e-6 or af_alt > 1 + 1e-6:
        return None, "bad_freq"

    candidates: list[tuple[str, str, str, bool, str]] = []
    for ld_ref, ld_alt, snp_id in ld_records:
        if source_ref == ld_ref and source_alt == ld_alt:
            candidates.append((ld_ref, ld_alt, snp_id, False, "exact_alt_effect"))
        elif source_ref == ld_alt and source_alt == ld_ref:
            candidates.append((ld_ref, ld_alt, snp_id, True, "exact_ref_effect_flipped"))
        elif comp(source_ref) == ld_ref and comp(source_alt) == ld_alt:
            candidates.append((ld_ref, ld_alt, snp_id, False, "strand_flip_alt_effect"))
        elif comp(source_ref) == ld_alt and comp(source_alt) == ld_ref:
            candidates.append((ld_ref, ld_alt, snp_id, True, "strand_flip_ref_effect_flipped"))

    if not candidates:
        return None, "allele_mismatch"
    if len(candidates) > 1:
        return None, "multiple_ld_allele_matches"

    ld_ref, ld_alt, snp_id, flip, reason = candidates[0]
    out_beta = -beta if flip else beta
    out_freq = 1.0 - af_alt if flip else af_alt
    out_freq = min(max(out_freq, 0.0), 1.0)

    out = {
        "SNP": snp_id,
        "CHR": row["chrom"],
        "BP": row["pos"],
        "A1": ld_alt,
        "A2": ld_ref,
        "freq": fmt(out_freq),
        "b": fmt(out_beta),
        "se": row["sebeta"],
        "p": row["pval"],
        "n": str(SAMPLE_N),
        "rsids": row["rsids"],
        "nearest_genes": row["nearest_genes"],
        "af_alt_cases": row["af_alt_cases"],
        "af_alt_controls": row["af_alt_controls"],
        "original_chrom": row["original_chrom"],
        "original_pos": row["original_pos"],
        "liftover_strand": row["liftover_strand"],
        "source_ref": source_ref,
        "source_alt": source_alt,
        "source_af_alt": row["af_alt"],
        "source_beta": row["beta"],
        "harmonization": reason,
    }
    return out, reason


def main() -> int:
    HARMONIZED_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    out_tsv = HARMONIZED_DIR / f"{PREFIX}.hg19.onek1k_ld_harmonized.tsv.gz"
    out_ma = HARMONIZED_DIR / f"{PREFIX}.hg19.onek1k_ld_harmonized.smr.ma"
    rejected_tsv = HARMONIZED_DIR / f"{PREFIX}.hg19.onek1k_ld_rejected.tsv.gz"

    ld_map, ld_counts = load_ld_map()
    seen_snps: set[str] = set()
    counts = {
        "input_rows": 0,
        "harmonized_rows": 0,
        "non_autosome": 0,
        "non_snv": 0,
        "bad_position": 0,
        "no_ld_position": 0,
        "allele_mismatch": 0,
        "multiple_ld_allele_matches": 0,
        "duplicate_harmonized_snp": 0,
        "bad_effect": 0,
        "bad_freq": 0,
        "exact_alt_effect": 0,
        "exact_ref_effect_flipped": 0,
        "strand_flip_alt_effect": 0,
        "strand_flip_ref_effect_flipped": 0,
    }

    with gzip.open(INPUT, "rt", newline="") as fin, gzip.open(out_tsv, "wt", newline="") as fout, out_ma.open(
        "w", newline=""
    ) as maout, gzip.open(rejected_tsv, "wt", newline="") as rejout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"No header found in {INPUT}")
        rejected_fields = list(reader.fieldnames) + ["reject_reason"]
        annotated_writer = csv.DictWriter(fout, delimiter="\t", fieldnames=ANNOTATED_COLUMNS)
        ma_writer = csv.DictWriter(maout, delimiter="\t", fieldnames=SMR_COLUMNS)
        rejected_writer = csv.DictWriter(rejout, delimiter="\t", fieldnames=rejected_fields)
        annotated_writer.writeheader()
        ma_writer.writeheader()
        rejected_writer.writeheader()

        for row in reader:
            counts["input_rows"] += 1
            chrom = row["chrom"].replace("chr", "")
            row["chrom"] = chrom
            if chrom not in {str(i) for i in range(1, 23)}:
                counts["non_autosome"] += 1
                row["reject_reason"] = "non_autosome"
                rejected_writer.writerow(row)
                continue
            if len(row["ref"]) != 1 or len(row["alt"]) != 1 or row["ref"].upper() not in "ACGT" or row["alt"].upper() not in "ACGT":
                counts["non_snv"] += 1
                row["reject_reason"] = "non_snv"
                rejected_writer.writerow(row)
                continue
            try:
                pos = int(row["pos"])
            except ValueError:
                counts["bad_position"] += 1
                row["reject_reason"] = "bad_position"
                rejected_writer.writerow(row)
                continue

            ld_records = ld_map.get((chrom, pos))
            if not ld_records:
                counts["no_ld_position"] += 1
                row["reject_reason"] = "no_ld_position"
                rejected_writer.writerow(row)
                continue

            harmonized, reason = harmonize_row(row, ld_records)
            if harmonized is None:
                counts[reason] += 1
                row["reject_reason"] = reason
                rejected_writer.writerow(row)
                continue
            if harmonized["SNP"] in seen_snps:
                counts["duplicate_harmonized_snp"] += 1
                row["reject_reason"] = "duplicate_harmonized_snp"
                rejected_writer.writerow(row)
                continue

            seen_snps.add(harmonized["SNP"])
            counts["harmonized_rows"] += 1
            counts[reason] += 1
            annotated_writer.writerow(harmonized)
            ma_writer.writerow({col: harmonized[col] for col in SMR_COLUMNS})

            if counts["input_rows"] % 1_000_000 == 0:
                print(
                    f"processed={counts['input_rows']} harmonized={counts['harmonized_rows']} "
                    f"no_ld={counts['no_ld_position']} mismatch={counts['allele_mismatch']}",
                    flush=True,
                )

    qc_path = QC_DIR / f"{PREFIX}.hg19_onek1k_ld_harmonization_qc.tsv"
    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in counts.items():
            writer.writerow([key, value])
        for key, value in ld_counts.items():
            writer.writerow([key, value])
        writer.writerow(["phenocode", PHENOCODE])
        writer.writerow(["sample_n", SAMPLE_N])
        writer.writerow(["input", INPUT.relative_to(PROJECT_DIR)])
        writer.writerow(["output_tsv_gz", out_tsv.relative_to(PROJECT_DIR)])
        writer.writerow(["output_smr_ma", out_ma.relative_to(PROJECT_DIR)])
        writer.writerow(["rejected_tsv_gz", rejected_tsv.relative_to(PROJECT_DIR)])

    report = REPORT_DIR / "FINNGEN_R13_I9_VTE_ONEK1K_LD_HARMONIZATION_QC.md"
    lines = [
        "# FinnGen R13 I9_VTE to OneK1K LD Harmonization QC",
        "",
        "Final FinnGen replication SMR GWAS output is harmonized to `SNP=chr:pos:REF:ALT`, `A1=ALT`, `A2=REF`, `freq=ALT frequency`, and `b=ALT effect`.",
        "",
        "## Inputs",
        "",
        f"- Lifted FinnGen GRCh37 input: `{INPUT.relative_to(PROJECT_DIR)}`",
        f"- OneK1K LD variants loaded: {ld_counts['ld_variants']:,}",
        f"- SMR sample size `n`: {SAMPLE_N:,}",
        "",
        "## Counts",
        "",
        f"- Input rows: {counts['input_rows']:,}",
        f"- Harmonized rows: {counts['harmonized_rows']:,}",
        f"- Non-autosome rows: {counts['non_autosome']:,}",
        f"- Non-SNV rows: {counts['non_snv']:,}",
        f"- No LD position: {counts['no_ld_position']:,}",
        f"- Allele mismatch: {counts['allele_mismatch']:,}",
        f"- Multiple LD allele matches: {counts['multiple_ld_allele_matches']:,}",
        f"- Duplicate harmonized SNP rows removed: {counts['duplicate_harmonized_snp']:,}",
        f"- Bad effect rows: {counts['bad_effect']:,}",
        f"- Bad frequency rows: {counts['bad_freq']:,}",
        f"- Exact ALT-effect rows: {counts['exact_alt_effect']:,}",
        f"- Exact REF-effect rows flipped: {counts['exact_ref_effect_flipped']:,}",
        f"- Strand-flip ALT-effect rows: {counts['strand_flip_alt_effect']:,}",
        f"- Strand-flip REF-effect rows flipped: {counts['strand_flip_ref_effect_flipped']:,}",
        "",
        "## Outputs",
        "",
        f"- Annotated TSV: `{out_tsv.relative_to(PROJECT_DIR)}`",
        f"- SMR `.ma`: `{out_ma.relative_to(PROJECT_DIR)}`",
        f"- Rejected rows: `{rejected_tsv.relative_to(PROJECT_DIR)}`",
        f"- QC table: `{qc_path.relative_to(PROJECT_DIR)}`",
    ]
    report.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_tsv}")
    print(f"Wrote {out_ma}")
    print(f"Wrote {rejected_tsv}")
    print(f"Wrote {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
