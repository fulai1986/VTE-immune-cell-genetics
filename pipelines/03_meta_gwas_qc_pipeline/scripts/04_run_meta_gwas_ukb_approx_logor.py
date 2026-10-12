#!/usr/bin/env python3
"""Run standard VTE meta-GWAS using UKB approximate log(OR) effects.

This pipeline preserves the original meta-GWAS outputs and creates a new
standard result set with UKB `ukb-d-I9_VTE` converted from Z scores to an
approximate log(OR) scale before METAL inverse-variance meta-analysis.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import subprocess
from collections import defaultdict
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
GWAS_DIR = PROJECT_DIR / "gwas"
FORMATTED_DIR = GWAS_DIR / "formatted"
META_DIR = GWAS_DIR / "meta"
METAL_DIR = META_DIR / "metal"
INPUT_DIR = META_DIR / "input_ukb_approx_logor"
INPUT_UNZIPPED_DIR = META_DIR / "input_unzipped_ukb_approx_logor"
QC_DIR = GWAS_DIR / "qc"
REPORT_DIR = GWAS_DIR / "reports"
HARMONIZED_DIR = GWAS_DIR / "harmonized" / "oneK1K_ld"
LD_DIR = PROJECT_DIR / "results" / "ld_reference" / "onek1k"
METAL_BIN = Path(os.environ.get("METAL_BIN", PROJECT_DIR / "tools" / "metal" / "METAL" / "build" / "metal" / "metal"))

PREFIX = "vte_discovery_finngen_independent_ivw_ukb_approx_logor"
OLD_PREFIX = "vte_discovery_finngen_independent_ivw"

DATASET_ORDER = [
    ("decode_chb_cvdc_dbds", 232072),
    ("decode_iceland", 323445),
    ("decode_intermountain", 58233),
    ("ukb_i9_vte_approx_logor", 361194),
]

DATASET_INPUTS = {
    "decode_chb_cvdc_dbds": META_DIR / "input" / "decode_chb_cvdc_dbds.metal.txt.gz",
    "decode_iceland": META_DIR / "input" / "decode_iceland.metal.txt.gz",
    "decode_intermountain": META_DIR / "input" / "decode_intermountain.metal.txt.gz",
    "ukb_i9_vte_approx_logor": INPUT_DIR / "ukb_i9_vte_approx_logor.metal.txt.gz",
}

METAL_COLUMNS = ["MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value", "Freq1", "Weight"]
FULL_COLUMNS = [
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
    "n_cohorts",
    "direction",
    "het_i2",
    "het_chisq",
    "het_df",
    "het_p",
    "metal_marker",
]
SMR_COLUMNS = ["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]
ANNOTATED_COLUMNS = FULL_COLUMNS + [
    "original_snp",
    "original_a1",
    "original_a2",
    "original_freq",
    "original_b",
    "harmonization",
]
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def marker_name(row: dict[str, str]) -> str:
    alleles = sorted([row["A1"].upper(), row["A2"].upper()])
    return f"{row['CHR']}:{row['BP']}:{alleles[0]}:{alleles[1]}"


def prepare_ukb_metal_input() -> dict[str, int | str]:
    """Deduplicate converted UKB approximate-logOR rows into METAL format."""
    source = FORMATTED_DIR / "ukb_i9_vte.approx_logor.standardized.tsv.gz"
    out_path = DATASET_INPUTS["ukb_i9_vte_approx_logor"]
    best: dict[str, dict[str, str]] = {}
    raw_rows = 0
    duplicate_marker_rows = 0

    with gzip.open(source, "rt", newline="") as handle:
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
        for marker in sorted(best, key=lambda x: (int(x.split(":")[0]), int(x.split(":")[1]), x)):
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
        "dataset_id": "ukb_i9_vte_approx_logor",
        "standardized_rows": raw_rows,
        "metal_rows": len(best),
        "duplicate_marker_rows_removed": duplicate_marker_rows,
        "metal_input": str(out_path.relative_to(PROJECT_DIR)),
    }


def write_input_qc(ukb_row: dict[str, int | str]) -> None:
    summary_path = QC_DIR / f"{PREFIX}.metal_input_prepare_summary.tsv"
    existing_summary = QC_DIR / "metal_input_prepare_summary.tsv"
    rows: list[dict[str, str]] = []
    if existing_summary.exists():
        with existing_summary.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row["dataset_id"].startswith("decode_"):
                    rows.append(row)
    rows.append({k: str(v) for k, v in ukb_row.items()})

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def prepare_unzipped_inputs() -> None:
    INPUT_UNZIPPED_DIR.mkdir(parents=True, exist_ok=True)
    for dataset, _ in DATASET_ORDER:
        source = DATASET_INPUTS[dataset]
        target = INPUT_UNZIPPED_DIR / f"{dataset}.metal.txt"
        with gzip.open(source, "rt") as fin, target.open("w") as fout:
            for chunk in iter(lambda: fin.read(1024 * 1024), ""):
                if not chunk:
                    break
                fout.write(chunk)


def write_metal_protocol() -> Path:
    protocol = METAL_DIR / f"{PREFIX}.metal"
    lines = [
        "# METAL protocol for FinnGen-independent VTE discovery meta-GWAS",
        "# UKB input uses approximate log(OR) converted from Z scores.",
        "",
        "SCHEME STDERR",
        "AVERAGEFREQ ON",
        "MINMAXFREQ ON",
        "GENOMICCONTROL OFF",
        "",
        "MARKER MarkerName",
        "ALLELE Allele1 Allele2",
        "EFFECT Effect",
        "STDERR StdErr",
        "PVALUE P-value",
        "FREQ Freq1",
        "WEIGHT Weight",
        "",
    ]
    for dataset, _ in DATASET_ORDER:
        lines.append(f"PROCESS {INPUT_UNZIPPED_DIR / f'{dataset}.metal.txt'}")
    lines += [
        "",
        f"OUTFILE {METAL_DIR / f'{PREFIX}.'} .tbl",
        "ANALYZE HETEROGENEITY",
        "QUIT",
        "",
    ]
    protocol.write_text("\n".join(lines))
    return protocol


def run_metal(protocol: Path) -> Path:
    log_path = GWAS_DIR / "logs" / f"{PREFIX}.metal.log"
    with log_path.open("w") as log:
        subprocess.run([str(METAL_BIN), str(protocol)], cwd=PROJECT_DIR, check=True, stdout=log, stderr=subprocess.STDOUT)
    return log_path


def parse_marker(marker: str) -> tuple[str, str, str, str]:
    chrom, bp, allele_low, allele_high = marker.split(":")
    return chrom, bp, allele_low.upper(), allele_high.upper()


def cohort_count(direction: str) -> int:
    return sum(1 for c in direction if c != "?")


def cohort_n(direction: str) -> int:
    return sum(n for (_, n), flag in zip(DATASET_ORDER, direction) if flag != "?")


def sortable_p(row: dict[str, str]) -> float:
    value = as_float(row["p"])
    return value if math.isfinite(value) else math.inf


def postprocess_metal() -> dict[str, int]:
    metal_tbl = METAL_DIR / f"{PREFIX}.1.tbl"
    full_gz = META_DIR / f"{PREFIX}.metal_full.tsv.gz"
    primary_gz = META_DIR / f"{PREFIX}.primary_ncohort2.tsv.gz"
    smr_ma = META_DIR / f"{PREFIX}.primary_ncohort2.smr.ma"
    top_tsv = QC_DIR / f"{PREFIX}.top_hits.tsv"
    qc_tsv = QC_DIR / f"{PREFIX}.postprocess_qc.tsv"
    report = REPORT_DIR / f"{PREFIX}.METAL_POSTPROCESS_QC.md"

    counts = {
        "metal_rows": 0,
        "primary_ncohort2_rows": 0,
        "ncohort1_rows": 0,
        "ncohort2_rows": 0,
        "ncohort3_rows": 0,
        "ncohort4_rows": 0,
        "genomewide_p5e-8_primary": 0,
        "heidi_usable_primary_p5e-8_hetp_gt0.01": 0,
    }
    top_rows: list[dict[str, str]] = []

    with metal_tbl.open() as fin, gzip.open(full_gz, "wt", newline="") as full_out, gzip.open(
        primary_gz, "wt", newline=""
    ) as primary_out, smr_ma.open("w", newline="") as smr_out:
        reader = csv.DictReader(fin, delimiter="\t")
        full_writer = csv.DictWriter(full_out, delimiter="\t", fieldnames=FULL_COLUMNS)
        primary_writer = csv.DictWriter(primary_out, delimiter="\t", fieldnames=FULL_COLUMNS)
        smr_writer = csv.DictWriter(smr_out, delimiter="\t", fieldnames=SMR_COLUMNS)
        full_writer.writeheader()
        primary_writer.writeheader()
        smr_writer.writeheader()

        for row in reader:
            counts["metal_rows"] += 1
            chrom, bp, _allele_low, _allele_high = parse_marker(row["MarkerName"])
            a1 = row["Allele1"].upper()
            a2 = row["Allele2"].upper()
            direction = row["Direction"]
            n_cohorts = cohort_count(direction)
            n = cohort_n(direction)
            p = as_float(row["P-value"])
            het_p = as_float(row["HetPVal"])
            out = {
                "SNP": f"{chrom}:{bp}:{a2}:{a1}",
                "CHR": chrom,
                "BP": bp,
                "A1": a1,
                "A2": a2,
                "freq": row["Freq1"],
                "b": row["Effect"],
                "se": row["StdErr"],
                "p": row["P-value"],
                "n": str(n),
                "n_cohorts": str(n_cohorts),
                "direction": direction,
                "het_i2": row["HetISq"],
                "het_chisq": row["HetChiSq"],
                "het_df": row["HetDf"],
                "het_p": row["HetPVal"],
                "metal_marker": row["MarkerName"],
            }
            full_writer.writerow(out)
            counts[f"ncohort{n_cohorts}_rows"] += 1
            if n_cohorts >= 2:
                counts["primary_ncohort2_rows"] += 1
                primary_writer.writerow(out)
                smr_writer.writerow({k: out[k] for k in SMR_COLUMNS})
                if math.isfinite(p) and p < 5e-8:
                    counts["genomewide_p5e-8_primary"] += 1
                    if math.isfinite(het_p) and het_p > 0.01:
                        counts["heidi_usable_primary_p5e-8_hetp_gt0.01"] += 1
                top_rows.append(out)
                if len(top_rows) > 10000:
                    top_rows = sorted(top_rows, key=sortable_p)[:1000]

    top_rows = sorted(top_rows, key=sortable_p)[:1000]
    with top_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FULL_COLUMNS)
        writer.writeheader()
        writer.writerows(top_rows)

    with qc_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in counts.items():
            writer.writerow([key, value])

    lines = [
        "# VTE Discovery METAL Postprocess QC: UKB approximate log(OR)",
        "",
        f"- METAL rows: {counts['metal_rows']:,}",
        f"- Primary rows with at least 2 cohorts: {counts['primary_ncohort2_rows']:,}",
        f"- Rows present in 1 cohort: {counts['ncohort1_rows']:,}",
        f"- Rows present in 2 cohorts: {counts['ncohort2_rows']:,}",
        f"- Rows present in 3 cohorts: {counts['ncohort3_rows']:,}",
        f"- Rows present in 4 cohorts: {counts['ncohort4_rows']:,}",
        f"- Primary rows with P < 5e-8: {counts['genomewide_p5e-8_primary']:,}",
        f"- P < 5e-8 and heterogeneity P > 0.01: {counts['heidi_usable_primary_p5e-8_hetp_gt0.01']:,}",
        "",
        "## Outputs",
        "",
        f"- Full METAL result gzip: `{full_gz.relative_to(PROJECT_DIR)}`",
        f"- Primary n_cohorts>=2 result: `{primary_gz.relative_to(PROJECT_DIR)}`",
        f"- SMR-ready GWAS summary: `{smr_ma.relative_to(PROJECT_DIR)}`",
        f"- Top hits: `{top_tsv.relative_to(PROJECT_DIR)}`",
        f"- QC table: `{qc_tsv.relative_to(PROJECT_DIR)}`",
        "",
        "## Top 10 Primary Signals",
        "",
        "| SNP | A1 | A2 | beta | se | P | n_cohorts | direction | HetP |",
        "|---|---|---|---:|---:|---:|---:|---|---:|",
    ]
    for row in top_rows[:10]:
        lines.append(
            f"| {row['SNP']} | {row['A1']} | {row['A2']} | {row['b']} | {row['se']} | "
            f"{row['p']} | {row['n_cohorts']} | {row['direction']} | {row['het_p']} |"
        )
    report.write_text("\n".join(lines) + "\n")
    return counts


def passes_filter(n_cohorts: int, stratum: str) -> bool:
    if stratum == "ncohort3plus":
        return n_cohorts >= 3
    if stratum == "ncohort4":
        return n_cohorts == 4
    raise ValueError(f"Unknown stratum: {stratum}")


def make_ncohort_strata() -> None:
    full_gz = META_DIR / f"{PREFIX}.metal_full.tsv.gz"
    strata = {
        "ncohort3plus": {
            "label": "n_cohorts>=3",
            "tsv": META_DIR / f"{PREFIX}.primary_ncohort3plus.tsv.gz",
            "ma": META_DIR / f"{PREFIX}.primary_ncohort3plus.smr.ma",
            "rows": 0,
            "p5e8": 0,
            "p5e8_hetp_gt0_01": 0,
        },
        "ncohort4": {
            "label": "n_cohorts=4",
            "tsv": META_DIR / f"{PREFIX}.primary_ncohort4.tsv.gz",
            "ma": META_DIR / f"{PREFIX}.primary_ncohort4.smr.ma",
            "rows": 0,
            "p5e8": 0,
            "p5e8_hetp_gt0_01": 0,
        },
    }
    writers = {}
    try:
        for key, spec in strata.items():
            tsv_handle = gzip.open(spec["tsv"], "wt", newline="")
            ma_handle = spec["ma"].open("w", newline="")
            tsv_writer = csv.DictWriter(tsv_handle, delimiter="\t", fieldnames=FULL_COLUMNS)
            ma_writer = csv.DictWriter(ma_handle, delimiter="\t", fieldnames=SMR_COLUMNS)
            tsv_writer.writeheader()
            ma_writer.writeheader()
            writers[key] = (tsv_handle, ma_handle, tsv_writer, ma_writer)

        with gzip.open(full_gz, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                n_cohorts = int(row["n_cohorts"])
                p = as_float(row["p"])
                het_p = as_float(row["het_p"])
                for key, spec in strata.items():
                    if not passes_filter(n_cohorts, key):
                        continue
                    _, _, tsv_writer, ma_writer = writers[key]
                    tsv_writer.writerow(row)
                    ma_writer.writerow({col: row[col] for col in SMR_COLUMNS})
                    spec["rows"] += 1
                    if math.isfinite(p) and p < 5e-8:
                        spec["p5e8"] += 1
                        if math.isfinite(het_p) and het_p > 0.01:
                            spec["p5e8_hetp_gt0_01"] += 1
    finally:
        for tsv_handle, ma_handle, _, _ in writers.values():
            tsv_handle.close()
            ma_handle.close()

    qc_path = QC_DIR / f"{PREFIX}.ncohort_strata_qc.tsv"
    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["stratum", "filter", "rows", "p_lt_5e-8", "p_lt_5e-8_hetp_gt_0.01", "tsv_gz", "smr_ma"])
        for key, spec in strata.items():
            writer.writerow(
                [
                    key,
                    spec["label"],
                    spec["rows"],
                    spec["p5e8"],
                    spec["p5e8_hetp_gt0_01"],
                    spec["tsv"].relative_to(PROJECT_DIR),
                    spec["ma"].relative_to(PROJECT_DIR),
                ]
            )


def comp(allele: str) -> str:
    return allele.translate(COMPLEMENT)


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
    counts = {"ld_variants": 0, "ld_bad_ids": 0, "ld_duplicate_positions": 0}
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


def format_float(value: float) -> str:
    return f"{value:.10g}"


def harmonize_row(row: dict[str, str], ld_records: list[tuple[str, str, str]]) -> tuple[dict[str, str] | None, str]:
    original_a1 = row["A1"].upper()
    original_a2 = row["A2"].upper()
    original_freq = as_float(row["freq"])
    original_b = as_float(row["b"])
    if not math.isfinite(original_freq) or not math.isfinite(original_b):
        return None, "bad_effect_or_freq"

    candidates: list[tuple[str, str, str, bool, str]] = []
    for ref, alt, snp_id in ld_records:
        if original_a1 == alt and original_a2 == ref:
            candidates.append((ref, alt, snp_id, False, "exact_alt_effect"))
        elif original_a1 == ref and original_a2 == alt:
            candidates.append((ref, alt, snp_id, True, "exact_ref_effect_flipped"))
        elif comp(original_a1) == alt and comp(original_a2) == ref:
            candidates.append((ref, alt, snp_id, False, "strand_flip_alt_effect"))
        elif comp(original_a1) == ref and comp(original_a2) == alt:
            candidates.append((ref, alt, snp_id, True, "strand_flip_ref_effect_flipped"))

    if len(candidates) == 0:
        return None, "allele_mismatch"
    if len(candidates) > 1:
        return None, "multiple_ld_allele_matches"

    ref, alt, snp_id, flip, match_type = candidates[0]
    b = -original_b if flip else original_b
    freq = 1.0 - original_freq if flip else original_freq
    if freq < -1e-6 or freq > 1 + 1e-6:
        return None, "bad_harmonized_freq"
    out = dict(row)
    out.update(
        {
            "SNP": snp_id,
            "A1": alt,
            "A2": ref,
            "freq": format_float(min(max(freq, 0.0), 1.0)),
            "b": format_float(b),
            "original_snp": row["SNP"],
            "original_a1": original_a1,
            "original_a2": original_a2,
            "original_freq": row["freq"],
            "original_b": row["b"],
            "harmonization": match_type,
        }
    )
    return out, match_type


def harmonize_one(label: str, input_path: Path, ld_map: dict[tuple[str, int], list[tuple[str, str, str]]]) -> dict[str, int | str]:
    out_tsv = HARMONIZED_DIR / f"{PREFIX}.primary_{label}.onek1k_ld_harmonized.tsv.gz"
    out_ma = HARMONIZED_DIR / f"{PREFIX}.primary_{label}.onek1k_ld_harmonized.smr.ma"
    rejected_tsv = HARMONIZED_DIR / f"{PREFIX}.primary_{label}.onek1k_ld_rejected.tsv.gz"
    counts: dict[str, int | str] = {
        "dataset": label,
        "input": str(input_path.relative_to(PROJECT_DIR)),
        "input_rows": 0,
        "harmonized_rows": 0,
        "no_ld_position": 0,
        "allele_mismatch": 0,
        "multiple_ld_allele_matches": 0,
        "bad_effect_or_freq": 0,
        "bad_harmonized_freq": 0,
        "exact_alt_effect": 0,
        "exact_ref_effect_flipped": 0,
        "strand_flip_alt_effect": 0,
        "strand_flip_ref_effect_flipped": 0,
        "output_tsv_gz": str(out_tsv.relative_to(PROJECT_DIR)),
        "output_smr_ma": str(out_ma.relative_to(PROJECT_DIR)),
        "rejected_tsv_gz": str(rejected_tsv.relative_to(PROJECT_DIR)),
    }
    with gzip.open(input_path, "rt", newline="") as fin, gzip.open(out_tsv, "wt", newline="") as fout, out_ma.open(
        "w", newline=""
    ) as maout, gzip.open(rejected_tsv, "wt", newline="") as rejout:
        reader = csv.DictReader(fin, delimiter="\t")
        annotated_writer = csv.DictWriter(fout, delimiter="\t", fieldnames=ANNOTATED_COLUMNS)
        ma_writer = csv.DictWriter(maout, delimiter="\t", fieldnames=SMR_COLUMNS)
        rejected_writer = csv.DictWriter(rejout, delimiter="\t", fieldnames=list(reader.fieldnames or []) + ["reject_reason"])
        annotated_writer.writeheader()
        ma_writer.writeheader()
        rejected_writer.writeheader()
        for row in reader:
            counts["input_rows"] = int(counts["input_rows"]) + 1
            try:
                key = (row["CHR"], int(row["BP"]))
            except ValueError:
                row["reject_reason"] = "bad_position"
                rejected_writer.writerow(row)
                continue
            ld_records = ld_map.get(key)
            if not ld_records:
                counts["no_ld_position"] = int(counts["no_ld_position"]) + 1
                row["reject_reason"] = "no_ld_position"
                rejected_writer.writerow(row)
                continue
            harmonized, reason = harmonize_row(row, ld_records)
            if harmonized is None:
                counts[reason] = int(counts.get(reason, 0)) + 1
                row["reject_reason"] = reason
                rejected_writer.writerow(row)
                continue
            counts["harmonized_rows"] = int(counts["harmonized_rows"]) + 1
            counts[reason] = int(counts[reason]) + 1
            annotated_writer.writerow({col: harmonized[col] for col in ANNOTATED_COLUMNS})
            ma_writer.writerow({col: harmonized[col] for col in SMR_COLUMNS})
    return counts


def harmonize_to_onek1k_ld() -> None:
    ld_map, ld_counts = load_ld_map()
    inputs = [
        ("ncohort3plus", META_DIR / f"{PREFIX}.primary_ncohort3plus.tsv.gz"),
        ("ncohort4", META_DIR / f"{PREFIX}.primary_ncohort4.tsv.gz"),
    ]
    dataset_counts = [harmonize_one(label, input_path, ld_map) for label, input_path in inputs]
    qc_path = QC_DIR / f"{PREFIX}.onek1k_ld_harmonization_qc.tsv"
    with qc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["metric", "value"])
        writer.writeheader()
        for key, value in ld_counts.items():
            writer.writerow({"metric": key, "value": value})
    dataset_qc_path = QC_DIR / f"{PREFIX}.onek1k_ld_harmonization_dataset_qc.tsv"
    with dataset_qc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(dataset_counts[0].keys()))
        writer.writeheader()
        writer.writerows(dataset_counts)


def write_run_summary(log_path: Path) -> None:
    report = REPORT_DIR / f"{PREFIX}.RUN_SUMMARY.md"
    lines = [
        "# VTE meta-GWAS standard run: UKB approximate log(OR)",
        "",
        "## Rationale",
        "",
        "The original OpenGWAS `ukb-d-I9_VTE` effect sizes were not on a log(OR)-comparable scale. In this standard workflow, UKB effects are first converted from signed Z scores to approximate log(OR), then combined with the deCODE cohorts by METAL inverse-variance meta-analysis.",
        "",
        "## Key outputs",
        "",
        f"- METAL protocol: `gwas/meta/metal/{PREFIX}.metal`",
        f"- METAL log: `{log_path.relative_to(PROJECT_DIR)}`",
        f"- Full METAL table: `gwas/meta/{PREFIX}.metal_full.tsv.gz`",
        f"- n_cohorts>=3 table: `gwas/meta/{PREFIX}.primary_ncohort3plus.tsv.gz`",
        f"- n_cohorts>=3 SMR `.ma`: `gwas/meta/{PREFIX}.primary_ncohort3plus.smr.ma`",
        f"- n_cohorts>=3 OneK1K-LD harmonized SMR `.ma`: `gwas/harmonized/oneK1K_ld/{PREFIX}.primary_ncohort3plus.onek1k_ld_harmonized.smr.ma`",
        f"- n_cohorts=4 table: `gwas/meta/{PREFIX}.primary_ncohort4.tsv.gz`",
        f"- QC tables prefix: `gwas/qc/{PREFIX}.*.tsv`",
        "",
        "## Caveat",
        "",
        "UKB remains an approximate log(OR) conversion rather than native logistic-regression summary statistics; this should be described in Methods.",
    ]
    report.write_text("\n".join(lines) + "\n")


def main() -> int:
    for path in [INPUT_DIR, INPUT_UNZIPPED_DIR, METAL_DIR, META_DIR, QC_DIR, REPORT_DIR, HARMONIZED_DIR, GWAS_DIR / "logs"]:
        path.mkdir(parents=True, exist_ok=True)
    if not METAL_BIN.exists():
        raise FileNotFoundError(f"METAL binary not found: {METAL_BIN}")
    ukb_row = prepare_ukb_metal_input()
    write_input_qc(ukb_row)
    prepare_unzipped_inputs()
    protocol = write_metal_protocol()
    log_path = run_metal(protocol)
    postprocess_metal()
    make_ncohort_strata()
    harmonize_to_onek1k_ld()
    write_run_summary(log_path)
    print(f"Completed {PREFIX}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
