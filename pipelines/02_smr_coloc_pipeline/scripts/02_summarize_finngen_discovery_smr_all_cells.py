#!/usr/bin/env python3
"""Summarize FinnGen discovery SMR across all OneK1K Science14 cell types."""

from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
CELL_FILE = PROJECT_DIR / "config" / "science14_celltypes.tsv"
FINNGEN_ROOT = PROJECT_DIR / "results" / "smr_analysis" / "finngen_r13_discovery_onek1k_params"
META_ROOT = PROJECT_DIR / "results" / "smr_analysis" / "discovery_ncohort3plus_onek1k_params"
QC_DIR = PROJECT_DIR / "results" / "smr_analysis" / "qc"
REPORT_DIR = PROJECT_DIR / "reports" / "smr_analysis"


def cell_to_id(cell: str) -> str:
    return cell.replace(" ", "_")


def parse_float(value: str) -> float:
    if value in {"", "NA", None}:
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def fmt(value: float) -> str:
    if math.isnan(value):
        return "NA"
    return f"{value:.6g}"


def bh_fdr(p_values: list[float]) -> list[float]:
    indexed = sorted((p, i) for i, p in enumerate(p_values) if not math.isnan(p))
    q_values = [math.nan] * len(p_values)
    previous = 1.0
    m = len(indexed)
    for reverse_rank, (p, i) in enumerate(reversed(indexed), start=1):
        rank = m - reverse_rank + 1
        q = min(previous, p * m / rank)
        previous = q
        q_values[i] = min(q, 1.0)
    return q_values


def read_cells() -> list[str]:
    cells: list[str] = []
    with CELL_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cells.append(cell_to_id(row["science14_cell_type"]))
    return cells


def read_smr_file(path: Path, cell: str, chrom: int) -> tuple[list[str], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            rows.append({"cell": cell, "chr": str(chrom), **row})
    return ["cell", "chr"] + fieldnames, rows


def read_meta_lookup(cell: str) -> dict[tuple[str, str], dict[str, str]]:
    path = META_ROOT / cell / f"{cell}.discovery_ncohort3plus_onek1k_params.allchr.smr.tsv.gz"
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    if not path.exists():
        return lookup
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            lookup[(row["Gene"], row["ProbeChr"])] = row
    return lookup


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    path.write_bytes(path.read_bytes().replace(b"\r\n", b"\n"))


def main() -> int:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    cells = read_cells()
    all_summary: list[dict[str, str]] = []
    all_fdr05: list[dict[str, str]] = []
    all_heidi001: list[dict[str, str]] = []
    all_top20: list[dict[str, str]] = []
    all_validation: list[dict[str, str]] = []
    output_fields: list[str] | None = None

    for cell in cells:
        cell_rows: list[dict[str, str]] = []
        chr_qc: list[dict[str, str]] = []
        for chrom in range(1, 23):
            smr_path = (
                FINNGEN_ROOT
                / cell
                / f"chr{chrom}"
                / f"{cell}.chr{chrom}.finngen_r13_discovery_onek1k_params.smr"
            )
            if not smr_path.exists():
                raise FileNotFoundError(smr_path)
            fields, rows = read_smr_file(smr_path, cell, chrom)
            output_fields = fields + ["q_SMR_cell_BH"]
            cell_rows.extend(rows)
            chr_qc.append(
                {
                    "cell": cell,
                    "chr": str(chrom),
                    "smr_rows": str(len(rows)),
                    "smr_path": str(smr_path.relative_to(PROJECT_DIR)),
                }
            )

        q_values = bh_fdr([parse_float(row["p_SMR"]) for row in cell_rows])
        for row, q_value in zip(cell_rows, q_values):
            row["q_SMR_cell_BH"] = fmt(q_value)

        combined_path = FINNGEN_ROOT / cell / f"{cell}.finngen_r13_discovery_onek1k_params.allchr.smr.tsv.gz"
        with gzip.open(combined_path, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=output_fields)
            writer.writeheader()
            writer.writerows(cell_rows)

        chr_qc_path = QC_DIR / f"{cell}.finngen_r13_discovery_onek1k_params.chromosome_qc.tsv"
        write_tsv(chr_qc_path, list(chr_qc[0].keys()), chr_qc)

        top_rows = sorted(cell_rows, key=lambda row: parse_float(row["p_SMR"]))
        fdr05_rows = [row for row in cell_rows if parse_float(row["q_SMR_cell_BH"]) < 0.05]
        heidi001_rows = [
            row
            for row in fdr05_rows
            if not math.isnan(parse_float(row["p_HEIDI"])) and parse_float(row["p_HEIDI"]) >= 0.01
        ]
        meta_lookup = read_meta_lookup(cell)
        for row in sorted(fdr05_rows, key=lambda item: (parse_float(item["q_SMR_cell_BH"]), parse_float(item["p_SMR"]))):
            meta = meta_lookup.get((row["Gene"], row["ProbeChr"]))
            if meta is None:
                continue
            b_finngen = parse_float(row["b_SMR"])
            b_meta = parse_float(meta["b_SMR"])
            p_meta = parse_float(meta["p_SMR"])
            p_heidi_meta = parse_float(meta["p_HEIDI"])
            all_validation.append(
                {
                    "cell": cell,
                    "Gene": row["Gene"],
                    "chr": row["chr"],
                    "topSNP_finngen": row["topSNP"],
                    "p_GWAS_finngen": row["p_GWAS"],
                    "p_eQTL_finngen": row["p_eQTL"],
                    "b_SMR_finngen": row["b_SMR"],
                    "p_SMR_finngen": row["p_SMR"],
                    "q_SMR_finngen_cell_BH": row["q_SMR_cell_BH"],
                    "p_HEIDI_finngen": row["p_HEIDI"],
                    "heidi001_finngen_pass": str(parse_float(row["p_HEIDI"]) >= 0.01 if row["p_HEIDI"] != "NA" else False),
                    "topSNP_meta": meta["topSNP"],
                    "p_GWAS_meta": meta["p_GWAS"],
                    "p_eQTL_meta": meta["p_eQTL"],
                    "b_SMR_meta": meta["b_SMR"],
                    "p_SMR_meta": meta["p_SMR"],
                    "p_HEIDI_meta": meta["p_HEIDI"],
                    "same_direction": str((b_finngen > 0 and b_meta > 0) or (b_finngen < 0 and b_meta < 0)),
                    "meta_nominal_p_SMR_lt_0.05": str(p_meta < 0.05),
                    "meta_heidi001_pass": str(not math.isnan(p_heidi_meta) and p_heidi_meta >= 0.01),
                    "meta_heidi005_pass": str(not math.isnan(p_heidi_meta) and p_heidi_meta >= 0.05),
                }
            )

        all_fdr05.extend(fdr05_rows)
        all_heidi001.extend(heidi001_rows)
        all_top20.extend(top_rows[:20])
        p_smr = [parse_float(row["p_SMR"]) for row in cell_rows]
        p_heidi = [parse_float(row["p_HEIDI"]) for row in cell_rows if row["p_HEIDI"] != "NA"]
        top = top_rows[0]
        all_summary.append(
            {
                "cell": cell,
                "n_smr_rows": str(len(cell_rows)),
                "min_p_SMR": fmt(min(p_smr)),
                "min_q_SMR_cell_BH": fmt(min(parse_float(row["q_SMR_cell_BH"]) for row in cell_rows)),
                "n_cell_BH_FDR_lt_0.05": str(len(fdr05_rows)),
                "n_cell_BH_FDR_lt_0.05_and_HEIDI_ge_0.01": str(len(heidi001_rows)),
                "n_cell_BH_FDR_lt_0.05_and_HEIDI_ge_0.05": str(
                    sum(
                        parse_float(row["q_SMR_cell_BH"]) < 0.05
                        and not math.isnan(parse_float(row["p_HEIDI"]))
                        and parse_float(row["p_HEIDI"]) >= 0.05
                        for row in cell_rows
                    )
                ),
                "n_HEIDI_available": str(len(p_heidi)),
                "n_HEIDI_ge_0.01": str(sum(value >= 0.01 for value in p_heidi)),
                "n_HEIDI_ge_0.05": str(sum(value >= 0.05 for value in p_heidi)),
                "top_gene": top["Gene"],
                "top_snp": top["topSNP"],
                "top_p_GWAS": top["p_GWAS"],
                "top_p_eQTL": top["p_eQTL"],
                "top_p_SMR": top["p_SMR"],
                "top_q_SMR_cell_BH": top["q_SMR_cell_BH"],
                "top_p_HEIDI": top["p_HEIDI"],
            }
        )

    if output_fields is None:
        raise RuntimeError("No SMR output was read")

    summary_path = QC_DIR / "finngen_r13_discovery_onek1k_params.cellwise_bh_fdr_heidi001_summary.tsv"
    fdr05_path = QC_DIR / "finngen_r13_discovery_onek1k_params.cellwise_bh_fdr_lt_0.05_hits.tsv"
    heidi001_path = QC_DIR / "finngen_r13_discovery_onek1k_params.cellwise_bh_fdr_lt_0.05_heidi_ge_0.01_hits.tsv"
    top20_path = QC_DIR / "finngen_r13_discovery_onek1k_params.top20_by_cell.with_cellwise_fdr.tsv"
    validation_path = QC_DIR / "finngen_r13_discovery_onek1k_params.fdr05_hits.meta_validation.tsv"

    write_tsv(summary_path, list(all_summary[0].keys()), sorted(all_summary, key=lambda row: row["cell"]))
    write_tsv(
        fdr05_path,
        output_fields,
        sorted(all_fdr05, key=lambda row: (row["cell"], parse_float(row["q_SMR_cell_BH"]), parse_float(row["p_SMR"]))),
    )
    write_tsv(
        heidi001_path,
        output_fields,
        sorted(all_heidi001, key=lambda row: (row["cell"], parse_float(row["q_SMR_cell_BH"]), parse_float(row["p_SMR"]))),
    )
    write_tsv(
        top20_path,
        output_fields,
        sorted(all_top20, key=lambda row: (row["cell"], parse_float(row["p_SMR"]))),
    )
    validation_fields = list(all_validation[0].keys()) if all_validation else ["cell", "Gene"]
    write_tsv(validation_path, validation_fields, all_validation)

    report_path = REPORT_DIR / "FINNGEN_R13_DISCOVERY_ONEK1K_PARAMS_CELLWISE_FDR_HEIDI001_SUMMARY.md"
    lines = [
        "# FinnGen R13 Discovery SMR Cell-wise FDR Summary",
        "",
        "SMR parameter set: `--smr --peqtl-smr 0.05 --diff-freq-prop 0.99 --thread-num 4`.",
        "Discovery GWAS: FinnGen R13 `I9_VTE`, GRCh37-lifted and harmonized to OneK1K LD.",
        "FDR rule: BH-FDR is calculated within each cell type across chr1-22 SMR rows.",
        "HEIDI sensitivity threshold for prioritized hits: `p_HEIDI >= 0.01`.",
        "",
        "## Cell Summary",
        "",
        "| Cell | Rows | min p_SMR | min cell FDR | FDR<0.05 | FDR<0.05 + HEIDI>=0.01 | FDR<0.05 + HEIDI>=0.05 | Top gene | Top p_HEIDI |",
        "|---|---:|---:|---:|---:|---:|---:|---|---:|",
    ]
    for row in sorted(all_summary, key=lambda item: item["cell"]):
        lines.append(
            f"| {row['cell']} | {int(row['n_smr_rows']):,} | {row['min_p_SMR']} | {row['min_q_SMR_cell_BH']} | "
            f"{row['n_cell_BH_FDR_lt_0.05']} | {row['n_cell_BH_FDR_lt_0.05_and_HEIDI_ge_0.01']} | "
            f"{row['n_cell_BH_FDR_lt_0.05_and_HEIDI_ge_0.05']} | {row['top_gene']} | {row['top_p_HEIDI']} |"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Summary: `{summary_path.relative_to(PROJECT_DIR)}`",
            f"- FDR<0.05 hits: `{fdr05_path.relative_to(PROJECT_DIR)}`",
            f"- FDR<0.05 and HEIDI>=0.01 hits: `{heidi001_path.relative_to(PROJECT_DIR)}`",
            f"- Top20 by cell: `{top20_path.relative_to(PROJECT_DIR)}`",
            f"- Meta validation for FDR<0.05 hits: `{validation_path.relative_to(PROJECT_DIR)}`",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n")

    print(summary_path)
    print(fdr05_path, len(all_fdr05))
    print(heidi001_path, len(all_heidi001))
    print(validation_path, len(all_validation))
    print(report_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
