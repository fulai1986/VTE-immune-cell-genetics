#!/usr/bin/env python3
"""Lift FinnGen R13 I9_VTE summary statistics from GRCh38 to GRCh37."""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path

from pyliftover import LiftOver


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
GWAS_DIR = PROJECT_DIR / "gwas"
FINNGEN_DIR = GWAS_DIR / "replication" / "finngen_r13"
OUT_DIR = FINNGEN_DIR / "liftover_grch37"
QC_DIR = GWAS_DIR / "qc"
REPORT_DIR = GWAS_DIR / "reports"

INPUT = FINNGEN_DIR / "finngen_R13_I9_VTE.gz"
CHAIN = GWAS_DIR / "external" / "ucsc" / "hg38ToHg19.over.chain.gz"
OUTPUT = OUT_DIR / "finngen_R13_I9_VTE.hg19.tsv.gz"
UNMAPPED = OUT_DIR / "finngen_R13_I9_VTE.hg38_unmapped.tsv.gz"
QC_TSV = QC_DIR / "finngen_R13_I9_VTE.hg38_to_hg19_liftover_qc.tsv"
REPORT = REPORT_DIR / "FINNGEN_R13_I9_VTE_LIFTOVER_GRCH37_QC.md"

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def normalize_header(fieldnames: list[str]) -> list[str]:
    return ["chrom" if field == "#chrom" else field for field in fieldnames]


def complement_allele(allele: str) -> str:
    return allele.upper().translate(COMPLEMENT)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    lifter = LiftOver(str(CHAIN))
    counts = {
        "input_rows": 0,
        "lifted_rows": 0,
        "unmapped_rows": 0,
        "multiple_mapping_rows": 0,
        "negative_strand_rows": 0,
        "non_autosome_lifted_rows": 0,
        "bad_input_rows": 0,
    }

    with gzip.open(INPUT, "rt", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"No header found in {INPUT}")
        input_fields = normalize_header(reader.fieldnames)
        output_fields = input_fields + ["original_chrom", "original_pos", "liftover_strand"]
        unmapped_fields = input_fields + ["unmapped_reason"]

        with gzip.open(OUTPUT, "wt", newline="") as fout, gzip.open(UNMAPPED, "wt", newline="") as funmapped:
            writer = csv.DictWriter(fout, delimiter="\t", fieldnames=output_fields)
            unmapped_writer = csv.DictWriter(funmapped, delimiter="\t", fieldnames=unmapped_fields)
            writer.writeheader()
            unmapped_writer.writeheader()

            for raw_row in reader:
                counts["input_rows"] += 1
                row = {("chrom" if key == "#chrom" else key): value for key, value in raw_row.items()}
                chrom = row["chrom"].replace("chr", "")
                try:
                    pos = int(row["pos"])
                except ValueError:
                    row["unmapped_reason"] = "bad_position"
                    unmapped_writer.writerow(row)
                    counts["bad_input_rows"] += 1
                    continue

                mappings = lifter.convert_coordinate(f"chr{chrom}", pos - 1)
                if not mappings:
                    row["unmapped_reason"] = "no_mapping"
                    unmapped_writer.writerow(row)
                    counts["unmapped_rows"] += 1
                    continue
                if len(mappings) > 1:
                    row["unmapped_reason"] = "multiple_mapping"
                    unmapped_writer.writerow(row)
                    counts["multiple_mapping_rows"] += 1
                    continue

                new_chrom_raw, new_pos0, strand, _score = mappings[0]
                new_chrom = new_chrom_raw.replace("chr", "")
                new_pos = int(new_pos0) + 1
                out = dict(row)
                out["original_chrom"] = chrom
                out["original_pos"] = str(pos)
                out["chrom"] = new_chrom
                out["pos"] = str(new_pos)
                out["liftover_strand"] = strand
                out["ref"] = out["ref"].upper()
                out["alt"] = out["alt"].upper()

                if strand == "-":
                    out["ref"] = complement_allele(out["ref"])
                    out["alt"] = complement_allele(out["alt"])
                    counts["negative_strand_rows"] += 1

                if new_chrom not in {str(i) for i in range(1, 23)}:
                    counts["non_autosome_lifted_rows"] += 1
                counts["lifted_rows"] += 1
                writer.writerow(out)

                if counts["input_rows"] % 1_000_000 == 0:
                    print(
                        f"processed={counts['input_rows']} lifted={counts['lifted_rows']} "
                        f"unmapped={counts['unmapped_rows']} multiple={counts['multiple_mapping_rows']}",
                        flush=True,
                    )

    with QC_TSV.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in counts.items():
            writer.writerow([key, value])
        writer.writerow(["input", INPUT.relative_to(PROJECT_DIR)])
        writer.writerow(["chain", CHAIN.relative_to(PROJECT_DIR)])
        writer.writerow(["output", OUTPUT.relative_to(PROJECT_DIR)])
        writer.writerow(["unmapped", UNMAPPED.relative_to(PROJECT_DIR)])

    lines = [
        "# FinnGen R13 I9_VTE GRCh38 to GRCh37 Liftover QC",
        "",
        f"- Input: `{INPUT.relative_to(PROJECT_DIR)}`",
        f"- Chain: `{CHAIN.relative_to(PROJECT_DIR)}`",
        f"- Output: `{OUTPUT.relative_to(PROJECT_DIR)}`",
        f"- Unmapped: `{UNMAPPED.relative_to(PROJECT_DIR)}`",
        "",
        "## Counts",
        "",
        f"- Input rows: {counts['input_rows']:,}",
        f"- Lifted rows: {counts['lifted_rows']:,}",
        f"- Unmapped rows: {counts['unmapped_rows']:,}",
        f"- Multiple mapping rows: {counts['multiple_mapping_rows']:,}",
        f"- Negative-strand rows requiring allele complement: {counts['negative_strand_rows']:,}",
        f"- Non-autosome lifted rows: {counts['non_autosome_lifted_rows']:,}",
        f"- Bad input rows: {counts['bad_input_rows']:,}",
        "",
        "## Notes",
        "",
        "- Input positions are one-based GRCh38 coordinates.",
        "- `pyliftover` uses zero-based coordinates internally, so positions are converted as `pos - 1` before liftover and `+ 1` after liftover.",
        "- For negative-strand mappings, `ref` and `alt` are complemented. `beta` and `af_alt` remain attached to the converted ALT allele label.",
        "- The lifted file still requires allele/ID harmonization to the OneK1K GRCh37 LD reference before SMR.",
        "",
        f"QC table: `{QC_TSV.relative_to(PROJECT_DIR)}`",
    ]
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUTPUT}")
    print(f"Wrote {UNMAPPED}")
    print(f"Wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
