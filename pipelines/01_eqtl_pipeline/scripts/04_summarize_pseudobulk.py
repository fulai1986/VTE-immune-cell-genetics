#!/usr/bin/env python3
"""Summarize all completed all-pools pseudobulk outputs."""

from __future__ import annotations

import csv
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
QC_DIR = PROJECT_DIR / "results" / "qc"
REPORT_DIR = PROJECT_DIR / "reports"


def read_summary(path: Path) -> dict[str, str]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["metric"]: row["value"] for row in reader}


def main() -> int:
    rows = []
    for path in sorted(QC_DIR.glob("*.all_pools.pseudobulk_summary.tsv")):
        summary = read_summary(path)
        rows.append(
            {
                "cell_type": summary.get("target_cell_type", path.name.split(".")[0]),
                "pools": summary.get("pools", ""),
                "cells_retained": int(summary.get("cells_retained", "0")),
                "donors_retained": int(summary.get("donors_retained", "0")),
                "genes_raw": int(summary.get("genes_raw", "0")),
                "genes_expression": int(summary.get("genes_expression", "0")),
                "summary_file": str(path.relative_to(PROJECT_DIR)),
            }
        )

    out = QC_DIR / "all_pseudobulk_summary.tsv"
    with out.open("w", newline="") as handle:
        fieldnames = [
            "cell_type",
            "pools",
            "cells_retained",
            "donors_retained",
            "genes_raw",
            "genes_expression",
            "summary_file",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    report = REPORT_DIR / "STEP2_PSEUDOBULK_SUMMARY.md"
    lines = [
        "# Step 2 Pseudobulk Summary",
        "",
        f"- Completed all-pools cell types: {len(rows)}",
        f"- Summary table: `{out.relative_to(PROJECT_DIR)}`",
        "",
        "| Cell type | Cells | Donors | Expression genes |",
        "|---|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['cell_type']} | {row['cells_retained']} | {row['donors_retained']} | {row['genes_expression']} |"
        )
    report.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out}")
    print(f"Wrote {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
