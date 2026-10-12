#!/usr/bin/env python3
"""Run targeted meta-GWAS SMR validation for FinnGen Tier 1-3 signals.

Uses the new standard VTE meta-GWAS with UKB approximate log(OR) conversion and
`n_cohorts >= 3`, harmonized to the OneK1K LD/eQTL allele convention.
"""

from __future__ import annotations

import csv
import gzip
import math
import subprocess
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
SMR = Path(os.environ.get("SMR_BIN", PROJECT_DIR / "tools" / "smr" / "source" / "SMR" / "build" / "RelWithDebInfo" / "smr"))
GWAS = (
    PROJECT_DIR
    / "gwas"
    / "harmonized"
    / "oneK1K_ld"
    / "vte_discovery_finngen_independent_ivw_ukb_approx_logor.primary_ncohort3plus.onek1k_ld_harmonized.smr.ma"
)
LD_ROOT = PROJECT_DIR / "results" / "ld_reference" / "onek1k"
EQTL_ROOT = PROJECT_DIR / "results" / "smr"
OUT_ROOT = PROJECT_DIR / "results" / "smr_analysis" / "meta_ncohort3plus_ukb_approx_logor_validation"
LOG_ROOT = PROJECT_DIR / "logs" / "smr_analysis" / "meta_ncohort3plus_ukb_approx_logor_validation"
QC_DIR = PROJECT_DIR / "results" / "smr_analysis" / "qc"
REPORT_DIR = PROJECT_DIR / "reports" / "smr_analysis"

CANDIDATE_TABLE = QC_DIR / "finngen_r13_discovery_tier1_3.with_smr_beta_or.tsv"

RUN_LABEL = "meta_ncohort3plus_ukb_approx_logor"


def parse_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def read_candidates() -> list[dict[str, str]]:
    path = CANDIDATE_TABLE
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tier = row.get("tier", "").strip('"')
            if tier not in {"Tier 1", "Tier 2", "Tier 3"}:
                continue
            rows.append(row)
    return rows


def output_prefix(cell: str, chrom: int) -> Path:
    return OUT_ROOT / cell / f"chr{chrom}" / f"{cell}.chr{chrom}.{RUN_LABEL}"


def smr_result_path(cell: str, chrom: int) -> Path:
    return Path(str(output_prefix(cell, chrom)) + ".smr")


def run_chr(cell: str, chrom: int, thread_num: int = 4) -> tuple[int, Path]:
    out_prefix = output_prefix(cell, chrom)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    log_path = LOG_ROOT / f"{cell}.chr{chrom}.{RUN_LABEL}.smr.log"
    result_path = smr_result_path(cell, chrom)
    if result_path.exists() and result_path.stat().st_size > 0:
        return 0, log_path

    bfile = LD_ROOT / f"chr{chrom}" / f"onek1k.chr{chrom}.r2_0.8.maf_0.01.allele_id"
    beqtl = EQTL_ROOT / cell / f"chr{chrom}" / f"{cell}.chr{chrom}.r2_0.8.maf_0.01.allele_id"
    cmd = [
        str(SMR),
        "--bfile",
        str(bfile),
        "--gwas-summary",
        str(GWAS),
        "--beqtl-summary",
        str(beqtl),
        "--out",
        str(out_prefix),
        "--thread-num",
        str(thread_num),
        "--smr",
        "--peqtl-smr",
        "0.05",
        "--diff-freq-prop",
        "0.99",
    ]
    with log_path.open("w") as log:
        completed = subprocess.run(cmd, cwd=PROJECT_DIR, stdout=log, stderr=subprocess.STDOUT, check=False)
    return completed.returncode, log_path


def read_smr(cell: str, chrom: int) -> list[dict[str, str]]:
    path = smr_result_path(cell, chrom)
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sign(value: float) -> int:
    if not math.isfinite(value) or value == 0:
        return 0
    return 1 if value > 0 else -1


def summarize(candidates: list[dict[str, str]], run_qc: list[dict[str, str]]) -> None:
    out_path = QC_DIR / "finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_validation.tsv"
    display_path = (
        PROJECT_DIR
        / "manuscript_results"
        / "completed_results"
        / "result2_smr_coloc"
        / "source_tables"
        / "tier1_3_meta_ncohort3plus_ukb_approx_logor_validation_display.tsv"
    )
    report_path = REPORT_DIR / "FINNGEN_R13_DISCOVERY_TIER1_3_META_NCOHORT3PLUS_UKB_APPROX_LOGOR_VALIDATION.md"
    display_path.parent.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    output_rows: list[dict[str, str]] = []
    display_rows: list[dict[str, str]] = []
    for cand in candidates:
        cell = cand["cell"].strip('"')
        gene = cand["Gene"].strip('"')
        tier = cand["tier"].strip('"')
        chrom = int(cand["chr"].strip('"'))
        top_snp = cand["topSNP"].strip('"')
        finngen_b = parse_float(cand.get("b_SMR_finngen", cand.get("b_SMR", "")))
        if not math.isfinite(finngen_b):
            finngen_or = parse_float(cand.get("OR_SMR_finngen", ""))
            finngen_b = math.log(finngen_or) if math.isfinite(finngen_or) and finngen_or > 0 else math.nan

        rows = read_smr(cell, chrom)
        gene_rows = [row for row in rows if row.get("Gene") == gene or row.get("ProbeID") == gene]
        if not gene_rows:
            out = {
                **cand,
                "meta_ukb_approx_found": "False",
                "meta_ukb_approx_b_SMR": "NA",
                "meta_ukb_approx_se_SMR": "NA",
                "meta_ukb_approx_OR": "NA",
                "meta_ukb_approx_OR_95CI": "NA",
                "meta_ukb_approx_p_SMR": "NA",
                "meta_ukb_approx_p_HEIDI": "NA",
                "meta_ukb_approx_same_direction": "NA",
                "meta_ukb_approx_nominal_p_SMR_lt_0.05": "False",
                "meta_ukb_approx_heidi001_pass": "False",
                "meta_ukb_approx_validation_pass": "False",
            }
        else:
            row = min(gene_rows, key=lambda x: parse_float(x.get("p_SMR", "nan")))
            b = parse_float(row.get("b_SMR", ""))
            se = parse_float(row.get("se_SMR", ""))
            p_smr = parse_float(row.get("p_SMR", ""))
            p_heidi = parse_float(row.get("p_HEIDI", ""))
            same_direction = sign(b) != 0 and sign(b) == sign(finngen_b)
            nominal = math.isfinite(p_smr) and p_smr < 0.05
            heidi001 = math.isfinite(p_heidi) and p_heidi >= 0.01
            or_value = math.exp(b) if math.isfinite(b) else math.nan
            lo = math.exp(b - 1.96 * se) if math.isfinite(b) and math.isfinite(se) else math.nan
            hi = math.exp(b + 1.96 * se) if math.isfinite(b) and math.isfinite(se) else math.nan
            out = {
                **cand,
                "meta_ukb_approx_found": "True",
                "meta_ukb_approx_probe": row.get("ProbeID", ""),
                "meta_ukb_approx_topSNP": row.get("topSNP", ""),
                "meta_ukb_approx_b_SMR": f"{b:.8g}",
                "meta_ukb_approx_se_SMR": f"{se:.8g}",
                "meta_ukb_approx_OR": f"{or_value:.8g}",
                "meta_ukb_approx_OR_95CI": f"{or_value:.2f} ({lo:.2f}-{hi:.2f})",
                "meta_ukb_approx_p_SMR": f"{p_smr:.8g}",
                "meta_ukb_approx_p_HEIDI": f"{p_heidi:.8g}" if math.isfinite(p_heidi) else "NA",
                "meta_ukb_approx_same_direction": str(same_direction),
                "meta_ukb_approx_nominal_p_SMR_lt_0.05": str(nominal),
                "meta_ukb_approx_heidi001_pass": str(heidi001),
                "meta_ukb_approx_validation_pass": str(same_direction and nominal and heidi001),
            }
        output_rows.append(out)
        display_rows.append(
            {
                "Tier": tier,
                "Cell-gene pair": f"{cell.replace('_', ' ')} - {gene}",
                "FinnGen top SNP": top_snp,
                "FinnGen direction": "Risk" if sign(finngen_b) > 0 else "Protective" if sign(finngen_b) < 0 else "NA",
                "Meta OR (95% CI)": out["meta_ukb_approx_OR_95CI"],
                "Meta p_SMR": out["meta_ukb_approx_p_SMR"],
                "Meta p_HEIDI": out["meta_ukb_approx_p_HEIDI"],
                "Direction consistent": out["meta_ukb_approx_same_direction"],
                "Meta validation": "Pass" if out["meta_ukb_approx_validation_pass"] == "True" else "Not pass",
            }
        )

    fields = list(output_rows[0].keys())
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)

    display_fields = list(display_rows[0].keys())
    with display_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=display_fields)
        writer.writeheader()
        writer.writerows(display_rows)

    counts: dict[str, tuple[int, int]] = {}
    for row in output_rows:
        tier = row["tier"].strip('"')
        n, p = counts.get(tier, (0, 0))
        counts[tier] = (n + 1, p + (1 if row["meta_ukb_approx_validation_pass"] == "True" else 0))
    total_pass = sum(1 for row in output_rows if row["meta_ukb_approx_validation_pass"] == "True")
    lines = [
        "# FinnGen Tier 1-3 validation in meta-GWAS with UKB approximate log(OR)",
        "",
        f"- GWAS: `{GWAS.relative_to(PROJECT_DIR)}`",
        "- SMR parameters: `--smr --peqtl-smr 0.05 --diff-freq-prop 0.99 --thread-num 4`.",
        "- Validation rule: same SMR direction as FinnGen, meta `p_SMR < 0.05`, and meta `p_HEIDI >= 0.01`.",
        "",
        "## Counts",
        "",
        f"- Total Tier 1-3 candidates: {len(output_rows)}",
        f"- Meta validated: {total_pass}",
        "",
        "| Tier | N | Pass |",
        "|---|---:|---:|",
    ]
    for tier in ["Tier 1", "Tier 2", "Tier 3"]:
        n, p = counts.get(tier, (0, 0))
        lines.append(f"| {tier} | {n} | {p} |")
    lines.extend(["", "## SMR run QC", "", "| Cell | Chr | Return code | Rows | Log |", "|---|---:|---:|---:|---|"])
    for row in run_qc:
        lines.append(f"| {row['cell']} | {row['chr']} | {row['return_code']} | {row['rows']} | `{row['log']}` |")
    lines.extend(["", "## Outputs", "", f"- Full validation table: `{out_path.relative_to(PROJECT_DIR)}`", f"- Display table: `{display_path.relative_to(PROJECT_DIR)}`"])
    report_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_path}")
    print(f"Wrote {display_path}")
    print(f"Wrote {report_path}")


def main() -> int:
    candidates = read_candidates()
    combos = sorted({(row["cell"].strip('"'), int(row["chr"].strip('"'))) for row in candidates})
    run_qc: list[dict[str, str]] = []
    failures = []
    for cell, chrom in combos:
        rc, log_path = run_chr(cell, chrom)
        result_path = smr_result_path(cell, chrom)
        rows = sum(1 for _ in result_path.open()) - 1 if result_path.exists() else -1
        run_qc.append(
            {
                "cell": cell,
                "chr": str(chrom),
                "return_code": str(rc),
                "rows": str(rows),
                "log": str(log_path.relative_to(PROJECT_DIR)),
            }
        )
        print(f"{cell}\tchr{chrom}\trc={rc}\trows={rows}", flush=True)
        if rc != 0:
            failures.append((cell, chrom, log_path))
    if failures:
        for cell, chrom, log_path in failures:
            print(f"FAIL\t{cell}\tchr{chrom}\t{log_path}")
        return 1
    summarize(candidates, run_qc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
