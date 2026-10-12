#!/usr/bin/env python3
"""Prepare per-cell-type donor covariates for cis-eQTL mapping."""

from __future__ import annotations

import csv
import re
import os
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
CONFIG_DIR = PROJECT_DIR / "config"
QC_DIR = PROJECT_DIR / "results" / "qc"
PSEUDOBULK_DIR = PROJECT_DIR / "results" / "pseudobulk"
METADATA_DIR = PROJECT_DIR / "metadata"
COV_DIR = METADATA_DIR / "covariates"
REPORT_DIR = PROJECT_DIR / "reports"


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def read_paths() -> dict[str, Path]:
    with (CONFIG_DIR / "project_paths.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: Path(row["path"]) for row in reader}


def zscore(values: pd.Series) -> pd.Series:
    arr = pd.to_numeric(values, errors="coerce").astype(float)
    sd = arr.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(arr)), index=values.index)
    return (arr - arr.mean()) / sd


def expression_pcs(expr_path: Path, donors: list[str], max_pcs: int) -> pd.DataFrame:
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr = expr.loc[:, donors]
    mat = expr.to_numpy(dtype=float).T
    mat = np.nan_to_num(mat, nan=0.0)
    mat = mat - mat.mean(axis=0, keepdims=True)
    n_pcs = min(max_pcs, mat.shape[0] - 3, mat.shape[1])
    if n_pcs <= 0:
        return pd.DataFrame({"donor_id": donors})
    u, s, _ = np.linalg.svd(mat, full_matrices=False)
    pc = u[:, :n_pcs] * s[:n_pcs]
    out = pd.DataFrame(pc, columns=[f"exprPC{i}" for i in range(1, n_pcs + 1)])
    out.insert(0, "donor_id", donors)
    return out


def main() -> int:
    COV_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    paths = read_paths()
    genotype = pd.read_csv(paths["genotype_pc_table"], sep="\t", dtype={"donor_id": str, "SampleID": str})
    genotype["sex_male"] = (genotype["sex"].astype(str).str.lower() == "male").astype(int)

    summaries = pd.read_csv(QC_DIR / "all_pseudobulk_summary.tsv", sep="\t")
    summary_rows = []

    for _, row in summaries.iterrows():
        cell_type = row["cell_type"]
        safe = safe_name(cell_type)
        donor_qc_path = QC_DIR / f"{safe}.all_pools.pseudobulk_donor_qc.tsv"
        expr_path = PSEUDOBULK_DIR / f"{safe}.all_pools.expression_rankint.tsv.gz"
        if not donor_qc_path.exists() or not expr_path.exists():
            summary_rows.append(
                {
                    "cell_type": cell_type,
                    "status": "MISSING_INPUT",
                    "n_donors": 0,
                    "n_covariates": 0,
                    "n_expr_pcs": 0,
                    "notes": f"missing {donor_qc_path.name} or {expr_path.name}",
                }
            )
            continue

        donor_qc = pd.read_csv(donor_qc_path, sep="\t", dtype={"individual_id": str}).rename(
            columns={"individual_id": "donor_id"}
        )
        cov = donor_qc.merge(genotype, on="donor_id", how="inner", validate="one_to_one")
        cov = cov.sort_values("donor_id").reset_index(drop=True)
        donors = cov["donor_id"].tolist()

        n_expr_pcs = min(10, len(donors) - 3)
        expr_pc = expression_pcs(expr_path, donors, max_pcs=10)
        cov = cov.merge(expr_pc, on="donor_id", how="left", validate="one_to_one")

        cov_out = pd.DataFrame(
            {
                "donor_id": cov["donor_id"],
                "SampleID": cov["SampleID"],
                "age": cov["age"],
                "sex": cov["sex"],
                "cells_in_celltype": cov["cells_in_celltype"],
                "library_size": cov["library_size"],
                "cov_intercept": 1.0,
                "cov_age_z": zscore(cov["age"]),
                "cov_sex_male": cov["sex_male"],
                "cov_log_cells_z": zscore(np.log1p(cov["cells_in_celltype"])),
                "cov_log_library_size_z": zscore(np.log1p(cov["library_size"])),
            }
        )
        for i in range(1, 11):
            cov_out[f"cov_genotype_PC{i}_z"] = zscore(cov[f"PC{i}"])
        for col in [c for c in cov.columns if c.startswith("exprPC")]:
            cov_out[f"cov_{col}_z"] = zscore(cov[col])

        long_path = COV_DIR / f"{safe}.all_pools.covariates_long.tsv"
        wide_path = COV_DIR / f"{safe}.all_pools.covariates_matrixeqtl.tsv"
        donor_list_path = COV_DIR / f"{safe}.all_pools.donor_list.txt"
        cov_out.to_csv(long_path, sep="\t", index=False)
        donor_list_path.write_text("\n".join(donors) + "\n")

        cov_cols = [c for c in cov_out.columns if c.startswith("cov_")]
        wide = cov_out[["donor_id"] + cov_cols].set_index("donor_id").T
        wide.insert(0, "Covariate", wide.index)
        wide.to_csv(wide_path, sep="\t", index=False)

        summary_rows.append(
            {
                "cell_type": cell_type,
                "status": "OK",
                "n_donors": len(cov_out),
                "n_covariates": len(cov_cols),
                "n_expr_pcs": len([c for c in cov_cols if c.startswith("cov_exprPC")]),
                "notes": "low_donor_count" if len(cov_out) < 200 else "",
            }
        )
        print(f"Wrote covariates for {cell_type}: donors={len(cov_out)} covariates={len(cov_cols)}")

    summary = pd.DataFrame(summary_rows)
    summary_path = QC_DIR / "all_covariate_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    report_path = REPORT_DIR / "STEP3_EQTL_COVARIATES.md"
    lines = [
        "# Step 3 eQTL Covariates",
        "",
        f"- Completed cell types: {(summary['status'] == 'OK').sum()}",
        f"- Summary table: `{summary_path.relative_to(PROJECT_DIR)}`",
        "- Covariate outputs: `metadata/covariates/*.covariates_long.tsv` and `*.covariates_matrixeqtl.tsv`",
        "",
        "| Cell type | Donors | Covariates | Expr PCs | Notes |",
        "|---|---:|---:|---:|---|",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"| {row['cell_type']} | {row['n_donors']} | {row['n_covariates']} | {row['n_expr_pcs']} | {row['notes']} |"
        )
    report_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
