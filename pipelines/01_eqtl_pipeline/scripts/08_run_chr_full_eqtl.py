#!/usr/bin/env python3
"""Run full chromosome cis-eQTL for one cell type with cell-type-scoped outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
CONFIG_DIR = PROJECT_DIR / "config"
PSEUDOBULK_DIR = PROJECT_DIR / "results" / "pseudobulk"
EQTL_DIR = PROJECT_DIR / "results" / "eqtl"
QC_DIR = PROJECT_DIR / "results" / "qc" / "eqtl"
REPORT_DIR = PROJECT_DIR / "reports" / "eqtl"
METADATA_DIR = PROJECT_DIR / "metadata"


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def read_paths() -> dict[str, Path]:
    with (CONFIG_DIR / "project_paths.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: Path(row["path"]) for row in reader}


def read_passwords() -> tuple[str, str]:
    try:
        return os.environ["IMPUTE_PASS_A"], os.environ["IMPUTE_PASS_B"]
    except KeyError as exc:
        raise RuntimeError("Set IMPUTE_PASS_A and IMPUTE_PASS_B before reading imputation archives.") from exc


def password_for_chrom(chrom: int) -> str:
    pass_a, pass_b = read_passwords()
    return pass_a if chrom in {1, 10, 11, 12} else pass_b


def parse_info(info: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for part in info.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        if key in {"AF", "MAF", "R2"}:
            try:
                out[key] = float(value)
            except ValueError:
                pass
    return out


def stream_vcf(chrom: int, imputation_dir: Path):
    cmd = [
        "unzip",
        "-P",
        password_for_chrom(chrom),
        "-p",
        str(imputation_dir / f"chr_{chrom}.zip"),
        f"chr{chrom}.dose.vcf.gz",
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdout is not None
    handle = gzip.GzipFile(fileobj=proc.stdout)
    try:
        for raw in handle:
            yield raw.decode("utf-8", errors="replace").rstrip("\n")
    finally:
        try:
            handle.close()
        except Exception:
            pass
        if proc.poll() is None:
            proc.terminate()


def dosage_from_sample(sample: str, ds_index: int) -> float:
    parts = sample.split(":")
    if ds_index >= len(parts) or parts[ds_index] in {".", ""}:
        return np.nan
    try:
        return float(parts[ds_index])
    except ValueError:
        return np.nan


def residualize_matrix(matrix: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    q, _ = np.linalg.qr(covariates)
    return matrix - q @ (q.T @ matrix)


def residualize_vector(vector: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    q, _ = np.linalg.qr(covariates)
    return vector - q @ (q.T @ vector)


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    n = len(pvalues)
    if n == 0:
        return pvalues
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def write_rows(path: Path, rows: list[dict], header_written: bool) -> bool:
    if not rows:
        return header_written
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False, mode="a", header=not header_written, compression="gzip")
    return True


def run(args: argparse.Namespace) -> int:
    paths = read_paths()
    cell_safe = safe_name(args.cell_type)
    cell_eqtl_dir = EQTL_DIR / cell_safe
    cell_qc_dir = QC_DIR / cell_safe
    cell_report_dir = REPORT_DIR / cell_safe
    cell_eqtl_dir.mkdir(parents=True, exist_ok=True)
    cell_qc_dir.mkdir(parents=True, exist_ok=True)
    cell_report_dir.mkdir(parents=True, exist_ok=True)

    prefix = f"chr{args.chrom}.full"
    out_path = cell_eqtl_dir / f"{prefix}.cis_eqtl.tsv.gz"
    top_path = cell_eqtl_dir / f"{prefix}.top_by_gene.tsv"
    summary_path = cell_qc_dir / f"{prefix}.summary.tsv"
    report_path = cell_report_dir / f"{prefix}.md"
    if out_path.exists() and not args.overwrite:
        raise SystemExit(f"Output exists; use --overwrite to replace: {out_path}")
    for path in [out_path, top_path, summary_path, report_path]:
        if path.exists() and args.overwrite:
            path.unlink()

    print(f"Reading expression for {args.cell_type}")
    expr = pd.read_csv(PSEUDOBULK_DIR / f"{cell_safe}.all_pools.expression_rankint.tsv.gz", sep="\t", index_col=0)
    cov = pd.read_csv(METADATA_DIR / "covariates" / f"{cell_safe}.all_pools.covariates_long.tsv", sep="\t", dtype={"donor_id": str})
    sample_map = pd.read_csv(METADATA_DIR / "imputed_vcf_sample_to_donor_map.tsv", sep="\t", dtype={"donor_id": str})
    sample_map = sample_map[sample_map["has_donor_id"].astype(str).str.lower().isin(["true", "1"])]

    gene_annot = pd.read_csv(paths["gene_annotation"], sep="\t")
    genes = gene_annot[(gene_annot["chr"].astype(str) == str(args.chrom)) & (gene_annot["Probe"].isin(expr.index))].copy()
    genes = genes.sort_values(["start", "Probe"]).reset_index(drop=True)
    expr = expr.loc[genes["Probe"]]

    common_donors = [d for d in cov["donor_id"].tolist() if d in expr.columns]
    cov = cov.set_index("donor_id").loc[common_donors].reset_index()
    expr = expr.loc[:, common_donors]
    cov_cols = [c for c in cov.columns if c.startswith("cov_")]
    x = cov[cov_cols].to_numpy(dtype=float)
    y = expr.to_numpy(dtype=float).T
    y_resid = residualize_matrix(y, x)
    y_norm = np.sqrt(np.sum(y_resid * y_resid, axis=0))

    donor_to_col = {d: i for i, d in enumerate(common_donors)}
    vcf_to_donor = dict(zip(sample_map["vcf_sample"], sample_map["donor_id"]))
    gene_starts = genes["start"].to_numpy(dtype=int)
    gene_ends = genes["end"].to_numpy(dtype=int)
    gene_ids = genes["Probe"].tolist()
    gene_symbols = genes["gene_symbol"].tolist()

    header_written = False
    buffer: list[dict] = []
    top_by_gene: dict[str, dict] = {}
    p_for_fdr = []
    variants_seen = variants_pass_qc = variants_tested = cis_pairs = 0
    header_samples = sample_indices = sample_donors = None

    print(f"Streaming chr{args.chrom}; genes={len(genes)} donors={len(common_donors)}")
    for line in stream_vcf(args.chrom, paths["imputation_dir"]):
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header_samples = line.split("\t")[9:]
            sample_indices = []
            sample_donors = []
            for i, sample in enumerate(header_samples):
                donor = vcf_to_donor.get(sample)
                if donor in donor_to_col:
                    sample_indices.append(i)
                    sample_donors.append(donor)
            print(f"Matched VCF samples to {len(sample_indices)} expression donors")
            continue
        if header_samples is None:
            continue

        parts = line.split("\t")
        variants_seen += 1
        pos = int(parts[1])
        info = parse_info(parts[7])
        r2 = info.get("R2", np.nan)
        maf_info = info.get("MAF", np.nan)
        if not np.isfinite(r2) or r2 < args.min_r2 or not np.isfinite(maf_info) or maf_info < args.min_maf:
            continue
        variants_pass_qc += 1
        cis_mask = (gene_starts <= pos + args.cis_window) & (gene_ends >= pos - args.cis_window)
        if not cis_mask.any():
            continue
        fmt = parts[8].split(":")
        if "DS" not in fmt:
            continue
        ds_index = fmt.index("DS")
        dosage = np.full(len(common_donors), np.nan)
        for sample_idx, donor in zip(sample_indices, sample_donors):
            dosage[donor_to_col[donor]] = dosage_from_sample(parts[9 + sample_idx], ds_index)
        ok = np.isfinite(dosage)
        if ok.sum() < max(50, len(common_donors) * 0.8):
            continue
        mean_dosage = np.nanmean(dosage)
        dosage[~ok] = mean_dosage
        af = mean_dosage / 2.0
        maf = min(af, 1.0 - af)
        if maf < args.min_maf:
            continue
        g_resid = residualize_vector(dosage, x)
        g_norm_sq = float(np.dot(g_resid, g_resid))
        if g_norm_sq <= 0 or not np.isfinite(g_norm_sq):
            continue
        g_norm = math.sqrt(g_norm_sq)
        gene_idx = np.where(cis_mask)[0]
        dots = g_resid @ y_resid[:, gene_idx]
        corr = np.clip(dots / (g_norm * y_norm[gene_idx]), -0.999999, 0.999999)
        df = len(common_donors) - x.shape[1] - 2
        t_stat = corr * np.sqrt(df / np.maximum(1e-12, 1 - corr * corr))
        pvals = 2 * stats.t.sf(np.abs(t_stat), df=df)
        beta = dots / g_norm_sq
        se = np.abs(beta / t_stat)
        se[~np.isfinite(se)] = np.nan
        snp_id = parts[2] if parts[2] != "." else f"{parts[0]}:{parts[1]}:{parts[3]}:{parts[4]}"

        for local_i, gene_i in enumerate(gene_idx):
            row = {
                "cell_type": args.cell_type,
                "chrom": args.chrom,
                "pos": pos,
                "snp_id": snp_id,
                "ref": parts[3],
                "alt": parts[4],
                "maf": maf,
                "imputation_r2": r2,
                "gene_id": gene_ids[gene_i],
                "gene_symbol": gene_symbols[gene_i],
                "gene_start": int(gene_starts[gene_i]),
                "gene_end": int(gene_ends[gene_i]),
                "n_donors": len(common_donors),
                "n_covariates": x.shape[1],
                "beta": beta[local_i],
                "se": se[local_i],
                "t": t_stat[local_i],
                "p": pvals[local_i],
            }
            buffer.append(row)
            p_for_fdr.append(float(pvals[local_i]))
            current_top = top_by_gene.get(row["gene_id"])
            if current_top is None or row["p"] < current_top["p"]:
                top_by_gene[row["gene_id"]] = row.copy()
        variants_tested += 1
        cis_pairs += len(gene_idx)
        if len(buffer) >= args.write_buffer:
            header_written = write_rows(out_path, buffer, header_written)
            buffer.clear()
        if variants_tested % args.progress_every == 0:
            print(f"tested_variants={variants_tested}; cis_pairs={cis_pairs}; variants_seen={variants_seen}")

    header_written = write_rows(out_path, buffer, header_written)
    top = pd.DataFrame(top_by_gene.values())
    if len(top):
        top["fdr_within_top_by_gene"] = bh_fdr(top["p"].to_numpy(dtype=float))
        top = top.sort_values("p")
    top.to_csv(top_path, sep="\t", index=False)
    min_p = min(p_for_fdr) if p_for_fdr else np.nan

    summary = pd.DataFrame(
        [
            {"metric": "cell_type", "value": args.cell_type},
            {"metric": "chrom", "value": args.chrom},
            {"metric": "donors", "value": len(common_donors)},
            {"metric": "genes_considered", "value": len(genes)},
            {"metric": "covariates", "value": x.shape[1]},
            {"metric": "variants_seen", "value": variants_seen},
            {"metric": "variants_pass_qc", "value": variants_pass_qc},
            {"metric": "variants_tested_with_cis_genes", "value": variants_tested},
            {"metric": "cis_pairs_tested", "value": cis_pairs},
            {"metric": "top_genes", "value": len(top)},
            {"metric": "min_p", "value": min_p},
            {"metric": "full_result", "value": str(out_path.relative_to(PROJECT_DIR))},
            {"metric": "top_by_gene", "value": str(top_path.relative_to(PROJECT_DIR))},
        ]
    )
    summary.to_csv(summary_path, sep="\t", index=False)
    report_path.write_text(
        "\n".join(
            [
                f"# Full cis-eQTL {args.cell_type} chr{args.chrom}",
                "",
                f"- Donors: {len(common_donors)}",
                f"- Genes considered: {len(genes)}",
                f"- Variants seen: {variants_seen}",
                f"- Variants passing QC: {variants_pass_qc}",
                f"- Variants tested with cis genes: {variants_tested}",
                f"- Cis pairs tested: {cis_pairs}",
                f"- Top genes: {len(top)}",
                f"- Minimum P: {min_p}",
                f"- Full result: `{out_path.relative_to(PROJECT_DIR)}`",
                f"- Top by gene: `{top_path.relative_to(PROJECT_DIR)}`",
                f"- Summary: `{summary_path.relative_to(PROJECT_DIR)}`",
            ]
        )
        + "\n"
    )
    print(f"Wrote {out_path}")
    print(f"Wrote {top_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--chrom", type=int, required=True)
    parser.add_argument("--cis-window", type=int, default=1_000_000)
    parser.add_argument("--min-r2", type=float, default=0.8)
    parser.add_argument("--min-maf", type=float, default=0.01)
    parser.add_argument("--write-buffer", type=int, default=100_000)
    parser.add_argument("--progress-every", type=int, default=10_000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(run(parse_args()))
