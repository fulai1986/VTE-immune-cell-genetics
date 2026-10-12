#!/usr/bin/env python3
"""Prepare SMR qfile inputs and sparse BESD files from project eQTL results.

The project eQTL beta is the effect of VCF DS, i.e. ALT allele dosage. SMR A1
is the effect allele, so qfile rows are exported as A1=ALT, A2=REF, Freq=ALT AF.
SNP IDs are allele-specific (chr:pos:REF:ALT) to match the cleaned LD reference.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import subprocess
import os
from pathlib import Path


PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
CONFIG_DIR = PROJECT_DIR / "config"
EQTL_DIR = PROJECT_DIR / "results" / "eqtl"
SMR_DIR = PROJECT_DIR / "results" / "smr"
QC_DIR = PROJECT_DIR / "results" / "qc" / "smr"
SMR_BIN = Path(os.environ.get("SMR_BIN", PROJECT_DIR / "tools" / "smr" / "source" / "SMR" / "build" / "RelWithDebInfo" / "smr"))


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def read_cell_types() -> list[str]:
    path = CONFIG_DIR / "science14_celltypes.tsv"
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row["science14_cell_type"] for row in reader]


def read_paths() -> dict[str, Path]:
    with (CONFIG_DIR / "project_paths.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: Path(row["path"]) for row in reader}


def read_passwords() -> tuple[str, str]:
    try:
        return os.environ["IMPUTE_PASS_A"], os.environ["IMPUTE_PASS_B"]
    except KeyError as exc:
        raise RuntimeError("Set IMPUTE_PASS_A and IMPUTE_PASS_B before reading imputation archives.") from exc


def password_for_chrom(chrom: int, pass_a: str, pass_b: str) -> str:
    return pass_a if chrom in {1, 10, 11, 12} else pass_b


def parse_info(info: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for part in info.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        if key not in {"AF", "MAF", "R2"}:
            continue
        try:
            out[key] = float(value.split(",")[0])
        except ValueError:
            pass
    return out


def finite_float(value: str, default: float = math.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    return f"{value:.12g}"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def load_chr_af(chrom: int, imputation_dir: Path, pass_a: str, pass_b: str) -> dict[str, tuple[str, str, float]]:
    """Return allele-specific SNP ID -> (REF, ALT, ALT AF)."""
    zip_path = imputation_dir / f"chr_{chrom}.zip"
    member = f"chr{chrom}.dose.vcf.gz"
    cmd = ["unzip", "-P", password_for_chrom(chrom, pass_a, pass_b), "-p", str(zip_path), member]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdout is not None
    handle = gzip.GzipFile(fileobj=proc.stdout)
    af: dict[str, tuple[str, str, float]] = {}
    try:
        for raw in handle:
            line = raw.decode("utf-8", errors="replace").rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 8:
                continue
            chrom_value, pos, _snp, ref, alt = parts[:5]
            if "," in alt or ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"}:
                continue
            info = parse_info(parts[7])
            alt_af = info.get("AF", math.nan)
            if not math.isfinite(alt_af):
                continue
            allele_id = f"{chrom_value}:{pos}:{ref}:{alt}"
            af[allele_id] = (ref, alt, alt_af)
    finally:
        try:
            handle.close()
        except Exception:
            pass
        if proc.poll() is None:
            proc.terminate()
    stderr = proc.stderr.read().decode(errors="replace") if proc.stderr else ""
    if proc.wait() not in {0, -15}:
        raise RuntimeError(f"Failed to stream {member} from {zip_path}: {stderr}")
    return af


def qfile_prefix(cell_safe: str, chrom: int) -> Path:
    return SMR_DIR / cell_safe / f"chr{chrom}" / f"{cell_safe}.chr{chrom}.r2_0.8.maf_0.01.allele_id"


def write_qfile(cell_type: str, chrom: int, af_by_id: dict[str, tuple[str, str, float]], min_p: float, overwrite: bool) -> dict[str, int | str]:
    cell_safe = safe_name(cell_type)
    in_path = EQTL_DIR / cell_safe / f"chr{chrom}.full.cis_eqtl.tsv.gz"
    if not in_path.exists():
        raise FileNotFoundError(f"Missing eQTL file: {in_path}")
    out_prefix = qfile_prefix(cell_safe, chrom)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    qfile = Path(str(out_prefix) + ".qfile.txt")
    manifest = Path(str(out_prefix) + ".qfile_manifest.tsv")
    if qfile.exists() and manifest.exists() and not overwrite:
        return {"cell_safe": cell_safe, "chrom": chrom, "status": "qfile_exists"}

    rows = missing_af = missing_se = missing_required = p_floored = 0
    unique_snps: set[str] = set()
    unique_probes: set[str] = set()
    with open_text(in_path) as inp, qfile.open("w", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        required = {
            "chrom",
            "pos",
            "snp_id",
            "ref",
            "alt",
            "gene_id",
            "gene_symbol",
            "gene_start",
            "beta",
            "se",
            "p",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{in_path} missing columns: {', '.join(sorted(missing))}")

        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp", "Gene", "Orientation", "b", "se", "p"])
        for row in reader:
            ref = row["ref"]
            alt = row["alt"]
            allele_id = f'{row["chrom"]}:{row["pos"]}:{ref}:{alt}'
            af_record = af_by_id.get(allele_id)
            if af_record is None:
                missing_af += 1
                continue
            beta = finite_float(row["beta"])
            se = finite_float(row["se"])
            p = finite_float(row["p"])
            if not math.isfinite(beta) or not math.isfinite(se) or se <= 0 or not math.isfinite(p):
                missing_required += 1
                if not math.isfinite(se) or se <= 0:
                    missing_se += 1
                continue
            if p <= 0 or p < min_p:
                p = min_p
                p_floored += 1
            _vcf_ref, _vcf_alt, alt_af = af_record
            probe = row["gene_id"]
            writer.writerow(
                [
                    allele_id,
                    row["chrom"],
                    row["pos"],
                    alt,
                    ref,
                    format_float(alt_af),
                    probe,
                    row["chrom"],
                    row["gene_start"],
                    row["gene_symbol"] or probe,
                    "+",
                    format_float(beta),
                    format_float(se),
                    format_float(p),
                ]
            )
            rows += 1
            unique_snps.add(allele_id)
            unique_probes.add(probe)

    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        writer.writerow(["cell_type", cell_type])
        writer.writerow(["cell_safe", cell_safe])
        writer.writerow(["chrom", chrom])
        writer.writerow(["qfile", qfile])
        writer.writerow(["qfile_rows", rows])
        writer.writerow(["unique_snps", len(unique_snps)])
        writer.writerow(["unique_probes", len(unique_probes)])
        writer.writerow(["rows_skipped_missing_af", missing_af])
        writer.writerow(["rows_skipped_missing_required", missing_required])
        writer.writerow(["rows_skipped_missing_se", missing_se])
        writer.writerow(["p_values_floored", p_floored])
        writer.writerow(["allele_direction", "beta_is_ALT_dosage_effect; A1=ALT; A2=REF; Freq=ALT_AF"])
        writer.writerow(["snp_id_format", "chr:pos:REF:ALT"])

    return {
        "cell_safe": cell_safe,
        "chrom": chrom,
        "status": "qfile_written",
        "qfile_rows": rows,
        "unique_snps": len(unique_snps),
        "unique_probes": len(unique_probes),
        "rows_skipped_missing_af": missing_af,
    }


def run_smr_make_besd(cell_safe: str, chrom: int, overwrite: bool) -> dict[str, int | str]:
    prefix = qfile_prefix(cell_safe, chrom)
    qfile = Path(str(prefix) + ".qfile.txt")
    besd = Path(str(prefix) + ".besd")
    esi = Path(str(prefix) + ".esi")
    epi = Path(str(prefix) + ".epi")
    log = Path(str(prefix) + ".make_besd.log")
    if besd.exists() and esi.exists() and epi.exists() and not overwrite:
        return {"cell_safe": cell_safe, "chrom": chrom, "besd_status": "besd_exists"}
    if not qfile.exists():
        raise FileNotFoundError(f"Missing qfile: {qfile}")
    cmd = [str(SMR_BIN), "--qfile", str(qfile), "--make-besd", "--out", str(prefix)]
    with log.open("w") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    return {"cell_safe": cell_safe, "chrom": chrom, "besd_status": "besd_written"}


def run_query_qc(cell_safe: str, chrom: int, threshold: float, overwrite: bool) -> dict[str, int | str]:
    prefix = qfile_prefix(cell_safe, chrom)
    out_prefix = Path(str(prefix) + f".query_p{threshold:g}")
    out_txt = Path(str(out_prefix) + ".txt")
    log = Path(str(out_prefix) + ".log")
    if out_txt.exists() and not overwrite:
        status = "query_exists"
    else:
        cmd = [str(SMR_BIN), "--beqtl-summary", str(prefix), "--query", str(threshold), "--out", str(out_prefix)]
        with log.open("w") as handle:
            proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)
        status = "query_written"
    n_rows = 0
    if out_txt.exists():
        with out_txt.open() as handle:
            n_rows = max(sum(1 for _ in handle) - 1, 0)
    return {"cell_safe": cell_safe, "chrom": chrom, "query_status": status, "query_rows": n_rows}


def write_summary(rows: list[dict[str, int | str]]) -> Path:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    path = QC_DIR / "smr_qfile_besd_prepare_summary.tsv"
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-type", action="append", help="Science14 cell type; may be repeated. Default: all.")
    parser.add_argument("--chrom", type=int, action="append", help="Chromosome to process; may be repeated. Default: 1-22.")
    parser.add_argument("--min-p", type=float, default=1e-300)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--make-besd", action="store_true", help="Run SMR --qfile --make-besd after qfile export.")
    parser.add_argument("--query-qc", action="store_true", help="Query generated BESD as a readability check.")
    parser.add_argument("--query-threshold", type=float, default=1.0)
    args = parser.parse_args()

    if args.make_besd and not SMR_BIN.exists():
        raise FileNotFoundError(f"Compiled SMR binary not found: {SMR_BIN}")

    paths = read_paths()
    pass_a, pass_b = read_passwords()
    cell_types = args.cell_type or read_cell_types()
    chroms = args.chrom or list(range(1, 23))

    summaries: list[dict[str, int | str]] = []
    for chrom in chroms:
        print(f"Loading chr{chrom} ALT AF from imputed VCF", flush=True)
        af_by_id = load_chr_af(chrom, paths["imputation_dir"], pass_a, pass_b)
        print(f"chr{chrom}: loaded {len(af_by_id)} allele-specific AF records", flush=True)
        for cell_type in cell_types:
            cell_safe = safe_name(cell_type)
            print(f"Preparing {cell_safe} chr{chrom}", flush=True)
            row = write_qfile(cell_type, chrom, af_by_id, args.min_p, args.overwrite)
            if args.make_besd:
                row.update(run_smr_make_besd(cell_safe, chrom, args.overwrite))
            if args.query_qc:
                row.update(run_query_qc(cell_safe, chrom, args.query_threshold, args.overwrite))
            summaries.append(row)
            print(row, flush=True)

    summary_path = write_summary(summaries)
    print(f"Wrote summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
