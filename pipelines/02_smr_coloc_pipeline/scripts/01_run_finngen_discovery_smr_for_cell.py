#!/usr/bin/env python3
"""Run FinnGen R13 discovery SMR for one OneK1K PBMC cell type across chr1-22."""
from __future__ import annotations
import argparse, subprocess
import os
from pathlib import Path
PROJECT_DIR = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
SMR = Path(os.environ.get("SMR_BIN", PROJECT_DIR / "tools/smr/source/SMR/build/RelWithDebInfo/smr"))
GWAS = PROJECT_DIR / "gwas/replication/finngen_r13/harmonized_oneK1K_ld/finngen_R13_I9_VTE.hg19.onek1k_ld_harmonized.smr.ma"
LD_ROOT = PROJECT_DIR / "results/ld_reference/onek1k"
EQTL_ROOT = PROJECT_DIR / "results/smr"
OUT_ROOT = PROJECT_DIR / "results/smr_analysis/finngen_r13_discovery_onek1k_params"
LOG_ROOT = PROJECT_DIR / "logs/smr_analysis/finngen_r13_discovery_onek1k_params"
def output_prefix(cell: str, chrom: int) -> Path:
    return OUT_ROOT / cell / f"chr{chrom}" / f"{cell}.chr{chrom}.finngen_r13_discovery_onek1k_params"
def smr_result_path(cell: str, chrom: int) -> Path:
    return Path(str(output_prefix(cell, chrom)) + ".smr")
def run_chr(cell: str, chrom: int, thread_num: int, force: bool) -> tuple[int, Path]:
    out_prefix = output_prefix(cell, chrom)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    log_path = LOG_ROOT / f"{cell}.chr{chrom}.finngen_r13_discovery_onek1k_params.smr.log"
    if smr_result_path(cell, chrom).exists() and smr_result_path(cell, chrom).stat().st_size > 0 and not force:
        return 0, log_path
    bfile = LD_ROOT / f"chr{chrom}" / f"onek1k.chr{chrom}.r2_0.8.maf_0.01.allele_id"
    beqtl = EQTL_ROOT / cell / f"chr{chrom}" / f"{cell}.chr{chrom}.r2_0.8.maf_0.01.allele_id"
    cmd = [str(SMR), "--bfile", str(bfile), "--gwas-summary", str(GWAS), "--beqtl-summary", str(beqtl), "--out", str(out_prefix), "--thread-num", str(thread_num), "--smr", "--peqtl-smr", "0.05", "--diff-freq-prop", "0.99"]
    with log_path.open("w") as log:
        done = subprocess.run(cmd, cwd=PROJECT_DIR, stdout=log, stderr=subprocess.STDOUT, check=False)
    return done.returncode, log_path
def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument('--cell', required=True); ap.add_argument('--thread-num', type=int, default=4); ap.add_argument('--force', action='store_true'); args=ap.parse_args()
    failures=[]
    for chrom in range(1,23):
        rc, log=run_chr(args.cell, chrom, args.thread_num, args.force)
        out=smr_result_path(args.cell, chrom)
        rows=sum(1 for _ in out.open())-1 if out.exists() else -1
        print(f"chr{chrom}	{'OK' if rc==0 else f'FAIL({rc})'}	rows={rows}	log={log.relative_to(PROJECT_DIR)}", flush=True)
        if rc != 0: failures.append((chrom, log))
    if failures:
        for chrom, log in failures: print(f"chr{chrom}: {log}")
        return 1
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
