#!/usr/bin/env python3

from pathlib import Path
import gzip
import os
import pandas as pd


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
TIER_FILE = PROJECT / "results/smr_analysis/qc/finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf_beta_or.tsv"
FINNGEN = PROJECT / "gwas/replication/finngen_r13/harmonized_oneK1K_ld/finngen_R13_I9_VTE.hg19.onek1k_ld_harmonized.tsv.gz"
OUT = PROJECT / "manuscript_results/tables/top_snp_rsid_mapping_hg19.tsv"


def main():
    tier = pd.read_csv(TIER_FILE, sep="\t")
    wanted = set(tier["topSNP"].astype(str))
    found = {}
    with gzip.open(FINNGEN, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        required = ["SNP", "CHR", "BP", "A1", "A2", "rsids", "original_chrom", "original_pos"]
        missing = [name for name in required if name not in idx]
        if missing:
            raise RuntimeError(f"Missing FinnGen columns: {missing}")
        for line in handle:
            row = line.rstrip("\n").split("\t")
            snp = row[idx["SNP"]]
            if snp not in wanted:
                continue
            found[snp] = {
                "topSNP": snp,
                "rsid": row[idx["rsids"]].split(",")[0],
                "CHR_hg19": row[idx["CHR"]],
                "BP_hg19": row[idx["BP"]],
                "A1": row[idx["A1"]],
                "A2": row[idx["A2"]],
                "original_chrom_GRCh38": row[idx["original_chrom"]],
                "original_pos_GRCh38": row[idx["original_pos"]],
            }
            if len(found) == len(wanted):
                break
    missing = sorted(wanted - set(found))
    if missing:
        raise RuntimeError(f"Missing rsID mappings for: {missing}")
    out = pd.DataFrame([found[snp] for snp in sorted(found)])
    out.to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
