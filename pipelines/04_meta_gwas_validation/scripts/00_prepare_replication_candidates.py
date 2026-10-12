#!/usr/bin/env python3
"""Export the nine FinnGen Tier 1-3 candidates for independent replication."""
import os
from pathlib import Path
import pandas as pd


def main():
    base = Path(os.environ["PROJECT_DIR"]).expanduser().resolve() / "results/smr_analysis/qc"
    frame = pd.read_csv(base / "finngen_r13_discovery_fdr05_hits.evidence_tiers_with_coloc_abf_beta_or.tsv", sep="\t")
    frame = frame[frame["tier"].isin(["Tier 1", "Tier 2", "Tier 3"])].copy()
    if len(frame) != 9 or frame.duplicated(["cell", "Gene"]).any():
        raise ValueError("Expected nine distinct prioritized cell-gene pairs")
    frame = frame.rename(columns={"smr_p_SMR": "p_SMR_finngen", "PP.H4.abf": "coloc_abf_PP4"})
    columns = ["tier", "cell", "Gene", "chr", "topSNP", "b_SMR_finngen", "OR_SMR_finngen", "p_SMR_finngen", "smr_q_SMR_cell_BH", "smr_p_HEIDI", "coloc_abf_PP4", "b_SMR_meta", "OR_SMR_meta", "meta_p_SMR", "meta_p_HEIDI", "same_direction"]
    frame = frame[[column for column in columns if column in frame]]
    frame = frame.rename(columns={"meta_p_SMR": "p_SMR_meta", "meta_p_HEIDI": "p_HEIDI_meta"})
    frame.to_csv(base / "finngen_r13_discovery_tier1_3.with_smr_beta_or.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
