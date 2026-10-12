#!/usr/bin/env python3
"""Apply candidate-set cell-wise BH-FDR to replication SMR results."""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


def main():
    base = Path(os.environ["PROJECT_DIR"]).expanduser().resolve() / "results/smr_analysis/qc"
    stem = "finngen_r13_discovery_tier1_3.meta_ncohort3plus_ukb_approx_logor_validation"
    frame = pd.read_csv(base / f"{stem}.tsv", sep="\t")
    p = pd.to_numeric(frame["meta_ukb_approx_p_SMR"], errors="coerce")
    frame["meta_ukb_approx_q_SMR_cell_BH"] = np.nan
    for indices in frame.groupby("cell", sort=False).groups.values():
        valid = p.loc[indices].dropna()
        if len(valid):
            frame.loc[valid.index, "meta_ukb_approx_q_SMR_cell_BH"] = false_discovery_control(valid.to_numpy())
    frame["meta_ukb_approx_cellwise_fdr_lt_0.05"] = frame["meta_ukb_approx_q_SMR_cell_BH"].lt(0.05)
    heidi = pd.to_numeric(frame["meta_ukb_approx_p_HEIDI"], errors="coerce")
    concordant = frame["meta_ukb_approx_same_direction"].astype(str).str.lower().isin(["true", "1"])
    frame["meta_ukb_approx_validation_pass_cellwise_fdr"] = frame["meta_ukb_approx_cellwise_fdr_lt_0.05"] & concordant & np.isfinite(heidi) & heidi.ge(0.01)
    frame.to_csv(base / f"{stem}.cellwise_fdr.tsv", sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
