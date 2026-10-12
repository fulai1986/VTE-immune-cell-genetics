#!/usr/bin/env python3
"""Export the fixed all-cell annotation and cluster-level COVID cell counts."""
import os
from pathlib import Path

import pandas as pd


def main():
    base = Path(os.environ["PROJECT_DIR"]).expanduser().resolve() / "covid_validation/hca_onek1k_mapping/results"
    frame = pd.read_csv(base / "hca_obs_with_onek1k_compatible_annotation.curated.tsv.gz", sep="\t")
    if frame["hca_cell_id"].duplicated().any():
        raise ValueError("Duplicate HCA cell identifiers")
    # Preserve H5AD observation order for downstream matrix indexing.
    frame.to_csv(base / "hca_all_cells_onek1k_compatible_annotation.curated.tsv.gz",
                 sep="\t", index=False, compression="gzip")
    counts = (frame.loc[frame["Status"].eq("Covid")]
              .groupby(["full_clustering", "curated_OneK1K_compatible", "curated_mapping_confidence"], dropna=False)
              .size().reset_index(name="covid_cells")
              .sort_values("covid_cells", ascending=False))
    counts.to_csv(base / "hca_covid_cluster_onek1k_compatible_counts.tsv", sep="\t", index=False)
    print(f"Exported {len(frame)} cells and {len(counts)} COVID clusters")


if __name__ == "__main__":
    main()
