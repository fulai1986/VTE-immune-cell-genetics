#!/usr/bin/env python3
"""Create an analysis workspace and copy editable configuration templates."""
import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path)
    args = parser.parse_args()
    workspace = args.workspace.expanduser().resolve()
    if workspace == REPO or REPO in workspace.parents:
        parser.error("Choose a workspace outside the source repository.")
    for name in ["config", "metadata", "results", "logs", "reports", "scratch", "reference/onek1k_official/esnp", "resources/msigdb", "gwas/config", "gwas/download", "gwas/replication/finngen_r13", "covid_validation/stephenson_hca"]:
        (workspace / name).mkdir(parents=True, exist_ok=True)
    templates = {
        "config/hallmark_biological_category_map.tsv": "config/hallmark_biological_category_map.tsv",
        "pipelines/01_eqtl_pipeline/config/project_paths.example.tsv": "config/project_paths.tsv",
        "pipelines/02_smr_coloc_pipeline/config/science14_celltypes.tsv": "config/science14_celltypes.tsv",
        "pipelines/03_meta_gwas_qc_pipeline/config/gwas_dataset_manifest.example.tsv": "gwas/config/gwas_dataset_manifest.tsv",
    }
    for source, target in templates.items():
        destination = workspace / target
        if destination.exists():
            print(f"Keep existing {destination}")
            continue
        with (REPO / source).open(newline="") as handle:
            rows = [row for row in csv.reader(handle, delimiter="\t") if row]
        for row in rows[1:]:
            if row[0] == "project_dir":
                row[1] = str(workspace)
        with destination.open("w", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
        print(destination)


if __name__ == "__main__":
    main()
