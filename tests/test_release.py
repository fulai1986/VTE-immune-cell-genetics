import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
PIPE = ROOT / "pipelines/05_covid_hca_single_cell_validation"


class ReleaseTests(unittest.TestCase):
    def test_workspace_initialization_preserves_local_configuration(self):
        with tempfile.TemporaryDirectory(prefix="vte workspace ") as directory:
            command = [sys.executable, str(ROOT / "tools/init_workspace.py"), directory]
            subprocess.run(command, check=True, capture_output=True)
            config = Path(directory) / "config/project_paths.tsv"
            data = pd.read_csv(config, sep="\t")
            self.assertEqual(data.loc[data.key.eq("project_dir"), "path"].iloc[0], str(Path(directory).resolve()))
            content = config.read_text()
            subprocess.run(command, check=True, capture_output=True)
            self.assertEqual(config.read_text(), content)
            subprocess.run([sys.executable, str(ROOT / "run.py"), "run", "hca-mapping", "--project", directory, "--dry-run"], check=True, capture_output=True)

    def test_bh_and_bootstrap(self):
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {"PROJECT_DIR": directory}):
            scope = runpy.run_path(str(PIPE / "analysis/07_summarize_liana_high_low_contrasts.py"))
            values = pd.Series([0.001, 0.04, 0.2, 0.009, 0.9, 0.07])
            np.testing.assert_allclose(scope["bh_fdr"](values), false_discovery_control(values.to_numpy()))
            ci = scope["bootstrap_median_difference"](np.ones(10), np.full(10, 3.0), 12)
            np.testing.assert_allclose(ci, [2, 2])

    def test_tertiles_and_workspace_relocation(self):
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {"PROJECT_DIR": directory}):
            scope = runpy.run_path(str(PIPE / "analysis/06_build_high_low_hallmark_contrasts.py"))
            path = scope["EXPOSURE"]
            path.parent.mkdir(parents=True)
            pd.DataFrame({"cell": ["CD8_ET"] * 84, "Gene": ["SLC44A2"] * 84,
                          "sample_id": [f"sample{i}" for i in range(84)], "candidate_mean_expr": np.arange(84)}).to_csv(path, sep="\t", index=False)
            counts = scope["exposure_groups"]()["group"].value_counts().to_dict()
            self.assertEqual(counts, {"Low": 28, "Middle": 28, "High": 28})

    def test_fixed_mapping_export_preserves_cell_order(self):
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {"PROJECT_DIR": directory}):
            base = Path(directory) / "covid_validation/hca_onek1k_mapping/results"
            base.mkdir(parents=True)
            frame = pd.DataFrame({"hca_cell_id": ["z", "a", "m"], "Status": ["Covid", "Healthy", "Covid"],
                                  "full_clustering": ["CD8.TE"] * 3, "curated_OneK1K_compatible": ["CD8 ET"] * 3,
                                  "curated_mapping_confidence": ["High"] * 3})
            frame.to_csv(base / "hca_obs_with_onek1k_compatible_annotation.curated.tsv.gz", sep="\t", index=False)
            runpy.run_path(str(PIPE / "mapping/06_export_fixed_mapping_and_covid_counts.py"))["main"]()
            result = pd.read_csv(base / "hca_all_cells_onek1k_compatible_annotation.curated.tsv.gz", sep="\t")
            self.assertEqual(result.hca_cell_id.tolist(), ["z", "a", "m"])
            counts = pd.read_csv(base / "hca_covid_cluster_onek1k_compatible_counts.tsv", sep="\t")
            self.assertEqual(counts.covid_cells.tolist(), [2])

    def test_discovery_uses_filesystem_cell_identifiers(self):
        stages = json.loads((ROOT / "config/workflows.json").read_text())
        commands = [x for x in stages["finngen-discovery"] if "--cell" in x]
        self.assertEqual(len(commands), 14)
        self.assertTrue(all(" " not in x[x.index("--cell") + 1] for x in commands))
        self.assertIn("--allele-specific-ids", stages["ld-reference"][0])


if __name__ == "__main__":
    unittest.main()
