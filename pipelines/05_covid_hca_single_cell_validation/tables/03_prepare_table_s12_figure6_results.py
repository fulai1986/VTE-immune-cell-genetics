#!/usr/bin/env python3
"""Prepare the complete long-format source data for Supplementary Table S12."""

import os
from pathlib import Path
import json

import pandas as pd


PROJECT = Path(os.environ["PROJECT_DIR"]).expanduser().resolve()
RESULTS = PROJECT / "covid_validation/stephenson_hca/result5_expression_context/results"
OUT = PROJECT / "manuscript_results/completed_results/result5_covid_single_cell_validation/supplementary"
TMP = PROJECT / "scratch/table_s12"

CD8_FILE = RESULTS / "cd8et_slc44a2_pseudobulk_mean_hallmark_high_low.tsv"
LIANA_FILE = RESULTS / "figure5c_cd8et_monocyte_platelet_bootstrap_contrasts.tsv"
PLATELET_FILE = RESULTS / "platelet_hallmark_high_low_cd8et_slc44a2_28v28.tsv"
CATEGORY_FILE = PROJECT / "config/hallmark_biological_category_map.tsv"

TSV_OUT = OUT / "Table_S12_complete_slc44a2_high_low_results_underlying_figure6.tsv"
JSON_OUT = TMP / "Table_S12_rows.json"

CD8_DISPLAYED = {
    "TNFA_SIGNALING_VIA_NFKB",
    "INFLAMMATORY_RESPONSE",
    "IL6_JAK_STAT3_SIGNALING",
    "COMPLEMENT",
    "COAGULATION",
    "HYPOXIA",
    "REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "APOPTOSIS",
    "P53_PATHWAY",
    "TGF_BETA_SIGNALING",
    "PI3K_AKT_MTOR_SIGNALING",
    "HEME_METABOLISM",
}


def clean_hallmark(value: str) -> str:
    value = value.removeprefix("HALLMARK_")
    replacements = {
        "TNFA_SIGNALING_VIA_NFKB": "TNF-NFkB signaling",
        "IL6_JAK_STAT3_SIGNALING": "IL6-JAK-STAT3 signaling",
        "TGF_BETA_SIGNALING": "TGF-beta signaling",
        "PI3K_AKT_MTOR_SIGNALING": "PI3K-AKT-mTOR signaling",
        "REACTIVE_OXYGEN_SPECIES_PATHWAY": "Reactive oxygen species",
        "P53_PATHWAY": "P53 pathway",
        "MTORC1_SIGNALING": "mTORC1 signaling",
        "MYC_TARGETS_V1": "MYC targets V1",
        "MYC_TARGETS_V2": "MYC targets V2",
        "E2F_TARGETS": "E2F targets",
        "G2M_CHECKPOINT": "G2M checkpoint",
        "KRAS_SIGNALING_UP": "KRAS signaling up",
        "KRAS_SIGNALING_DN": "KRAS signaling down",
        "UV_RESPONSE_UP": "UV response up",
        "UV_RESPONSE_DN": "UV response down",
    }
    return replacements.get(value, value.replace("_", " ").title())


def nullable_number(value):
    return None if pd.isna(value) else float(value)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)

    category = pd.read_csv(CATEGORY_FILE, sep="\t").set_index("hallmark")["category"].to_dict()
    records = []

    cd8 = pd.read_csv(CD8_FILE, sep="\t")
    assert len(cd8) == 50 and cd8["hallmark"].nunique() == 50
    cd8 = cd8.sort_values(["FDR", "p", "hallmark"], kind="stable")
    for _, row in cd8.iterrows():
        raw_name = str(row["hallmark"])
        records.append(
            {
                "Section": "CD8 ET Hallmark programs",
                "Figure panel": "Figure 6C",
                "Outcome class": "Hallmark score",
                "Biological category": category.get(raw_name, ""),
                "Communication axis": "",
                "Outcome / metric": clean_hallmark(raw_name),
                "Source identifier": raw_name,
                "Low donors (n)": int(row["n_low"]),
                "High donors (n)": int(row["n_high"]),
                "Median Low": float(row["median_low"]),
                "Median High": float(row["median_high"]),
                "Delta (High-Low)": float(row["delta_high_minus_low"]),
                "Bootstrap 95% CI lower": None,
                "Bootstrap 95% CI upper": None,
                "Two-sided Mann-Whitney P": float(row["p"]),
                "BH-FDR": float(row["FDR"]),
                "FDR family (tests)": 50,
                "Genes present (n)": int(row["n_genes_present"]),
                "Displayed in Figure 6": "Yes" if raw_name in CD8_DISPLAYED else "No",
            }
        )

    liana = pd.read_csv(LIANA_FILE, sep="\t")
    assert len(liana) == 6 and not liana.duplicated(["axis", "metric"]).any()
    axis_order = {"CD8_ET_to_monocyte": 1, "CD8_ET_to_platelet": 2, "monocyte_to_platelet": 3}
    metric_order = {"n_cellphone_p_lt_0.05": 1, "sum_rank_strength_sig": 2}
    liana = liana.assign(
        axis_order=liana["axis"].map(axis_order), metric_order=liana["metric"].map(metric_order)
    ).sort_values(["axis_order", "metric_order"])
    for _, row in liana.iterrows():
        records.append(
            {
                "Section": "Donor-level LIANA communication",
                "Figure panel": "Figure 6D-E",
                "Outcome class": "Communication metric",
                "Biological category": "",
                "Communication axis": str(row["axis_label"]).replace("->", "→"),
                "Outcome / metric": str(row["metric_label"]).replace("Significant LR pairs", "Number of significant LR pairs").replace("Summed rank evidence", "Summed LIANA rank evidence"),
                "Source identifier": str(row["metric"]),
                "Low donors (n)": int(row["n_low"]),
                "High donors (n)": int(row["n_high"]),
                "Median Low": float(row["median_low"]),
                "Median High": float(row["median_high"]),
                "Delta (High-Low)": float(row["delta_high_minus_low"]),
                "Bootstrap 95% CI lower": nullable_number(row["ci_low"]),
                "Bootstrap 95% CI upper": nullable_number(row["ci_high"]),
                "Two-sided Mann-Whitney P": float(row["p"]),
                "BH-FDR": float(row["FDR"]),
                "FDR family (tests)": 6,
                "Genes present (n)": None,
                "Displayed in Figure 6": "Yes",
            }
        )

    platelet = pd.read_csv(PLATELET_FILE, sep="\t")
    assert len(platelet) == 9 and platelet["hallmark"].nunique() == 9
    platelet = platelet.sort_values(["FDR", "p", "hallmark"], kind="stable")
    for _, row in platelet.iterrows():
        source_name = str(row["hallmark"])
        raw_name = source_name.removeprefix("HALLMARK_")
        records.append(
            {
                "Section": "Platelet Hallmark programs",
                "Figure panel": "Figure 6F",
                "Outcome class": "Hallmark score",
                "Biological category": category.get(raw_name, ""),
                "Communication axis": "",
                "Outcome / metric": clean_hallmark(source_name),
                "Source identifier": source_name,
                "Low donors (n)": int(row["n_low"]),
                "High donors (n)": int(row["n_high"]),
                "Median Low": float(row["median_low"]),
                "Median High": float(row["median_high"]),
                "Delta (High-Low)": float(row["delta_high_minus_low"]),
                "Bootstrap 95% CI lower": None,
                "Bootstrap 95% CI upper": None,
                "Two-sided Mann-Whitney P": float(row["p"]),
                "BH-FDR": float(row["FDR"]),
                "FDR family (tests)": 9,
                "Genes present (n)": int(row["n_genes_present"]),
                "Displayed in Figure 6": "Yes",
            }
        )

    out = pd.DataFrame(records)
    assert len(out) == 65
    assert out["Section"].value_counts().to_dict() == {
        "CD8 ET Hallmark programs": 50,
        "Platelet Hallmark programs": 9,
        "Donor-level LIANA communication": 6,
    }
    out.to_csv(TSV_OUT, sep="\t", index=False, na_rep="")
    json_records = json.loads(out.to_json(orient="records"))
    JSON_OUT.write_text(json.dumps(json_records, indent=2), encoding="utf-8")
    print(out.groupby("Section").size().to_string())
    print(f"Wrote {TSV_OUT}")
    print(f"Wrote {JSON_OUT}")


if __name__ == "__main__":
    main()
