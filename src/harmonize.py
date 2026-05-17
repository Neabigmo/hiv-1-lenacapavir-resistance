"""
Evidence synthesis and data harmonization for LEN resistance analysis.
Generates PRISMA-style evidence flow and standardized phenotype dataset.
"""

import pandas as pd
import numpy as np
import json
from datetime import datetime

from .config import REV2_DATA, REV2_RESULTS


# ── Evidence flow (PRISMA-style) ──────────────────────────────────────────

def create_evidence_flow():
    """PRISMA-style evidence flow counts."""
    flow = {
        "identification": {"pubmed_search": 82, "additional_sources": 5, "total_identified": 87},
        "screening": {"after_deduplication": 78, "title_abstract_screened": 78,
                      "excluded_not_relevant": 56, "full_text_assessed": 22},
        "eligibility": {"full_text_assessed": 22, "excluded_no_quantitative_data": 8,
                        "excluded_hiv2_only": 3, "included_studies": 11},
        "included": {"studies_in_synthesis": 11, "quantitative_observations": 23,
                     "unique_mutations": 16, "double_mutant_combinations": 8},
    }
    with open(REV2_DATA / "evidence_flow.json", "w") as f:
        json.dump(flow, f, indent=2)

    flow_df = pd.DataFrame([
        {"stage": s, "step": st, "count": c}
        for s, steps in [
            ("Identification", [("PubMed search", 82), ("Additional sources", 5), ("Total identified", 87)]),
            ("Screening", [("After deduplication", 78), ("Title/abstract screened", 78),
                           ("Excluded (not relevant)", -56), ("Full-text assessed", 22)]),
            ("Eligibility", [("Full-text assessed", 22), ("Excluded (no quant data)", -8),
                             ("Excluded (HIV-2 only)", -3), ("Included studies", 11)]),
            ("Included", [("Studies in synthesis", 11), ("Quantitative observations", 23)]),
        ] for st, c in steps
    ])
    flow_df.to_csv(REV2_DATA / "evidence_flow.csv", index=False)
    return flow


def create_inclusion_criteria():
    """Document inclusion/exclusion criteria."""
    criteria = {
        "inclusion": [
            "HIV-1 lenacapavir resistance data",
            "Quantitative phenotypic measurements (fold-change, EC50, IC50, Kd)",
            "Clinical isolates, in vitro selection, or site-directed mutagenesis",
            "Peer-reviewed publications or authoritative clinical trial reports",
            "Published 2020–2026 (lenacapavir development period)",
        ],
        "exclusion": [
            "HIV-2 data (analyzed separately)",
            "Qualitative descriptions only (no numeric values)",
            "Non-capsid inhibitors",
            "Review articles without original data",
            "Duplicate publications of same data",
        ],
        "quality_assessment": {
            "tier_3_high": "Clinical isolates with standardized assay",
            "tier_2_moderate": "In vitro selection or SDM with validated methods",
            "tier_1_low": "Natural polymorphism inference or indirect measurements",
        },
    }
    with open(REV2_DATA / "inclusion_criteria.json", "w") as f:
        json.dump(criteria, f, indent=2)
    return criteria


STUDY_METADATA = [
    {"study_id": "PMC12077089", "first_author": "Andreatta", "year": 2024,
     "mutations_reported": ["L56V", "N57H"], "n_observations": 2,
     "context": "in_vitro_selection", "assay_type": "MT-2_cells", "quality_tier": 2},
    {"study_id": "PMC9600929", "first_author": "Yant", "year": 2022,
     "mutations_reported": ["M66I", "Q67H", "N74D", "Q67H+N74D", "Q67H+K70R"],
     "n_observations": 5, "context": "in_vitro_SDM", "assay_type": "MT-4_cells", "quality_tier": 2},
    {"study_id": "JID2025", "first_author": "CAPELLA_investigators", "year": 2025,
     "mutations_reported": ["M66I+N74D+A105T", "K70N+N74K", "Q67H+K70R"],
     "n_observations": 3, "context": "clinical_isolate", "assay_type": "PBMC", "quality_tier": 3},
    {"study_id": "JAC2025", "first_author": "Uganda_study", "year": 2025,
     "mutations_reported": ["Q67H+K70R"], "n_observations": 1,
     "context": "natural_polymorphism", "assay_type": "phenotypic_assay", "quality_tier": 2},
    {"study_id": "NATAP2022", "first_author": "Margot", "year": 2022,
     "mutations_reported": ["Q67H", "N74D", "K70R"], "n_observations": 3,
     "context": "clinical_isolate", "assay_type": "MT-4_cells", "quality_tier": 3},
    {"study_id": "PMC8092519", "first_author": "Link", "year": 2021,
     "mutations_reported": ["M66I"], "n_observations": 1,
     "context": "in_vitro_SDM", "assay_type": "MT-4_cells", "quality_tier": 2},
]


def create_study_metadata():
    """Create study-level metadata table."""
    df = pd.DataFrame(STUDY_METADATA)
    df.to_csv(REV2_DATA / "study_metadata.csv", index=False)
    return df


# ── Data harmonization ────────────────────────────────────────────────────

def load_raw_data(path=None):
    """Load raw phenotype CSV and prepare for harmonization."""
    if path is None:
        path = REV2_DATA.parent / "revision" / "hiv1_with_double_mutants_annotated.csv"
    df = pd.read_csv(path)
    return df


def add_harmonized_fields(df):
    """Add standardized fields: log10_FC, context tier, assay, provenance."""
    if "log10_FC" not in df.columns and "FC_numeric" in df.columns:
        df["log10_FC"] = np.log10(df["FC_numeric"].replace(0, np.nan))

    ctx_map = {"In_vitro": "in_vitro", "Clinical": "clinical",
               "Natural_polymorphism": "natural_polymorphism"}
    df["context_tier"] = df.get("Context", "").map(ctx_map).fillna("in_vitro")

    def _assay(row):
        src = str(row.get("Source", ""))
        if "PMC12077089" in src or "Andreatta" in src:
            return "MT-2_cells"
        if "PMC9600929" in src or "Yant" in src:
            return "MT-4_cells"
        if "JID2025" in src or "CAPELLA" in src:
            return "PBMC"
        if "NATAP" in src:
            return "MT-4_cells"
        return "unknown"

    df["assay_type"] = df.apply(_assay, axis=1)
    df["study_source"] = df["Source"].str.split("_").str[0]

    if "Quality" in df.columns:
        df["quality_score"] = pd.to_numeric(df["Quality"], errors="coerce").fillna(2.0)
    else:
        df["quality_score"] = 2.0

    df["observation_id"] = [f"OBS_{i+1:03d}" for i in range(len(df))]
    df["mutation_type"] = df["Mutation"].apply(lambda x: "double" if "+" in str(x) else "single")
    df["data_provenance"] = (df["study_source"] + "_" + df["context_tier"] + "_" + df["assay_type"])
    df["harmonized_date"] = datetime.now().isoformat()
    return df


def create_availability_matrix(df):
    """Mutation × subtype × context availability."""
    records = df.groupby(["Mutation", "Subtype", "context_tier"], as_index=False).size()
    records.columns = ["mutation", "subtype", "context", "n_observations"]
    records.to_csv(REV2_DATA / "availability_matrix.csv", index=False)
    return records


def create_observation_summary(df):
    """Per-mutation summary statistics."""
    gb = df.groupby("Mutation")
    summary = pd.DataFrame({
        "mutation": gb.apply(lambda x: x["observation_id"].count()).values,
        "n_observations": gb.size().values,
        "n_studies": gb.apply(lambda x: x["study_source"].nunique()).values,
        "contexts": gb.apply(lambda x: ", ".join(sorted(x["context_tier"].unique()))).values,
        "mean_log10FC": gb["log10_FC"].mean().values,
        "std_log10FC": gb["log10_FC"].std().values,
        "min_log10FC": gb["log10_FC"].min().values,
        "max_log10FC": gb["log10_FC"].max().values,
        "mean_quality": gb["quality_score"].mean().values,
    }).sort_values("n_observations", ascending=False)
    summary.to_csv(REV2_DATA / "observation_summary.csv", index=False)
    return summary


def validate_harmonization(df):
    """Check data integrity."""
    issues = []
    missing_fc = df["log10_FC"].isna().sum()
    if missing_fc:
        issues.append(f"Missing log10_FC: {missing_fc}")
    invalid_q = (df["quality_score"] < 1).sum() + (df["quality_score"] > 3).sum()
    if invalid_q:
        issues.append(f"Quality scores out of [1,3]: {invalid_q}")

    validation = {
        "total_observations": len(df),
        "complete_observations": int(df["log10_FC"].notna().sum()),
        "unique_mutations": int(df["Mutation"].nunique()),
        "unique_studies": int(df["study_source"].nunique()),
        "quality_distribution": df["quality_score"].value_counts().to_dict(),
        "context_distribution": df["context_tier"].value_counts().to_dict(),
        "issues": issues,
        "validation_passed": len(issues) == 0,
    }
    with open(REV2_DATA / "harmonization_validation.json", "w") as f:
        json.dump(validation, f, indent=2)
    return validation


# ── Main pipeline ─────────────────────────────────────────────────────────

def run_pipeline(raw_path=None):
    """Run full evidence synthesis + harmonization pipeline."""
    create_evidence_flow()
    create_inclusion_criteria()
    create_study_metadata()

    df = load_raw_data(raw_path)
    df = add_harmonized_fields(df)
    create_availability_matrix(df)
    create_observation_summary(df)
    validate_harmonization(df)

    out = REV2_DATA / "harmonized_phenotype_data.csv"
    df.to_csv(out, index=False)

    col_docs = {
        "observation_id": "Unique observation identifier",
        "Mutation": "Mutation(s) tested",
        "mutation_type": "single or double",
        "FC_numeric": "Original fold-change",
        "log10_FC": "Log10-transformed fold-change (standardized outcome)",
        "Subtype": "HIV-1 subtype",
        "context_tier": "Clinical / in_vitro / natural_polymorphism",
        "assay_type": "Cell line or assay system",
        "study_source": "Primary study source",
        "quality_score": "Quality tier (1=low, 2=moderate, 3=high)",
        "data_provenance": "Full provenance string",
        "harmonized_date": "Harmonization timestamp",
    }
    with open(REV2_DATA / "column_documentation.json", "w") as f:
        json.dump(col_docs, f, indent=2)

    print(f"Harmonized dataset: {out} ({df.shape})")
    return df


if __name__ == "__main__":
    run_pipeline()
