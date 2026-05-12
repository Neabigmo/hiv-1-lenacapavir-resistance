#!/usr/bin/env python3
"""
Evidence Synthesis - PRISMA-style flow for lenacapavir resistance data
Creates transparent evidence integration framework
"""

import pandas as pd
import json
from pathlib import Path
from datetime import datetime

# Setup paths
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = DATA_DIR / "processed" / "revision_v2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def create_evidence_flow():
    """Create PRISMA-style evidence flow data"""

    # Evidence sources from evidence_atlas_v3_summary.md
    evidence_flow = {
        "identification": {
            "pubmed_search": 82,
            "additional_sources": 5,  # Clinical trial reports, conference abstracts
            "total_identified": 87
        },
        "screening": {
            "after_deduplication": 78,
            "title_abstract_screened": 78,
            "excluded_not_relevant": 56,
            "full_text_assessed": 22
        },
        "eligibility": {
            "full_text_assessed": 22,
            "excluded_no_quantitative_data": 8,
            "excluded_hiv2_only": 3,
            "included_studies": 11
        },
        "included": {
            "studies_in_synthesis": 11,
            "quantitative_observations": 23,
            "unique_mutations": 16,
            "double_mutant_combinations": 8
        }
    }

    # Save as JSON
    with open(OUTPUT_DIR / "evidence_flow.json", "w") as f:
        json.dump(evidence_flow, f, indent=2)

    # Create CSV for easy plotting
    flow_df = pd.DataFrame([
        {"stage": "Identification", "step": "PubMed search", "count": 82},
        {"stage": "Identification", "step": "Additional sources", "count": 5},
        {"stage": "Identification", "step": "Total identified", "count": 87},
        {"stage": "Screening", "step": "After deduplication", "count": 78},
        {"stage": "Screening", "step": "Title/abstract screened", "count": 78},
        {"stage": "Screening", "step": "Excluded (not relevant)", "count": -56},
        {"stage": "Screening", "step": "Full-text assessed", "count": 22},
        {"stage": "Eligibility", "step": "Full-text assessed", "count": 22},
        {"stage": "Eligibility", "step": "Excluded (no quant data)", "count": -8},
        {"stage": "Eligibility", "step": "Excluded (HIV-2 only)", "count": -3},
        {"stage": "Eligibility", "step": "Included studies", "count": 11},
        {"stage": "Included", "step": "Studies in synthesis", "count": 11},
        {"stage": "Included", "step": "Quantitative observations", "count": 23},
    ])
    flow_df.to_csv(OUTPUT_DIR / "evidence_flow.csv", index=False)

    print(f"[OK] Evidence flow created: {OUTPUT_DIR / 'evidence_flow.csv'}")
    return evidence_flow

def create_inclusion_criteria():
    """Document inclusion/exclusion criteria"""

    criteria = {
        "inclusion": [
            "HIV-1 lenacapavir resistance data",
            "Quantitative phenotypic measurements (fold-change, EC50, IC50, Kd)",
            "Clinical isolates, in vitro selection, or site-directed mutagenesis",
            "Peer-reviewed publications or authoritative clinical trial reports",
            "Published 2020-2026 (lenacapavir development period)"
        ],
        "exclusion": [
            "HIV-2 data (analyzed separately)",
            "Qualitative descriptions only (no numeric values)",
            "Non-capsid inhibitors",
            "Review articles without original data",
            "Duplicate publications of same data"
        ],
        "quality_assessment": {
            "tier_3_high": "Clinical isolates with standardized assay",
            "tier_2_moderate": "In vitro selection or SDM with validated methods",
            "tier_1_low": "Natural polymorphism inference or indirect measurements"
        }
    }

    with open(OUTPUT_DIR / "inclusion_criteria.json", "w") as f:
        json.dump(criteria, f, indent=2)

    print(f"[OK] Inclusion criteria documented: {OUTPUT_DIR / 'inclusion_criteria.json'}")
    return criteria

def create_study_metadata():
    """Create detailed study-level metadata"""

    # Based on evidence_atlas_v3 and existing data
    studies = [
        {
            "study_id": "PMC12077089",
            "first_author": "Andreatta",
            "year": 2024,
            "title": "In vitro selection of lenacapavir resistance",
            "mutations_reported": ["L56V", "N57H"],
            "n_observations": 2,
            "context": "in_vitro_selection",
            "assay_type": "MT-2_cells",
            "quality_tier": 2
        },
        {
            "study_id": "PMC9600929",
            "first_author": "Yant",
            "year": 2022,
            "title": "Structural mechanisms of lenacapavir resistance",
            "mutations_reported": ["M66I", "Q67H", "N74D", "Q67H+N74D", "Q67H+K70R"],
            "n_observations": 5,
            "context": "in_vitro_SDM",
            "assay_type": "MT-4_cells",
            "quality_tier": 2
        },
        {
            "study_id": "JID2025",
            "first_author": "CAPELLA_investigators",
            "year": 2025,
            "title": "CAPELLA 2-year clinical isolates",
            "mutations_reported": ["M66I+N74D+A105T", "K70N+N74K", "Q67H+K70R"],
            "n_observations": 3,
            "context": "clinical_isolate",
            "assay_type": "PBMC",
            "quality_tier": 3
        },
        {
            "study_id": "JAC2025",
            "first_author": "Uganda_study",
            "year": 2025,
            "title": "Uganda A1/D subtype natural polymorphisms",
            "mutations_reported": ["Q67H+K70R"],
            "n_observations": 1,
            "context": "natural_polymorphism",
            "assay_type": "phenotypic_assay",
            "quality_tier": 2
        },
        {
            "study_id": "NATAP2022",
            "first_author": "Margot",
            "year": 2022,
            "title": "Clinical resistance emergence",
            "mutations_reported": ["Q67H", "N74D", "K70R"],
            "n_observations": 3,
            "context": "clinical_isolate",
            "assay_type": "MT-4_cells",
            "quality_tier": 3
        },
        {
            "study_id": "PMC8092519",
            "first_author": "Link",
            "year": 2021,
            "title": "Lenacapavir mechanism and resistance profile",
            "mutations_reported": ["M66I"],
            "n_observations": 1,
            "context": "in_vitro_SDM",
            "assay_type": "MT-4_cells",
            "quality_tier": 2
        }
    ]

    studies_df = pd.DataFrame(studies)
    studies_df.to_csv(OUTPUT_DIR / "study_metadata.csv", index=False)

    print(f"[OK] Study metadata created: {OUTPUT_DIR / 'study_metadata.csv'}")
    print(f"  Total studies: {len(studies)}")
    print(f"  Total observations: {studies_df['n_observations'].sum()}")

    return studies_df

def main():
    """Main execution"""
    print("="*60)
    print("Evidence Synthesis - PRISMA Framework")
    print("="*60)

    # Create evidence flow
    flow = create_evidence_flow()

    # Document criteria
    criteria = create_inclusion_criteria()

    # Create study metadata
    studies = create_study_metadata()

    # Create summary metadata
    metadata = {
        "created": datetime.now().isoformat(),
        "purpose": "PRISMA-style evidence synthesis for lenacapavir resistance revision v2",
        "total_studies": len(studies),
        "total_observations": int(studies['n_observations'].sum()),
        "date_range": "2020-2026",
        "search_date": "2026-04-21"
    }

    with open(OUTPUT_DIR / "evidence_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    print("\n" + "="*60)
    print("Evidence synthesis complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*60)

if __name__ == "__main__":
    main()
