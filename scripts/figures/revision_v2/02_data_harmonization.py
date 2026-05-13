#!/usr/bin/env python3
"""
Data Harmonization - Create fully annotated phenotype dataset
Adds provenance, quality tiers, and standardized fields
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime

# Setup paths
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "data"
INPUT_FILE = DATA_DIR / "processed" / "revision" / "hiv1_with_double_mutants_annotated.csv"
OUTPUT_DIR = DATA_DIR / "processed" / "revision_v2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_and_clean_data():
    """Load existing data and prepare for harmonization"""
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} observations from {INPUT_FILE.name}")
    print(f"Columns: {df.columns.tolist()}")
    return df

def add_harmonized_fields(df):
    """Add standardized fields for evidence synthesis"""

    # Ensure log10_FC exists
    if 'log10_FC' not in df.columns and 'FC_numeric' in df.columns:
        df['log10_FC'] = np.log10(df['FC_numeric'].replace(0, np.nan))

    # Map context to standardized tiers
    context_mapping = {
        'In_vitro': 'in_vitro',
        'Clinical': 'clinical',
        'Natural_polymorphism': 'natural_polymorphism'
    }
    df['context_tier'] = df['Context'].map(context_mapping).fillna('in_vitro')

    # Assay type from source
    def infer_assay_type(row):
        source = str(row.get('Source', ''))
        if 'PMC12077089' in source or 'Andreatta' in source:
            return 'MT-2_cells'
        elif 'PMC9600929' in source or 'Yant' in source:
            return 'MT-4_cells'
        elif 'JID2025' in source or 'CAPELLA' in source:
            return 'PBMC'
        elif 'NATAP' in source:
            return 'MT-4_cells'
        else:
            return 'unknown'

    df['assay_type'] = df.apply(infer_assay_type, axis=1)

    # Study source (clean)
    df['study_source'] = df['Source'].str.split('_').str[0]

    # Quality score (already exists, ensure it's numeric)
    if 'Quality' in df.columns:
        df['quality_score'] = pd.to_numeric(df['Quality'], errors='coerce').fillna(2.0)
    else:
        df['quality_score'] = 2.0  # Default moderate quality

    # Observation ID
    df['observation_id'] = [f"OBS_{i+1:03d}" for i in range(len(df))]

    # Mutation type (single vs double)
    df['mutation_type'] = df['Mutation'].apply(
        lambda x: 'double' if '+' in str(x) else 'single'
    )

    # Data provenance
    df['data_provenance'] = df.apply(
        lambda row: f"{row['study_source']}_{row['context_tier']}_{row['assay_type']}",
        axis=1
    )

    # Harmonization timestamp
    df['harmonized_date'] = datetime.now().isoformat()

    return df

def create_availability_matrix(df):
    """Create mutation × subtype × context availability matrix"""

    # Get unique values
    mutations = df['Mutation'].unique()
    subtypes = df['Subtype'].unique()
    contexts = df['context_tier'].unique()

    # Create matrix
    availability = []
    for mut in mutations:
        for subtype in subtypes:
            for context in contexts:
                count = len(df[
                    (df['Mutation'] == mut) &
                    (df['Subtype'] == subtype) &
                    (df['context_tier'] == context)
                ])
                if count > 0:
                    availability.append({
                        'mutation': mut,
                        'subtype': subtype,
                        'context': context,
                        'n_observations': count
                    })

    avail_df = pd.DataFrame(availability)
    avail_df.to_csv(OUTPUT_DIR / "availability_matrix.csv", index=False)
    print(f"[OK] Availability matrix: {len(avail_df)} mutation-subtype-context combinations")

    return avail_df

def create_observation_summary(df):
    """Create per-mutation observation summary"""

    summary = df.groupby('Mutation').agg({
        'observation_id': 'count',
        'study_source': lambda x: len(x.unique()),
        'context_tier': lambda x: ', '.join(sorted(x.unique())),
        'log10_FC': ['mean', 'std', 'min', 'max'],
        'quality_score': 'mean'
    }).reset_index()

    summary.columns = [
        'mutation', 'n_observations', 'n_studies', 'contexts',
        'mean_log10FC', 'std_log10FC', 'min_log10FC', 'max_log10FC',
        'mean_quality'
    ]

    summary = summary.sort_values('n_observations', ascending=False)
    summary.to_csv(OUTPUT_DIR / "observation_summary.csv", index=False)
    print(f"[OK] Observation summary: {len(summary)} unique mutations")

    return summary

def validate_harmonization(df):
    """Validate harmonized dataset"""

    issues = []

    # Check for missing log10_FC
    missing_fc = df['log10_FC'].isna().sum()
    if missing_fc > 0:
        issues.append(f"Missing log10_FC: {missing_fc} observations")

    # Check for missing provenance
    missing_prov = df['data_provenance'].isna().sum()
    if missing_prov > 0:
        issues.append(f"Missing provenance: {missing_prov} observations")

    # Check quality scores
    invalid_quality = df[~df['quality_score'].between(1, 3)].shape[0]
    if invalid_quality > 0:
        issues.append(f"Invalid quality scores: {invalid_quality} observations")

    validation = {
        "total_observations": len(df),
        "complete_observations": len(df[df['log10_FC'].notna()]),
        "unique_mutations": df['Mutation'].nunique(),
        "unique_studies": df['study_source'].nunique(),
        "quality_distribution": df['quality_score'].value_counts().to_dict(),
        "context_distribution": df['context_tier'].value_counts().to_dict(),
        "issues": issues,
        "validation_passed": len(issues) == 0
    }

    with open(OUTPUT_DIR / "harmonization_validation.json", "w") as f:
        json.dump(validation, f, indent=2)

    print(f"\n[OK] Validation complete:")
    print(f"  Total observations: {validation['total_observations']}")
    print(f"  Complete observations: {validation['complete_observations']}")
    print(f"  Unique mutations: {validation['unique_mutations']}")
    print(f"  Issues found: {len(issues)}")
    if issues:
        for issue in issues:
            print(f"    - {issue}")

    return validation

def main():
    """Main execution"""
    print("="*60)
    print("Data Harmonization - Phenotype Dataset")
    print("="*60)

    # Load data
    df = load_and_clean_data()

    # Add harmonized fields
    print("\nAdding harmonized fields...")
    df_harmonized = add_harmonized_fields(df)

    # Create availability matrix
    print("\nCreating availability matrix...")
    avail = create_availability_matrix(df_harmonized)

    # Create observation summary
    print("\nCreating observation summary...")
    summary = create_observation_summary(df_harmonized)

    # Validate
    print("\nValidating harmonization...")
    validation = validate_harmonization(df_harmonized)

    # Save harmonized dataset
    output_file = OUTPUT_DIR / "harmonized_phenotype_data.csv"
    df_harmonized.to_csv(output_file, index=False)
    print(f"\n[OK] Harmonized dataset saved: {output_file}")
    print(f"  Shape: {df_harmonized.shape}")
    print(f"  Columns: {len(df_harmonized.columns)}")

    # Save column documentation
    column_docs = {
        "observation_id": "Unique identifier for each observation",
        "Mutation": "Mutation(s) tested (single or combination)",
        "mutation_type": "single or double",
        "FC_numeric": "Original fold-change value",
        "log10_FC": "Log10-transformed fold-change (standardized outcome)",
        "Subtype": "HIV-1 subtype",
        "context_tier": "Experimental context: clinical/in_vitro/natural_polymorphism",
        "assay_type": "Cell line or assay system used",
        "study_source": "Primary study identifier",
        "quality_score": "Quality tier (1=low, 2=moderate, 3=high)",
        "data_provenance": "Full provenance string",
        "harmonized_date": "Date of harmonization"
    }

    with open(OUTPUT_DIR / "column_documentation.json", "w") as f:
        json.dump(column_docs, f, indent=2)

    print("\n" + "="*60)
    print("Data harmonization complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*60)

if __name__ == "__main__":
    main()
