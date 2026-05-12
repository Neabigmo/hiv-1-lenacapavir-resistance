#!/usr/bin/env python3
"""
Epistasis and Compensatory Mutation Analysis
Proper baseline calculation and context-dependent analysis
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json

# Setup paths
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "results" / "revision_v2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load harmonized data"""
    df = pd.read_csv(DATA_DIR / "harmonized_phenotype_data.csv")
    df_complete = df[df['log10_FC'].notna()].copy()
    print(f"Loaded {len(df_complete)} complete observations")
    return df_complete

def identify_double_mutants(df):
    """Identify double mutant combinations"""
    df['is_double'] = df['Mutation'].str.contains(r'\+', regex=True)
    double_mutants = df[df['is_double']].copy()

    print(f"\nFound {len(double_mutants)} double mutant observations:")
    for _, row in double_mutants.iterrows():
        print(f"  {row['Mutation']}: {row['FC_numeric']:.1f}-fold (log10={row['log10_FC']:.2f})")

    return double_mutants

def calculate_epistasis(df, double_mutants):
    """Calculate epistasis with proper baseline"""

    epistasis_results = []

    for _, dm_row in double_mutants.iterrows():
        combination = dm_row['Mutation']
        observed_log10fc = dm_row['log10_FC']
        observed_fc = dm_row['FC_numeric']

        # Parse combination
        mutations = combination.split('+')
        if len(mutations) != 2:
            continue

        mut1, mut2 = mutations

        # Find single mutant effects
        mut1_data = df[df['Mutation'] == mut1]
        mut2_data = df[df['Mutation'] == mut2]

        if len(mut1_data) == 0 or len(mut2_data) == 0:
            # Single mutants not measured
            epistasis_results.append({
                'combination': combination,
                'mut1': mut1,
                'mut2': mut2,
                'observed_log10fc': observed_log10fc,
                'observed_fc': observed_fc,
                'expected_log10fc': np.nan,
                'expected_fc': np.nan,
                'interaction_residual': np.nan,
                'interaction_type': 'unknown_singles',
                'note': 'Single mutants not measured'
            })
            continue

        # Calculate expected (additive in log space)
        mut1_log10fc = mut1_data['log10_FC'].mean()
        mut2_log10fc = mut2_data['log10_FC'].mean()

        expected_log10fc = mut1_log10fc + mut2_log10fc
        expected_fc = 10 ** expected_log10fc

        # Interaction residual
        interaction_residual = observed_log10fc - expected_log10fc

        # Classify interaction
        if abs(interaction_residual) < 0.3:  # <2-fold difference
            interaction_type = 'additive'
        elif interaction_residual > 0.3:
            interaction_type = 'positive_synergy'
        else:
            interaction_type = 'negative_synergy'

        epistasis_results.append({
            'combination': combination,
            'mut1': mut1,
            'mut2': mut2,
            'observed_log10fc': observed_log10fc,
            'observed_fc': observed_fc,
            'mut1_log10fc': mut1_log10fc,
            'mut2_log10fc': mut2_log10fc,
            'expected_log10fc': expected_log10fc,
            'expected_fc': expected_fc,
            'interaction_residual': interaction_residual,
            'fold_amplification': observed_fc / expected_fc if expected_fc > 0 else np.nan,
            'interaction_type': interaction_type
        })

    epistasis_df = pd.DataFrame(epistasis_results)
    return epistasis_df

def analyze_context_dependent_combinations(df):
    """Identify combinations observed in multiple contexts"""

    # Find combinations with multiple observations
    combo_counts = df[df['Mutation'].str.contains(r'\+', regex=True)].groupby('Mutation').size()
    recurrent_combos = combo_counts[combo_counts > 1].index.tolist()

    print(f"\nRecurrent combinations (observed >1 time): {len(recurrent_combos)}")

    context_results = []

    for combo in recurrent_combos:
        combo_data = df[df['Mutation'] == combo]

        print(f"\n{combo}: {len(combo_data)} observations")

        for _, row in combo_data.iterrows():
            context_results.append({
                'combination': combo,
                'observation_id': row['observation_id'],
                'log10_FC': row['log10_FC'],
                'FC': row['FC_numeric'],
                'context': row['context_tier'],
                'study': row['study_source'],
                'subtype': row['Subtype']
            })
            print(f"  {row['study_source']} ({row['context_tier']}): {row['FC_numeric']:.1f}-fold")

    context_df = pd.DataFrame(context_results)
    return context_df

def analyze_compensatory_patterns(df):
    """Analyze M66I-centered combinations for compensatory patterns"""

    # Find M66I alone
    m66i_alone = df[df['Mutation'] == 'M66I']

    if len(m66i_alone) == 0:
        print("\nWarning: M66I single mutant not found")
        m66i_log10fc = 3.505  # From literature
    else:
        m66i_log10fc = m66i_alone['log10_FC'].mean()

    print(f"\nM66I alone: log10_FC = {m66i_log10fc:.3f}")

    # Find M66I combinations
    m66i_combos = df[df['Mutation'].str.contains('M66I', regex=False) &
                     df['Mutation'].str.contains(r'\+', regex=True)]

    print(f"M66I combinations found: {len(m66i_combos)}")

    compensatory_results = []

    for _, row in m66i_combos.iterrows():
        combo = row['Mutation']
        observed_log10fc = row['log10_FC']
        observed_fc = row['FC_numeric']

        # Check if lower than M66I alone
        if observed_log10fc < m66i_log10fc:
            pattern = 'putative_compensatory'
            note = f"Lower than M66I alone ({10**m66i_log10fc:.0f}-fold)"
        else:
            pattern = 'additive_or_synergistic'
            note = f"Higher than or equal to M66I alone"

        compensatory_results.append({
            'combination': combo,
            'observed_log10fc': observed_log10fc,
            'observed_fc': observed_fc,
            'm66i_alone_log10fc': m66i_log10fc,
            'm66i_alone_fc': 10 ** m66i_log10fc,
            'difference_from_m66i': observed_log10fc - m66i_log10fc,
            'pattern': pattern,
            'note': note
        })

        print(f"  {combo}: {observed_fc:.0f}-fold ({pattern})")

    compensatory_df = pd.DataFrame(compensatory_results)
    return compensatory_df

def classify_mutations(df, epistasis_df, compensatory_df):
    """Classify all mutations by evidence strength"""

    classifications = []

    # Single mutants
    single_mutants = df[~df['Mutation'].str.contains(r'\+', regex=True)]

    for mutation in single_mutants['Mutation'].unique():
        mut_data = single_mutants[single_mutants['Mutation'] == mutation]
        mean_log10fc = mut_data['log10_FC'].mean()
        n_obs = len(mut_data)

        # Classify by resistance level
        if mean_log10fc > 3.0:  # >1000-fold
            category = 'primary_resistant'
            evidence = 'strong'
        elif mean_log10fc > 1.0:  # >10-fold
            category = 'moderate_resistant'
            evidence = 'moderate'
        else:
            category = 'low_resistant'
            evidence = 'moderate'

        classifications.append({
            'mutation': mutation,
            'type': 'single',
            'mean_log10fc': mean_log10fc,
            'n_observations': n_obs,
            'category': category,
            'evidence_strength': evidence
        })

    # Double mutants
    for _, row in epistasis_df.iterrows():
        combo = row['combination']

        if row['interaction_type'] == 'positive_synergy':
            category = 'synergistic'
            evidence = 'moderate'
        elif row['interaction_type'] == 'additive':
            category = 'additive'
            evidence = 'moderate'
        else:
            category = 'sub_additive'
            evidence = 'moderate'

        # Check if compensatory
        if combo in compensatory_df['combination'].values:
            comp_row = compensatory_df[compensatory_df['combination'] == combo].iloc[0]
            if comp_row['pattern'] == 'putative_compensatory':
                category = 'putative_compensatory'
                evidence = 'hypothesis_generating'

        classifications.append({
            'mutation': combo,
            'type': 'double',
            'mean_log10fc': row['observed_log10fc'],
            'n_observations': 1,
            'category': category,
            'evidence_strength': evidence
        })

    classification_df = pd.DataFrame(classifications)
    return classification_df

def main():
    """Main execution"""
    print("="*60)
    print("Epistasis and Compensatory Mutation Analysis")
    print("="*60)

    # Load data
    df = load_data()

    # Identify double mutants
    double_mutants = identify_double_mutants(df)

    # Calculate epistasis
    print("\nCalculating epistasis...")
    epistasis = calculate_epistasis(df, double_mutants)
    epistasis.to_csv(OUTPUT_DIR / "epistasis_matrix.csv", index=False)
    print(f"[OK] Epistasis analysis saved: {len(epistasis)} combinations")

    # Context-dependent combinations
    print("\nAnalyzing context-dependent combinations...")
    context_specific = analyze_context_dependent_combinations(df)
    context_specific.to_csv(OUTPUT_DIR / "context_specific_combinations.csv", index=False)
    print(f"[OK] Context-specific analysis saved")

    # Compensatory patterns
    print("\nAnalyzing compensatory patterns...")
    compensatory = analyze_compensatory_patterns(df)
    compensatory.to_csv(OUTPUT_DIR / "compensatory_patterns.csv", index=False)
    print(f"[OK] Compensatory analysis saved")

    # Classify all mutations
    print("\nClassifying mutations...")
    classifications = classify_mutations(df, epistasis, compensatory)
    classifications.to_csv(OUTPUT_DIR / "mutation_classifications.csv", index=False)
    print(f"[OK] Classifications saved: {len(classifications)} mutations")

    # Summary
    print("\n" + "="*60)
    print("Summary")
    print("="*60)

    print("\nEpistasis patterns:")
    for itype in epistasis['interaction_type'].unique():
        count = len(epistasis[epistasis['interaction_type'] == itype])
        print(f"  {itype}: {count}")

    print("\nCompensatory patterns:")
    for pattern in compensatory['pattern'].unique():
        count = len(compensatory[compensatory['pattern'] == pattern])
        print(f"  {pattern}: {count}")

    print("\nMutation categories:")
    for cat in classifications['category'].unique():
        count = len(classifications[classifications['category'] == cat])
        print(f"  {cat}: {count}")

    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()
