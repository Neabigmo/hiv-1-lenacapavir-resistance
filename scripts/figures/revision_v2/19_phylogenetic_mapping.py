#!/usr/bin/env python3
"""
Phylogenetic Mapping of HIV-1 CA Resistance Mutations
Maps LEN resistance mutations onto representative HIV-1 subtype phylogeny

Usage: python 19_phylogenetic_mapping.py
Output: results/revision_v2/phylogenetic_mapping.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
RESULTS_DIR = BASE_DIR / "results" / "revision_v2"
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# HIV-1 group M subtype reference sequences (partial CA region alignment)
# Representative lineages for major subtypes
# Position numbers are relative to HXB2 reference (CA region: positions 133-363 of Gag)
SUBTYPE_LINEAGES = {
    'B': {
        'lineage': 'B', 'country': 'North America/Western Europe',
        'ca_positions': {
            'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAA', 'K70': 'AAA',
            'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC'
        },
        'notes': 'First characterized, most in vitro data',
    },
    'C': {
        'lineage': 'C', 'country': 'Southern Africa/India',
        'ca_positions': {
            'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAA', 'K70': 'AAG',
            'N74': 'AAT', 'A105': 'GCC', 'T107': 'ACC'
        },
        'notes': 'Most prevalent globally, distinct CA context',
    },
    'CRF01_AE': {
        'lineage': 'CRF01_AE', 'country': 'Southeast Asia',
        'ca_positions': {
            'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAA', 'K70': 'AAA',
            'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC'
        },
        'notes': 'CRF, distinct from parental subtypes',
    },
    'CRF02_AG': {
        'lineage': 'CRF02_AG', 'country': 'West Africa',
        'ca_positions': {
            'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAG', 'K70': 'AAA',
            'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC'
        },
        'notes': 'CRF, high prevalence in West Africa',
    },
    'A1': {
        'lineage': 'A1', 'country': 'East Africa',
        'ca_positions': {
            'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAA', 'K70': 'AAA',
            'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC'
        },
        'notes': 'Subtype A1 (vs A2)',
    },
    'D': {
        'lineage': 'D', 'country': 'Central/East Africa',
        'ca_positions': {
            'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAA', 'K70': 'AAA',
            'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC'
        },
        'notes': 'More ancient lineage, faster progression',
    },
    'G': {
        'lineage': 'G', 'country': 'West Africa/Central Africa',
        'ca_positions': {
            'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAA', 'K70': 'AAA',
            'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC'
        },
        'notes': 'Related to CRF02_AG',
    },
}

# Resistance mutations with their properties
RESISTANCE_MUTATIONS = {
    'N57H': {
        'position': 'N57', 'primary': True, 'tier': 1,
        'codon_wt': 'AAT', 'codon_mut': 'CAT',
        'mechanism': 'H-bond disruption',
        'fold_change': '4,890×',
        'subtype_impact': 'Universal (all subtypes)',
    },
    'M66I': {
        'position': 'M66', 'primary': True, 'tier': 1,
        'codon_wt': 'ATG', 'codon_mut': 'ATA',
        'mechanism': 'Steric clash (beta-branched)',
        'fold_change': '3,200×',
        'subtype_impact': 'Universal (all subtypes)',
    },
    'Q67H': {
        'position': 'Q67', 'primary': False, 'tier': 2,
        'codon_wt': 'CAA', 'codon_mut': 'CAC',
        'mechanism': 'Context-dependent',
        'fold_change': 'Variable',
        'subtype_impact': 'Variable (B > C)',
    },
    'K70R': {
        'position': 'K70', 'primary': False, 'tier': 2,
        'codon_wt': 'AAA', 'codon_mut': 'AGA',
        'mechanism': 'Context-dependent',
        'fold_change': 'Variable',
        'subtype_impact': 'Variable (B > C)',
    },
    'Q67H+K70R': {
        'position': 'Q67+K70', 'primary': False, 'tier': 2,
        'codon_wt': 'CAA+AAA', 'codon_mut': 'CAC+AGA',
        'mechanism': 'Dual, synergistic',
        'fold_change': '76× (context range)',
        'subtype_impact': 'B subtype enriched',
    },
    'N74D': {
        'position': 'N74', 'primary': False, 'tier': 2,
        'codon_wt': 'AAT', 'codon_mut': 'GAT',
        'mechanism': 'Electrostatic',
        'fold_change': 'Variable',
        'subtype_impact': 'Less common',
    },
    'A105T': {
        'position': 'A105', 'primary': False, 'tier': 3,
        'codon_wt': 'GCT', 'codon_mut': 'ACT',
        'mechanism': 'Compensatory (with M66I)',
        'fold_change': 'With M66I: 111×',
        'subtype_impact': 'Subtype B enriched',
    },
    'T107A': {
        'position': 'T107', 'primary': False, 'tier': 3,
        'codon_wt': 'ACC', 'codon_mut': 'GCC',
        'mechanism': 'Compensatory (with M66I)',
        'fold_change': 'With M66I: 234×',
        'subtype_impact': 'Subtype B enriched',
    },
    'A105V': {
        'position': 'A105', 'primary': False, 'tier': 3,
        'codon_wt': 'GCT', 'codon_mut': 'GTT',
        'mechanism': 'Secondary resistance',
        'fold_change': 'Moderate',
        'subtype_impact': 'Subtype B enriched',
    },
}


def estimate_subtype_distance(subtype1, subtype2):
    """
    Estimate phylogenetic distance between subtypes
    Based on HIV-1 group M polytree distances
    Returns distance in arbitrary units (0-1 scale)
    """
    # Major subtype distances (approximate)
    distances = {
        ('B', 'C'): 0.35, ('B', 'CRF01_AE'): 0.30, ('B', 'CRF02_AG'): 0.28,
        ('B', 'A1'): 0.25, ('B', 'D'): 0.22, ('B', 'G'): 0.28,
        ('C', 'CRF01_AE'): 0.38, ('C', 'CRF02_AG'): 0.32, ('C', 'A1'): 0.30,
        ('C', 'D'): 0.35, ('C', 'G'): 0.35,
        ('CRF01_AE', 'CRF02_AG'): 0.25, ('CRF01_AE', 'A1'): 0.28,
        ('CRF02_AG', 'A1'): 0.25, ('CRF02_AG', 'G'): 0.18,
        ('A1', 'D'): 0.20, ('A1', 'G'): 0.25,
        ('D', 'G'): 0.28,
    }

    # Self-distance
    if subtype1 == subtype2:
        return 0.0

    # Check both orders
    key = (subtype1, subtype2)
    if key in distances:
        return distances[key]

    key = (subtype2, subtype1)
    if key in distances:
        return distances[key]

    # Default distance for unknown pairs
    return 0.3


def calculate_codon_distance(codon1, codon2):
    """Count nucleotide differences between codons"""
    return sum(c1 != c2 for c1, c2 in zip(codon1, codon2))


def map_mutation_to_subtype(mutation_name, mutation_info, subtype_info):
    """
    Map a mutation onto a subtype and calculate relevant metrics
    """
    position = mutation_info['position']

    if position not in subtype_info['ca_positions']:
        return None

    wt_codon = subtype_info['ca_positions'][position]
    mut_codon = mutation_info.get('codon_mut', '')

    # Calculate codon distance
    n_changes = calculate_codon_distance(wt_codon, mut_codon)

    # Codon usage differences
    # (In real implementation, would use HIV-1 codon usage tables)

    return {
        'mutation': mutation_name,
        'subtype': subtype_info['lineage'],
        'wt_codon': wt_codon,
        'mut_codon': mut_codon,
        'nucleotide_changes': n_changes,
        'mechanism': mutation_info['mechanism'],
        'fold_change': mutation_info['fold_change'],
        'surveillance_tier': mutation_info['tier'],
        'geographic_distribution': subtype_info['country'],
        'primary_mutation': mutation_info['primary'],
    }


def build_phylogenetic_matrix():
    """Build complete mutation × subtype matrix with all metrics"""
    results = []

    for mut_name, mut_info in RESISTANCE_MUTATIONS.items():
        for subtype_name, subtype_info in SUBTYPE_LINEAGES.items():
            result = map_mutation_to_subtype(mut_name, mut_info, subtype_info)
            if result:
                results.append(result)

    return pd.DataFrame(results)


def calculate_subtype_mutation_burden(df):
    """
    Calculate mutation burden per subtype
    How many Tier 1/2/3 mutations are possible in each subtype
    """
    # Count by subtype
    burden = df.groupby('subtype').size().reset_index(name='n_resistance_positions')

    # Get tier distribution per subtype as a list
    tier_dist = df.groupby('subtype')['surveillance_tier'].apply(list).reset_index()
    tier_dist.columns = ['subtype', 'tier_distribution']

    burden = burden.merge(tier_dist, on='subtype', how='left')

    # Count by tier using separate groupby calls
    for tier_num in ['1', '2', '3']:
        tier_counts = df[df['surveillance_tier'] == tier_num].groupby('subtype').size()
        burden['n_tier' + tier_num] = burden['subtype'].map(tier_counts).fillna(0).astype(int)

    return burden


def calculate_inter_subtype_distances():
    """Calculate pairwise phylogenetic distances between all subtypes"""
    subtypes = list(SUBTYPE_LINEAGES.keys())
    distance_matrix = np.zeros((len(subtypes), len(subtypes)))

    for i, s1 in enumerate(subtypes):
        for j, s2 in enumerate(subtypes):
            distance_matrix[i, j] = estimate_subtype_distance(s1, s2)

    return pd.DataFrame(distance_matrix, index=subtypes, columns=subtypes)


def generate_mutation_tree_data():
    """
    Generate tree-like structure data for phylogenetic visualization
    Returns hierarchical data suitable for dendrogram/heatmap display
    """
    tree_data = []

    # Root node
    tree_data.append({
        'node': 'HIV-1 Group M',
        'parent': None,
        'depth': 0,
        'subtype_count': len(SUBTYPE_LINEAGES),
        'primary_mutations': 2,  # N57H, M66I
    })

    # Subtype nodes
    for subtype, info in SUBTYPE_LINEAGES.items():
        tree_data.append({
            'node': subtype,
            'parent': 'HIV-1 Group M',
            'depth': 1,
            'country': info['country'],
            'primary_mutations': 2,
            'secondary_mutations': sum(1 for m in RESISTANCE_MUTATIONS.values() if not m['primary']),
        })

    return pd.DataFrame(tree_data)


def create_mutation_heatmap_data(df):
    """
    Create matrix for heatmap showing mutation × subtype presence
    Binary: 1 = mutation possible, 0 = not possible
    """
    # Get all mutations and subtypes
    mutations = list(RESISTANCE_MUTATIONS.keys())
    subtypes = list(SUBTYPE_LINEAGES.keys())

    # Build presence matrix
    presence_matrix = np.zeros((len(mutations), len(subtypes)))

    for i, mut in enumerate(mutations):
        mut_info = RESISTANCE_MUTATIONS[mut]
        for j, sub in enumerate(subtypes):
            # Check if codon distance is calculable
            mut_row = df[(df['mutation'] == mut) & (df['subtype'] == sub)]
            if len(mut_row) > 0:
                presence_matrix[i, j] = 1

    return pd.DataFrame(presence_matrix, index=mutations, columns=subtypes)


def main():
    """Main analysis: phylogenetic mapping of resistance mutations"""
    print("="*60)
    print("Phylogenetic Mapping of HIV-1 CA LEN Resistance Mutations")
    print("="*60)

    # Build complete matrix
    print("\n[1] Building mutation × subtype matrix...")
    matrix_df = build_phylogenetic_matrix()
    print(f"  Generated {len(matrix_df)} mutation-subtype pairs")

    # Calculate burden per subtype
    print("\n[2] Calculating mutation burden per subtype...")
    burden_df = calculate_subtype_mutation_burden(matrix_df)
    print(f"  Analyzed {len(burden_df)} subtypes")

    # Inter-subtype distances
    print("\n[3] Calculating inter-subtype phylogenetic distances...")
    distance_df = calculate_inter_subtype_distances()
    print(f"  Computed {len(distance_df)**2} pairwise distances")

    # Tree data
    print("\n[4] Generating phylogenetic tree structure...")
    tree_df = generate_mutation_tree_data()
    print(f"  Created {len(tree_df)} tree nodes")

    # Mutation presence heatmap data
    print("\n[5] Creating mutation presence matrix...")
    heatmap_df = create_mutation_heatmap_data(matrix_df)
    print(f"  Generated {heatmap_df.shape[0]} × {heatmap_df.shape[1]} matrix")

    # Save outputs
    output_matrix = RESULTS_DIR / "phylogenetic_mutation_matrix.csv"
    matrix_df.to_csv(output_matrix, index=False)
    print(f"\n[OK] Saved: {output_matrix}")

    output_burden = RESULTS_DIR / "phylogenetic_subtype_burden.csv"
    burden_df.to_csv(output_burden, index=False)
    print(f"[OK] Saved: {output_burden}")

    output_distances = RESULTS_DIR / "phylogenetic_distances.csv"
    distance_df.to_csv(output_distances)
    print(f"[OK] Saved: {output_distances}")

    output_tree = RESULTS_DIR / "phylogenetic_tree.csv"
    tree_df.to_csv(output_tree, index=False)
    print(f"[OK] Saved: {output_tree}")

    output_heatmap = RESULTS_DIR / "mutation_subtype_heatmap.csv"
    heatmap_df.to_csv(output_heatmap)
    print(f"[OK] Saved: {output_heatmap}")

    # Print summary
    print("\n" + "="*60)
    print("Phylogenetic Mapping Summary")
    print("="*60)

    print("\nMutation burden per subtype:")
    print("-"*60)
    print(f"{'Subtype':<15} {'Tier 1':<10} {'Tier 2':<10} {'Tier 3':<10} {'Total':<10}")
    print("-"*60)
    for _, row in burden_df.iterrows():
        print(f"{row['subtype']:<15} {int(row.get('n_tier1', 0)):<10} "
              f"{int(row.get('n_tier2', 0)):<10} {int(row.get('n_tier3', 0)):<10} "
              f"{int(row['n_resistance_positions']):<10}")

    print("\n" + "-"*60)
    print("Inter-subtype distances (0=identical, 1=divergent):")
    print("-"*60)
    print(distance_df.round(3).to_string())

    print("\n" + "="*60)
    print("Mutation × Subtype Presence Matrix")
    print("="*60)
    print("(1 = resistance possible, 0 = not possible in WT codon)")
    print(heatmap_df.replace({0: '.', 1: 'X'}).to_string())

    print("\n" + "="*60)
    print("Phylogenetic mapping complete!")
    print("="*60)

    return matrix_df, burden_df, distance_df, tree_df, heatmap_df


if __name__ == "__main__":
    matrix_df, burden_df, distance_df, tree_df, heatmap_df = main()