#!/usr/bin/env python3
"""
Genetic Barrier Analysis for HIV-1 Capsid LEN Resistance Mutations
Calculates nucleotide substitution requirements for each RAM

Usage: python 18_genetic_barrier.py
Output: results/revision_v2/genetic_barrier_scores.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
RESULTS_DIR = BASE_DIR / "results" / "revision_v2"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# HIV-1 codon tables (standard genetic code)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Amino acid properties for classification
AA_PROPERTY = {
    'A': 'hydrophobic', 'C': 'special', 'D': 'negative', 'E': 'negative',
    'F': 'hydrophobic', 'G': 'special', 'H': 'positive', 'I': 'hydrophobic',
    'K': 'positive', 'L': 'hydrophobic', 'M': 'hydrophobic', 'N': 'polar',
    'P': 'special', 'Q': 'polar', 'R': 'positive', 'S': 'polar',
    'T': 'polar', 'V': 'hydrophobic', 'W': 'hydrophobic', 'Y': 'hydrophobic',
    '*': 'stop'
}

# HIV-1 consensus B subtype codons for key CA residues
HIV_CONSENSUS_CODONS = {
    'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAA', 'K70': 'AAA',
    'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC',
    'L56': 'CTG', 'P1': 'CCA',
}

# Resistance mutations and their WT/target codons
RESISTANCE_MUTATIONS = {
    'N57H': {'wt_codon': 'AAT', 'mut_codon': 'CAT', 'position': 'N57', 'fitness_cost': 'low'},
    'M66I': {'wt_codon': 'ATG', 'mut_codon': 'ATA', 'position': 'M66', 'fitness_cost': 'high'},
    'Q67H': {'wt_codon': 'CAA', 'mut_codon': 'CAC', 'position': 'Q67', 'fitness_cost': 'low'},
    'K70R': {'wt_codon': 'AAA', 'mut_codon': 'AGA', 'position': 'K70', 'fitness_cost': 'low'},
    'K70N': {'wt_codon': 'AAA', 'mut_codon': 'AAT', 'position': 'K70', 'fitness_cost': 'moderate'},
    'N74D': {'wt_codon': 'AAT', 'mut_codon': 'GAT', 'position': 'N74', 'fitness_cost': 'moderate'},
    'A105T': {'wt_codon': 'GCT', 'mut_codon': 'ACT', 'position': 'A105', 'fitness_cost': 'low'},
    'T107A': {'wt_codon': 'ACC', 'mut_codon': 'GCC', 'position': 'T107', 'fitness_cost': 'low'},
    'A105V': {'wt_codon': 'GCT', 'mut_codon': 'GTT', 'position': 'A105', 'fitness_cost': 'moderate'},
}

# Subtype-specific codon usage (representative)
SUBTYPE_CODONS = {
    'B': {
        'N57': 'AAT', 'M66': 'ATG', 'Q67': ['CAA', 'CAG'], 'K70': 'AAA',
        'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC',
    },
    'C': {
        'N57': 'AAT', 'M66': 'ATG', 'Q67': ['CAA'], 'K70': ['AAA', 'AAG'],
        'N74': 'AAT', 'A105': ['GCT', 'GCC'], 'T107': 'ACC',
    },
    'CRF01_AE': {
        'N57': 'AAT', 'M66': 'ATG', 'Q67': 'CAA', 'K70': 'AAA',
        'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC',
    },
    'CRF02_AG': {
        'N57': 'AAT', 'M66': 'ATG', 'Q67': ['CAA', 'CAG'], 'K70': 'AAA',
        'N74': 'AAT', 'A105': 'GCT', 'T107': 'ACC',
    },
}


def count_nucleotide_changes(codon1, codon2):
    """Count nucleotide differences between two codons"""
    return sum(c1 != c2 for c1, c2 in zip(codon1, codon2))


def get_synonymous_codons(aa, codon_table=CODON_TABLE):
    """Get all codons that code for the same amino acid"""
    return [codon for codon, aa2 in codon_table.items() if aa2 == aa and aa != '*']


def calculate_barrier_score(wt_codon, mut_codon, aa_property_change=True):
    """
    Calculate genetic barrier score (1-10 scale)
    Higher = harder to acquire (more nucleotide changes, property change)
    """
    n_changes = count_nucleotide_changes(wt_codon, mut_codon)
    wt_aa = CODON_TABLE.get(wt_codon, 'X')
    mut_aa = CODON_TABLE.get(mut_codon, 'X')

    # Base score from nucleotide distance
    score = n_changes * 2

    # Property change penalty
    if aa_property_change and wt_aa in AA_PROPERTY and mut_aa in AA_PROPERTY:
        if AA_PROPERTY[wt_aa] != AA_PROPERTY[mut_aa]:
            score += 1.5

    # Synonymous option penalty (if synonymous paths exist, barrier is lower)
    syn_wt = get_synonymous_codons(wt_aa)
    syn_mut = get_synonymous_codons(mut_aa)
    if len(syn_wt) > 1:
        score -= 0.5  # More WT options means easier transition
    if len(syn_mut) > 1:
        score -= 0.5  # More target options means easier to hit

    return min(10, max(1, score))


def classify_barrier(score):
    """Classify barrier level from score"""
    if score >= 7:
        return 'High'
    elif score >= 4:
        return 'Moderate'
    else:
        return 'Low'


def calculate_escape_probability(n_changes, fitness_cost=0.0):
    """
    Estimate relative escape probability based on barrier and fitness cost
    Uses empirical approximation based on HIV-1 within-patient evolution rates
    """
    # 1 substitution ≈ 10^-4 probability per replication cycle
    # 2 substitutions ≈ 10^-8
    # 3 substitutions ≈ 10^-12

    prob_per_cycle = 10 ** (-4 * n_changes)

    # Fitness cost reduces effective probability
    if fitness_cost == 'high':
        prob_per_cycle *= 0.01
    elif fitness_cost == 'moderate':
        prob_per_cycle *= 0.1
    else:
        prob_per_cycle *= 0.5

    return prob_per_cycle


def analyze_mutation(mut_name, mut_info, subtypes=None):
    """Analyze genetic barrier for a single mutation across subtypes"""
    if subtypes is None:
        subtypes = list(SUBTYPE_CODONS.keys())

    wt_codon_raw = mut_info['wt_codon']
    mut_codon_raw = mut_info['mut_codon']
    position = mut_info['position']
    fitness = mut_info.get('fitness_cost', 'moderate')

    results = []

    for subtype in subtypes:
        if subtype not in SUBTYPE_CODONS:
            continue

        sub_codons = SUBTYPE_CODONS[subtype]
        if position not in sub_codons:
            continue

        wt_codon_list = sub_codons[position]
        if isinstance(wt_codon_list, str):
            wt_codon_list = [wt_codon_list]

        for wt_codon in wt_codon_list:
            n_changes = count_nucleotide_changes(wt_codon, mut_codon_raw)
            score = calculate_barrier_score(wt_codon, mut_codon_raw, aa_property_change=True)
            barrier_class = classify_barrier(score)
            escape_prob = calculate_escape_probability(n_changes, fitness_cost=fitness)

            results.append({
                'mutation': mut_name,
                'subtype': subtype,
                'wt_codon': wt_codon,
                'mut_codon': mut_codon_raw,
                'n_nucleotide_changes': n_changes,
                'barrier_score': round(score, 2),
                'barrier_class': barrier_class,
                'fitness_cost': fitness,
                'relative_escape_prob': escape_prob,
            })

    return results


def analyze_combination_mutation(mut_name, mutations, subtypes=None):
    """Analyze genetic barrier for combination mutations (double mutants)"""
    if subtypes is None:
        subtypes = list(SUBTYPE_CODONS.keys())

    results = []

    for subtype in subtypes:
        total_changes = 0
        for mut in mutations:
            if mut in RESISTANCE_MUTATIONS:
                mut_info = RESISTANCE_MUTATIONS[mut]
                position = mut_info['position']
                if subtype in SUBTYPE_CODONS and position in SUBTYPE_CODONS[subtype]:
                    wt_codon_list = SUBTYPE_CODONS[subtype][position]
                    if isinstance(wt_codon_list, str):
                        wt_codon_list = [wt_codon_list]
                    # Take first codon as representative
                    wt_codon = wt_codon_list[0]
                    n_changes = count_nucleotide_changes(wt_codon, mut_info['mut_codon'])
                    total_changes += n_changes

        # Combined barrier
        combined_score = min(10, total_changes * 1.5)
        barrier_class = classify_barrier(combined_score)

        results.append({
            'mutation': mut_name,
            'subtype': subtype,
            'wt_codon': 'multi',
            'mut_codon': 'multi',
            'n_nucleotide_changes': total_changes,
            'barrier_score': round(combined_score, 2),
            'barrier_class': barrier_class,
            'fitness_cost': 'high',
            'relative_escape_prob': 10 ** (-4 * total_changes),
        })

    return results


def main():
    """Main analysis: calculate genetic barriers for all LEN resistance mutations"""
    print("="*60)
    print("Genetic Barrier Analysis for HIV-1 CA LEN Resistance")
    print("="*60)

    all_results = []

    # Single mutations
    print("\n[1] Analyzing single resistance mutations...")
    for mut_name, mut_info in RESISTANCE_MUTATIONS.items():
        results = analyze_mutation(mut_name, mut_info)
        all_results.extend(results)
        print(f"  {mut_name}: {len(results)} subtype-codon combinations analyzed")

    # Double mutations (combinations)
    print("\n[2] Analyzing combination mutations...")
    combinations = [
        ('Q67H+K70R', ['Q67H', 'K70R']),
        ('Q67H+N74D', ['Q67H', 'N74D']),
        ('M66I+A105T', ['M66I', 'A105T']),
        ('M66I+T107A', ['M66I', 'T107A']),
        ('A105T+T107A', ['A105T', 'T107A']),
    ]

    for combo_name, mutations in combinations:
        results = analyze_combination_mutation(combo_name, mutations)
        all_results.extend(results)
        print(f"  {combo_name}: {len(results)} subtype combinations analyzed")

    # Create DataFrame
    df = pd.DataFrame(all_results)

    # Summary by mutation (using B subtype as representative)
    print("\n[3] Generating summary...")
    summary = df[df['subtype'] == 'B'].groupby('mutation').agg({
        'n_nucleotide_changes': 'first',
        'barrier_score': 'first',
        'barrier_class': 'first',
        'fitness_cost': 'first',
        'relative_escape_prob': 'first',
    }).reset_index()

    summary = summary.sort_values('barrier_score', ascending=False)

    # Save outputs
    output_csv = RESULTS_DIR / "genetic_barrier_scores.csv"
    df.to_csv(output_csv, index=False)
    print(f"\n[OK] Saved: {output_csv}")

    summary_csv = RESULTS_DIR / "genetic_barrier_summary.csv"
    summary.to_csv(summary_csv, index=False)
    print(f"[OK] Saved: {summary_csv}")

    # Print summary table
    print("\n" + "="*60)
    print("Genetic Barrier Summary (B subtype)")
    print("="*60)
    print(f"{'Mutation':<20} {'Changes':<10} {'Score':<8} {'Class':<12} {'Fitness':<10}")
    print("-"*60)
    for _, row in summary.iterrows():
        print(f"{row['mutation']:<20} {row['n_nucleotide_changes']:<10} "
              f"{row['barrier_score']:<8.2f} {row['barrier_class']:<12} {row['fitness_cost']:<10}")

    print("\n" + "="*60)
    print("Barriers by Class:")
    print("-"*60)
    for barrier_class in ['High', 'Moderate', 'Low']:
        muts = summary[summary['barrier_class'] == barrier_class]['mutation'].tolist()
        print(f"\n{barrier_class} barrier ({len(muts)} mutations):")
        for m in muts:
            print(f"  - {m}")

    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

    return df, summary


if __name__ == "__main__":
    df, summary = main()