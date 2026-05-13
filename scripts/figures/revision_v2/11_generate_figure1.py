#!/usr/bin/env python3
"""
Figure 1: Evidence Landscape - PRISMA + Availability Matrix (revision_v2)
2 panels: PRISMA flow diagram, data availability matrix
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
from pathlib import Path

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
RESULTS_DIR = BASE_DIR / "results" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "manuscript" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Visual standardization for JMV submission
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['lines.linewidth'] = 1.2
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

# Muted color palette
SINGLE_COLOR = '#6B8BA4'   # blue-gray for single mutations
COMBO_COLOR = '#B85450'    # muted red for combination mutations
MATRIX_CMAP = 'Blues'      # sequential palette for availability matrix
PRISMA_COLOR = '#6B8BA4'   # blue-gray for PRISMA boxes

plt.style.use('seaborn-v0_8-paper')


def create_prisma_flow(ax):
    """Panel A: PRISMA 2020 flow diagram"""

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    ax.text(5, 9.6, 'A. PRISMA 2020 Flow Diagram',
            ha='center', fontsize=14, weight='bold')

    # PRISMA boxes with counts (larger boxes, better spacing)
    boxes = [
        (5, 8.0, 'Records identified\n(n = 87)', '82 PubMed + 5 additional sources'),
        (5, 5.8, 'Records after deduplication\n(n = 78)', ''),
        (5, 3.6, 'Full-text articles assessed\n(n = 22)', ''),
        (5, 1.4, 'Sources included in synthesis\n(n = 11)', '24 quantitative observations'),
    ]

    # Draw boxes
    for x, y, title, note in boxes:
        rect = mpatches.FancyBboxPatch((x - 2.0, y - 0.75), 4.0, 1.5,
                                        boxstyle="round,pad=0.08",
                                        facecolor='white',
                                        edgecolor=PRISMA_COLOR, linewidth=1.8)
        ax.add_patch(rect)
        ax.text(x, y, title, ha='center', va='center', fontsize=11, weight='bold')
        if note:
            ax.text(x, y - 0.85, note, ha='center', fontsize=9, style='italic', color='#555555')

    # Exclusion labels
    exclusion_labels = [
        (7.5, 7.2, 'Title/abstract screening'),
        (7.5, 6.6, 'Excluded (n = 56)'),
        (7.5, 4.6, 'Full-text eligibility'),
        (7.5, 4.0, 'Excluded (n = 11):'),
        (7.5, 3.4, '8 no quantitative data'),
        (7.5, 2.8, '3 HIV-2 only'),
    ]

    for x, y, label in exclusion_labels:
        ax.text(x, y, label, ha='left', fontsize=9, color='#555555')

    # Exclusion box (dashed)
    excl_rect = mpatches.FancyBboxPatch((7.0, 2.3), 3.0, 2.6,
                                         boxstyle="round,pad=0.05",
                                         facecolor='#F8F8F8',
                                         edgecolor='#B85450', linewidth=1.3,
                                         linestyle='--')
    ax.add_patch(excl_rect)

    # Arrows between main boxes
    arrow_props = dict(arrowstyle='->', color=PRISMA_COLOR, lw=2.0)
    ax.annotate('', xy=(5, 6.65), xytext=(5, 7.15), arrowprops=arrow_props)
    ax.annotate('', xy=(5, 4.45), xytext=(5, 4.95), arrowprops=arrow_props)
    ax.annotate('', xy=(5, 2.25), xytext=(5, 2.75), arrowprops=arrow_props)

    # Arrow from screening to exclusion
    ax.annotate('', xy=(7.2, 6.7), xytext=(6.85, 6.7),
                arrowprops=dict(arrowstyle='->', color='#B85450', lw=1.2, linestyle='dashed'))

    # Search date annotation
    ax.text(5, 0.15, 'Search: PubMed (2025-12) | 82 PubMed + 5 additional sources',
            ha='center', fontsize=8, color='#888888', style='italic')


def create_availability_matrix(ax):
    """Panel B: Mutation × Subtype × Context availability heatmap"""

    avail_df = pd.read_csv(DATA_DIR / "availability_matrix.csv")

    # Remove aggregate/control rows
    aggregate_mutations = ['GCSMs_median', 'SDM_GCSMs', 'Multiclass_resistant']
    avail_df = avail_df[~avail_df['mutation'].isin(aggregate_mutations)].copy()

    # Create pivot table
    avail_df['subtype_context'] = avail_df['subtype'] + '_' + avail_df['context']

    pivot = avail_df.pivot_table(
        index='mutation',
        columns='subtype_context',
        values='n_observations',
        fill_value=0
    )

    # Sort by single vs combination (add horizontal separator)
    mutation_order = pivot.index.tolist()
    single_mutations = [m for m in mutation_order if '+' not in m]
    combo_mutations = [m for m in mutation_order if '+' in m]
    sorted_mutation_order = single_mutations + combo_mutations
    pivot = pivot.reindex(sorted_mutation_order)

    # Plot heatmap with sequential palette (Blues)
    sns.heatmap(pivot, annot=True, fmt='g', cmap=MATRIX_CMAP,
                cbar_kws={'label': 'N observations'},
                linewidths=0.5, linecolor='gray',
                ax=ax, vmin=0, vmax=5)

    # Add horizontal separator between single and combination mutations
    if single_mutations and combo_mutations:
        n_single = len(single_mutations)
        ax.axhline(y=n_single, color='black', linewidth=1.5)

    ax.set_xlabel('Subtype × Context', fontsize=9)
    ax.set_ylabel('Mutation', fontsize=9)
    ax.set_title('B. Data Availability Matrix', fontsize=11, weight='bold', loc='left')
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)

    # Add mutation grouping labels
    if single_mutations and combo_mutations:
        n_single = len(single_mutations)
        ax.text(-0.5, n_single / 2 - 0.5, 'Single\nmutants',
                ha='center', va='center', fontsize=7, color=SINGLE_COLOR, weight='bold',
                transform=ax.get_yaxis_transform(), rotation=0)
        ax.text(-0.5, n_single + len(combo_mutations) / 2 - 0.5, 'Combination\nmutants',
                ha='center', va='center', fontsize=7, color=COMBO_COLOR, weight='bold',
                transform=ax.get_yaxis_transform(), rotation=0)


def main():
    """Generate Figure 1 - Evidence Landscape (2 panels: PRISMA + matrix)"""
    print("="*60)
    print("Generating Figure 1: Evidence Landscape (PRISMA + Availability Matrix)")
    print("="*60)

    # Create figure with 2 panels (1x2 grid)
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(1, 2, figure=fig, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Generate panels
    print("\nPanel A: PRISMA flow diagram...")
    create_prisma_flow(ax1)

    print("Panel B: Availability matrix...")
    create_availability_matrix(ax2)

    # Save
    fig.subplots_adjust(wspace=0.35)

    output_pdf = OUTPUT_DIR / "figure1_evidence_landscape.pdf"
    output_png = OUTPUT_DIR / "figure1_evidence_landscape.png"

    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)

    print(f"\n[OK] Saved: {output_pdf}")
    print(f"[OK] Saved: {output_png}")
    print("\n" + "="*60)
    print("Figure 1 complete!")
    print("="*60)

    plt.close()


if __name__ == "__main__":
    main()
