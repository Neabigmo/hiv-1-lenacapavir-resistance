#!/usr/bin/env python3
"""
Figure 1: Evidence Landscape - Redesigned UpSet Study-Centric View (revision_v2)
Shows 11 independent studies as columns with 3 data source types colored
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "manuscript" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Color palette for 3 data source types (light red, blue, green)
COLORS = {
    'Clinical Trial': '#FFB3B3',      # Light red
    'In Vitro': '#B3D9FF',            # Light blue
    'Population': '#B3FFB3'           # Light green
}

COLOR_HEX = {
    'Clinical Trial': '#FFB3B3',
    'In Vitro': '#B3D9FF',
    'Population': '#B3FFB3'
}


def categorize_source(study_source):
    """Categorize study into one of 3 data source types"""
    clinical_trials = ['CAPELLA', 'JID2025']
    population_studies = ['EMJ2025']
    # Everything else is In Vitro

    if study_source in clinical_trials:
        return 'Clinical Trial'
    elif study_source in population_studies:
        return 'Population'
    else:
        return 'In Vitro'


def load_and_prepare_data():
    """Load and prepare the data with source categorization"""
    pheno_df = pd.read_csv(DATA_DIR / "harmonized_phenotype_data.csv")
    pheno_df['source_category'] = pheno_df['study_source'].apply(categorize_source)

    # Get unique studies in consistent order
    studies_order = ['CAPELLA', 'JID2025', 'NATAP2022', 'PMC8092519',
                     'PMC12077089', 'PMC9600929', 'PMID40231850', 'Structural',
                     'AAC', 'JAC2025', 'EMJ2025']

    # Filter to only studies that exist in data
    studies_order = [s for s in studies_order if s in pheno_df['study_source'].unique()]

    # Get unique mutations (excluding non-resistance entries)
    resistance_mutations = [
        'Q67H', 'Q67K', 'Q67H+K70R', 'Q67H+N74D', 'Q67H+T107N', 'Q67H+N74S',
        'K70N', 'K70N+N74K',
        'N74D', 'L56V', 'N57H', 'L56V+N57H',
        'M66I', 'M66I+A105T', 'M66I+T107A', 'M66I+N74D', 'M66I+N74D+A105T',
        'Q67K+K70H', 'K436E+I437T',
        'GCSMs_median', 'SDM_GCSMs', 'Multiclass_resistant'
    ]

    mutations_in_data = [m for m in resistance_mutations if m in pheno_df['Mutation'].unique()]

    return pheno_df, studies_order, mutations_in_data


def create_study_centric_upset(ax_matrix, ax_bars, ax_legend, df, studies_order, mutations):
    """Create study-centric UpSet visualization

    Args:
        ax_matrix: Axes for the study × mutation matrix
        ax_bars: Axes for the left stacked bar chart
        ax_legend: Axes for legend
        df: harmonized phenotype dataframe
        studies_order: ordered list of 11 studies
        mutations: list of mutations to display
    """

    # Build presence matrix: which study reported which mutation
    presence_matrix = pd.DataFrame(0, index=studies_order, columns=mutations)

    # Also track the count per cell and the category
    for _, row in df.iterrows():
        study = row['study_source']
        mutation = row['Mutation']
        category = row['source_category']
        if study in presence_matrix.index and mutation in presence_matrix.columns:
            presence_matrix.loc[study, mutation] += 1

    # Calculate stacked bar values for each mutation (left side)
    # For each mutation, count how many observations from each category
    mutation_category_counts = {}
    for mut in mutations:
        mut_data = df[df['Mutation'] == mut]
        counts = {
            'Clinical Trial': len(mut_data[mut_data['source_category'] == 'Clinical Trial']),
            'In Vitro': len(mut_data[mut_data['source_category'] == 'In Vitro']),
            'Population': len(mut_data[mut_data['source_category'] == 'Population'])
        }
        mutation_category_counts[mut] = counts

    # Sort mutations by total count (descending)
    total_counts = {mut: sum(mutation_category_counts[mut].values()) for mut in mutations}
    sorted_mutations = sorted(mutations, key=lambda x: total_counts.get(x, 0), reverse=True)

    # Draw the left stacked bar chart
    bar_height = 0.8
    y_positions = np.arange(len(sorted_mutations))

    left = np.zeros(len(sorted_mutations))
    for category in ['Clinical Trial', 'In Vitro', 'Population']:
        heights = [mutation_category_counts[mut].get(category, 0) for mut in sorted_mutations]
        ax_bars.barh(y_positions, heights, left=left, height=bar_height,
                     color=COLORS[category], edgecolor='white', linewidth=0.5)
        left = left + np.array(heights)

    ax_bars.set_yticks(y_positions)
    ax_bars.set_yticklabels(sorted_mutations, fontsize=8)
    ax_bars.set_xlabel('N observations', fontsize=9)
    ax_bars.set_xlim(0, max(6, left.max() * 1.1))
    ax_bars.grid(axis='x', alpha=0.3)

    # Highlight combination mutations in yticklabels
    for i, mut in enumerate(sorted_mutations):
        if '+' in mut:
            ax_bars.get_yticklabels()[i].set_color('#A23B72')
            ax_bars.get_yticklabels()[i].set_fontweight('bold')

    # Draw the matrix (studies as columns, mutations as rows)
    for i, study in enumerate(studies_order):
        category = categorize_source(study)
        for j, mut in enumerate(sorted_mutations):
            count = presence_matrix.loc[study, mut] if study in presence_matrix.index else 0
            if count > 0:
                # Draw a dot/marker
                ax_matrix.scatter(i, j, s=count * 100, c=COLORS[category],
                                  edgecolors='black', linewidth=0.5, zorder=3)
            else:
                # Draw empty circle for presence indicator
                ax_matrix.scatter(i, j, s=30, c='lightgray',
                                  edgecolors='gray', linewidth=0.3, zorder=2, alpha=0.3)

    ax_matrix.set_xticks(range(len(studies_order)))
    ax_matrix.set_xticklabels(studies_order, rotation=45, ha='right', fontsize=7)
    ax_matrix.set_yticks(range(len(sorted_mutations)))
    ax_matrix.set_yticklabels([])  # No need for y-labels on matrix
    ax_matrix.set_xlim(-0.5, len(studies_order) - 0.5)
    ax_matrix.set_ylim(-0.5, len(sorted_mutations) - 0.5)
    ax_matrix.invert_yaxis()
    ax_matrix.tick_params(axis='x', labelsize=7)

    # Add grid
    ax_matrix.set_axisbelow(True)
    for i in range(len(studies_order)):
        ax_matrix.axvline(i - 0.5, color='lightgray', linewidth=0.5, zorder=1)

    # Create legend
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['Clinical Trial'], edgecolor='black',
                       label=f'Clinical Trial (n={len([s for s in studies_order if categorize_source(s) == "Clinical Trial"])})'),
        mpatches.Patch(facecolor=COLORS['In Vitro'], edgecolor='black',
                       label=f'In Vitro (n={len([s for s in studies_order if categorize_source(s) == "In Vitro"])})'),
        mpatches.Patch(facecolor=COLORS['Population'], edgecolor='black',
                       label=f'Population (n={len([s for s in studies_order if categorize_source(s) == "Population"])})')
    ]
    ax_legend.legend(handles=legend_elements, loc='center', frameon=True,
                     fontsize=8, title='Data Source Type', title_fontsize=9)
    ax_legend.axis('off')


def create_stacked_bar_detail(ax, df, mutations):
    """Create detailed stacked bar showing category breakdown per mutation"""

    # Prepare data
    resistance_mutations = [
        'Q67H', 'Q67K', 'Q67H+K70R', 'Q67H+N74D', 'Q67H+T107N', 'Q67H+N74S',
        'K70N', 'K70N+N74K', 'N74D', 'L56V', 'N57H', 'L56V+N57H',
        'M66I', 'M66I+A105T', 'M66I+T107A', 'M66I+N74D', 'M66I+N74D+A105T',
        'Q67K+K70H', 'K436E+I437T'
    ]
    mutations_in_data = [m for m in resistance_mutations if m in df['Mutation'].unique()]

    # Calculate category counts per mutation
    category_data = []
    for mut in mutations_in_data:
        mut_df = df[df['Mutation'] == mut]
        counts = {
            'mutation': mut,
            'total': len(mut_df),
            'Clinical Trial': len(mut_df[mut_df['source_category'] == 'Clinical Trial']),
            'In Vitro': len(mut_df[mut_df['source_category'] == 'In Vitro']),
            'Population': len(mut_df[mut_df['source_category'] == 'Population'])
        }
        category_data.append(counts)

    category_df = pd.DataFrame(category_data)
    category_df = category_df.sort_values('total', ascending=True)

    # Draw stacked bar
    y_pos = np.arange(len(category_df))
    left = np.zeros(len(category_df))

    for cat, color in [('Population', COLORS['Population']),
                        ('In Vitro', COLORS['In Vitro']),
                        ('Clinical Trial', COLORS['Clinical Trial'])]:
        heights = category_df[cat].values
        ax.barh(y_pos, heights, left=left, height=0.7,
                color=color, edgecolor='white', linewidth=0.5, label=cat)
        left = left + heights

    ax.set_yticks(y_pos)
    ax.set_yticklabels(category_df['mutation'], fontsize=8)
    ax.set_xlabel('N observations', fontsize=9)
    ax.grid(axis='x', alpha=0.3)

    # Color-code mutation labels
    for i, mut in enumerate(category_df['mutation']):
        if '+' in mut:
            ax.get_yticklabels()[i].set_color('#A23B72')


def create_availability_summary(ax, df):
    """Create a summary heatmap showing study × context tier availability"""

    # Create pivot
    pivot = df.pivot_table(
        index='study_source',
        columns='context_tier',
        values='FC',
        aggfunc='count',
        fill_value=0
    )

    # Reorder by data source category
    clinical = ['CAPELLA', 'JID2025']
    in_vitro = ['NATAP2022', 'PMC8092519', 'PMC12077089', 'PMC9600929',
                'PMID40231850', 'Structural', 'AAC', 'JAC2025']
    population = ['EMJ2025']

    study_order = clinical + in_vitro + population
    study_order = [s for s in study_order if s in pivot.index]

    pivot = pivot.reindex(study_order)

    # Create color mapping
    row_colors = pd.Series([categorize_source(s) for s in pivot.index], index=pivot.index).map(COLORS)

    sns.heatmap(pivot, annot=True, fmt='g', cmap='YlGnBu',
                cbar_kws={'label': 'N observations'}, ax=ax,
                linewidths=0.5, linecolor='gray')

    ax.set_xlabel('Context Tier', fontsize=9)
    ax.set_ylabel('Study', fontsize=9)
    ax.set_title('B. Study × Context Tier Availability', fontsize=11, weight='bold', loc='left')
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=7)
    plt.setp(ax.get_xticklabels(), rotation=0, fontsize=8)


def create_abundance_bar(ax, df):
    """Create top horizontal bar showing abundance of 3 data source types"""
    # Calculate abundance per source type
    counts = {
        'Clinical Trial': len(df[df['source_category'] == 'Clinical Trial']),
        'In Vitro': len(df[df['source_category'] == 'In Vitro']),
        'Population': len(df[df['source_category'] == 'Population'])
    }

    categories = ['Clinical Trial', 'In Vitro', 'Population']
    values = [counts[cat] for cat in categories]
    colors = [COLORS[cat] for cat in categories]
    total = sum(values)

    # Draw horizontal stacked bar
    x_pos = 0
    for i, (cat, val) in enumerate(zip(categories, values)):
        ax.barh(0, val, left=x_pos, height=0.5, color=colors[i],
                edgecolor='black', linewidth=0.5)
        # Add percentage label
        pct = val / total * 100
        if pct > 5:  # Only show if large enough
            ax.text(x_pos + val/2, 0, f'{pct:.1f}%',
                   ha='center', va='center', fontsize=9, fontweight='bold')
        x_pos += val

    # Set axis properties
    ax.set_ylim(-0.3, 0.3)
    ax.set_yticks([0])
    ax.set_yticklabels([f'N={total}'])
    ax.set_xlabel('N observations', fontsize=10)
    ax.set_title('Data Source Abundance', fontsize=11, weight='bold', loc='left')

    # Add legend
    legend_patches = [mpatches.Patch(color=colors[i], edgecolor='black', label=f'{cat} ({counts[cat]})')
                      for i, cat in enumerate(categories)]
    ax.legend(handles=legend_patches, loc='upper right', frameon=True, fontsize=8,
              title='Source Type', title_fontsize=9)

    # Remove y-axis spine for cleaner look
    ax.spines['left'].set_visible(False)


def create_figure1(df, studies_order, mutations, sorted_mutations, presence_matrix,
                   mutation_category_counts, study_category_counts, OUTPUT_DIR):
    """Generate Figure 1 with properly aligned study-centric layout

    Layout semantics:
    - Left stacked bars: each row = one mutation, bars show count by source type
    - Middle column: mutation names (text labels)
    - Bottom matrix: studies as columns (x), mutations as rows (y), dots/circles
    - Top bar: study-level summary aligned with matrix columns
    - Group indicator: colored blocks above matrix showing study categories

    CRITICAL: All columns must share the same x-coordinate system, and
    all rows must share the same y-coordinate system.
    """
    n_mutations = len(sorted_mutations)
    n_studies = len(studies_order)

    # Calculate max observation count for scaling
    max_obs = max(left.max() for left in [
        np.array([sum(mutation_category_counts[mut].values()) for mut in sorted_mutations])
    ])

    fig = plt.figure(figsize=(16, 10))

    # Define grid: 3 rows × 3 columns
    # Row 0: title (left), empty (middle), group indicator + study bars (right)
    # Row 1: left bars | mutation names | study matrix
    # Row 2: empty (for spacing)
    gs = fig.add_gridspec(3, 3,
                          width_ratios=[0.12, 0.04, 0.84],
                          height_ratios=[0.08, 0.84, 0.08],
                          left=0.06, right=0.97, top=0.95, bottom=0.08,
                          wspace=0.08, hspace=0.02)

    # ===== Get actual bounding boxes for alignment =====
    # We need to ensure the study columns in top and bottom align exactly

    # Top group indicator will be drawn inside the matrix axes using spans
    # Left bars will be drawn in their own axes but with matching y-coordinates

    # ===== Title area =====
    ax_title = fig.add_subplot(gs[0, 0])
    ax_title.text(0.5, 0.5, 'Evidence\nLandscape', fontsize=14, fontweight='bold',
                  ha='center', va='center', transform=ax_title.transAxes)
    ax_title.axis('off')

    # ===== Left: Stacked bars for each mutation =====
    ax_left = fig.add_subplot(gs[1, 0])

    y_positions = np.arange(n_mutations)
    bar_height = 0.65  # Slightly less than 1 to show gap between rows

    left = np.zeros(n_mutations)
    for category in ['Clinical Trial', 'In Vitro', 'Population']:
        heights = [mutation_category_counts[mut].get(category, 0) for mut in sorted_mutations]
        ax_left.barh(y_positions, heights, left=left, height=bar_height,
                    color=COLORS[category], edgecolor='white', linewidth=0.5)

        # Add count labels inside bars
        for idx, h in enumerate(heights):
            if h > 0 and h >= 2:  # Only label if segment is visible
                ax_left.text(left[idx] + h/2, y_positions[idx],
                            str(int(h)), ha='center', va='center',
                            fontsize=6, fontweight='bold', color='black')
        left = left + np.array(heights)

    ax_left.set_yticks([])
    ax_left.set_xlabel('N', fontsize=8)
    ax_left.set_xlim(0, max(4, left.max() * 1.3))
    ax_left.set_ylim(-0.5, n_mutations - 0.5)
    # NO invert_yaxis - match mutation name order (top to bottom = sorted_mutations order)

    # ===== CRITICAL: Add horizontal grid lines matching matrix rows =====
    # This creates visual continuity between left bars and right matrix
    for j in range(n_mutations):
        ax_left.axhline(j + 0.5, color='#e0e0e0', linewidth=0.3, zorder=0)

    # ===== Middle: Mutation names (text labels matching matrix rows) =====
    # Hide all y-axis labels from middle panel - keep only left panel labels
    ax_mid = fig.add_subplot(gs[1, 1])
    for i, mut in enumerate(sorted_mutations):
        ax_mid.text(0.5, i, mut, fontsize=7, ha='center', va='center',
                   fontweight='bold' if '+' in mut else 'normal',
                   color='#A23B72' if '+' in mut else 'black')
    ax_mid.set_ylim(-0.5, n_mutations - 0.5)
    # NO invert_yaxis - match mutation name order
    ax_mid.set_yticks([])
    ax_mid.set_yticklabels([])
    ax_mid.tick_params(axis='y', labelleft=False, labelright=False)
    ax_mid.tick_params(axis='x', labelbottom=False, labeltop=False)
    ax_mid.spines['top'].set_visible(False)
    ax_mid.spines['right'].set_visible(False)
    ax_mid.spines['bottom'].set_visible(False)
    ax_mid.spines['left'].set_visible(False)

    # ===== Right: Study × Mutation matrix =====
    ax_matrix = fig.add_subplot(gs[1, 2])

    # Draw matrix - studies as columns (x), mutations as rows (y)
    for i, study in enumerate(studies_order):
        study_category = categorize_source(study)
        for j, mut in enumerate(sorted_mutations):
            count = presence_matrix.loc[study, mut] if study in presence_matrix.index else 0
            if count > 0:
                # Draw filled circle (●) - size based on observation count
                size = 80 + count * 50
                ax_matrix.scatter(i, j, s=size, c=COLORS[study_category],
                                edgecolors='black', linewidth=0.8, zorder=3, marker='o')
            else:
                # Draw empty circle (○) for absence
                ax_matrix.scatter(i, j, s=25, c='white',
                                edgecolors='gray', linewidth=0.5, zorder=2, marker='o')

    ax_matrix.set_xticks(range(n_studies))
    ax_matrix.set_xticklabels(studies_order, rotation=45, ha='right', fontsize=8)
    ax_matrix.set_xlim(-0.5, n_studies - 0.5)
    ax_matrix.set_ylim(-0.5, n_mutations - 0.5)
    # NO invert_yaxis - match mutation name order (top to bottom = sorted_mutations order)
    ax_matrix.tick_params(axis='x', labelsize=8)
    ax_matrix.tick_params(axis='y', labelleft=False)
    ax_matrix.set_yticks([])  # Completely remove y-axis ticks and labels from right matrix panel

    # Add vertical grid lines
    for i in range(n_studies):
        ax_matrix.axvline(i + 0.5, color='lightgray', linewidth=0.5, zorder=1)

    # Reset ylim to only matrix area (no top bar chart)
    ax_matrix.set_ylim(-0.5, n_mutations - 0.5)

    # ===== Draw spanning background strip across ALL columns (left to right) =====
    # This connects the top bars visually to the entire width
    # Only draw horizontal lines within the MATRIX area (not in bar region)
    matrix_bottom_y = n_mutations - 0.5
    for j in range(n_mutations):
        row_center_y = j  # Row center positions are 0, 1, 2, ... n_mutations-1
        # Only draw grid line if it's within matrix bounds
        if -0.5 <= row_center_y + 0.5 <= matrix_bottom_y:
            for ax in [ax_left, ax_mid, ax_matrix]:
                ax.axhline(j + 0.5, color='#e0e0e0', linewidth=0.4, zorder=0)

    # Extend vertical grid lines from matrix through middle to left panel
    for i in range(n_studies):
        ax_matrix.axvline(i + 0.5, color='#d0d0d0', linewidth=0.5, zorder=0)
        ax_mid.axvline(0.5, color='#e0e0e0', linewidth=0.4, zorder=0)  # Middle column boundary

    # ===== Legend =====
    ax_legend = fig.add_axes([0.72, 0.04, 0.25, 0.05])
    ax_legend.axis('off')

    legend_elements = [
        mpatches.Patch(facecolor=COLORS['Clinical Trial'], edgecolor='black', linewidth=0.5,
                       label=f'Clinical Trial (n={len([s for s in studies_order if categorize_source(s) == "Clinical Trial"])})'),
        mpatches.Patch(facecolor=COLORS['In Vitro'], edgecolor='black', linewidth=0.5,
                       label=f'In Vitro (n={len([s for s in studies_order if categorize_source(s) == "In Vitro"])})'),
        mpatches.Patch(facecolor=COLORS['Population'], edgecolor='black', linewidth=0.5,
                       label=f'Population (n={len([s for s in studies_order if categorize_source(s) == "Population"])})')
    ]
    ax_legend.legend(handles=legend_elements, loc='center', frameon=True,
                    fontsize=7, ncol=3, title='Data Source Type', title_fontsize=8)

    # ===== Figure caption =====
    caption_text = (
        "a | Each dot represents a study reporting a mutation (● = present, ○ = absent). "
        "Dot color indicates data source type. Left bars show mutation observation counts by source. "
        "Combination mutations marked in purple. CT = Clinical Trial, Pop = Population."
    )
    fig.text(0.06, 0.03, caption_text, fontsize=7, ha='left', va='top',
             wrap=True, style='italic')

    # Save figure (use lower DPI first to avoid MemoryError, then save high-res separately)
    import gc
    gc.collect()
    plt.savefig(OUTPUT_DIR / "figure1_study_centric_upset.png", dpi=150, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "figure1_study_centric_upset.pdf", bbox_inches='tight')
    print(f"\n[OK] Saved: {OUTPUT_DIR / 'figure1_study_centric_upset.png'}")
    print(f"[OK] Saved: {OUTPUT_DIR / 'figure1_study_centric_upset.pdf'}")

    plt.close()


if __name__ == "__main__":
    print("="*60)
    print("Generating Figure 1: Evidence Landscape (Study-Centric)")
    print("="*60)

    # Load data
    df, studies_order, mutations = load_and_prepare_data()

    print(f"\nData loaded:")
    print(f"  - {len(studies_order)} studies")
    print(f"  - {len(mutations)} mutations")
    print(f"\nStudy categories:")
    for cat in ['Clinical Trial', 'In Vitro', 'Population']:
        n_studies = len([s for s in studies_order if categorize_source(s) == cat])
        print(f"  - {cat}: {n_studies} studies")

    # Build presence matrix
    presence_matrix = pd.DataFrame(0, index=studies_order, columns=mutations)
    for _, row in df.iterrows():
        study = row['study_source']
        mutation = row['Mutation']
        if study in presence_matrix.index and mutation in presence_matrix.columns:
            presence_matrix.loc[study, mutation] += 1

    # Sort mutations by total count (descending)
    total_counts = presence_matrix.sum(axis=0).sort_values(ascending=False)
    sorted_mutations = total_counts.index.tolist()

    # Calculate mutation category counts (UNIQUE studies per category)
    # This aligns semantically with right matrix which shows presence per study
    mutation_category_counts = {}
    for mut in sorted_mutations:
        mut_data = df[df['Mutation'] == mut]
        unique_studies_ct = mut_data['study_source'].nunique()
        counts = {
            'Clinical Trial': mut_data[mut_data['source_category'] == 'Clinical Trial']['study_source'].nunique(),
            'In Vitro': mut_data[mut_data['source_category'] == 'In Vitro']['study_source'].nunique(),
            'Population': mut_data[mut_data['source_category'] == 'Population']['study_source'].nunique()
        }
        mutation_category_counts[mut] = counts

    # Calculate study-level category counts
    study_category_counts = {}
    for study in studies_order:
        study_data = df[df['study_source'] == study]
        counts = {
            'Clinical Trial': len(study_data[study_data['source_category'] == 'Clinical Trial']),
            'In Vitro': len(study_data[study_data['source_category'] == 'In Vitro']),
            'Population': len(study_data[study_data['source_category'] == 'Population'])
        }
        study_category_counts[study] = counts

    # Generate figure with aligned layout
    create_figure1(df, studies_order, mutations, sorted_mutations, presence_matrix,
                   mutation_category_counts, study_category_counts, OUTPUT_DIR)
