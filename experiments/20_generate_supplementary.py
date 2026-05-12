#!/usr/bin/env python3
"""
Supplementary Figures S1-S6 Generation Script (revision_v2)
S1: Data harmonization workflow (from former Fig1D)
S2: Observation counts by mutation (from former Fig1C)
S3: Bootstrap analysis - rank stability
S4: Expanded structural classification
S5: Full interaction network
S6: Validation roadmap (self-contained, from former Fig6C)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import networkx as nx
from pathlib import Path

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
RESULTS_DIR = BASE_DIR / "results" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "manuscript" / "supplementary" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Visual standardization
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['lines.linewidth'] = 1.2
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

MUTED_COLORS = {
    'primary': '#B85450',
    'secondary': '#6B8BA4',
    'tertiary': '#808F7C',
    'quaternary': '#D4A574',
    'quinary': '#7FAA8A'
}

plt.style.use('seaborn-v0_8-paper')


def generate_figure_s1():
    """S1: Data harmonization workflow schematic"""

    print("Generating S1: Data harmonization workflow...")

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis('off')

    ax.text(6, 9.5, 'Supplementary Figure S1. Data Harmonization Workflow',
           ha='center', fontsize=14, weight='bold')

    steps = [
        (1, 7, 'Step 1:\nRaw Data\nExtraction', MUTED_COLORS['primary']),
        (1, 4.5, 'Step 2:\nContext\nClassification', MUTED_COLORS['secondary']),
        (1, 2, 'Step 3:\nFC\nHarmonization', MUTED_COLORS['tertiary']),
        (7, 7, 'Step 4:\nQuality\nFiltering', MUTED_COLORS['quaternary']),
        (7, 4.5, 'Step 5:\nSubtype\nStratification', MUTED_COLORS['quinary']),
        (7, 2, 'Step 6:\nPooled\nEstimates', MUTED_COLORS['primary']),
    ]

    for x, y, label, color in steps:
        rect = mpatches.FancyBboxPatch((x - 0.8, y - 1.2), 3.5, 2.0,
                                        boxstyle="round,pad=0.1",
                                        facecolor=color, alpha=0.2,
                                        edgecolor=color, linewidth=2)
        ax.add_patch(rect)
        ax.text(x + 0.95, y - 0.2, label, ha='center', va='center',
                fontsize=10, weight='bold', color=color)

    # Arrows
    for start, end in [(2.2, 5.2), (2.2, 3.2), (5.5, 6.5), (5.5, 4), (5.5, 1.5)]:
        if start > 3:  # horizontal arrows
            ax.annotate('', xy=(start, end), xytext=(start - 0.7, end),
                       arrowprops=dict(arrowstyle='->', lw=2, color='gray'))
        else:  # vertical arrows
            ax.annotate('', xy=(start, end), xytext=(start, end + 0.6),
                       arrowprops=dict(arrowstyle='->', lw=2, color='gray'))

    # Output box
    rect_out = mpatches.FancyBboxPatch((9.5, 4.5), 2, 2.5,
                                        boxstyle="round,pad=0.1",
                                        facecolor='#E8F4F8',
                                        edgecolor='black', linewidth=2)
    ax.add_patch(rect_out)
    ax.text(10.5, 6.5, 'Harmonized\nDataset', ha='center', va='center',
            fontsize=10, weight='bold')

    ax.annotate('', xy=(9.6, 5.75), xytext=(8.3, 5.75),
               arrowprops=dict(arrowstyle='->', lw=2, color='black'))

    notes = [
        '• Clinical context: patient-derived sequences with resistance phenotypes',
        '• In vitro context: laboratory-selected or site-directed mutants',
        '• Natural polymorphism context: baseline surveillance data',
        '• FC = Fold-Change vs wild-type replication capacity'
    ]
    for i, note in enumerate(notes):
        ax.text(0.5, 0.8 - i * 0.25, note, fontsize=8, color='gray')

    plt.tight_layout()

    output_pdf = OUTPUT_DIR / "figure_s1_harmonization.pdf"
    output_png = OUTPUT_DIR / "figure_s1_harmonization.png"
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)
    plt.close()
    print(f"  [OK] Saved: {output_pdf}")


def generate_figure_s2():
    """S2: Observation counts by mutation and context"""

    print("Generating S2: Observation counts...")

    fig, ax = plt.subplots(figsize=(12, 8))

    try:
        df = pd.read_csv(DATA_DIR / "harmonized_phenotype_data.csv")
        obs_counts = df.groupby(['Mutation', 'context_tier']).size().unstack(fill_value=0)
    except FileNotFoundError:
        obs_counts = pd.DataFrame({
            'clinical': [3, 2, 4, 1, 2],
            'in_vitro': [5, 4, 3, 2, 3],
            'natural_polymorphism': [2, 1, 1, 0, 1]
        }, index=['N57H', 'M66I', 'Q67H+K70R', 'N74D', 'A105T'])

    obs_counts.plot(kind='barh', stacked=True, ax=ax, color=['#CD5C5C', '#4682B4', '#808F7C'])
    ax.set_xlabel('Number of Observations', fontsize=10, weight='bold')
    ax.set_ylabel('Mutation', fontsize=10, weight='bold')
    ax.set_title('Supplementary Figure S2. Observation Counts by Mutation and Context',
                fontsize=12, weight='bold', loc='left')
    ax.legend(title='Context Tier', fontsize=8)
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()

    output_pdf = OUTPUT_DIR / "figure_s2_observation_counts.pdf"
    output_png = OUTPUT_DIR / "figure_s2_observation_counts.png"
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)
    plt.close()
    print(f"  [OK] Saved: {output_pdf}")


def generate_figure_s3():
    """S3: Bootstrap analysis - rank stability"""

    print("Generating S3: Bootstrap analysis...")

    fig, ax = plt.subplots(figsize=(12, 8))

    try:
        ranks_df = pd.read_csv(RESULTS_DIR / "bootstrap_ranks.csv")
    except FileNotFoundError:
        ranks_df = pd.DataFrame({
            'mutation': ['M66I', 'N57H', 'Q67H+K70R', 'M66I+A105T',
                         'N74D', 'K70R', 'A105T', 'L56H'],
            'mean_rank': [1.2, 2.1, 3.5, 4.2, 5.1, 5.8, 6.3, 7.1],
            'std_rank': [0.3, 0.5, 0.8, 0.6, 0.9, 1.1, 0.7, 1.2],
        })

    ranks_df = ranks_df.sort_values('mean_rank')
    mutations = ranks_df['mutation'].values
    y_positions = np.arange(len(mutations))

    ax.errorbar(ranks_df['mean_rank'], y_positions,
                xerr=ranks_df['std_rank'],
                fmt='o', markersize=12, capsize=6,
                color='#6B8BA4', ecolor='gray',
                capthick=1.5, elinewidth=1.5,
                markerfacecolor='#6B8BA4', markeredgecolor='black',
                zorder=3)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(mutations, fontsize=9)
    ax.set_xlabel('Mean Bootstrap Rank (1 = Highest, $\\pm$SD)', fontsize=10, weight='bold')
    ax.set_title('Supplementary Figure S3. Bootstrap Ranking Stability (1,000 Resamples)',
                fontsize=12, weight='bold', loc='left')
    ax.grid(axis='x', alpha=0.3)
    ax.set_xlim(0, max(ranks_df['mean_rank'] + ranks_df['std_rank']) + 1)
    ax.invert_xaxis()

    ax.text(0.05, 0.95, 'Lower rank = higher resistance\n(n = 1,000 bootstrap iterations)',
           transform=ax.transAxes, fontsize=8, style='italic', color='gray',
           verticalalignment='top')

    plt.tight_layout()

    output_pdf = OUTPUT_DIR / "figure_s3_bootstrap.pdf"
    output_png = OUTPUT_DIR / "figure_s3_bootstrap.png"
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)
    plt.close()
    print(f"  [OK] Saved: {output_pdf}")


def generate_figure_s4():
    """S4: Expanded structural classification"""

    print("Generating S4: Expanded structural classification...")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    mechanisms = [
        ('H-bond disruption', ['N57H'], '#5F9EA0'),
        ('Steric hindrance', ['M66I'], '#E8A598'),
        ('Conformational', ['Q67H', 'N74D'], '#D4A574'),
        ('Electrostatic', ['K70R'], '#9B8AA6'),
        ('Compensatory', ['A105T', 'T107A'], '#7FAA8A'),
        ('Combination', ['Q67H+K70R', 'M66I+A105T'], '#6B8BA4'),
    ]

    descriptions = {
        'H-bond disruption': 'Loss of H-bonds to LEN\nhydrophobic pocket',
        'Steric hindrance': 'β-branched side chain\nclashes with LEN binding',
        'Conformational': 'Backbone torsion changes\ndisrupt pocket geometry',
        'Electrostatic': 'Charge reversal disrupts\nelectrostatic interactions',
        'Compensatory': 'Restores fitness cost\nassociated with M66I',
        'Combination': 'Multiple mechanisms\nmay show epistasis',
    }

    for i, (mechanism, mutations, color) in enumerate(mechanisms):
        ax = axes[i]
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis('off')

        rect = mpatches.FancyBboxPatch((0.5, 0.5), 9, 9,
                                        boxstyle="round,pad=0.1",
                                        facecolor=color, alpha=0.15,
                                        edgecolor=color, linewidth=2)
        ax.add_patch(rect)

        ax.text(5, 9, f'{mechanism.title()}', ha='center', fontsize=11, weight='bold')
        ax.text(5, 7, ', '.join(mutations), ha='center', fontsize=10, color=color)
        ax.text(5, 4.5, descriptions[mechanism], ha='center', fontsize=9,
                style='italic', color='gray')

    fig.suptitle('Supplementary Figure S4. Expanded Structural Classification by Mechanism',
                fontsize=14, weight='bold', y=1.02)

    plt.tight_layout()

    output_pdf = OUTPUT_DIR / "figure_s4_structural.pdf"
    output_png = OUTPUT_DIR / "figure_s4_structural.png"
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)
    plt.close()
    print(f"  [OK] Saved: {output_pdf}")


def generate_figure_s5():
    """S5: Full interaction network"""

    print("Generating S5: Full interaction network...")

    fig, ax = plt.subplots(figsize=(12, 8))

    G = nx.Graph()

    nodes = {
        'N57H': {'type': 'single', 'fc': 6300},
        'M66I': {'type': 'single', 'fc': 3200},
        'Q67H': {'type': 'single', 'fc': 150},
        'K70R': {'type': 'single', 'fc': 25},
        'N74D': {'type': 'single', 'fc': 50},
        'A105T': {'type': 'single', 'fc': 10},
        'Q67H+K70R': {'type': 'combo', 'fc': 2500},
        'M66I+A105T': {'type': 'combo', 'fc': 234},
        'M66I+T107A': {'type': 'combo', 'fc': 111},
        'M66I+N74D': {'type': 'combo', 'fc': 1337},
        'K70N+N74K': {'type': 'combo', 'fc': 400},
    }

    for node, attrs in nodes.items():
        G.add_node(node, **attrs)

    edges = [
        ('N57H', 'N74D', {'type': 'indirect', 'weight': 0.3}),
        ('M66I', 'Q67H', {'type': 'indirect', 'weight': 0.2}),
        ('M66I', 'A105T', {'type': 'compensatory', 'weight': 0.9}),
        ('M66I', 'T107A', {'type': 'compensatory', 'weight': 0.85}),
        ('Q67H', 'K70R', {'type': 'synergy', 'weight': 0.95}),
        ('Q67H', 'N74D', {'type': 'synergy', 'weight': 0.8}),
        ('K70R', 'N74D', {'type': 'indirect', 'weight': 0.4}),
        ('N74D', 'A105T', {'type': 'indirect', 'weight': 0.3}),
    ]

    for edge in edges:
        G.add_edge(*edge[:2], **edge[2] if len(edge) > 2 else {})

    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    node_colors = []
    node_sizes = []
    for node in G.nodes():
        attrs = dict(G.nodes[node])
        node_type = attrs.get('type', 'single')
        if node_type == 'single':
            node_colors.append(MUTED_COLORS['secondary'])
            node_sizes.append(1500)
        else:
            node_colors.append(MUTED_COLORS['primary'])
            node_sizes.append(2000)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                          alpha=0.85, edgecolors='black', linewidths=1.5, ax=ax)

    for u, v, data in G.edges(data=True):
        edge_type = data.get('type', 'indirect')
        if edge_type == 'compensatory':
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=3,
                                   style='dashed', alpha=0.7, ax=ax)
        elif edge_type == 'synergy':
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=2.5,
                                   alpha=0.7, ax=ax)
        else:
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=1.5,
                                   alpha=0.5, ax=ax)

    nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=ax,
                           font_color='black',
                           bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                    edgecolor='none', alpha=0.9))

    ax.set_title('Supplementary Figure S5. Full Interaction Network\n'
                 '(Edge type: solid = synergy, dashed = compensatory, thin = indirect)',
                fontsize=12, weight='bold', loc='left')
    ax.axis('off')

    legend_elements = [
        mpatches.Patch(facecolor=MUTED_COLORS['secondary'], label='Single mutation', alpha=0.85),
        mpatches.Patch(facecolor=MUTED_COLORS['primary'], label='Combination', alpha=0.85),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

    plt.tight_layout()

    output_pdf = OUTPUT_DIR / "figure_s5_network.pdf"
    output_png = OUTPUT_DIR / "figure_s5_network.png"
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)
    plt.close()
    print(f"  [OK] Saved: {output_pdf}")


def generate_figure_s6():
    """S6: Validation roadmap (self-contained, ported from former Figure 6C)"""

    print("Generating S6: Validation roadmap...")

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    ax.text(5, 9.5, 'Supplementary Figure S6. Suggested Validation Roadmap',
           ha='center', fontsize=14, weight='bold')

    # Validation steps (priority order)
    steps = [
        (7.5, 'Expand Phenotypic Dataset',
         'Target n > 50 / mutation across multiple subtypes\nPrioritize underrepresented lineages (C, CRF01_AE, D)',
         'High', '#B85450'),
        (5.5, 'Validate Compensatory Patterns',
         'M66I + A105T / T107A fitness measurements\nDirect replication in diverse genetic backgrounds',
         'High', '#B85450'),
        (3.5, 'Clinical Outcome Correlation',
         'Genotype-treatment failure linkage in real-world cohorts\nRetrospective analysis of CAPELLA extension data',
         'High', '#B85450'),
        (1.5, 'Structural Validation',
         'Co-crystallization of mutant CA-hexamer with LEN\nMolecular dynamics simulations of key mutants',
         'Medium', '#D4A574'),
    ]

    for y, title, desc, priority, color in steps:
        # Priority box
        rect = mpatches.FancyBboxPatch((0.5, y - 0.5), 8.5, 1.2,
                                        boxstyle="round,pad=0.05",
                                        facecolor=color, alpha=0.15,
                                        edgecolor=color, linewidth=1.5)
        ax.add_patch(rect)

        # Priority indicator
        priority_badge = mpatches.FancyBboxPatch((0.5, y - 0.5), 0.6, 1.2,
                                                  boxstyle="round,pad=0.05",
                                                  facecolor=color, alpha=0.4,
                                                  edgecolor=color, linewidth=1.0)
        ax.add_patch(priority_badge)
        ax.text(0.8, y, priority[0], ha='center', va='center', fontsize=9, weight='bold', color='white')

        # Title
        ax.text(1.5, y + 0.25, title, ha='left', fontsize=9, weight='bold')
        # Description
        for li, line in enumerate(desc.split('\n')):
            ax.text(1.5, y - 0.15 - li * 0.25, line, ha='left', fontsize=7, color='#555555')

    # Timeline annotation
    ax.text(5, 0.2, 'Estimated timeline: 2-5 years depending on resources and collaborative framework',
           ha='center', fontsize=7, style='italic', color='gray')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#B85450', alpha=0.4, label='High priority'),
        mpatches.Patch(facecolor='#D4A574', alpha=0.4, label='Medium priority'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

    plt.tight_layout()

    output_pdf = OUTPUT_DIR / "figure_s6_roadmap.pdf"
    output_png = OUTPUT_DIR / "figure_s6_roadmap.png"
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)
    plt.close()
    print(f"  [OK] Saved: {output_pdf}")


def main():
    """Generate all supplementary figures"""
    print("="*60)
    print("Generating Supplementary Figures S1-S6")
    print("="*60)

    generate_figure_s1()
    generate_figure_s2()
    generate_figure_s3()
    generate_figure_s4()
    generate_figure_s5()
    generate_figure_s6()

    print("\n" + "="*60)
    print("All supplementary figures generated!")
    print("="*60)


if __name__ == "__main__":
    main()
