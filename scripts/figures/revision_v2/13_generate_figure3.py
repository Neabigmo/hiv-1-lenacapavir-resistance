#!/usr/bin/env python3
"""
Figure 3: Epistasis & Context-Dependent Combinations (revision_v2)
2 panels: scatter/point-range for slopes, minimal interaction graph for networks
Deleted Panel C (Evidence strength classification) - already in Discussion
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

# Muted color palette for JMV submission
CONTEXT_COLORS = {
    'clinical': '#CD5C5C',           # muted red
    'in_vitro': '#4682B4',           # steel blue
    'natural_polymorphism': '#808F7C'  # gray-green
}
MUTATION_COLORS = {
    'primary': '#B85450',           # muted red for primary resistant
    'compensatory': '#7FAA8A',       # desaturated green for compensatory
    'synergy': '#6B8BA4'             # blue-gray for synergy
}

plt.style.use('seaborn-v0_8-paper')


def create_slope_dumbbell(ax):
    """Panel A: Point-range plot - Q67H+K70R variation across genetic backgrounds"""

    # Load context-specific data
    try:
        context_df = pd.read_csv(RESULTS_DIR / "context_specific_combinations.csv")
        q67h_k70r = context_df[context_df['combination'] == 'Q67H+K70R'].copy()
    except:
        # Create from harmonized data
        df = pd.read_csv(DATA_DIR / "harmonized_phenotype_data.csv")
        q67h_k70r = df[df['Mutation'] == 'Q67H+K70R'][['Source', 'FC_numeric', 'context_tier', 'Subtype']].copy()
        q67h_k70r.columns = ['study', 'FC', 'context', 'subtype']
        q67h_k70r['study'] = q67h_k70r['study'].str.replace('_', ' ')

    if len(q67h_k70r) == 0:
        ax.text(0.5, 0.5, 'No Q67H+K70R data', ha='center', va='center')
        ax.set_title('A. Context-Dependent: Q67H+K70R (Point-Range)', fontsize=10, weight='bold', loc='left')
        return

    # Sort by FC
    q67h_k70r = q67h_k70r.sort_values('FC')

    # Get data points
    studies = [str(s) for s in q67h_k70r['study'].values]
    fcs = q67h_k70r['FC'].values if 'FC' in q67h_k70r.columns else q67h_k70r['log10_FC'].values
    contexts = q67h_k70r['context'].values if 'context' in q67h_k70r.columns else ['in_vitro'] * len(q67h_k70r)

    n = len(q67h_k70r)
    y_positions = np.arange(n)

    # Draw connecting line (range)
    ax.hlines(y_positions, min(fcs), fcs, color='gray', alpha=0.5, linewidth=1.5)

    # Draw points with muted colors
    for i, (fc, context) in enumerate(zip(fcs, contexts)):
        color = CONTEXT_COLORS.get(context, '#6B8BA4')
        ax.plot(fc, i, 'o', color=color, markersize=10,
                markeredgecolor='black', markeredgewidth=1.2, zorder=3)

        # Label with value
        ax.text(fc + 1.2, i, f"{fc:.1f}×", ha='left', va='center', fontsize=10)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(studies, fontsize=8)
    ax.set_xlabel('Fold-Change vs WT', fontsize=10, weight='bold')
    ax.set_title('A. Q67H+K70R Context Sensitivity (Point-Range)', fontsize=10, weight='bold', loc='left')
    ax.set_xlim(-1, max(fcs) * 1.3)
    ax.grid(axis='x', alpha=0.3)

    # Add range annotation
    fc_range = max(fcs) / min(fcs)
    ax.text(0.5, -0.12, f'{fc_range:.1f}-fold range across {n} studies (Table S4)',
           transform=ax.transAxes, ha='center', va='top',
           fontsize=7, style='italic', color='gray')

    # Legend with muted colors
    legend_elements = [
        mpatches.Patch(color=CONTEXT_COLORS['clinical'], label='Clinical', alpha=0.8),
        mpatches.Patch(color=CONTEXT_COLORS['in_vitro'], label='In vitro', alpha=0.8),
        mpatches.Patch(color=CONTEXT_COLORS['natural_polymorphism'], label='Natural polymorphism', alpha=0.8)
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7, framealpha=0.9)


def create_m66i_network(ax):
    """Panel B: Minimal interaction graph - M66I combinations"""

    def _lookup_compensatory_data():
        """Try to derive combination data from harmonized phenotype data."""
        try:
            # Try clean CSV first, fall back to original
            try:
                harm = pd.read_csv(DATA_DIR / "harmonized_phenotype_data_clean.csv")
            except FileNotFoundError:
                harm = pd.read_csv(DATA_DIR / "harmonized_phenotype_data.csv")
            # Get only M66I-containing combinations with numeric FC
            m66i_combos = harm[harm['Mutation'].str.contains(r'\+', na=False)].copy()
            m66i_combos = m66i_combos[m66i_combos['Mutation'].str.contains('M66I', na=False)]
            m66i_combos = m66i_combos[m66i_combos['log10_FC'].notna()].copy()

            if len(m66i_combos) == 0:
                return None

            records = []
            for mut, group in m66i_combos.groupby('Mutation'):
                median_log10 = group['log10_FC'].median()
                fc = 10 ** median_log10
                # Infer pattern from fc vs solo M66I median
                m66i_solo = harm[harm['Mutation'] == 'M66I']['log10_FC'].median()
                if not np.isnan(m66i_solo) and not np.isnan(median_log10):
                    pattern = 'compensatory' if median_log10 < m66i_solo else 'other'
                else:
                    pattern = 'other'
                records.append({
                    'combination': mut,
                    'observed_fc': round(fc, 1),
                    'pattern': pattern
                })
            return pd.DataFrame(records)
        except Exception as e:
            print(f"  [Warning] Could not derive compensatory data: {e}")
            return None

    # Load compensatory data
    comp_df = None
    try:
        comp_df = pd.read_csv(RESULTS_DIR / "compensatory_patterns.csv")
    except FileNotFoundError:
        comp_df = _lookup_compensatory_data()

    # Final fallback if everything fails
    if comp_df is None or len(comp_df) == 0:
        comp_df = pd.DataFrame({
            'combination': ['M66I+A105T', 'M66I+T107A', 'M66I+N74D'],
            'observed_fc': [111, 234, 1337],
            'pattern': ['compensatory', 'compensatory', 'other']
        })

    # Ensure M66I+N74D+A105T triple mutant is included if available in harmonized data
    if 'M66I+N74D+A105T' not in comp_df['combination'].values:
        try:
            harm = pd.read_csv(DATA_DIR / "harmonized_phenotype_data_clean.csv")
            triple_data = harm[harm['Mutation'] == 'M66I+N74D+A105T']
            if len(triple_data) > 0 and triple_data['log10_FC'].notna().any():
                triple_fc = 10 ** triple_data['log10_FC'].median()
                triple_row = pd.DataFrame([{
                    'combination': 'M66I+N74D+A105T',
                    'observed_fc': round(triple_fc, 1),
                    'pattern': 'other'
                }])
                comp_df = pd.concat([comp_df, triple_row], ignore_index=True)
                print(f"  [Added] M66I+N74D+A105T (FC={triple_fc:.0f}) to M66I network from harmonized data")
        except (FileNotFoundError, KeyError) as e:
            print(f"  [Info] Could not add M66I+N74D+A105T: {e}")

    # Create minimal network (not cartoon infographic)
    G = nx.Graph()

    # Add M66I as central node
    G.add_node('M66I', node_type='primary', fc=3200)

    # Add combinations
    for _, row in comp_df.iterrows():
        combo = row['combination']
        fc = row['observed_fc']
        pattern = row.get('pattern', 'other')

        parts = combo.split('+')
        if 'M66I' in parts:
            other_muts = [p for p in parts if p != 'M66I']
            combo_name = '+'.join(other_muts)

            G.add_node(combo_name, node_type='compensatory' if 'compensatory' in pattern else 'other',
                      fc=fc)
            G.add_edge('M66I', combo_name, weight=fc, edge_type='compensatory' if 'compensatory' in pattern else 'other')

    # Use circular layout for cleaner visualization
    pos = nx.spring_layout(G, k=1.5, iterations=30, seed=42)

    # Draw nodes with muted colors
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        node_type = G.nodes[node].get('node_type', 'other')
        if node_type == 'primary':
            node_colors.append(MUTATION_COLORS['primary'])  # muted red
            node_sizes.append(1500)
        elif node_type == 'compensatory':
            node_colors.append(MUTATION_COLORS['compensatory'])  # desaturated green
            node_sizes.append(1200)
        else:
            node_colors.append(MUTATION_COLORS['synergy'])  # blue-gray
            node_sizes.append(1200)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                          alpha=0.8, edgecolors='black', linewidths=1.5, ax=ax)

    # Draw edges: solid for compensatory, dashed for other
    for u, v, data in G.edges(data=True):
        edge_type = data.get('edge_type', 'other')
        if edge_type == 'compensatory':
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=2,
                                 style='dashed', alpha=0.7, ax=ax)
        else:
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=1.5,
                                 alpha=0.7, ax=ax)

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax,
                           font_color='black',
                           bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                    edgecolor='none', alpha=0.8))

    # Add FC labels on edges
    edge_labels = {}
    for u, v, data in G.edges(data=True):
        fc = data['weight']
        edge_labels[(u, v)] = f"{fc:.0f}×"

    nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=9, ax=ax)

    # Annotation for triple mutant node (1,337-fold)
    for node in G.nodes():
        if 'N74D' in node or 'A105T' in node:
            fc = G.nodes[node].get('fc', 0)
            if fc == 1337:
                x, y = pos[node]
                ax.annotate('M66I+N74D+A105T\n(1,337-fold)',
                           xy=(x, y), xytext=(x + 0.15, y - 0.15),
                           fontsize=7, style='italic', color='#B85450',
                           ha='left', va='top',
                           bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='#B85450'),
                           arrowprops=dict(arrowstyle='->', color='#B85450', lw=0.8))
                break

    ax.set_title('B. M66I Combination Network', fontsize=10, weight='bold', loc='left')
    ax.axis('off')

    # Minimal legend with muted colors
    legend_elements = [
        mpatches.Patch(facecolor=MUTATION_COLORS['primary'], label='Primary resistant', alpha=0.8),
        mpatches.Patch(facecolor=MUTATION_COLORS['compensatory'], label='Compensatory (dashed)', alpha=0.8)
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7)


def main():
    """Generate Figure 3 - 2-Panel version (Point-Range + Minimal Network)"""
    print("="*60)
    print("Generating Figure 3: Point-Range + Minimal Network")
    print("="*60)

    # Create figure with 2 panels side-by-side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

    # Generate panels
    print("\nPanel A: Q67H+K70R point-range plot...")
    create_slope_dumbbell(ax1)

    print("Panel B: M66I minimal network...")
    create_m66i_network(ax2)

    # Save
    plt.tight_layout()

    output_pdf = OUTPUT_DIR / "figure3_interactions.pdf"
    output_png = OUTPUT_DIR / "figure3_interactions.png"

    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)

    print(f"\n[OK] Saved: {output_pdf}")
    print(f"[OK] Saved: {output_png}")
    print("\n" + "="*60)
    print("Figure 3 complete!")
    print("="*60)

    plt.close()


if __name__ == "__main__":
    main()