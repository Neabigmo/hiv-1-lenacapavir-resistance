#!/usr/bin/env python3
"""
Supplementary Figure S7: Sensitivity Analysis of Mutation Ranking Stability
3 panels:
  (left)   Full dataset (baseline ranking)
  (center) Remove single-observation mutations
  (right)  Subtype B only

Each panel shows mean log10FC rank per mutation; N57H and M66I
are highlighted to demonstrate top-rank robustness.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "manuscript" / "supplementary" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Visual standardization (matching main figures)
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['lines.linewidth'] = 1.2
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

HIGHLIGHT_COLOR = '#B85450'     # muted red for N57H/M66I
DEFAULT_COLOR = '#6B8BA4'       # blue-gray for others
STABLE_COLOR = '#7FAA8A'        # muted green for stable lower ranks

plt.style.use('seaborn-v0_8-paper')


def compute_ranking(df):
    """Compute per-mutation mean log10FC ranking."""
    grouped = df.groupby('Mutation')['log10_FC'].agg(['mean', 'count'])
    grouped = grouped[grouped['mean'].notna()]
    if len(grouped) == 0:
        return grouped
    grouped['rank'] = grouped['mean'].rank(ascending=False)
    grouped = grouped.sort_values('rank')
    return grouped


def draw_rank_panel(ax, ranking, title, highlight_mutations=None):
    """Draw a Cleveland-style dot plot of mutation ranks."""
    if highlight_mutations is None:
        highlight_mutations = {'N57H', 'M66I'}

    mutations = ranking.index.tolist()
    ranks = ranking['rank'].values
    y_pos = np.arange(len(mutations))

    # Colors: highlight key mutations
    colors = []
    for m in mutations:
        if m in highlight_mutations:
            colors.append(HIGHLIGHT_COLOR)
        else:
            colors.append(DEFAULT_COLOR)

    ax.scatter(ranks, y_pos, c=colors, s=100, edgecolors='black',
               linewidths=0.8, zorder=3, alpha=0.9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(mutations, fontsize=8)
    ax.set_xlabel('Mean Rank', fontsize=9)
    ax.set_title(title, fontsize=10, weight='bold', loc='left')

    # Invert x-axis so rank 1 is on the right (top resistance)
    ax.invert_xaxis()
    ax.set_xlim(len(mutations) + 1, 0.5)
    ax.grid(axis='x', alpha=0.3, zorder=1)

    # Annotate top rank value
    top_mut = mutations[0]
    top_rank = ranks[0]
    ax.annotate(f'{top_mut}: rank {int(top_rank)}',
                xy=(top_rank, 0), xytext=(top_rank - 0.3, 0.3),
                fontsize=6, color=HIGHLIGHT_COLOR, fontweight='bold',
                ha='center')


def main():
    print("=" * 60)
    print("Generating Figure S7: Sensitivity Analysis")
    print("=" * 60)

    df = pd.read_csv(DATA_DIR / "harmonized_phenotype_data.csv")
    df = df[df['log10_FC'].notna()].copy()
    print(f"  Total observations with valid log10FC: {len(df)}")

    # --- Compute rankings for each perturbation ---

    # 1. Full dataset (baseline)
    baseline = compute_ranking(df)
    print(f"  Baseline: {len(baseline)} mutations ranked")

    # 2. Remove single-observation mutations
    obs_counts = df.groupby('Mutation')['log10_FC'].count()
    mutations_with_multiple = obs_counts[obs_counts >= 2].index
    df_multi = df[df['Mutation'].isin(mutations_with_multiple)]
    multi_obs = compute_ranking(df_multi)
    n_removed_single = len(baseline) - len(multi_obs)
    print(f"  Single-obs removed: {len(multi_obs)} mutations ranked "
          f"({n_removed_single} single-obs mutations excluded)")

    # 3. Subtype B only
    # Subtype column contains values like 'B', 'Subtype_B', 'Mixed', 'Clinical_isolate'
    subtype_b_mask = (
        df['Subtype'].astype(str).str.upper().isin(['B', 'SUBTYPE_B'])
    ) | (
        df['Subtype_original'].astype(str).str.upper().isin(['B', 'SUBTYPE_B'])
    ) if 'Subtype_original' in df.columns else (
        df['Subtype'].astype(str).str.upper().isin(['B', 'SUBTYPE_B'])
    )

    df_subtype_b = df[subtype_b_mask].copy()

    subtype_b_rank = compute_ranking(df_subtype_b)
    print(f"  Subtype B only: {len(subtype_b_rank)} mutations from "
          f"{len(df_subtype_b)} observations")

    # --- Create figure ---
    fig = plt.figure(figsize=(16, 7))
    gs = GridSpec(1, 3, figure=fig, wspace=0.4)

    panels = [
        (fig.add_subplot(gs[0, 0]), baseline,
         'Full dataset\n(baseline)'),
        (fig.add_subplot(gs[0, 1]), multi_obs,
         f'Remove single-observation\nmutations (n={n_removed_single} removed)'),
        (fig.add_subplot(gs[0, 2]), subtype_b_rank,
         f'Subtype B only\n(n={len(df_subtype_b)} observations)'),
    ]

    for ax, ranking, title in panels:
        draw_rank_panel(ax, ranking, title)

    fig.suptitle('Supplementary Figure S7. Sensitivity Analysis of Mutation Ranking Stability',
                 fontsize=13, weight='bold', y=1.01)

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=HIGHLIGHT_COLOR,
               markersize=10, label='N57H / M66I (top-ranked)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=DEFAULT_COLOR,
               markersize=10, label='Other mutations'),
    ]
    fig.legend(handles=legend_elements, loc='lower center',
               ncol=2, fontsize=9, frameon=True)

    fig.subplots_adjust(top=0.88, bottom=0.18, wspace=0.4)

    output_pdf = OUTPUT_DIR / "figure_s7_sensitivity.pdf"
    output_png = OUTPUT_DIR / "figure_s7_sensitivity.png"

    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)

    print(f"\n[OK] Saved: {output_pdf}")
    print(f"[OK] Saved: {output_png}")
    print("\n" + "=" * 60)
    print("Figure S7 complete!")
    print("=" * 60)

    plt.close()


if __name__ == "__main__":
    main()
