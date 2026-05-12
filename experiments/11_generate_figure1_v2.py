#!/usr/bin/env python3
"""
Figure 1: Evidence Landscape — PRISMA 2020 Flow Diagram + Data Availability Matrix
Panel A: PRISMA 2020 flow diagram (vertical layout, professional academic format)
Panel B: Mutation x Subtype x Context availability heatmap
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from pathlib import Path
import os, shutil, tempfile

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
RESULTS_DIR = BASE_DIR / "results" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "manuscript" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Global style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['axes.linewidth'] = 0.8

# Color palette
PRISMA_BLUE = '#6B8BA4'       # blue-gray for main flow boxes
PRISMA_BG = '#F4F7FA'         # very light blue-gray for main box fill
EXCL_RED = '#B85450'           # muted red for exclusion boxes
EXCL_BG = '#FFF8F8'            # very light pink for exclusion box fill
NOTE_GRAY = '#666666'          # gray for secondary notes
ARROW_MAIN = '#6B8BA4'         # arrow color for main flow
ARROW_EXCL = '#B85450'         # arrow color for exclusions


def create_prisma_flow(ax):
    """Panel A: PRISMA 2020 flow diagram with professional vertical layout.

    4 main boxes (Identification -> Screening -> Eligibility -> Included)
    are arranged vertically on the left. 3 exclusion boxes sit on the right
    at corresponding vertical levels with dashed borders and dashed arrows.
    """
    ax.set_xlim(0, 12.5)
    ax.set_ylim(0, 10.5)
    ax.axis('off')

    # Panel title
    ax.text(0.3, 10.1, 'A. PRISMA 2020 Flow Diagram',
            fontsize=14, weight='bold', ha='left', va='center')

    # ================================================================
    # MAIN FLOW BOXES (vertical, left side)
    # Box: center_x=3.5, width=4.4, height=1.6
    # ================================================================
    main_cx = 3.5
    main_w = 4.4
    main_h = 1.6

    # Each box: (center_y, title, count, note_under_count)
    main_boxes = [
        (8.2, 'Records identified',  'n = 87',  '82 PubMed + 5 other'),
        (6.0, 'Records screened',    'n = 78',  None),
        (3.8, 'Full-text assessed',  'n = 22',  None),
        (1.6, 'Studies included',    'n = 11',  '26 quantitative observations'),
    ]

    for cy, title, count, note in main_boxes:
        # Draw rounded rectangle
        box = FancyBboxPatch(
            (main_cx - main_w / 2, cy - main_h / 2),
            main_w, main_h,
            boxstyle="round,pad=0.1",
            facecolor=PRISMA_BG,
            edgecolor=PRISMA_BLUE,
            linewidth=2.2,
            zorder=2,
        )
        ax.add_patch(box)

        # Title line (smaller, bold)
        ax.text(main_cx, cy + 0.25, title,
                ha='center', va='center',
                fontsize=10, weight='bold', zorder=3)

        # Count line (larger, bold)
        ax.text(main_cx, cy - 0.18, count,
                ha='center', va='center',
                fontsize=14, weight='bold', zorder=3)

        # Note line (italic, gray)
        if note:
            ax.text(main_cx, cy - 0.55, note,
                    ha='center', va='center',
                    fontsize=8, style='italic', color=NOTE_GRAY, zorder=3)

    # ================================================================
    # VERTICAL ARROWS connecting main boxes
    # ================================================================
    main_bottoms = [7.4, 5.2, 3.0]
    main_tops    = [6.8, 4.6, 2.4]

    for bot, top in zip(main_bottoms, main_tops):
        ax.annotate(
            '', xy=(main_cx, top), xytext=(main_cx, bot),
            arrowprops=dict(
                arrowstyle='->', color=ARROW_MAIN,
                lw=2.8, mutation_scale=35,
            ),
        )

    # ================================================================
    # EXCLUSION BOXES (right side, dashed borders)
    # Box: center_x=9.2, width=3.0
    # ================================================================
    excl_cx = 9.5
    excl_w = 3.2

    # ---- E1: Duplicates ----
    e1_cy = 8.2
    e1_h = 1.2
    e1_box = FancyBboxPatch(
        (excl_cx - excl_w / 2, e1_cy - e1_h / 2),
        excl_w, e1_h,
        boxstyle="round,pad=0.06",
        facecolor=EXCL_BG,
        edgecolor=EXCL_RED,
        linewidth=1.5,
        linestyle='--',
        zorder=2,
    )
    ax.add_patch(e1_box)
    ax.text(excl_cx, e1_cy + 0.15, 'Duplicates removed',
            ha='center', va='center', fontsize=9, zorder=3)
    ax.text(excl_cx, e1_cy - 0.25, 'n = 9',
            ha='center', va='center', fontsize=12, weight='bold',
            color=EXCL_RED, zorder=3)

    # ---- E2: Title/abstract excluded ----
    e2_cy = 6.0
    e2_h = 1.2
    e2_box = FancyBboxPatch(
        (excl_cx - excl_w / 2, e2_cy - e2_h / 2),
        excl_w, e2_h,
        boxstyle="round,pad=0.06",
        facecolor=EXCL_BG,
        edgecolor=EXCL_RED,
        linewidth=1.5,
        linestyle='--',
        zorder=2,
    )
    ax.add_patch(e2_box)
    ax.text(excl_cx, e2_cy + 0.15, 'Title/abstract excluded',
            ha='center', va='center', fontsize=9, zorder=3)
    ax.text(excl_cx, e2_cy - 0.25, 'n = 56',
            ha='center', va='center', fontsize=12, weight='bold',
            color=EXCL_RED, zorder=3)

    # ---- E3: Full-text excluded (taller to accommodate reasons) ----
    e3_cy = 3.8
    e3_h = 1.8
    e3_box = FancyBboxPatch(
        (excl_cx - excl_w / 2, e3_cy - e3_h / 2),
        excl_w, e3_h,
        boxstyle="round,pad=0.06",
        facecolor=EXCL_BG,
        edgecolor=EXCL_RED,
        linewidth=1.5,
        linestyle='--',
        zorder=2,
    )
    ax.add_patch(e3_box)
    ax.text(excl_cx, e3_cy + 0.55, 'Full-text excluded',
            ha='center', va='center', fontsize=9, zorder=3)
    ax.text(excl_cx, e3_cy + 0.10, 'n = 11',
            ha='center', va='center', fontsize=12, weight='bold',
            color=EXCL_RED, zorder=3)
    ax.text(excl_cx, e3_cy - 0.35, '8 no quantitative data',
            ha='center', va='center', fontsize=8.5,
            color=NOTE_GRAY, zorder=3)
    ax.text(excl_cx, e3_cy - 0.65, '3 HIV-2 only',
            ha='center', va='center', fontsize=8.5,
            color=NOTE_GRAY, zorder=3)

    # ================================================================
    # DASHED ARROWS from main boxes to exclusion boxes
    # ================================================================
    # From right edge of main box to left edge of exclusion box
    excl_arrow_head = excl_cx - excl_w / 2  # 7.9
    main_right = main_cx + main_w / 2        # 5.7

    excl_centers = [8.2, 6.0, 3.8]
    for cy in excl_centers:
        ax.annotate(
            '', xy=(excl_arrow_head, cy), xytext=(main_right + 0.05, cy),
            arrowprops=dict(
                arrowstyle='->', color=ARROW_EXCL,
                lw=1.5, linestyle='dashed', mutation_scale=20,
            ),
        )


def create_availability_matrix(ax):
    """Panel B: Mutation x Subtype x Context availability heatmap."""
    avail_df = pd.read_csv(DATA_DIR / "availability_matrix.csv")

    # Create pivot table
    avail_df['subtype_context'] = avail_df['subtype'] + '_' + avail_df['context']

    pivot = avail_df.pivot_table(
        index='mutation',
        columns='subtype_context',
        values='n_observations',
        fill_value=0,
    )

    # Sort: single mutations first, then combination mutations
    mutation_order = pivot.index.tolist()
    single_muts = [m for m in mutation_order if '+' not in m]
    combo_muts = [m for m in mutation_order if '+' in m]
    sorted_order = single_muts + combo_muts
    pivot = pivot.reindex(sorted_order)

    # Add row sum column and move it to leftmost
    pivot['Sum'] = pivot.sum(axis=1).astype(int)
    sum_col = pivot.pop('Sum')
    pivot.insert(0, 'Sum', sum_col)

    # Heatmap using Blues colormap
    sns.heatmap(
        pivot, annot=True, fmt='g', cmap='Blues',
        cbar_kws={'label': 'N observations'},
        linewidths=0.5, linecolor='gray',
        ax=ax, vmin=0, vmax=None,
    )

    # Separator between single and combination mutations
    if single_muts and combo_muts:
        n_single = len(single_muts)
        ax.axhline(y=n_single, color='black', linewidth=1.5)

    # Separator after Sum column (first column)
    ax.axvline(x=1, color='#333333', linewidth=1.5, linestyle='-')

    ax.set_xlabel('Subtype x Context', fontsize=9)
    ax.set_ylabel('Mutation', fontsize=9)
    ax.set_title('B. Data Availability Matrix',
                 fontsize=12, weight='bold', loc='left', pad=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)

    # Mutation group labels
    if single_muts and combo_muts:
        n_single = len(single_muts)
        n_combo = len(combo_muts)
        ax.text(
            -0.35, n_single / 2 - 0.5, 'Single\nmutants',
            ha='center', va='center', fontsize=7, color='#2E6B9E',
            weight='bold', transform=ax.get_yaxis_transform(),
        )
        ax.text(
            -0.35, n_single + n_combo / 2 - 0.5, 'Combination\nmutants',
            ha='center', va='center', fontsize=7, color='#2E6B9E',
            weight='bold', transform=ax.get_yaxis_transform(),
        )


def main():
    """Generate Figure 1 — 2 panels: PRISMA flow + Availability Matrix."""
    print('=' * 60)
    print('Generating Figure 1: Evidence Landscape')
    print('  Panel A: PRISMA 2020 Flow Diagram')
    print('  Panel B: Data Availability Matrix')
    print('=' * 60)

    # Figure layout: 1 x 2 grid
    fig = plt.figure(figsize=(16, 8.5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.30)

    # Panel A: PRISMA flow diagram
    print('\n  [Panel A] Drawing PRISMA 2020 flow diagram...')
    ax_a = fig.add_subplot(gs[0, 0])
    create_prisma_flow(ax_a)

    # Panel B: Availability matrix
    print('  [Panel B] Drawing data availability matrix...')
    ax_b = fig.add_subplot(gs[0, 1])
    create_availability_matrix(ax_b)

    # Save (use temp to avoid permission conflicts with open file viewers)
    tmp = Path(tempfile.mkdtemp())
    tmp_pdf = tmp / 'figure1.pdf'
    tmp_png = tmp / 'figure1.png'

    fig.savefig(str(tmp_pdf), dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(str(tmp_png), dpi=300, bbox_inches='tight', pad_inches=0.3)

    output_pdf = OUTPUT_DIR / 'figure1_evidence_landscape_v4.pdf'
    output_png = OUTPUT_DIR / 'figure1_evidence_landscape_v4.png'

    shutil.move(str(tmp_pdf), str(output_pdf))
    shutil.move(str(tmp_png), str(output_png))
    shutil.rmtree(str(tmp), ignore_errors=True)

    print(f'\n  [OK] Saved: {output_pdf}')
    print(f'  [OK] Saved: {output_png}')
    print('=' * 60)
    print('Figure 1 complete!')
    print('=' * 60)

    plt.close()


if __name__ == '__main__':
    main()
