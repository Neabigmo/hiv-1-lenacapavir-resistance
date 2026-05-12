#!/usr/bin/env python3
"""
Figure 4: Structural Coherence of LEN Resistance Mutations (Revision)
3-panel layout:
  A (top, full width): Hexamer structural localization (muted colors)
  B (bottom left): 2x2 mechanism mini-panels
  C (bottom right): Structure-phenotype scatter correlation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from pathlib import Path
from scipy import stats

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
DATA_DIR = BASE_DIR / "data" / "processed" / "revision_v2"
RESULTS_DIR = BASE_DIR / "results" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "manuscript" / "figures"
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

# Muted colors
RESIDUE_COLORS = {
    'N57': '#5F9EA0',     # muted teal
    'M66': '#E8A598',     # muted salmon
    'Q67': '#D4A574',     # muted ochre
    'N74': '#6A7B8B',     # slate blue
    'K70': '#9B8AA6',     # muted purple
    'A105': '#7FAA8A',    # desaturated green
    'T107': '#7FAA8A',    # desaturated green
}
MECHANISM_COLORS = {
    'hbond': '#5F9EA0',
    'steric': '#E8A598',
    'synergy': '#D4A574',
    'compensatory': '#7FAA8A'
}

plt.style.use('seaborn-v0_8-paper')


def create_panel_a(ax):
    """Panel A: PyMOL structure render of CA hexamer with mutation sites"""

    # Load user's PyMOL render
    img_path = BASE_DIR / "manuscript" / "figures" / "pymol- structure.png"
    if not img_path.exists():
        raise FileNotFoundError(
            f"PyMOL render not found at {img_path}. "
            "Generate it via PyMOL before running this script."
        )

    img = plt.imread(img_path)
    ax.imshow(img)
    ax.axis('off')
    ax.set_title('A. Structural Localization of Major LEN Resistance Positions',
                  fontsize=12, weight='bold', loc='left')


def create_panel_b(ax):
    """Panel B: 2x2 mechanism mini-panels"""

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis('off')

    ax.text(5, 5.85, 'B. Resistance Mechanisms', ha='center', fontsize=12, weight='bold')

    # 2x2 grid positions
    mech_positions = [
        (2.7, 3.05, 2.1, 2.5, 'hbond', 'N57H', 'H-bond\ndisruption', '2 H-bonds lost\nHydrophobic pocket'),
        (5.2, 3.05, 2.1, 2.5, 'steric', 'M66I', 'Steric\nhindrance', 'β-branched clash\nBinding interference'),
        (2.7, 0.45, 2.1, 2.4, 'synergy', 'Q67H+N74D', 'Conform./electrostatic\nsynergy', 'Dual mechanism\nNo directional epistasis'),
        (5.2, 0.45, 2.1, 2.4, 'compensatory', 'A105T/T107A', 'Putative\ncompensation', 'Putative fitness-restoring pattern\nPending validation'),
    ]

    for x, y, w, h, mech, label, title, desc in mech_positions:
        color = MECHANISM_COLORS[mech]
        region = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                               facecolor=color, alpha=0.15,
                               edgecolor=color, linewidth=1.5)
        ax.add_patch(region)
        ax.text(x + w / 2, y + h - 0.3, label, ha='center', fontsize=10, weight='bold', color=color)
        # Mechanism title
        for li, line in enumerate(title.split('\n')):
            ax.text(x + w / 2, y + h - 1.0 - li * 0.3, line, ha='center', fontsize=7, style='italic', color='#333333')
        # Description
        for li, line in enumerate(desc.split('\n')):
            ax.text(x + w / 2, y + 0.5 - li * 0.25, line, ha='center', fontsize=6.5, color='#555555')


def create_panel_c(ax):
    """Panel C: Structure-phenotype scatter correlation"""

    # Data
    mutations_data = {
        'N57H': {'perturbation': 7.2, 'log10fc': 3.69},
        'M66I': {'perturbation': 8.5, 'log10fc': 3.51},
        'Q67H': {'perturbation': 6.0, 'log10fc': 2.18},
        'A105T': {'perturbation': 4.0, 'log10fc': 1.0},
        'Q67H+N74D': {'perturbation': 5.5, 'log10fc': 2.18},
        'N74D': {'perturbation': 3.5, 'log10fc': 1.70},
        'K70R': {'perturbation': 3.0, 'log10fc': 1.40},
    }

    color_map = {
        'N57H': MECHANISM_COLORS['hbond'],
        'M66I': MECHANISM_COLORS['steric'],
        'Q67H': MECHANISM_COLORS['synergy'],
        'A105T': MECHANISM_COLORS['compensatory'],
        'Q67H+N74D': MECHANISM_COLORS['synergy'],
        'N74D': '#6A7B8B',
        'K70R': '#9B8AA6',
    }

    perturbations = [d['perturbation'] for d in mutations_data.values()]
    log10fcs = [d['log10fc'] for d in mutations_data.values()]
    labels = list(mutations_data.keys())
    colors = [color_map.get(l, '#6B8BA4') for l in labels]

    # Scatter
    ax.scatter(perturbations, log10fcs, s=120, c=colors, edgecolors='black',
               linewidth=1.2, alpha=0.85, zorder=3)

    # Labels
    for i, (x, y, label) in enumerate(zip(perturbations, log10fcs, labels)):
        ax.annotate(label, (x, y), xytext=(5, 5),
                   textcoords='offset points', fontsize=7, fontweight='bold')

    # Regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(perturbations, log10fcs)
    x_line = np.array([min(perturbations) - 0.5, max(perturbations) + 0.5])
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, '--', color='gray', linewidth=1.2, alpha=0.7, zorder=2)

    # Statistics box
    ax.text(0.05, 0.95, f'r = {r_value:.3f}\np = {p_value:.3f}',
            transform=ax.transAxes, fontsize=9,
            verticalalignment='top', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.9, edgecolor='gray'))

    ax.set_xlabel('Structural Perturbation Score', fontsize=9, weight='bold')
    ax.set_ylabel('log$_{10}$(FC)', fontsize=9, weight='bold')
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3, zorder=1)
    ax.set_title('C. Structure-Phenotype Correlation', fontsize=11, weight='bold', loc='left')

    ax.axhline(y=3, color='#CD5C5C', linestyle='--', alpha=0.3, linewidth=0.8)
    ax.axhline(y=2, color='gray', linestyle='--', alpha=0.3, linewidth=0.8)


def main():
    """Generate Figure 4 - Structural coherence (3 panels)"""
    print("="*60)
    print("Generating Figure 4: Structural Coherence (3-panel layout)")
    print("="*60)

    fig = plt.figure(figsize=(16, 14))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.8, 1.0], hspace=0.2, wspace=0.3)

    ax1 = fig.add_subplot(gs[0, :])    # A: hexamer (full width)
    ax2 = fig.add_subplot(gs[1, 0])    # B: 2x2 mechanisms
    ax3 = fig.add_subplot(gs[1, 1])    # C: scatter correlation

    print("\nPanel A: Hexamer structural localization...")
    create_panel_a(ax1)

    print("Panel B: 2x2 mechanism panels...")
    create_panel_b(ax2)

    print("Panel C: Structure-phenotype correlation...")
    create_panel_c(ax3)

    output_pdf = OUTPUT_DIR / "figure4_structure.pdf"
    output_png = OUTPUT_DIR / "figure4_structure.png"

    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)

    print(f"\n[OK] Saved: {output_pdf}")
    print(f"[OK] Saved: {output_png}")
    print("\n" + "="*60)
    print("Figure 4 complete!")
    print("="*60)

    plt.close()


if __name__ == "__main__":
    main()
