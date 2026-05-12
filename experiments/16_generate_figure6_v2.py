#!/usr/bin/env python3
"""
Figure 6: Claim Grading & Surveillance Framework (Improved Academic Design)
3 panels: claim matrix, surveillance framework schematic, validation roadmap
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches
from pathlib import Path

# Setup
BASE_DIR = Path(__file__).parent.parent.parent
RESULTS_DIR = BASE_DIR / "results" / "revision_v2"
OUTPUT_DIR = BASE_DIR / "manuscript" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

plt.style.use('seaborn-v0_8-paper')

def create_claim_matrix(ax):
    """Panel A: Claim grading matrix"""

    # Define claims with evidence strength
    claims = [
        ('Mutation-level ranking is stable', 'strong', 'Phenotype', 'Bootstrap CV, n=1000'),
        ('M66I/N57H are high resistance', 'strong', 'Phenotype', 'Multiple studies, n>5'),
        ('Subtype adds limited incremental value', 'moderate', 'Phenotype', 'Model comparison, AIC'),
        ('M66I causes steric hindrance', 'strong', 'Structure', 'Crystal structure, ΔΔG'),
        ('N57H loses H-bonds', 'strong', 'Structure', 'Crystal structure, ΔΔG'),
        ('Q67H+K70R is context-sensitive', 'moderate', 'Phenotype', '3 observations, 46x variation'),
        ('M66I+A105T is compensatory', 'hypothesis', 'Phenotype', '1 observation, needs validation'),
        ('Resistance positions are conserved', 'strong', 'Evolution', 'Conservation score >0.95'),
        ('Subtype baseline frequencies differ', 'moderate', 'Evolution', 'Los Alamos data, n>1000'),
        ('Fitness cost correlates with resistance', 'hypothesis', 'Evolution', 'Limited data, n=3'),
        ('Mutation-first surveillance is actionable', 'strong', 'Translation', 'Stable ranking, clear tiers'),
        ('Compensatory evolution is plausible', 'hypothesis', 'Translation', 'Indirect evidence only'),
    ]

    # Create dataframe
    df = pd.DataFrame(claims, columns=['claim', 'evidence', 'domain', 'basis'])

    # Create matrix for heatmap
    evidence_order = ['strong', 'moderate', 'hypothesis']
    domain_order = ['Phenotype', 'Structure', 'Evolution', 'Translation']

    # Count claims by domain and evidence
    matrix = np.zeros((len(domain_order), len(evidence_order)))
    for domain_idx, domain in enumerate(domain_order):
        for ev_idx, ev in enumerate(evidence_order):
            count = len(df[(df['domain'] == domain) & (df['evidence'] == ev)])
            matrix[domain_idx, ev_idx] = count

    # Plot heatmap
    sns.heatmap(matrix, annot=True, fmt='.0f', cmap='YlGn',
                cbar_kws={'label': 'Number of Claims'},
                linewidths=2, linecolor='white',
                xticklabels=['Strong', 'Moderate', 'Hypothesis-\nGenerating'],
                yticklabels=domain_order,
                ax=ax, vmin=0, vmax=5)

    ax.set_xlabel('Evidence Strength', fontsize=9, weight='bold')
    ax.set_ylabel('Domain', fontsize=9, weight='bold')
    ax.set_title('A. Claim Grading Matrix', fontsize=10, weight='bold', loc='left')
    plt.setp(ax.get_xticklabels(), rotation=0, fontsize=7)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=7)

    # Add summary text
    total_strong = int(matrix[:, 0].sum())
    total_moderate = int(matrix[:, 1].sum())
    total_hypothesis = int(matrix[:, 2].sum())

    summary_text = f'Total: {total_strong} strong, {total_moderate} moderate, {total_hypothesis} hypothesis-generating'
    ax.text(0.5, -0.12, summary_text, transform=ax.transAxes,
           ha='center', fontsize=6.5, weight='bold',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

def create_surveillance_framework(ax):
    """Panel B: Mutation-first surveillance framework schematic"""

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(5, 9.7, 'B. Mutation-First Surveillance Framework', ha='center',
           fontsize=10, weight='bold')

    # Color scheme
    colors = {
        'layer1': ('#E3F2FD', '#1976D2'),
        'layer2': ('#FFF3E0', '#F57C00'),
        'layer3': ('#E8F5E9', '#388E3C'),
        'layer4': ('#FFEBEE', '#D32F2F')
    }

    # Layer 1: Data collection
    box1 = FancyBboxPatch((0.5, 8.2), 9, 1.2, boxstyle="round,pad=0.08",
                          facecolor=colors['layer1'][0], edgecolor=colors['layer1'][1],
                          linewidth=2.5, zorder=2)
    ax.add_patch(box1)
    ax.text(5, 9.1, 'Layer 1: Clinical Sequence Surveillance', ha='center',
           fontsize=9, weight='bold', zorder=3)
    ax.text(5, 8.65, 'Genotype failures → Flag CA mutations', ha='center',
           fontsize=7, zorder=3)

    # Arrow
    arrow1 = FancyArrowPatch((5, 8.1), (5, 7.5), arrowstyle='->', mutation_scale=25,
                            linewidth=3, color='#333', zorder=1)
    ax.add_patch(arrow1)

    # Layer 2: Mutation-level triage
    box2 = FancyBboxPatch((0.5, 5.5), 9, 2, boxstyle="round,pad=0.08",
                          facecolor=colors['layer2'][0], edgecolor=colors['layer2'][1],
                          linewidth=2.5, zorder=2)
    ax.add_patch(box2)
    ax.text(5, 7.2, 'Layer 2: Mutation-Level Triage', ha='center',
           fontsize=9, weight='bold', zorder=3)

    # Three priority boxes (centered as a group around x=5.0)
    priority_box_width = 1.8
    priority_data = [
        (2.5, 6.1, 'HIGH', 'M66I, N57H', '(>1000×)', '#FFCDD2', '#C62828'),
        (5.0, 6.1, 'MEDIUM', 'Q67H, K70R', '(100-1000×)', '#FFE0B2', '#EF6C00'),
        (7.5, 6.1, 'LOW', 'Other', '(<100×)', '#BBDEFB', '#1565C0'),
    ]

    for x, y, level, muts, fc, bg_color, edge_color in priority_data:
        pbox = FancyBboxPatch((x - priority_box_width / 2, y), priority_box_width, 0.9, boxstyle="round,pad=0.05",
                             facecolor=bg_color, edgecolor=edge_color,
                             linewidth=2, alpha=0.9, zorder=3)
        ax.add_patch(pbox)
        ax.text(x, y+0.7, level, ha='center', fontsize=8, weight='bold',
               color=edge_color, zorder=4)
        ax.text(x, y+0.45, muts, ha='center', fontsize=6.5, zorder=4)
        ax.text(x, y+0.2, fc, ha='center', fontsize=6, zorder=4)

    # Arrow
    arrow2 = FancyArrowPatch((5, 5.4), (5, 4.8), arrowstyle='->', mutation_scale=25,
                            linewidth=3, color='#333', zorder=1)
    ax.add_patch(arrow2)

    # Layer 3: Context refinement
    box3 = FancyBboxPatch((0.5, 3.5), 9, 1.3, boxstyle="round,pad=0.08",
                          facecolor=colors['layer3'][0], edgecolor=colors['layer3'][1],
                          linewidth=2.5, zorder=2)
    ax.add_patch(box3)
    ax.text(5, 4.5, 'Layer 3: Context Refinement (Optional)', ha='center',
           fontsize=9, weight='bold', zorder=3)
    ax.text(5, 4.1, 'Check subtype | combinations | compensatory', ha='center',
           fontsize=7, zorder=3)
    ax.text(5, 3.75, 'Note: Subtype limited incremental value', ha='center',
           fontsize=6, style='italic', color='#666', zorder=3)

    # Arrow
    arrow3 = FancyArrowPatch((5, 3.4), (5, 2.8), arrowstyle='->', mutation_scale=25,
                            linewidth=3, color='#333', zorder=1)
    ax.add_patch(arrow3)

    # Layer 4: Action
    box4 = FancyBboxPatch((0.5, 1.5), 9, 1.3, boxstyle="round,pad=0.08",
                          facecolor=colors['layer4'][0], edgecolor=colors['layer4'][1],
                          linewidth=2.5, zorder=2)
    ax.add_patch(box4)
    ax.text(5, 2.5, 'Layer 4: Clinical Action', ha='center',
           fontsize=9, weight='bold', zorder=3)
    ax.text(5, 2.05, 'High: regimen change | Medium: monitor | Low: follow-up',
           ha='center', fontsize=7, zorder=3)

    # Removed side callout to keep the panel clean and avoid visual distraction.

def create_validation_roadmap(ax):
    """Panel C: Next-step validation roadmap"""

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(5, 9.7, 'C. Validation Roadmap', ha='center', fontsize=10, weight='bold')

    # Define validation steps
    steps = [
        (7.8, 'Expand Phenotypic Dataset', 'n>50/mutation, more subtypes', 'High', 'Ongoing'),
        (5.9, 'Validate Compensatory Patterns', 'M66I+A105T/T107A fitness', 'High', 'Needed'),
        (4.0, 'Clinical Outcome Correlation', 'Genotype-failure link', 'High', 'Needed'),
        (2.1, 'Structural Validation', 'Co-crystals, MD simulations', 'Medium', 'Partial'),
    ]

    for y, title, desc, priority, status in steps:
        # Priority color
        if priority == 'High':
            bg_color = '#FFCDD2'
            edge_color = '#C62828'
        elif priority == 'Medium':
            bg_color = '#FFE0B2'
            edge_color = '#EF6C00'
        else:
            bg_color = '#BBDEFB'
            edge_color = '#1565C0'

        # Status marker
        if status == 'Ongoing':
            marker = '✓'
            marker_color = '#388E3C'
        elif status == 'Partial':
            marker = '~'
            marker_color = '#F57C00'
        else:
            marker = '○'
            marker_color = '#D32F2F'

        # Draw box
        box = FancyBboxPatch((0.5, y-0.5), 7.5, 1.1, boxstyle="round,pad=0.08",
                            facecolor=bg_color, edgecolor=edge_color,
                            linewidth=2, alpha=0.4, zorder=2)
        ax.add_patch(box)

        # Add text
        ax.text(0.7, y+0.35, f'{marker} {title}', ha='left', fontsize=8,
               weight='bold', color=marker_color, zorder=3)
        ax.text(0.7, y-0.05, desc, ha='left', fontsize=6.5, zorder=3)

        # Priority label
        ax.text(8.3, y-0.05, priority, ha='center', va='center', fontsize=7,
               weight='bold', color=edge_color,
               bbox=dict(boxstyle='round', facecolor='white', edgecolor=edge_color,
                        linewidth=2, pad=0.3), zorder=3)

    # Add legend
    legend_y = 9.2
    ax.text(0.5, legend_y, 'Status:', ha='left', fontsize=7, weight='bold')
    ax.text(1.3, legend_y, '✓=Ongoing', ha='left', fontsize=6, color='#388E3C')
    ax.text(2.3, legend_y, '~=Partial', ha='left', fontsize=6, color='#F57C00')
    ax.text(3.2, legend_y, '○=Needed', ha='left', fontsize=6, color='#D32F2F')

    # Timeline annotation
    ax.text(5, 0.5, 'Timeline: 2-5 years', ha='center', fontsize=7,
           style='italic', weight='600',
           bbox=dict(boxstyle='round', facecolor='#FFF9C4', edgecolor='#F57F17',
                    linewidth=1.5, pad=0.4))

def main():
    """Generate Figure 6"""
    print("="*60)
    print("Generating Figure 6: Claim Grading & Framework (Improved Design)")
    print("="*60)

    # Layout: A/B/C in one horizontal row
    fig = plt.figure(figsize=(24, 7), constrained_layout=True)
    gs = GridSpec(1, 3, figure=fig, wspace=0.22)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # Generate panels
    print("\nPanel A: Claim matrix...")
    create_claim_matrix(ax1)

    print("Panel B: Surveillance framework...")
    create_surveillance_framework(ax2)

    print("Panel C: Validation roadmap...")
    create_validation_roadmap(ax3)

    # Save
    output_pdf = OUTPUT_DIR / "figure6_framework.pdf"
    output_png = OUTPUT_DIR / "figure6_framework.png"

    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)

    print(f"\n[OK] Saved: {output_pdf}")
    print(f"\n[OK] Saved: {output_png}")
    print("\n" + "="*60)
    print("Figure 6 complete!")
    print("="*60)

    plt.close()

if __name__ == "__main__":
    main()
