#!/usr/bin/env python3
"""
Figure 6: Claim Grading & Surveillance Framework
3 panels: claim matrix, surveillance framework schematic, validation roadmap
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from pathlib import Path

# Setup
BASE_DIR = Path(__file__).parent.parent.parent.parent
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
    ax.text(0.5, -0.18, summary_text, transform=ax.transAxes,
           ha='center', fontsize=7, weight='bold',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

def create_surveillance_framework(ax):
    """Panel B: Mutation-first surveillance framework schematic"""

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(5, 9.5, 'B. Mutation-First Surveillance Framework', ha='center', fontsize=10, weight='bold')

    # Layer 1: Data collection
    rect1 = mpatches.FancyBboxPatch((0.5, 7.5), 9, 1.2,
                                     boxstyle="round,pad=0.1",
                                     facecolor='#E8F4F8', edgecolor='black', linewidth=2)
    ax.add_patch(rect1)
    ax.text(5, 8.5, 'Layer 1: Clinical Sequence Surveillance', ha='center', fontsize=8, weight='bold')
    ax.text(5, 8, 'Genotype failures | Flag CA mutations',
           ha='center', fontsize=6)

    # Arrow
    ax.arrow(5, 7.4, 0, -0.5, head_width=0.3, head_length=0.2, fc='black', ec='black', linewidth=2)

    # Layer 2: Mutation-level triage
    rect2 = mpatches.FancyBboxPatch((0.5, 5), 9, 1.8,
                                     boxstyle="round,pad=0.1",
                                     facecolor='#FFF4E6', edgecolor='black', linewidth=2)
    ax.add_patch(rect2)
    ax.text(5, 6.5, 'Layer 2: Mutation-Level Triage', ha='center', fontsize=8, weight='bold')

    # Three priority boxes
    priority_boxes = [
        (1.5, 5.3, 'High\nM66I, N57H\n(>1000x)', '#E74C3C'),
        (4.5, 5.3, 'Medium\nQ67H, K70R\n(100-1000x)', '#F39C12'),
        (7.5, 5.3, 'Low\nOther\n(<100x)', '#3498DB'),
    ]

    for x, y, label, color in priority_boxes:
        box = mpatches.FancyBboxPatch((x-0.8, y), 1.6, 0.9,
                                       boxstyle="round,pad=0.05",
                                       facecolor=color, edgecolor='black',
                                       linewidth=1.5, alpha=0.7)
        ax.add_patch(box)
        ax.text(x, y+0.45, label, ha='center', va='center', fontsize=6, weight='bold')

    # Arrow
    ax.arrow(5, 4.9, 0, -0.5, head_width=0.3, head_length=0.2, fc='black', ec='black', linewidth=2)

    # Layer 3: Context refinement
    rect3 = mpatches.FancyBboxPatch((0.5, 2.5), 9, 1.3,
                                     boxstyle="round,pad=0.1",
                                     facecolor='#E8F5E9', edgecolor='black', linewidth=2)
    ax.add_patch(rect3)
    ax.text(5, 3.6, 'Layer 3: Context Refinement (Optional)', ha='center', fontsize=8, weight='bold')
    ax.text(5, 3.1, 'Check subtype | combinations | compensatory',
           ha='center', fontsize=6)
    ax.text(5, 2.7, 'Note: Subtype limited incremental value',
           ha='center', fontsize=5, style='italic', color='gray')

    # Arrow
    ax.arrow(5, 2.4, 0, -0.5, head_width=0.3, head_length=0.2, fc='black', ec='black', linewidth=2)

    # Layer 4: Action
    rect4 = mpatches.FancyBboxPatch((0.5, 0.5), 9, 1.2,
                                     boxstyle="round,pad=0.1",
                                     facecolor='#FFE6E6', edgecolor='red', linewidth=2)
    ax.add_patch(rect4)
    ax.text(5, 1.5, 'Layer 4: Clinical Action', ha='center', fontsize=8, weight='bold')
    ax.text(5, 1, 'High: regimen change | Medium: monitor | Low: follow-up',
           ha='center', fontsize=6)

    # Add side annotation
    ax.text(9.8, 5, 'Mutation signal\nstable across\ncontexts',
           ha='center', va='center', fontsize=5, rotation=-90,
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

def create_validation_roadmap(ax):
    """Panel C: Next-step validation roadmap"""

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(5, 9.5, 'C. Validation Roadmap', ha='center', fontsize=10, weight='bold')

    # Define validation steps (reduced to 4 most critical)
    steps = [
        # (y_pos, title, description, priority, status)
        (7.5, 'Expand Phenotypic Dataset', 'n>50/mutation, more subtypes', 'High', 'Ongoing'),
        (5.5, 'Validate Compensatory Patterns', 'M66I+A105T/T107A fitness', 'High', 'Needed'),
        (3.5, 'Clinical Outcome Correlation', 'Genotype-failure link', 'High', 'Needed'),
        (1.5, 'Structural Validation', 'Co-crystals, MD sims', 'Medium', 'Partial'),
    ]

    for y, title, desc, priority, status in steps:
        # Priority color
        if priority == 'High':
            color = '#E74C3C'
        elif priority == 'Medium':
            color = '#F39C12'
        else:
            color = '#3498DB'

        # Status marker
        if status == 'Ongoing':
            marker = '[OK]'
            marker_color = 'green'
        elif status == 'Partial':
            marker = '~'
            marker_color = 'orange'
        else:
            marker = '[ ]'
            marker_color = 'red'

        # Draw box
        box = mpatches.FancyBboxPatch((0.5, y-0.4), 8, 1.0,
                                       boxstyle="round,pad=0.05",
                                       facecolor=color, edgecolor='black',
                                       linewidth=1.5, alpha=0.3)
        ax.add_patch(box)

        # Add text
        ax.text(0.7, y+0.25, f'{marker} {title}', ha='left', fontsize=7, weight='bold', color=marker_color)
        ax.text(0.7, y-0.15, desc, ha='left', fontsize=6)

        # Priority label
        ax.text(8.8, y, priority, ha='center', va='center', fontsize=6, weight='bold',
               bbox=dict(boxstyle='round', facecolor='white', edgecolor=color, linewidth=2))

    # Add legend
    legend_y = 9
    ax.text(0.5, legend_y, 'Status:', ha='left', fontsize=6, weight='bold')
    ax.text(1.3, legend_y, '[OK]=Ongoing', ha='left', fontsize=5, color='green')
    ax.text(2.5, legend_y, '~=Partial', ha='left', fontsize=5, color='orange')
    ax.text(3.5, legend_y, '[ ]=Needed', ha='left', fontsize=5, color='red')

    # Add timeline annotation
    ax.text(5, 0.3, 'Timeline: 2-5 years',
           ha='center', fontsize=6, style='italic',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

def main():
    """Generate Figure 6"""
    print("="*60)
    print("Generating Figure 6: Claim Grading & Framework")
    print("="*60)

    # Create figure with increased size and spacing
    fig = plt.figure(figsize=(22, 17))
    gs = GridSpec(3, 1, figure=fig, hspace=0.6)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])

    # Generate panels
    print("\nPanel A: Claim matrix...")
    create_claim_matrix(ax1)

    print("Panel B: Surveillance framework...")
    create_surveillance_framework(ax2)

    print("Panel C: Validation roadmap...")
    create_validation_roadmap(ax3)

    # Save
    plt.tight_layout()

    output_pdf = OUTPUT_DIR / "figure6_framework.pdf"
    output_png = OUTPUT_DIR / "figure6_framework.png"

    fig.savefig(output_pdf, dpi=300, bbox_inches='tight', pad_inches=0.3)
    fig.savefig(output_png, dpi=150, bbox_inches='tight', pad_inches=0.3)

    print(f"\n[OK] Saved: {output_pdf}")
    print(f"\n[OK] Saved: {output_png}")
    print("\n" + "="*60)
    print("Figure 6 complete!")
    print("="*60)

    plt.close()

if __name__ == "__main__":
    main()
