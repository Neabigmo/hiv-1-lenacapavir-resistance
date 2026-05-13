#!/usr/bin/env python3
"""
Graphical Abstract: Study Logic Flow
Single comprehensive diagram showing the complete study workflow
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

# Setup
BASE_DIR = Path(__file__).parent.parent.parent
OUTPUT_DIR = BASE_DIR / "manuscript" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

plt.style.use('seaborn-v0_8-paper')

def create_graphical_abstract():
    """Create comprehensive graphical abstract"""

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(7, 9.5, 'Mutation-First Resistance Architecture of HIV-1 Lenacapavir',
           ha='center', fontsize=14, weight='bold')
    ax.text(7, 9, 'From Compiled Evidence to Graded Claims',
           ha='center', fontsize=11, style='italic', color='gray')

    # ===== SECTION 1: Question (Top) =====
    question_box = FancyBboxPatch((0.5, 7.8), 13, 0.8,
                                   boxstyle="round,pad=0.1",
                                   facecolor='#FFE6E6', edgecolor='red', linewidth=2)
    ax.add_patch(question_box)
    ax.text(7, 8.4, 'Research Question', ha='center', fontsize=10, weight='bold')
    ax.text(7, 8, 'Does HIV-1 subtype modulate lenacapavir resistance, or is mutation identity dominant?',
           ha='center', fontsize=9)

    # Arrow down
    ax.arrow(7, 7.7, 0, -0.4, head_width=0.3, head_length=0.15, fc='black', ec='black', linewidth=2)

    # ===== SECTION 2: Evidence Collection =====
    evidence_box = FancyBboxPatch((0.5, 6.2), 4, 1,
                                   boxstyle="round,pad=0.08",
                                   facecolor='#E8F4F8', edgecolor='blue', linewidth=2)
    ax.add_patch(evidence_box)
    ax.text(2.5, 6.9, 'Evidence Collection', ha='center', fontsize=9, weight='bold')
    ax.text(2.5, 6.6, 'PRISMA synthesis', ha='center', fontsize=7)
    ax.text(2.5, 6.4, '87 identified → 11 included', ha='center', fontsize=7)

    # Arrow down
    ax.arrow(2.5, 6.1, 0, -0.3, head_width=0.2, head_length=0.1, fc='blue', ec='blue', linewidth=1.5)

    # ===== SECTION 3: Harmonization =====
    harmonize_box = FancyBboxPatch((0.5, 4.8), 4, 0.9,
                                    boxstyle="round,pad=0.08",
                                    facecolor='#E8F4F8', edgecolor='blue', linewidth=2)
    ax.add_patch(harmonize_box)
    ax.text(2.5, 5.4, 'Data Harmonization', ha='center', fontsize=9, weight='bold')
    ax.text(2.5, 5.1, '24 observations, 23 complete', ha='center', fontsize=7)
    ax.text(2.5, 4.9, 'Context tiers + quality scores', ha='center', fontsize=7)

    # Arrow to three analysis branches
    ax.arrow(2.5, 4.7, 0, -0.3, head_width=0.2, head_length=0.1, fc='black', ec='black', linewidth=2)
    ax.arrow(2.5, 4.3, 2, -0.5, head_width=0.2, head_length=0.1, fc='black', ec='black', linewidth=1.5)
    ax.arrow(2.5, 4.3, 7, -0.5, head_width=0.2, head_length=0.1, fc='black', ec='black', linewidth=1.5)

    # ===== SECTION 4: Three Analysis Lines =====

    # Line 1: Phenotypic Analysis
    pheno_box = FancyBboxPatch((0.3, 2.5), 3.5, 1.2,
                                boxstyle="round,pad=0.08",
                                facecolor='#E8F5E9', edgecolor='green', linewidth=2)
    ax.add_patch(pheno_box)
    ax.text(2, 3.5, 'Phenotypic Analysis', ha='center', fontsize=9, weight='bold', color='green')
    ax.text(2, 3.2, '• Model comparison (M0-M3)', ha='left', fontsize=7)
    ax.text(2, 3, '• Bootstrap ranking (n=1000)', ha='left', fontsize=7)
    ax.text(2, 2.8, '• Epistasis & context effects', ha='left', fontsize=7)
    ax.text(2, 2.6, 'Result: Mutation-only model best', ha='center', fontsize=7, weight='bold', style='italic')

    # Line 2: Structural Analysis
    struct_box = FancyBboxPatch((4.5, 2.5), 3.5, 1.2,
                                 boxstyle="round,pad=0.08",
                                 facecolor='#F3E5F5', edgecolor='purple', linewidth=2)
    ax.add_patch(struct_box)
    ax.text(6.2, 3.5, 'Structural Analysis', ha='center', fontsize=9, weight='bold', color='purple')
    ax.text(6.2, 3.2, '• ΔΔG calculations', ha='left', fontsize=7)
    ax.text(6.2, 3, '• H-bond & steric metrics', ha='left', fontsize=7)
    ax.text(6.2, 2.8, '• Perturbation scores', ha='left', fontsize=7)
    ax.text(6.2, 2.6, 'Result: Mechanism validation', ha='center', fontsize=7, weight='bold', style='italic')

    # Line 3: Evolutionary Analysis
    evol_box = FancyBboxPatch((8.7, 2.5), 3.5, 1.2,
                               boxstyle="round,pad=0.08",
                               facecolor='#FFF3E0', edgecolor='orange', linewidth=2)
    ax.add_patch(evol_box)
    ax.text(10.4, 3.5, 'Evolutionary Analysis', ha='center', fontsize=9, weight='bold', color='orange')
    ax.text(10.4, 3.2, '• Conservation scores', ha='left', fontsize=7)
    ax.text(10.4, 3, '• Subtype frequencies', ha='left', fontsize=7)
    ax.text(10.4, 2.8, '• Natural polymorphisms', ha='left', fontsize=7)
    ax.text(10.4, 2.6, 'Result: High conservation', ha='center', fontsize=7, weight='bold', style='italic')

    # Arrows converging to graded claims
    ax.arrow(2, 2.4, 3, -0.6, head_width=0.2, head_length=0.1, fc='black', ec='black', linewidth=1.5)
    ax.arrow(6.2, 2.4, 0.8, -0.6, head_width=0.2, head_length=0.1, fc='black', ec='black', linewidth=1.5)
    ax.arrow(10.4, 2.4, -3.4, -0.6, head_width=0.2, head_length=0.1, fc='black', ec='black', linewidth=1.5)

    # ===== SECTION 5: Graded Claims =====
    claims_box = FancyBboxPatch((4.5, 0.8), 5, 0.9,
                                 boxstyle="round,pad=0.08",
                                 facecolor='#FFFDE7', edgecolor='#F57F17', linewidth=3)
    ax.add_patch(claims_box)
    ax.text(7, 1.5, 'Evidence-Graded Claims', ha='center', fontsize=10, weight='bold')
    ax.text(5, 1.2, 'Strong:', ha='left', fontsize=8, weight='bold', color='green')
    ax.text(5.8, 1.2, 'Mutation ranking stable', ha='left', fontsize=7)
    ax.text(5, 1, 'Moderate:', ha='left', fontsize=8, weight='bold', color='orange')
    ax.text(5.8, 1, 'Subtype adds limited value', ha='left', fontsize=7)
    ax.text(5, 0.85, 'Hypothesis:', ha='left', fontsize=8, weight='bold', color='red')
    ax.text(5.8, 0.85, 'Compensatory evolution', ha='left', fontsize=7)

    # Arrow to output
    ax.arrow(7, 0.7, 0, -0.3, head_width=0.3, head_length=0.1, fc='black', ec='black', linewidth=2)

    # ===== SECTION 6: Output =====
    output_box = FancyBboxPatch((3.5, 0.05), 7, 0.3,
                                 boxstyle="round,pad=0.05",
                                 facecolor='#C8E6C9', edgecolor='green', linewidth=2)
    ax.add_patch(output_box)
    ax.text(7, 0.2, 'Mutation-First Surveillance Framework for Clinical Implementation',
           ha='center', fontsize=9, weight='bold')

    # ===== Side annotations =====
    # Left side: Data flow
    ax.text(0.2, 6, 'Data\nFlow', ha='center', va='center', fontsize=8, rotation=90,
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    # Right side: Evidence strength
    ax.text(13.5, 4, 'Evidence\nStrength\nIncreases', ha='center', va='center', fontsize=8, rotation=-90,
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

    # Add key numbers
    ax.text(12.5, 7, 'Key Numbers:', ha='left', fontsize=8, weight='bold',
           bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', linewidth=1))
    ax.text(12.5, 6.7, '• 11 sources', ha='left', fontsize=7)
    ax.text(12.5, 6.5, '• 24 observations', ha='left', fontsize=7)
    ax.text(12.5, 6.3, '• 16 mutations', ha='left', fontsize=7)
    ax.text(12.5, 6.1, '• 3 analysis lines', ha='left', fontsize=7)
    ax.text(12.5, 5.9, '• 12 graded claims', ha='left', fontsize=7)

    return fig

def main():
    """Generate graphical abstract"""
    print("="*60)
    print("Generating Graphical Abstract")
    print("="*60)

    fig = create_graphical_abstract()

    # Save
    output_pdf = OUTPUT_DIR / "graphical_abstract.pdf"
    output_png = OUTPUT_DIR / "graphical_abstract.png"

    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')

    print(f"\n[OK] Saved: {output_pdf}")
    print(f"\n[OK] Saved: {output_png}")
    print("\n" + "="*60)
    print("Graphical Abstract complete!")
    print("="*60)

    plt.close()

if __name__ == "__main__":
    main()
