#!/usr/bin/env pymol
"""
Figure 4: LEN-CA Hexamer Structure (PDB 6VKV)
Correct biological assembly rendering.

Key difference from asymmetric unit:
  - fetch 6vkv        -> 3 chains (asymmetric unit)
  - fetch 6vkv, type=pdb1 -> 6 chains (true hexamer)
"""

import os
import sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUTPUT_DIR = os.path.join(BASE_DIR, "manuscript", "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_FIGURE = os.path.join(OUTPUT_DIR, "figure4a_pymol_depth.png")


def create_pymol_figure():
    """Create publication-quality LEN-CA hexamer structure figure."""
    from pymol import cmd

    print("="*60)
    print("Generating Figure 4: LEN-CA Hexamer Structure")
    print("PDB: 6VKV (Biological Assembly)")
    print("="*60)

    # Load biological assembly directly
    print("\n[Step 1] Fetching biological assembly (type=pdb1)...")
    cmd.fetch("6vkv", type="pdb1", async_=0)

    # Remove solvent
    cmd.remove("solvent")

    # Setup environment
    print("[Step 2] Setting up rendering environment...")
    cmd.bg_color("white")
    cmd.hide("everything", "all")
    cmd.set("ray_opaque_background", "off")

    # Show cartoon representation
    print("[Step 3] Rendering protein cartoon with chain colors...")
    cmd.show("cartoon")

    # Apply chainbow coloring for distinct subunit colors
    cmd.util.chainbow("polymer")

    # LEN ligand (resn QNG)
    print("[Step 4] Rendering LEN ligand (QNG)...")
    cmd.select("lig", "resn QNG")
    cmd.show("sticks", "lig")
    cmd.color("red", "lig")
    cmd.set("stick_radius", 0.22, "lig")

    # Mutation sites as spheres
    print("[Step 5] Highlighting mutation sites...")

    mutations = [
        ("n57", "resid 57", "green"),
        ("m66", "resid 66", "red"),
        ("q67", "resid 67", "orange"),
        ("n74", "resid 74", "blue"),
        ("a105", "resid 105", "cyan"),
        ("t107", "resid 107", "cyan"),
    ]

    for name, selection, color in mutations:
        cmd.select(name, selection)
        cmd.show("spheres", name)
        cmd.color(color, name)
        cmd.set("sphere_scale", 0.4, name)
        print(f"  Added {name} - {color}")

    # Transparency for depth
    cmd.set("cartoon_transparency", 0.15)

    # Rendering quality
    print("[Step 6] Setting up high-quality rendering...")
    cmd.set("antialias", 2)
    cmd.set("specular", 0.3)
    cmd.set("orthoscopic", "on")
    cmd.set("depth_cue", "off")

    # Camera setup
    print("[Step 7] Setting view angle...")
    cmd.orient()
    cmd.zoom(5)

    # Ray tracing and export
    print(f"\n[Export] Rendering at 2400x2400 pixels...")
    cmd.png(OUTPUT_FIGURE, width=2400, height=2400, dpi=300, ray=1)

    # Cleanup
    cmd.delete("all")

    print(f"\n[OK] Figure saved: {OUTPUT_FIGURE}")
    return True


def main():
    try:
        success = create_pymol_figure()
        if success:
            print("\n" + "="*60)
            print("Figure 4 complete!")
            print("="*60)
    except Exception as e:
        print(f"\n[ERROR] {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()