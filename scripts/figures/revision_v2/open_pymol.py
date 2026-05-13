#!/usr/bin/env python3
"""
Launch PyMOL with pre-configured LEN-CA hexamer structure
Run this script to open PyMOL GUI with the structure loaded
"""

from pymol import cmd
import sys

# Fetch the structure
print("Fetching 6VKV structure...")
cmd.fetch("6VKV", async_=0)

print("Setting up visualization...")
cmd.bg_color("white")

# Hide everything, show cartoon
cmd.hide("everything", "all")
cmd.show("cartoon", "all")
cmd.color("white", "all")

# Color chains distinctly
chains = ['A', 'B', 'C', 'D', 'E', 'F']
chain_colors = ['blue', 'green', 'yellow', 'orange', 'purple', 'pink']
for i, chain in enumerate(chains):
    cmd.color(chain_colors[i], f"chain {chain}")

# Show surface
cmd.show("surface", "all")
cmd.set("surface_color", "white")
cmd.set("transparency", 0.2)

# LEN ligand
cmd.hide("everything", "organic")
cmd.show("sticks", "organic")
cmd.color("red", "organic")
cmd.set("stick_radius", 0.15, "organic")

# Highlight key resistance residues
cmd.select("res_n57", "resn ASP and resid 57")
cmd.show("spheres", "res_n57")
cmd.color("green", "res_n57")
cmd.set("sphere_scale", 0.8, "res_n57")

cmd.select("res_m66", "resn MET and resid 66")
cmd.show("spheres", "res_m66")
cmd.color("firebrick", "res_m66")
cmd.set("sphere_scale", 0.9, "res_m66")

cmd.select("res_q67", "resn GLN and resid 67")
cmd.show("spheres", "res_q67")
cmd.color("orange", "res_q67")
cmd.set("sphere_scale", 0.8, "res_q67")

cmd.select("res_n74", "resn ASP and resid 74")
cmd.show("spheres", "res_n74")
cmd.color("blue", "res_n74")
cmd.set("sphere_scale", 0.8, "res_n74")

cmd.select("res_a105", "resn ALA and resid 105")
cmd.show("spheres", "res_a105")
cmd.color("cyan", "res_a105")
cmd.set("sphere_scale", 0.7, "res_a105")

cmd.select("res_t107", "resn THR and resid 107")
cmd.show("spheres", "res_t107")
cmd.color("cyan", "res_t107")
cmd.set("sphere_scale", 0.7, "res_t107")

# Camera setup
cmd.zoom("all", 100, 0)
cmd.center("organic")
cmd.rotate("z", 20)
cmd.rotate("y", 60)
cmd.rotate("x", 15)

# Set labels
cmd.set("label_font_id", 18)
cmd.set("label_size", 20)
cmd.set("label_color", "black")
cmd.set("label_bg_color", "white")
cmd.set("label_bg_alpha", 0.9)

# Add labels at residue positions
cmd.label("resn ASP and resid 57", '"N57"')
cmd.label("resn MET and resid 66", '"M66"')
cmd.label("resn GLN and resid 67", '"Q67"')
cmd.label("resn ASP and resid 74", '"N74"')
cmd.label("resn ALA and resid 105", '"A105"')
cmd.label("resn THR and resid 107", '"T107"')
cmd.label("resn LEN", '"LEN"')

# Save session file
session_path = "H:/2026try/4.20/JMV/manuscript/figures/lenca_hexamer.pse"
cmd.save(session_path)
print(f"\nSession saved: {session_path}")
print("You can load this file later to continue from this view.")

print("\n" + "="*60)
print("PyMOL session ready!")
print("Adjust the view as you like, then:")
print("1. File > Save Session As... to save your adjustments")
print("2. Use File > Export Image for high-res PNG")
print("="*60)

# Launch GUI
cmd.extend("show_gui", __import__('pymol').launch)