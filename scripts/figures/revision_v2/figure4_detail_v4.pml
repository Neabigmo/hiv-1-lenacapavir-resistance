#!/usr/bin/env pymol
# Figure 4: Binding Site Detail - HIV-1 Capsid LEN-CA with Mutations
# Ligand: yellow, Pocket: gray, Mutations: colored spheres, H-bond: green dashed

reinitialize
fetch 6vkv, async=0
remove solvent
hide everything

############################################
# 1. Chain A - gray cartoon
############################################
show cartoon, chain A
color gray85, chain A
set cartoon_fancy_sheets, 1
set cartoon_smooth_segments, 3
set cartoon_flatness, 0.3

############################################
# 2. Surface - semi-transparent gray
############################################
show surface, chain A
set transparency, 0.6
color gray85, chain A

############################################
# 3. Mutation spheres - 6 sites
############################################
# N57 - green
select n57_ca, resid 57 and name CA and chain A
show spheres, n57_ca
color green, n57_ca
set sphere_scale, 0.65

# M66 - red
select m66_ca, resid 66 and name CA and chain A
show spheres, m66_ca
color red, m66_ca
set sphere_scale, 0.65

# Q67 - orange
select q67_ca, resid 67 and name CA and chain A
show spheres, q67_ca
color orange, q67_ca
set sphere_scale, 0.65

# N74 - blue
select n74_ca, resid 74 and name CA and chain A
show spheres, n74_ca
color blue, n74_ca
set sphere_scale, 0.65

# A105 - cyan
select a105_ca, resid 105 and name CA and chain A
show spheres, a105_ca
color cyan, a105_ca
set sphere_scale, 0.65

# T107 - lightblue (浅青色)
select t107_ca, resid 107 and name CA and chain A
show spheres, t107_ca
color lightblue, t107_ca
set sphere_scale, 0.65

############################################
# 4. Surface coloring near mutations
############################################
select surf57, (chain A and resid 57 around 5 and elem C)
color green, surf57

select surf66, (chain A and resid 66 around 5 and elem C)
color red, surf66

select surf67, (chain A and resid 67 around 5 and elem C)
color orange, surf67

select surf74, (chain A and resid 74 around 5 and elem C)
color blue, surf74

select surf105, (chain A and resid 105 around 5 and elem C)
color cyan, surf105

select surf107, (chain A and resid 107 around 5 and elem C)
color lightblue, surf107

############################################
# 5. BINDING POCKET - chain A only, GRAY
# CRITICAL: Only select chain A residues, NOT ligand atoms
############################################
# Select chain A residues within 5A of the ligand
select pocket_residues, (chain A within 5 of resn QNG) and chain A
show sticks, pocket_residues
color gray85, pocket_residues
set stick_radius, 0.12

############################################
# 6. LIGAND - YELLOW, only near chain A
# CRITICAL: Select ligand atoms ONLY (resn QNG), color yellow
############################################
# Select ligand atoms within 8A of chain A residue 57
select ligand_atoms, (resn QNG around 6 and chain A around 8) or (resn QNG and chain A around 8)
show sticks, ligand_atoms
color yellow, ligand_atoms
set stick_radius, 0.25
set valence, 0

# Also show as spheres for visibility
show spheres, ligand_atoms
color yellow, ligand_atoms
set sphere_scale, 0.3

############################################
# 7. H-bond: N57 backbone N to ligand
############################################
select n57_n, resid 57 and name N and chain A
show spheres, n57_n
color green, n57_n
set sphere_scale, 0.35

# Select ligand oxygen atoms
select lig_o, resn QNG and elem O

# Create H-bond dashed line
distance hb57, n57_n, lig_o, 3.5
set dash_color, green
set dash_gap, 0.0
set dash_width, 12.0
set dash_radius, 0.08
hide label, hb57

############################################
# 8. Render settings
############################################
bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3
set shininess, 0
set sphere_quality, 2
set stick_quality, 2

############################################
# 9. Zoom and orient - detail view
############################################
zoom (chain A and resid 50-115) or resn QNG, 6
orient (chain A and resid 50-115) or resn QNG

############################################
# TO SAVE IMAGE:
# png figure4_detail.png, width=2400, height=2400, ray=1
############################################
