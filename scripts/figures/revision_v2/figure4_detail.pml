#!/usr/bin/env pymol
# Figure 4 Detail: Binding Site with Mutations + Ligand
# Chain A + ligand on chain D + surface with mutation colors

reinitialize

fetch 6vkv, async=0

remove solvent

hide everything

############################################
# Chain A - protein cartoon
############################################

show cartoon, chain A
color gray85, chain A

set cartoon_fancy_sheets, 1
set cartoon_smooth_segments, 3
set cartoon_flatness, 0.3

############################################
# Ligand on chain D - yellow sticks
# MUST show chain D explicitly after hide everything
############################################

select ligA, chain D
show sticks, ligA
color yellow, ligA
set stick_radius, 0.25

############################################
# Surface - silver-white with transparency
# Must show after hide everything
############################################

show surface, chain A
set transparency, 0.7
color gray70, chain A

# Color surface atoms near mutations with mutation colors
color green, (chain A and resid 57 around 4 and elem C)
color red, (chain A and resid 66 around 4 and elem C)
color orange, (chain A and resid 67 around 4 and elem C)
color blue, (chain A and resid 74 around 4 and elem C)
color cyan, (chain A and resid 105 around 4 and elem C)
color cyan, (chain A and resid 107 around 4 and elem C)

# Hide other chains completely (only keep A and D)
hide cartoon, (not chain A)
hide spheres, (not chain A)

############################################
# 5 Mutation sites - colored spheres
# Including T107 for cyan (compensation region)
############################################

# N57 - green
select n57_ca, resid 57 and name CA and chain A
show spheres, n57_ca
color green, n57_ca
set sphere_scale, 0.7

# M66 - red
select m66_ca, resid 66 and name CA and chain A
show spheres, m66_ca
color red, m66_ca
set sphere_scale, 0.7

# Q67 - orange
select q67_ca, resid 67 and name CA and chain A
show spheres, q67_ca
color orange, q67_ca
set sphere_scale, 0.7

# N74 - blue
select n74_ca, resid 74 and name CA and chain A
show spheres, n74_ca
color blue, n74_ca
set sphere_scale, 0.7

# A105 - cyan
select a105_ca, resid 105 and name CA and chain A
show spheres, a105_ca
color cyan, a105_ca
set sphere_scale, 0.7

# T107 - cyan (compensation region, second point)
select t107_ca, resid 107 and name CA and chain A
show spheres, t107_ca
color cyan, t107_ca
set sphere_scale, 0.7

############################################
# Hydrogen bond: N57 backbone N to ligand O29
# PROMINENT green dashed line
############################################

# Select N57 backbone N
select n57_n, resid 57 and name N and chain A

# Select ligand oxygens near N57
select lig_o, (chain D around 3.5 and elem O)

# Create distance measurement for H-bond visualization
distance hbond_n57_lig, n57_n, lig_o, 3.5

# Style the hydrogen bond dashed line
set dash_color, green
set dash_gap, 0.15
set dash_width, 5.0
set dash_radius, 0.12
hide label, hbond_n57_lig

############################################
# Zoom on binding site interface
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3
set shininess, 0

# Zoom to show mutation-ligand interface
zoom (chain A and resid 55-110) or chain D, 12

############################################
# TO SAVE IMAGE:
# png figure4_detail.png, width=2400, height=2400, ray=1
############################################