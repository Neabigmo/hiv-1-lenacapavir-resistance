#!/usr/bin/env pymol
# Figure 4: True CA Hexamer - Biological Assembly
# Use actual PDB 6VKV biological assembly (pdb1 = hexamer)

reinitialize

# Fetch biological assembly 1 - actual 6-subunit hexamer
fetch 6vkv, type=pdb1, async=0

remove solvent

hide everything

############################################
# Show all 6 chains with cartoon
############################################

show cartoon, all
color gray85, all

set cartoon_fancy_sheets, 1

############################################
# Light surface overlay
############################################

show surface, all
set transparency, 0.85, all

############################################
# Mutation sites - 5 spheres per chain
# Position: 57=green, 66=red, 67=orange, 74=blue, 105=cyan
############################################

select pos57, resid 57
show spheres, pos57
color green, pos57
set sphere_scale, 0.4, pos57

select pos66, resid 66
show spheres, pos66
color red, pos66
set sphere_scale, 0.4, pos66

select pos67, resid 67
show spheres, pos67
color orange, pos67
set sphere_scale, 0.4, pos67

select pos74, resid 74
show spheres, pos74
color blue, pos74
set sphere_scale, 0.4, pos74

select pos105, resid 105
show spheres, pos105
color cyan, pos105
set sphere_scale, 0.4, pos105

############################################
# Settings - top down view
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3

zoom
turn x, 85
zoom

############################################
# TO SAVE IMAGE:
# png figure4_hexamer_true.png, width=2400, height=2400, ray=1
############################################