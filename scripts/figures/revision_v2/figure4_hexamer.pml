#!/usr/bin/env pymol
# Figure 4: True CA Hexamer - Biological Assembly
# Use actual PDB 6VKV biological assembly (hexamer) with light surface

reinitialize

# Fetch biological assembly 1 - this should be the actual hexamer
fetch 6vkv, type=pdb1, async=0

remove solvent

hide everything

############################################
# Show all chains in hexamer with cartoon
############################################

show cartoon, all
color gray85, all

set cartoon_fancy_sheets, 1

############################################
# Light surface overlay (very transparent)
############################################

show surface, all
set transparency, 0.88

############################################
# Mutation sites as colored spheres on surface
# Surface will show through mutations
############################################

# N57 - green (hydrogen bond loss, 4890x fold change)
select pos57, resid 57
show spheres, pos57
color green, pos57
set sphere_scale, 0.5, pos57

# M66 - red (steric hindrance, 3200x fold change)
select pos66, resid 66
show spheres, pos66
color red, pos66
set sphere_scale, 0.5, pos66

# Q67 - orange (environmental sensitivity, 76x fold change)
select pos67, resid 67
show spheres, pos67
color orange, pos67
set sphere_scale, 0.5, pos67

# N74 - blue (minimal effect, 5x fold change)
select pos74, resid 74
show spheres, pos74
color blue, pos74
set sphere_scale, 0.5, pos74

# A105/T107 - cyan (compensation region)
select pos105, resid 105
show spheres, pos105
color cyan, pos105
set sphere_scale, 0.5, pos105

select pos107, resid 107
show spheres, pos107
color cyan, pos107
set sphere_scale, 0.5, pos107

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
# png figure4_hexamer.png, width=2400, height=2400, ray=1
############################################