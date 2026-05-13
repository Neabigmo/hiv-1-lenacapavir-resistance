#!/usr/bin/env pymol
# Figure 4: True CA Hexamer - Biological Assembly
# Real hexamer structure (C6 symmetry, A6 stoichiometry)

reinitialize

# Fetch the actual biological assembly - this is the true hexamer
fetch 6vkv, type=pdb1, async=0

remove solvent

hide everything

############################################
# Show all 6 chains in hexamer with cartoon
############################################

show cartoon, all
color gray85, all

set cartoon_fancy_sheets, 1
set cartoon_smooth_segments, 3
set cartoon_flatness, 0.3

############################################
# Light surface overlay (very transparent)
############################################

show surface, all
set transparency, 0.85, all

############################################
# 5 Mutation sites as colored spheres on each chain
# (These will appear on all 6 chains due to symmetry)
############################################

# N57 - green (hydrogen bond loss, 4890x fold change)
select pos57, resid 57
show spheres, pos57
color green, pos57
set sphere_scale, 0.45, pos57

# M66 - red (steric hindrance, 3200x fold change)
select pos66, resid 66
show spheres, pos66
color red, pos66
set sphere_scale, 0.45, pos66

# Q67 - orange (environmental sensitivity, 76x fold change)
select pos67, resid 67
show spheres, pos67
color orange, pos67
set sphere_scale, 0.45, pos67

# N74 - blue (minimal effect, 5x fold change)
select pos74, resid 74
show spheres, pos74
color blue, pos74
set sphere_scale, 0.45, pos74

# A105 - cyan (compensation region)
select pos105, resid 105
show spheres, pos105
color cyan, pos105
set sphere_scale, 0.45, pos105

############################################
# Top-down view of hexamer
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3
set shininess, 0

# Rotate to top-down view (C6 axis along Z)
zoom
turn x, 85
zoom

############################################
# TO SAVE IMAGE:
# png figure4_hexamer_true.png, width=2400, height=2400, ray=1
############################################