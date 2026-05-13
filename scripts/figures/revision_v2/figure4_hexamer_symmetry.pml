#!/usr/bin/env pymol
# Figure 4: True Hexamer with Crystal Symmetry
# Generate complete 6-fold symmetric hexamer from asymmetric unit

reinitialize

# Fetch asymmetric unit
fetch 6vkv, async=0

remove solvent

hide everything

############################################
# Generate crystal symmetry to get full hexamer
# Space group P 6 has 6-fold rotational symmetry
############################################

# Apply crystal symmetry to generate symmetry-related copies
# This creates the complete hexamer
symexp hex, 6vkv, 6vkv, 1

# Now we should have the full hexamer with 6 copies
# The naming will be: 6vkv, hex_1, hex_2, hex_3, hex_4, hex_5 (6 total)

############################################
# Show all hexamer chains
############################################

show cartoon, hex
color gray85, hex

set cartoon_fancy_sheets, 1

############################################
# LEN ligand (resn QNG) - yellow sticks
# Select all ligand copies in the hexamer
############################################

select lig, resn QNG
show sticks, lig
color yellow, lig
set stick_radius, 0.3

############################################
# All 5 Mutation sites - spheres on all chains
# Position: 57=green, 66=red, 67=blue, 74=lightblue, 105=cyan
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
color blue, pos67
set sphere_scale, 0.4, pos67

select pos74, resid 74
show spheres, pos74
color lightblue, pos74
set sphere_scale, 0.4, pos74

select pos105, resid 105
show spheres, pos105
color cyan, pos105
set sphere_scale, 0.4, pos105

############################################
# Settings
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3

# Top-down view of hexamer
orient hex
turn x, 85
zoom 2

############################################
# Count how many chains we have
# print cmd.get_model("polymer").get_chains()
############################################

############################################
# TO SAVE IMAGE:
# png figure4_hexamer_true.png, width=2400, height=2400, ray=1
############################################
