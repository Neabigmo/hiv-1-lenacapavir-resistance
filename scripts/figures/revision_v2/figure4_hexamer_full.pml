#!/usr/bin/env pymol
# Figure 4: Complete Hexamer with All 6 Subunits
# Shows all chains A-F with real mutation distribution and LEN drug

reinitialize

# Fetch biological assembly (hexamer - 6 chains)
fetch 6vkv, type=pdb1, async=0

remove solvent

hide everything

############################################
# Protein cartoon - all 6 chains with different colors
############################################

show cartoon, chain A+B+C+D+E+F

# Each chain a different color for clarity
color gray70, chain A
color gray80, chain B
color gray70, chain C
color gray80, chain D
color gray70, chain E
color gray80, chain F

set cartoon_fancy_sheets, 1

############################################
# LEN ligand (resn QNG) - yellow sticks
# All ligand copies in the hexamer
############################################

select lig, resn QNG
show sticks, lig
color yellow, lig
set stick_radius, 0.3

############################################
# Mutation sites on all chains
# Each chain has mutations at positions: 57, 66, 67, 74, 105
# Shown as colored spheres
############################################

# Position 57 - green (N57)
select pos57, resid 57
show spheres, pos57
color green, pos57
set sphere_scale, 0.4, pos57

# Position 66 - red (M66)
select pos66, resid 66
show spheres, pos66
color red, pos66
set sphere_scale, 0.4, pos66

# Position 67 - blue (Q67)
select pos67, resid 67
show spheres, pos67
color blue, pos67
set sphere_scale, 0.4, pos67

# Position 74 - light blue (N74)
select pos74, resid 74
show spheres, pos74
color lightblue, pos74
set sphere_scale, 0.4, pos74

# Position 105 - cyan (A105)
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
orient chain A+B+C+D+E+F
turn x, 85
zoom 2

############################################
# TO SAVE IMAGE:
# png figure4_hexamer_full.png, width=2400, height=2400, ray=1
############################################
