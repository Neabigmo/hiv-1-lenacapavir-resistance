#!/usr/bin/env pymol
# Figure 4 Detail: Single Chain - 5 Mutations
# Clean, simplified view with one sphere per mutation site
# NO individual hydrogen bond lines - just show proximity

reinitialize

fetch 6vkv, async=0

remove solvent

hide everything

############################################
# Protein - silver-gray cartoon
############################################

show cartoon, chain A
color gray85, chain A

set cartoon_fancy_sheets, 1

############################################
# LEN ligand (resn QNG) - yellow sticks
############################################

select lig, resn QNG and chain A
show sticks, lig
color yellow, lig
set stick_radius, 0.3

############################################
# 5 Mutation sites - ONE sphere each
# No labels, no multiple selections
############################################

# N57 - green
select n57, resid 57 and name CA and chain A
show spheres, n57
color green, n57
set sphere_scale, 0.6, n57

# M66 - red
select m66, resid 66 and name CA and chain A
show spheres, m66
color red, m66
set sphere_scale, 0.6, m66

# Q67 - blue
select q67, resid 67 and name CA and chain A
show spheres, q67
color blue, q67
set sphere_scale, 0.6, q67

# N74 - light blue
select n74, resid 74 and name CA and chain A
show spheres, n74
color lightblue, n74
set sphere_scale, 0.6, n74

# A105 - cyan
select a105, resid 105 and name CA and chain A
show spheres, a105
color cyan, a105
set sphere_scale, 0.6, a105

############################################
# NO hydrogen bond lines - just show spheres
# The spheres at correct positions indicate interactions
############################################

############################################
# Settings
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3

# Zoom into binding site
center lig
zoom lig, 8

############################################
# TO SAVE IMAGE:
# png figure4_single_chain_simple.png, width=2400, height=2400, ray=1
############################################