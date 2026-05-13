#!/usr/bin/env pymol
# Figure 4 Detail: Single Chain - 5 Mutations with Key Interactions
# Simplified: 1 sphere per mutation, only 2 key hydrogen bonds

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
# 5 Mutation sites - ONE sphere each (larger)
############################################

# N57 - green
select n57_ca, resid 57 and name CA and chain A
show spheres, n57_ca
color green, n57_ca
set sphere_scale, 0.7, n57_ca

# M66 - red
select m66_ca, resid 66 and name CA and chain A
show spheres, m66_ca
color red, m66_ca
set sphere_scale, 0.7, m66_ca

# Q67 - blue
select q67_ca, resid 67 and name CA and chain A
show spheres, q67_ca
color blue, q67_ca
set sphere_scale, 0.7, q67_ca

# N74 - light blue
select n74_ca, resid 74 and name CA and chain A
show spheres, n74_ca
color lightblue, n74_ca
set sphere_scale, 0.7, n74_ca

# A105 - cyan
select a105_ca, resid 105 and name CA and chain A
show spheres, a105_ca
color cyan, a105_ca
set sphere_scale, 0.7, a105_ca

############################################
# ONLY 2 key hydrogen bonds - simplified
############################################

select lig_o1, lig and name O1
select lig_o2, lig and name O2

# N57 to ligand - main H-bond
distance hbond1, n57_ca, lig_o1
set dash_color, green, hbond1
set dash_gap, 0.4, hbond1
set dash_width, 2.0, hbond1

# N74 to ligand - main H-bond
distance hbond2, n74_ca, lig_o2
set dash_color, lightblue, hbond2
set dash_gap, 0.4, hbond2
set dash_width, 2.0, hbond2

############################################
# Settings
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3

center lig
zoom lig, 8

############################################
# TO SAVE IMAGE:
# png figure4_single_chain_interactions.png, width=2400, height=2400, ray=1
############################################