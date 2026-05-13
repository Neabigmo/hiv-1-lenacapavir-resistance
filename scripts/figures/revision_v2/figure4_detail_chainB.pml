#!/usr/bin/env pymol
# Figure 4 Detail: Chain B Binding Site
# Detailed view of LEN binding pocket - Chain B

reinitialize

# Fetch asymmetric unit
fetch 6vkv, async=0

remove solvent

hide everything

############################################
# Protein - silver-gray cartoon for chain B
############################################

show cartoon, chain B
color gray85, chain B

# Cartoon settings
set cartoon_fancy_sheets, 1

############################################
# LEN ligand (resn QNG) - yellow sticks
# Select from both chains to show at interface
############################################

select lig, resn QNG and chain A+B
show sticks, lig
color yellow, lig
set stick_radius, 0.25

############################################
# Mutation sites - spheres (no labels)
############################################

# N74 - light blue
select n74, resid 74 and chain B
show spheres, n74
color lightblue, n74
set sphere_scale, 0.5, n74

# A105 - cyan
select a105, resid 105 and chain B
show spheres, a105
color cyan, a105
set sphere_scale, 0.5, a105

############################################
# Binding pocket residues - sticks
############################################

select pocket, byres (chain B within 5 of lig)
show sticks, pocket
set stick_radius, 0.1

############################################
# Hydrogen bonds to ligand
############################################

# N74 to ligand
select n74_ca, resid 74 and name CA and chain B
select lig_o, lig and name O1
distance hbond_n74, n74_ca, lig_o
set dash_color, lightblue, hbond_n74
set dash_gap, 0.3, hbond_n74
set dash_length, 0.1, hbond_n74

############################################
# Hydrophobic contacts
############################################

# A105 hydrophobic interaction with ligand
select a105_cb, resid 105 and name CB and chain B
distance hydro_a105, a105_cb, lig
set dash_color, cyan, hydro_a105
set dash_gap, 0.5, hydro_a105
set dash_length, 0.15, hydro_a105

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
# png figure4_detail_chainB.png, width=2400, height=2400, ray=1
############################################
