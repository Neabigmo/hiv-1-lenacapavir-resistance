#!/usr/bin/env pymol
# Figure 4 Detail: Chains A and B - All 5 mutation sites each
# Shows both chains with complete mutation sites visible

reinitialize

# Fetch asymmetric unit
fetch 6vkv, async=0

remove solvent

hide everything

############################################
# Protein - silver-gray cartoon for chains A and B
############################################

show cartoon, chain A+B
color gray85, chain A+B

# Cartoon settings
set cartoon_fancy_sheets, 1

############################################
# LEN ligand (resn QNG) - yellow sticks
############################################

select lig, resn QNG and chain A+B
show sticks, lig
color yellow, lig
set stick_radius, 0.25

############################################
# Mutation sites - spheres on chains A and B
# All 5 sites shown on EACH chain where they occur
############################################

# Chain A mutation sites
# N57 - green
select n57, resid 57 and chain A
show spheres, n57
color green, n57
set sphere_scale, 0.5, n57

# M66 - red
select m66, resid 66 and chain A
show spheres, m66
color red, m66
set sphere_scale, 0.5, m66

# Q67 - blue
select q67, resid 67 and chain A
show spheres, q67
color blue, q67
set sphere_scale, 0.5, q67

# N74 - light blue
select n74, resid 74 and chain A
show spheres, n74
color lightblue, n74
set sphere_scale, 0.5, n74

# A105 - cyan
select a105, resid 105 and chain A
show spheres, a105
color cyan, a105
set sphere_scale, 0.5, a105

# Chain B mutation sites
# N57 - green
select n57_b, resid 57 and chain B
show spheres, n57_b
color green, n57_b
set sphere_scale, 0.5, n57_b

# M66 - red
select m66_b, resid 66 and chain B
show spheres, m66_b
color red, m66_b
set sphere_scale, 0.5, m66_b

# Q67 - blue
select q67_b, resid 67 and chain B
show spheres, q67_b
color blue, q67_b
set sphere_scale, 0.5, q67_b

# N74 - light blue
select n74_b, resid 74 and chain B
show spheres, n74_b
color lightblue, n74_b
set sphere_scale, 0.5, n74_b

# A105 - cyan
select a105_b, resid 105 and chain B
show spheres, a105_b
color cyan, a105_b
set sphere_scale, 0.5, a105_b

############################################
# Binding pocket residues - sticks
############################################

select pocket, byres (chain A+B within 5 of lig)
show sticks, pocket
set stick_radius, 0.1

############################################
# Hydrogen bonds to ligand
############################################

# N57 to ligand (chain A)
select n57_ca, resid 57 and name CA and chain A
select lig_o, lig and name O1
distance hbond_n57, n57_ca, lig_o
set dash_color, green, hbond_n57
set dash_gap, 0.3, hbond_n57
set dash_length, 0.1, hbond_n57

# N74 to ligand (chain B)
select n74_ca, resid 74 and name CA and chain B
distance hbond_n74, n74_ca, lig_o
set dash_color, lightblue, hbond_n74
set dash_gap, 0.3, hbond_n74
set dash_length, 0.1, hbond_n74

# Q67 to ligand (chain A)
select q67_cb, resid 67 and name CB and chain A
distance hbond_q67, q67_cb, lig
set dash_color, blue, hbond_q67
set dash_gap, 0.3, hbond_q67
set dash_length, 0.1, hbond_q67

############################################
# Hydrophobic contacts
############################################

# M66 hydrophobic interaction with ligand (chain A)
select m66_cb, resid 66 and name CB and chain A
distance hydro_m66, m66_cb, lig
set dash_color, red, hydro_m66
set dash_gap, 0.5, hydro_m66
set dash_length, 0.15, hydro_m66

# A105 hydrophobic interaction (chain B)
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
zoom lig, 6

############################################
# TO SAVE IMAGE:
# png figure4_detail_AB.png, width=2400, height=2400, ray=1
############################################
