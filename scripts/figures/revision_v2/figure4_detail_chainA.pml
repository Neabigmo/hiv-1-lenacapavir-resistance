#!/usr/bin/env pymol
# Figure 4 Detail: Chain A with All 5 Mutations
# Shows all 5 mutation sites in the binding pocket region (chains A+B)

reinitialize

fetch 6vkv, async=0

remove solvent

hide everything

############################################
# Protein - silver-gray cartoon for chain A only
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
set stick_radius, 0.25

############################################
# All 5 Mutation sites - spheres
# N57, M66, Q67 on chain A
# N74, A105 on chain B (at binding interface)
############################################

# N57 - green (chain A)
select n57, resid 57 and chain A
show spheres, n57
color green, n57
set sphere_scale, 0.5, n57

# M66 - red (chain A)
select m66, resid 66 and chain A
show spheres, m66
color red, m66
set sphere_scale, 0.5, m66

# Q67 - blue (chain A)
select q67, resid 67 and chain A
show spheres, q67
color blue, q67
set sphere_scale, 0.5, q67

# N74 - light blue (chain B - at binding interface)
# Small transparent spheres - minimal visual presence
select n74, resid 74 and chain B
show spheres, n74
color lightblue, n74
set sphere_scale, 0.3, n74
set sphere_transparency, 0.5, n74

# A105 - cyan (chain B - at binding interface)
# Small transparent spheres - minimal visual presence
select a105, resid 105 and chain B
show spheres, a105
color cyan, a105
set sphere_scale, 0.3, a105
set sphere_transparency, 0.5, a105

############################################
# Binding pocket residues - sticks
############################################

select pocket, byres (chain A within 5 of lig)
show sticks, pocket
set stick_radius, 0.1

############################################
# Hydrogen bonds to ligand
############################################

select n57_ca, resid 57 and name CA and chain A
select lig_o, lig and name O1
distance hbond_n57, n57_ca, lig_o
set dash_color, green, hbond_n57
set dash_gap, 0.3, hbond_n57
set dash_length, 0.1, hbond_n57

select n74_ca, resid 74 and name CA and chain B
distance hbond_n74, n74_ca, lig_o
set dash_color, lightblue, hbond_n74
set dash_gap, 0.3, hbond_n74
set dash_length, 0.1, hbond_n74

select q67_cb, resid 67 and name CB and chain A
distance hbond_q67, q67_cb, lig
set dash_color, blue, hbond_q67
set dash_gap, 0.3, hbond_q67
set dash_length, 0.1, hbond_q67

############################################
# Hydrophobic contacts
############################################

select m66_cb, resid 66 and name CB and chain A
distance hydro_m66, m66_cb, lig
set dash_color, red, hydro_m66
set dash_gap, 0.5, hydro_m66
set dash_length, 0.15, hydro_m66

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

center lig
zoom lig, 7

############################################
# TO SAVE IMAGE:
# png figure4_detail_chainA.png, width=2400, height=2400, ray=1
############################################
