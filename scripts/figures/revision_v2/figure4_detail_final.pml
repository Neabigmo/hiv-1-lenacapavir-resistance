#!/usr/bin/env pymol
# Figure 4 Detail: Binding Site with Mutations + Ligand

reinitialize

fetch 6vkv, async=0
remove solvent
hide everything

# Chain A cartoon
show cartoon, chain A
color gray85, chain A
set cartoon_fancy_sheets, 1
set cartoon_smooth_segments, 3
set cartoon_flatness, 0.3

# Ligand - light yellow sticks
select lig, resn QNG
show sticks, lig
color yelloworange, lig
set stick_radius, 0.25
set valence, 0

# Surface on chain A
show surface, chain A
set transparency, 0.65
color gray80, chain A

# Color surface near mutations
color green, (chain A and resid 57 around 5 and elem C)
color red, (chain A and resid 66 around 5 and elem C)
color orange, (chain A and resid 67 around 5 and elem C)
color blue, (chain A and resid 74 around 5 and elem C)
color cyan, (chain A and resid 105 around 5 and elem C)
color cyan, (chain A and resid 107 around 5 and elem C)

# Mutation spheres
select n57_ca, resid 57 and name CA and chain A
show spheres, n57_ca
color green, n57_ca
set sphere_scale, 0.7, n57_ca

select m66_ca, resid 66 and name CA and chain A
show spheres, m66_ca
color red, m66_ca
set sphere_scale, 0.7, m66_ca

select q67_ca, resid 67 and name CA and chain A
show spheres, q67_ca
color orange, q67_ca
set sphere_scale, 0.7, q67_ca

select n74_ca, resid 74 and name CA and chain A
show spheres, n74_ca
color blue, n74_ca
set sphere_scale, 0.7, n74_ca

select a105_ca, resid 105 and name CA and chain A
show spheres, a105_ca
color cyan, a105_ca
set sphere_scale, 0.7, a105_ca

select t107_ca, resid 107 and name CA and chain A
show spheres, t107_ca
color cyan, t107_ca
set sphere_scale, 0.7, t107_ca

# H-bond
select n57_n, resid 57 and name N and chain A
show spheres, n57_n
color green, n57_n
set sphere_scale, 0.4, n57_n

select lig_o, (resn QNG and elem O)
distance hbond_n57_lig, n57_n, lig_o, 4.0

set dash_color, green
set dash_gap, 0.0
set dash_width, 8.0
set dash_radius, 0.1
hide label, hbond_n57_lig

# Settings
bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3
set shininess, 0
set sphere_quality, 3
set stick_quality, 3

# Zoom - detail view
zoom (chain A and resid 55-110) or lig, 5
orient (chain A and resid 50-115) or lig
turn y, 15

# TO SAVE: png figure4_detail.png, width=2400, height=2400, ray=1