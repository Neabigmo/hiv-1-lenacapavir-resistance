#!/usr/bin/env pymol
# Figure 4 Detail: Single Chain - 5 Mutations
# Only show ONE chain (chain A) with 5 mutation sites
# NO ligand, clean view showing mutation positions only

reinitialize

# Fetch asymmetric unit to get chain A
fetch 6vkv, async=0

remove solvent

hide everything

############################################
# Show only chain A - gray cartoon
############################################

show cartoon, chain A
color gray85, chain A

set cartoon_fancy_sheets, 1
set cartoon_smooth_segments, 3
set cartoon_flatness, 0.3

############################################
# 5 Mutation sites as colored spheres
# Only on chain A - clean single chain view
############################################

# N57 - green (hydrogen bond loss, 4890x)
select n57_ca, resid 57 and name CA and chain A
show spheres, n57_ca
color green, n57_ca
set sphere_scale, 0.7, n57_ca

# M66 - red (steric hindrance, 3200x)
select m66_ca, resid 66 and name CA and chain A
show spheres, m66_ca
color red, m66_ca
set sphere_scale, 0.7, m66_ca

# Q67 - orange (environmental sensitivity, 76x)
select q67_ca, resid 67 and name CA and chain A
show spheres, q67_ca
color orange, q67_ca
set sphere_scale, 0.7, q67_ca

# N74 - blue (minimal effect, 5x)
select n74_ca, resid 74 and name CA and chain A
show spheres, n74_ca
color blue, n74_ca
set sphere_scale, 0.7, n74_ca

# A105 - cyan (compensation region, first point)
select a105_ca, resid 105 and name CA and chain A
show spheres, a105_ca
color cyan, a105_ca
set sphere_scale, 0.7, a105_ca

# T107 - cyan (compensation region, second point)
select t107_ca, resid 107 and name CA and chain A
show spheres, t107_ca
color cyan, t107_ca
set sphere_scale, 0.7, t107_ca

############################################
# Settings - zoom into mutation region
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3
set shininess, 0

# Center and zoom on the mutation cluster
center (chain A and resid 57-110)
zoom (chain A and resid 57-110), 8

############################################
# TO SAVE IMAGE:
# png figure4_single_chain_5mut.png, width=2400, height=2400, ray=1
############################################