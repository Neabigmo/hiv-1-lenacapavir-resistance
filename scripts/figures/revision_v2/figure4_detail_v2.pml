#!/usr/bin/env pymol
# Figure 4 Detail: Binding Site with Mutations + Ligand
# FIXED: Use resn QNG for ligand, proper H-bond selection

reinitialize

fetch 6vkv, async=0

remove solvent

hide everything

############################################
# Chain A - protein cartoon
############################################

show cartoon, chain A
color gray85, chain A

set cartoon_fancy_sheets, 1
set cartoon_smooth_segments, 3
set cartoon_flatness, 0.3

############################################
# Ligand - select by residue name QNG, NOT chain
# This is the correct way to select the ligand
############################################

select lig, resn QNG
show sticks, lig
color yellow, lig
set stick_radius, 0.3
set valence, 0

# Show ligand as spheres too for visibility
show spheres, lig
color yellow, lig
set sphere_scale, 0.4

############################################
# Surface - silver-white with transparency
# Show ONLY on chain A
############################################

show surface, chain A
set transparency, 0.65
color gray80, chain A

# Color surface near mutations with mutation colors
select surf57, (chain A and resid 57 around 5 and elem C)
color green, surf57

select surf66, (chain A and resid 66 around 5 and elem C)
color red, surf66

select surf67, (chain A and resid 67 around 5 and elem C)
color orange, surf67

select surf74, (chain A and resid 74 around 5 and elem C)
color blue, surf74

select surf105, (chain A and resid 105 around 5 and elem C)
color cyan, surf105

select surf107, (chain A and resid 107 around 5 and elem C)
color cyan, surf107

############################################
# 6 Mutation sites - colored spheres
# Green=57, Red=66, Orange=67, Blue=74, Cyan=105,107
############################################

# N57 - green
select n57_ca, resid 57 and name CA and chain A
show spheres, n57_ca
color green, n57_ca
set sphere_scale, 0.6, n57_ca

# M66 - red
select m66_ca, resid 66 and name CA and chain A
show spheres, m66_ca
color red, m66_ca
set sphere_scale, 0.6, m66_ca

# Q67 - orange
select q67_ca, resid 67 and name CA and chain A
show spheres, q67_ca
color orange, q67_ca
set sphere_scale, 0.6, q67_ca

# N74 - blue
select n74_ca, resid 74 and name CA and chain A
show spheres, n74_ca
color blue, n74_ca
set sphere_scale, 0.6, n74_ca

# A105 - cyan (compensation region, point 1)
select a105_ca, resid 105 and name CA and chain A
show spheres, a105_ca
color cyan, a105_ca
set sphere_scale, 0.6, a105_ca

# T107 - cyan (compensation region, point 2)
select t107_ca, resid 107 and name CA and chain A
show spheres, t107_ca
color cyan, t107_ca
set sphere_scale, 0.6, t107_ca

############################################
# Hydrogen bond: N57 backbone N to ligand oxygen
# PROMINENT green dashed line
############################################

# First, show N57 backbone N as sphere to ensure it's visible
select n57_n, resid 57 and name N and chain A
show spheres, n57_n
color green, n57_n
set sphere_scale, 0.35, n57_n

# Select ligand oxygens
select lig_o, (resn QNG and elem O)

# Create distance between N57 N and nearest ligand O
distance hbond, n57_n, lig_o, 3.5

# Style the hydrogen bond - make it VERY prominent
set dash_color, green
set dash_gap, 0.0
set dash_width, 10.0
set dash_radius, 0.08

# Hide text label on the distance
hide label, hbond

############################################
# Settings - zoom on binding site
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3
set shininess, 0
set sphere_quality, 2
set stick_quality, 2

# Show everything in a good view
zoom (chain A and resid 50-115) or lig, 8

# Orient to show the binding site well
orient (chain A and resid 50-115) or lig

############################################
# TO SAVE IMAGE:
# png figure4_detail.png, width=2400, height=2400, ray=1
############################################
