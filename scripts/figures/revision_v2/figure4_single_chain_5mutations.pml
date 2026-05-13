#!/usr/bin/env pymol

reinitialize
fetch 6vkv, async=0
remove solvent
hide everything

# Chain A cartoon
show cartoon, chain A
color gray85, chain A
set cartoon_fancy_sheets, 1

# Surface
show surface, chain A
color gray85, chain A
set transparency, 0.5

# Mutation selections
select mut57, chain A and resid 57
select mut66, chain A and resid 66
select mut67, chain A and resid 67
select mut74, chain A and resid 74
select mut105, chain A and resid 105
select mut107, chain A and resid 107
select n57_n, chain A and resid 57 and name N

# Mutation spheres
show spheres, mut57+mut66+mut67+mut74+mut105+mut107
set sphere_scale, 0.6
color green, mut57
color red, mut66
color orange, mut67
color blue, mut74
color cyan, mut105
color lightblue, mut107

# Surface coloring
color green, (chain A and resid 57 around 4 and elem C)
color red, (chain A and resid 66 around 4 and elem C)
color orange, (chain A and resid 67 around 4 and elem C)
color blue, (chain A and resid 74 around 4 and elem C)
color cyan, (chain A and resid 105 around 4 and elem C)
color lightblue, (chain A and resid 107 around 4 and elem C)

# Binding pocket - chain A ONLY, gray
select pocket, byres (chain A within 5 of (resn QNG and chain A around 6))
show sticks, pocket
set stick_radius, 0.1
color gray85, pocket

# N57 backbone N
color green, n57_n
set sphere_scale, 0.35, n57_n

# H-bond
select lig_o, (resn QNG and elem O and chain A around 6)
distance hb57, n57_n, lig_o, 3.5
set dash_color, green
set dash_gap, 0.0
set dash_width, 10.0
set dash_radius, 0.1
hide label, hb57

# LIGAND - separate, colored YELLOWORANGE
select this_ligand, (resn QNG and chain A around 6)
show sticks, this_ligand
color yelloworange, this_ligand
set stick_radius, 0.3

# Render settings
bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3
set shininess, 0
set sphere_quality, 3
set stick_quality, 3

zoom pocket, 4
orient pocket

# TO SAVE: png fig4.png, width=2400, height=2400, ray=1