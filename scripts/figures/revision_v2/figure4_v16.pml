#!/usr/bin/env pymol
# Figure 4: Binding Site Detail - NO SURFACE VERSION

reinitialize
fetch 6vkv, async=0
remove solvent
hide everything

# 1. Chain A cartoon - gray (NO SURFACE)
show cartoon, chain A
color gray85, chain A

# 2. ALL 6 Mutations - colored CA spheres
select m57, resid 57 and name CA
show spheres, m57
color green, m57
set sphere_scale, 0.8

select m66, resid 66 and name CA
show spheres, m66
color red, m66
set sphere_scale, 0.8

select m67, resid 67 and name CA
show spheres, m67
color orange, m67
set sphere_scale, 0.8

select m74, resid 74 and name CA
show spheres, m74
color blue, m74
set sphere_scale, 0.8

select m105, resid 105 and name CA
show spheres, m105
color cyan, m105
set sphere_scale, 0.8

select m107, resid 107 and name CA
show spheres, m107
color lightblue, m107
set sphere_scale, 0.8

# 3. Pocket - chain A residue 55-75 sticks, GRAY
select pocket, chain A and resid 55-75
show sticks, pocket
color gray85, pocket
set stick_radius, 0.12

# 4. Ligand - YELLOW
select ligand, resn QNG
show sticks, ligand
color yelloworange, ligand
set stick_radius, 0.25

# 5. H-bond N57 to ligand
select n57_n, resid 57 and name N
show spheres, n57_n
color green, n57_n
set sphere_scale, 0.6

select lig_o, resn QNG and name O
distance hb1, n57_n, lig_o, 4.0
set dash_color, red
set dash_width, 20
set dash_gap, 0.4
set dash_radius, 0.15
hide label, hb1
show dashes, hb1

# 6. Settings
bg_color white
set ray_opaque_background, on

# 7. Zoom - show ALL mutations
zoom (chain A and resid 57-110), 4

# TO SAVE: png figure4.png, width=2400, height=2400, ray=1