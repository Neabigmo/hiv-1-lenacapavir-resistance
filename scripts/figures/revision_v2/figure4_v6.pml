#!/usr/bin/env pymol
# Figure 4: Binding Site Detail
# Order: 1.gray 2.mutations 3.ligand YELLOW

reinitialize
fetch 6vkv, async=0
remove solvent
hide everything

# 1. Chain A cartoon - gray
show cartoon, chain A
color gray85, chain A
set cartoon_fancy_sheets, 1

# 2. Surface - gray, semi-transparent
show surface, chain A
set transparency, 0.6
color gray85, chain A

# 3. Pocket - only chain A atoms, GRAY
select pocket, chain A and resid 55-75
show sticks, pocket
color gray85, pocket
set stick_radius, 0.12

# 4. Mutations - colored spheres
# 57=green, 66=red, 67=orange, 74=blue, 105=cyan, 107=lightblue
select m57, resid 57 and chain A and name CA
show spheres, m57
color green, m57
set sphere_scale, 0.6

select m66, resid 66 and chain A and name CA
show spheres, m66
color red, m66
set sphere_scale, 0.6

select m67, resid 67 and chain A and name CA
show spheres, m67
color orange, m67
set sphere_scale, 0.6

select m74, resid 74 and chain A and name CA
show spheres, m74
color blue, m74
set sphere_scale, 0.6

select m105, resid 105 and chain A and name CA
show spheres, m105
color cyan, m105
set sphere_scale, 0.6

select m107, resid 107 and chain A and name CA
show spheres, m107
color lightblue, m107
set sphere_scale, 0.6

# 5. Ligand - YELLOW (resn QNG is the drug, on chains D/L/R)
select ligand_only, resn QNG
show sticks, ligand_only
set stick_radius, 0.25
set valence, 0
color yelloworange, ligand_only

# 6. H-bond N57 to ligand
select n57_backbone, resid 57 and chain A and name N
show spheres, n57_backbone
color green, n57_backbone
set sphere_scale, 0.3

select ligand_oxy, resn QNG and elem O
distance hb1, n57_backbone, ligand_oxy, 3.5
set dash_color, green
set dash_width, 10
set dash_gap, 0
set dash_radius, 0.06
hide label, hb1

# 7. Settings
bg_color white
set ray_opaque_background, on
set antialias, 2

# 8. Zoom
zoom (chain A and resid 55-110), 5
orient (chain A and resid 55-110)

# TO SAVE: png figure4.png, width=2400, height=2400, ray=1
