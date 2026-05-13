#!/usr/bin/env pymol
# Figure 4: Binding Site Detail - ENHANCED VISIBILITY
# 改进：氢键更粗、配体更突出

reinitialize
fetch 6vkv, async=0
remove solvent
hide everything

# 1. Chain A cartoon - gray
show cartoon, chain A
color gray85, chain A

# 2. Surface - gray with transparency
show surface, chain A
set transparency, 0.7
color gray85, chain A

# 3. ALL 6 Mutations - LARGER colored CA spheres
select m57, resid 57 and name CA
show spheres, m57
color green, m57
set sphere_scale, 1.2

select m66, resid 66 and name CA
show spheres, m66
color red, m66
set sphere_scale, 1.2

select m67, resid 67 and name CA
show spheres, m67
color orange, m67
set sphere_scale, 1.2

select m74, resid 74 and name CA
show spheres, m74
color blue, m74
set sphere_scale, 1.2

select m105, resid 105 and name CA
show spheres, m105
color cyan, m105
set sphere_scale, 1.2

select m107, resid 107 and name CA
show spheres, m107
color lightblue, m107
set sphere_scale, 1.2

# 4. Pocket - chain A residue 55-75 sticks, GRAY
select pocket, chain A and resid 55-75
show sticks, pocket
color gray85, pocket
set stick_radius, 0.15

# 5. Ligand - YELLOW, LARGER STICKS
select ligand, resn QNG
show sticks, ligand
color yelloworange, ligand
set stick_radius, 0.35
util.colorrainbow ligand, _ ligand

# 5b. Ligand surface - yellow transparent for better visibility
show surface, ligand
color yelloworange, ligand
set transparency, 0.5

# 6. H-bond N57 to ligand - ENHANCED RED DASH
select n57_n, resid 57 and name N
show spheres, n57_n
color green, n57_n
set sphere_scale, 0.8

select lig_o, resn QNG and name O
distance hb1, n57_n, lig_o, 3.5
set dash_color, red
set dash_width, 25
set dash_gap, 0.2
set dash_radius, 0.25
set label_distance_digits, 1
set label_position, (0, -1, 0)
hide label, hb1
show dashes, hb1

# 7. Settings
bg_color white
set ray_opaque_background, on

# 8. Zoom - show ALL mutations
zoom (chain A and resid 57-110), 4
orient (chain A and resid 57-110)

# TO SAVE: png figure4.png, width=2400, height=2400, ray=1