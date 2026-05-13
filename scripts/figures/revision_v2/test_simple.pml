#!/usr/bin/env pymol
# TEST: 验证配体和口袋选择

reinitialize
fetch 6vkv, async=0
remove solvent
hide everything

# 1. Chain A cartoon - 灰色
show cartoon, chain A
color gray85, chain A

# 2. 突变球
select m57, resid 57 and chain A and name CA
show spheres, m57
color green, m57
set sphere_scale, 0.7

select m66, resid 66 and chain A and name CA
show spheres, m66
color red, m66
set sphere_scale, 0.7

# 3. 口袋 - 只选chain A的原子
select pocket_only_chainA, chain A and resid 55-75
show sticks, pocket_only_chainA
color gray85, pocket_only_chainA
set stick_radius, 0.12

# 4. 配体 - 只选chain D上的QNG (第一个配体)
# 关键：配体在chain D上，不是chain A
select ligand_chainD, chain D and resn QNG
show sticks, ligand_chainD
color yelloworange, ligand_chainD
set stick_radius, 0.25

# 5. 氢键
select n57_n, resid 57 and chain A and name N
show spheres, n57_n
color green, n57_n
set sphere_scale, 0.4

select lig_o, chain D and resn QNG and elem O
distance hb, n57_n, lig_o, 3.5
set dash_color, green
set dash_width, 12
hide label, hb

# 6. 放大
zoom (chain A and resid 50-80), 4

bg_color white
set ray_opaque_background, on