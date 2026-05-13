#!/usr/bin/env pymol
# Figure 4: Binding Site Detail - HIV-1 Capsid LEN-CA with Mutations
# Order: 1.银灰 2.突变 3.药物黄色

reinitialize
fetch 6vkv, async=0
remove solvent
hide everything

############################################
# 1. Chain A - 全部银灰色卡通
############################################
show cartoon, chain A
color gray85, chain A
set cartoon_fancy_sheets, 1

############################################
# 2. Surface - 银灰色半透明
############################################
show surface, chain A
set transparency, 0.6
color gray85, chain A

############################################
# 3. 口袋 - chain A上的残基，银灰色
# 注意：只选chain A，不选配体
############################################
select pocket_residues, chain A and (resn QNG around 6)
show sticks, pocket_residues
color gray85, pocket_residues
set stick_radius, 0.12

############################################
# 4. 突变球 - 六个位点
############################################
# N57 - 绿色
select n57_ca, resid 57 and name CA and chain A
show spheres, n57_ca
color green, n57_ca
set sphere_scale, 0.65

# M66 - 红色
select m66_ca, resid 66 and name CA and chain A
show spheres, m66_ca
color red, m66_ca
set sphere_scale, 0.65

# Q67 - 橙色
select q67_ca, resid 67 and name CA and chain A
show spheres, q67_ca
color orange, q67_ca
set sphere_scale, 0.65

# N74 - 蓝色
select n74_ca, resid 74 and name CA and chain A
show spheres, n74_ca
color blue, n74_ca
set sphere_scale, 0.65

# A105 - 青色
select a105_ca, resid 105 and name CA and chain A
show spheres, a105_ca
color cyan, a105_ca
set sphere_scale, 0.65

# T107 - 浅青色
select t107_ca, resid 107 and name CA and chain A
show spheres, t107_ca
color lightblue, t107_ca
set sphere_scale, 0.65

############################################
# 5. 配体/药物 - 只选QNG残基本身，黄色
# 关键：resn QNG是配体，只在chains D,L,R上
############################################
select this_ligand, resn QNG
show sticks, this_ligand
set stick_radius, 0.25
set valence, 0
color yelloworange, this_ligand

############################################
# 6. N57氢键 - 绿色
############################################
select n57_n, resid 57 and name N and chain A
show spheres, n57_n
color green, n57_n
set sphere_scale, 0.35

select lig_o, resn QNG and elem O
distance hb57, n57_n, lig_o, 3.5
set dash_color, green
set dash_gap, 0.0
set dash_width, 12.0
set dash_radius, 0.08
hide label, hb57

############################################
# 7. 渲染设置
############################################
bg_color white
set ray_opaque_background, on
set antialias, 2
set specular, 0.3
set shininess, 0
set sphere_quality, 2
set stick_quality, 2

############################################
# 8. 放大局部
############################################
zoom (chain A and resid 55-110), 5
orient (chain A and resid 55-110)

############################################
# TO SAVE: png figure4.png, width=2400, height=2400, ray=1
############################################
