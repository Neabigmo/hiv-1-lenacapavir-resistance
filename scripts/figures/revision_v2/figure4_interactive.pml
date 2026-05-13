#!/usr/bin/env pymol
# Figure 4: LEN-CA Hexamer Structure (PDB 6VKV)
# Semi-transparent surface with visible cartoon underneath
# Both layers visible at any zoom level

reinitialize

# Fetch biological assembly
fetch 6vkv, type=pdb1, async=0

remove solvent

hide everything

############################################
# Protein cartoon -条带结构 (shown underneath)
############################################

show cartoon, polymer

# Each chain different color
util.chainbow

# Cartoon settings
set cartoon_fancy_sheets, 1
set cartoon_smooth_segements, 1

############################################
# Protein surface - semi-transparent overlay
############################################

# Create surface on top of cartoon
show surface, polymer

# Silver-gray surface
color gray85, polymer

# Constant transparency - NOT zoom dependent
set surface_mode, 3
set transparency, 0.25, polymer

# Disable depth_cue so transparency is constant
set depth_cue, 0

# Metallic specular
set specular, 0.6
set spec_reflect, 0.3
set shininess, 80
set ambient, 0.3

############################################
# LEN ligand (resn QNG) - light yellow sticks
############################################

select lig, resn QNG
show sticks, lig
color 0xeedd88, lig
set stick_radius, 0.3
set sphere_scale, 0.25, lig

############################################
# Mutation sites - spheres visible through surface
############################################

# N57 - green
select n57, resid 57
show spheres, n57
color 0x66cc66, n57
set sphere_scale, 0.6, n57

# M66 - red
select m66, resid 66
show spheres, m66
color 0xcc6666, m66
set sphere_scale, 0.6, m66

# Q67 - orange
select q67, resid 67
show spheres, q67
color 0xcc9966, q67
set sphere_scale, 0.6, q67

# N74 - blue
select n74, resid 74
show spheres, n74
color 0x6699cc, n74
set sphere_scale, 0.6, n74

# A105 - lighter cyan
select a105, resid 105
show spheres, a105
color 0x99dddd, a105
set sphere_scale, 0.6, a105

# T107 - darker cyan
select t107, resid 107
show spheres, t107
color 0x55aaaa, t107
set sphere_scale, 0.6, t107

############################################
# Settings
############################################

bg_color white
set ray_opaque_background, on
set antialias, 2

# Camera - top down view
orient polymer
turn x, 85
zoom 3

############################################
# TO SAVE IMAGE:
# png figure4.png, width=2400, height=2400, ray=1
############################################