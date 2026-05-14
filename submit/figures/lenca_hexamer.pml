# PyMOL Macro File for LEN-CA Hexamer (6VKV)
# Double-click this file to open in PyMOL

fetch 6VKV

bg_color white

# Hide all, show cartoon
hide everything, all
show cartoon, all
color white, all

# Color chains
color blue, chain A
color green, chain B
color yellow, chain C
color orange, chain D
color purple, chain E
color pink, chain F

# Show surface
show surface, all
set surface_color, white
set transparency, 0.2

# LEN ligand
hide everything, organic
show sticks, organic
color red, organic
set stick_radius, 0.15, organic

# Key resistance residues
select res_n57, resn ASP and resid 57
show spheres, res_n57
color green, res_n57
set sphere_scale, 0.8, res_n57

select res_m66, resn MET and resid 66
show spheres, res_m66
color firebrick, res_m66
set sphere_scale, 0.9, res_m66

select res_q67, resn GLN and resid 67
show spheres, res_q67
color orange, res_q67
set sphere_scale, 0.8, res_q67

select res_n74, resn ASP and resid 74
show spheres, res_n74
color blue, res_n74
set sphere_scale, 0.8, res_n74

select res_a105, resn ALA and resid 105
show spheres, res_a105
color cyan, res_a105
set sphere_scale, 0.7, res_a105

select res_t107, resn THR and resid 107
show spheres, res_t107
color cyan, res_t107
set sphere_scale, 0.7, res_t107

# Camera setup - zoom in on hexamer
zoom all, 100
center organic

# Rotate for good view
rotate z, 20
rotate y, 60
rotate x, 15

# Labels
set label_font_id, 18
set label_size, 24
set label_color, black
set label_bg_color, white
set label_bg_alpha, 0.9

label resn ASP and resid 57, "N57"
label resn MET and resid 66, "M66"
label resn GLN and resid 67, "Q67"
label resn ASP and resid 74, "N74"
label resn ALA and resid 105, "A105"
label resn THR and resid 107, "T107"
label resn LEN, "LEN"

# Ray tracing settings for high quality
set antialias, 2
set ray_trace_mode, 0

# Zoom to fit
zoom
