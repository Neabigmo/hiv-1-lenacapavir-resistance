#!/usr/bin/env pymol
"""
Figure 4: High-quality ray-traced rendering
Run this AFTER adjusting view in figure4_interactive.pml
This will save a 2400x2400 publication-quality image.
"""

# Output settings
png_width = 2400
png_height = 2400
output_file = "figure4_final.png"

# High-quality ray tracing settings
print("Setting up high-quality ray tracing...")
set ray_trace_mode, 0
set antialias, 2
set valence, 0
set ambient_occlusion_mode, 1
set ambient_occlusion_scale, 1.0
set ambient_occlusion_smooth, 3
set ray_shadow, 1

# Ensure view is optimized
zoom all, 0.8

# Ray trace and save
print(f"Rendering at {png_width}x{png_height} pixels...")
png output_file, width=png_width, height=png_height, dpi=300, ray=1

print(f"Done! Image saved as: {output_file}")
