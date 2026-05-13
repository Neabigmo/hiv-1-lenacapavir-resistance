#!/usr/bin/env pymol
# Check if Chain A and Chain B are identical
# If yes, one chain represents all

reinitialize

fetch 6vkv, async=0

remove solvent

# Align chains A and B to check similarity
align chain B, chain A

# Print RMSD - if 0, chains are identical
print "RMSD between chain A and B:", cmd.align("chain B", "chain A")[0]

# If RMSD is ~0, chains are identical
# Let's also check CA atoms only
align chain B and name CA, chain A and name CA
print "RMSD (CA atoms only):", cmd.align("chain B and name CA", "chain A and name CA")[0]

# Show both chains colored differently
show cartoon, chain A+B
color red, chain A
color blue, chain B
set cartoon_fancy_sheets, 1

# Zoom to see both
orient chain A+B
