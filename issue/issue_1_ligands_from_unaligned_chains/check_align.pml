# Load structures
load /tmp/alignment_check/4HJO.pdb, query_4HJO
load /tmp/alignment_check/5CNN.pdb, target_5CNN

# Align target to query
align target_5CNN, query_4HJO

# Show as cartoon
hide all
show cartoon, query_4HJO
show cartoon, target_5CNN

# Color differently
color cyan, query_4HJO
color orange, target_5CNN

# Show ligands as sticks
select query_lig, query_4HJO and hetatm and not resn HOH
select target_lig, target_5CNN and hetatm and not resn HOH
show sticks, query_lig
show sticks, target_lig

# Highlight ANP specifically
select anp, resn ANP
show spheres, anp
color magenta, anp

# Center view
zoom

# Print alignment info
print "Alignment complete. Check if ANP (magenta spheres) is near the binding site."
