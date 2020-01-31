## Calculate the Solvent-Accessible Surface Area (SASA) of individual residues at a time or as a group in PyMOL,
## then print to an output file.

from pymol import cmd
import pymol

# start with a fresh pymol session
cmd.reinitialize

## test
"""
print("Hello, so far this is NOT working when invoked in PyMOL")
## NOTE - (this works) to run a sript in PyMOL, it needs to be saved in Pymol's working directory.
## NOTE - to find a PyMOL session's working directory, simply use 'pwd' in PyMOL's command line
## NOTE - 'reinitialize' removes all objects from this pymol session's memory (it's a reset as if the program was restarted) 
"""
# load the pymol session containing the objects this script will work on
cmd.load('Full-length_AHA2 mixed_xlinks_FIGURES_11-06-2019.pse','all') 

# SASA settings
cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

# Bulk SASA calculations are done and then results are immediately typecast from float to string to allow string concatenation
areaProtein = str(cmd.get_area("all", 1, 0))
areaMixedLysines = str(cmd.get_area("mixedPlusTrimer", 1, 0))

# Print SASA results to an output text file TODO: expand this to include a list: SASA for each individual lysine!
Output = open('SASA.txt', 'w')
Output.write("Solvent-Accessible Surface Area Calculations in Angstrom:" + "\n\n")
Output.write("Protein's solvent-accessible surface area is " + areaProtein + "\n")
Output.write("Total solvent-accessible surface area of crosslinked lysines is " + areaMixedLysines + "\n")
Output.close()