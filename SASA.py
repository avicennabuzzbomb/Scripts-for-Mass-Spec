## Calculate the Solvent-Accessible Surface Area (SASA) of individual residues at a time or as a group in PyMOL,
## then print to an output file.

"""
## NOTE - To run a sript in PyMOL, it needs to be saved in Pymol's working directory.
## NOTE - to find a PyMOL session's working directory, simply use 'pwd' in PyMOL's command line; use cd to change directories.
## NOTE - 'reinitialize' removes all objects from this pymol session's memory (it's a reset as if the program was restarted) 
"""

from pymol import cmd  # PyMOL's cmd and pymol modules carry the nomenclature the python shell needs
import pymol

# SASA settings
cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

###############################################################################################################################
# METHOD: Iterate through residues of interest, calculate SASA one at a time
def eachLys(x):

    SASA_each_lys = []
    #cmd.select("thisLys", "resi x")
    #eachArea = cmd.get_area("thisLys", 1, 0)  # only if this residue is a lysine, get the SASA for this aa position.
    SASA_each_lys = [x] # fill list

    return SASA_each_lys
###############################################################################################################################

# start with a fresh pymol session
cmd.reinitialize

myspace = {'eachLys': eachLys}  # pass this function (eachLys) into the expression namespace, 
                                # instead of evaluating the expression in the global pymol module namespace


# load the pymol session containing the objects this script will work on
cmd.load('5KSD_C-terminal.pse','all') 

# Create a selection of all lysines, get the full aa sequence and count all lysines in the selection using this string
cmd.select("Lysines", "resn Lys and Chain A")
fasta = cmd.get_fastastr('all')
num_aa = len(fasta)
numLys = fasta.count('K')
print("FASTA sequence is: \n" + fasta + "\n\nNumber of lysines in this selection is " + str(numLys))




# iterate through the selection "Lysines" and store in a list; retrieve only unique resi from this list
# Goal - calculate and print SASA one residue at a time but NOT for every atom.

#store iteration directly as a list; from this, resi can be extracted. This sends iterate to eachLys(), which sends back the list of resi.
aggregate = [cmd.iterate('(Lysines)', 'eachLys(resi)', space=myspace)]

# move elements to a new list, exclude duplicates (eliminates connection to the method call)
unique_resi = ["Residue Positions"]
aggSize = len(aggregate)
print("Size of the Aggregate list is " + str(aggSize))
size = len(unique_resi)
agg_index = 0
"""
# Keep only unique lysine positions
while numLys > size - 1:                  # as long as the list of unique positions (minus header) is smaller than number of lysines, continue
    cycles = 0                            # keep track of this while loop
    for i in aggregate:
        if i <= aggSize:                   # values outside the amino acid range are excluded
            if unique_resi[size-1] == aggregate[agg_index]:   # duplicate values are not saved in the new list
                agg_index += 1            # the index of the next value in aggregate must update for either check
                pass
            else:
                agg_index += 1            # the index of the next value in aggregate must update for either check
                unique_resi[size] = i
                size += 1                 # keep track of new list size, so its indices can be used to check new values
    cycles += 1
    if cycles >= 1000:                    # if this loop does not complete on its own, break it
        break

print(size)    # to show that the list can store it. Goal is to extract unique resi and use to get each unique lys' SASA in a loop.
"""

"""###
# Bulk SASA calculations are done and then results are immediately typecast from float to string to allow string concatenation
areaProtein = str(cmd.get_area("all", 1, 0))
areaLysines = str(cmd.get_area("Lysines", 1, 0))

# Print SASA results to an output text file TODO: expand this to include a list: SASA for each individual lysine!
Output = open('SASA.txt', 'w')
Output.write("Solvent-Accessible Surface Area Calculations in Square Angstroms:" + "\n\n")
Output.write("Protein's solvent-accessible surface area is " + areaProtein + "\n")
Output.write("Total solvent-accessible surface area of crosslinked lysines is " + areaLysines + "\n")
Output.close()
"""###