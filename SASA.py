## Calculate the Solvent-Accessible Surface Area (SASA) of individual residues at a time or as a group in PyMOL,
## then print to an output file.

"""
## NOTE - To run a sript in PyMOL, it needs to be saved in Pymol's working directory.
## NOTE - to find a PyMOL session's working directory, simply use 'pwd' in PyMOL's command line; use cd to change directories.
## NOTE - 'reinitialize' removes all objects from this pymol session's memory (it's a reset as if the program was restarted) 
"""

from pymol import cmd  # PyMOL's methods and commands
import pymol
import re       # methods for string editing
import decimal  # methods for correct rounding

# SASA settings
cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
cmd.set('dot_density', 3)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

# start with a fresh pymol session
cmd.reinitialize

# load the pymol session containing the objects this script will work on
cmd.load('5KSD_C-terminal.pse','all') 

# Create a selection of all lysines, get the full aa sequence and count all lysines in the selection using this string
cmd.select("Lysines", "resn Lys and Chain A")
fasta = cmd.get_fastastr('Chain A')
fasta = re.sub("[^A-Z]+", "", fasta)     # remove everything except uppercase letters; Python exports fasta with aa's capitalized; meta data is lower cased and contains whitespace.
num_aa = str(len(fasta))
numLys = fasta.count('K')
print("FASTA sequence is: \n" + fasta + "\n\nNumber of lysines in this selection is " + str(numLys))


def find(s, ch):
    positions = []
    for count, ltr in enumerate(s, 1):     # syntax like  [pos for pos, char in enumerate(s) if char == c]  is a cleaner way...
        if ltr == ch:
            count = str(count)
            cmd.select("sele", "resn lys and resi " + count)     # string concatenation after typecasting count to a string allows the script to do the selections for me.
            sasa = (cmd.get_area("sele", 1, 0))                  # TODO: SASA should be rounded to three decimal places
            thisLys = [ltr + " " + count, sasa]    
            positions.append(thisLys)
    
    return positions

# get a list of lysines and their positions just from the string, fasta.
ch = 'K'
positions = find(fasta, ch)

print("Length of fasta string (characters) is " + num_aa)
print(positions)

# Bulk SASA calculations are done and then results are immediately typecast from float to string to allow string concatenation
areaProtein = str(cmd.get_area("all", 1, 0))
areaLysines = str(cmd.get_area("Lysines", 1, 0))

# Print SASA results to an output text file
Output = open('SASA.txt', 'w')
Output.write("FASTA sequence is: \n" + fasta)
Output.write("\n\nProtein's solvent-accessible surface area is " + areaProtein + "\n")
Output.write("Total solvent-accessible surface area of crosslinked lysines is " + areaLysines + "\n")
Output.write("SASA for each residue:\n\n" + "Here is where I need to print from the list [positions] to this output stream! TODO-implement print list\n")
Output.write("SASA of each lysine:\n")

Output.close()
