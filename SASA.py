## Calculate the Solvent-Accessible Surface Area (SASA) of individual residues at a time or as a group in PyMOL,
## then print to an output file.

## This script works by starting a fresh PyMOL session, then importing a lone residue (lysine in this case) and getting that
## residue's total solvent accessible surface area. This result is written to a new output file. Then the PyMol session is 
## reset and the protein pdb file is loaded. To access individual lysines one at a time, a string is made from PyMOL's fasta
## sequence command, and this string is copied and "cleaned" so that a version containing only the amino acids is left.
## Then, the string is used to get the position (index) of each lysine from python's enumerator. Whenever 'K' is found in the
## string, its index is recorded, where the first index of the string is set to '1'. 

"""
## NOTE - To run a sript in PyMOL, it needs to be saved in Pymol's working directory.
## NOTE - to find a PyMOL session's working directory, simply use 'pwd' in PyMOL's command line; use cd to change directories.
## NOTE - 'reinitialize' removes all objects from this pymol session's memory (it's a reset as if the program was restarted) 
"""

from pymol import cmd  # PyMOL's methods and commands; this is required for the script to use PyMOL's functions.
import pymol
import re       # methods for string editing
import decimal  # methods for correct rounding
import csv      # methods for handling csv file i/o

# SASA settings
cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

#open a file writer object to save results
Output = open('SASA.csv', 'w')   #TODO note - csv is preferable to text. Need to implement writing into separate cells.

# start with a fresh pymol session
cmd.reinitialize

# do a calculation of the residue when completely free in solution (not in a tripeptide)
cmd.fragment("lys")
maximumSASA = cmd.get_area("lys", 1, 0)
maximumSASA_str = str(maximumSASA)
Output.write("A free lysine has a maximum SASA of " + maximumSASA_str + "\n")

# get rid of the free lysine so that only the protein of interest is now considered
cmd.reinitialize

# load the pymol session containing the objects this script will work on
cmd.load('5KSD_C-terminal.pse','all') 

# Create a selection of all lysines, get the full aa sequence and count all lysines in the selection using this string
cmd.select("Lysines", "resn Lys and Chain A")
fasta = cmd.get_fastastr('Chain A')
clean_fasta = re.sub("[^A-Z]+", "", fasta)     # remove everything except uppercase letters; Python exports fasta with aa's capitalized; meta data is lower cased and contains whitespace.
num_aa = str(len(clean_fasta))
numLys = clean_fasta.count('K')
print("FASTA sequence is: \n" + fasta + "\n\nNumber of lysines in this selection is " + str(numLys))

## Method for getting each lysine position and using it to make a separate residue selection, and calculate and write SASA to output.
def find(s, ch, thr):
    assert type(thr) == float  #TODO using the threshold, residues can be annotated as exposed or buried with respect to a cutoff value.
    print("Length of fasta string (characters) is " + num_aa)
    Output.write("\n\nResidue Position" + "     " + "Absolute SASA" + "     " + "Relative SASA")
    for count, ltr in enumerate(s, 1):     # setting [0] of the string = [1], find each 'K' and get its index (aa position).
        if ltr == ch:

            # Typecasting count (index or aa position) to a string allows the script to use count in a PyMOL selection-expression.
            count = str(count)

            # selection algebra for picking one specific residue at a time and performing operation(s) on it.
            # note that this specific line of code can be universally useful for making discrete selections.
            cmd.select("sele", "resn lys and resi " + count)

            # SASA and relative SASA are calculated for each K in the string fasta, rounded to 3 decimal places, 
            # and then typecast to string so it can be printed easily. In the future, make rounding optional.
            sasa = cmd.get_area("sele", 1, 0)
            rel_sasa = sasa / maximumSASA

            sasa = str(round(sasa, 3))
            rel_sasa = str(round(rel_sasa, 3))

            # Each SASA printed to console and then to output
            print("\n" + ltr + count + " | " + sasa + " | " + rel_sasa)
            Output.write("\n" + ltr + count + "     " + sasa + "        " + rel_sasa)  
   
    return

# Bulk SASA calculations are done and then results are immediately typecast from float to string to allow string concatenation
areaProtein = str(cmd.get_area("all", 1, 0))
areaLysines = str(cmd.get_area("Lysines", 1, 0))

# Calculate SASA results to an output text file
Output.write("FASTA sequence is: \n" + fasta)
Output.write("\n\nProtein's solvent-accessible surface area is " + areaProtein + "\n")
Output.write("Total SASA of all lysines is " + areaLysines + "\n")
Output.write("\nSASA for each lysine:")

# get a list of lysines and their positions just from the string, fasta; the method call also writes to the output.
ch = 'K'
threshold = 0.25
find(clean_fasta, ch, threshold)
Output.close()
