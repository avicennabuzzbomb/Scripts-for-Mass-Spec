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

## TODO TODO TODO Add arguments and a docstring.

from pymol import cmd  # PyMOL's methods and commands; this is required for the script to use PyMOL's functions.
import pymol    # this one might be unnecessary since we specifically need pymol.cmd()
import re       # methods for string editing
import decimal  # methods for correct rounding
import csv      # methods for handling csv file i/o
import time     # methods for tracking efficiency of the code (CPU time)

## "Stopwatch" starts now
start = time.time()

## start with a fresh pymol session
cmd.reinitialize

## SASA settings
cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU


################################################################################################################################################
## Method for getting each lysine position and using it to make a separate residue selection, and calculate and write SASA to output.
def find(s, ch, thr):
    assert type(thr) == float  #TODO using the threshold, residues can be annotated as exposed or buried with respect to a cutoff value.
    print("Length of fasta string (characters) is " + num_aa)

    header = ["Residue", "Absolute SASA", "Relative SASA", "Exposed or Buried?"]
    threshold = 0.25           #TODO make this a user-input value

    with open('SASA.csv', 'w') as file:                         # correct implementation for writing values into csv, in the same row, vals in separate cells
        writer = csv.writer(file, delimiter = ',')
        writer.writerow(header)

        #for count, ltr in enumerate(s, 1):     # setting [0] of the string = [1], find each 'K' and get its index (aa position). TODO - make the 1 in enumerate(s, 1) a user input value.
        for count, ltr in enumerate(s, 10):     # fasta of 5KSD is truncated in 2 directions, and "start index" for 5ksd is actually 12, not 1.
                                                # Adjusting enumerate index by the difference in amino acids (N-terminal) corrects this. Note, 11 (12 - 1) are missing from the front,
                                                # but python is a 0-index language, therefore make [10] the start index (0 - 10 = 11 index positions).
                                                # TODO given this index problem, this script can be made generalizable by requiring the user to input either a fasta
                                                # of the full length protein and/or the position where the N-terminal portion of the enzyme begins. (ex., user here would
                                                # enter '12' as the first amino acid position, code would subtract 1 to adjust for 0 index, then it would be saved in a variable
                                                # used by the enumerator). TODO TODO TODO This script works very well but could really use an argparser for these inputs.
            if ltr == ch:
                ## Typecasting count (index or aa position) to a string allows the script to use count in a PyMOL selection-expression.
                count = str(count)
                residue = "" + ltr + count
                ## selection algebra for picking one specific residue at a time and performing operation(s) on it.
                ## note that this specific line of code can be universally useful for making discrete selections.
                cmd.select("sele", "resn lys and resi " + count)

                ## SASA and relative SASA are calculated for each K in the string fasta, rounded to 3 decimal places, 
                ## and then typecast to string so it can be printed easily. In the future, make rounding optional.
                sasa = cmd.get_area("sele", 1, 0)
                rel_sasa = sasa / maximumSASA
                burial = "Exposed"
                if rel_sasa <= threshold:
                    burial = "Buried"
                sasa = str(round(sasa, 3))
                rel_sasa = str(round(rel_sasa, 3))
                ## Each SASA printed to console and then to output
                print("\n" + residue + " | " + sasa + " | " + rel_sasa + " | " + burial)
                current_row = [residue, sasa, rel_sasa, burial]
                writer.writerow(current_row)

    return
################################################################################################################################################

## do a calculation of the residue when completely free in solution (not in a tripeptide)
cmd.fragment("lys")  # Summon a disembodied lysine
cmd.select("Rgroup", "sidechain and lys")    # Select only the sidechain (non-backbone atoms)
maximumSASA = cmd.get_area("Rgroup", 1, 0)  #get solvent exposed area of the lysine
maximumSASA_str = str(maximumSASA)  #typecast to a string type (easier to deal with if printing into a string)

## get rid of the free lysine so that only the protein of interest is now considered
cmd.reinitialize

## load the pymol session, or fetch a PDB entry, containing the objects this script will work on.
## When loading:
# cmd.load('5KSD_C-terminal.pse','all') 

## When fetching:
cmd.fetch("5KSD", "Chain A")  #TODO need to figure out how to correctly handle fetched objects; the script currently cannot
# make valid selections and returns a list of nonsense SASA values instead of getting the real calculation.

## Create a selection of all lysines, get the full aa sequence and count all lysines in the selection using this string
cmd.select("Lysines", "resn Lys and Chain A")
fasta = cmd.get_fastastr('Chain A')
clean_fasta = re.sub("[^A-Z]+", "", fasta)     # remove everything except uppercase letters; Python exports fasta with aa's capitalized; meta data is lower cased and contains whitespace.
num_aa = str(len(clean_fasta))
numLys = clean_fasta.count('K')
print("FASTA sequence is: \n" + fasta + "\n\nNumber of lysines in this selection is " + str(numLys))

## Bulk SASA calculations are done and then results are immediately typecast from float to string to allow string concatenation
areaProtein = str(cmd.get_area("Chain A", 1, 0))
areaLysines = str(cmd.get_area("Lysines", 1, 0))

## get a list of lysines and their positions just from the string, fasta; the method call also writes to the output.
ch = 'K'
threshold = 0.25   # in the future, this can be an argument (for user input)

## one row at a time, comma is delimiter (for cells); call method and print row by row
find(clean_fasta, ch, threshold)

## "Stopwatch" stops now; print runtime
stop = time.time()
print("Relative SASA is calculated using absolute SASA of the sidechain of a free lysine = " + maximumSASA_str + ";\nthe threshold for classifying lysines as buried vs. exposed is " + str(threshold))
print("Time (seconds) taken for SASA calculations: " + str(stop - start))