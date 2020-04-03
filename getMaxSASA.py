import pymol    # this one might be unnecessary since we specifically need pymol.cmd()
from pymol import cmd  # PyMOL's methods and commands; this is required for the script to use PyMOL's functions.
import csv

## start with a fresh pymol session
cmd.reinitialize

## SASA settings
cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

# List of biological amino acids
aa = ["arg", "his", "lys", "asp", "glu", "ser", "thr", "asn", "gln", "cys", "gly", "pro", "ala", "val", "ile", "leu", "met", "phe", "tyr", "trp"]

## a header for the table of residue sidechain SASA values
header = ["Table of residue sidechain SASA values"]
subheader = ["Residue", "SASA value"]

with open('MaxSASAs.csv', 'w') as file: # write output values into csv row by row, vals in separate columns
    writer = csv.writer(file, delimiter = ',')
    writer.writerow(header)
    writer.writerow(subheader)
    for residue in aa:
        ## do a calculation of the residue when completely free in solution (not in a tripeptide)
        cmd.fragment(residue)  # Summon a disembodied lysine
        cmd.select("Rgroup", "sidechain and " + residue)    # Select only the sidechain (non-backbone atoms)
        maximumSASA = cmd.get_area("Rgroup", 1, 0)  # get solvent exposed area of the residue sidechain
        maximumSASA_str = str(maximumSASA)  #typecast to a string type (easier to deal with if printing into a string)

        ## banish the free lysine (clear the session for loading only the protein of interest)
        cmd.reinitialize

        current_row = [residue, maximumSASA]
        writer.writerow(current_row)


cmd.reinitialize